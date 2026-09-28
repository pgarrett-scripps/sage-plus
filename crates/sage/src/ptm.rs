//! PTM (post-translational modification) site localization.
//!
//! Sage pre-enumerates every variable-modification combination as a separate
//! [`Peptide`] in the database, so a scored PSM already points at a *single*
//! arrangement of its modifications. This module takes that winning peptide
//! and, for each variable modification it carries, asks the question MaxQuant /
//! MSFragger answer with their "site" reports: **which residue actually carries
//! the modification, and with what confidence?**
//!
//! For each distinct variable-mod delta mass on the peptide we
//!   1. recover the set of candidate residues from the search's
//!      [`ModificationSpecificity`] rules (e.g. all S/T/Y for Phospho),
//!   2. enumerate every way to distribute the `k` copies of that mass across
//!      the candidate sites (keeping all other modifications pinned in place),
//!   3. re-score each arrangement against the experimental spectrum using only
//!      *site-determining ions* — fragments whose mass differs between
//!      arrangements,
//!   4. convert the per-arrangement scores into an **AScore**-style delta
//!      between the two best arrangements and a per-site **localization
//!      probability** (the Andromeda / MaxQuant convention).

use itertools::Itertools;
use serde::Serialize;

use crate::ion_series::{IonSeries, Kind};
use crate::mass::Tolerance;
use crate::modification::{ModificationDefinition, ModificationSpecificity};
use crate::peptide::Peptide;
use crate::peptide::Site;
use crate::spectrum::{select_most_intense_peak, ProcessedSpectrum};
use std::sync::Arc;

/// Two modifications are considered the same delta mass if their masses agree
/// to within this tolerance (modification deltas are stored as `f32`).
const MASS_EPS: f32 = 1e-3;

/// Maximum number of site arrangements to enumerate for a single modification.
/// Peptides with many candidate residues (e.g. long, S/T-rich phosphopeptides)
/// can otherwise generate a combinatorial explosion; when the count exceeds
/// this cap the modification is reported as un-localized.
const MAX_ARRANGEMENTS: usize = 4096;

/// Localization confidence for a single candidate site.
#[derive(Serialize, Clone, Debug, PartialEq)]
pub struct SiteScore {
    pub attachment: crate::ptm_library::Attachment,
    /// 0-based residue index within the peptide.
    pub position: usize,
    /// Amino acid residue at this position.
    pub residue: u8,
    /// Marginal localization probability for this site (0..=1).
    pub probability: f32,
}

/// Localization result for one variable modification (one delta mass) on a
/// peptide.
#[derive(Serialize, Clone, Debug, PartialEq)]
pub struct ModLocalization {
    /// Delta mass that was localized.
    pub mass: f32,
    /// Registered Unimod label for `mass`, if any (e.g. `"Phospho"`).
    pub label: Option<String>,
    /// Number of copies (`k`) of this modification on the peptide.
    pub site_count: usize,
    /// Number of candidate residues the modification could occupy.
    pub candidate_sites: usize,
    /// Number of site-determining theoretical ion slots (the binomial trials).
    pub site_determining_ions: u32,
    /// Number of site-determining ions matched by the best arrangement.
    pub site_determining_matched: u32,
    /// AScore: best arrangement score minus the second-best (0 if unambiguous).
    pub delta_score: f32,
    /// Absolute score difference from target/decoy arrangement competition.
    pub target_decoy_score: f32,
    /// Whether an impossible-site decoy arrangement won the competition.
    pub decoy_winner: bool,
    /// Whether a balanced target/decoy competition could be constructed.
    pub competition_eligible: bool,
    /// Dataset-level false-localization-rate q-value, assigned after all PSMs
    /// have been localized.
    pub localization_q_value: f32,
    /// The `k` highest-probability sites, sorted by position. These are the
    /// "localized" positions reported in the site table.
    pub best_sites: Vec<SiteScore>,
    /// Marginal localization probability for *every* candidate site.
    pub all_sites: Vec<SiteScore>,
}

impl ModLocalization {
    /// A target-decoy separation cannot resolve equally scoring target arrangements.
    /// Keep such arrangements in the competition population but do not accept one
    /// arbitrary site assignment as a confidently localized result.
    pub fn set_competition_q_value(&mut self, q_value: f32) {
        self.localization_q_value = if self.candidate_sites > self.site_count
            && (!self.delta_score.is_finite() || self.delta_score <= 0.0)
        {
            1.0
        } else {
            q_value
        };
    }
}

/// Identity of a modification type for the false-localization-rate
/// competition: the reported name (Unimod label, configured name, or signed
/// mass) and the delta mass on a 0.001 Da grid, the tolerance at which two
/// localized masses are treated as the same modification.
#[derive(Clone, Debug, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct ModificationType {
    pub name: String,
    mass_milli_da: i64,
}

impl ModificationType {
    pub fn new(name: &str, mass: f32) -> Self {
        Self {
            name: name.to_owned(),
            mass_milli_da: (mass as f64 / MASS_EPS as f64).round() as i64,
        }
    }
}

impl ModLocalization {
    /// Name written to the `modification` column of the site reports.
    pub fn reported_name(&self) -> String {
        self.label
            .clone()
            .unwrap_or_else(|| format!("{:+}", self.mass))
    }

    /// Modification type whose competition this localization belongs to.
    pub fn modification_type(&self) -> ModificationType {
        ModificationType::new(&self.reported_name(), self.mass)
    }
}

/// All per-modification localization results for a single PSM.
#[derive(Serialize, Clone, Debug, Default, PartialEq)]
pub struct Localization {
    pub mods: Vec<ModLocalization>,
}

/// Mirror of [`crate::scoring`]'s private `max_fragment_charge`, so the
/// localization search considers the same fragment charge range as scoring.
fn max_fragment_charge(max_fragment_charge: Option<u8>, precursor_charge: u8) -> u8 {
    precursor_charge
        .min(
            max_fragment_charge
                .map(|c| c + 1)
                .unwrap_or(precursor_charge),
        )
        .max(2)
}

/// Localize every variable modification carried by `peptide` against `spectrum`.
///
/// * `potential_mods` accepts full definitions from
///   [`crate::database::IndexedDatabase::localization_mods`] to retain identity.
///   Legacy `(specificity, mass)` rules are also accepted with mass-based identity.
/// * `ion_kinds` should be the same set of fragment ion kinds used for scoring.
pub fn localize<R: LocalizationRule>(
    peptide: &Peptide,
    spectrum: &ProcessedSpectrum,
    ion_kinds: &[Kind],
    potential_mods: &[R],
    fragment_tol: Tolerance,
    user_max_fragment_charge: Option<u8>,
    precursor_charge: u8,
) -> Localization {
    let max_charge = max_fragment_charge(user_max_fragment_charge, precursor_charge);

    // Group physical attachment sites by modification identity.
    let mut mods = Vec::new();
    for group in modification_groups(potential_mods, peptide) {
        if let Some(loc) = localize_mass(
            peptide,
            spectrum,
            ion_kinds,
            &group,
            fragment_tol,
            max_charge,
        ) {
            mods.push(loc);
        }
    }
    Localization { mods }
}

/// Return whether a peptide carries at least one residue-specific variable
/// modification that the localizer can move between candidate sites.
pub fn has_localizable_modification<R: LocalizationRule>(
    peptide: &Peptide,
    potential_mods: &[R],
) -> bool {
    modification_groups(potential_mods, peptide)
        .iter()
        .any(|group| {
            group
                .candidates(peptide)
                .iter()
                .any(|&index| group.is_placed(peptide, index))
        })
}

/// Maximum number of positional isomers scored for one PSM's `isomer_delta`.
/// Beyond it the isomers are not scored and the delta is null.
pub const MAX_POSITIONAL_ISOMERS: usize = MAX_ARRANGEMENTS;

/// The positional isomers of a peptidoform: the same sequence and
/// modification composition with the variable modifications at other sites.
#[derive(Debug)]
pub enum PositionalIsomers {
    /// The peptidoform carries no variable modification from the rules.
    Unmodified,
    /// Its variable modifications have exactly one possible placement.
    SinglePlacement,
    /// More than [`MAX_POSITIONAL_ISOMERS`] placements; none were built.
    TooMany,
    /// Every other placement, excluding `peptide` itself.
    Isomers(Vec<Peptide>),
}

/// Enumerate every positional isomer of `peptide` under `potential_mods`.
///
/// All variable modifications move jointly: each modification type keeps its
/// copy count and may occupy any site its rules allow that no other
/// modification holds. Modifications not described by `potential_mods`
/// (static modifications, search-time mass offsets) stay pinned.
pub fn positional_isomers<R: LocalizationRule>(
    peptide: &Peptide,
    potential_mods: &[R],
    max_isomers: usize,
) -> PositionalIsomers {
    use crate::peptide::AppliedModification;

    let length = peptide.sequence.len();
    let groups = modification_groups(potential_mods, peptide)
        .into_iter()
        .filter_map(|group| {
            let mut sites = group
                .sites
                .iter()
                .map(|&site| encode_site(site, length))
                .collect::<Vec<_>>();
            sites.sort_unstable();
            sites.dedup();
            let placed = sites
                .iter()
                .copied()
                .filter(|&index| group.is_placed(peptide, index))
                .collect::<Vec<_>>();
            (!placed.is_empty()).then_some((group, sites, placed))
        })
        .collect::<Vec<_>>();
    if groups.is_empty() {
        return PositionalIsomers::Unmodified;
    }

    // Split the applied modifications into movable placements and pinned ones.
    let mut pinned = Vec::new();
    let mut templates: Vec<Option<AppliedModification>> = vec![None; groups.len()];
    for applied in peptide.applied_modifications() {
        let index = encode_site(applied.site, length);
        let owner = groups.iter().position(|(group, _, placed)| {
            placed.contains(&index)
                && (applied.modification.mass - group.mass).abs() < MASS_EPS
                && group
                    .definition
                    .as_ref()
                    .is_none_or(|definition| applied.modification == definition.as_ref())
        });
        let owned = AppliedModification {
            site: applied.site,
            modification: Arc::new(applied.modification.clone()),
            kind: applied.kind,
        };
        match owner {
            Some(group) if templates[group].is_none() => templates[group] = Some(owned),
            Some(_) => {}
            None => pinned.push(owned),
        }
    }
    let pinned_sites = pinned
        .iter()
        .map(|applied| encode_site(applied.site, length))
        .collect::<Vec<_>>();

    let mut total = 1usize;
    let mut choices = Vec::with_capacity(groups.len());
    for (_, sites, placed) in &groups {
        let available = sites
            .iter()
            .copied()
            .filter(|index| !pinned_sites.contains(index) || placed.contains(index))
            .collect::<Vec<_>>();
        total = total.saturating_mul(num_combinations(available.len(), placed.len()));
        choices.push(available);
    }
    if total > max_isomers {
        return PositionalIsomers::TooMany;
    }

    let mut isomers = Vec::new();
    for arrangement in groups
        .iter()
        .zip(&choices)
        .map(|((_, _, placed), available)| {
            available.iter().copied().combinations(placed.len())
        })
        .multi_cartesian_product()
    {
        if arrangement
            .iter()
            .zip(&groups)
            .all(|(chosen, (_, _, placed))| chosen == placed)
        {
            continue;
        }
        let mut occupied = arrangement.iter().flatten().copied().collect::<Vec<_>>();
        let count = occupied.len();
        occupied.sort_unstable();
        occupied.dedup();
        if occupied.len() != count {
            continue;
        }
        let mut applied = pinned.clone();
        for (chosen, template) in arrangement.iter().zip(&templates) {
            let template = template.as_ref().expect("every group has a placement");
            applied.extend(chosen.iter().map(|&index| AppliedModification {
                site: decode_site(index, length),
                modification: template.modification.clone(),
                kind: template.kind,
            }));
        }
        isomers.push(peptide.with_applied_modifications(applied));
    }
    if isomers.is_empty() {
        PositionalIsomers::SinglePlacement
    } else {
        PositionalIsomers::Isomers(isomers)
    }
}

/// Whether `a` and `b` are distinct placements of the same sequence and
/// modification composition.
pub fn is_positional_isomer(a: &Peptide, b: &Peptide) -> bool {
    if a.decoy != b.decoy || a.sequence != b.sequence || a.label_channel != b.label_channel {
        return false;
    }
    let mut left = a.applied_modifications().collect::<Vec<_>>();
    let mut right = b.applied_modifications().collect::<Vec<_>>();
    if left.len() != right.len() {
        return false;
    }
    left.sort_unstable();
    right.sort_unstable();
    if left == right {
        return false;
    }
    let mut left = left
        .iter()
        .map(|applied| (applied.modification, applied.kind))
        .collect::<Vec<_>>();
    let mut right = right
        .iter()
        .map(|applied| (applied.modification, applied.kind))
        .collect::<Vec<_>>();
    left.sort_unstable();
    right.sort_unstable();
    left == right
}

/// Full definitions preserve identity. Mass-only rules remain supported for callers
/// that do not retain modification metadata.
pub trait LocalizationRule {
    fn specificity(&self) -> ModificationSpecificity;
    fn mass(&self) -> f32;
    fn definition(&self) -> Option<Arc<ModificationDefinition>>;
    fn sites(&self, peptide: &Peptide) -> Vec<Site> {
        peptide.rule_sites(self.specificity())
    }
}

pub struct ResolvedLocalizationRule {
    pub specificity: ModificationSpecificity,
    pub definition: Arc<ModificationDefinition>,
    pub sites: Vec<Site>,
}

impl LocalizationRule for ResolvedLocalizationRule {
    fn specificity(&self) -> ModificationSpecificity {
        self.specificity
    }
    fn mass(&self) -> f32 {
        self.definition.mass
    }
    fn definition(&self) -> Option<Arc<ModificationDefinition>> {
        Some(self.definition.clone())
    }
    fn sites(&self, _: &Peptide) -> Vec<Site> {
        self.sites.clone()
    }
}

impl LocalizationRule for (ModificationSpecificity, f32) {
    fn specificity(&self) -> ModificationSpecificity {
        self.0
    }
    fn mass(&self) -> f32 {
        self.1
    }
    fn definition(&self) -> Option<Arc<ModificationDefinition>> {
        None
    }
}

impl LocalizationRule for (ModificationSpecificity, Arc<ModificationDefinition>) {
    fn specificity(&self) -> ModificationSpecificity {
        self.0
    }
    fn mass(&self) -> f32 {
        self.1.mass
    }
    fn definition(&self) -> Option<Arc<ModificationDefinition>> {
        Some(self.1.clone())
    }
}

struct ModificationGroup {
    mass: f32,
    definition: Option<Arc<ModificationDefinition>>,
    specificities: Vec<ModificationSpecificity>,
    sites: Vec<Site>,
}

impl ModificationGroup {
    fn is_placed(&self, peptide: &Peptide, index: usize) -> bool {
        let site = decode_site(index, peptide.sequence.len());
        peptide.applied_modifications().any(|applied| {
            applied.site == site
                && (applied.modification.mass - self.mass).abs() < MASS_EPS
                && self
                    .definition
                    .as_ref()
                    .is_none_or(|definition| applied.modification == definition.as_ref())
        })
    }

    fn candidates(&self, peptide: &Peptide) -> Vec<usize> {
        let mut candidates = self
            .sites
            .iter()
            .copied()
            .filter_map(|site| {
                let index = encode_site(site, peptide.sequence.len());
                let occupied = peptide
                    .applied_modifications()
                    .any(|applied| applied.site == site);
                (!occupied || self.is_placed(peptide, index)).then_some(index)
            })
            .collect::<Vec<_>>();
        candidates.sort_unstable();
        candidates.dedup();
        candidates
    }

    fn residues(&self) -> Vec<u8> {
        self.specificities
            .iter()
            .flat_map(|specificity| match specificity {
                ModificationSpecificity::Residue(r)
                | ModificationSpecificity::Internal(r)
                | ModificationSpecificity::PeptideN(Some(r))
                | ModificationSpecificity::PeptideC(Some(r))
                | ModificationSpecificity::ProteinN(Some(r))
                | ModificationSpecificity::ProteinC(Some(r)) => vec![*r],
                ModificationSpecificity::Motif(motif) => motif.site_residues(),
                _ => Vec::new(),
            })
            .collect()
    }

    fn label(&self) -> Option<String> {
        self.definition
            .as_ref()
            .and_then(|definition| definition.name.as_deref().map(str::to_owned))
            .or_else(|| crate::unimod::label_for(self.mass))
    }
}

fn encode_site(site: Site, length: usize) -> usize {
    match site {
        Site::Nterm => length,
        Site::Cterm => length + 1,
        Site::Sequence(i) => i as usize,
    }
}
fn decode_site(index: usize, length: usize) -> Site {
    if index == length {
        Site::Nterm
    } else if index == length + 1 {
        Site::Cterm
    } else {
        Site::Sequence(index as u32)
    }
}
fn site_score(index: usize, peptide: &Peptide, probability: f32) -> SiteScore {
    let site = decode_site(index, peptide.sequence.len());
    let position = match site {
        Site::Nterm => 0,
        Site::Cterm => peptide.sequence.len() - 1,
        Site::Sequence(i) => i as usize,
    };
    SiteScore {
        attachment: crate::ptm_library::Attachment::from_site(site),
        position,
        residue: peptide.sequence[position],
        probability,
    }
}

fn modification_groups<R: LocalizationRule>(
    potential_mods: &[R],
    peptide: &Peptide,
) -> Vec<ModificationGroup> {
    let mut groups: Vec<ModificationGroup> = Vec::new();
    for rule in potential_mods {
        let definition = rule.definition();
        if let Some(group) = groups.iter_mut().find(|group| {
            group.definition == definition && (group.mass - rule.mass()).abs() < MASS_EPS
        }) {
            if !group.specificities.contains(&rule.specificity()) {
                group.specificities.push(rule.specificity());
            }
            group.sites.extend(rule.sites(peptide));
        } else {
            groups.push(ModificationGroup {
                mass: rule.mass(),
                definition,
                specificities: vec![rule.specificity()],
                sites: rule.sites(peptide),
            });
        }
    }
    groups
}

fn localize_mass(
    peptide: &Peptide,
    spectrum: &ProcessedSpectrum,
    ion_kinds: &[Kind],
    group: &ModificationGroup,
    fragment_tol: Tolerance,
    max_charge: u8,
) -> Option<ModLocalization> {
    let mass = group.mass;
    let candidates = group.candidates(peptide);
    let k = candidates
        .iter()
        .filter(|&&index| group.is_placed(peptide, index))
        .count();

    if k == 0 || candidates.is_empty() {
        return None;
    }

    let total_c = candidates.len();
    let n_arrangements = num_combinations(total_c, k);
    if n_arrangements > MAX_ARRANGEMENTS {
        // Too ambiguous to enumerate; report the existing placement with no
        // confidence so the modification still appears in the report.
        let placed: Vec<usize> = candidates
            .iter()
            .copied()
            .filter(|&idx| group.is_placed(peptide, idx))
            .collect();
        return Some(ModLocalization {
            mass,
            label: group.label(),
            site_count: k,
            candidate_sites: total_c,
            site_determining_ions: 0,
            site_determining_matched: 0,
            delta_score: 0.0,
            target_decoy_score: 0.0,
            decoy_winner: true,
            competition_eligible: false,
            localization_q_value: 1.0,
            best_sites: placed
                .iter()
                .map(|&p| site_score(p, peptide, f32::NAN))
                .collect(),
            all_sites: candidates
                .iter()
                .map(|&p| site_score(p, peptide, f32::NAN))
                .collect(),
        });
    }

    // Estimate the per-ion random match probability from the spectrum, used as
    // the success probability of the binomial site-determining-ion model.
    let p_random = random_match_probability(spectrum, peptide.monoisotopic, fragment_tol);

    // Use the same number of impossible-site decoy candidates as valid target
    // candidates. Equal target/decoy search spaces make direct competition and
    // dataset-level FLR counting interpretable without a size correction.
    let decoy_candidates = balanced_decoy_candidates(peptide, &group.residues(), total_c);
    let mut scoring_candidates = candidates.clone();
    if let Some(decoys) = &decoy_candidates {
        scoring_candidates.extend(decoys.iter().copied());
        scoring_candidates.sort_unstable();
    }

    // Score every arrangement.
    struct Arrangement {
        sites: Vec<usize>,
        score: f64,
        matched: u32,
    }

    let mut total_trials = 0u32;
    let mut arrangements: Vec<Arrangement> = Vec::with_capacity(n_arrangements);
    for combo in candidates.iter().copied().combinations(k) {
        let variant = build_variant(peptide, mass, &candidates, &combo);
        let (matched, trials) = score_arrangement(
            &variant,
            spectrum,
            ion_kinds,
            &scoring_candidates,
            k,
            fragment_tol,
            max_charge,
        );
        total_trials = trials; // identical across arrangements
        let pvalue = binomial_tail(matched, trials, p_random);
        let score = -10.0 * pvalue.log10();
        arrangements.push(Arrangement {
            sites: combo,
            score,
            matched,
        });
    }

    // AScore delta: best - second-best arrangement score.
    arrangements.sort_by(|a, b| b.score.total_cmp(&a.score));
    let mut delta_score = if arrangements.len() >= 2 {
        (arrangements[0].score - arrangements[1].score) as f32
    } else {
        0.0
    };
    if let Some(best) = arrangements.first() {
        let signature = |sites: &[usize]| {
            let mut mapped = sites
                .iter()
                .map(|&index| match decode_site(index, peptide.sequence.len()) {
                    Site::Nterm => 0,
                    Site::Cterm => peptide.sequence.len() - 1,
                    Site::Sequence(i) => i as usize,
                })
                .collect::<Vec<_>>();
            mapped.sort_unstable();
            mapped
        };
        if arrangements
            .iter()
            .skip(1)
            .any(|alternative| signature(&alternative.sites) == signature(&best.sites))
        {
            delta_score = 0.0;
        }
    }
    let best_matched = arrangements.first().map(|a| a.matched).unwrap_or(0);
    let best_target_score = arrangements.first().map(|a| a.score).unwrap_or(0.0);

    let best_decoy_score = decoy_candidates.as_ref().and_then(|decoys| {
        decoys
            .iter()
            .copied()
            .combinations(k)
            .map(|combo| {
                let variant = build_variant(peptide, mass, &candidates, &combo);
                let (matched, trials) = score_arrangement(
                    &variant,
                    spectrum,
                    ion_kinds,
                    &scoring_candidates,
                    k,
                    fragment_tol,
                    max_charge,
                );
                let pvalue = binomial_tail(matched, trials, p_random);
                -10.0 * pvalue.log10()
            })
            .max_by(|a, b| a.total_cmp(b))
    });
    let (target_decoy_score, decoy_winner) = match best_decoy_score {
        Some(decoy_score) => (
            (best_target_score - decoy_score).abs() as f32,
            decoy_score >= best_target_score,
        ),
        None => (0.0, true),
    };

    // Per-arrangement posterior weights via a numerically stable softmax over
    // `score / 10 * ln(10)` (equivalent to normalizing 10^(score/10) = 1/pvalue).
    let log_weights: Vec<f64> = arrangements
        .iter()
        .map(|a| a.score / 10.0 * std::f64::consts::LN_10)
        .collect();
    let max_lw = log_weights
        .iter()
        .cloned()
        .fold(f64::NEG_INFINITY, f64::max);
    let denom: f64 = log_weights.iter().map(|lw| (lw - max_lw).exp()).sum();

    // Marginal probability per candidate site: sum the posterior of every
    // arrangement that places the modification on that site.
    let mut marginals: Vec<f64> = vec![0.0; total_c];
    let position_index: std::collections::HashMap<usize, usize> = candidates
        .iter()
        .enumerate()
        .map(|(i, &pos)| (pos, i))
        .collect();
    for (arr, lw) in arrangements.iter().zip(log_weights.iter()) {
        let w = (lw - max_lw).exp() / denom;
        for site in &arr.sites {
            marginals[position_index[site]] += w;
        }
    }

    let all_sites: Vec<SiteScore> = candidates
        .iter()
        .enumerate()
        .map(|(i, &pos)| site_score(pos, peptide, marginals[i] as f32))
        .collect();

    // Choose the `k` sites with the highest marginal probability as the
    // localized positions.
    let mut ranked: Vec<usize> = (0..total_c).collect();
    ranked.sort_by(|&a, &b| marginals[b].total_cmp(&marginals[a]));
    let mut best_sites: Vec<SiteScore> = ranked
        .into_iter()
        .take(k)
        .map(|i| all_sites[i].clone())
        .collect();
    best_sites.sort_by_key(|s| s.position);

    Some(ModLocalization {
        mass,
        label: group.label(),
        site_count: k,
        candidate_sites: total_c,
        site_determining_ions: total_trials,
        site_determining_matched: best_matched,
        delta_score,
        target_decoy_score,
        decoy_winner,
        competition_eligible: decoy_candidates.is_some(),
        localization_q_value: 1.0,
        best_sites,
        all_sites,
    })
}

/// Select evenly distributed, unmodified residues that are impossible target
/// sites. Returning the same number as target candidates balances the two
/// arrangement search spaces.
fn balanced_decoy_candidates(
    peptide: &Peptide,
    target_residues: &[u8],
    target_count: usize,
) -> Option<Vec<usize>> {
    let pool = peptide
        .sequence
        .iter()
        .enumerate()
        .filter(|(idx, residue)| {
            !target_residues.contains(residue) && peptide.modification_at(*idx) == 0.0
        })
        .map(|(idx, _)| idx)
        .collect::<Vec<_>>();
    if pool.len() < target_count || target_count == 0 {
        return None;
    }
    Some(
        (0..target_count)
            .map(|i| pool[(2 * i + 1) * pool.len() / (2 * target_count)])
            .collect(),
    )
}

/// Convert target/decoy localization competitions into monotonic q-values.
/// Input is `(competition_score, decoy_winner)` in caller order.
pub fn target_decoy_q_values(evidence: &[(f32, bool)]) -> Vec<f32> {
    let mut order = (0..evidence.len()).collect::<Vec<_>>();
    order.sort_by(|&a, &b| evidence[b].0.total_cmp(&evidence[a].0));

    let mut targets = 0usize;
    let mut decoys = 0usize;
    let mut prefix_fdr = vec![1.0f32; order.len()];
    let mut start = 0usize;
    while start < order.len() {
        let score = evidence[order[start]].0;
        let mut end = start + 1;
        while end < order.len() && evidence[order[end]].0 == score {
            end += 1;
        }
        for &idx in &order[start..end] {
            if evidence[idx].1 {
                decoys += 1;
            } else {
                targets += 1;
            }
        }
        let fdr = ((decoys + 1) as f32 / targets.max(1) as f32).min(1.0);
        prefix_fdr[start..end].fill(fdr);
        start = end;
    }

    let mut minimum = 1.0f32;
    for fdr in prefix_fdr.iter_mut().rev() {
        minimum = minimum.min(*fdr);
        *fdr = minimum;
    }

    let mut q_values = vec![1.0; evidence.len()];
    for (rank, &original) in order.iter().enumerate() {
        q_values[original] = prefix_fdr[rank];
    }
    q_values
}

/// Localization q-values with a separate target/decoy competition for each
/// modification type.
///
/// `targets` holds `(modification type, competition_score, decoy_winner)` for
/// every target-PSM localization in the population; `probes` holds
/// `(modification type, competition_score)` for decoy-PSM localizations, which
/// read their q-value off the curve of their own type (1.0 when that type has
/// no target population). Both outputs are in caller order.
///
/// A false localization rate is specific to a modification: residue
/// specificity, the number of candidate sites and the chance that the PSM
/// itself is a wrong peptidoform all differ between, for example, phospho
/// S/T/Y and oxidation of a single Met. Pooling them lets one type's decoy
/// wins set another type's q-values. Following LuciPHOr (Fermin et al. 2013,
/// Mol Cell Proteomics; LuciPHOr2, Fermin et al. 2015, Bioinformatics) and the
/// decoy-amino-acid FLR of Ramsbottom et al. 2022 (J Proteome Res), the FLR is
/// estimated per modification type.
///
/// A PSM carrying several modification types (phospho plus oxidation, say)
/// contributes one localization per type, each to its own competition. Each
/// was scored with the other types held at their placed sites, so the rows are
/// the same as in a pooled competition; only who they compete with changes.
pub fn target_decoy_q_values_by_type<K: Eq + std::hash::Hash + Clone>(
    targets: &[(K, f32, bool)],
    probes: &[(K, f32)],
) -> (Vec<f32>, Vec<f32>) {
    let mut groups: std::collections::HashMap<K, (Vec<usize>, Vec<usize>)> =
        std::collections::HashMap::new();
    for (ix, (key, _, _)) in targets.iter().enumerate() {
        groups.entry(key.clone()).or_default().0.push(ix);
    }
    for (ix, (key, _)) in probes.iter().enumerate() {
        groups.entry(key.clone()).or_default().1.push(ix);
    }
    let mut target_q = vec![1.0f32; targets.len()];
    let mut probe_q = vec![1.0f32; probes.len()];
    for (target_ix, probe_ix) in groups.into_values() {
        let evidence = target_ix
            .iter()
            .map(|&ix| (targets[ix].1, targets[ix].2))
            .collect::<Vec<_>>();
        let q_values = target_decoy_q_values(&evidence);
        let scores = probe_ix.iter().map(|&ix| probes[ix].1).collect::<Vec<_>>();
        let probe_values = q_values_at_scores(&evidence, &q_values, &scores);
        for (&ix, q) in target_ix.iter().zip(q_values) {
            target_q[ix] = q;
        }
        for (&ix, q) in probe_ix.iter().zip(probe_values) {
            probe_q[ix] = q;
        }
    }
    (target_q, probe_q)
}

/// Read q-values for new competition scores off an existing q-value curve.
///
/// `evidence` and `q_values` are the population and output of
/// [`target_decoy_q_values`]. Each probe score receives the smallest q-value
/// of any population entry scoring at or below it, or 1.0 when none does, so
/// `q <= cutoff` admits a probe exactly when its score clears the score
/// threshold that the cutoff sets on the population.
pub fn q_values_at_scores(evidence: &[(f32, bool)], q_values: &[f32], probes: &[f32]) -> Vec<f32> {
    let mut curve = evidence
        .iter()
        .zip(q_values)
        .map(|(&(score, _), &q)| (score, q))
        .collect::<Vec<_>>();
    curve.sort_by(|a, b| a.0.total_cmp(&b.0));
    let mut minimum = 1.0f32;
    for point in curve.iter_mut() {
        minimum = minimum.min(point.1);
        point.1 = minimum;
    }
    probes
        .iter()
        .map(|&score| {
            let at_or_below = curve.partition_point(|point| point.0.total_cmp(&score).is_le());
            at_or_below
                .checked_sub(1)
                .map(|ix| curve[ix].1)
                .unwrap_or(1.0)
        })
        .collect()
}

/// Clone `peptide` and relocate the target `mass`: clear it from every
/// candidate position, then place it on the chosen positions. The total mass is
/// invariant, so `monoisotopic` does not change.
fn build_variant(peptide: &Peptide, mass: f32, candidates: &[usize], chosen: &[usize]) -> Peptide {
    let mut variant = peptide.clone();
    variant.relocate_modification_mass(mass, candidates, chosen, MASS_EPS);
    variant
}

/// Count matched site-determining ions for `variant` and the (arrangement
/// independent) number of site-determining theoretical ion slots.
fn score_arrangement(
    variant: &Peptide,
    spectrum: &ProcessedSpectrum,
    ion_kinds: &[Kind],
    candidates: &[usize],
    k: usize,
    fragment_tol: Tolerance,
    max_charge: u8,
) -> (u32, u32) {
    let total_c = candidates.len();
    let mut matched = 0u32;
    let mut trials = 0u32;

    for kind in ion_kinds {
        for (idx, ion) in IonSeries::new(variant, *kind).enumerate() {
            // Candidate positions covered by this fragment.
            let c_in_region = match kind {
                Kind::A | Kind::B | Kind::C => {
                    // prefix [0, idx]
                    candidates
                        .iter()
                        .filter(|&&p| p <= idx || p == variant.sequence.len())
                        .count()
                }
                Kind::X | Kind::Y | Kind::Z | Kind::ZDot => {
                    // suffix [idx + 1, len - 1]
                    candidates
                        .iter()
                        .filter(|&&p| {
                            (p > idx && p < variant.sequence.len())
                                || p == variant.sequence.len() + 1
                        })
                        .count()
                }
            };
            if !is_site_determining(c_in_region, total_c, k) {
                continue;
            }
            for charge in 1..max_charge {
                trials += 1;
                let mz = ion.monoisotopic_mass / charge as f32;
                if select_most_intense_peak(
                    &spectrum.masses,
                    &spectrum.intensities,
                    mz,
                    fragment_tol,
                    None,
                )
                .is_some()
                {
                    matched += 1;
                }
            }
        }
    }
    (matched, trials)
}

/// A fragment covering `c_in_region` of the `total_c` candidate sites is
/// site-determining iff the number of modifications it contains is *not*
/// constant across all `k`-subsets of the candidate sites.
fn is_site_determining(c_in_region: usize, total_c: usize, k: usize) -> bool {
    let max_count = k.min(c_in_region);
    let min_count = k.saturating_sub(total_c - c_in_region);
    max_count != min_count
}

/// Estimate the probability that a single theoretical ion matches an
/// experimental peak by chance: (#peaks) * (tolerance window width) / (m/z
/// range), evaluated at a representative fragment m/z.
fn random_match_probability(
    spectrum: &ProcessedSpectrum,
    peptide_mono: f32,
    fragment_tol: Tolerance,
) -> f64 {
    let n_peaks = spectrum.masses.len();
    if n_peaks == 0 {
        return 1e-3;
    }
    let lo_mass = spectrum.masses.first().copied().unwrap_or(0.0);
    let hi_mass = spectrum.masses.last().copied().unwrap_or(lo_mass + 1.0);
    let range = (hi_mass - lo_mass).max(1.0);

    // Representative fragment m/z ~ half the peptide mass.
    let center = (peptide_mono / 2.0).max(lo_mass);
    let (lo, hi) = fragment_tol.bounds(center);
    let window = (hi - lo).abs().max(1e-4);

    let p = n_peaks as f64 * window as f64 / range as f64;
    p.clamp(1e-6, 0.5)
}

/// Cumulative binomial upper tail: P(X >= `successes`) for `trials` trials with
/// success probability `p`. Computed in log space to avoid overflow.
fn binomial_tail(successes: u32, trials: u32, p: f64) -> f64 {
    if trials == 0 {
        return 1.0;
    }
    let successes = successes.min(trials);
    if successes == 0 {
        return 1.0;
    }
    let ln_p = p.ln();
    let ln_q = (1.0 - p).ln();
    let mut sum = 0.0f64;
    for x in successes..=trials {
        let ln_pmf = ln_choose(trials, x) + x as f64 * ln_p + (trials - x) as f64 * ln_q;
        sum += ln_pmf.exp();
    }
    sum.clamp(1e-300, 1.0)
}

fn ln_choose(n: u32, k: u32) -> f64 {
    ln_factorial(n) - ln_factorial(k) - ln_factorial(n - k)
}

fn ln_factorial(n: u32) -> f64 {
    // Exact within f64 precision; `n` here is bounded by the number of
    // theoretical fragment ions, which is small.
    (1..=n).map(|i| (i as f64).ln()).sum()
}

/// `n choose k`, saturating at `usize::MAX` to avoid overflow when checking the
/// arrangement cap.
fn num_combinations(n: usize, k: usize) -> usize {
    if k > n {
        return 0;
    }
    let k = k.min(n - k);
    let mut result: u128 = 1;
    for i in 0..k {
        result = result.saturating_mul((n - i) as u128) / (i as u128 + 1);
        if result > usize::MAX as u128 {
            return usize::MAX;
        }
    }
    result as usize
}

#[cfg(test)]
#[path = "../tests/unit/ptm.rs"]
mod test;
