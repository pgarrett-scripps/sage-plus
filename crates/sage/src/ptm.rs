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

/// Default fewest matched separating ions by which the best arrangement must
/// beat the runner-up before its sites are accepted as localized. Separating
/// ions are the fragments whose mass differs between those two arrangements;
/// with a margin of 0, arrangements tied on those ions also compete.
pub const DEFAULT_MIN_SEPARATING_MARGIN: u32 = 1;

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
    /// Whether a balanced target/decoy competition could be constructed. False
    /// when every candidate site is occupied: those sites are certain for the
    /// peptide, so they take no part in the false-localization-rate estimate.
    pub competition_eligible: bool,
    /// Matched separating ions of the best arrangement minus those of the
    /// runner-up (0 when there is no runner-up).
    pub separating_margin: i32,
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

    /// Whether this localization enters the false-localization-rate competition:
    /// a balanced decoy competition exists and the best arrangement beats the
    /// runner-up by at least `min_separating_margin` matched separating ions.
    /// Localizations that fail the margin keep q-value 1.0 and do not count as
    /// targets or decoys.
    pub fn competes(&self, min_separating_margin: u32) -> bool {
        self.competition_eligible && self.separating_margin >= min_separating_margin as i32
    }
}

/// Largest gap between the delta masses of two unnamed modifications that
/// belong to one modification type.
pub const TYPE_MASS_TOLERANCE: f32 = 2e-3;

/// Identity of a modification type for the false-localization-rate and
/// site-level FDR competitions.
///
/// A named modification (configured name or Unimod label) is its name,
/// whatever its delta mass. Unnamed, mass-only modifications are grouped by
/// [`modification_types`], which clusters their masses so that masses closer
/// than [`TYPE_MASS_TOLERANCE`] always share a type: there is no fixed grid
/// whose boundaries could split one modification in two.
#[derive(Clone, Debug, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub enum ModificationType {
    Named(String),
    /// Index of the mass cluster among the population's unnamed masses.
    Unnamed(usize),
}

/// Modification types for a population of `(name, delta mass)` identities,
/// in caller order. Every localization or site that competes together must
/// be passed in one call, so that targets and decoys get the same keys.
///
/// Named identities are keyed by name. Unnamed masses are sorted and split
/// wherever two neighbours are more than [`TYPE_MASS_TOLERANCE`] apart
/// (single linkage), so the grouping depends only on the gaps between the
/// masses present.
pub fn modification_types(identities: &[(Option<&str>, f32)]) -> Vec<ModificationType> {
    let mut unnamed = identities
        .iter()
        .enumerate()
        .filter(|(_, (name, _))| name.is_none())
        .map(|(ix, &(_, mass))| (mass, ix))
        .collect::<Vec<_>>();
    unnamed.sort_by(|a, b| a.0.total_cmp(&b.0));
    let mut cluster = vec![0usize; identities.len()];
    let mut current = 0usize;
    for (position, &(mass, ix)) in unnamed.iter().enumerate() {
        if position > 0 && mass - unnamed[position - 1].0 > TYPE_MASS_TOLERANCE {
            current += 1;
        }
        cluster[ix] = current;
    }
    identities
        .iter()
        .enumerate()
        .map(|(ix, (name, _))| match name {
            Some(name) => ModificationType::Named((*name).to_owned()),
            None => ModificationType::Unnamed(cluster[ix]),
        })
        .collect()
}

impl ModLocalization {
    /// Name written to the `modification` column of the site reports.
    pub fn reported_name(&self) -> String {
        self.label
            .clone()
            .unwrap_or_else(|| format!("{:+}", self.mass))
    }

    /// Input to [`modification_types`]: the label (`None` when unnamed) and
    /// the delta mass.
    pub fn type_identity(&self) -> (Option<&str>, f32) {
        (self.label.as_deref(), self.mass)
    }
}

/// All per-modification localization results for a single PSM.
#[derive(Serialize, Clone, Debug, Default, PartialEq)]
pub struct Localization {
    pub mods: Vec<ModLocalization>,
}

impl Localization {
    /// The peptide in ProForma 2.0 notation with every localized modification
    /// moved to its best sites.
    ///
    /// A modification whose localization q-value is at most `q_cutoff` (and
    /// below 1) is written plainly on its best sites. Any other localization
    /// with more candidates than copies is written as a scored position group
    /// (ProForma 2.0, LeDuc et al. 2022): the modification on each best site
    /// with its label and score, `S[Phospho#g1(0.90)]`, and every other
    /// candidate scoring at least 0.005 as `T[#g1(0.09)]`. Scores are the
    /// site probabilities, rounded to two decimals. Groups are not written for
    /// terminal candidates, which ProForma position groups do not cover; those
    /// stay on their best sites.
    ///
    /// A ProForma group label stands for one modification instance with one
    /// preferred location (section 4.4.2), so `k` copies get `k` groups, one
    /// per best site in sequence order: `S[Phospho#g1(0.95)]` and
    /// `S[Phospho#g2(0.90)]`. Every other candidate could hold any copy, so it
    /// carries one tag per group, `T[#g1(0.15)][#g2(0.15)]` (several bracketed
    /// tags on one residue are allowed by section 4.5). Each tag shows the
    /// site's own marginal probability, so the per-site probabilities read the
    /// same from any group and the full mass is written once per copy. A best
    /// site is not listed in the other copies' groups: one residue holds at
    /// most one copy. One copy gives the single-group output unchanged.
    pub fn proforma(&self, peptide: &Peptide, q_cutoff: f32) -> String {
        let length = peptide.sequence.len();
        let encode = |site: &SiteScore| match site.attachment {
            crate::ptm_library::Attachment::PeptideNTerm
            | crate::ptm_library::Attachment::ProteinNTerm => length,
            crate::ptm_library::Attachment::PeptideCTerm
            | crate::ptm_library::Attachment::ProteinCTerm => length + 1,
            crate::ptm_library::Attachment::Residue => site.position,
        };

        let mut variant = peptide.clone();
        // Residue index -> (group label on the modification, extra group tags).
        let mut annotations: std::collections::BTreeMap<usize, (Option<String>, Vec<String>)> =
            std::collections::BTreeMap::new();
        let mut group = 0usize;
        for modification in &self.mods {
            if modification.best_sites.is_empty()
                || modification.candidate_sites <= modification.site_count
            {
                continue;
            }
            let candidates = modification
                .all_sites
                .iter()
                .map(encode)
                .collect::<Vec<_>>();
            let chosen = modification
                .best_sites
                .iter()
                .map(encode)
                .collect::<Vec<_>>();
            variant.relocate_modification_mass(modification.mass, &candidates, &chosen, MASS_EPS);

            let confident = modification.localization_q_value <= q_cutoff
                && modification.localization_q_value < 1.0;
            if confident || candidates.iter().any(|&index| index >= length) {
                continue;
            }
            // One group per copy, numbered by best site in sequence order.
            let mut copies = chosen.clone();
            copies.sort_unstable();
            copies.dedup();
            let first = group + 1;
            group += copies.len();
            for (site, &index) in modification.all_sites.iter().zip(&candidates) {
                let score = format!("({:.2})", site.probability);
                let entry = annotations.entry(index).or_default();
                if let Some(copy) = copies.iter().position(|&c| c == index) {
                    entry.0 = Some(format!("#g{}{score}", first + copy));
                } else if site.probability >= 0.005 {
                    entry
                        .1
                        .extend((first..=group).map(|label| format!("[#g{label}{score}]")));
                }
            }
        }

        let mut out = String::new();
        if let Some(mass) = variant.nterm {
            out.push_str(&variant.modification_tag(Site::Nterm, mass));
            out.push('-');
        }
        for (index, &residue) in variant.sequence.iter().enumerate() {
            out.push(residue as char);
            let mass = variant.modification_at(index);
            let annotation = annotations.get(&index);
            if mass != 0.0 {
                let mut tag = variant.modification_tag(Site::Sequence(index as u32), mass);
                if let Some(label) = annotation.and_then(|a| a.0.as_deref()) {
                    tag.insert_str(tag.len() - 1, label);
                }
                out.push_str(&tag);
            }
            if let Some((_, extra)) = annotation {
                extra.iter().for_each(|tag| out.push_str(tag));
            }
        }
        if let Some(mass) = variant.cterm {
            out.push('-');
            out.push_str(&variant.modification_tag(Site::Cterm, mass));
        }
        out
    }
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
            separating_margin: 0,
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
    let separating_margin = match arrangements.as_slice() {
        [best, runner_up, ..] => {
            let best_variant = build_variant(peptide, mass, &candidates, &best.sites);
            let runner_variant = build_variant(peptide, mass, &candidates, &runner_up.sites);
            let (best_hits, runner_hits) = separating_matches(
                &best_variant,
                &runner_variant,
                spectrum,
                ion_kinds,
                fragment_tol,
                max_charge,
            );
            best_hits as i32 - runner_hits as i32
        }
        _ => 0,
    };
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
        competition_eligible: total_c > k && decoy_candidates.is_some(),
        separating_margin,
        // Every candidate occupied: the placement is fixed by the peptide, and
        // only a wrong peptide (covered by the PSM q-value) could misplace it.
        localization_q_value: if total_c == k { 0.0 } else { 1.0 },
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

/// Matched ions that separate two arrangements of one peptide: for every
/// fragment whose mass differs between `a` and `b`, count matches of `a`'s ion
/// and of `b`'s ion, at the charges [`score_arrangement`] uses.
fn separating_matches(
    a: &Peptide,
    b: &Peptide,
    spectrum: &ProcessedSpectrum,
    ion_kinds: &[Kind],
    fragment_tol: Tolerance,
    max_charge: u8,
) -> (u32, u32) {
    let matches = |mass: f32| {
        (1..max_charge)
            .filter(|&charge| {
                select_most_intense_peak(
                    &spectrum.masses,
                    &spectrum.intensities,
                    mass / charge as f32,
                    fragment_tol,
                    None,
                )
                .is_some()
            })
            .count() as u32
    };
    let (mut a_hits, mut b_hits) = (0u32, 0u32);
    for kind in ion_kinds {
        for (ion_a, ion_b) in IonSeries::new(a, *kind).zip(IonSeries::new(b, *kind)) {
            if (ion_a.monoisotopic_mass - ion_b.monoisotopic_mass).abs() > MASS_EPS {
                a_hits += matches(ion_a.monoisotopic_mass);
                b_hits += matches(ion_b.monoisotopic_mass);
            }
        }
    }
    (a_hits, b_hits)
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
