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
//!
//! Optionally, [`LossRollup`] folds neutral-loss peaks into their parent
//! fragments before step 3 (our own design, off by default; see
//! [`LossRollup`] for how each mode handles the phospho/water collision).

use itertools::Itertools;
use serde::{Deserialize, Serialize};
use smallvec::{smallvec, SmallVec};

use crate::ion_series::{IonSeries, Kind};
use crate::mass::{Tolerance, H2O, NH3};
use crate::modification::{ModificationDefinition, ModificationSpecificity, NeutralLossMode};
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
    /// Neutral-loss roll-up evidence, `None` when roll-up is off.
    /// [`LossRollup::Observed`]: peaks folded into an observed parent.
    /// [`LossRollup::Placement`]: site-determining ion slots of the best
    /// arrangements matched only through a loss form, summed over modifications.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub loss_rollup_peaks: Option<u32>,
}

/// Neutral-loss roll-up before site-determining ions are scored. This is our
/// own design, not a published localization method.
///
/// Generic losses are water and ammonia; modification losses are the
/// `neutral_losses` of the modifications on the peptide. A phosphorylated
/// fragment minus H3PO4 has the mass of the unmodified fragment minus water
/// (79.966 + 18.011 = 97.977 Da), so the same peak can be the loss of either
/// placement. The modes treat that collision differently:
///
/// * `Observed` folds a peak into whichever parent (+H2O, +NH3 or +loss mass)
///   is present in the spectrum, choosing the most intense parent. The folded
///   peak is removed and its intensity added to the parent. Loss peaks never
///   become matches, so a collision peak cannot credit either placement; it
///   only leaves the spectrum, which lowers the random-match rate.
/// * `Placement` scores each arrangement against its own fragments: a
///   site-determining ion slot matches if the intact ion or one of its loss
///   forms is observed. A loss form that has the mass of any form of the same
///   fragment with a different modification count (for phospho, the -H3PO4
///   form of the phospho fragment and the -H2O form of the unmodified one) is
///   dropped from both, so a collision peak credits neither placement. Each
///   slot's chance probability grows with its number of forms, and the score
///   uses the Poisson-binomial tail.
#[derive(
    Serialize, Deserialize, Clone, Copy, Debug, Default, PartialEq, Eq, schemars::JsonSchema,
)]
#[serde(rename_all = "snake_case")]
pub enum LossRollup {
    #[default]
    Off,
    Observed,
    Placement,
}

/// Generic neutral losses rolled up for every fragment.
const GENERIC_LOSSES: [f32; 2] = [H2O, NH3];

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
    localize_with_rollup(
        peptide,
        spectrum,
        ion_kinds,
        potential_mods,
        fragment_tol,
        user_max_fragment_charge,
        precursor_charge,
        LossRollup::Off,
    )
}

/// [`localize`] with an optional neutral-loss roll-up ([`LossRollup`]).
#[allow(clippy::too_many_arguments)]
pub fn localize_with_rollup<R: LocalizationRule>(
    peptide: &Peptide,
    spectrum: &ProcessedSpectrum,
    ion_kinds: &[Kind],
    potential_mods: &[R],
    fragment_tol: Tolerance,
    user_max_fragment_charge: Option<u8>,
    precursor_charge: u8,
    rollup: LossRollup,
) -> Localization {
    let max_charge = max_fragment_charge(user_max_fragment_charge, precursor_charge);
    let groups = modification_groups(potential_mods, peptide);

    let mut rolled = None;
    if rollup == LossRollup::Observed && !groups.is_empty() {
        let mut losses: Vec<f32> = GENERIC_LOSSES.to_vec();
        for applied in peptide.applied_modifications() {
            losses.extend(applied.modification.neutral_losses.iter().copied());
        }
        losses.sort_unstable_by(f32::total_cmp);
        losses.dedup_by(|a, b| (*a - *b).abs() < MASS_EPS);
        rolled = Some(roll_up_observed_losses(
            spectrum,
            &losses,
            fragment_tol,
            max_charge,
        ));
    }
    let (scored_spectrum, folded) = match &rolled {
        Some((rolled, folded)) => (rolled, Some(*folded)),
        None => (spectrum, None),
    };

    // Group physical attachment sites by modification identity.
    let mut mods = Vec::new();
    let mut loss_only = 0u32;
    for group in groups {
        if let Some((loc, loss_matches)) = localize_mass(
            peptide,
            scored_spectrum,
            ion_kinds,
            &group,
            fragment_tol,
            max_charge,
            rollup,
        ) {
            loss_only += loss_matches;
            mods.push(loc);
        }
    }
    let loss_rollup_peaks = match rollup {
        LossRollup::Off => None,
        LossRollup::Observed => folded.filter(|_| !mods.is_empty()),
        LossRollup::Placement => (!mods.is_empty()).then_some(loss_only),
    };
    Localization {
        mods,
        loss_rollup_peaks,
    }
}

/// Placement-independent roll-up ([`LossRollup::Observed`]). Peaks are visited
/// from the highest mass down, so a parent is resolved before its losses. A
/// peak whose mass plus a loss (at fragment charge 1..`max_charge`, in sage's
/// peak convention) lands on an observed, heavier peak is folded into that
/// peak's root parent; when several parents are observed the most intense one
/// wins. Returns the simplified spectrum and the number of folded peaks.
fn roll_up_observed_losses(
    spectrum: &ProcessedSpectrum,
    losses: &[f32],
    fragment_tol: Tolerance,
    max_charge: u8,
) -> (ProcessedSpectrum, u32) {
    let masses = &spectrum.masses;
    let intensities = &spectrum.intensities;
    let n = masses.len();
    let mut root: Vec<usize> = (0..n).collect();
    let mut summed: Vec<f32> = intensities.clone();
    let mut folded = 0u32;
    for i in (0..n).rev() {
        let mut best: Option<usize> = None;
        for &loss in losses {
            for charge in 1..max_charge {
                let target = masses[i] + loss / charge as f32;
                if let Some(j) =
                    select_most_intense_peak(masses, intensities, target, fragment_tol, None)
                {
                    if j > i && best.is_none_or(|b| intensities[j] > intensities[b]) {
                        best = Some(j);
                    }
                }
            }
        }
        if let Some(parent) = best {
            let r = root[parent];
            root[i] = r;
            summed[r] += summed[i];
            folded += 1;
        }
    }
    let kept = |values: &[f32]| -> Vec<f32> {
        (0..n)
            .filter(|&idx| root[idx] == idx)
            .map(|idx| values[idx])
            .collect()
    };
    let mut rolled = spectrum.clone();
    rolled.masses = kept(masses);
    rolled.intensities = kept(&summed);
    if spectrum.mobilities.len() == n {
        rolled.mobilities = kept(&spectrum.mobilities);
    }
    if spectrum.charges.len() == n {
        rolled.charges = (0..n)
            .filter(|&idx| root[idx] == idx)
            .map(|idx| spectrum.charges[idx])
            .collect();
    }
    if spectrum.charge_is_known.len() == n {
        rolled.charge_is_known = (0..n)
            .filter(|&idx| root[idx] == idx)
            .map(|idx| spectrum.charge_is_known[idx])
            .collect();
    }
    (rolled, folded)
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
    rollup: LossRollup,
) -> Option<(ModLocalization, u32)> {
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
        return Some((
            ModLocalization {
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
            },
            0,
        ));
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
        loss_only: u32,
    }

    // (score, matched, trials, slots matched only through a loss form)
    let score_variant = |variant: &Peptide| -> (f64, u32, u32, u32) {
        if rollup == LossRollup::Placement {
            let slots = score_arrangement_with_losses(
                variant,
                spectrum,
                ion_kinds,
                &scoring_candidates,
                k,
                group,
                fragment_tol,
                max_charge,
            );
            let probabilities = slots
                .iter()
                .map(|slot| 1.0 - (1.0 - p_random).powi(slot.forms as i32))
                .collect::<Vec<_>>();
            let matched = slots.iter().filter(|slot| slot.matched).count() as u32;
            let loss_only = slots.iter().filter(|slot| slot.loss_only).count() as u32;
            let pvalue = poisson_binomial_tail(matched, &probabilities);
            (
                -10.0 * pvalue.log10(),
                matched,
                slots.len() as u32,
                loss_only,
            )
        } else {
            let (matched, trials) = score_arrangement(
                variant,
                spectrum,
                ion_kinds,
                &scoring_candidates,
                k,
                fragment_tol,
                max_charge,
            );
            let pvalue = binomial_tail(matched, trials, p_random);
            (-10.0 * pvalue.log10(), matched, trials, 0)
        }
    };

    let mut total_trials = 0u32;
    let mut arrangements: Vec<Arrangement> = Vec::with_capacity(n_arrangements);
    for combo in candidates.iter().copied().combinations(k) {
        let variant = build_variant(peptide, mass, &candidates, &combo);
        let (score, matched, trials, loss_only) = score_variant(&variant);
        total_trials = trials; // identical across arrangements
        arrangements.push(Arrangement {
            sites: combo,
            score,
            matched,
            loss_only,
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
    let best_loss_only = arrangements.first().map(|a| a.loss_only).unwrap_or(0);
    let best_target_score = arrangements.first().map(|a| a.score).unwrap_or(0.0);

    let best_decoy_score = decoy_candidates.as_ref().and_then(|decoys| {
        decoys
            .iter()
            .copied()
            .combinations(k)
            .map(|combo| {
                let variant = build_variant(peptide, mass, &candidates, &combo);
                score_variant(&variant).0
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

    Some((
        ModLocalization {
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
        },
        best_loss_only,
    ))
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

/// One site-determining ion slot (fragment and charge) under
/// [`LossRollup::Placement`].
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
struct RollupSlot {
    /// Intact form plus the loss forms kept after collision removal.
    forms: u32,
    /// Any form is observed.
    matched: bool,
    /// Only a loss form is observed.
    loss_only: bool,
}

/// Whether fragment `series_index` of `kind` contains `site`.
fn fragment_contains(kind: Kind, series_index: usize, site: Site) -> bool {
    match (kind, site) {
        (Kind::A | Kind::B | Kind::C, Site::Nterm) => true,
        (Kind::A | Kind::B | Kind::C, Site::Cterm) => false,
        (Kind::A | Kind::B | Kind::C, Site::Sequence(index)) => index as usize <= series_index,
        (Kind::X | Kind::Y | Kind::Z | Kind::ZDot, Site::Nterm) => false,
        (Kind::X | Kind::Y | Kind::Z | Kind::ZDot, Site::Cterm) => true,
        (Kind::X | Kind::Y | Kind::Z | Kind::ZDot, Site::Sequence(index)) => {
            index as usize > series_index
        }
    }
}

/// Every total neutral loss from `copies` modifications, each given by its
/// loss masses and whether the retained form is kept (`optional`). Always
/// includes 0.0.
fn combine_losses<'a>(copies: impl Iterator<Item = (&'a [f32], bool)>) -> SmallVec<[f32; 8]> {
    let mut totals: SmallVec<[f32; 8]> = smallvec![0.0];
    for (losses, optional) in copies {
        if losses.is_empty() {
            continue;
        }
        let mut next: SmallVec<[f32; 8]> = SmallVec::new();
        for total in &totals {
            if optional {
                next.push(*total);
            }
            for loss in losses {
                next.push(total + loss);
            }
        }
        next.sort_unstable_by(f32::total_cmp);
        next.dedup_by(|a, b| (*a - *b).abs() < MASS_EPS);
        totals = next;
    }
    if !totals.contains(&0.0) {
        totals.insert(0, 0.0);
    }
    totals
}

/// Placement-aware roll-up ([`LossRollup::Placement`]): score `variant`'s own
/// site-determining fragments, counting a slot as matched when the intact ion
/// or one of its loss forms is observed. Loss forms are the generic water and
/// ammonia losses plus the fragment's modification losses. A loss form whose
/// mass equals any form of the same fragment carrying a different number of
/// the localized modification is dropped, so that peak credits no placement.
#[allow(clippy::too_many_arguments)]
fn score_arrangement_with_losses(
    variant: &Peptide,
    spectrum: &ProcessedSpectrum,
    ion_kinds: &[Kind],
    candidates: &[usize],
    k: usize,
    group: &ModificationGroup,
    fragment_tol: Tolerance,
    max_charge: u8,
) -> Vec<RollupSlot> {
    let total_c = candidates.len();
    let group_mass = group.mass;
    let is_group =
        |modification: &ModificationDefinition| (modification.mass - group_mass).abs() < MASS_EPS;
    let loss_spec = |modification: &ModificationDefinition| {
        (
            modification.neutral_losses.to_vec(),
            modification.neutral_loss_mode == NeutralLossMode::Optional,
        )
    };
    let group_losses = variant
        .applied_modifications()
        .find(|applied| is_group(applied.modification))
        .map(|applied| loss_spec(applied.modification))
        .unwrap_or_else(|| (Vec::new(), true));

    let mut slots = Vec::new();
    for kind in ion_kinds {
        for (idx, ion) in IonSeries::new(variant, *kind).enumerate() {
            let c_in_region = match kind {
                Kind::A | Kind::B | Kind::C => candidates
                    .iter()
                    .filter(|&&p| p <= idx || p == variant.sequence.len())
                    .count(),
                Kind::X | Kind::Y | Kind::Z | Kind::ZDot => candidates
                    .iter()
                    .filter(|&&p| {
                        (p > idx && p < variant.sequence.len()) || p == variant.sequence.len() + 1
                    })
                    .count(),
            };
            if !is_site_determining(c_in_region, total_c, k) {
                continue;
            }
            let min_count = k.saturating_sub(total_c - c_in_region);
            let max_count = k.min(c_in_region);

            let mut copies = 0usize;
            let mut pinned: Vec<(Vec<f32>, bool)> = Vec::new();
            for applied in variant.applied_modifications() {
                if !fragment_contains(*kind, idx, applied.site) {
                    continue;
                }
                if is_group(applied.modification) {
                    copies += 1;
                } else if !applied.modification.neutral_losses.is_empty() {
                    pinned.push(loss_spec(applied.modification));
                }
            }
            // Loss masses (0.0 = intact) of this fragment with `count` copies.
            let losses_with =
                |count: usize| -> SmallVec<[f32; 8]> {
                    let mut losses =
                        combine_losses(pinned.iter().map(|(l, o)| (l.as_slice(), *o)).chain(
                            (0..count).map(|_| (group_losses.0.as_slice(), group_losses.1)),
                        ));
                    losses.extend(GENERIC_LOSSES);
                    losses.sort_unstable_by(f32::total_cmp);
                    losses.dedup_by(|a, b| (*a - *b).abs() < MASS_EPS);
                    losses
                };
            let other_masses = (min_count..=max_count)
                .filter(|&count| count != copies)
                .flat_map(|count| {
                    let shift = (count as f32 - copies as f32) * group_mass;
                    losses_with(count)
                        .into_iter()
                        .map(move |loss| ion.monoisotopic_mass + shift - loss)
                })
                .collect::<Vec<_>>();
            let forms = losses_with(copies)
                .into_iter()
                .map(|loss| ion.monoisotopic_mass - loss)
                .enumerate()
                .filter(|&(form, form_mass)| {
                    form == 0
                        || (form_mass > 0.0
                            && !other_masses.iter().any(|&other| {
                                let (lo, hi) = fragment_tol.bounds(form_mass);
                                other >= lo && other <= hi
                            }))
                })
                .map(|(_, form_mass)| form_mass)
                .collect::<SmallVec<[f32; 8]>>();

            for charge in 1..max_charge {
                let observed = |mass: f32| {
                    select_most_intense_peak(
                        &spectrum.masses,
                        &spectrum.intensities,
                        mass / charge as f32,
                        fragment_tol,
                        None,
                    )
                    .is_some()
                };
                let intact = observed(forms[0]);
                let matched = intact || forms[1..].iter().any(|&mass| observed(mass));
                slots.push(RollupSlot {
                    forms: forms.len() as u32,
                    matched,
                    loss_only: matched && !intact,
                });
            }
        }
    }
    slots
}

/// Poisson-binomial upper tail: P(X >= `successes`) for independent trials
/// with the given success probabilities.
fn poisson_binomial_tail(successes: u32, probabilities: &[f64]) -> f64 {
    if probabilities.is_empty() || successes == 0 {
        return 1.0;
    }
    let mut dist = vec![0.0f64; probabilities.len() + 1];
    dist[0] = 1.0;
    for (n, &p) in probabilities.iter().enumerate() {
        for x in (0..=n + 1).rev() {
            let stay = dist[x] * (1.0 - p);
            let step = if x > 0 { dist[x - 1] * p } else { 0.0 };
            dist[x] = stay + step;
        }
    }
    let s = (successes as usize).min(probabilities.len());
    dist[s..].iter().sum::<f64>().clamp(1e-300, 1.0)
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
