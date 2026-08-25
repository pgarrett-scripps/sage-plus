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
use crate::modification::ModificationSpecificity;
use crate::peptide::Peptide;
use crate::spectrum::{select_most_intense_peak, ProcessedSpectrum};

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
/// * `potential_mods` is the search's `(specificity, mass)` list (available as
///   [`crate::database::IndexedDatabase::potential_mods`]); it is used to map a
///   delta mass back onto its candidate residues.
/// * `ion_kinds` should be the same set of fragment ion kinds used for scoring.
pub fn localize(
    peptide: &Peptide,
    spectrum: &ProcessedSpectrum,
    ion_kinds: &[Kind],
    potential_mods: &[(ModificationSpecificity, f32)],
    fragment_tol: Tolerance,
    user_max_fragment_charge: Option<u8>,
    precursor_charge: u8,
) -> Localization {
    let max_charge = max_fragment_charge(user_max_fragment_charge, precursor_charge);

    // Map each variable delta mass -> set of residues that may carry it. Only
    // residue-specificity mods are localizable (terminal mods have a single
    // possible position and are not relocated).
    let mut mods = Vec::new();
    for (mass, residues) in residue_specificities(potential_mods) {
        if let Some(loc) = localize_mass(
            peptide,
            spectrum,
            ion_kinds,
            mass,
            &residues,
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
pub fn has_localizable_modification(
    peptide: &Peptide,
    potential_mods: &[(ModificationSpecificity, f32)],
) -> bool {
    residue_specificities(potential_mods)
        .iter()
        .any(|(mass, residues)| {
            peptide.sequence.iter().enumerate().any(|(index, residue)| {
                residues.contains(residue)
                    && (peptide.modification_at(index) - mass).abs() < MASS_EPS
            })
        })
}

/// Collapse the `(specificity, mass)` list into `(mass, residues)` groups,
/// unioning e.g. `Residue(S)`, `Residue(T)`, `Residue(Y)` that share the
/// phospho delta mass.
fn residue_specificities(potential_mods: &[(ModificationSpecificity, f32)]) -> Vec<(f32, Vec<u8>)> {
    let mut groups: Vec<(f32, Vec<u8>)> = Vec::new();
    for (spec, mass) in potential_mods {
        let residue = match spec {
            ModificationSpecificity::Residue(r) => *r,
            // Terminal-specificity mods are not relocated.
            _ => continue,
        };
        match groups.iter_mut().find(|(m, _)| (m - mass).abs() < MASS_EPS) {
            Some((_, residues)) => {
                if !residues.contains(&residue) {
                    residues.push(residue);
                }
            }
            None => groups.push((*mass, vec![residue])),
        }
    }
    groups
}

fn localize_mass(
    peptide: &Peptide,
    spectrum: &ProcessedSpectrum,
    ion_kinds: &[Kind],
    mass: f32,
    residues: &[u8],
    fragment_tol: Tolerance,
    max_charge: u8,
) -> Option<ModLocalization> {
    // Candidate positions: residues matching the specificity that either
    // already carry this mass or are currently unmodified (so we never displace
    // a *different* modification when relocating this one).
    let candidates: Vec<usize> = peptide
        .sequence
        .iter()
        .enumerate()
        .filter(|(idx, residue)| {
            residues.contains(residue) && {
                let m = peptide.modification_at(*idx);
                m == 0.0 || (m - mass).abs() < MASS_EPS
            }
        })
        .map(|(idx, _)| idx)
        .collect();

    // Number of copies currently placed on a candidate site.
    let k = candidates
        .iter()
        .filter(|&&idx| (peptide.modification_at(idx) - mass).abs() < MASS_EPS)
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
            .filter(|&idx| (peptide.modification_at(idx) - mass).abs() < MASS_EPS)
            .collect();
        return Some(ModLocalization {
            mass,
            label: crate::unimod::label_for(mass),
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
                .map(|&p| SiteScore {
                    position: p,
                    residue: peptide.sequence[p],
                    probability: f32::NAN,
                })
                .collect(),
            all_sites: candidates
                .iter()
                .map(|&p| SiteScore {
                    position: p,
                    residue: peptide.sequence[p],
                    probability: f32::NAN,
                })
                .collect(),
        });
    }

    // Estimate the per-ion random match probability from the spectrum, used as
    // the success probability of the binomial site-determining-ion model.
    let p_random = random_match_probability(spectrum, peptide.monoisotopic, fragment_tol);

    // Use the same number of impossible-site decoy candidates as valid target
    // candidates. Equal target/decoy search spaces make direct competition and
    // dataset-level FLR counting interpretable without a size correction.
    let decoy_candidates = balanced_decoy_candidates(peptide, residues, total_c);
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
    let delta_score = if arrangements.len() >= 2 {
        (arrangements[0].score - arrangements[1].score) as f32
    } else {
        0.0
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
        .map(|(i, &pos)| SiteScore {
            position: pos,
            residue: peptide.sequence[pos],
            probability: marginals[i] as f32,
        })
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
        label: crate::unimod::label_for(mass),
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

/// Clone `peptide` and relocate the target `mass`: clear it from every
/// candidate position, then place it on the chosen positions. The total mass is
/// invariant, so `monoisotopic` does not change.
fn build_variant(peptide: &Peptide, mass: f32, candidates: &[usize], chosen: &[usize]) -> Peptide {
    let mut variant = peptide.clone();
    for &pos in candidates {
        if (variant.modifications[pos] - mass).abs() < MASS_EPS {
            variant.modifications[pos] = 0.0;
        }
    }
    for &pos in chosen {
        variant.modifications[pos] = mass;
    }
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
                    candidates.iter().filter(|&&p| p <= idx).count()
                }
                Kind::X | Kind::Y | Kind::Z => {
                    // suffix [idx + 1, len - 1]
                    candidates.iter().filter(|&&p| p > idx).count()
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
