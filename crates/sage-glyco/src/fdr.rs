//! Two-level FDR for glycopeptides.
//!
//! 1. **Peptide level.** A linear discriminant over peptide features plus
//!    glyco features separates target from decoy peptides; target-decoy
//!    competition gives a PSM-level `peptide_q`. If the fit fails, it is
//!    retried with increasing ridge regularization and finally falls back
//!    to the hyperscore alone, so a run always produces q-values.
//! 2. **Glycan level.** Every explanation of the precursor delta
//!    (composition, isotope error, ammonium adducts) and its decoy twin get
//!    a log-likelihood-ratio score from their per-class Y-ion and oxonium
//!    counts (see [`crate::evidence`]). The best of all targets and decoys
//!    wins (exact ties by a fair coin), and target-decoy competition on the
//!    winner's score plus its lead over the other label gives `glycan_q`. This follows PTM-Shepherd's glycan assignment for
//!    MSFragger-Glyco, with three changes: per-spectrum random match rates
//!    instead of a global one, ion classes whose match rates are refitted
//!    from confident PSMs (one EM step), and a decoy twin per candidate
//!    explanation, so competition happens between compositions of the same
//!    mass rather than against a global decoy list.

use sage_core::ml::linear_discriminant::LinearDiscriminantAnalysis;

use crate::evidence::{
    ClassCounts, CLASSES, NEUAC_ABSENT, NEUAC_PRESENT, NEUGC_ABSENT, NEUGC_PRESENT,
};
use crate::search::{Explanation, GlycoCandidate};

/// Initial per-class match probability of a true explanation's ions.
const INITIAL_RATES: [f64; CLASSES] = [0.6, 0.4, 0.3, 0.2, 0.8, 0.1, 0.7, 0.1];
/// Pseudo-count weight of the initial rates in the EM refit.
const PRIOR_WEIGHT: f64 = 20.0;
/// Rates are kept inside this range so every log term stays finite.
const RATE_BOUNDS: (f64, f64) = (0.005, 0.995);

/// Glycan-level score model.
#[derive(Clone, Debug, PartialEq)]
pub struct GlycanModel {
    /// Match probability of each ion class for a correct explanation.
    pub rates: [f64; CLASSES],
    /// Log prior of isotope errors, indexed from `isotope_min`.
    pub isotope_prior: Vec<f64>,
    pub isotope_min: i8,
    /// Log prior of the number of ammonium adducts.
    pub adduct_prior: Vec<f64>,
    /// Standard deviation of the precursor mass error, ppm.
    pub ppm_sigma: f64,
}

impl GlycanModel {
    pub fn new(isotopes: (i8, i8), max_adducts: u8, explain_ppm: f32) -> Self {
        let isotope_prior = (isotopes.0..=isotopes.1)
            .map(|isotope| match isotope {
                0 => 0.7f64.ln(),
                1 => 0.25f64.ln(),
                _ => 0.05f64.ln(),
            })
            .collect();
        let adduct_prior = (0..=max_adducts)
            .map(|adducts| match adducts {
                0 => 0.75f64.ln(),
                1 => 0.2f64.ln(),
                _ => 0.05f64.ln(),
            })
            .collect();
        GlycanModel {
            rates: INITIAL_RATES,
            isotope_prior,
            isotope_min: isotopes.0,
            adduct_prior,
            ppm_sigma: (explain_ppm as f64 / 2.0).max(1.0),
        }
    }

    /// Log-likelihood ratio of `counts` under "this explanation is right"
    /// against "its ions match at random", plus mass and prior terms.
    pub fn score(
        &self,
        explanation: &Explanation,
        counts: &ClassCounts,
        random_y: f32,
        random_oxonium: f32,
    ) -> f64 {
        let mut score = 0.0;
        for class in 0..CLASSES {
            let generated = counts.generated[class] as f64;
            if generated == 0.0 {
                continue;
            }
            let matched = counts.matched[class] as f64;
            let random = match class {
                NEUAC_PRESENT | NEUAC_ABSENT | NEUGC_PRESENT | NEUGC_ABSENT => random_oxonium,
                _ => random_y,
            } as f64;
            let random = random.clamp(RATE_BOUNDS.0, RATE_BOUNDS.1);
            let p = self.rates[class];
            score += matched * (p / random).ln()
                + (generated - matched) * ((1.0 - p) / (1.0 - random)).ln();
        }
        let isotope = (explanation.isotope - self.isotope_min) as usize;
        score += self.isotope_prior.get(isotope).copied().unwrap_or(-5.0);
        score += self
            .adduct_prior
            .get(explanation.adducts as usize)
            .copied()
            .unwrap_or(-5.0);
        let z = explanation.error_ppm as f64 / self.ppm_sigma;
        score - 0.5 * z * z
    }

    /// The best explanation of a candidate over targets and decoy twins.
    ///
    /// An exact tie between the best target and the best decoy means the
    /// fragments cannot tell them apart (typically no ion beyond the shared
    /// core matched), so it is settled by a coin flip seeded from the
    /// spectrum. Awarding ties to the target would let every uninformative
    /// explanation pass as a target and hide its error from target-decoy
    /// competition.
    pub fn assign(&self, candidate: &GlycoCandidate) -> Assignment {
        let mut best = [(f64::NEG_INFINITY, 0usize); 2];
        let mut second = f64::NEG_INFINITY;
        for (index, explanation) in candidate.explanations.iter().enumerate() {
            for (slot, counts) in [&explanation.target, &explanation.decoy]
                .into_iter()
                .enumerate()
            {
                let score = self.score(
                    explanation,
                    counts,
                    candidate.random_y,
                    candidate.random_oxonium,
                );
                if score > best[slot].0 {
                    second = second.max(best[slot].0);
                    best[slot] = (score, index);
                } else {
                    second = second.max(score);
                }
            }
        }
        let decoy = match best[0].0.total_cmp(&best[1].0) {
            std::cmp::Ordering::Greater => false,
            std::cmp::Ordering::Less => true,
            std::cmp::Ordering::Equal => coin(&candidate.feature.spec_id, best[0].1),
        };
        let (winner, loser) = if decoy {
            (best[1], best[0])
        } else {
            (best[0], best[1])
        };
        Assignment {
            explanation: winner.1,
            decoy,
            score: winner.0,
            runner_up: second.max(loser.0),
            opposite: loser.0,
        }
    }

    /// One EM step: re-estimate class rates, isotope and adduct priors and
    /// the mass error spread from confidently assigned explanations.
    pub fn refit<'a>(&mut self, confident: impl Iterator<Item = &'a Explanation>) {
        let mut generated = [0.0f64; CLASSES];
        let mut matched = [0.0f64; CLASSES];
        let mut isotopes = vec![1.0f64; self.isotope_prior.len()];
        let mut adducts = vec![1.0f64; self.adduct_prior.len()];
        let mut ppm = Vec::new();
        for explanation in confident {
            for class in 0..CLASSES {
                generated[class] += explanation.target.generated[class] as f64;
                matched[class] += explanation.target.matched[class] as f64;
            }
            if let Some(n) = isotopes.get_mut((explanation.isotope - self.isotope_min) as usize) {
                *n += 1.0;
            }
            if let Some(n) = adducts.get_mut(explanation.adducts as usize) {
                *n += 1.0;
            }
            ppm.push(explanation.error_ppm as f64);
        }
        if ppm.len() < 50 {
            log::info!(
                "glyco: {} confident glycan assignments, keeping the default model",
                ppm.len()
            );
            return;
        }
        for class in 0..CLASSES {
            let prior = INITIAL_RATES[class];
            self.rates[class] = ((PRIOR_WEIGHT * prior + matched[class])
                / (PRIOR_WEIGHT + generated[class]))
                .clamp(RATE_BOUNDS.0, RATE_BOUNDS.1);
        }
        let total: f64 = isotopes.iter().sum();
        self.isotope_prior = isotopes.iter().map(|n| (n / total).ln()).collect();
        let total: f64 = adducts.iter().sum();
        self.adduct_prior = adducts.iter().map(|n| (n / total).ln()).collect();
        let mean = ppm.iter().sum::<f64>() / ppm.len() as f64;
        let var = ppm.iter().map(|x| (x - mean).powi(2)).sum::<f64>() / ppm.len() as f64;
        self.ppm_sigma = var.sqrt().max(1.0);
    }
}

impl Assignment {
    /// Target-decoy competition statistic: the winning score plus its lead
    /// over the best explanation of the other label. The score measures how
    /// well the precursor delta and fragments fit any glycan; the lead is
    /// the evidence for this composition over its decoy twin. Both are
    /// symmetric under swapping a wrong composition with its twin, so
    /// competition stays fair. On the yeast entrapment the sum passes 10%
    /// more glycoPSMs than the score alone (781 against 708), with 2.4%
    /// against 1.3% non-yeast glycans, both below MSFragger-Glyco's 3.8%.
    pub fn competition(&self) -> f64 {
        let lead = self.score - self.opposite;
        // A candidate without a decoy twin (impossible in practice) counts
        // as a large lead.
        self.score + if lead.is_finite() { lead } else { 1e3 }
    }
}

/// Deterministic fair coin from a spectrum identifier.
fn coin(spec_id: &str, salt: usize) -> bool {
    let mut hash = 0xcbf2_9ce4_8422_2325u64 ^ salt as u64;
    for byte in spec_id.bytes() {
        hash = (hash ^ byte as u64).wrapping_mul(0x0100_0000_01b3);
    }
    // Finalizer, so the low bit depends on every input byte.
    hash ^= hash >> 33;
    hash = hash.wrapping_mul(0xff51_afd7_ed55_8ccd);
    hash ^= hash >> 33;
    hash & 1 == 1
}

/// The winning explanation of one candidate.
#[derive(Copy, Clone, Debug, PartialEq)]
pub struct Assignment {
    pub explanation: usize,
    pub decoy: bool,
    pub score: f64,
    /// Best score among the other target and decoy explanations.
    pub runner_up: f64,
    /// Best score of the other label (the best decoy for a target winner).
    pub opposite: f64,
}

/// A scored glycopeptide-spectrum match.
#[derive(Clone, Debug)]
pub struct GlycoPsm {
    pub candidate: GlycoCandidate,
    pub discriminant: f64,
    pub peptide_q: f32,
    pub assignment: Assignment,
    pub glycan_q: f32,
    pub passes_fdr: bool,
}

/// FDR summary numbers.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct FdrSummary {
    /// Spectra with a selected candidate.
    pub candidates: usize,
    /// Explained peptide candidates before selection.
    pub explained: usize,
    /// How the peptide discriminant was obtained.
    pub model: &'static str,
    pub peptide_passing: usize,
    pub passing: usize,
    pub glycan_decoys_passing_peptide: usize,
    pub rates: [f64; CLASSES],
}

const PEPTIDE_FEATURES: usize = 32;

/// One spectrum's selected candidate: its index, and the lead of its glycan
/// score over the next explained peptide candidate of the same spectrum
/// (over zero when it is the only one).
#[derive(Copy, Clone, Debug, PartialEq)]
pub struct Selection {
    pub index: usize,
    pub lead: f64,
}

/// Pick, for each spectrum, the explained peptide candidate whose best
/// glycan explanation scores highest. Candidates of one spectrum are
/// contiguous (see [`crate::GlycoSearch::search`]); ties keep the better
/// peptide rank. Target and decoy peptides are treated alike, so peptide
/// target-decoy competition stays fair.
pub fn select(candidates: &[GlycoCandidate], assignments: &[Assignment]) -> Vec<Selection> {
    let mut selected = Vec::new();
    let mut start = 0;
    while start < candidates.len() {
        let same = |c: &GlycoCandidate| {
            c.feature.file_id == candidates[start].feature.file_id
                && c.feature.spec_id == candidates[start].feature.spec_id
        };
        let mut end = start + 1;
        while end < candidates.len() && same(&candidates[end]) {
            end += 1;
        }
        let score = |i: usize| assignments[i].score;
        // First maximum, so ties keep the better peptide rank.
        let best = (start + 1..end).fold(start, |b, i| if score(i) > score(b) { i } else { b });
        let second = (start..end)
            .filter(|&i| i != best)
            .map(score)
            .fold(None, |m: Option<f64>, s| Some(m.map_or(s, |m| m.max(s))))
            .unwrap_or(0.0);
        let lead = score(best) - second;
        selected.push(Selection {
            index: best,
            lead: if lead.is_finite() { lead } else { 0.0 },
        });
        start = end;
    }
    selected
}

fn peptide_features(
    candidate: &GlycoCandidate,
    best_glycan: f64,
    lead: f64,
) -> [f64; PEPTIDE_FEATURES] {
    let feature = &candidate.feature;
    let best = candidate
        .explanations
        .iter()
        .map(|e| e.target.y_matched() as f64 / e.target.y_generated().max(1) as f64)
        .fold(0.0, f64::max);
    let core = candidate
        .explanations
        .iter()
        .map(|e| e.target.matched[0] as f64 / e.target.generated[0].max(1) as f64)
        .fold(0.0, f64::max);
    let y_intensity = candidate
        .explanations
        .iter()
        .map(|e| e.y_intensity as f64)
        .fold(0.0, f64::max);
    let y_matched = candidate
        .explanations
        .iter()
        .map(|e| e.target.y_matched() as f64)
        .fold(0.0, f64::max);
    let ppm = candidate
        .explanations
        .iter()
        .map(|e| e.error_ppm.abs() as f64)
        .fold(f64::INFINITY, f64::min);
    [
        feature.hyperscore.ln_1p(),
        feature.delta_next.ln_1p(),
        feature.delta_best.ln_1p(),
        (-feature.poisson).max(0.0).ln_1p(),
        feature.matched_intensity_pct as f64,
        (feature.matched_peaks as f64).ln_1p(),
        feature.longest_b as f64,
        feature.longest_y as f64,
        feature.peptide_len as f64,
        feature.missed_cleavages as f64,
        feature.charge as f64,
        core,
        best,
        f64::from(u8::from(candidate.anchored)),
        y_intensity,
        candidate.oxonium_ions as f64,
        candidate.oxonium_fraction as f64,
        best_glycan.clamp(-50.0, 200.0),
        ppm.min(100.0),
        (candidate.explanations.len() as f64).ln(),
        (feature.rank as f64).ln(),
        lead.clamp(-50.0, 50.0),
        y_matched.ln_1p(),
        (candidate.explanations[0].target.matched[0] as f64).ln_1p(),
        f64::from(u8::from(candidate.y1)),
        (candidate.hex_ratio as f64).clamp(-10.0, 10.0),
        (candidate.site.matched as f64).ln_1p(),
        candidate.site.matched as f64 / candidate.site.possible.max(1) as f64,
        (candidate.site.bare as f64).ln_1p(),
        (candidate.site.high_charge as f64).ln_1p(),
        f64::from(u8::from(candidate.twin_hyperscore.is_some())),
        candidate.twin_hyperscore.map_or(0.0, |twin| {
            feature.hyperscore.ln_1p() - twin.max(0.0).ln_1p()
        }),
    ]
}

/// Target-decoy q-values of items sorted by descending score. `decoy[i]`
/// marks decoys; tied scores share one q-value. Returns q in input order.
pub fn q_values(scores: &[f64], decoy: &[bool]) -> Vec<f32> {
    let mut order: Vec<usize> = (0..scores.len()).collect();
    order.sort_by(|&a, &b| scores[b].total_cmp(&scores[a]));
    let mut q = vec![1.0f32; scores.len()];
    let (mut targets, mut decoys) = (0usize, 1usize);
    let mut start = 0;
    while start < order.len() {
        let mut end = start + 1;
        while end < order.len() && scores[order[end]] == scores[order[start]] {
            end += 1;
        }
        for &index in &order[start..end] {
            if decoy[index] {
                decoys += 1;
            } else {
                targets += 1;
            }
        }
        let fdr = decoys as f32 / targets.max(1) as f32;
        for &index in &order[start..end] {
            q[index] = fdr;
        }
        start = end;
    }
    let mut min = 1.0f32;
    for &index in order.iter().rev() {
        min = min.min(q[index]);
        q[index] = min;
    }
    q
}

/// Fit the peptide discriminant: LDA over standardized features, retried
/// with ridge regularization, then hyperscore as the last resort.
fn discriminant(
    rows: &[[f64; PEPTIDE_FEATURES]],
    decoy: &[bool],
    hyperscore: &[f64],
) -> (Vec<f64>, &'static str) {
    let n = rows.len().max(1) as f64;
    let mut mean = [0.0; PEPTIDE_FEATURES];
    let mut sd = [0.0; PEPTIDE_FEATURES];
    for row in rows {
        for j in 0..PEPTIDE_FEATURES {
            mean[j] += row[j] / n;
        }
    }
    for row in rows {
        for j in 0..PEPTIDE_FEATURES {
            sd[j] += (row[j] - mean[j]).powi(2) / n;
        }
    }
    // Constant features carry no information and would make the scatter
    // matrix singular; standardize them to zero.
    let scale: Vec<f64> = sd
        .iter()
        .map(|v| if v.sqrt() > 1e-9 { 1.0 / v.sqrt() } else { 0.0 })
        .collect();
    let standardize = |row: &[f64; PEPTIDE_FEATURES]| -> [f64; PEPTIDE_FEATURES] {
        std::array::from_fn(|j| (row[j] - mean[j]) * scale[j])
    };
    for (regularization, name) in [
        (0.0, "lda"),
        (1e-3, "lda_ridge_1e-3"),
        (1e-2, "lda_ridge_1e-2"),
        (1e-1, "lda_ridge_1e-1"),
    ] {
        let Some(lda) =
            LinearDiscriminantAnalysis::train_regularized(rows, decoy, standardize, regularization)
        else {
            continue;
        };
        let scores: Vec<f64> = rows
            .iter()
            .map(|row| lda.score(&standardize(row)))
            .collect();
        if scores.iter().all(|s| s.is_finite()) {
            return (scores, name);
        }
    }
    log::warn!("glyco: peptide LDA failed, ranking by hyperscore");
    (hyperscore.iter().map(|h| h.ln_1p()).collect(), "hyperscore")
}

/// Score candidates at both levels. `explain_ppm` sets the initial mass
/// error spread, `peptide_fdr` and `glycan_fdr` the reporting thresholds.
pub fn score(
    candidates: Vec<GlycoCandidate>,
    isotopes: (i8, i8),
    max_adducts: u8,
    explain_ppm: f32,
    peptide_fdr: f32,
    glycan_fdr: f32,
) -> (Vec<GlycoPsm>, GlycanModel, FdrSummary) {
    let mut model = GlycanModel::new(isotopes, max_adducts, explain_ppm);
    let explained = candidates.len();
    let assignments: Vec<Assignment> = candidates.iter().map(|c| model.assign(c)).collect();

    // One candidate per spectrum: the best-supported explained peptide.
    let selection = select(&candidates, &assignments);
    let mut keep = vec![false; candidates.len()];
    for s in &selection {
        keep[s.index] = true;
    }
    let assignments: Vec<Assignment> = selection.iter().map(|s| assignments[s.index]).collect();
    let leads: Vec<f64> = selection.iter().map(|s| s.lead).collect();
    let candidates: Vec<GlycoCandidate> = candidates
        .into_iter()
        .zip(keep)
        .filter_map(|(c, keep)| keep.then_some(c))
        .collect();

    // Peptide level.
    let rows: Vec<_> = candidates
        .iter()
        .zip(&assignments)
        .zip(&leads)
        .map(|((candidate, assignment), lead)| peptide_features(candidate, assignment.score, *lead))
        .collect();
    let peptide_decoy: Vec<bool> = candidates.iter().map(|c| c.feature.label == -1).collect();
    let hyperscore: Vec<f64> = candidates.iter().map(|c| c.feature.hyperscore).collect();
    let (discriminant, model_name) = discriminant(&rows, &peptide_decoy, &hyperscore);
    drop(rows);
    let peptide_q = q_values(&discriminant, &peptide_decoy);

    // Glycan level: refit on confident target peptides whose winner is a
    // target explanation clearly ahead of every alternative.
    let confident = candidates
        .iter()
        .zip(&assignments)
        .zip(&peptide_q)
        .filter(|((c, a), q)| {
            c.feature.label != -1 && **q <= peptide_fdr && !a.decoy && a.score - a.runner_up > 2.0
        })
        .map(|((c, a), _)| &c.explanations[a.explanation]);
    model.refit(confident);
    let assignments: Vec<Assignment> = candidates.iter().map(|c| model.assign(c)).collect();

    let eligible: Vec<usize> = (0..candidates.len())
        .filter(|&i| !peptide_decoy[i] && peptide_q[i] <= peptide_fdr)
        .collect();
    let glycan_decoy: Vec<bool> = eligible.iter().map(|&i| assignments[i].decoy).collect();
    let mut glycan_q = vec![1.0f32; candidates.len()];
    let stat: Vec<f64> = eligible
        .iter()
        .map(|&i| assignments[i].competition())
        .collect();
    for (&index, q) in eligible.iter().zip(q_values(&stat, &glycan_decoy)) {
        glycan_q[index] = q;
    }

    let mut summary = FdrSummary {
        candidates: candidates.len(),
        explained,
        model: model_name,
        peptide_passing: eligible.len(),
        glycan_decoys_passing_peptide: glycan_decoy.iter().filter(|d| **d).count(),
        rates: model.rates,
        ..Default::default()
    };
    let psms: Vec<GlycoPsm> = candidates
        .into_iter()
        .enumerate()
        .map(|(i, candidate)| {
            let passes_fdr = !peptide_decoy[i]
                && peptide_q[i] <= peptide_fdr
                && !assignments[i].decoy
                && glycan_q[i] <= glycan_fdr;
            GlycoPsm {
                candidate,
                discriminant: discriminant[i],
                peptide_q: peptide_q[i],
                assignment: assignments[i],
                glycan_q: glycan_q[i],
                passes_fdr,
            }
        })
        .collect();
    summary.passing = psms.iter().filter(|p| p.passes_fdr).count();
    (psms, model, summary)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn q_values_follow_target_decoy_competition() {
        let scores = [5.0, 4.0, 3.0, 2.0, 1.0];
        let decoy = [false, false, true, false, true];
        let q = q_values(&scores, &decoy);
        // (1 + decoys) / targets, then cumulative minimum from the bottom.
        assert_eq!(q, vec![0.5, 0.5, 0.6666667, 0.6666667, 1.0]);
    }

    #[test]
    fn tied_scores_share_a_q_value() {
        let q = q_values(&[2.0, 2.0, 1.0], &[false, true, false]);
        assert_eq!(q[0], q[1]);
    }

    fn explanation(matched: u16, decoy_matched: u16) -> Explanation {
        let mut target = ClassCounts::default();
        target.generated[2] = 10;
        target.matched[2] = matched;
        let mut decoy = target;
        decoy.matched[2] = decoy_matched;
        Explanation {
            composition: 0,
            isotope: 0,
            adducts: 0,
            error_ppm: 1.0,
            target,
            decoy,
            y_intensity: 0.0,
        }
    }

    #[test]
    fn more_matched_ions_win_and_ties_are_a_fair_coin() {
        let model = GlycanModel::new((0, 1), 1, 20.0);
        let score = |e: &Explanation, c: &ClassCounts| model.score(e, c, 0.1, 0.1);
        let e = explanation(6, 1);
        assert!(score(&e, &e.target) > score(&e, &e.decoy));

        let candidate = |explanations| GlycoCandidate {
            feature: sage_core::scoring::Feature {
                spec_id: "scan=1".into(),
                ..Default::default()
            },
            peptide_mass: 1000.0,
            oxonium_ions: 3,
            oxonium_fraction: 0.1,
            anchored: true,
            y1: true,
            hex_ratio: 0.0,
            random_y: 0.1,
            random_oxonium: 0.1,
            explanations,
            site: Default::default(),
            twin_hyperscore: None,
        };
        let a = model.assign(&candidate(vec![explanation(1, 6)]));
        assert!(a.decoy);
        // Ties: deterministic, and roughly half go to the decoy.
        let tie = candidate(vec![explanation(3, 3)]);
        assert_eq!(model.assign(&tie), model.assign(&tie));
        let decoys = (0..1000)
            .filter(|scan| {
                let mut tie = tie.clone();
                tie.feature.spec_id = format!("scan={scan}");
                model.assign(&tie).decoy
            })
            .count();
        assert!((400..600).contains(&decoys), "{decoys}");
        let a = model.assign(&candidate(vec![explanation(1, 1), explanation(7, 2)]));
        assert_eq!((a.explanation, a.decoy), (1, false));
        // The competition statistic adds the lead over the best decoy.
        let (low, high) = (explanation(1, 1), explanation(7, 2));
        let best_decoy = score(&low, &low.decoy).max(score(&high, &high.decoy));
        assert_eq!(a.opposite, best_decoy);
        assert_eq!(a.competition(), 2.0 * a.score - a.opposite);
    }

    #[test]
    fn selection_keeps_the_best_supported_candidate_per_spectrum() {
        let candidate = |scan: &str| GlycoCandidate {
            feature: sage_core::scoring::Feature {
                spec_id: scan.into(),
                ..Default::default()
            },
            peptide_mass: 1000.0,
            oxonium_ions: 3,
            oxonium_fraction: 0.1,
            anchored: true,
            y1: true,
            hex_ratio: 0.0,
            random_y: 0.1,
            random_oxonium: 0.1,
            explanations: vec![explanation(1, 0)],
            site: Default::default(),
            twin_hyperscore: None,
        };
        let assignment = |score: f64| Assignment {
            explanation: 0,
            decoy: false,
            score,
            runner_up: f64::NEG_INFINITY,
            opposite: f64::NEG_INFINITY,
        };
        let candidates = [
            candidate("a"),
            candidate("a"),
            candidate("a"),
            candidate("b"),
            candidate("c"),
            candidate("c"),
        ];
        let assignments: Vec<_> = [1.0, 5.0, 3.0, 2.0, 4.0, 4.0]
            .into_iter()
            .map(assignment)
            .collect();
        let selected = select(&candidates, &assignments);
        assert_eq!(
            selected,
            vec![
                Selection {
                    index: 1,
                    lead: 2.0
                },
                // A lone candidate leads over zero.
                Selection {
                    index: 3,
                    lead: 2.0
                },
                // Ties keep the better-ranked (earlier) peptide.
                Selection {
                    index: 4,
                    lead: 0.0
                },
            ]
        );
    }
}
