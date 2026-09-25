//! Crosslink rescoring and false discovery rates.
//!
//! A linear discriminant combines the CSM features, trained with TT matches as
//! targets and TD/DD matches as decoys. FDR is `(TD - DD) / TT`, computed
//! separately for intra- and inter-protein links (the two have very different
//! prior rates of false matches), at the CSM and at the residue-pair level.

use crate::search::{protein_position, Class, Csm};
use sage_core::database::IndexedDatabase;
use sage_core::ml::linear_discriminant::LinearDiscriminantAnalysis;
use std::collections::HashMap;
use std::sync::Arc;

const FEATURES: usize = 12;
const REGULARIZATION: f64 = 1e-3;

/// Target counts passing 1% FDR.
#[derive(Copy, Clone, Debug, Default, PartialEq, Eq)]
pub struct FdrSummary {
    pub intra_csms: usize,
    pub inter_csms: usize,
    pub intra_residue_pairs: usize,
    pub inter_residue_pairs: usize,
    /// Whether the discriminant was fitted (else the combined hyperscore is
    /// used).
    pub discriminant_fitted: bool,
}

fn features(csm: &Csm) -> [f64; FEATURES] {
    let (a, b) = (&csm.alpha, &csm.beta);
    [
        csm.hyperscore.ln_1p(),
        a.hyperscore.min(b.hyperscore).max(0.0).ln_1p(),
        a.hyperscore.max(b.hyperscore).max(0.0).ln_1p(),
        csm.delta_next.max(0.0).ln_1p(),
        (a.doublet as u8 + b.doublet as u8) as f64,
        (csm.precursor_ppm.abs() as f64).ln_1p(),
        csm.isotope_error.unsigned_abs() as f64,
        (a.matched_peaks.min(b.matched_peaks) as f64).ln_1p(),
        (a.matched_peaks.max(b.matched_peaks) as f64).ln_1p(),
        (csm.matched_intensity_pct.max(0.0) as f64).ln_1p(),
        csm.charge as f64,
        (a.length.min(b.length) as f64).ln(),
    ]
}

/// Score CSMs with a discriminant and assign CSM and residue-pair q-values.
pub fn assign_q_values(csms: &mut [Csm], db: &IndexedDatabase) -> FdrSummary {
    let decoy: Vec<bool> = csms.iter().map(|c| c.class != Class::TT).collect();
    let lda = LinearDiscriminantAnalysis::train_regularized(csms, &decoy, features, REGULARIZATION)
        .filter(|lda| {
            csms.iter()
                .take(1)
                .all(|c| lda.score(&features(c)).is_finite())
        });
    let fitted = lda.is_some();
    for csm in csms.iter_mut() {
        csm.discriminant_score = match &lda {
            Some(lda) => lda.score(&features(csm)),
            None => csm.hyperscore,
        };
    }
    if !fitted {
        log::warn!("crosslink discriminant could not be fitted; using the combined hyperscore");
    }

    let mut summary = FdrSummary {
        discriminant_fitted: fitted,
        ..Default::default()
    };

    // CSM level.
    let scored: Vec<(f64, Class, bool)> = csms
        .iter()
        .map(|c| (c.discriminant_score, c.class, c.intra))
        .collect();
    for (csm, q) in csms.iter_mut().zip(grouped_q_values(&scored)) {
        csm.csm_q = q;
    }
    for csm in csms.iter() {
        if csm.class == Class::TT && csm.csm_q <= 0.01 {
            if csm.intra {
                summary.intra_csms += 1;
            } else {
                summary.inter_csms += 1;
            }
        }
    }

    // Residue-pair level: aggregate the CSMs of each linked residue pair and
    // score the pair with its own discriminant (see `PairEvidence`).
    let keys: Vec<Option<ResiduePair>> = csms.iter().map(|c| residue_pair(db, c)).collect();
    let mut grouped: HashMap<&ResiduePair, PairEvidence> = HashMap::new();
    for (csm, key) in csms.iter().zip(&keys) {
        if let Some(key) = key {
            grouped
                .entry(key)
                .and_modify(|e| e.add(csm))
                .or_insert_with(|| PairEvidence::new(csm));
        }
    }
    let mut pairs: Vec<(&ResiduePair, PairEvidence)> = grouped.into_iter().collect();
    pairs.sort_unstable_by(|a, b| a.0.cmp(b.0));
    let evidence: Vec<&PairEvidence> = pairs.iter().map(|(_, e)| e).collect();
    let pair_decoy: Vec<bool> = evidence.iter().map(|e| e.class != Class::TT).collect();
    let pair_lda = LinearDiscriminantAnalysis::train_regularized(
        &evidence,
        &pair_decoy,
        |e: &&PairEvidence| e.features(),
        REGULARIZATION,
    )
    .filter(|lda| {
        evidence
            .iter()
            .all(|e| lda.score(&e.features()).is_finite())
    });
    if pair_lda.is_none() && !pairs.is_empty() {
        log::warn!("residue-pair discriminant could not be fitted; using the best CSM score");
    }
    let scored: Vec<(f64, Class, bool)> = evidence
        .iter()
        .map(|e| {
            let score = match &pair_lda {
                Some(lda) => lda.score(&e.features()),
                None => e.best,
            };
            (score, e.class, e.intra)
        })
        .collect();
    let pair_q: HashMap<ResiduePair, f32> = pairs
        .iter()
        .zip(grouped_q_values(&scored))
        .map(|((key, _), q)| ((*key).clone(), q))
        .collect();
    for ((key, _), (_, class, intra)) in pairs.iter().zip(&scored) {
        if *class == Class::TT && pair_q[*key] <= 0.01 {
            if *intra {
                summary.intra_residue_pairs += 1;
            } else {
                summary.inter_residue_pairs += 1;
            }
        }
    }
    for (csm, key) in csms.iter_mut().zip(&keys) {
        csm.residue_pair_q = key.as_ref().map_or(1.0, |key| pair_q[key]);
    }
    summary
}

/// CSMs at or below this CSM q-value count as support for their residue pair.
const SUPPORT_Q: f32 = 0.05;

/// Evidence for one residue pair, pooled over its CSMs: the best CSM score
/// and the number of supporting CSMs (CSM q <= `SUPPORT_Q`). A pair seen in a
/// single spectrum is still scored, by its best CSM. Repeated support raises
/// the score, because random wrong pairs, like decoy pairs, rarely collect
/// several good CSMs. Counting every CSM, or adding distinct charge states or
/// peptide forms, did not help consistently on the two ground-truth sets.
struct PairEvidence {
    best: f64,
    supporting: usize,
    class: Class,
    intra: bool,
}

impl PairEvidence {
    fn new(csm: &Csm) -> Self {
        let mut evidence = PairEvidence {
            best: csm.discriminant_score,
            supporting: 0,
            class: csm.class,
            intra: csm.intra,
        };
        evidence.add(csm);
        evidence
    }

    fn add(&mut self, csm: &Csm) {
        if csm.discriminant_score > self.best {
            self.best = csm.discriminant_score;
            self.class = csm.class;
            self.intra = csm.intra;
        }
        if csm.csm_q <= SUPPORT_Q {
            self.supporting += 1;
        }
    }

    fn features(&self) -> [f64; 2] {
        [self.best, (self.supporting as f64).ln_1p()]
    }
}

/// Unordered pair of (protein, 1-based position), decoy tags kept.
type ResiduePair = [(Arc<str>, u32); 2];

fn residue_pair(db: &IndexedDatabase, csm: &Csm) -> Option<ResiduePair> {
    let a = protein_position(&db[csm.alpha.peptide], csm.alpha.site)?;
    let b = protein_position(&db[csm.beta.peptide], csm.beta.site)?;
    Some(if a <= b { [a, b] } else { [b, a] })
}

/// q-values with FDR `(TD - DD) / TT`, computed separately for each
/// `intra` group and made monotone. Input order is preserved.
pub fn grouped_q_values(items: &[(f64, Class, bool)]) -> Vec<f32> {
    let mut q = vec![1.0f32; items.len()];
    for group in [true, false] {
        let indices: Vec<usize> = (0..items.len()).filter(|&i| items[i].2 == group).collect();
        let scored: Vec<(f64, Class)> = indices.iter().map(|&i| (items[i].0, items[i].1)).collect();
        for (index, value) in indices.into_iter().zip(q_values(&scored)) {
            q[index] = value;
        }
    }
    q
}

/// q-values with FDR `(TD - DD) / TT` for one group. Input order is
/// preserved; tied scores share a q-value.
pub fn q_values(items: &[(f64, Class)]) -> Vec<f32> {
    let mut order: Vec<usize> = (0..items.len()).collect();
    order.sort_by(|&a, &b| items[b].0.total_cmp(&items[a].0));
    let mut fdr = vec![1.0f32; items.len()];
    let (mut tt, mut td, mut dd) = (0usize, 0usize, 0usize);
    let mut start = 0;
    while start < order.len() {
        let mut end = start;
        while end < order.len() && items[order[end]].0 == items[order[start]].0 {
            match items[order[end]].1 {
                Class::TT => tt += 1,
                Class::TD => td += 1,
                Class::DD => dd += 1,
            }
            end += 1;
        }
        let value = (td.saturating_sub(dd)) as f32 / tt.max(1) as f32;
        for &index in &order[start..end] {
            fdr[index] = value.min(1.0);
        }
        start = end;
    }
    // Monotone: q at a score is the lowest FDR at that score or below.
    let mut running = 1.0f32;
    for &index in order.iter().rev() {
        running = running.min(fdr[index]);
        fdr[index] = running;
    }
    fdr
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn fdr_subtracts_dd_from_td() {
        use Class::*;
        let items = [
            (10.0, TT),
            (9.0, TT),
            (8.0, TD),
            (7.0, TT),
            (6.0, TT),
            (5.0, DD),
            (4.0, TD),
            (3.0, TT),
        ];
        let q = q_values(&items);
        assert_eq!(q[0], 0.0);
        assert_eq!(q[1], 0.0);
        // At score 6: 4 TT, 1 TD, 0 DD -> 0.25; at score 5 still 1 TD - 1 DD.
        assert_eq!(q[4], 0.0);
        assert_eq!(q[2], 0.0);
        // At score 3: 5 TT, 2 TD, 1 DD -> 0.2.
        assert!((q[7] - 0.2).abs() < 1e-6);
    }

    #[test]
    fn q_values_are_monotone_and_tie_aware() {
        use Class::*;
        let items = [(5.0, TT), (5.0, TD), (4.0, TT), (3.0, TT)];
        let q = q_values(&items);
        assert_eq!(q[0], q[1]);
        assert!((q[3] - 1.0 / 3.0).abs() < 1e-6);
        assert!((q[0] - 1.0 / 3.0).abs() < 1e-6);
    }

    fn csm(score: f64, q: f32, class: Class) -> Csm {
        let chain = crate::search::ChainMatch {
            peptide: sage_core::database::PeptideIx(0),
            site: sage_core::peptide::Site::Sequence(0),
            hyperscore: 0.0,
            matched_peaks: 0,
            doublet: false,
            length: 8,
        };
        Csm {
            file_id: 0,
            spectrum_id: String::new(),
            rt: 0.0,
            charge: 3,
            expmass: 0.0,
            calcmass: 0.0,
            isotope_error: 0,
            precursor_ppm: 0.0,
            alpha: chain.clone(),
            beta: chain,
            hyperscore: 0.0,
            delta_next: 0.0,
            matched_peaks: 0,
            matched_intensity_pct: 0.0,
            doublets: 0,
            class,
            intra: true,
            discriminant_score: score,
            csm_q: q,
            residue_pair_q: 1.0,
        }
    }

    #[test]
    fn pair_evidence_counts_only_supporting_csms() {
        let mut pair = PairEvidence::new(&csm(2.0, 0.2, Class::TT));
        // A single CSM that does not pass support is still scored by its best score.
        assert_eq!(pair.features(), [2.0, 0.0]);
        pair.add(&csm(5.0, 0.01, Class::TT));
        pair.add(&csm(1.0, 0.04, Class::TT));
        assert_eq!(pair.best, 5.0);
        assert_eq!(pair.supporting, 2);
        assert!((pair.features()[1] - 3f64.ln()).abs() < 1e-12);
    }

    #[test]
    fn best_csm_sets_pair_class() {
        let mut pair = PairEvidence::new(&csm(1.0, 0.5, Class::TD));
        pair.add(&csm(3.0, 0.5, Class::TT));
        assert_eq!(pair.class, Class::TT);
    }

    #[test]
    fn groups_are_independent() {
        use Class::*;
        let items = [(10.0, TT, true), (9.0, TD, false), (8.0, TT, false)];
        let q = grouped_q_values(&items);
        assert_eq!(q[0], 0.0);
        assert_eq!(q[2], 1.0);
    }
}
