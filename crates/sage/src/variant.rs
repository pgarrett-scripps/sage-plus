//! Variant peptides: peptides of `database.variant_fasta` proteins only.
//!
//! A variant peptide is one whose every protein is a variant protein (see
//! [`ProteinClasses::is_variant`]); a sequence that a reference protein also
//! explains is a reference peptide. Variant peptides compete in the search
//! like any other peptide, so the global q-values keep their meaning, and get
//! a second, class-specific q-value as well: target-decoy FDR estimated within
//! the variant class alone (a separate FDR), at the PSM and the peptide level.
//! A variant class is much smaller than the reference proteome and its false
//! discoveries need not be in proportion to the whole run's, which is why the
//! global q-value alone does not control the error rate among variant
//! identifications.
//!
//! The estimator is the run's own: `(decoys + 1) / targets` over complete
//! tied-score groups with a running minimum ([`grouped_q_values`]) for PSMs,
//! and the picked-peptide competition ([`picked_peptide_q_values`]) for
//! peptides, each over variant-class PSMs only. A class with few decoys is
//! conservative, not optimistic: with the +1, fewer than about 100 target PSMs
//! above a threshold can never reach 1%.

use crate::database::{IndexedDatabase, PeptideIx};
use crate::fasta::ProteinClasses;
use crate::fdr::picked_peptide_q_values;
use crate::ml::qvalue::grouped_q_values;
use crate::scoring::Feature;
use fnv::{FnvHashMap, FnvHashSet};
use itertools::{Either, Itertools};
use serde::{Deserialize, Serialize};

/// The q-value threshold the counts are reported at.
pub const Q_VALUE: f32 = 0.01;

/// Class-specific q-values of the variant PSMs of a run, keyed by `psm_id`.
/// PSMs that are not variant-class have no entry.
#[derive(Clone, Debug, Default)]
pub struct VariantQValues {
    spectrum: FnvHashMap<usize, f32>,
    peptide: FnvHashMap<usize, f32>,
}

impl VariantQValues {
    /// The variant spectrum-level q-value of PSM `psm_id`, `None` when the
    /// PSM's peptide is not a variant peptide.
    pub fn spectrum_q(&self, psm_id: usize) -> Option<f32> {
        self.spectrum.get(&psm_id).copied()
    }

    /// The variant peptide-level q-value of PSM `psm_id`, `None` when the
    /// PSM's peptide is not a variant peptide.
    pub fn peptide_q(&self, psm_id: usize) -> Option<f32> {
        self.peptide.get(&psm_id).copied()
    }
}

/// Targets and decoys of one level of the variant class, and how many targets
/// pass at [`Q_VALUE`] by the class q-value and by the global one. These are
/// what explain an empty variant result: a class needs about 100 targets above
/// a threshold before anything can pass at 1%.
#[derive(
    Clone, Copy, Debug, Default, PartialEq, Eq, Serialize, Deserialize, schemars::JsonSchema,
)]
pub struct VariantClassCounts {
    /// Target PSMs (or peptides) of the class that took part.
    pub targets: usize,
    /// Decoy PSMs (or peptides) of the class that took part.
    pub decoys: usize,
    /// Targets at or below 1% by the class-specific q-value.
    pub passing_variant_q: usize,
    /// Targets at or below 1% by the run's global q-value.
    pub passing_global_q: usize,
}

/// Is this PSM's peptide a variant peptide? Decoys keep the accessions of the
/// targets they were reversed from, so a variant peptide's decoy is in the
/// class too and the class has its own decoys.
pub fn is_variant(db: &IndexedDatabase, classes: &ProteinClasses, feature: &Feature) -> bool {
    classes.is_variant(db[feature.peptide_idx].proteins.iter())
}

/// Variant spectrum- and peptide-level q-values of `features`, which must
/// already carry their discriminant score and global q-values, with the
/// class counts at both levels.
pub fn variant_q_values(
    db: &IndexedDatabase,
    features: &[Feature],
    classes: &ProteinClasses,
) -> (VariantQValues, VariantClassCounts, VariantClassCounts) {
    let class = features
        .iter()
        .filter(|feature| is_variant(db, classes, feature))
        .collect::<Vec<_>>();
    let spectrum_q = grouped_q_values(
        &class,
        |feature| f64::from(feature.discriminant_score),
        |feature| feature.label == -1,
        |_| (),
    );
    let (peptide_q, _) =
        picked_peptide_q_values(db, features, |feature| is_variant(db, classes, feature));

    let mut q_values = VariantQValues::default();
    let mut psms = VariantClassCounts::default();
    for (feature, &q) in class.iter().zip(&spectrum_q) {
        q_values.spectrum.insert(feature.psm_id, q);
        match feature.label == -1 {
            true => psms.decoys += 1,
            false => {
                psms.targets += 1;
                psms.passing_variant_q += usize::from(q <= Q_VALUE);
                psms.passing_global_q += usize::from(feature.spectrum_q <= Q_VALUE);
            }
        }
    }

    let mut peptides = VariantClassCounts::default();
    let mut seen: FnvHashSet<PeptideIx> = FnvHashSet::default();
    let mut passing: FnvHashSet<(PeptideIx, bool)> = FnvHashSet::default();
    for (feature, q) in features.iter().zip(peptide_q) {
        let Some(q) = q else { continue };
        q_values.peptide.insert(feature.psm_id, q);
        let decoy = db[feature.peptide_idx].decoy;
        if seen.insert(feature.peptide_idx) {
            match decoy {
                true => peptides.decoys += 1,
                false => peptides.targets += 1,
            }
        }
        if !decoy && q <= Q_VALUE && passing.insert((feature.peptide_idx, true)) {
            peptides.passing_variant_q += 1;
        }
        if !decoy && feature.peptide_q <= Q_VALUE && passing.insert((feature.peptide_idx, false)) {
            peptides.passing_global_q += 1;
        }
    }
    (q_values, psms, peptides)
}

/// Run protein inference (`infer`: picked protein FDR, grouping, picked
/// protein-group FDR) on every PSM except the variant ones, then give each
/// variant PSM a protein q-value and protein-group q-value of 1 and its raw
/// protein list as its group, as an ungrouped PSM gets. Variant peptides are
/// kept out because their error rate is the class's, not the proteome's: a
/// variant peptide counted as protein evidence would let variant false
/// discoveries move the protein results. The order of `features` is kept.
pub fn infer_proteins_without_variants<R>(
    db: &IndexedDatabase,
    classes: &ProteinClasses,
    features: &mut Vec<Feature>,
    infer: impl FnOnce(&mut [Feature]) -> R,
) -> R {
    let held = features
        .iter()
        .map(|feature| is_variant(db, classes, feature))
        .collect::<Vec<_>>();
    if !held.contains(&true) {
        return infer(features);
    }
    let (mut variant, mut reference): (Vec<Feature>, Vec<Feature>) = std::mem::take(features)
        .into_iter()
        .zip(&held)
        .partition_map(|(feature, &held)| match held {
            true => Either::Left(feature),
            false => Either::Right(feature),
        });
    let result = infer(&mut reference);
    for feature in &mut variant {
        let peptide = &db[feature.peptide_idx];
        feature.protein_q = 1.0;
        feature.protein_group_q = 1.0;
        feature.protein_groups = Some(peptide.proteins(&db.decoy_tag, db.generate_decoys));
        feature.num_protein_groups = peptide.proteins.len() as u32;
    }
    let mut variant = variant.into_iter();
    let mut reference = reference.into_iter();
    features.extend(held.into_iter().map(|held| {
        match held {
            true => variant.next(),
            false => reference.next(),
        }
        .expect("every feature was split into exactly one side")
    }));
    result
}

#[cfg(test)]
#[path = "../tests/unit/variant.rs"]
mod test;
