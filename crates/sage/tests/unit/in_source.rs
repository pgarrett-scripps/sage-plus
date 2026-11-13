use super::*;
use crate::database::{IndexedDatabase, PeptideIx};
use crate::enzyme::{Digest, ProteinOccurrence};
use crate::modification::{ModificationDefinition, NeutralLossMode};
use crate::peptide::{AppliedModification, CompactModifications, ModificationKind, Peptide, Site};
use crate::scoring::Feature;
use std::sync::Arc;

/// A peptide at `start` in protein `protein`, as the digest would place it.
fn peptide(sequence: &str, protein: &str, start: u32) -> Peptide {
    let mut peptide = Peptide::try_from(Digest {
        sequence: sequence.into(),
        ..Default::default()
    })
    .unwrap();
    peptide.protein_sites = Arc::from([ProteinOccurrence {
        protein: protein.into(),
        start: Some(start),
        prev_aa: None,
        next_aa: None,
        source: None,
        met_clipped: false,
    }]);
    peptide
}

/// Place one modification of `mass` on residue `index`, naming it `name`.
fn modify(peptide: &mut Peptide, index: u32, mass: f32, name: Option<&str>) {
    peptide.modifications = CompactModifications::from_applied([AppliedModification {
        site: Site::Sequence(index),
        modification: Arc::new(ModificationDefinition {
            mass,
            name: name.map(Arc::from),
            neutral_losses: Arc::from([]),
            site_losses: None,
            neutral_loss_mode: NeutralLossMode::Optional,
            channel_offsets: Arc::default(),
        }),
        kind: ModificationKind::Ordinary,
    }])
    .unwrap();
    peptide.monoisotopic += mass;
}

fn psm(psm_id: usize, peptide_idx: u32, file_id: usize, rt: f32) -> Feature {
    Feature {
        psm_id,
        peptide_idx: PeptideIx(peptide_idx),
        file_id,
        rt,
        label: 1,
        spectrum_q: 0.001,
        peptide_q: 0.001,
        ..Default::default()
    }
}

fn database(peptides: Vec<Peptide>) -> IndexedDatabase {
    IndexedDatabase {
        peptides,
        generate_decoys: false,
        ..Default::default()
    }
}

/// A fixed tolerance keeps the unit tests independent of the derivation.
fn settings(rt_tolerance: f32) -> InSourceSettings {
    InSourceSettings {
        rt_tolerance: Some(rt_tolerance),
    }
}

#[test]
fn flags_a_co_eluting_piece_and_names_the_root() {
    // PEPTIDEK at protein position 10; PEPTIDE is its first seven residues.
    let db = database(vec![
        peptide("PEPTIDEK", "P1", 10),
        peptide("PEPTIDE", "P1", 10),
    ]);
    let features = vec![psm(0, 0, 0, 30.0), psm(1, 1, 0, 30.05)];

    let annotations = annotate(&features, &db, &settings(0.2));

    let flagged = annotations
        .psm(1)
        .expect("the piece is an in-source product");
    assert_eq!(flagged.root, PeptideIx(0));
    assert_eq!(flagged.kind, InSourceKind::Fragment);
    assert_eq!(flagged.kind.label(), "fragment");
    // The intact peptide is not a product of anything.
    assert!(annotations.psm(0).is_none());
    assert_eq!(annotations.flagged_psms(), 1);
}

#[test]
fn a_piece_that_does_not_co_elute_is_not_flagged() {
    let db = database(vec![
        peptide("PEPTIDEK", "P1", 10),
        peptide("PEPTIDE", "P1", 10),
    ]);
    // The same two peptides, 4 minutes apart: a separate experimental peptide.
    let features = vec![psm(0, 0, 0, 30.0), psm(1, 1, 0, 34.0)];

    let annotations = annotate(&features, &db, &settings(0.2));

    assert!(annotations.psm(1).is_none());
    // It is still a structural candidate, which is what the check counts.
    let stats = &annotations.files()[0];
    assert_eq!(stats.candidate_pairs, 1);
    assert_eq!(stats.coeluting_pairs, 0);
}

#[test]
fn the_same_peptide_is_a_product_at_one_retention_time_and_not_at_another() {
    let db = database(vec![
        peptide("PEPTIDEK", "P1", 10),
        peptide("PEPTIDE", "P1", 10),
    ]);
    let features = vec![
        psm(0, 0, 0, 30.0),
        psm(1, 1, 0, 30.1),
        psm(2, 1, 0, 44.0),
        psm(3, 1, 0, 44.1),
    ];

    let annotations = annotate(&features, &db, &settings(0.2));

    assert!(annotations.psm(1).is_some());
    assert!(annotations.psm(2).is_none());
    assert!(annotations.psm(3).is_none());
}

#[test]
fn flags_a_water_loss_of_the_whole_peptide() {
    // The same span, 18.010565 Da lighter: a neutral loss in the source.
    let mut lighter = peptide("PEPTIDE", "P1", 10);
    modify(&mut lighter, 2, -crate::mass::H2O, None);
    let db = database(vec![peptide("PEPTIDE", "P1", 10), lighter]);
    let features = vec![psm(0, 0, 0, 20.0), psm(1, 1, 0, 20.02)];

    let annotations = annotate(&features, &db, &settings(0.2));

    let flagged = annotations.psm(1).expect("the lighter form is a product");
    assert_eq!(flagged.root, PeptideIx(0));
    assert_eq!(flagged.kind, InSourceKind::Loss("water".into()));
    assert_eq!(flagged.kind.label(), "loss:water");
}

#[test]
fn flags_the_loss_of_a_whole_modification_by_its_name() {
    // The root carries phosphorylation; the product is the same span without
    // it, which is the labile modification lost in the source.
    let mut root = peptide("PEPSTIDE", "P1", 0);
    modify(&mut root, 3, 79.96633, Some("Phospho"));
    let db = database(vec![root, peptide("PEPSTIDE", "P1", 0)]);
    let features = vec![psm(0, 0, 0, 12.0), psm(1, 1, 0, 12.03)];

    let annotations = annotate(&features, &db, &settings(0.2));

    let flagged = annotations
        .psm(1)
        .expect("the unmodified form is a product");
    assert_eq!(flagged.kind, InSourceKind::Loss("Phospho".into()));
    assert_eq!(flagged.kind.label(), "loss:Phospho");
}

#[test]
fn flags_a_piece_that_also_lost_water() {
    let mut piece = peptide("PEPTIDE", "P1", 10);
    modify(&mut piece, 1, -crate::mass::H2O, None);
    let db = database(vec![peptide("PEPTIDEKRK", "P1", 10), piece]);
    let features = vec![psm(0, 0, 0, 55.0), psm(1, 1, 0, 55.1)];

    let annotations = annotate(&features, &db, &settings(0.2));

    let flagged = annotations.psm(1).expect("a piece can carry one loss");
    assert_eq!(flagged.kind, InSourceKind::FragmentWithLoss("water".into()));
    assert_eq!(flagged.kind.label(), "fragment+loss:water");
}

#[test]
fn a_heavier_piece_is_not_a_product() {
    // The piece carries a modification its span of the root does not: it
    // cannot have come from the root by losing something.
    let mut piece = peptide("PEPTIDE", "P1", 10);
    modify(&mut piece, 1, 42.010565, Some("Acetyl"));
    let db = database(vec![peptide("PEPTIDEK", "P1", 10), piece]);
    let features = vec![psm(0, 0, 0, 30.0), psm(1, 1, 0, 30.02)];

    let annotations = annotate(&features, &db, &settings(0.2));

    assert!(annotations.psm(1).is_none());
    assert_eq!(annotations.files()[0].candidate_pairs, 0);
}

#[test]
fn a_mass_difference_that_is_no_listed_loss_is_not_flagged() {
    let mut piece = peptide("PEPTIDE", "P1", 10);
    modify(&mut piece, 1, -13.37, None);
    let db = database(vec![peptide("PEPTIDEK", "P1", 10), piece]);
    let features = vec![psm(0, 0, 0, 30.0), psm(1, 1, 0, 30.02)];

    let annotations = annotate(&features, &db, &settings(0.2));

    assert!(annotations.psm(1).is_none());
}

#[test]
fn the_root_is_the_longest_co_eluting_peptide_that_contains_the_product() {
    // A piece of a piece points at the intact peptide, not at the middle one.
    let db = database(vec![
        peptide("PEPTIDEKRSTK", "P1", 4),
        peptide("PEPTIDEK", "P1", 4),
        peptide("PEPTIDE", "P1", 4),
    ]);
    let features = vec![psm(0, 0, 0, 40.0), psm(1, 1, 0, 40.02), psm(2, 2, 0, 40.04)];

    let annotations = annotate(&features, &db, &settings(0.2));

    assert_eq!(annotations.psm(2).unwrap().root, PeptideIx(0));
    assert_eq!(annotations.psm(1).unwrap().root, PeptideIx(0));
}

#[test]
fn only_the_longest_root_that_co_elutes_is_named() {
    // The intact peptide elutes elsewhere; the middle peptide co-elutes, so it
    // is the root of the shortest piece.
    let db = database(vec![
        peptide("PEPTIDEKRSTK", "P1", 4),
        peptide("PEPTIDEK", "P1", 4),
        peptide("PEPTIDE", "P1", 4),
    ]);
    let features = vec![psm(0, 0, 0, 10.0), psm(1, 1, 0, 40.02), psm(2, 2, 0, 40.04)];

    let annotations = annotate(&features, &db, &settings(0.2));

    assert_eq!(annotations.psm(2).unwrap().root, PeptideIx(1));
    assert!(annotations.psm(1).is_none());
}

#[test]
fn a_piece_in_another_protein_is_not_a_product_of_this_one() {
    let db = database(vec![
        peptide("PEPTIDEK", "P1", 10),
        peptide("PEPTIDE", "P2", 300),
    ]);
    let features = vec![psm(0, 0, 0, 30.0), psm(1, 1, 0, 30.02)];

    let annotations = annotate(&features, &db, &settings(0.2));

    assert!(annotations.psm(1).is_none());
}

#[test]
fn a_product_and_its_root_must_be_in_the_same_file() {
    let db = database(vec![
        peptide("PEPTIDEK", "P1", 10),
        peptide("PEPTIDE", "P1", 10),
    ]);
    let features = vec![psm(0, 0, 0, 30.0), psm(1, 1, 1, 30.02)];

    let annotations = annotate(&features, &db, &settings(0.2));

    assert!(annotations.psm(1).is_none());
    assert_eq!(annotations.files().len(), 2);
}

#[test]
fn unconfident_and_decoy_psms_take_no_part() {
    let mut decoy_root = peptide("PEPTIDEK", "P1", 10);
    decoy_root.decoy = true;
    let db = database(vec![
        peptide("PEPTIDEK", "P1", 10),
        peptide("PEPTIDE", "P1", 10),
        decoy_root,
    ]);

    // A root above the q-value cutoff cannot be a parent.
    let mut unconfident = psm(0, 0, 0, 30.0);
    unconfident.spectrum_q = 0.5;
    let annotations = annotate(&[unconfident, psm(1, 1, 0, 30.02)], &db, &settings(0.2));
    assert!(annotations.psm(1).is_none());

    // Neither can a decoy, nor is a decoy product flagged.
    let mut decoy_psm = psm(2, 2, 0, 30.0);
    decoy_psm.label = -1;
    let annotations = annotate(&[decoy_psm, psm(1, 1, 0, 30.02)], &db, &settings(0.2));
    assert!(annotations.psm(1).is_none());
}

#[test]
fn a_quantified_row_is_judged_on_its_own_apex_retention_time() {
    let db = database(vec![
        peptide("PEPTIDEK", "P1", 10),
        peptide("PEPTIDE", "P1", 10),
    ]);
    // The root's PSMs sit at 30.0; the product's PSM at 34.0 is not flagged,
    // but an LFQ peak of the product apexing at 30.05 is.
    let features = vec![psm(0, 0, 0, 30.0), psm(1, 1, 0, 34.0)];

    let annotations = annotate(&features, &db, &settings(0.2));

    assert!(annotations.psm(1).is_none());
    let quantified = annotations
        .precursor(PeptideIx(1), 0, 30.05)
        .expect("the quantified peak co-elutes with the root");
    assert_eq!(quantified.root, PeptideIx(0));
    assert_eq!(quantified.kind, InSourceKind::Fragment);
    assert!(annotations.precursor(PeptideIx(1), 0, 34.0).is_none());
    // A file the annotation never saw has no tolerance, so no verdict.
    assert!(annotations.precursor(PeptideIx(1), 7, 30.05).is_none());
}

#[test]
fn the_check_reports_the_histogram_and_what_chance_would_give() {
    let db = database(vec![
        peptide("PEPTIDEK", "P1", 10),
        peptide("PEPTIDE", "P1", 10),
    ]);
    // One co-eluting product, and twelve unrelated PSMs of the same peptide
    // spread over the run, which is what the background is estimated from.
    let mut features = vec![psm(0, 0, 0, 30.0), psm(1, 1, 0, 30.02)];
    for step in 0..12 {
        features.push(psm(2 + step, 1, 0, 31.0 + step as f32 * 0.2));
    }

    let annotations = annotate(&features, &db, &settings(0.1));
    let stats = &annotations.files()[0];

    assert_eq!(stats.file_id, 0);
    assert_eq!(stats.rt_tolerance, 0.1);
    assert_eq!(stats.rt_tolerance_source, "configured");
    assert_eq!(stats.confident_psms, 14);
    assert_eq!(stats.candidate_pairs, 13);
    assert_eq!(stats.coeluting_pairs, 1);
    assert_eq!(stats.flagged_psms, 1);
    // Bins are one tolerance wide and sorted, and the co-eluting pair sits in
    // the bin that starts at zero.
    assert!(stats
        .rt_difference_histogram
        .windows(2)
        .all(|pair| pair[0].bin < pair[1].bin));
    let zero = stats
        .rt_difference_histogram
        .iter()
        .find(|entry| entry.bin == 0)
        .unwrap();
    assert_eq!(zero.pairs, 1);
    assert_eq!(zero.start, 0.0);
    assert!((zero.end - 0.1).abs() < 1e-6);
    // The estimate is the background bins' mean count over the two bins of
    // the co-elution window, and the unrelated PSMs put it far below what
    // actually co-eluted.
    let background: usize = stats
        .rt_difference_histogram
        .iter()
        .filter(|entry| BACKGROUND_BIN_RANGE.contains(&entry.bin.abs()))
        .map(|entry| entry.pairs)
        .sum();
    assert!(background >= 7, "background bins hold {background} pairs");
    let bins = BACKGROUND_BIN_RANGE.len() * 2;
    assert!((stats.expected_by_chance - 2.0 * background as f32 / bins as f32).abs() < 1e-6);
    assert!(stats.expected_by_chance < stats.coeluting_pairs as f32);
}

#[test]
fn the_tolerance_comes_from_the_spread_of_repeated_peptides() {
    // Twenty observations, 0.00 to 0.19: the 95th percentile is the value at
    // rank round(0.95 * 19) = 18.
    let mut spreads = (0..20).map(|step| step as f32 / 100.0).collect::<Vec<_>>();
    let (tolerance, source) = derive_rt_tolerance(&mut spreads, 60.0);
    assert_eq!(source, "derived");
    assert!((tolerance - 0.18).abs() < 1e-6);

    // Too few observations: a fraction of the run instead.
    let mut thin = vec![0.5; MIN_SPREAD_OBSERVATIONS - 1];
    let (tolerance, source) = derive_rt_tolerance(&mut thin, 60.0);
    assert_eq!(source, "fallback");
    assert!((tolerance - 60.0 * FALLBACK_RT_FRACTION).abs() < 1e-6);

    // Repeats all at one retention time, and no usable span.
    let mut flat = vec![0.0; MIN_SPREAD_OBSERVATIONS];
    let (tolerance, source) = derive_rt_tolerance(&mut flat, 0.0);
    assert_eq!(source, "fallback");
    assert_eq!(tolerance, FALLBACK_RT_TOLERANCE);
}

#[test]
fn a_derived_tolerance_is_recorded_per_file() {
    let db = database(vec![
        peptide("PEPTIDEK", "P1", 10),
        peptide("PEPTIDE", "P1", 10),
    ]);
    // Enough repeats of the root to derive a tolerance from the file itself.
    let mut features = vec![psm(0, 1, 0, 30.0)];
    for step in 0..MIN_SPREAD_OBSERVATIONS + 1 {
        features.push(psm(1 + step, 0, 0, 30.0 + step as f32 * 0.01));
    }

    let annotations = annotate(&features, &db, &InSourceSettings::default());
    let stats = &annotations.files()[0];

    assert_eq!(stats.rt_tolerance_source, "derived");
    assert!(stats.rt_tolerance > 0.0);
}
