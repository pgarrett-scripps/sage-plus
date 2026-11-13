use super::*;
use crate::database::{IndexedDatabase, PeptideIx};
use crate::enzyme::{Digest, ProteinOccurrence};
use crate::mass::{H2O, NH3};
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

fn definition(mass: f32, name: Option<&str>, losses: &[f32]) -> Arc<ModificationDefinition> {
    Arc::new(ModificationDefinition {
        mass,
        name: name.map(Arc::from),
        neutral_losses: Arc::from(losses),
        site_losses: None,
        neutral_loss_mode: NeutralLossMode::Optional,
        channel_offsets: Arc::default(),
    })
}

/// Place modifications on a peptide, keeping its mass in step.
fn modify(peptide: &mut Peptide, sites: &[(u32, Arc<ModificationDefinition>, ModificationKind)]) {
    peptide.modifications =
        CompactModifications::from_applied(sites.iter().map(|(index, modification, kind)| {
            AppliedModification {
                site: Site::Sequence(*index),
                modification: modification.clone(),
                kind: *kind,
            }
        }))
        .unwrap();
    peptide.monoisotopic += sites
        .iter()
        .map(|(_, modification, _)| modification.mass)
        .sum::<f32>();
}

/// A peptide carrying one ordinary modification of `mass` on residue `index`.
fn modified(sequence: &str, protein: &str, start: u32, index: u32, mass: f32) -> Peptide {
    let mut peptide = peptide(sequence, protein, start);
    modify(
        &mut peptide,
        &[(
            index,
            definition(mass, None, &[]),
            ModificationKind::Ordinary,
        )],
    );
    peptide
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
fn flags_a_co_eluting_suffix_and_names_the_parent() {
    // KPEPTIDE at protein position 10; PEPTIDE is its C-terminal seven
    // residues, which have the ordinary peptide mass of that sequence.
    let db = database(vec![
        peptide("KPEPTIDE", "P1", 10),
        peptide("PEPTIDE", "P1", 11),
    ]);
    let features = vec![psm(0, 0, 0, 30.0), psm(1, 1, 0, 30.05)];

    let annotations = annotate(&features, &db, &settings(0.2));

    let flagged = annotations
        .psm(1)
        .expect("the suffix is an in-source product");
    assert_eq!(flagged.parent, PeptideIx(0));
    assert_eq!(flagged.kind.placement, Placement::Suffix);
    assert_eq!(flagged.kind.loss, None);
    assert_eq!(flagged.kind.label(), "fragment:suffix");
    // The intact peptide is not a product of anything.
    assert!(annotations.psm(0).is_none());
    assert_eq!(annotations.flagged_psms(), 1);
}

#[test]
fn a_prefix_at_the_intact_peptide_mass_is_not_a_product() {
    // Cleavage in the source leaves the N-terminal piece 18.010565 Da lighter
    // than the peptide of the same sequence, so a peptide at that sequence's
    // ordinary mass did not come from this parent.
    let db = database(vec![
        peptide("PEPTIDEK", "P1", 10),
        peptide("PEPTIDE", "P1", 10),
    ]);
    let features = vec![psm(0, 0, 0, 30.0), psm(1, 1, 0, 30.02)];

    let annotations = annotate(&features, &db, &settings(0.2));

    assert!(annotations.psm(1).is_none());
    assert_eq!(annotations.files()[0].totals.candidate_pairs, 0);
}

#[test]
fn flags_a_water_light_prefix_and_internal_product() {
    // Both carry the water deficit of their new C-terminus, as a b-type
    // species does.
    let prefix = modified("PEPTIDE", "P1", 10, 1, -H2O);
    let internal = modified("EPTID", "P1", 11, 1, -H2O);
    let db = database(vec![peptide("PEPTIDEKR", "P1", 10), prefix, internal]);
    let features = vec![psm(0, 0, 0, 44.0), psm(1, 1, 0, 44.02), psm(2, 2, 0, 44.04)];

    let annotations = annotate(&features, &db, &settings(0.2));

    let flagged = annotations.psm(1).expect("the prefix is water-light");
    assert_eq!(flagged.kind.placement, Placement::Prefix);
    assert_eq!(flagged.kind.loss, None);
    assert_eq!(flagged.kind.label(), "fragment:prefix");
    let flagged = annotations
        .psm(2)
        .expect("the internal piece is water-light");
    assert_eq!(flagged.kind.placement, Placement::Internal);
    assert_eq!(flagged.kind.label(), "fragment:internal");
}

#[test]
fn flags_a_water_loss_of_the_whole_peptide_by_its_mass() {
    let lighter = modified("PEPTIDE", "P1", 10, 2, -H2O);
    let db = database(vec![peptide("PEPTIDE", "P1", 10), lighter]);
    let features = vec![psm(0, 0, 0, 20.0), psm(1, 1, 0, 20.02)];

    let annotations = annotate(&features, &db, &settings(0.2));

    let flagged = annotations.psm(1).expect("the lighter form is a product");
    assert_eq!(flagged.parent, PeptideIx(0));
    assert_eq!(flagged.kind.placement, Placement::SameSpan);
    assert_eq!(flagged.kind.loss, Some(H2O));
    assert_eq!(flagged.kind.label(), "loss:18.0106");
}

#[test]
fn flags_a_suffix_that_also_lost_ammonia() {
    let mut product = peptide("PEPTIDE", "P1", 13);
    modify(
        &mut product,
        &[(1, definition(-NH3, None, &[]), ModificationKind::Ordinary)],
    );
    let db = database(vec![peptide("KRTPEPTIDE", "P1", 10), product]);
    let features = vec![psm(0, 0, 0, 55.0), psm(1, 1, 0, 55.1)];

    let annotations = annotate(&features, &db, &settings(0.2));

    let flagged = annotations.psm(1).expect("a product may carry one loss");
    assert_eq!(flagged.kind.placement, Placement::Suffix);
    assert_eq!(flagged.kind.loss, Some(NH3));
    assert_eq!(flagged.kind.label(), "fragment:suffix+loss:17.0265");
}

#[test]
fn a_whole_modification_is_a_loss_only_where_the_modification_declares_it() {
    // The parent carries phosphorylation that declares the 79.96633 loss, and
    // the product is the same span without it: the labile modification lost in
    // the source.
    let mut declared_parent = peptide("PEPSTIDE", "P1", 0);
    modify(
        &mut declared_parent,
        &[(
            3,
            definition(79.96633, Some("Phospho"), &[79.96633]),
            ModificationKind::Ordinary,
        )],
    );
    let db = database(vec![declared_parent, peptide("PEPSTIDE", "P1", 0)]);
    let annotations = annotate(
        &[psm(0, 0, 0, 12.0), psm(1, 1, 0, 12.03)],
        &db,
        &settings(0.2),
    );
    let flagged = annotations
        .psm(1)
        .expect("the declared loss makes this a product");
    assert_eq!(flagged.kind.placement, Placement::SameSpan);
    assert_eq!(flagged.kind.label(), "loss:79.9663");

    // The same shape with a modification that declares no loss is not a
    // product: an unmodified peptide co-eluting with its modified form is not
    // evidence of in-source fragmentation.
    let mut silent_parent = peptide("PEPSTIDE", "P1", 0);
    modify(
        &mut silent_parent,
        &[(
            3,
            definition(79.96633, Some("Phospho"), &[]),
            ModificationKind::Ordinary,
        )],
    );
    let db = database(vec![silent_parent, peptide("PEPSTIDE", "P1", 0)]);
    let annotations = annotate(
        &[psm(0, 0, 0, 12.0), psm(1, 1, 0, 12.03)],
        &db,
        &settings(0.2),
    );
    assert!(annotations.psm(1).is_none());
}

#[test]
fn a_labelling_modification_never_makes_a_product() {
    // A labelling pair: the heavy channel declares its own mass as a loss, and
    // it still must not turn the light channel into its product.
    let mut heavy = peptide("PEPTIDEK", "P1", 0);
    modify(
        &mut heavy,
        &[(
            7,
            definition(8.014199, Some("Lys8"), &[8.014199]),
            ModificationKind::Label,
        )],
    );
    let db = database(vec![heavy, peptide("PEPTIDEK", "P1", 0)]);

    let annotations = annotate(
        &[psm(0, 0, 0, 25.0), psm(1, 1, 0, 25.01)],
        &db,
        &settings(0.2),
    );

    assert!(annotations.psm(1).is_none());
}

#[test]
fn a_modification_on_another_residue_is_a_difference_not_an_equal_total() {
    // Same total mass, phosphate on a different residue: two differing sites,
    // which one neutral loss cannot explain.
    let parent = modified("PEPSTIDESK", "P1", 0, 3, 79.96633);
    let product = modified("PEPSTIDESK", "P1", 0, 8, 79.96633);
    let db = database(vec![parent, product]);

    let annotations = annotate(
        &[psm(0, 0, 0, 18.0), psm(1, 1, 0, 18.02)],
        &db,
        &settings(0.2),
    );

    assert!(annotations.psm(1).is_none());
    assert_eq!(annotations.files()[0].totals.candidate_pairs, 0);
}

#[test]
fn a_heavier_product_is_not_a_product() {
    let product = modified("PEPTIDE", "P1", 13, 1, 42.010565);
    let db = database(vec![peptide("KRTPEPTIDE", "P1", 10), product]);

    let annotations = annotate(
        &[psm(0, 0, 0, 30.0), psm(1, 1, 0, 30.02)],
        &db,
        &settings(0.2),
    );

    assert!(annotations.psm(1).is_none());
    assert_eq!(annotations.files()[0].totals.candidate_pairs, 0);
}

#[test]
fn a_mass_difference_that_is_no_listed_loss_is_not_flagged() {
    let product = modified("PEPTIDE", "P1", 13, 1, -13.37);
    let db = database(vec![peptide("KRTPEPTIDE", "P1", 10), product]);

    let annotations = annotate(
        &[psm(0, 0, 0, 30.0), psm(1, 1, 0, 30.02)],
        &db,
        &settings(0.2),
    );

    assert!(annotations.psm(1).is_none());
}

#[test]
fn a_product_that_does_not_co_elute_is_not_flagged() {
    let db = database(vec![
        peptide("KPEPTIDE", "P1", 10),
        peptide("PEPTIDE", "P1", 11),
    ]);
    let features = vec![psm(0, 0, 0, 30.0), psm(1, 1, 0, 34.0)];

    let annotations = annotate(&features, &db, &settings(0.2));

    assert!(annotations.psm(1).is_none());
    // It is still a candidate pair, which is what the check counts.
    let stats = &annotations.files()[0].totals;
    assert_eq!(stats.candidate_pairs, 1);
    assert_eq!(stats.coeluting_pairs, 0);
}

#[test]
fn the_same_peptide_is_a_product_at_one_retention_time_and_not_at_another() {
    let db = database(vec![
        peptide("KPEPTIDE", "P1", 10),
        peptide("PEPTIDE", "P1", 11),
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
fn the_parent_is_the_longest_co_eluting_peptide_that_contains_the_product() {
    // A product of a product points at the intact peptide, not at the middle
    // peptide. All three share the C-terminus, so all are suffix products.
    let db = database(vec![
        peptide("KRSTPEPTIDEK", "P1", 4),
        peptide("PEPTIDEK", "P1", 8),
        peptide("TIDEK", "P1", 11),
    ]);
    let features = vec![psm(0, 0, 0, 40.0), psm(1, 1, 0, 40.02), psm(2, 2, 0, 40.04)];

    let annotations = annotate(&features, &db, &settings(0.2));

    assert_eq!(annotations.psm(2).unwrap().parent, PeptideIx(0));
    assert_eq!(annotations.psm(1).unwrap().parent, PeptideIx(0));
}

#[test]
fn only_a_parent_that_co_elutes_is_named() {
    let db = database(vec![
        peptide("KRSTPEPTIDEK", "P1", 4),
        peptide("PEPTIDEK", "P1", 8),
        peptide("TIDEK", "P1", 11),
    ]);
    // The intact peptide elutes elsewhere, so the middle peptide is the parent
    // of the shortest product.
    let features = vec![psm(0, 0, 0, 10.0), psm(1, 1, 0, 40.02), psm(2, 2, 0, 40.04)];

    let annotations = annotate(&features, &db, &settings(0.2));

    assert_eq!(annotations.psm(2).unwrap().parent, PeptideIx(1));
    assert!(annotations.psm(1).is_none());
}

#[test]
fn a_product_in_another_protein_is_not_a_product_of_this_one() {
    let db = database(vec![
        peptide("KPEPTIDE", "P1", 10),
        peptide("PEPTIDE", "P2", 300),
    ]);

    let annotations = annotate(
        &[psm(0, 0, 0, 30.0), psm(1, 1, 0, 30.02)],
        &db,
        &settings(0.2),
    );

    assert!(annotations.psm(1).is_none());
}

#[test]
fn a_product_and_its_parent_must_be_in_the_same_file() {
    let db = database(vec![
        peptide("KPEPTIDE", "P1", 10),
        peptide("PEPTIDE", "P1", 11),
    ]);

    let annotations = annotate(
        &[psm(0, 0, 0, 30.0), psm(1, 1, 1, 30.02)],
        &db,
        &settings(0.2),
    );

    assert!(annotations.psm(1).is_none());
    assert_eq!(annotations.files().len(), 2);
}

#[test]
fn a_candidate_pair_is_found_even_when_the_parent_is_confident_in_another_file() {
    // The parent is identified in file 1 only. The pair is still a candidate,
    // because sequence and mass do not depend on the file, but the product in
    // file 0 has no parent there to co-elute with.
    let db = database(vec![
        peptide("KPEPTIDE", "P1", 10),
        peptide("PEPTIDE", "P1", 11),
    ]);
    let features = vec![psm(0, 1, 0, 30.0), psm(1, 0, 1, 30.0), psm(2, 1, 1, 30.02)];

    let annotations = annotate(&features, &db, &settings(0.2));

    assert_eq!(annotations.candidates(PeptideIx(1)).len(), 1);
    assert!(annotations.psm(0).is_none());
    assert_eq!(annotations.files()[0].totals.candidate_pairs, 0);
    assert!(annotations.psm(2).is_some());
    assert_eq!(annotations.files()[1].totals.candidate_pairs, 1);
}

#[test]
fn unconfident_and_decoy_psms_take_no_part() {
    let mut decoy_parent = peptide("KPEPTIDE", "P1", 10);
    decoy_parent.decoy = true;
    let db = database(vec![
        peptide("KPEPTIDE", "P1", 10),
        peptide("PEPTIDE", "P1", 11),
        decoy_parent,
    ]);

    // A parent above the q-value cutoff is no parent.
    let mut unconfident = psm(0, 0, 0, 30.0);
    unconfident.spectrum_q = 0.5;
    let annotations = annotate(&[unconfident, psm(1, 1, 0, 30.02)], &db, &settings(0.2));
    assert!(annotations.psm(1).is_none());

    // Neither is a decoy, and a decoy PSM is never flagged.
    let mut decoy_psm = psm(2, 2, 0, 30.0);
    decoy_psm.label = -1;
    let annotations = annotate(&[decoy_psm, psm(1, 1, 0, 30.02)], &db, &settings(0.2));
    assert!(annotations.psm(1).is_none());
}

#[test]
fn a_quantified_row_is_judged_on_the_retention_times_its_caller_has() {
    let db = database(vec![
        peptide("KPEPTIDE", "P1", 10),
        peptide("PEPTIDE", "P1", 11),
    ]);
    // The product has no PSM in file 1 at all, which is the shape of a row
    // transferred between runs; its peak there is still judged.
    let features = vec![psm(0, 0, 0, 30.0), psm(1, 1, 0, 34.0), psm(2, 0, 1, 31.0)];

    let annotations = annotate(&features, &db, &settings(0.2));
    let candidates = annotations.candidates(PeptideIx(1));
    assert_eq!(candidates.len(), 1);
    assert_eq!(candidates[0].parent, PeptideIx(0));

    // The PSM at 34.0 in file 0 is not flagged, but a peak apexing at 30.05 is.
    assert!(annotations.psm(1).is_none());
    let quantified = annotations
        .coeluting(0, 30.05, candidates, |_| Some(30.0))
        .expect("the quantified peak co-elutes with the parent's peak");
    assert_eq!(quantified.parent, PeptideIx(0));
    assert!(annotations
        .coeluting(0, 34.0, candidates, |_| Some(30.0))
        .is_none());
    // In file 1 the product was never identified; the parent's apex decides.
    assert!(annotations
        .coeluting(1, 31.02, candidates, |_| Some(31.0))
        .is_some());
    // A parent without a peak in that file cannot be named.
    assert!(annotations
        .coeluting(1, 31.02, candidates, |_| None)
        .is_none());
    // A file the annotation never saw has no tolerance, so no verdict.
    assert!(annotations
        .coeluting(7, 30.05, candidates, |_| Some(30.0))
        .is_none());
}

#[test]
fn the_check_splits_by_placement_and_estimates_what_chance_would_give() {
    // One co-eluting suffix product, plus twelve unrelated PSMs of the same
    // peptide spread over the run, which is what the background comes from.
    let db = database(vec![
        peptide("KPEPTIDE", "P1", 10),
        peptide("PEPTIDE", "P1", 11),
    ]);
    let mut features = vec![psm(0, 0, 0, 30.0), psm(1, 1, 0, 30.02)];
    for step in 0..12 {
        features.push(psm(2 + step, 1, 0, 31.0 + step as f32 * 0.25));
    }

    let annotations = annotate(&features, &db, &settings(0.1));
    let file = &annotations.files()[0];

    assert_eq!(file.file_id, 0);
    assert_eq!(file.rt_tolerance, 0.1);
    assert_eq!(file.rt_tolerance_source, "configured");
    assert_eq!(file.confident_psms, 14);
    assert_eq!(file.flagged_psms, 1);
    assert_eq!(file.totals.placement, "all");
    assert_eq!(file.totals.candidate_pairs, 13);
    assert_eq!(file.totals.coeluting_pairs, 1);
    // Every candidate here is a suffix, so the split repeats the totals.
    assert_eq!(file.placements.len(), 1);
    assert_eq!(file.placements[0].placement, "suffix");
    assert_eq!(file.placements[0].candidate_pairs, 13);
    assert_eq!(file.placements[0].flagged_psms, 1);

    // The histogram grid is fixed and symmetric about zero, and the pair at
    // +0.02 sits in the bin that starts at zero.
    let histogram = &file.totals.rt_difference_histogram;
    assert!(histogram.windows(2).all(|pair| pair[0].bin < pair[1].bin));
    let zero = histogram.iter().find(|entry| entry.bin == 0).unwrap();
    assert_eq!(zero.pairs, 1);
    assert_eq!(zero.start, 0.0);
    assert!((zero.end - HISTOGRAM_BIN).abs() < 1e-6);

    // The estimate is the background band's density over the window's width,
    // and the unrelated PSMs put it far below what co-eluted.
    let background: usize = histogram
        .iter()
        .filter(|entry| {
            let distance = entry.start.abs();
            (BACKGROUND_BAND.start..BACKGROUND_BAND.end).contains(&distance)
        })
        .map(|entry| entry.pairs)
        .sum();
    assert!(background >= 8, "background band holds {background} pairs");
    let density = background as f32 / (2.0 * (BACKGROUND_BAND.end - BACKGROUND_BAND.start));
    assert!((file.totals.expected_by_chance - density * 2.0 * 0.1).abs() < 1e-6);
    assert!(file.totals.expected_by_chance < file.totals.coeluting_pairs as f32);

    // The run-level record is the files' records summed on the same grid.
    let totals = annotations.run_totals();
    assert_eq!(totals.candidate_pairs, file.totals.candidate_pairs);
    assert_eq!(
        totals.rt_difference_histogram,
        file.totals.rt_difference_histogram
    );
    let placements = annotations.run_placements();
    assert_eq!(placements.len(), 1);
    assert_eq!(placements[0].placement, "suffix");
}

#[test]
fn a_difference_beyond_the_histogram_is_counted_apart() {
    let db = database(vec![
        peptide("KPEPTIDE", "P1", 10),
        peptide("PEPTIDE", "P1", 11),
    ]);
    // Six minutes apart, past the histogram's half-width of five.
    let annotations = annotate(
        &[psm(0, 0, 0, 10.0), psm(1, 1, 0, 16.0)],
        &db,
        &settings(0.1),
    );
    let stats = &annotations.files()[0].totals;

    assert_eq!(stats.candidate_pairs, 1);
    assert_eq!(stats.outside_pairs, 1);
    assert!(stats.rt_difference_histogram.is_empty());
    assert_eq!(stats.expected_by_chance, 0.0);
}

#[test]
fn the_tolerance_comes_from_the_spread_of_repeated_peptides() {
    // Twenty observations, 0.00 to 0.19: the median is the value at rank
    // round(0.5 * 19) = 10.
    let mut spreads = (0..20).map(|step| step as f32 / 100.0).collect::<Vec<_>>();
    let (tolerance, source) = derive_rt_tolerance(&mut spreads, 60.0);
    assert_eq!(source, "derived");
    assert!((tolerance - 0.10).abs() < 1e-6, "{tolerance}");

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
        peptide("KPEPTIDE", "P1", 10),
        peptide("PEPTIDE", "P1", 11),
    ]);
    // Enough repeats of the parent to derive a tolerance from the file itself.
    let mut features = vec![psm(0, 1, 0, 30.0)];
    for step in 0..MIN_SPREAD_OBSERVATIONS + 1 {
        features.push(psm(1 + step, 0, 0, 30.0 + step as f32 * 0.01));
    }

    let annotations = annotate(&features, &db, &InSourceSettings::default());
    let stats = &annotations.files()[0];

    assert_eq!(stats.rt_tolerance_source, "derived");
    assert!(stats.rt_tolerance > 0.0);
}

#[test]
fn an_unset_option_follows_the_digest_and_true_or_false_overrides_it() {
    assert_eq!(
        InSourceConfig::resolve(None, true),
        Some(InSourceSettings::default())
    );
    assert_eq!(InSourceConfig::resolve(None, false), None);
    assert_eq!(
        InSourceConfig::resolve(Some(InSourceConfig::Enabled(true)), false),
        Some(InSourceSettings::default())
    );
    assert_eq!(
        InSourceConfig::resolve(Some(InSourceConfig::Enabled(false)), true),
        None
    );
    let configured = InSourceSettings {
        rt_tolerance: Some(0.5),
    };
    assert_eq!(
        InSourceConfig::resolve(Some(InSourceConfig::Settings(configured)), false),
        Some(configured)
    );
}

#[test]
fn a_tolerance_of_zero_or_worse_is_rejected() {
    assert!(InSourceSettings::default().validate().is_ok());
    assert!(InSourceSettings {
        rt_tolerance: Some(0.25)
    }
    .validate()
    .is_ok());
    for bad in [0.0, -1.0, f32::NAN, f32::INFINITY] {
        assert!(
            InSourceSettings {
                rt_tolerance: Some(bad)
            }
            .validate()
            .is_err(),
            "{bad} should be rejected"
        );
    }
}
