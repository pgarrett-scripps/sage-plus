use super::*;
use crate::database::Builder;
use crate::enzyme::Digest;
use crate::ion_series::IonSeries;
use crate::modification::{ModificationDefinition, ModificationSpecificity, NeutralLossMode};
use crate::peptide::Peptide;
use std::{collections::HashMap, sync::Arc};

#[test]
fn score_ordering_uses_hyperscore() {
    let database_first = Score {
        peptide: PeptideIx(0),
        hyperscore: 1.0,
        ..Default::default()
    };
    let database_last = Score {
        peptide: PeptideIx(10_000),
        hyperscore: 100.0,
        ..Default::default()
    };

    assert!(database_last > database_first);
    assert_eq!(
        database_last.partial_cmp(&database_first),
        Some(database_last.cmp(&database_first))
    );

    let mut candidates = vec![database_first, database_last];
    bounded_min_heapify(&mut candidates, 1);
    assert_eq!(candidates[0].peptide, database_last.peptide);
}

#[test]
fn exact_prefilter_preserves_tied_isobaric_scoring() {
    let builder = Builder {
        generate_decoys: Some(false),
        ..Builder::default()
    };
    let parameters = builder.make_parameters();
    let peptides = parameters.peptides_from_tsv("sequence\nLEKPI\nPEILK\nPELIK\n");
    let expected = peptides
        .iter()
        .find(|peptide| peptide.sequence.as_ref() == b"PEILK")
        .unwrap()
        .clone();
    let mut masses = [Kind::B, Kind::Y]
        .into_iter()
        .flat_map(|kind| IonSeries::new(&expected, kind))
        .map(|ion| ion.monoisotopic_mass)
        .collect::<Vec<_>>();
    masses.sort_by(f32::total_cmp);
    let query = ProcessedSpectrum {
        level: 2,
        id: "tied-isobaric".into(),
        precursors: vec![Precursor {
            mz: expected.monoisotopic / 2.0 + PROTON,
            charge: Some(2),
            ..Precursor::default()
        }],
        intensities: vec![1.0; masses.len()],
        charges: vec![1; masses.len()],
        total_ion_current: masses.len() as f32,
        masses,
        ..ProcessedSpectrum::default()
    };

    let full = parameters.clone().build_from_peptides(peptides);
    let make_scorer = |database| Scorer {
        db: database,
        precursor_tol: Tolerance::Da(-0.01, 0.01),
        fragment_tol: Tolerance::Da(-0.01, 0.01),
        min_matched_peaks: 1,
        min_isotope_err: 0,
        max_isotope_err: 0,
        min_precursor_charge: 2,
        max_precursor_charge: 2,
        override_precursor_charge: false,
        max_fragment_charge: Some(1),
        chimera: false,
        report_psms: 2,
        wide_window: false,
        annotate_matches: false,
        mass_shift_ppm: crate::ambiguity::DEFAULT_MASS_SHIFT_PPM,
        score_type: ScoreType::SageHyperScore,
    };
    let full_scorer = make_scorer(&full);
    let full_features = full_scorer.score(&query);
    assert_eq!(
        full[full_features[0].peptide_idx].sequence.as_ref(),
        b"PEILK"
    );

    let keep = AtomicBitSet::new(full.peptides.len());
    full_scorer.exact_prefilter(&query, &keep);
    let survivors = full
        .peptides
        .iter()
        .enumerate()
        .filter(|(index, _)| keep.contains(*index))
        .map(|(_, peptide)| peptide.clone())
        .collect::<Vec<_>>();
    assert!(survivors.len() < full.peptides.len());
    let reduced = parameters.build_from_peptides(survivors);
    let reduced_features = make_scorer(&reduced).score(&query);

    assert_eq!(
        full[full_features[0].peptide_idx].sequence,
        reduced[reduced_features[0].peptide_idx].sequence
    );
    assert_eq!(full_features[0].hyperscore, reduced_features[0].hyperscore);
    assert_eq!(
        full_features[0].scored_candidates,
        reduced_features[0].scored_candidates
    );
    assert_eq!(
        full_features[0].matched_peaks,
        reduced_features[0].matched_peaks
    );
}

#[test]
fn longest_series() {
    let mut run = Run::default();

    run.matched(1);
    run.matched(2);
    run.matched(3);
    run.matched(3);
    run.matched(3);

    assert_eq!(run.length, 3);
    assert_eq!(run.longest, 3);

    run.matched(5);
    run.matched(5);
    assert_eq!(run.length, 1);
    assert_eq!(run.longest, 3);
    run.matched(6);
    assert_eq!(run.length, 2);
}

#[test]
fn test_max_fragment_charge() {
    assert_eq!(max_fragment_charge(None, 1), 2);
    assert_eq!(max_fragment_charge(None, 2), 2);
    assert_eq!(max_fragment_charge(None, 3), 3);
    assert_eq!(max_fragment_charge(None, 4), 4);
    assert_eq!(max_fragment_charge(Some(1), 2), 2);
    assert_eq!(max_fragment_charge(Some(1), 3), 2);
    assert_eq!(max_fragment_charge(Some(2), 4), 3);
    assert_eq!(max_fragment_charge(Some(4), 1), 2);
}

#[test]
fn fragment_match_index_enforces_known_charge_and_infers_unknown_charge() {
    let query = ProcessedSpectrum {
        masses: vec![500.0, 1_000.0],
        intensities: vec![10.0, 20.0],
        charges: vec![1, 2],
        charge_is_known: vec![false, true],
        ..ProcessedSpectrum::default()
    };
    let index = FragmentMatchIndex::new(&query, 3);

    assert_eq!(
        index.select_peak(&query, 1_000.0, 2, Tolerance::Da(-0.01, 0.01)),
        Some(1)
    );
    assert_eq!(
        index.select_peak(&query, 1_000.0, 1, Tolerance::Da(-0.01, 0.01)),
        None
    );

    let two_known_charges = ProcessedSpectrum {
        masses: vec![1_000.0, 1_000.0],
        intensities: vec![15.0, 20.0],
        charges: vec![1, 2],
        charge_is_known: vec![true, true],
        ..ProcessedSpectrum::default()
    };
    let index = FragmentMatchIndex::new(&two_known_charges, 3);
    assert_eq!(
        index.select_peak(&two_known_charges, 1_000.0, 1, Tolerance::Da(-0.01, 0.01)),
        Some(0)
    );
    assert_eq!(
        index.select_peak(&two_known_charges, 1_000.0, 2, Tolerance::Da(-0.01, 0.01)),
        Some(1)
    );

    let unknown_only = ProcessedSpectrum {
        masses: vec![500.0],
        intensities: vec![10.0],
        charges: vec![1],
        charge_is_known: vec![false],
        ..ProcessedSpectrum::default()
    };
    let index = FragmentMatchIndex::new(&unknown_only, 3);
    assert_eq!(
        index.select_peak(&unknown_only, 1_000.0, 2, Tolerance::Da(-0.01, 0.01)),
        Some(0)
    );

    let da_scaled_unknown = ProcessedSpectrum {
        masses: vec![500.0075],
        intensities: vec![10.0],
        charges: vec![1],
        charge_is_known: vec![false],
        ..ProcessedSpectrum::default()
    };
    let index = FragmentMatchIndex::new(&da_scaled_unknown, 3);
    assert_eq!(
        index.select_peak(&da_scaled_unknown, 1_000.0, 2, Tolerance::Da(-0.01, 0.01)),
        Some(0)
    );

    let da_unscaled_known = ProcessedSpectrum {
        masses: vec![1_000.015],
        intensities: vec![10.0],
        charges: vec![2],
        charge_is_known: vec![true],
        ..ProcessedSpectrum::default()
    };
    let index = FragmentMatchIndex::new(&da_unscaled_known, 3);
    assert_eq!(
        index.select_peak(&da_unscaled_known, 1_000.0, 2, Tolerance::Da(-0.01, 0.01)),
        None
    );
}

#[test]
fn equal_nonzero_isotope_bounds_are_honored() {
    let peptide = crate::peptide::Peptide::try_from(Digest {
        sequence: "PEPTIDER".into(),
        protein: Arc::from("protein"),
        ..Digest::default()
    })
    .unwrap();
    let fragment_masses = [Kind::B, Kind::Y]
        .into_iter()
        .flat_map(|kind| IonSeries::new(&peptide, kind))
        .map(|ion| ion.monoisotopic_mass)
        .collect::<Vec<_>>();

    let parameters = Builder::default().make_parameters();
    let database = parameters.build_from_peptides(vec![peptide.clone()]);
    let precursor_charge = 2;
    let precursor = Precursor {
        mz: (peptide.monoisotopic + NEUTRON) / precursor_charge as f32 + PROTON,
        charge: Some(precursor_charge),
        ..Precursor::default()
    };
    let mut query = ProcessedSpectrum {
        level: 2,
        id: "isotope-test".into(),
        precursors: vec![precursor],
        masses: fragment_masses.clone(),
        intensities: vec![1.0; fragment_masses.len()],
        charges: vec![1; fragment_masses.len()],
        total_ion_current: fragment_masses.len() as f32,
        ..ProcessedSpectrum::default()
    };
    query.masses.sort_by(f32::total_cmp);

    let scorer = Scorer {
        db: &database,
        precursor_tol: Tolerance::Da(-0.01, 0.01),
        fragment_tol: Tolerance::Da(-0.01, 0.01),
        min_matched_peaks: 1,
        min_isotope_err: 1,
        max_isotope_err: 1,
        min_precursor_charge: 2,
        max_precursor_charge: 2,
        override_precursor_charge: false,
        max_fragment_charge: Some(1),
        chimera: false,
        report_psms: 1,
        wide_window: false,
        annotate_matches: false,
        mass_shift_ppm: crate::ambiguity::DEFAULT_MASS_SHIFT_PPM,
        score_type: ScoreType::SageHyperScore,
    };

    let features = scorer.score(&query);
    assert_eq!(features.len(), 1);
    assert_eq!(features[0].isotope_error, NEUTRON);
}

#[test]
fn isotope_offsets_are_honored_for_labeled_precursors() {
    let builder: Builder = serde_json::from_value(serde_json::json!({
        "generate_decoys": false,
        "static_mods": {
            "R": {
                "mass": 0.0,
                "name": "SILAC-R",
                "channel_offsets": {"light": 0.0, "heavy": 10.008269}
            }
        }
    }))
    .unwrap();
    let parameters = builder.make_parameters();
    parameters.validate_channels().unwrap();
    let peptides = parameters.peptides_from_tsv("sequence\nPEPTIDER\n");
    let heavy = peptides
        .iter()
        .find(|peptide| peptide.label_channel.as_deref() == Some("heavy"))
        .unwrap()
        .clone();
    let fragment_masses = [Kind::B, Kind::Y]
        .into_iter()
        .flat_map(|kind| IonSeries::new(&heavy, kind))
        .map(|ion| ion.monoisotopic_mass)
        .collect::<Vec<_>>();
    let database = parameters.build_from_peptides(peptides);
    let precursor_charge = 2;
    let precursor = Precursor {
        mz: (heavy.monoisotopic + NEUTRON) / precursor_charge as f32 + PROTON,
        charge: Some(precursor_charge),
        ..Precursor::default()
    };
    let mut query = ProcessedSpectrum {
        level: 2,
        id: "labeled-isotope-test".into(),
        precursors: vec![precursor],
        masses: fragment_masses.clone(),
        intensities: vec![1.0; fragment_masses.len()],
        charges: vec![1; fragment_masses.len()],
        total_ion_current: fragment_masses.len() as f32,
        ..ProcessedSpectrum::default()
    };
    query.masses.sort_by(f32::total_cmp);
    let scorer = Scorer {
        db: &database,
        precursor_tol: Tolerance::Da(-0.01, 0.01),
        fragment_tol: Tolerance::Da(-0.01, 0.01),
        min_matched_peaks: 1,
        min_isotope_err: 1,
        max_isotope_err: 1,
        min_precursor_charge: 2,
        max_precursor_charge: 2,
        override_precursor_charge: false,
        max_fragment_charge: Some(1),
        chimera: false,
        report_psms: 1,
        wide_window: false,
        annotate_matches: false,
        mass_shift_ppm: crate::ambiguity::DEFAULT_MASS_SHIFT_PPM,
        score_type: ScoreType::SageHyperScore,
    };

    let features = scorer.score(&query);
    assert_eq!(features.len(), 1);
    assert_eq!(features[0].isotope_error, NEUTRON);
    assert_eq!(
        database[features[0].peptide_idx].label_channel.as_deref(),
        Some("heavy")
    );
}

#[test]
fn neutral_loss_alternatives_count_once_per_cleavage_and_charge() {
    let modification = Arc::new(ModificationDefinition {
        mass: 20.0,
        name: Some(Arc::from("TestMod")),
        neutral_losses: Arc::from([10.0]),
        neutral_loss_mode: NeutralLossMode::Optional,
        channel_offsets: Arc::default(),
    });
    let peptide = Peptide::try_from(Digest {
        sequence: "AMK".into(),
        ..Default::default()
    })
    .unwrap()
    .apply(
        &[(
            ModificationSpecificity::Residue(b'M'),
            modification,
            Some(1),
        )],
        &HashMap::default(),
        1,
        None,
    )
    .into_iter()
    .find(|peptide| peptide.to_string().contains("TestMod"))
    .unwrap();

    let group = IonGroupSeries::new(&peptide, Kind::B).nth(1).unwrap();
    assert_eq!(group.variants.len(), 2);
    let mut variants = group.variants;
    variants.sort_by(|a, b| a.monoisotopic_mass.total_cmp(&b.monoisotopic_mass));

    let db = IndexedDatabase {
        peptides: vec![peptide],
        ion_kinds: vec![Kind::B],
        ..Default::default()
    };
    let scorer = Scorer {
        db: &db,
        precursor_tol: Tolerance::Da(-0.01, 0.01),
        fragment_tol: Tolerance::Da(-0.01, 0.01),
        min_matched_peaks: 1,
        min_isotope_err: 0,
        max_isotope_err: 0,
        min_precursor_charge: 2,
        max_precursor_charge: 2,
        override_precursor_charge: false,
        max_fragment_charge: Some(1),
        chimera: false,
        report_psms: 1,
        wide_window: false,
        annotate_matches: true,
        mass_shift_ppm: crate::ambiguity::DEFAULT_MASS_SHIFT_PPM,
        score_type: ScoreType::SageHyperScore,
    };
    let query = ProcessedSpectrum {
        masses: variants
            .iter()
            .map(|variant| variant.monoisotopic_mass)
            .collect(),
        intensities: vec![100.0, 10.0],
        charges: vec![1, 1],
        total_ion_current: 110.0,
        ..Default::default()
    };
    let pre_score = PreScore {
        peptide: PeptideIx(0),
        precursor_charge: 2,
        ..Default::default()
    };

    let (score, fragments, _) = scorer.score_candidate(&query, &pre_score, true);
    assert_eq!(score.matched_b, 1);
    assert_eq!(score.summed_b, 100.0);
    let fragments = fragments.unwrap();
    assert_eq!(fragments.fragment_ordinals.len(), 1);
    assert_eq!(fragments.neutral_losses, vec![10.0]);

    let deferred = scorer.annotate_candidate(
        &query,
        &Feature {
            peptide_idx: PeptideIx(0),
            charge: 2,
            ..Default::default()
        },
    );
    assert_eq!(deferred.kinds, fragments.kinds);
    assert_eq!(deferred.charges, fragments.charges);
    assert_eq!(deferred.fragment_ordinals, fragments.fragment_ordinals);
    assert_eq!(deferred.intensities, fragments.intensities);
    assert_eq!(deferred.mz_calculated, fragments.mz_calculated);
    assert_eq!(deferred.mz_experimental, fragments.mz_experimental);
    assert_eq!(deferred.neutral_losses, fragments.neutral_losses);
}

#[test]
fn deferred_chimera_annotation_replays_filtered_preceding_ranks() {
    let peptide = Peptide::try_from(Digest {
        sequence: "PEPTIDER".into(),
        ..Default::default()
    })
    .unwrap();
    let mut peaks = [Kind::B, Kind::Y]
        .into_iter()
        .flat_map(|kind| IonSeries::new(&peptide, kind))
        .map(|ion| (ion.monoisotopic_mass, 100.0, 1))
        .collect::<Vec<_>>();
    peaks.sort_unstable_by(|a, b| a.0.total_cmp(&b.0));
    let query = ProcessedSpectrum {
        level: 2,
        masses: peaks.iter().map(|peak| peak.0).collect(),
        intensities: peaks.iter().map(|peak| peak.1).collect(),
        charges: peaks.iter().map(|peak| peak.2).collect(),
        total_ion_current: peaks.iter().map(|peak| peak.1).sum(),
        ..Default::default()
    };
    let database = IndexedDatabase {
        peptides: vec![peptide],
        ion_kinds: vec![Kind::B, Kind::Y],
        ..Default::default()
    };
    let scorer = |chimera| Scorer {
        db: &database,
        precursor_tol: Tolerance::Da(-0.01, 0.01),
        fragment_tol: Tolerance::Da(-0.01, 0.01),
        min_matched_peaks: 1,
        min_isotope_err: 0,
        max_isotope_err: 0,
        min_precursor_charge: 2,
        max_precursor_charge: 2,
        override_precursor_charge: false,
        max_fragment_charge: Some(1),
        chimera,
        report_psms: 2,
        wide_window: false,
        annotate_matches: false,
        mass_shift_ppm: crate::ambiguity::DEFAULT_MASS_SHIFT_PPM,
        score_type: ScoreType::SageHyperScore,
    };
    let rank_one = Feature {
        peptide_idx: PeptideIx(0),
        charge: 2,
        rank: 1,
        ..Default::default()
    };
    let rank_two = Feature {
        rank: 2,
        ..rank_one.clone()
    };
    let features = [&rank_one, &rank_two];

    let replayed = scorer(true).annotate_ranked_candidates(&query, &features, &[false, true]);
    assert!(replayed[0].is_none());
    assert_eq!(
        replayed[1].as_ref().unwrap().fragment_ordinals.len(),
        0,
        "rank one must remove its peaks even when it is filtered from output"
    );

    let independent = scorer(false).annotate_ranked_candidates(&query, &features, &[false, true]);
    assert!(!independent[1]
        .as_ref()
        .unwrap()
        .fragment_ordinals
        .is_empty());
}
