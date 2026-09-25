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
        mass_recalibration: None,
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
        mass_recalibration: None,
    };

    let features = scorer.score(&query);
    assert_eq!(features.len(), 1);
    assert_eq!(features[0].isotope_error, NEUTRON);
    // The matched isotope error is not a residual modification mass.
    assert_eq!(features[0].mass_shift, 0.0);
    assert_eq!(features[0].ambiguity_sequence, "PEPTIDER");
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
        mass_recalibration: None,
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
        mass_recalibration: None,
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
        mass_recalibration: None,
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

mod mass_offsets {
    use super::*;
    use crate::database::{MassOffsetAssignment, Parameters};
    use crate::fasta::Fasta;
    use crate::modification::{SearchMode, SiteMode, VarModEntry, VariableModification};
    use crate::peptide::Site;
    use crate::ptm_library::{PtmLibrary, PtmLibrarySite};

    const PHOSPHO: f32 = 79.966_33;
    const FASTA: &str = ">P1\nMAGSPEPTSIDEKLLSAYGNRWTTPEGSAR\n>P2\nGGSTVLAPEDKAAAAR\n";

    fn phospho(search_mode: SearchMode, site_mode: SiteMode) -> Vec<VarModEntry> {
        vec![VarModEntry::Detailed(VariableModification {
            mass: PHOSPHO,
            max_count: Some(1),
            name: Some("Phospho".into()),
            neutral_losses: Vec::new(),
            neutral_loss_mode: NeutralLossMode::Optional,
            site_mode,
            search_mode,
            channel_offsets: Default::default(),
        })]
    }

    fn parameters(search_mode: SearchMode, site_mode: SiteMode) -> Parameters {
        let mut variable_mods = HashMap::new();
        for residue in ["S", "T", "Y"] {
            variable_mods.insert(residue.to_string(), phospho(search_mode, site_mode));
        }
        Builder {
            max_variable_mods: Some(1),
            variable_mods: Some(variable_mods),
            static_mods: Some(HashMap::from([(
                "C".to_string(),
                crate::modification::StaticModEntry::Mass(57.021_464),
            )])),
            peptide_min_mass: Some(300.0),
            enzyme: Some(crate::database::EnzymeBuilder {
                min_len: Some(5),
                ..Default::default()
            }),
            ..Builder::default()
        }
        .make_parameters()
    }

    fn database(search_mode: SearchMode) -> IndexedDatabase {
        let parameters = parameters(search_mode, SiteMode::Exhaustive);
        let fasta = Fasta::parse(FASTA.into(), "rev_", true).unwrap();
        parameters.clone().build(fasta)
    }

    fn spectrum(peptide: &Peptide) -> ProcessedSpectrum {
        let mut masses = [Kind::B, Kind::Y]
            .into_iter()
            .flat_map(|kind| IonSeries::new(peptide, kind))
            .map(|ion| ion.monoisotopic_mass)
            .collect::<Vec<_>>();
        masses.sort_by(f32::total_cmp);
        ProcessedSpectrum {
            level: 2,
            id: "offset".into(),
            precursors: vec![Precursor {
                mz: peptide.monoisotopic / 2.0 + PROTON,
                charge: Some(2),
                ..Precursor::default()
            }],
            intensities: vec![10.0; masses.len()],
            charges: vec![1; masses.len()],
            total_ion_current: 10.0 * masses.len() as f32,
            masses,
            ..ProcessedSpectrum::default()
        }
    }

    fn scorer(database: &IndexedDatabase, chimera: bool) -> Scorer<'_> {
        Scorer {
            db: database,
            precursor_tol: Tolerance::Ppm(-10.0, 10.0),
            fragment_tol: Tolerance::Ppm(-10.0, 10.0),
            min_matched_peaks: 4,
            min_isotope_err: -1,
            max_isotope_err: 2,
            min_precursor_charge: 2,
            max_precursor_charge: 3,
            override_precursor_charge: false,
            max_fragment_charge: Some(1),
            chimera,
            report_psms: 2,
            wide_window: false,
            annotate_matches: false,
            mass_shift_ppm: crate::ambiguity::DEFAULT_MASS_SHIFT_PPM,
            score_type: ScoreType::SageHyperScore,
            mass_recalibration: None,
        }
    }

    /// The target peptidoform from the expanded database, phosphorylated at
    /// the second S/T/Y candidate (PEPTS*IDEK).
    fn expanded_target(database: &IndexedDatabase) -> Peptide {
        database
            .peptides
            .iter()
            .find(|peptide| peptide.to_string() == "MAGSPEPTS[Phospho]IDEK")
            .expect("expanded database contains the phosphopeptide")
            .clone()
    }

    #[test]
    fn offset_modifications_are_not_indexed() {
        let expanded = database(SearchMode::Database);
        let offset = database(SearchMode::MassOffset);
        let unmodified = parameters(SearchMode::Database, SiteMode::Exhaustive);
        let mut plain = unmodified.clone();
        plain.variable_mods.clear();
        let plain = plain.build(Fasta::parse(FASTA.into(), "rev_", true).unwrap());

        assert_eq!(offset.peptides.len(), plain.peptides.len());
        assert_eq!(offset.fragments.len(), plain.fragments.len());
        assert!(expanded.peptides.len() > offset.peptides.len());
        assert_eq!(offset.mass_offsets.len(), 1);
        assert_eq!(offset.mass_offsets[0].specificities.len(), 3);
        // Localization still knows the offset site rules.
        assert_eq!(offset.potential_mods, expanded.potential_mods);
    }

    #[test]
    fn offset_search_matches_expanded_search() {
        let expanded = database(SearchMode::Database);
        let offset = database(SearchMode::MassOffset);
        let target = expanded_target(&expanded);
        let query = spectrum(&target);

        let expanded_hits = scorer(&expanded, false).score(&query);
        let mut offset_hits = scorer(&offset, false).score(&query);
        assert_eq!(
            expanded[expanded_hits[0].peptide_idx].to_string(),
            "MAGSPEPTS[Phospho]IDEK"
        );
        let resolved = offset.resolve_peptide(&offset_hits[0]).into_owned();
        assert_eq!(resolved.to_string(), "MAGSPEPTS[Phospho]IDEK");
        assert_eq!(
            offset_hits[0].mass_offset,
            Some(MassOffsetAssignment {
                offset: 0,
                site: Site::Sequence(8),
            })
        );
        assert_eq!(offset_hits[0].hyperscore, expanded_hits[0].hyperscore);
        assert_eq!(offset_hits[0].matched_peaks, expanded_hits[0].matched_peaks);
        assert_eq!(offset_hits[0].calcmass, expanded_hits[0].calcmass);
        // Tied placement isomers compete exactly as in the expanded database.
        assert_eq!(offset_hits[0].delta_next, expanded_hits[0].delta_next);
        let decoy = expanded
            .peptides
            .iter()
            .find(|peptide| peptide.decoy && peptide.modification_at(4) != 0.0)
            .unwrap();
        let base = offset
            .peptides
            .iter()
            .find(|peptide| peptide.decoy && peptide.sequence == decoy.sequence)
            .unwrap();
        let placed = base.with_mass_offset(Site::Sequence(4), &offset.mass_offsets[0].definition);
        assert_eq!(placed.monoisotopic.to_bits(), decoy.monoisotopic.to_bits());
        // The assigned offset is not residual precursor error.
        assert!(offset_hits[0].delta_mass.abs() < 1.0);
        assert!((offset_hits[0].expmass - offset_hits[0].calcmass).abs() < 0.01);

        // Materialized peptides behave like indexed ones for downstream steps.
        let mut offset = offset;
        offset.materialize_mass_offsets(&mut offset_hits);
        let peptide = &offset[offset_hits[0].peptide_idx];
        assert!(offset_hits[0].peptide_idx.0 as usize >= offset.peptides.len());
        assert_eq!(peptide.to_string(), "MAGSPEPTS[Phospho]IDEK");
        assert_eq!(peptide.proteins.as_slice(), &[Arc::<str>::from("P1")]);
        assert_eq!(peptide.protein_sites[0].start, Some(0));
        assert_eq!(peptide.monoisotopic, resolved.monoisotopic);
        let materialized = offset.offset_peptides.len();
        offset.materialize_mass_offsets(&mut offset_hits);
        assert_eq!(offset.offset_peptides.len(), materialized);
    }

    #[test]
    fn offset_placements_are_localized_to_the_same_site() {
        let expanded = database(SearchMode::Database);
        let mut offset = database(SearchMode::MassOffset);
        let query = spectrum(&expanded_target(&expanded));
        let mut hits = scorer(&offset, false).score(&query);
        offset.materialize_mass_offsets(&mut hits);

        let localize = |database: &IndexedDatabase, feature: &Feature| {
            crate::ptm::localize(
                &database[feature.peptide_idx],
                &query,
                &database.ion_kinds,
                &database.potential_mods,
                Tolerance::Ppm(-10.0, 10.0),
                Some(1),
                feature.charge,
            )
        };
        let expanded_hits = scorer(&expanded, false).score(&query);
        let expected = localize(&expanded, &expanded_hits[0]);
        let observed = localize(&offset, &hits[0]);
        assert_eq!(observed, expected);
        assert_eq!(observed.mods[0].best_sites[0].position, 8);
        assert_eq!(observed.mods[0].label.as_deref(), Some("Phospho"));
    }

    #[test]
    fn unshifted_peptides_still_compete_with_offsets() {
        let expanded = database(SearchMode::Database);
        let offset = database(SearchMode::MassOffset);
        let target = expanded
            .peptides
            .iter()
            .find(|peptide| peptide.to_string() == "GGSTVLAPEDK")
            .unwrap()
            .clone();
        let hits = scorer(&offset, false).score(&spectrum(&target));
        assert_eq!(hits[0].mass_offset, None);
        assert_eq!(offset[hits[0].peptide_idx].to_string(), "GGSTVLAPEDK");
    }

    #[test]
    fn offsets_only_skips_the_unshifted_hypothesis() {
        let expanded = database(SearchMode::Database);
        let mut offset = database(SearchMode::MassOffset);
        offset.offsets_only = true;
        // An unmodified precursor is no longer searched...
        let plain = expanded
            .peptides
            .iter()
            .find(|peptide| peptide.to_string() == "GGSTVLAPEDK")
            .unwrap()
            .clone();
        let hits = scorer(&offset, false).score(&spectrum(&plain));
        assert!(hits.iter().all(|hit| hit.mass_offset.is_some()));
        // ...while offset precursors are found as before.
        let hits = scorer(&offset, false).score(&spectrum(&expanded_target(&expanded)));
        assert_eq!(
            offset.resolve_peptide(&hits[0]).to_string(),
            "MAGSPEPTS[Phospho]IDEK"
        );
    }

    #[test]
    fn exact_prefilter_keeps_offset_base_peptides() {
        let expanded = database(SearchMode::Database);
        let offset = database(SearchMode::MassOffset);
        let query = spectrum(&expanded_target(&expanded));
        let keep = AtomicBitSet::new(offset.peptides.len());
        scorer(&offset, false).exact_prefilter(&query, &keep);
        let base = offset
            .peptides
            .iter()
            .position(|peptide| peptide.to_string() == "MAGSPEPTSIDEK")
            .unwrap();
        assert!(keep.contains(base));

        let survivors = offset
            .peptides
            .iter()
            .enumerate()
            .filter(|(index, _)| keep.contains(*index))
            .map(|(_, peptide)| peptide.clone())
            .collect::<Vec<_>>();
        assert!(survivors.len() < offset.peptides.len());
        let reduced =
            parameters(SearchMode::MassOffset, SiteMode::Exhaustive).build_from_peptides(survivors);
        let full = scorer(&offset, false).score(&query);
        let filtered = scorer(&reduced, false).score(&query);
        assert_eq!(
            offset.resolve_peptide(&full[0]).to_string(),
            reduced.resolve_peptide(&filtered[0]).to_string()
        );
        assert_eq!(full[0].hyperscore, filtered[0].hyperscore);
        assert_eq!(full[0].matched_peaks, filtered[0].matched_peaks);
    }

    #[test]
    fn shifted_fragment_lookups_retrieve_offset_candidates() {
        let expanded = database(SearchMode::Database);
        let offset = database(SearchMode::MassOffset);
        let target = expanded_target(&expanded);
        let base = offset
            .peptides
            .iter()
            .position(|peptide| peptide.to_string() == "MAGSPEPTSIDEK")
            .unwrap();
        let unshifted = [Kind::B, Kind::Y]
            .into_iter()
            .flat_map(|kind| IonSeries::new(&offset.peptides[base], kind))
            .map(|ion| ion.monoisotopic_mass)
            .collect::<Vec<_>>();
        // Keep only fragments that carry the modification.
        let mut query = spectrum(&target);
        let retained = query
            .masses
            .iter()
            .map(|mass| !unshifted.iter().any(|base| (base - mass).abs() < 0.01))
            .collect::<Vec<_>>();
        let mut keep_flags = retained.iter();
        query.masses.retain(|_| *keep_flags.next().unwrap());
        query.intensities.truncate(query.masses.len());
        query.charges.truncate(query.masses.len());
        assert!(query.masses.len() >= 4);

        let keep = AtomicBitSet::new(offset.peptides.len());
        scorer(&offset, false).exact_prefilter(&query, &keep);
        assert!(keep.contains(base));
        let hits = scorer(&offset, false).score(&query);
        assert_eq!(
            offset.resolve_peptide(&hits[0]).to_string(),
            "MAGSPEPTS[Phospho]IDEK"
        );
    }

    #[test]
    fn chimera_removes_peaks_of_the_placed_peptidoform() {
        let expanded = database(SearchMode::Database);
        let offset = database(SearchMode::MassOffset);
        let query = spectrum(&expanded_target(&expanded));
        let hits = scorer(&offset, true).score(&query);
        assert_eq!(
            offset.resolve_peptide(&hits[0]).to_string(),
            "MAGSPEPTS[Phospho]IDEK"
        );
        // All peaks were explained by rank one, so nothing remains to match.
        assert_eq!(hits.len(), 1);
    }

    #[test]
    fn occupied_sites_do_not_receive_offsets() {
        let offset = database(SearchMode::MassOffset);
        let peptide = Peptide::try_from(Digest {
            sequence: "SCTK".into(),
            ..Default::default()
        })
        .unwrap()
        .with_mass_offset(Site::Sequence(0), &offset.mass_offsets[0].definition);
        let sites = offset.mass_offset_sites(&peptide, &offset.mass_offsets[0]);
        assert_eq!(sites, vec![Site::Sequence(2)]);
    }

    #[test]
    fn library_offsets_use_library_sites_and_mirror_them_on_decoys() {
        let mut parameters = parameters(SearchMode::MassOffset, SiteMode::Library);
        parameters.loaded_ptm_library = Some(Arc::new(PtmLibrary::new(vec![PtmLibrarySite {
            attachment: Default::default(),
            protein: "P1".into(),
            position: 8,
            residue: b'S',
            modification: "Phospho".into(),
        }])));
        let database = parameters.build(Fasta::parse(FASTA.into(), "rev_", true).unwrap());
        let offset = &database.mass_offsets[0];
        let target = database
            .peptides
            .iter()
            .find(|peptide| peptide.to_string() == "MAGSPEPTSIDEK")
            .unwrap();
        assert_eq!(
            database.mass_offset_sites(target, offset),
            vec![Site::Sequence(8)]
        );
        let decoy = database
            .peptides
            .iter()
            .find(|peptide| peptide.decoy && peptide.sequence == target.reverse().sequence)
            .unwrap();
        // Internal residues are mirrored: 12 - 8 = 4.
        assert_eq!(decoy.sequence[4], b'S');
        assert_eq!(
            database.mass_offset_sites(decoy, offset),
            vec![Site::Sequence(4)]
        );
        let other = database
            .peptides
            .iter()
            .find(|peptide| peptide.to_string() == "GGSTVLAPEDK")
            .unwrap();
        assert!(database.mass_offset_sites(other, offset).is_empty());
    }

    #[test]
    fn duplicate_offset_and_indexed_peptidoforms_are_reported_once() {
        let mut parameters = parameters(SearchMode::Database, SiteMode::Exhaustive);
        let mut unnamed = phospho(SearchMode::MassOffset, SiteMode::Exhaustive);
        if let VarModEntry::Detailed(modification) = &mut unnamed[0] {
            modification.name = None;
        }
        parameters
            .variable_mods
            .insert(ModificationSpecificity::Residue(b'S'), {
                let mut entries =
                    parameters.variable_mods[&ModificationSpecificity::Residue(b'S')].clone();
                if let VarModEntry::Detailed(modification) = &mut entries[0] {
                    modification.name = None;
                }
                entries.extend(unnamed);
                entries
            });
        let database = parameters.build(Fasta::parse(FASTA.into(), "rev_", true).unwrap());
        assert_eq!(database.mass_offsets.len(), 1);
        let target = database
            .peptides
            .iter()
            .find(|peptide| {
                peptide.sequence.as_ref() == b"MAGSPEPTSIDEK" && peptide.modification_at(8) != 0.0
            })
            .unwrap()
            .clone();
        let hits = scorer(&database, false).score(&spectrum(&target));
        assert_eq!(hits[0].mass_offset, None, "indexed form wins exact ties");
        assert!(hits.len() < 2 || hits[1].hyperscore < hits[0].hyperscore);
    }

    /// Spectra shifted beyond the search tolerance are recovered by a
    /// search-time correction, while reported errors stay raw.
    #[test]
    fn recalibration_recovers_shifted_spectrum_and_keeps_raw_errors() {
        use crate::mass_recalibration::{
            FileMassCorrection, GroupMassCorrection, MassErrorModel, MassModelAxes, MassModelKind,
            MassRecalibration,
        };
        use crate::spectrum::{AcquisitionGroup, Activation, MassAnalyzer};
        let expanded = database(SearchMode::Database);
        let target = expanded_target(&expanded);
        let mut query = spectrum(&target);
        query.acquisition = AcquisitionGroup {
            analyzer: MassAnalyzer::Orbitrap,
            activation: Activation::Hcd,
        };
        let ion_trap = AcquisitionGroup {
            analyzer: MassAnalyzer::IonTrap,
            activation: Activation::Cid,
        };
        let shift = 1.0 + 15e-6;
        for precursor in &mut query.precursors {
            precursor.mz *= shift;
        }
        for mass in &mut query.masses {
            *mass = (*mass + PROTON) * shift - PROTON;
        }

        let plain = scorer(&expanded, false).score(&query);
        assert!(plain
            .first()
            .is_none_or(|hit| expanded[hit.peptide_idx].to_string() != "MAGSPEPTS[Phospho]IDEK"));

        let offset = MassErrorModel {
            kind: MassModelKind::Static,
            axes: MassModelAxes::None,
            intercept_ppm: 15.0,
            rt: None,
            mz: None,
            max_abs_ppm: 20.0,
        };
        let recalibration = MassRecalibration {
            files: vec![FileMassCorrection {
                precursor: Some(offset.clone()),
                fragment: vec![
                    GroupMassCorrection {
                        group: query.acquisition,
                        model: Some(offset.clone()),
                    },
                    // Another analyzer's model must not touch this spectrum.
                    GroupMassCorrection {
                        group: ion_trap,
                        model: Some(MassErrorModel {
                            intercept_ppm: -15.0,
                            ..offset
                        }),
                    },
                ],
            }],
        };
        let corrected = Scorer {
            mass_recalibration: Some(Arc::new(recalibration)),
            ..scorer(&expanded, false)
        };
        let hits = corrected.score(&query);
        let hit = &hits[0];
        assert_eq!(
            expanded[hit.peptide_idx].to_string(),
            "MAGSPEPTS[Phospho]IDEK"
        );
        assert!((hit.delta_mass - 15.0).abs() < 0.2, "{}", hit.delta_mass);
        assert!(hit.aligned_delta_mass.abs() < 0.2);
        assert!((hit.expmass - query.precursors[0].mz * 2.0 + 2.0 * PROTON).abs() < 1e-3);
        assert!((hit.average_ppm - 15.0).abs() < 0.2, "{}", hit.average_ppm);
        assert!((hit.signed_fragment_ppm - 15.0).abs() < 0.2);
        assert!(hit.aligned_average_ppm < 0.2);

        // Annotation reports observed peak m/z, not corrected m/z.
        let fragments = corrected.annotate_candidate(&query, hit);
        let observed = query
            .masses
            .iter()
            .map(|mass| mass + PROTON)
            .collect::<Vec<_>>();
        for mz in &fragments.mz_experimental {
            assert!(observed.iter().any(|o| (o - mz).abs() < 1e-3), "{mz}");
        }

        // The same spectrum recorded by an analyzer without a fragment model
        // keeps its fragment error; only the per-file precursor model applies.
        let mut other = query.clone();
        other.acquisition = AcquisitionGroup {
            analyzer: MassAnalyzer::Tof,
            activation: Activation::Hcd,
        };
        // Widen the fragment window so the uncorrected 15 ppm peaks still match.
        let wide = Scorer {
            fragment_tol: Tolerance::Ppm(-30.0, 30.0),
            mass_recalibration: corrected.mass_recalibration.clone(),
            ..scorer(&expanded, false)
        };
        let hits = wide.score(&other);
        let hit = &hits[0];
        assert!(hit.aligned_delta_mass.abs() < 0.2);
        assert!(
            (hit.aligned_average_ppm - 15.0).abs() < 0.5,
            "{}",
            hit.aligned_average_ppm
        );
    }
}
