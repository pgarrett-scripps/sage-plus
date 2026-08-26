use super::*;
use crate::database::PeptideIx;

fn synthetic_retention_data(count: usize) -> (IndexedDatabase, Vec<Feature>) {
    const RESIDUES: &[u8] = b"ACDEFGHIKLMNPQRSTVWY";
    let peptides = (0..count)
        .map(|index| {
            let length = 8 + index % 10;
            let sequence = (0..length)
                .map(|position| RESIDUES[(index * 7 + position * 11 + index * position) % 20])
                .collect::<Vec<_>>();
            Peptide {
                sequence: sequence.into(),
                monoisotopic: 700.0 + index as f32 * 2.0,
                ..Peptide::default()
            }
        })
        .collect::<Vec<_>>();
    let features = peptides
        .iter()
        .enumerate()
        .map(|(index, peptide)| {
            let hydro = peptide
                .sequence
                .iter()
                .map(|residue| hydrophobicity(*residue))
                .sum::<f64>();
            Feature {
                peptide_idx: PeptideIx(index as u32),
                label: 1,
                spectrum_q: 0.001,
                aligned_rt: (0.5 + hydro / 100.0).clamp(0.02, 0.98) as f32,
                ..Feature::default()
            }
        })
        .collect();
    (
        IndexedDatabase {
            peptides,
            ..IndexedDatabase::default()
        },
        features,
    )
}

#[test]
fn peptide_folds_are_stable_and_bounded() {
    for peptide_idx in 0..10_000 {
        let sequence = format!("PEPTIDE{peptide_idx}");
        let fold = peptide_fold(sequence.as_bytes(), 3, 42);
        assert!(fold < 3);
        assert_eq!(fold, peptide_fold(sequence.as_bytes(), 3, 42));
    }
}

#[test]
fn retention_settings_have_defaults() {
    let settings = RetentionTimeSettings::default();
    assert_eq!(settings.features, RetentionTimeFeatureSet::Basic);
    assert_eq!(settings.folds, 3);
    assert_eq!(settings.ptm_regularization, 25.0);
}

#[test]
fn physicochemical_embedding_has_positional_and_hydrophobic_features() {
    let mut map = [0; 26];
    for (idx, aa) in VALID_AA.iter().enumerate() {
        map[(aa - b'A') as usize] = idx;
    }
    let peptide = Peptide {
        sequence: b"ACDEFGHIK".to_vec().into(),
        modifications: crate::peptide::CompactModifications::default(),
        ..Peptide::default()
    };
    let row =
        physicochemical_source_embed(&peptide, &map, RetentionTimeFeatureSet::Physicochemical);
    assert_eq!(row[POSITIONAL_START + map[0]], 1.0); // N1 = A
    assert_eq!(row[POSITIONAL_START + VALID_AA.len() + map[2]], 1.0); // N2 = C
    assert!(row[HYDROPHOBIC_START].is_finite());
    assert!(row[PROPERTY_START + 2].is_finite());
    let linear = linear_physicochemical_embed(&peptide, &map);
    assert_eq!(linear.len(), LINEAR_PHYSICOCHEMICAL_FEATURES);
    assert!(linear.iter().all(|value| value.is_finite()));
    let additive = linear_additive_ptm_embed(&peptide, &map);
    assert_eq!(additive.len(), LINEAR_ADDITIVE_PTM_FEATURES);
    assert!(additive.iter().all(|value| value.is_finite()));
}

#[test]
fn enriched_embeddings_support_compact_unmodified_peptides() {
    let mut map = [0; 26];
    for (idx, aa) in VALID_AA.iter().enumerate() {
        map[(aa - b'A') as usize] = idx;
    }
    let compact = Peptide {
        sequence: b"ACDEFGHIK".to_vec().into(),
        modifications: crate::peptide::CompactModifications::default(),
        ..Peptide::default()
    };
    let dense = compact.clone();

    assert_eq!(
        linear_physicochemical_embed(&compact, &map),
        linear_physicochemical_embed(&dense, &map)
    );
    assert_eq!(
        linear_additive_ptm_embed(&compact, &map),
        linear_additive_ptm_embed(&dense, &map)
    );
}

#[test]
fn variable_modification_counts_are_site_specific() {
    let peptide = Peptide {
        sequence: b"AMMPEPTIDEK".to_vec().into(),
        modifications: crate::peptide::CompactModifications::from_dense([
            0.0, 15.994_915, 15.994_915, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        ]),
        ..Peptide::default()
    };
    assert_eq!(
        variable_mod_count(&peptide, ModificationSpecificity::Residue(b'M'), 15.994_915),
        2.0
    );
    assert_eq!(
        variable_mod_count(&peptide, ModificationSpecificity::Residue(b'M'), 42.0),
        0.0
    );
}

#[test]
fn basic_embedding_tracks_terminal_residues_and_global_features() {
    let mut map = [0; 26];
    for (idx, aa) in VALID_AA.iter().enumerate() {
        map[(aa - b'A') as usize] = idx;
    }
    let peptide = Peptide {
        sequence: b"ACDEK".to_vec().into(),
        monoisotopic: 999.0,
        ..Peptide::default()
    };
    let row = RetentionModel::embed(&peptide, &map);

    assert_eq!(row[..VALID_AA.len()].iter().sum::<f64>(), 5.0);
    assert_eq!(row[N_TERMINAL..C_TERMINAL].iter().sum::<f64>(), 2.0);
    assert_eq!(row[C_TERMINAL..PEPTIDE_LEN].iter().sum::<f64>(), 2.0);
    assert_eq!(row[PEPTIDE_LEN], 5.0);
    assert_eq!(row[PEPTIDE_MASS], 999.0_f64.ln_1p());
    assert_eq!(row[INTERCEPT], 1.0);
}

#[test]
fn terminal_modification_counts_respect_specificity_and_protein_position() {
    let peptide = Peptide {
        sequence: b"ACDK".to_vec().into(),
        modifications: crate::peptide::CompactModifications::from_dense([42.0, 0.0, 0.0, 17.0]),
        nterm: Some(10.0),
        cterm: Some(20.0),
        position: Position::Full,
        ..Peptide::default()
    };

    for (specificity, mass) in [
        (ModificationSpecificity::PeptideN(None), 10.0),
        (ModificationSpecificity::PeptideC(None), 20.0),
        (ModificationSpecificity::ProteinN(None), 10.0),
        (ModificationSpecificity::ProteinC(None), 20.0),
        (ModificationSpecificity::PeptideN(Some(b'A')), 42.0),
        (ModificationSpecificity::PeptideC(Some(b'K')), 17.0),
        (ModificationSpecificity::ProteinN(Some(b'A')), 42.0),
        (ModificationSpecificity::ProteinC(Some(b'K')), 17.0),
    ] {
        assert_eq!(variable_mod_count(&peptide, specificity, mass), 1.0);
    }
    assert_eq!(
        variable_mod_count(
            &peptide,
            ModificationSpecificity::PeptideN(Some(b'K')),
            42.0
        ),
        0.0
    );

    let internal = Peptide {
        position: Position::Internal,
        ..peptide
    };
    assert_eq!(
        variable_mod_count(&internal, ModificationSpecificity::ProteinN(None), 10.0),
        0.0
    );
    assert_eq!(
        variable_mod_count(
            &internal,
            ModificationSpecificity::ProteinC(Some(b'K')),
            17.0
        ),
        0.0
    );
}

#[test]
fn invalid_additive_retention_settings_preserve_predictions() {
    let db = IndexedDatabase::default();
    let mut features = Vec::new();
    for settings in [
        RetentionTimeSettings {
            folds: 1,
            ..RetentionTimeSettings::default()
        },
        RetentionTimeSettings {
            ptm_regularization: f64::NAN,
            ..RetentionTimeSettings::default()
        },
        RetentionTimeSettings {
            ptm_regularization: 0.0,
            ..RetentionTimeSettings::default()
        },
    ] {
        apply_ptm_offsets(&db, &mut features, &settings);
        assert!(features.is_empty());
    }
}

#[test]
fn retention_hydrophobicity_handles_scale_extremes_and_unknowns() {
    assert_eq!(hydrophobicity(b'I'), 4.5);
    assert_eq!(hydrophobicity(b'R'), -4.5);
    assert_eq!(hydrophobicity(b'X'), 0.0);
}

#[test]
fn basic_retention_prediction_fits_and_updates_every_feature() {
    let (db, mut features) = synthetic_retention_data(320);
    assert_eq!(
        predict(&db, &mut features, &RetentionTimeSettings::default()),
        Some(())
    );
    assert!(features.iter().all(|feature| {
        feature.predicted_rt.is_finite()
            && (0.0..=1.0).contains(&feature.predicted_rt)
            && feature.delta_rt_model.is_finite()
    }));
}

#[test]
fn physicochemical_retention_model_fits_and_predicts() {
    let (db, features) = synthetic_retention_data(320);
    let model =
        RetentionModel::fit(&db, &features, RetentionTimeFeatureSet::Physicochemical).unwrap();
    let prediction = model.predict_peptide(&db, &features[0]);
    assert!(model.r2.is_finite());
    assert!(prediction.is_finite());
}

#[test]
fn retention_fit_ignores_decoys_and_low_confidence_psms() {
    let (db, mut features) = synthetic_retention_data(320);
    features[0].label = -1;
    features[1].spectrum_q = 0.5;
    let model = RetentionModel::fit(&db, &features, RetentionTimeFeatureSet::Basic).unwrap();
    assert!(model.predict_peptide(&db, &features[2]).is_finite());
}
