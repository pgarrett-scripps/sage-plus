use super::*;
use crate::database::PeptideIx;

fn synthetic_mobility_data(count: usize) -> (IndexedDatabase, Vec<Feature>) {
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
            let charge = 2 + (index % 3) as u8;
            let hydro = peptide
                .sequence
                .iter()
                .map(|residue| hydrophobicity(*residue))
                .sum::<f64>();
            Feature {
                peptide_idx: PeptideIx(index as u32),
                label: 1,
                spectrum_q: 0.001,
                charge,
                ims: (0.8 + charge as f64 * 0.06 + hydro / 500.0) as f32,
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
fn ion_mobility_settings_have_safe_defaults() {
    let settings = IonMobilitySettings::default();
    assert!(settings.enabled);
    assert_eq!(settings.features, IonMobilityFeatureSet::Basic);
    assert_eq!(settings.folds, 3);
    assert_eq!(settings.ptm_regularization, 25.0);
    assert_eq!(settings.min_training_psms, 200);
}

#[test]
fn enriched_embedding_is_finite_and_charge_aware() {
    let map = amino_acid_map();
    let peptide = Peptide {
        sequence: b"ACDEFGHIK".to_vec().into(),
        modifications: crate::peptide::CompactModifications::default(),
        monoisotopic: 1000.0,
        ..Peptide::default()
    };
    let charge_two = enriched_embed(&peptide, 2, &map);
    let charge_three = enriched_embed(&peptide, 3, &map);
    assert!(charge_two.iter().all(|value| value.is_finite()));
    assert_ne!(charge_two, charge_three);
    assert_eq!(charge_two.len(), ENRICHED_FEATURES);
}

#[test]
fn zero_mobility_is_not_a_training_observation() {
    let mut feature = Feature::default();
    assert!(!valid_mobility(&feature));
    feature.ims = 1.1;
    assert!(valid_mobility(&feature));
    feature.ims = f32::NAN;
    assert!(!valid_mobility(&feature));
}

#[test]
fn ptm_row_has_global_and_charge_specific_effects() {
    let peptide = Peptide {
        sequence: b"AMPEPTIDEK".to_vec().into(),
        modifications: crate::peptide::CompactModifications::from_dense([
            0.0, 15.994_915, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        ]),
        ..Peptide::default()
    };
    let keys = [(ModificationSpecificity::Residue(b'M'), 15.994_915)];
    assert_eq!(
        MobilityPtmOffsetModel::row(&peptide, 2, &keys, &[2, 3]),
        vec![1.0, 1.0, 0.0]
    );
}

#[test]
fn basic_embedding_tracks_composition_charge_and_mass() {
    let map = amino_acid_map();
    let peptide = Peptide {
        sequence: b"ACDEK".to_vec().into(),
        modifications: crate::peptide::CompactModifications::default(),
        monoisotopic: 1_000.0,
        ..Peptide::default()
    };
    let row = basic_embed(&peptide, 2, &map);

    assert_eq!(row[..VALID_AA.len()].iter().sum::<f64>(), 5.0);
    assert!((row[BASIC_PCT_START..BASIC_N_TERMINAL].iter().sum::<f64>() - 1.0).abs() < 1e-12);
    assert_eq!(row[BASIC_CHARGE], 2.0);
    assert_eq!(row[BASIC_INV_CHARGE], 0.5);
    assert_eq!(row[BASIC_LEN], 5.0);
    assert_eq!(row[BASIC_MASS], 1.0);
    assert_eq!(row[BASIC_MZ], 0.5);
    assert_eq!(row[BASIC_INTERCEPT], 1.0);
}

#[test]
fn mobility_metrics_report_known_r_squared_and_error() {
    let features = [1.0, 2.0, 3.0].map(|ims| Feature {
        ims,
        ..Feature::default()
    });
    let (r2, mae) = prediction_metrics(&features, &[1.0, 2.5, 2.5], &[0, 1, 2]);

    assert!((r2 - 0.75).abs() < 1e-12);
    assert!((mae - (1.0 / 3.0)).abs() < 1e-12);
}

#[test]
fn mobility_prediction_guardrails_leave_features_untouched() {
    let db = IndexedDatabase::default();
    let original = Feature {
        label: 1,
        spectrum_q: 0.001,
        ims: 1.2,
        predicted_ims: 0.4,
        ..Feature::default()
    };

    let mut features = vec![original.clone()];
    let mut settings = IonMobilitySettings {
        enabled: false,
        ..IonMobilitySettings::default()
    };
    assert_eq!(predict(&db, &mut features, &settings), None);
    assert_eq!(features[0].predicted_ims, original.predicted_ims);

    settings.enabled = true;
    settings.folds = 1;
    assert_eq!(predict(&db, &mut features, &settings), None);

    settings.folds = 3;
    settings.min_training_psms = 2;
    assert_eq!(predict(&db, &mut features, &settings), None);
    assert_eq!(features[0].predicted_ims, original.predicted_ims);
}

#[test]
fn additive_mobility_offsets_validate_configuration_first() {
    let db = IndexedDatabase::default();
    let mut predictions = Vec::new();
    let invalid = IonMobilitySettings {
        ptm_regularization: 0.0,
        ..IonMobilitySettings::default()
    };
    assert!(apply_ptm_offsets(&db, &[], &invalid, &mut predictions)
        .unwrap_err()
        .contains("greater than zero"));

    let error =
        apply_ptm_offsets(&db, &[], &IonMobilitySettings::default(), &mut predictions).unwrap_err();
    assert!(error.contains("no variable modifications"));
}

#[test]
fn mobility_hydrophobicity_handles_scale_extremes_and_unknowns() {
    assert_eq!(hydrophobicity(b'I'), 4.5);
    assert_eq!(hydrophobicity(b'R'), -4.5);
    assert_eq!(hydrophobicity(b'X'), 0.0);
}

#[test]
fn basic_mobility_prediction_runs_cross_fitted_end_to_end() {
    let (db, mut features) = synthetic_mobility_data(420);
    let settings = IonMobilitySettings {
        folds: 3,
        min_training_psms: 300,
        ..IonMobilitySettings::default()
    };

    assert_eq!(predict(&db, &mut features, &settings), Some(()));
    assert!(features.iter().all(|feature| {
        feature.predicted_ims.is_finite()
            && (0.0..=2.0).contains(&feature.predicted_ims)
            && feature.delta_ims_model.is_finite()
    }));
}

#[test]
fn physicochemical_mobility_prediction_uses_enriched_features() {
    let (db, mut features) = synthetic_mobility_data(420);
    let settings = IonMobilitySettings {
        features: IonMobilityFeatureSet::Physicochemical,
        folds: 2,
        min_training_psms: 300,
        ..IonMobilitySettings::default()
    };

    assert_eq!(predict(&db, &mut features, &settings), Some(()));
    assert!(features
        .iter()
        .all(|feature| feature.predicted_ims.is_finite()));
}

#[test]
fn cross_fit_rejects_folds_without_enough_training_rows() {
    let (db, features) = synthetic_mobility_data(4);
    let settings = IonMobilitySettings {
        folds: 2,
        min_training_psms: 1,
        ..IonMobilitySettings::default()
    };
    let error = cross_fit::<BASIC_FEATURES>(&db, &features, &settings, basic_embed).unwrap_err();
    assert!(error.contains("training observations"));
}
