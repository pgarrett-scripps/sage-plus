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
fn residue_class_features_count_the_intended_residues() {
    let map = amino_acid_map();
    // 4 bulky (LVIW), 2 polar (ST), 3 positive (RKH), 1 negative (D),
    // 3 tiny (GAS), 3 branched (LVI); C, M and P belong to no class
    let peptide = Peptide {
        sequence: b"LVWIRKHDGACMSTP".to_vec().into(),
        modifications: crate::peptide::CompactModifications::default(),
        monoisotopic: 1_000.0,
        ..Peptide::default()
    };
    let expected = [4.0, 2.0, 3.0, 1.0, 3.0, 3.0];

    let row = basic_embed(&peptide, 2, &map);
    assert_eq!(
        [
            row[BASIC_NUM_BULKY],
            row[BASIC_NUM_UC_POLAR],
            row[BASIC_NUM_POSITIVE],
            row[BASIC_NUM_NEGATIVE],
            row[BASIC_NUM_TINY],
            row[BASIC_NUM_BRANCHED],
        ],
        expected
    );

    let row = enriched_embed(&peptide, 2, &map);
    let start = ENRICHED_FEATURES - ENRICHED_GLOBAL_FEATURES - ENRICHED_PROPERTY_FEATURES;
    let length = peptide.sequence.len() as f64;
    assert_eq!(
        row[start..start + ENRICHED_PROPERTY_FEATURES],
        expected.map(|count| count / length)
    );
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

fn noiseless_peptides(count: usize) -> Vec<(Peptide, u8)> {
    const RESIDUES: &[u8] = b"ACDEFGHIKLMNPQRSTVWY";
    let mut state = 0x9e37_79b9_7f4a_7c15u64;
    let mut next = move || {
        state ^= state << 13;
        state ^= state >> 7;
        state ^= state << 17;
        state
    };
    (0..count)
        .map(|_| {
            let length = 7 + (next() % 24) as usize;
            let mut sequence = (0..length - 1)
                .map(|_| RESIDUES[(next() % 20) as usize])
                .collect::<Vec<_>>();
            sequence.push(if next() % 2 == 0 { b'K' } else { b'R' });
            let oxidized = next() % 5 == 0;
            let monoisotopic = sequence
                .iter()
                .map(|&residue| crate::mass::monoisotopic(residue))
                .sum::<f32>()
                + crate::mass::H2O
                + if oxidized { 15.9949 } else { 0.0 };
            let charge = 2 + (next() % 3) as u8;
            (
                Peptide {
                    sequence: sequence.into(),
                    monoisotopic,
                    ..Peptide::default()
                },
                charge,
            )
        })
        .collect()
}

fn noiseless_weights<const D: usize>(seed: u64) -> [f64; D] {
    let mut state = seed;
    std::array::from_fn(|_| {
        state = state
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1_442_695_040_888_963_407);
        ((state >> 11) as f64 / (1u64 << 53) as f64 - 0.5) * 0.2
    })
}

fn dot<const D: usize>(row: &[f64; D], beta: &[f64]) -> f64 {
    row.iter().zip(beta).map(|(x, w)| x * w).sum()
}

#[test]
fn basic_mobility_embedding_fits_noiseless_linear_data() {
    let map = amino_acid_map();
    let weights = noiseless_weights::<BASIC_FEATURES>(3);
    let items = noiseless_peptides(20_000)
        .iter()
        .map(|(peptide, charge)| {
            let row = basic_embed(peptide, *charge, &map);
            (row, dot(&row, &weights))
        })
        .collect::<Vec<_>>();
    let scale = items.iter().map(|(_, y)| y.abs()).fold(0.0, f64::max);
    let lr = LinearRegression::fit::<_, BASIC_FEATURES>(&items, |_| true, |(x, _)| *x, |(_, y)| *y)
        .unwrap();
    let error = items
        .iter()
        .map(|(x, y)| (dot(x, &lr.beta) - y).abs())
        .fold(0.0, f64::max);
    eprintln!("basic mobility noiseless max |error| = {error:e} (max |y| = {scale:.3})");
    assert!(error < 1e-9, "max |error| = {error:e}");
}

#[test]
fn enriched_mobility_embedding_fits_noiseless_linear_data() {
    let map = amino_acid_map();
    let weights = noiseless_weights::<ENRICHED_FEATURES>(5);
    let items = noiseless_peptides(20_000)
        .iter()
        .map(|(peptide, charge)| {
            let row = enriched_embed(peptide, *charge, &map);
            (row, dot(&row, &weights))
        })
        .collect::<Vec<_>>();
    let lr =
        LinearRegression::fit::<_, ENRICHED_FEATURES>(&items, |_| true, |(x, _)| *x, |(_, y)| *y)
            .unwrap();
    let error = items
        .iter()
        .map(|(x, y)| (dot(x, &lr.beta) - y).abs())
        .fold(0.0, f64::max);
    eprintln!("enriched mobility noiseless max |error| = {error:e}");
    assert!(error < 1e-9, "max |error| = {error:e}");
}

#[test]
fn basic_mobility_predictions_ignore_redundant_class_columns() {
    let map = amino_acid_map();
    let weights = noiseless_weights::<BASIC_FEATURES>(9);
    let items = noiseless_peptides(8_000)
        .iter()
        .enumerate()
        .map(|(idx, (peptide, charge))| {
            let row = basic_embed(peptide, *charge, &map);
            (
                row,
                dot(&row, &weights) + ((idx as f64) * 0.37).sin() * 0.01,
            )
        })
        .collect::<Vec<_>>();
    let base =
        LinearRegression::fit::<_, BASIC_FEATURES>(&items, |_| true, |(x, _)| *x, |(_, y)| *y)
            .unwrap();
    // Residue-class counts are sums of residue counts. Replacing them with other
    // linear combinations of the residue counts must not change the fit.
    let relabel = |x: &[f64; BASIC_FEATURES]| -> [f64; BASIC_FEATURES] {
        let mut row = *x;
        for (offset, column) in (BASIC_NUM_BRANCHED..=BASIC_NUM_NEGATIVE).enumerate() {
            row[column] = x[offset] + 2.0 * x[offset + 3] - x[BASIC_LEN] + 0.5 * x[BASIC_INTERCEPT];
        }
        row
    };
    let relabeled = LinearRegression::fit::<_, BASIC_FEATURES>(
        &items,
        |_| true,
        |(x, _)| relabel(x),
        |(_, y)| *y,
    )
    .unwrap();
    let shift = items
        .iter()
        .map(|(x, _)| (dot(&relabel(x), &relabeled.beta) - dot(x, &base.beta)).abs())
        .fold(0.0, f64::max);
    eprintln!("basic mobility redundant-class-column shift = {shift:e}");
    assert!(
        shift < 1e-9,
        "redundant class columns moved predictions by {shift:e}"
    );
}

#[test]
fn mobility_fit_is_insensitive_to_the_standardized_ridge_penalty() {
    let map = amino_acid_map();
    let weights = noiseless_weights::<BASIC_FEATURES>(3);
    let clean = noiseless_peptides(50_000)
        .iter()
        .map(|(peptide, charge)| {
            let row = basic_embed(peptide, *charge, &map);
            (row, dot(&row, &weights))
        })
        .collect::<Vec<_>>();
    let noisy = clean
        .iter()
        .enumerate()
        .map(|(idx, (x, y))| (*x, y + ((idx as f64) * 0.37).sin() * 0.02))
        .collect::<Vec<_>>();
    let fit = |items: &[([f64; BASIC_FEATURES], f64)], ridge: f64| {
        LinearRegression::fit_with_ridge::<_, BASIC_FEATURES>(
            items,
            |_| true,
            |(x, _)| *x,
            |(_, y)| *y,
            ridge,
        )
        .unwrap()
    };
    let reference = fit(&noisy, crate::ml::regression::RIDGE_PER_ROW);
    for (ridge, tolerance) in [
        (1e-12, 1e-9),
        (1e-10, 1e-9),
        (1e-8, 1e-9),
        (1e-6, 1e-8),
        (1e-5, 1e-5),
        (1e-4, 1e-3),
        (1e-3, 1e-2),
    ] {
        let exact = fit(&clean, ridge);
        let error = clean
            .iter()
            .map(|(x, y)| (dot(x, &exact.beta) - y).abs())
            .fold(0.0, f64::max);
        let perturbed = fit(&noisy, ridge);
        let shift = noisy
            .iter()
            .map(|(x, _)| (dot(x, &perturbed.beta) - dot(x, &reference.beta)).abs())
            .fold(0.0, f64::max);
        eprintln!("ridge {ridge:e}: noiseless max |error| = {error:e}, noisy-fit shift vs default = {shift:e}");
        assert!(error < tolerance, "ridge {ridge:e}: error {error:e}");
        assert!(shift < tolerance, "ridge {ridge:e}: shift {shift:e}");
    }
}

fn oxidized(sequence: &[u8]) -> Peptide {
    let modifications = sequence
        .iter()
        .map(|&residue| if residue == b'M' { 15.994_915 } else { 0.0 })
        .collect::<Vec<f32>>();
    Peptide {
        sequence: sequence.to_vec().into(),
        modifications: crate::peptide::CompactModifications::from_dense(modifications),
        ..Peptide::default()
    }
}

/// Charge 2 PSMs with 1, 2 and 0 oxidations, all predicted at 0.5.
fn offset_data() -> (IndexedDatabase, Vec<Feature>, Vec<f64>) {
    let oxidation = (ModificationSpecificity::Residue(b'M'), 15.994_915);
    let db = IndexedDatabase {
        peptides: vec![oxidized(b"PEMK"), oxidized(b"MPMK"), oxidized(b"PEPK")],
        model_mods: vec![oxidation, oxidation],
        ..IndexedDatabase::default()
    };
    let features = [0.625f32, 0.75, 0.5625]
        .into_iter()
        .enumerate()
        .map(|(index, ims)| Feature {
            peptide_idx: PeptideIx(index as u32),
            label: 1,
            spectrum_q: 0.001,
            charge: 2,
            ims,
            ..Feature::default()
        })
        .collect();
    (db, features, vec![0.5; 3])
}

#[test]
fn mobility_ptm_offsets_split_global_and_charge_effects() {
    let (db, features, predictions) = offset_data();
    let model = MobilityPtmOffsetModel::fit(&db, &features, &predictions, &[0, 1, 2], 1.0).unwrap();
    assert_eq!(model.keys.len(), 1);
    assert_eq!(model.charges, vec![2]);
    // Rows are [count, count]; S = sum(count^2) = 5, R = sum(count * residual)
    // = 0.625. Ridge 1 on the global and 4 on the charge column:
    // [[6, 5], [5, 9]] x = [R, R] gives x = [4R, R] / 29.
    let r = 0.625;
    assert!((model.offsets[0] - 4.0 * r / 29.0).abs() < 1e-12);
    assert!((model.offsets[1] - r / 29.0).abs() < 1e-12);
    assert!((model.predict(&db[PeptideIx(1)], 2) - 10.0 * r / 29.0).abs() < 1e-12);
    // An unseen charge only gets the global effect.
    assert!((model.predict(&db[PeptideIx(1)], 3) - 8.0 * r / 29.0).abs() < 1e-12);
    assert_eq!(model.predict(&db[PeptideIx(2)], 2), 0.0);
}

#[test]
fn mobility_ptm_offsets_skip_invalid_rows_and_unmodified_peptides() {
    // A seed that puts the two oxidized peptides in different folds, so each
    // is corrected by a model trained on the other. Two folds never split
    // them: the hash parity depends only on the residues.
    let seed = (0..1000)
        .find(|&seed| peptide_fold(b"PEMK", 3, seed) != peptide_fold(b"MPMK", 3, seed))
        .unwrap();
    let settings = IonMobilitySettings {
        folds: 3,
        seed,
        ptm_regularization: 1.0,
        ..IonMobilitySettings::default()
    };
    let (db, mut features, mut predictions) = offset_data();
    let mut valid_only = predictions.clone();
    apply_ptm_offsets(&db, &features, &settings, &mut valid_only).unwrap();
    // Missing mobility: never trained on and never corrected.
    features.push(Feature {
        peptide_idx: PeptideIx(1),
        label: 1,
        spectrum_q: 0.001,
        charge: 2,
        ims: 0.0,
        ..Feature::default()
    });
    predictions.push(0.5);
    // A non-finite prediction stays non-finite.
    features.push(features[3].clone());
    predictions.push(f64::NAN);
    let (r2, _) = apply_ptm_offsets(&db, &features, &settings, &mut predictions).unwrap();
    assert!(r2.is_finite());
    // The invalid rows do not change the valid rows' corrections.
    assert_eq!(&predictions[..3], &valid_only[..]);
    assert_eq!(predictions[2], 0.5);
    assert_eq!(predictions[3], 0.5);
    assert!(predictions[4].is_nan());
    // Oxidation residuals are positive, so both oxidized peptides move up.
    assert!(predictions[0] > 0.5, "{}", predictions[0]);
    assert!(predictions[1] > predictions[0]);
}
