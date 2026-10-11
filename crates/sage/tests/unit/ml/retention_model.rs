use super::*;
use crate::database::PeptideIx;
use crate::enzyme::Position;

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

/// Deterministic tryptic-like peptides over the 20 canonical residues, with
/// realistic masses (so the mass column is nearly a linear combination of the
/// residue counts, as in real data).
fn noiseless_peptides(count: usize) -> Vec<Peptide> {
    const RESIDUES: &[u8] = b"ACDEFGHIKLMNPQRSTVWY";
    let mut state = 0x2545_f491_4f6c_dd1du64;
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
            Peptide {
                sequence: sequence.into(),
                monoisotopic,
                ..Peptide::default()
            }
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

fn map() -> [usize; 26] {
    let mut map = [0; 26];
    for (idx, aa) in VALID_AA.iter().enumerate() {
        map[(aa - b'A') as usize] = idx;
    }
    map
}

#[test]
fn basic_retention_embedding_fits_noiseless_linear_data() {
    let map = map();
    let peptides = noiseless_peptides(20_000);
    let weights = noiseless_weights::<BASIC_FEATURES>(7);
    let items = peptides
        .iter()
        .map(|peptide| {
            let row = RetentionModel::embed(peptide, &map);
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
    eprintln!("basic RT noiseless max |error| = {error:e} (max |y| = {scale:.3})");
    assert!(error < 1e-9, "max |error| = {error:e}");
    assert!((lr.r2 - 1.0).abs() < 1e-12, "r2 = {}", lr.r2);
}

#[test]
fn basic_retention_predictions_ignore_redundant_and_permuted_columns() {
    const AUGMENTED: usize = BASIC_FEATURES + 1;
    let map = map();
    let peptides = noiseless_peptides(8_000);
    let weights = noiseless_weights::<BASIC_FEATURES>(11);
    // A realistic, noisy target so the redundant directions are not trivially
    // pinned by an exact solution.
    let items = peptides
        .iter()
        .enumerate()
        .map(|(idx, peptide)| {
            let row = RetentionModel::embed(peptide, &map);
            (
                row,
                dot(&row, &weights) + ((idx as f64) * 0.37).sin() * 0.05,
            )
        })
        .collect::<Vec<_>>();
    let base =
        LinearRegression::fit::<_, BASIC_FEATURES>(&items, |_| true, |(x, _)| *x, |(_, y)| *y)
            .unwrap();
    // Append 2 * length - 0.5 * count(A) + 3 * intercept.
    let augment = |x: &[f64; BASIC_FEATURES]| -> [f64; AUGMENTED] {
        let mut row = [0.0; AUGMENTED];
        row[..BASIC_FEATURES].copy_from_slice(x);
        row[BASIC_FEATURES] = 2.0 * x[PEPTIDE_LEN] - 0.5 * x[0] + 3.0 * x[INTERCEPT];
        row
    };
    let augmented =
        LinearRegression::fit::<_, AUGMENTED>(&items, |_| true, |(x, _)| augment(x), |(_, y)| *y)
            .unwrap();
    let permute = |x: &[f64; BASIC_FEATURES]| -> [f64; BASIC_FEATURES] {
        std::array::from_fn(|j| x[BASIC_FEATURES - 1 - j])
    };
    let permuted = LinearRegression::fit::<_, BASIC_FEATURES>(
        &items,
        |_| true,
        |(x, _)| permute(x),
        |(_, y)| *y,
    )
    .unwrap();
    let mut redundant = 0.0f64;
    let mut reordered = 0.0f64;
    for (x, _) in &items {
        let reference = dot(x, &base.beta);
        redundant = redundant.max((dot(&augment(x), &augmented.beta) - reference).abs());
        reordered = reordered.max((dot(&permute(x), &permuted.beta) - reference).abs());
    }
    eprintln!("basic RT redundant-column shift = {redundant:e}, permutation shift = {reordered:e}");
    assert!(
        redundant < 1e-9,
        "redundant column moved predictions by {redundant:e}"
    );
    assert!(
        reordered < 1e-9,
        "column permutation moved predictions by {reordered:e}"
    );
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

fn offset_data() -> (IndexedDatabase, Vec<Feature>) {
    let oxidation = (ModificationSpecificity::Residue(b'M'), 15.994_915);
    let db = IndexedDatabase {
        peptides: vec![oxidized(b"PEMK"), oxidized(b"MPMK"), oxidized(b"PEPK")],
        // A duplicate key collapses; a key no peptide carries gets no offset.
        model_mods: vec![
            oxidation,
            oxidation,
            (ModificationSpecificity::Residue(b'S'), 79.966_33),
        ],
        ..IndexedDatabase::default()
    };
    let features = [0.625f32, 0.75, 0.5625]
        .into_iter()
        .enumerate()
        .map(|(index, aligned_rt)| Feature {
            peptide_idx: PeptideIx(index as u32),
            label: 1,
            spectrum_q: 0.001,
            aligned_rt,
            predicted_rt: 0.5,
            ..Feature::default()
        })
        .collect();
    (db, features)
}

#[test]
fn ptm_offsets_are_ridge_regressed_residuals_per_modification() {
    let (db, features) = offset_data();
    let model = PtmOffsetModel::fit(&db, &features, &[0, 1, 2], 1.0).unwrap();
    assert_eq!(model.keys.len(), 2);
    // Oxidation counts [1, 2, 0], residuals [0.125, 0.25, 0.0625]:
    // offset = (1 * 0.125 + 2 * 0.25) / (1 + 4 + ridge 1).
    let oxidation = 0.625 / 6.0;
    for (key, offset) in model.keys.iter().zip(&model.offsets) {
        let want = if key.0 == ModificationSpecificity::Residue(b'M') {
            oxidation
        } else {
            0.0
        };
        assert!((offset - want).abs() < 1e-12, "{offset}");
    }
    assert!((model.predict(&db[PeptideIx(1)]) - 2.0 * oxidation).abs() < 1e-12);
    assert_eq!(model.predict(&db[PeptideIx(2)]), 0.0);

    // Without modification keys there is nothing to fit.
    let bare = IndexedDatabase {
        peptides: db.peptides.clone(),
        ..IndexedDatabase::default()
    };
    assert!(PtmOffsetModel::fit(&bare, &features, &[0, 1, 2], 1.0).is_none());
}

/// A seed that puts the two oxidized peptides in different folds, so each is
/// corrected by a model trained on the other. Two folds never split them: the
/// hash parity depends only on the residues.
fn split_seed(folds: usize) -> u64 {
    (0..1000)
        .find(|&seed| peptide_fold(b"PEMK", folds, seed) != peptide_fold(b"MPMK", folds, seed))
        .unwrap()
}

#[test]
fn ptm_offsets_leave_unmodified_peptides_and_refresh_deltas() {
    let settings = RetentionTimeSettings {
        folds: 3,
        seed: split_seed(3),
        ptm_regularization: 1.0,
        ..RetentionTimeSettings::default()
    };
    let (db, mut features) = offset_data();
    let mut without_decoy = features.clone();
    apply_ptm_offsets(&db, &mut without_decoy, &settings);
    // A decoy never trains the offsets: a late decoy changes nothing.
    features.push(Feature {
        peptide_idx: PeptideIx(1),
        label: -1,
        aligned_rt: 0.9,
        predicted_rt: 0.5,
        ..Feature::default()
    });
    apply_ptm_offsets(&db, &mut features, &settings);
    for (with, without) in features.iter().zip(&without_decoy) {
        assert_eq!(with.predicted_rt, without.predicted_rt);
    }
    // Oxidation residuals are all positive, so both oxidized targets move up.
    assert!(
        features[0].predicted_rt > 0.5,
        "{}",
        features[0].predicted_rt
    );
    assert!(features[1].predicted_rt > features[0].predicted_rt);
    // The unmodified peptide keeps its prediction.
    assert_eq!(features[2].predicted_rt, 0.5);
    assert_eq!(features[2].delta_rt_model, 0.0625);
    for feature in &features {
        assert!((0.0..=1.0).contains(&feature.predicted_rt));
        assert_eq!(
            feature.delta_rt_model,
            (feature.aligned_rt - feature.predicted_rt).abs()
        );
    }
}
