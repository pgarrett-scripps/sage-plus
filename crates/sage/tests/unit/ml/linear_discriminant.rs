use super::*;
use crate::ml::*;

#[test]
fn linear_discriminant() {
    let a = Matrix::new([1., 2., 3., 4.], 2, 2);
    let eigenvector = [0.4159736, 0.90937671];
    assert!(all_close(
        &a.power_method(&[0.54, 0.34]),
        &eigenvector,
        1E-5
    ));

    #[rustfmt::skip]
        let feats: [[f64; 4]; 8] = [
            [5., 4., 3., 2.],
            [4., 5., 4., 3.],
            [6., 3., 4., 5.],
            [1., 0., 2., 9.],
            [5., 4., 4., 3.],
            [2., 1., 1., 9.5],
            [1., 0., 2., 8.],
            [3., 2., -2., 10.],
        ];

    let lda = LinearDiscriminantAnalysis::train::<_, 4>(
        &feats,
        &[false, false, false, true, false, true, true, true],
        |row| *row,
    )
    .expect("error training LDA");

    let mut scores: Vec<f64> = feats.iter().map(|row| lda.score(row)).collect();
    let norm = norm(&scores);
    scores = scores.into_iter().map(|s| s / norm).collect();

    let expected = [
        0.49706043,
        0.48920177,
        0.48920177,
        -0.07209359,
        0.51204672,
        -0.02849527,
        -0.04924864,
        -0.06055943,
    ];

    assert!(
        all_close(&scores, &expected, 1E-8),
        "{:?} {:?}",
        scores,
        expected
    );
}

/// Deterministic noise in [-1, 1).
fn noise(state: &mut u64) -> f64 {
    *state = state
        .wrapping_mul(6_364_136_223_846_793_005)
        .wrapping_add(1_442_695_040_888_963_407);
    ((*state >> 11) as f64 / (1u64 << 53) as f64) * 2.0 - 1.0
}

/// Targets score higher than decoys on hyperscore and related features, with
/// enough noise that every class has within-class variance.
fn synthetic_psms(targets: usize, decoys: usize) -> Vec<Feature> {
    let mut state = 42;
    (0..targets + decoys)
        .map(|index| {
            let decoy = index >= targets;
            let shift = if decoy { 0.0 } else { 1.0 };
            let mut n = || noise(&mut state);
            Feature {
                label: if decoy { -1 } else { 1 },
                rank: 1,
                charge: 2 + (index % 2) as u8,
                peptide_len: 10 + (index % 7),
                hyperscore: 20.0 + 10.0 * shift + 4.0 * n(),
                delta_next: 3.0 + 2.0 * shift + n(),
                delta_best: 1.0 + n().abs(),
                aligned_delta_mass: (3.0 - 2.0 * shift + 2.0 * n()) as f32,
                expmass: 1000.0,
                calcmass: 1000.0,
                aligned_average_ppm: (4.0 - shift + n()) as f32,
                poisson: -5.0 - 5.0 * shift + n(),
                matched_intensity_pct: (30.0 + 20.0 * shift + 10.0 * n()) as f32,
                matched_peaks: (8.0 + 4.0 * shift + 3.0 * n()).round() as u32,
                longest_b: (3.0 + 2.0 * shift + 2.0 * n()).abs().round() as u32,
                longest_y: (4.0 + 3.0 * shift + 2.0 * n()).abs().round() as u32,
                missed_cleavages: (index % 3) as u8,
                aligned_rt: (0.5 + 0.4 * n()) as f32,
                ..Default::default()
            }
        })
        .collect()
}

fn mean_discriminant(scores: &[Feature], decoy: bool) -> f64 {
    let selected = scores
        .iter()
        .filter(|feature| (feature.label == -1) == decoy)
        .map(|feature| feature.discriminant_score as f64)
        .collect::<Vec<_>>();
    selected.iter().sum::<f64>() / selected.len() as f64
}

#[test]
fn separable_psms_fit_a_model() {
    let mut scores = synthetic_psms(300, 300);
    score_psms(&mut scores, Tolerance::Ppm(-10.0, 10.0)).expect("model fits");
    assert!(mean_discriminant(&scores, false) > mean_discriminant(&scores, true));
    assert!(scores.iter().all(
        |feature| feature.discriminant_score.is_finite() && feature.posterior_error.is_finite()
    ));
}

#[test]
fn too_few_psms_fail_without_touching_scores() {
    let mut scores = synthetic_psms(300, MIN_CLASS_PSMS - 1);
    assert_eq!(
        score_psms(&mut scores, Tolerance::Ppm(-10.0, 10.0)),
        Err(LdaFailure::TooFewPsms {
            targets: 300,
            decoys: MIN_CLASS_PSMS - 1
        })
    );
    assert!(scores
        .iter()
        .all(|feature| feature.discriminant_score == 0.0));

    let mut no_decoys = synthetic_psms(300, 0);
    assert_eq!(
        score_psms(&mut no_decoys, Tolerance::Ppm(-10.0, 10.0)),
        Err(LdaFailure::TooFewPsms {
            targets: 300,
            decoys: 0
        })
    );
    // The generic trainer rejects an empty class too.
    assert_eq!(
        LinearDiscriminantAnalysis::train::<_, 1>(&[[1.0], [2.0]], &[false, false], |row| *row)
            .err(),
        Some(LdaFailure::TooFewPsms {
            targets: 2,
            decoys: 0
        })
    );
}

#[test]
fn non_finite_features_fail() {
    let mut scores = synthetic_psms(300, 300);
    scores[7].hyperscore = f64::NAN;
    assert_eq!(
        score_psms(&mut scores, Tolerance::Ppm(-10.0, 10.0)),
        Err(LdaFailure::NonFinite)
    );

    let mut infinite = synthetic_psms(300, 300);
    infinite[310].aligned_rt = f32::INFINITY;
    assert_eq!(
        score_psms(&mut infinite, Tolerance::Ppm(-10.0, 10.0)),
        Err(LdaFailure::NonFinite)
    );
}

#[test]
fn singular_scatter_fails() {
    // The first feature is constant within each class but differs between
    // them: its within-class scatter is zero, so only the solver's ridge
    // could set its weight.
    let rows = [
        [1.0, 0.3],
        [1.0, -0.2],
        [1.0, 0.5],
        [0.0, 0.1],
        [0.0, -0.4],
        [0.0, 0.2],
    ];
    let decoy = [false, false, false, true, true, true];
    assert_eq!(
        LinearDiscriminantAnalysis::train::<_, 2>(&rows, &decoy, |row| *row).err(),
        Some(LdaFailure::SingularScatter)
    );
}

#[test]
fn indistinguishable_classes_are_degenerate() {
    // Decoys are exact copies of the targets, so the class means match.
    let targets = synthetic_psms(100, 0);
    let mut scores = targets.clone();
    scores.extend(targets.into_iter().map(|feature| Feature {
        label: -1,
        ..feature
    }));
    assert_eq!(
        score_psms(&mut scores, Tolerance::Ppm(-10.0, 10.0)),
        Err(LdaFailure::Degenerate)
    );
}

#[test]
fn fallback_scores_are_finite_and_ordered() {
    let mut scores = synthetic_psms(300, 300);
    scores[0].poisson = f64::NEG_INFINITY;
    scores[1].poisson = f64::NAN;
    score_psms_fallback(&mut scores);
    assert!(scores.iter().all(
        |feature| feature.discriminant_score.is_finite() && feature.posterior_error.is_finite()
    ));
    assert!(mean_discriminant(&scores, false) > mean_discriminant(&scores, true));
    // An underflowed match probability ranks at the top.
    let best = scores
        .iter()
        .map(|feature| feature.discriminant_score)
        .fold(f32::NEG_INFINITY, f32::max);
    assert_eq!(scores[0].discriminant_score, best);
    // The best targets get a PEP below 1 (log10 PEP below 0).
    let top_target = scores
        .iter()
        .filter(|feature| feature.label == 1)
        .map(|feature| feature.posterior_error)
        .fold(f32::INFINITY, f32::min);
    assert!(top_target < 0.0);

    // Without decoys no PEP can be estimated: every PSM gets PEP 1.
    let mut targets_only = synthetic_psms(50, 0);
    score_psms_fallback(&mut targets_only);
    assert!(targets_only
        .iter()
        .all(|feature| feature.posterior_error == 0.0 && feature.discriminant_score > 0.0));
}

#[test]
fn constant_feature_columns_are_ignored() {
    // Pseudo-spectra and runs without an RT model hold several constant
    // columns (rank, delta_best, ion mobility, model deltas). They carry no
    // information and must not stop the fit.
    let mut scores = synthetic_psms(300, 300)
        .into_iter()
        .map(|feature| Feature {
            charge: 2,
            peptide_len: 12,
            delta_best: 0.1,
            missed_cleavages: 0,
            aligned_rt: 0.37,
            ..feature
        })
        .collect::<Vec<_>>();
    score_psms(&mut scores, Tolerance::Ppm(-10.0, 10.0)).expect("model fits");
    assert!(mean_discriminant(&scores, false) > mean_discriminant(&scores, true));

    // A constant column gets no weight, and the other weights match a fit
    // without it.
    let rows = [
        [3.0, 0.1, 0.3],
        [2.5, 0.1, -0.2],
        [3.5, 0.1, 0.4],
        [1.0, 0.1, 0.1],
        [1.4, 0.1, -0.4],
        [0.7, 0.1, 0.2],
    ];
    let decoy = [false, false, false, true, true, true];
    let with_constant =
        LinearDiscriminantAnalysis::train::<_, 3>(&rows, &decoy, |row| *row).expect("fits");
    let without_constant =
        LinearDiscriminantAnalysis::train::<_, 2>(&rows, &decoy, |row| [row[0], row[2]])
            .expect("fits");
    assert_eq!(with_constant.coef[1], 0.0);
    for row in &rows {
        let reduced = [row[0], row[2]];
        assert!((with_constant.score(row) - without_constant.score(&reduced)).abs() < 1e-9);
    }

    // The regularized fit detects the constant column before adding its ridge.
    let regularized =
        LinearDiscriminantAnalysis::train_regularized::<_, 3>(&rows, &decoy, |row| *row, 1e-3)
            .expect("fits");
    assert_eq!(regularized.coef[1], 0.0);
}
