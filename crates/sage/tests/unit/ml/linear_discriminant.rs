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

#[test]
fn library_rescoring_separates_synthetic_targets_and_decoys() {
    let mut features = (0..80)
        .map(|index| {
            let decoy = index >= 40;
            let variation = (index % 7) as f32 / 100.0;
            Feature {
                label: if decoy { -1 } else { 1 },
                spectral_angle: if decoy {
                    0.25 + variation
                } else {
                    0.80 + variation
                },
                delta_next: if decoy { 0.02 } else { 0.25 } + f64::from(variation),
                explained_library_intensity: if decoy { 0.30 } else { 0.85 },
                explained_query_intensity: if decoy { 0.25 } else { 0.80 },
                matched_peaks: if decoy { 5 } else { 14 } + index % 3,
                aligned_delta_mass: if decoy { 8.0 } else { 0.5 } + variation,
                aligned_average_ppm: if decoy { 7.0 } else { 0.7 } + variation,
                isotope_error: if decoy { 1.003 } else { 0.0 },
                predicted_rt: 0.5,
                delta_rt_model: if decoy { 0.35 } else { 0.02 } + variation,
                ims: 1.0,
                predicted_ims: 1.0,
                delta_ims_model: if decoy { 0.25 } else { 0.01 } + variation,
                ..Feature::default()
            }
        })
        .collect::<Vec<_>>();

    score_library_psms(&mut features).expect("library model should fit");
    let target_mean = features[..40]
        .iter()
        .map(|feature| feature.discriminant_score)
        .sum::<f32>()
        / 40.0;
    let decoy_mean = features[40..]
        .iter()
        .map(|feature| feature.discriminant_score)
        .sum::<f32>()
        / 40.0;
    assert!(target_mean > decoy_mean);
    assert!(features
        .iter()
        .all(|feature| feature.posterior_error.is_finite()));
}
