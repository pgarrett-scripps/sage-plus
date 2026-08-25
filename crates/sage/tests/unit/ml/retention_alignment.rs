use super::*;

fn transform_normalized(alignment: &Alignment, x: f64) -> f64 {
    alignment.transform((x * alignment.max_rt as f64) as f32) as f64
}

#[test]
fn nonlinear_warp_beats_affine_fit_with_outliers() {
    let points: Vec<_> = (0..160)
        .map(|i| {
            let x = i as f64 / 159.0;
            let expected = 0.04 + 0.9 * x + 0.055 * (std::f64::consts::PI * x).sin();
            let y = if i % 29 == 0 {
                1.0 - expected
            } else {
                expected
            };
            (x, y)
        })
        .collect();
    let alignment = fit_alignment(0, 100.0, points);
    assert!(alignment.knots.len() >= 4);

    let (nonlinear_error, affine_error) = (0..100).fold((0.0, 0.0), |(nl, linear), i| {
        let x = i as f64 / 99.0;
        let expected = 0.04 + 0.9 * x + 0.055 * (std::f64::consts::PI * x).sin();
        (
            nl + (transform_normalized(&alignment, x) - expected).abs(),
            linear + (alignment.slope as f64 * x + alignment.intercept as f64 - expected).abs(),
        )
    });
    assert!(nonlinear_error < affine_error * 0.45);
}

#[test]
fn fitted_warp_is_monotone() {
    let points: Vec<_> = (0..80)
        .map(|i| {
            let x = i as f64 / 79.0;
            let noise = if i % 7 == 0 { -0.04 } else { 0.0 };
            (x, x.powf(1.35) + noise)
        })
        .collect();
    let alignment = fit_alignment(0, 60.0, points);
    let transformed: Vec<_> = (0..200)
        .map(|i| transform_normalized(&alignment, i as f64 / 199.0))
        .collect();
    assert!(transformed.windows(2).all(|pair| pair[0] <= pair[1]));
}

#[test]
fn sparse_landmarks_use_affine_fallback() {
    let points = vec![(0.1, 0.25), (0.4, 0.46), (0.8, 0.74), (0.95, 0.845)];
    let alignment = fit_alignment(2, 80.0, points);
    assert!(alignment.knots.is_empty());
    assert!((transform_normalized(&alignment, 0.6) - 0.6).abs() < 0.01);
}

#[test]
fn robust_fit_handles_large_shift_and_bad_landmarks() {
    let points: Vec<_> = (0..100)
        .map(|i| {
            let x = i as f64 / 99.0;
            let expected = 0.2 + 0.65 * x;
            (
                x,
                if i % 17 == 0 {
                    1.0 - expected
                } else {
                    expected
                },
            )
        })
        .collect();
    let alignment = fit_alignment(0, 120.0, points);
    assert!((transform_normalized(&alignment, 0.5) - 0.525).abs() < 0.01);
}

#[test]
fn alignment_method_is_explicit_and_defaults_to_linear() {
    assert_eq!(AlignmentMethod::default(), AlignmentMethod::Linear);
}

#[test]
fn reference_alignment_recovers_affine_shift_with_outliers() {
    let points = (0..80)
        .map(|i| {
            let observed = i as f32 + 1.0;
            let reference = if i % 19 == 0 {
                200.0 - observed
            } else {
                observed * 0.85 + 3.0
            };
            (observed, reference)
        })
        .collect::<Vec<_>>();
    let alignment = fit_reference_alignment(&points, 16).unwrap();
    assert!((alignment.slope - 0.85).abs() < 0.01);
    assert!((alignment.intercept - 3.0).abs() < 0.1);
    assert!((alignment.transform(40.0) - 37.0).abs() < 0.1);
    assert!(alignment.inliers < alignment.points);
}
