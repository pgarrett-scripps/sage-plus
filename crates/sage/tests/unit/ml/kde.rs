use super::*;

#[test]
fn gaussian_kde_is_symmetric_and_peaks_near_the_sample() {
    let sample = [-1.0, 0.0, 1.0];
    let kde = Kde::new(&sample, std::convert::identity);

    assert!(kde.bandwidth.is_finite());
    assert!(kde.bandwidth > 0.0);
    assert!((kde.pdf(-0.5) - kde.pdf(0.5)).abs() < 1e-12);
    assert!(kde.pdf(0.0) > kde.pdf(5.0));
}

#[test]
fn bandwidth_adjustment_is_applied() {
    let sample = [-2.0, -1.0, 1.0, 2.0];
    let baseline = Kde::new(&sample, std::convert::identity);
    let adjusted = Kde::new(&sample, |bandwidth| bandwidth * 1.5);

    assert!((adjusted.bandwidth - baseline.bandwidth * 1.5).abs() < 1e-12);
}

#[test]
fn fitted_posterior_is_finite_bounded_and_improves_with_score() {
    let scores = [-4.0, -3.0, -2.0, 2.0, 3.0, 4.0];
    let decoys = [true, true, true, false, false, false];
    let estimator = Builder::default().bins(101).build(&scores, &decoys);
    let probabilities = (-4..=4)
        .map(|score| estimator.posterior_error(score as f64))
        .collect::<Vec<_>>();

    assert!(probabilities.iter().all(|value| value.is_finite()));
    assert!(probabilities
        .iter()
        .all(|value| (0.0..=1.0).contains(value)));
    assert!(probabilities
        .windows(2)
        .all(|window| window[0] >= window[1]));
    assert!(probabilities[0] > probabilities[probabilities.len() - 1]);
}

#[test]
fn posterior_interpolates_and_clamps_to_the_fitted_range() {
    let estimator = Estimator {
        bins: vec![1.0, 0.5, 0.0],
        min_score: 0.0,
        score_step: 1.0,
    };

    assert_eq!(estimator.posterior_error(0.5), 0.75);
    assert_eq!(estimator.posterior_error(-100.0), 1.0);
    assert_eq!(estimator.posterior_error(100.0), 0.0);
}
