use super::{average_finite, labeled_finite_values, median_finite, normalize_finite};
use sage_core::scoring::Feature;

#[test]
fn report_statistics_ignore_nonfinite_values() {
    assert_eq!(
        median_finite([f32::NAN, 3.0, 1.0, f32::INFINITY, 2.0]),
        Some(2.0)
    );
    assert_eq!(average_finite([f32::NAN, 2.0, 4.0]), Some(3.0));
    assert_eq!(median_finite([f32::NAN]), None);
    assert_eq!(average_finite([f32::INFINITY]), None);
}

#[test]
fn constant_report_values_normalize_without_nan() {
    assert_eq!(normalize_finite(vec![5.0]), vec![0.0]);
    assert_eq!(normalize_finite(vec![5.0, 5.0]), vec![0.0, 0.0]);
    assert_eq!(normalize_finite(vec![2.0, 4.0]), vec![0.0, 1.0]);
}

#[test]
fn report_score_series_filters_invalid_values_and_labels() {
    let features = vec![
        Feature {
            label: 1,
            discriminant_score: 2.0,
            ..Feature::default()
        },
        Feature {
            label: -1,
            discriminant_score: f32::NAN,
            ..Feature::default()
        },
        Feature {
            label: 0,
            discriminant_score: 1.0,
            ..Feature::default()
        },
    ];
    let values = labeled_finite_values(&features, |feature| feature.discriminant_score as f64);
    assert_eq!(values, (vec![2.0], vec![1]));
}
