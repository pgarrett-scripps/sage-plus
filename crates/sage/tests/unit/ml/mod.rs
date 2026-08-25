use super::*;

#[test]
fn vector_statistics_match_known_values() {
    let values = [1.0, 2.0, 3.0];

    assert_eq!(mean(&values), 2.0);
    assert!((std(&values) - (2.0_f64 / 3.0).sqrt()).abs() < 1e-12);
    assert_eq!(norm(&[3.0, 4.0]), 5.0);
}

#[test]
fn all_close_requires_equal_lengths_and_respects_tolerance() {
    assert!(all_close(&[1.0, 2.0], &[1.0001, 2.0001], 0.001));
    assert!(!all_close(&[1.0, 2.0], &[1.0], 0.001));
    assert!(!all_close(&[1.0], &[1.01], 0.001));
}
