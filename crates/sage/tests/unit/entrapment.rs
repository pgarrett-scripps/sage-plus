use super::*;

#[test]
fn the_combined_estimate_is_the_papers_equation_one() {
    // 100 entrapment discoveries beside 900 sample ones, with an entrapment
    // database five times the original: 100 * 1.2 / 1000.
    let estimate = EntrapmentEstimate::combined(900, 100, Some(5.0));
    assert_eq!(estimate.target_discoveries, 900);
    assert_eq!(estimate.entrapment_discoveries, 100);
    let fdp = estimate.estimated_fdp.expect("an estimate");
    assert!((fdp - 0.12).abs() < 1e-6, "{fdp}");

    // The 1/r term is what separates the estimate from the lower bound, and
    // it shrinks as the entrapment database grows.
    let lower_bound = 100.0 / 1000.0;
    assert!(fdp > lower_bound);
    let wide = combined_fdp(900, 100, Some(100.0)).unwrap();
    assert!(wide > lower_bound && wide < fdp);

    // A proportion cannot exceed one, however small the entrapment database.
    assert_eq!(combined_fdp(0, 10, Some(0.01)), Some(1.0));
}

#[test]
fn an_estimate_needs_a_ratio_and_some_discoveries() {
    assert_eq!(combined_fdp(100, 10, None), None);
    assert_eq!(combined_fdp(100, 10, Some(0.0)), None);
    assert_eq!(combined_fdp(100, 10, Some(-1.0)), None);
    assert_eq!(combined_fdp(100, 10, Some(f32::NAN)), None);
    assert_eq!(combined_fdp(0, 0, Some(5.0)), None);
    // No entrapment discovery at all is an estimate of zero, not a missing
    // one: the search reported nothing it should not have.
    assert_eq!(combined_fdp(100, 0, Some(5.0)), Some(0.0));
}

#[test]
fn the_database_ratio_is_entrapment_over_original() {
    assert_eq!(database_ratio(1000, 5000), Some(5.0));
    assert_eq!(database_ratio(0, 5000), None);
    assert_eq!(database_ratio(1000, 0), Some(0.0));
}
