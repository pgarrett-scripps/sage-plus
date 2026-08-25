use crate::mass::monoisotopic;

use super::{Tolerance, VALID_AA};

#[test]
fn smoke() {
    for ch in VALID_AA {
        assert!(monoisotopic(ch) > 0.0);
    }
}

#[test]
fn tolerances() {
    assert_eq!(
        Tolerance::Ppm(-10.0, 20.0).bounds(1000.0),
        (999.99, 1000.02)
    );
    assert_eq!(
        Tolerance::Ppm(-10.0, 10.0).bounds(487.0),
        (486.99513, 487.00487)
    );
    assert_eq!(
        Tolerance::Ppm(-50.0, 50.0).bounds(1000.0),
        (999.95, 1000.05)
    );
}
