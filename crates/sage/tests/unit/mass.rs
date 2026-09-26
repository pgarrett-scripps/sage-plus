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

fn assert_close(actual: (f32, f32), expected: (f32, f32)) {
    assert!(
        (actual.0 - expected.0).abs() < 1e-4 && (actual.1 - expected.1).abs() < 1e-4,
        "{actual:?} != {expected:?}"
    );
}

#[test]
fn pct_and_da_bounds() {
    assert_close(Tolerance::Pct(-1.0, 2.0).bounds(500.0), (495.0, 510.0));
    assert_close(Tolerance::Da(-0.5, 0.25).bounds(500.0), (499.5, 500.25));
    // Asymmetric windows that do not contain the center are allowed.
    assert_close(Tolerance::Da(1.0, 2.0).bounds(100.0), (101.0, 102.0));
}

#[test]
fn contains_is_inclusive_on_both_ends() {
    let t = Tolerance::Da(-1.0, 1.0);
    assert!(t.contains(100.0, 99.0));
    assert!(t.contains(100.0, 101.0));
    assert!(t.contains(100.0, 100.0));
    assert!(!t.contains(100.0, 98.99));
    assert!(!t.contains(100.0, 101.01));

    let t = Tolerance::Ppm(-10.0, 10.0);
    // 10 ppm of 1000 Da = 0.01 Da
    assert!(t.contains(1000.0, 1000.009));
    assert!(!t.contains(1000.0, 1000.011));
    assert!(!t.contains(1000.0, 999.989));
}

#[test]
fn ppm_to_delta_mass() {
    assert!((Tolerance::ppm_to_delta_mass(1000.0, 10.0) - 0.01).abs() < 1e-6);
    assert!((Tolerance::ppm_to_delta_mass(2500.0, -4.0) + 0.01).abs() < 1e-6);
    assert_eq!(Tolerance::ppm_to_delta_mass(1000.0, 0.0), 0.0);
}

#[test]
fn tolerance_scales_both_bounds_and_keeps_unit() {
    assert_eq!(
        Tolerance::Ppm(-10.0, 20.0) * 2.0,
        Tolerance::Ppm(-20.0, 40.0)
    );
    assert_eq!(Tolerance::Pct(-1.0, 3.0) * 0.5, Tolerance::Pct(-0.5, 1.5));
    match Tolerance::Da(-0.02, 0.04) * 3.0 {
        Tolerance::Da(lo, hi) => assert_close((lo, hi), (-0.06, 0.12)),
        other => panic!("unit changed: {other:?}"),
    }
}

#[test]
fn monoisotopic_rejects_non_residues() {
    // Non-uppercase bytes have no residue mass.
    assert_eq!(monoisotopic(b'a'), 0.0);
    assert_eq!(monoisotopic(b'1'), 0.0);
    assert_eq!(monoisotopic(b'-'), 0.0);
    // Uppercase letters without an amino acid (B, J, X, Z) map to 0.
    for ch in *b"BJXZ" {
        assert_eq!(monoisotopic(ch), 0.0, "{}", ch as char);
    }
    assert!((monoisotopic(b'G') - 57.02146).abs() < 1e-5);
    assert!((monoisotopic(b'W') - 186.07932).abs() < 1e-5);
    // Leucine and isoleucine are isobaric.
    assert_eq!(monoisotopic(b'L'), monoisotopic(b'I'));
}

#[test]
fn composition_sums_carbon_and_sulfur() {
    use super::{composition, Composition};
    // P(5) E(5) P(5) T(4) I(6) D(4) E(5) = 34 carbons, no sulfur
    let c: Composition = b"PEPTIDE".iter().map(|&aa| composition(aa)).sum();
    assert_eq!((c.carbon, c.sulfur), (34, 0));
    // M(5, S1) C(3, S1) M(5, S1)
    let c: Composition = b"MCM".iter().map(|&aa| composition(aa)).sum();
    assert_eq!((c.carbon, c.sulfur), (13, 3));
    // Unknown residues contribute nothing.
    let c: Composition = b"XBa".iter().map(|&aa| composition(aa)).sum();
    assert_eq!((c.carbon, c.sulfur), (0, 0));
    // Every valid residue has at least two carbons (glycine).
    for aa in VALID_AA {
        assert!(composition(aa).carbon >= 2, "{}", aa as char);
    }
    let empty: Composition = std::iter::empty().sum();
    assert_eq!((empty.carbon, empty.sulfur), (0, 0));
}
