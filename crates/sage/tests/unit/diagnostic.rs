use super::*;

#[test]
fn finds_ions_within_tolerance() {
    let ions = default_ions();
    // HexNAc 2 ppm off, Hex 50 ppm off, two peaks near NeuAc.
    let mz = [138.2, 163.0683, 204.0871, 292.1024, 292.103, 500.0];
    let intensity = [10.0, 10.0, 20.0, 5.0, 15.0, 40.0];
    let hits = find_ions(&ions, &mz, &intensity);
    let names = hits
        .iter()
        .map(|hit| ions[hit.ion].name.as_str())
        .collect::<Vec<_>>();
    assert_eq!(names, ["HexNAc", "NeuAc"]);
    assert_eq!(hits[0].mz, 204.0871);
    assert!((hits[0].relative_intensity - 0.2).abs() < 1e-6);
    // The more intense of the two NeuAc peaks is reported.
    assert_eq!(hits[1].mz, 292.103);
}

#[test]
fn user_tolerance_overrides_default() {
    let ions = vec![DiagnosticIon {
        name: "wide".into(),
        mz: 163.0601,
        tolerance: Some(Tolerance::Da(-0.01, 0.01)),
    }];
    let hits = find_ions(&ions, &[163.0683], &[1.0]);
    assert_eq!(hits.len(), 1);
    assert_eq!(hits[0].relative_intensity, 1.0);
}

#[test]
fn empty_spectrum_has_no_hits() {
    assert!(find_ions(&default_ions(), &[], &[]).is_empty());
    assert!(find_ions(&default_ions(), &[204.0867], &[0.0]).is_empty());
}
