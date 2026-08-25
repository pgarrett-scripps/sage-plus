use super::peptide_isotopes;

#[test]
fn smoke_isotopes() {
    let iso = peptide_isotopes(60, 5);
    let mut expected = [0.3972, 0.2824, 0.1869, 0.0846];
    expected.iter_mut().for_each(|val| *val /= 0.3972);

    let matched = iso.iter().zip(expected).all(|(a, b)| (a - b).abs() <= 0.02);

    assert!(matched, "{:?} {:?}", iso, expected);
}
