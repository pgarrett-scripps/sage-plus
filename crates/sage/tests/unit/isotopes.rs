use super::{averagine_isotopes, peptide_isotopes};

#[test]
fn smoke_isotopes() {
    let iso = peptide_isotopes(60, 5);
    let mut expected = [0.3972, 0.2824, 0.1869, 0.0846];
    expected.iter_mut().for_each(|val| *val /= 0.3972);

    let matched = iso.iter().zip(expected).all(|(a, b)| (a - b).abs() <= 0.02);

    assert!(matched, "{:?} {:?}", iso, expected);
}

#[test]
fn averagine_pattern_is_normalized_and_mass_dependent() {
    let low = averagine_isotopes(500.0);
    let high = averagine_isotopes(2500.0);

    assert!((low.iter().sum::<f32>() - 1.0).abs() < 1e-6);
    assert!((high.iter().sum::<f32>() - 1.0).abs() < 1e-6);
    assert!(low[0] > low[1]);
    assert!(high[1] > high[0]);
}
