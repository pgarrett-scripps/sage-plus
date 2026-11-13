use super::*;
use sage_core::spectrum::Precursor;

#[test]
fn twin_is_the_best_other_same_mass_peptide() {
    let feature = |label: i32, hyperscore: f64| Feature {
        label,
        hyperscore,
        ..Default::default()
    };
    // A target, its same-mass decoy, a same-mass target and an
    // unrelated peptide that scores highest.
    let features = [
        feature(1, 30.0),
        feature(-1, 20.0),
        feature(1, 25.0),
        feature(-1, 40.0),
    ];
    let masses = [1000.0, 1000.005, 999.995, 1100.0];
    assert_eq!(twin_hyperscore(&features, &masses, 0, "any"), Some(25.0));
    assert_eq!(
        twin_hyperscore(&features, &masses, 0, "opposite"),
        Some(20.0)
    );
    assert_eq!(
        twin_hyperscore(&features, &masses, 1, "opposite"),
        Some(30.0)
    );
    assert_eq!(twin_hyperscore(&features, &masses, 3, "any"), None);
    assert_eq!(twin_hyperscore(&features, &masses, 0, "off"), None);
}

#[test]
fn glycan_peaks_are_excluded_and_peptide_peaks_kept() {
    let mut placed = Peptide::try_from(sage_core::enzyme::Digest {
        sequence: "AANGTK".into(),
        ..Default::default()
    })
    .unwrap();
    let bare = placed.monoisotopic;
    placed.monoisotopic += HEXNAC as f32;
    let hex = crate::composition::Monosaccharide::Hex.mass() as f32;
    let glycan = 2.0 * HEXNAC as f32 + 5.0 * hex;
    // Peaks of undetermined charge, stored as singly charged neutral masses.
    let masses = [
        204.086_6 - PROTON,                     // HexNAc oxonium
        300.0,                                  // unrelated
        (bare + HEXNAC as f32) / 2.0,           // Y1 at 2+
        bare + 2.0 * HEXNAC as f32 + 3.0 * hex, // HexNAc(2)Hex(3) Y ion at 1+
    ];
    let query = ProcessedSpectrum {
        level: 2,
        masses: masses.to_vec(),
        charges: vec![1; masses.len()],
        charge_is_known: vec![false; masses.len()],
        intensities: vec![1.0; masses.len()],
        precursors: vec![Precursor {
            mz: (bare + glycan) / 3.0 + PROTON,
            charge: Some(3),
            ..Default::default()
        }],
        ..Default::default()
    };
    let exclusion = glycan_peak_exclusion(Tolerance::Ppm(-10.0, 10.0));
    let mut excluded = vec![false; masses.len()];
    exclusion(&placed, 3, &query, &mut excluded);
    assert_eq!(excluded, vec![true, false, true, true]);
}
