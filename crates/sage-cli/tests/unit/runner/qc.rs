use super::*;

fn ms1(mz: Vec<f32>, intensity: Vec<f32>, representation: Representation) -> RawSpectrum {
    RawSpectrum {
        ms_level: 1,
        representation,
        mz,
        intensity,
        ..Default::default()
    }
}

#[test]
fn polymer_scan_counts_centroided_ms1_only() {
    // Singly protonated PEG n = 10..=15, listed out of order, plus one
    // unrelated peak of the same total intensity.
    let mut mz = (10..=15)
        .rev()
        .map(|n| (18.010_565 + n as f64 * 44.026_215 + 1.007_276) as f32)
        .collect::<Vec<_>>();
    mz.push(300.0);
    let mut intensity = vec![1.0; 6];
    intensity.push(6.0);
    let spectra = vec![
        ms1(mz.clone(), intensity.clone(), Representation::Centroid),
        ms1(mz, intensity, Representation::Profile),
        RawSpectrum {
            ms_level: 2,
            ..Default::default()
        },
    ];
    let stats = scan_polymers(&spectra).expect("MS1 present");
    assert_eq!(stats.ms1_spectra, 1);
    assert_eq!(stats.skipped_profile_spectra, 1);
    let peg = &stats.polymers[0];
    assert_eq!(peg.name, "peg");
    assert!((peg.tic_pct - 50.0).abs() < 1e-9, "{peg:?}");
    assert!(stats.polymers[1..].iter().all(|share| share.tic_pct == 0.0));
}

#[test]
fn polymer_scan_needs_ms1() {
    let spectra = vec![RawSpectrum {
        ms_level: 2,
        ..Default::default()
    }];
    assert!(scan_polymers(&spectra).is_none());
}
