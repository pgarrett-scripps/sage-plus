use super::*;

#[test]
fn test_deisotope() {
    let mz = [
        800.9,
        800.9 + NEUTRON * 1.0,
        800.9 + NEUTRON * 2.0,
        803.4080,
        804.4108,
        805.4106,
        806.4116,
        810.0,
        812.0,
        812.0 + NEUTRON / 2.0,
    ];
    let int = [2., 1.5, 1., 4., 3., 2., 1., 1., 9.0, 4.5];
    let settings = DeisotopeSettings {
        min_score: 0.0,
        max_isotope_log2_ratio: 10.0,
        ..DeisotopeSettings::default()
    };
    let peaks = deisotope(&mz, &int, None, 2, settings, 800.91);

    assert_eq!(
        peaks,
        vec![
            Deisotoped {
                mz: 800.9,
                intensity: 2.0,
                charge: None,
                envelope: None,
            },
            Deisotoped {
                mz: 800.9 + NEUTRON * 1.0,
                intensity: 2.5,
                charge: Some(1),
                envelope: None,
            },
            Deisotoped {
                mz: 800.9 + NEUTRON * 2.0,
                intensity: 1.0,
                charge: Some(1),
                envelope: Some(1),
            },
            Deisotoped {
                mz: 803.4080,
                intensity: 10.0,
                charge: Some(1),
                envelope: None,
            },
            Deisotoped {
                mz: 804.4108,
                intensity: 3.0,
                charge: Some(1),
                envelope: Some(3),
            },
            Deisotoped {
                mz: 805.4106,
                intensity: 2.0,
                charge: Some(1),
                envelope: Some(3),
            },
            Deisotoped {
                mz: 806.4116,
                intensity: 1.0,
                charge: Some(1),
                envelope: Some(3),
            },
            Deisotoped {
                mz: 810.0,
                intensity: 1.0,
                charge: None,
                envelope: None,
            },
            Deisotoped {
                mz: 812.0,
                intensity: 13.5,
                charge: Some(2),
                envelope: None,
            },
            Deisotoped {
                mz: 812.0 + NEUTRON / 2.0,
                intensity: 4.5,
                charge: Some(2),
                envelope: Some(8),
            }
        ]
    );
}

#[test]
fn averagine_deisotope_accepts_a_rising_envelope() {
    let charge = 2;
    let neutral_mass = 2500.0;
    let mono_mz = neutral_mass / charge as f32 + PROTON;
    let pattern = crate::isotopes::averagine_isotopes(neutral_mass);
    let mz = [
        mono_mz + 2.0 * NEUTRON / charge as f32,
        mono_mz,
        mono_mz + NEUTRON / charge as f32,
    ];
    let intensity = [
        pattern[2] * 1000.0,
        pattern[0] * 1000.0,
        pattern[1] * 1000.0,
    ];

    let peaks = deisotope(
        &mz,
        &intensity,
        None,
        charge,
        DeisotopeSettings::default(),
        0.0,
    );
    let roots = peaks
        .iter()
        .filter(|peak| peak.envelope.is_none())
        .collect::<Vec<_>>();

    assert_eq!(roots.len(), 1);
    assert_eq!(roots[0].charge, Some(charge));
    assert!((roots[0].mz - mono_mz).abs() < 1e-4);
    assert!((roots[0].intensity - intensity.iter().sum::<f32>()).abs() < 1e-3);
}

#[test]
fn averagine_deisotope_preserves_intensity_with_competing_candidates() {
    let settings = DeisotopeSettings {
        min_score: 0.0,
        max_isotope_log2_ratio: 10.0,
        ..DeisotopeSettings::default()
    };
    let mz = [500.0, 500.0004, 500.0 + NEUTRON];
    let intensity = [100.0, 90.0, 50.0];

    let peaks = deisotope(&mz, &intensity, None, 1, settings, 0.0);
    let retained_intensity = peaks
        .iter()
        .filter(|peak| peak.envelope.is_none())
        .map(|peak| peak.intensity)
        .sum::<f32>();

    assert!((retained_intensity - intensity.iter().sum::<f32>()).abs() < 1e-5);
}

#[test]
fn averagine_deisotope_rejects_an_implausible_intensity_ratio() {
    let mono_mz = 500.0;
    let mz = [mono_mz, mono_mz + NEUTRON];
    let intensity = [1000.0, 0.01];

    let peaks = deisotope(&mz, &intensity, None, 1, DeisotopeSettings::default(), 0.0);

    assert!(peaks.iter().all(|peak| peak.charge.is_none()));
    assert!(peaks.iter().all(|peak| peak.envelope.is_none()));
}

#[test]
fn fragment_charge_array_constrains_deisotope_candidates() {
    let charge = 2;
    let neutral_mass = 1600.0;
    let mono_mz = neutral_mass / charge as f32 + PROTON;
    let pattern = crate::isotopes::averagine_isotopes(neutral_mass);
    let mz = [mono_mz, mono_mz + NEUTRON / charge as f32];
    let intensity = [pattern[0] * 1000.0, pattern[1] * 1000.0];

    let accepted = deisotope(
        &mz,
        &intensity,
        Some(&[charge, charge]),
        3,
        DeisotopeSettings::default(),
        0.0,
    );
    assert_eq!(accepted[0].charge, Some(charge));
    assert_eq!(accepted[1].envelope, Some(0));

    let conflicting = deisotope(
        &mz,
        &intensity,
        Some(&[1, 1]),
        3,
        DeisotopeSettings::default(),
        0.0,
    );
    assert!(conflicting.iter().all(|peak| peak.charge == Some(1)));
    assert!(conflicting.iter().all(|peak| peak.envelope.is_none()));
}

#[test]
fn process_ms2_uses_fragment_charge_without_isotope_collapsing() {
    let processor = SpectrumProcessor::new(10, false, 0.0);
    let spectrum = RawSpectrum {
        ms_level: 2,
        representation: Representation::Centroid,
        mz: vec![400.0],
        intensity: vec![25.0],
        fragment_charges: Some(vec![2]),
        ..RawSpectrum::default_with_file_id(7)
    };

    let processed = processor.process(spectrum);

    assert_eq!(processed.charges, vec![2]);
    assert_eq!(processed.charge_is_known, vec![true]);
    assert!((processed.masses[0] - (400.0 - PROTON) * 2.0).abs() < 0.001);
    assert!((processed.peak_mz(0) - 400.0).abs() < 0.001);
}

#[test]
fn select_most_intense_peak_uses_parallel_columns() {
    let masses = vec![99.0, 100.0, 100.01, 100.02, 101.0];
    let intensities = vec![10.0, 20.0, 50.0, 30.0, 100.0];

    let idx = select_most_intense_peak(
        &masses,
        &intensities,
        100.01,
        Tolerance::Da(-0.02, 0.02),
        None,
    )
    .expect("peak in tolerance");

    assert_eq!(idx, 2);
    assert_eq!(masses[idx], 100.01);
    assert_eq!(intensities[idx], 50.0);
}

#[test]
fn select_most_intense_peak_applies_offset() {
    let label = 126.127726;
    let masses = vec![label - PROTON - 0.01, label - PROTON, label - PROTON + 0.01];
    let intensities = vec![10.0, 100.0, 50.0];

    let idx = select_most_intense_peak(
        &masses,
        &intensities,
        label,
        Tolerance::Da(-0.005, 0.005),
        Some(-PROTON),
    )
    .expect("offset peak in tolerance");

    assert_eq!(idx, 1);
}

#[test]
fn select_most_intense_peak_includes_bounds_and_uses_later_tie() {
    let masses = vec![99.99, 100.0, 100.01, 100.02, 100.03];
    let intensities = vec![100.0, 40.0, 50.0, 50.0, 100.0];

    let idx = select_most_intense_peak(
        &masses,
        &intensities,
        100.01,
        Tolerance::Da(-0.01, 0.01),
        None,
    );

    assert_eq!(idx, Some(3));
}

#[test]
fn select_most_intense_peak_matches_two_bound_search_reference() {
    fn reference(
        masses: &[f32],
        intensities: &[f32],
        center: f32,
        tolerance: Tolerance,
        offset: Option<f32>,
    ) -> Option<usize> {
        let (lo, hi) = tolerance.bounds(center);
        let offset = offset.unwrap_or_default();
        let lo = lo + offset;
        let hi = hi + offset;
        let left = masses
            .partition_point(|mass| mass.total_cmp(&lo).is_lt())
            .saturating_sub(1);
        let right = masses[left..].partition_point(|mass| !mass.total_cmp(&hi).is_gt()) + left;

        let mut best_peak = None;
        let mut max_int = 0.0;
        for idx in (left..right).filter(|&idx| masses[idx] >= lo && masses[idx] <= hi) {
            if intensities[idx] >= max_int {
                max_int = intensities[idx];
                best_peak = Some(idx);
            }
        }
        best_peak
    }

    let signed_zero_masses = vec![-0.0, 0.0];
    let signed_zero_intensities = vec![100.0, 10.0];
    assert_eq!(
        select_most_intense_peak(
            &signed_zero_masses,
            &signed_zero_intensities,
            0.0,
            Tolerance::Da(0.0, 0.0),
            None,
        ),
        reference(
            &signed_zero_masses,
            &signed_zero_intensities,
            0.0,
            Tolerance::Da(0.0, 0.0),
            None,
        )
    );

    let masses = (0..400)
        .map(|idx| 50.0 + idx as f32 * 0.037)
        .collect::<Vec<_>>();
    let intensities = (0..masses.len())
        .map(|idx| ((idx * 37) % 101) as f32)
        .collect::<Vec<_>>();
    let tolerances = [
        Tolerance::Da(-0.02, 0.02),
        Tolerance::Ppm(-20.0, 20.0),
        Tolerance::Pct(-0.01, 0.01),
    ];

    for center_step in 0..600 {
        let center = 49.5 + center_step as f32 * 0.027;
        for tolerance in tolerances {
            for offset in [None, Some(-PROTON), Some(0.013)] {
                assert_eq!(
                    select_most_intense_peak(&masses, &intensities, center, tolerance, offset,),
                    reference(&masses, &intensities, center, tolerance, offset,),
                    "center={center}, tolerance={tolerance:?}, offset={offset:?}"
                );
            }
        }
    }
}

#[test]
fn process_ms1_without_mobility_builds_empty_mobility_column() {
    let processor = SpectrumProcessor::new(10, false, 0.0);
    let spectrum = RawSpectrum {
        ms_level: 1,
        mz: vec![102.0, 100.0, 101.0],
        intensity: vec![30.0, 10.0, 20.0],
        ..RawSpectrum::default_with_file_id(7)
    };

    let processed = processor.process(spectrum);

    assert_eq!(processed.file_id, 7);
    assert_eq!(
        processed.masses,
        vec![100.0 - PROTON, 101.0 - PROTON, 102.0 - PROTON]
    );
    assert_eq!(processed.intensities, vec![10.0, 20.0, 30.0]);
    assert_eq!(processed.charges, vec![1, 1, 1]);
    assert!(processed.mobilities.is_empty());
    assert_eq!(processed.total_ion_current, 60.0);
}

#[test]
fn sorted_ms1_reuses_raw_peak_allocations() {
    let processor = SpectrumProcessor::new(10, false, 0.0);
    let spectrum = RawSpectrum {
        ms_level: 1,
        mz: vec![100.0, 101.0, 102.0],
        intensity: vec![10.0, 20.0, 30.0],
        ..RawSpectrum::default()
    };
    let mass_ptr = spectrum.mz.as_ptr();
    let intensity_ptr = spectrum.intensity.as_ptr();

    let processed = processor.process(spectrum);

    assert_eq!(processed.masses.as_ptr(), mass_ptr);
    assert_eq!(processed.intensities.as_ptr(), intensity_ptr);
}

#[test]
fn process_ms1_with_mobility_sorts_all_columns_by_mass() {
    let processor = SpectrumProcessor::new(10, false, 0.0);
    let spectrum = RawSpectrum {
        ms_level: 1,
        mz: vec![102.0, 100.0, 101.0],
        intensity: vec![30.0, 10.0, 20.0],
        mobility: Some(vec![3.0, 1.0, 2.0]),
        ..RawSpectrum::default_with_file_id(7)
    };

    let processed = processor.process(spectrum);

    assert_eq!(
        processed.masses,
        vec![100.0 - PROTON, 101.0 - PROTON, 102.0 - PROTON]
    );
    assert_eq!(processed.intensities, vec![10.0, 20.0, 30.0]);
    assert_eq!(processed.charges, vec![1, 1, 1]);
    assert_eq!(processed.mobilities, vec![1.0, 2.0, 3.0]);
    assert_eq!(processed.masses.len(), processed.intensities.len());
    assert_eq!(processed.masses.len(), processed.charges.len());
    assert_eq!(processed.masses.len(), processed.mobilities.len());
}

#[test]
fn process_ms2_without_deisotoping_defaults_charges_to_one() {
    let processor = SpectrumProcessor::new(10, false, 0.0);
    let spectrum = RawSpectrum {
        ms_level: 2,
        representation: Representation::Centroid,
        mz: vec![102.0, 100.0, 101.0],
        intensity: vec![30.0, 10.0, 20.0],
        ..RawSpectrum::default_with_file_id(7)
    };

    let processed = processor.process(spectrum);

    assert_eq!(
        processed.masses,
        vec![100.0 - PROTON, 101.0 - PROTON, 102.0 - PROTON]
    );
    assert_eq!(processed.intensities, vec![10.0, 20.0, 30.0]);
    assert_eq!(processed.charges, vec![1, 1, 1]);
    assert_eq!(processed.charge_is_known, vec![false, false, false]);
    assert_eq!(processed.peak_mz(1), 101.0);
}

#[test]
fn sorted_ms2_without_deisotoping_reuses_raw_peak_allocations() {
    let processor = SpectrumProcessor::new(10, false, 0.0);
    let spectrum = RawSpectrum {
        ms_level: 2,
        representation: Representation::Centroid,
        mz: vec![100.0, 101.0, 102.0],
        intensity: vec![10.0, 20.0, 30.0],
        ..RawSpectrum::default()
    };
    let mass_ptr = spectrum.mz.as_ptr();
    let intensity_ptr = spectrum.intensity.as_ptr();

    let processed = processor.process(spectrum);

    assert_eq!(processed.masses.as_ptr(), mass_ptr);
    assert_eq!(processed.intensities.as_ptr(), intensity_ptr);
}

#[test]
fn process_ms2_with_deisotoping_tracks_reassigned_charge() {
    let processor = SpectrumProcessor::new(10, true, 0.0);
    let spectrum = RawSpectrum {
        ms_level: 2,
        representation: Representation::Centroid,
        mz: vec![812.0, 812.0 + NEUTRON / 2.0],
        intensity: vec![9.0, 4.5],
        ..RawSpectrum::default_with_file_id(7)
    };

    let processed = processor.process(spectrum);

    assert_eq!(processed.masses.len(), 1);
    assert_eq!(processed.intensities, vec![13.5]);
    assert_eq!(processed.charges, vec![2]);
    assert_eq!(processed.charge_is_known, vec![true]);
    assert!((processed.masses[0] - (812.0 - PROTON) * 2.0).abs() < 0.001);
    assert!((processed.peak_mz(0) - 812.0).abs() < 0.001);
}

#[test]
fn process_ms2_with_averagine_deisotoping_tracks_reassigned_charge() {
    let processor =
        SpectrumProcessor::with_deisotope_settings(10, DeisotopeSettings::default(), 0.0);
    let neutral_mass = 1600.0;
    let charge = 2;
    let mono_mz = neutral_mass / charge as f32 + PROTON;
    let pattern = crate::isotopes::averagine_isotopes(neutral_mass);
    let spectrum = RawSpectrum {
        ms_level: 2,
        representation: Representation::Centroid,
        precursors: vec![Precursor {
            charge: Some(3),
            ..Precursor::default()
        }],
        mz: vec![mono_mz, mono_mz + NEUTRON / charge as f32],
        intensity: vec![pattern[0] * 1000.0, pattern[1] * 1000.0],
        ..RawSpectrum::default_with_file_id(7)
    };

    let processed = processor.process(spectrum);

    assert_eq!(processed.masses.len(), 1);
    assert_eq!(processed.charges, vec![charge]);
    assert_eq!(processed.charge_is_known, vec![true]);
    assert!((processed.masses[0] - neutral_mass).abs() < 0.001);
}
