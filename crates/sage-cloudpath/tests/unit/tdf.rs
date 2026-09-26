use super::*;
use crate::tims_mobility::{BrukerMobilityScale, MobilityCalibration};

fn buffer(peaks: Vec<ImsPeak>) -> PeakBuffer {
    let mut order = (0..peaks.len()).collect::<Vec<_>>();
    order.sort_unstable_by(|left, right| {
        peaks[*right]
            .intensity
            .partial_cmp(&peaks[*left].intensity)
            .unwrap()
    });
    PeakBuffer {
        peaks,
        order,
        agg_buff: Vec::new(),
    }
}

#[test]
fn parses_precursor_metadata_without_losing_optional_values() {
    use timsrust::core::{Charge, FrameIndex, Im, Mz, Rt, ScanIndex};

    let precursor = timsrust::core::Precursor::new(
        Mz::from(500.25),
        Im::from(1.15),
        Rt::from(12.0),
        ScanIndex::try_from(5).unwrap(),
        Some(Charge::try_from(3).unwrap()),
        Some(1234.5),
        7,
        FrameIndex::try_from(42).unwrap(),
    );
    let parsed = TdfReader::parse_precursor(&precursor);

    assert_eq!(parsed.mz, 500.25);
    assert_eq!(parsed.charge, Some(3));
    assert_eq!(parsed.intensity, Some(1234.5));
    assert_eq!(parsed.spectrum_ref.as_deref(), Some("42"));
    assert_eq!(parsed.inverse_ion_mobility, Some(1.15));
}

#[test]
fn parses_real_bruker_directory() {
    let path = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/bruker/example_dia.d");
    let url = crate::Url::from_file_path(path).unwrap();
    let spectra =
        crate::util::read_spectra(&url, 3, None, BrukerProcessingConfig::default(), true).unwrap();

    assert!(!spectra.is_empty());
    assert!(spectra.iter().all(|spectrum| {
        spectrum.file_id == 3 && spectrum.mz.len() == spectrum.intensity.len()
    }));
    assert!(spectra.iter().any(|spectrum| spectrum.ms_level == 1));
    assert!(spectra.iter().any(|spectrum| spectrum.ms_level > 1));
    assert!(spectra
        .iter()
        .filter(|spectrum| spectrum.ms_level == 1)
        .all(|spectrum| spectrum
            .mobility
            .as_ref()
            .is_some_and(|values| { !values.is_empty() && values.len() == spectrum.mz.len() })));
    assert!(spectra
        .iter()
        .filter(|spectrum| spectrum.ms_level > 1)
        .all(|spectrum| !spectrum.precursors.is_empty()));
}

#[test]
fn mobility_offsets_expand_run_lengths_and_skip_empty_scans() {
    let scan_to_im = |scan: usize| 2.0 - scan as f32 * 0.25;
    let mobility = PeakBuffer::expand_mobility_iter(&[0, 2, 2, 5], &scan_to_im).collect::<Vec<_>>();

    assert_eq!(mobility, vec![2.0, 2.0, 1.5, 1.5, 1.5]);
}

#[test]
fn mobility_scale_can_be_set_alone() {
    let config: BrukerProcessingConfig =
        serde_json::from_str(r#"{"ion_mobility_scale": "linear"}"#).unwrap();
    assert_eq!(config.ion_mobility_scale, BrukerMobilityScale::Linear);
    let config: BrukerProcessingConfig = serde_json::from_str("{}").unwrap();
    assert_eq!(config.ion_mobility_scale, BrukerMobilityScale::Calibrated);
}

#[test]
fn calibrated_mobility_differs_from_linear_scale() {
    let path = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/bruker/example_dia.d");
    let url = crate::Url::from_file_path(&path).unwrap();
    let read = |scale| {
        let config = BrukerProcessingConfig {
            ion_mobility_scale: scale,
            ..BrukerProcessingConfig::default()
        };
        crate::util::read_spectra(&url, 0, None, config, true).unwrap()
    };
    let calibrated = read(BrukerMobilityScale::Calibrated);
    let linear = read(BrukerMobilityScale::Linear);
    assert_eq!(calibrated.len(), linear.len());

    let model = *MobilityCalibration::from_path(&path).unwrap().dominant();
    // Acquisition range limits used by the linear scale.
    let (lower, upper) = (0.6f32, 1.6f32);
    // Spectrum order is not stable between reads, so pair by level and id.
    let linear = linear
        .iter()
        .map(|spectrum| ((spectrum.ms_level, spectrum.id.clone()), spectrum))
        .collect::<std::collections::HashMap<_, _>>();
    let mut largest_shift = 0.0f32;
    for calibrated in &calibrated {
        let linear = linear[&(calibrated.ms_level, calibrated.id.clone())];
        if calibrated.ms_level == 1 {
            // MS1 centroiding merges peaks within a mobility tolerance, so the
            // peak lists can differ. Compare mean mobility instead.
            let mean = |values: &Vec<f32>| values.iter().sum::<f32>() / values.len().max(1) as f32;
            let (a, b) = (
                calibrated.mobility.as_ref().unwrap(),
                linear.mobility.as_ref().unwrap(),
            );
            largest_shift = largest_shift.max((mean(a) - mean(b)).abs());
            continue;
        }
        // MS2 spectra do not depend on the mobility scale.
        assert_eq!(calibrated.mz, linear.mz);
        assert_eq!(calibrated.intensity, linear.intensity);
        for (a, b) in calibrated.precursors.iter().zip(&linear.precursors) {
            let (a, b) = (
                a.inverse_ion_mobility.unwrap(),
                b.inverse_ion_mobility.unwrap(),
            );
            assert!((lower..=upper).contains(&b));
            // Invert timsrust's scale, linear in sqrt(1/K0) over 918 scans, to
            // recover the window-center scan. DIA centers use the dominant row.
            let slope = (lower.sqrt() - upper.sqrt()) / 918.0;
            let scan = ((b.sqrt() - upper.sqrt()) / slope).round() as f64;
            assert!((a - model.one_over_k0(scan) as f32).abs() < 1e-4);
        }
    }
    // On this run the scales differ by up to 0.087 1/K0 per scan (at scan 395).
    assert!(largest_shift > 1e-3, "largest shift {largest_shift}");
    assert!(largest_shift < 0.15, "largest shift {largest_shift}");
}

#[test]
fn centroiding_combines_nearby_mass_and_mobility_peaks() {
    let mut buffer = buffer(vec![
        ImsPeak {
            mz: 100.0,
            intensity: 10.0,
            im: 1.0,
        },
        ImsPeak {
            mz: 100.0004,
            intensity: 5.0,
            im: 1.01,
        },
        ImsPeak {
            mz: 101.0,
            intensity: 7.0,
            im: 1.2,
        },
    ]);

    let (mz, (intensity, mobility)) = buffer.fastcentroid_frame(5.0, 2.0);

    assert_eq!(mz, vec![100.0, 101.0]);
    assert_eq!(intensity, vec![15.0, 7.0]);
    assert_eq!(mobility, vec![1.0, 1.2]);
}

#[test]
fn centroiding_keeps_close_masses_separate_when_mobility_differs() {
    let mut buffer = buffer(vec![
        ImsPeak {
            mz: 100.0,
            intensity: 10.0,
            im: 1.0,
        },
        ImsPeak {
            mz: 100.0004,
            intensity: 5.0,
            im: 1.2,
        },
    ]);

    let (mz, (intensity, mobility)) = buffer.fastcentroid_frame(5.0, 2.0);

    assert_eq!(mz, vec![100.0, 100.0004]);
    assert_eq!(intensity, vec![10.0, 5.0]);
    assert_eq!(mobility, vec![1.0, 1.2]);
}

#[test]
fn clear_resets_all_reusable_storage() {
    let mut buffer = buffer(vec![ImsPeak {
        mz: 100.0,
        intensity: 10.0,
        im: 1.0,
    }]);
    buffer.agg_buff.push(buffer.peaks[0]);

    buffer.clear();

    assert!(buffer.peaks.is_empty());
    assert!(buffer.order.is_empty());
    assert!(buffer.agg_buff.is_empty());
}

#[test]
fn calibrated_mobility_reads_inputs_that_point_inside_the_directory() {
    let directory = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/bruker/example_dia.d");
    let read = |path: &Path| {
        let url = crate::Url::from_file_path(path).unwrap();
        let mut spectra =
            crate::util::read_spectra(&url, 0, None, BrukerProcessingConfig::default(), false)
                .unwrap();
        spectra.sort_by(|a, b| (a.ms_level, &a.id).cmp(&(b.ms_level, &b.id)));
        spectra
    };
    let expected = read(&directory);
    for file in ["analysis.tdf_bin", "analysis.tdf"] {
        let spectra = read(&directory.join(file));
        assert_eq!(spectra.len(), expected.len(), "{file}");
        for (actual, expected) in spectra.iter().zip(&expected) {
            assert_eq!(actual.id, expected.id, "{file}");
            let mobility = |spectrum: &RawSpectrum| {
                spectrum
                    .precursors
                    .iter()
                    .map(|precursor| precursor.inverse_ion_mobility)
                    .collect::<Vec<_>>()
            };
            assert_eq!(mobility(actual), mobility(expected), "{file}");
        }
    }
}

#[test]
fn calibrated_mobility_falls_back_to_linear_without_analysis_tdf() {
    let directory = tempfile::tempdir().unwrap();
    let path = directory.path().join("run.ms2");
    std::fs::create_dir(&path).unwrap();
    assert!(matches!(
        MobilityScale::new(&path, BrukerMobilityScale::Calibrated).unwrap(),
        MobilityScale::Linear(None)
    ));
    assert_eq!(
        BrukerMobilityScale::Calibrated.effective_for(&path),
        BrukerMobilityScale::Linear
    );
    let tdf = Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data/bruker/example_dia.d/analysis.tdf_bin");
    assert_eq!(
        BrukerMobilityScale::Calibrated.effective_for(&tdf),
        BrukerMobilityScale::Calibrated
    );
    assert_eq!(
        BrukerMobilityScale::Linear.effective_for(&tdf),
        BrukerMobilityScale::Linear
    );
}

#[test]
fn linear_scale_keeps_the_beta6_conversion() {
    // timsrust 0.6.6 changed its uncalibrated conversion; the linear scale
    // keeps the sqrt(1/K0) interpolation of timsrust 0.6.5 and Beta 6.
    let path = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/bruker/example_dia.d");
    let linear = LinearMobilityScale::from_path(&path).unwrap();
    assert_eq!(linear, LinearMobilityScale::new(0.6, 1.6, 918));
    let (lower, upper) = (0.6f64, 1.6f64);
    for scan in [0u32, 1, 395, 917, 918] {
        let expected = upper.sqrt() + (lower.sqrt() - upper.sqrt()) / 918.0 * f64::from(scan);
        assert_eq!(linear.one_over_k0(scan), expected * expected);
    }
    assert!(matches!(
        MobilityScale::new(&path, BrukerMobilityScale::Linear).unwrap(),
        MobilityScale::Linear(Some(scale)) if scale == linear
    ));
}

#[test]
fn diapasef_reads_one_spectrum_per_window_split() {
    // Regression test: timsrust 0.6 (Beta 2) read diaPASEF as precursor-anchored spectra
    // and ignored `ms2.frame_splitting_params`.
    let path = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/bruker/example_dia.d");
    let url = crate::Url::from_file_path(path).unwrap();
    let read = |frame_splitting_params| {
        let config = BrukerProcessingConfig {
            ms2: BrukerSpectrumConfig {
                frame_splitting_params,
                ..BrukerSpectrumConfig::default()
            },
            ..BrukerProcessingConfig::default()
        };
        crate::util::read_spectra(&url, 0, None, config, false).unwrap()
    };

    // 470 MS2 frames, each with two 26 m/z windows centered at 413 + 25k.
    let spectra = read(BrukerFrameWindowSplittingConfig::default());
    assert_eq!(spectra.len(), 940);
    for spectrum in &spectra {
        assert_eq!(spectrum.ms_level, 2);
        let precursor = &spectrum.precursors[0];
        let step = (precursor.mz - 413.0) / 25.0;
        assert!((step - step.round()).abs() < 1e-3, "{}", precursor.mz);
        assert_eq!(precursor.isolation_window, Some(Tolerance::Da(-13.0, 13.0)));
        assert_eq!(precursor.charge, None);
        assert!(precursor.inverse_ion_mobility.is_some());
    }

    // The splitting config is honored: two mobility splits per window double the count.
    let split = read(BrukerFrameWindowSplittingConfig::Quadrupole(
        BrukerQuadWindowExpansionStrategy::Even(2),
    ));
    assert_eq!(split.len(), 2 * spectra.len());
}
