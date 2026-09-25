use super::*;
use crate::tdf::BrukerProcessingConfig;

fn fixture() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/bruker/example_dia.d")
}

fn meta(id: usize, num_scans: usize, ms_ms_type: i64) -> FrameMeta {
    FrameMeta {
        id,
        num_scans,
        num_peaks: 1,
        ms_ms_type,
        rt: id as f64,
    }
}

/// A denoiser over synthetic frames, with no gates.
fn synthetic(params: &DenoiseParams, meta: Vec<FrameMeta>) -> Ms1Denoiser<'_> {
    Ms1Denoiser {
        denoiser: Denoiser::new(
            &params.filter,
            &params.stages(),
            acquisition(&meta),
            meta,
            false,
        )
        .unwrap(),
        path: PathBuf::from("synthetic.d/analysis.tdf"),
    }
}

#[test]
fn defaults_are_off_and_match_dnoise() {
    let config = BrukerDenoiseConfig::default();
    assert!(!config.enabled);
    let params = config.params();
    let filter = FilterParams::default();
    assert_eq!(params.filter.mz_half_width, filter.mz_half_width);
    assert_eq!(params.filter.min_feature_length, filter.min_feature_length);
    assert_eq!(params.filter.max_internal_gap, filter.max_internal_gap);
    assert_eq!(params.filter.num_iterations, filter.num_iterations);
    let halo = params.halo.expect("halo is on by default");
    assert_eq!(halo.peak_fraction, HaloParams::default().peak_fraction);
    assert_eq!(halo.mz_idx_half_width, 80);
    let polygon = params.polygon.expect("polygon gate is on by default");
    assert_eq!(
        (polygon.mz_pad, polygon.im_pad, polygon.overlap),
        (3.0, 0.015, true)
    );
    let dia = params.dia_ms1.expect("DIA MS1 gate is on by default");
    assert_eq!((dia.mz_pad, dia.im_pad, dia.overlap), (3.0, 0.015, true));
    assert_eq!(config.mobility_scale, BrukerMobilityScale::Calibrated);
    config.validate().unwrap();
}

#[test]
fn deserializes_partial_config_and_rejects_unknown_keys() {
    let config: BrukerDenoiseConfig =
        serde_json::from_str(r#"{"enabled": true, "halo": false, "iterations": 1}"#).unwrap();
    assert!(config.enabled);
    assert!(config.params().halo.is_none());
    assert_eq!(config.params().filter.num_iterations, 1);
    assert_eq!(config.mz_half_width, 3);

    assert!(serde_json::from_str::<BrukerDenoiseConfig>(r#"{"enable": true}"#).is_err());

    let bruker: BrukerProcessingConfig = serde_json::from_str("{}").unwrap();
    assert_eq!(bruker.denoise, BrukerDenoiseConfig::default());
    let bruker: BrukerProcessingConfig =
        serde_json::from_str(r#"{"denoise": {"enabled": true, "mobility_scale": "linear"}}"#)
            .unwrap();
    assert!(bruker.denoise.enabled);
    assert_eq!(bruker.denoise.mobility_scale, BrukerMobilityScale::Linear);
}

#[test]
fn validation_matches_dnoise() {
    let bad_halo = BrukerDenoiseConfig {
        halo_peak_fraction: 1.5,
        ..Default::default()
    };
    assert!(matches!(bad_halo.validate(), Err(DenoiseError::Config(_))));

    let bad_pad = BrukerDenoiseConfig {
        ms1_polygon_mz_pad: -1.0,
        ..Default::default()
    };
    assert!(bad_pad.validate().is_err());
    // A disabled stage is not checked, as in dnoise.
    let off = BrukerDenoiseConfig {
        ms1_polygon: false,
        ..bad_pad
    };
    off.validate().unwrap();
}

#[test]
fn classifies_acquisition_like_dnoise() {
    assert_eq!(acquisition(&[meta(1, 10, 0)]), Acquisition::Ms1Only);
    assert_eq!(
        acquisition(&[meta(1, 10, 0), meta(2, 10, 8)]),
        Acquisition::DdaPasef
    );
    assert_eq!(
        acquisition(&[meta(1, 10, 0), meta(2, 10, 9)]),
        Acquisition::DiaPasef
    );
    assert_eq!(acquisition(&[meta(1, 10, 10)]), Acquisition::PrmPasef);
    assert_eq!(
        acquisition(&[meta(1, 10, 8), meta(2, 10, 9)]),
        Acquisition::Mixed
    );
    assert_eq!(acquisition(&[meta(1, 10, 2)]), Acquisition::Unknown);
}

#[test]
fn keeps_mobility_streaks_and_drops_isolated_points() {
    let params = BrukerDenoiseConfig {
        halo: false,
        ..Default::default()
    }
    .params();
    // Frame ids need not start at 1 or be contiguous.
    let denoiser = synthetic(&params, vec![meta(5, 12, 0), meta(9, 12, 0)]);

    // An 8-scan streak at TOF 1000 and a lone point at TOF 5000 on scan 2.
    // The CSR row pointer stops at scan 9, before Frames.NumScans (12).
    let mut offsets = vec![0];
    let mut tof = Vec::new();
    for scan in 0..9 {
        if scan < 8 {
            tof.push(1000);
        }
        if scan == 2 {
            tof.push(5000);
        }
        offsets.push(tof.len());
    }
    let intensity = vec![100; tof.len()];
    let survivors = denoiser
        .denoise(9, &offsets, tof.clone(), intensity.clone())
        .unwrap();
    assert_eq!(survivors.len(), 8);
    assert!(survivors.iter().all(|&(_, tof, _)| tof == 1000));
    assert!(survivors.windows(2).all(|w| w[0].0 < w[1].0));

    assert!(matches!(
        denoiser.denoise(7, &offsets, tof.clone(), intensity.clone()),
        Err(DenoiseError::Metadata { .. })
    ));
    let too_many_scans = vec![0; 14];
    assert!(denoiser
        .denoise(5, &too_many_scans, Vec::new(), Vec::new())
        .is_err());
}

#[test]
fn builds_the_dia_window_gate_for_the_fixture() {
    let tdf = fixture().join("analysis.tdf");
    let params = BrukerDenoiseConfig::default().params();
    let denoiser = params.denoiser(&tdf).unwrap();
    assert!(denoiser.denoiser.dia_ms1_gate().is_some());
    assert!(denoiser.denoiser.polygon_gate().is_none());

    let linear = BrukerDenoiseConfig {
        mobility_scale: BrukerMobilityScale::Linear,
        ..Default::default()
    }
    .params();
    assert!(linear
        .denoiser(&tdf)
        .unwrap()
        .denoiser
        .dia_ms1_gate()
        .is_some());

    let no_gate = BrukerDenoiseConfig {
        dia_ms1_window: false,
        ..Default::default()
    }
    .params();
    assert!(no_gate
        .denoiser(&tdf)
        .unwrap()
        .denoiser
        .dia_ms1_gate()
        .is_none());
}

#[test]
fn denoises_only_ms1_spectra_when_reading_the_fixture() {
    let url = crate::Url::from_file_path(fixture()).unwrap();
    let plain =
        crate::util::read_spectra(&url, 0, None, BrukerProcessingConfig::default(), true).unwrap();
    let config = BrukerProcessingConfig {
        denoise: BrukerDenoiseConfig {
            enabled: true,
            ..Default::default()
        },
        ..Default::default()
    };
    let denoised = crate::util::read_spectra(&url, 0, None, config, true).unwrap();

    let peaks = |spectra: &[sage_core::spectrum::RawSpectrum], level: u8| {
        let mut peaks = spectra
            .iter()
            .filter(|s| s.ms_level == level)
            .map(|s| (s.id.clone(), s.mz.len()))
            .collect::<Vec<_>>();
        peaks.sort();
        peaks
    };
    assert_eq!(peaks(&plain, 2), peaks(&denoised, 2));
    let (plain_ms1, denoised_ms1) = (peaks(&plain, 1), peaks(&denoised, 1));
    assert_eq!(plain_ms1.len(), denoised_ms1.len());
    let total = |peaks: &[(String, usize)]| peaks.iter().map(|p| p.1).sum::<usize>();
    assert!(total(&denoised_ms1) < total(&plain_ms1));
    assert!(total(&denoised_ms1) > 0);
    for spectrum in denoised.iter().filter(|s| s.ms_level == 1) {
        assert!(spectrum.mz.windows(2).all(|w| w[0] <= w[1]));
        let mobility = spectrum.mobility.as_ref().unwrap();
        assert_eq!(mobility.len(), spectrum.mz.len());
        assert!(mobility.iter().all(|&im| (0.5..1.7).contains(&im)));
    }
}
