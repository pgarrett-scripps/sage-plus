use super::*;
use opentfraw::{PrecursorInfo, ScanMode, SpectrumRecord};

#[test]
fn converts_open_tf_raw_record() {
    let record = SpectrumRecord {
        index: 41,
        scan_number: 42,
        ms_level: 2,
        is_ms1: false,
        is_dia: false,
        is_wideband: false,
        polarity: None,
        // OpenTFRaw retains the nominal instrument mode even when asked
        // to return the centroid list.
        scan_mode: Some(ScanMode::Profile),
        filter: None,
        retention_time_min: 12.5,
        total_ion_current: 1234.0,
        base_peak_mz: 200.0,
        base_peak_intensity: 1000.0,
        low_mz: 100.0,
        high_mz: 1000.0,
        ion_injection_time_ms: Some(8.25),
        faims_cv: None,
        precursor: Some(PrecursorInfo {
            target_mz: Some(500.2),
            selected_mz: Some(500.25),
            isolation_width: Some(1.6),
            charge: Some(2),
            master_scan_number: Some(40),
            ..Default::default()
        }),
        mz: vec![100.1, 200.2],
        intensity: vec![10.0, 20.0],
    };

    let spectrum = ThermoRawReader::with_file_id(7).convert(record);
    assert_eq!(spectrum.file_id, 7);
    assert_eq!(spectrum.id, "controllerType=0 controllerNumber=1 scan=42");
    assert_eq!(spectrum.ms_level, 2);
    assert_eq!(spectrum.representation, Representation::Centroid);
    assert_eq!(spectrum.precursors[0].mz, 500.25);
    assert_eq!(spectrum.precursors[0].charge, Some(2));
    assert_eq!(
        spectrum.precursors[0].isolation_window,
        Some(Tolerance::Da(-0.8, 0.8))
    );
    assert_eq!(spectrum.mz, vec![100.1, 200.2]);
    assert_eq!(spectrum.intensity, vec![10.0, 20.0]);
}

#[test]
fn parses_real_raw_file() {
    let path =
        Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/thermo/Angiotensin_325-CID.raw");
    let url = crate::Url::from_file_path(path).unwrap();
    let spectra = crate::util::read_spectra(
        &url,
        3,
        None,
        crate::tdf::BrukerProcessingConfig::default(),
        false,
    )
    .unwrap();

    assert_eq!(spectra.len(), 10);
    assert!(spectra.iter().all(|spectrum| {
        spectrum.file_id == 3
            && spectrum.ms_level == 2
            && spectrum.representation == Representation::Centroid
            && !spectrum.mz.is_empty()
            && spectrum.mz.len() == spectrum.intensity.len()
    }));

    let first = &spectra[0];
    assert_eq!(first.id, "controllerType=0 controllerNumber=1 scan=1");
    assert!((first.total_ion_current - 37_687_076.0).abs() < 100.0);
    assert!((first.ion_injection_time - 7.422).abs() < 0.001);
    assert_eq!(first.precursors.len(), 1);
    assert_eq!(first.precursors[0].mz, 325.0);
    assert_eq!(first.precursors[0].charge, Some(1));
}
