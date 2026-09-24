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

fn record(scan_number: u32, ms_level: u32) -> SpectrumRecord {
    SpectrumRecord {
        index: scan_number as usize - 1,
        scan_number,
        ms_level,
        is_ms1: ms_level == 1,
        is_dia: false,
        is_wideband: false,
        polarity: None,
        scan_mode: None,
        filter: None,
        retention_time_min: 0.0,
        total_ion_current: 0.0,
        base_peak_mz: 0.0,
        base_peak_intensity: 0.0,
        low_mz: 0.0,
        high_mz: 0.0,
        ion_injection_time_ms: None,
        faims_cv: None,
        precursor: None,
        mz: Vec::new(),
        intensity: Vec::new(),
    }
}

#[test]
fn trailer_levels_follow_master_scans() {
    // MS1, MS2 of scan 1, MS3 of scan 2, MS2 of scan 1, missing trailer,
    // and a master that is not an earlier scan.
    let masters = [Some(0), Some(1), Some(2), Some(1), None, Some(9)];
    assert_eq!(
        trailer_levels(1, &masters),
        vec![Some(1), Some(2), Some(3), Some(2), None, None]
    );
}

#[test]
fn misaligned_events_defer_to_trailers() {
    // Events shifted by one scan, as decoded from some Orbitrap Fusion files.
    let masters = (0..1000)
        .map(|idx| Some(if idx % 4 == 0 { 0 } else { idx / 4 * 4 + 1 }))
        .collect::<Vec<_>>();
    let levels = trailer_levels(1, &masters);
    let shifted = (1..=1000)
        .map(|scan| record(scan, if scan % 4 == 2 { 1 } else { 2 }))
        .collect::<Vec<_>>();
    let disagreements = event_level_disagreements(1, &shifted, &levels);
    assert_eq!(disagreements, 500);
    assert!(events_misaligned(disagreements, shifted.len()));

    let aligned = (1..=1000)
        .map(|scan| record(scan, levels[scan as usize - 1].unwrap()))
        .collect::<Vec<_>>();
    assert_eq!(event_level_disagreements(1, &aligned, &levels), 0);
    assert!(!events_misaligned(10, 1000));
}

#[test]
fn ms3_without_precursor_mz_keeps_its_parent_scan() {
    let mut ms3 = record(9, 3);
    ms3.precursor = Some(PrecursorInfo {
        master_scan_number: Some(8),
        ..Default::default()
    });
    let spectrum = ThermoRawReader::with_file_id(0).convert(ms3);
    assert_eq!(
        spectrum.precursors[0].spectrum_ref.as_deref(),
        Some("controllerType=0 controllerNumber=1 scan=8")
    );

    let mut ms2 = record(8, 2);
    ms2.precursor = Some(PrecursorInfo {
        master_scan_number: Some(7),
        ..Default::default()
    });
    assert!(ThermoRawReader::with_file_id(0)
        .convert(ms2)
        .precursors
        .is_empty());
}

#[test]
fn trailer_level_replaces_event_precursor() {
    let mut ms1 = record(5, 2);
    ms1.precursor = Some(PrecursorInfo {
        selected_mz: Some(600.0),
        ..Default::default()
    });
    apply_trailer_level(&mut ms1, 1, None);
    assert_eq!(ms1.ms_level, 1);
    assert!(ms1.is_ms1);
    assert!(ms1.precursor.is_none());
}

/// Checks MS level counts on local files that are too large for the
/// repository. `SAGE_THERMO_LEVEL_CHECK` holds `path=ms1,ms2,ms3` entries
/// separated by `;`, with counts taken from a vendor mzML conversion.
#[test]
#[ignore]
fn raw_ms_levels_match_vendor_conversion() {
    let spec = std::env::var("SAGE_THERMO_LEVEL_CHECK").expect("SAGE_THERMO_LEVEL_CHECK");
    for entry in spec.split(';').filter(|entry| !entry.is_empty()) {
        let (path, counts) = entry.split_once('=').unwrap();
        let expected = counts
            .split(',')
            .map(|count| count.parse::<usize>().unwrap())
            .collect::<Vec<_>>();
        let spectra = ThermoRawReader::with_file_id(0).parse(path).unwrap();
        let observed = (1..=expected.len())
            .map(|level| {
                spectra
                    .iter()
                    .filter(|spectrum| spectrum.ms_level as usize == level)
                    .count()
            })
            .collect::<Vec<_>>();
        assert_eq!(observed, expected, "{path}");
        let missing = spectra
            .iter()
            .filter(|spectrum| spectrum.ms_level == 2 && spectrum.precursors.is_empty())
            .count();
        let unlinked = spectra
            .iter()
            .filter(|spectrum| {
                spectrum.ms_level > 2
                    && spectrum
                        .precursors
                        .first()
                        .and_then(|precursor| precursor.spectrum_ref.as_ref())
                        .is_none()
            })
            .count();
        println!(
            "{path}: levels {observed:?}, {missing} MS2 scans without a precursor, \
             {unlinked} MSn scans above MS2 without a parent scan"
        );
    }
}
