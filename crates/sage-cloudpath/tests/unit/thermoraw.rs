use super::*;
use sage_plus_raw::{PrecursorInfo, ScanMode, SpectrumRecord};

/// `ThermoRawReader::convert` for a record that already owns its arrays,
/// which is how these tests build their scans.
fn converted(reader: &ThermoRawReader, record: &SpectrumRecord) -> RawSpectrum {
    reader.convert(record, &record.mz, record.intensity.clone())
}

#[test]
fn converts_raw_record() {
    let record = SpectrumRecord {
        index: 41,
        scan_number: 42,
        ms_level: 2,
        is_ms1: false,
        is_dia: false,
        is_wideband: false,
        polarity: None,
        // The reader retains the nominal instrument mode even when asked
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
        dropped_peaks: 0,
    };

    let spectrum =
        ThermoRawReader::with_file_id(7).convert(&record, &record.mz, record.intensity.clone());
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
fn implausible_precursor_mz_is_not_searched() {
    let record = |selected_mz: f64, target_mz: Option<f64>| SpectrumRecord {
        index: 0,
        scan_number: 217,
        ms_level: 2,
        is_ms1: false,
        is_dia: false,
        is_wideband: false,
        polarity: None,
        scan_mode: None,
        filter: None,
        retention_time_min: 1.0,
        total_ion_current: 1.0,
        base_peak_mz: 200.0,
        base_peak_intensity: 1.0,
        low_mz: 100.0,
        high_mz: 1000.0,
        ion_injection_time_ms: None,
        faims_cv: None,
        precursor: Some(PrecursorInfo {
            selected_mz: Some(selected_mz),
            target_mz,
            ..Default::default()
        }),
        mz: vec![200.0],
        intensity: vec![1.0],
        dropped_peaks: 0,
    };
    let reader = ThermoRawReader::with_file_id(0);
    // A denormal precursor m/z is not searched.
    let implausible = record(2.1e-314, None);
    assert!(converted(&reader, &implausible).precursors.is_empty());
    // A plausible isolation target is used when the selected m/z is not.
    let target_only = record(0.0, Some(810.16));
    let spectrum = converted(&reader, &target_only);
    assert_eq!(spectrum.precursors[0].mz, 810.16);
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
        dropped_peaks: 0,
    }
}

#[test]
fn undecodable_scans_are_skipped() {
    // The read loop keeps going past a scan it cannot decode, counting the
    // failures and keeping the first error for the warning it logs.
    let mut skipped = None;
    assert!(skipped.is_none());

    note_skipped(&mut skipped, sage_plus_raw::Error::CorruptData("scan 2"));
    note_skipped(&mut skipped, sage_plus_raw::Error::CorruptData("scan 4"));
    let (failed, error) = skipped.unwrap();
    assert_eq!(failed, 2);
    assert_eq!(error.to_string(), "corrupt data: scan 2");
}

#[test]
fn plausible_precursor_bounds() {
    assert!(plausible_precursor_mz(806.96));
    assert!(!plausible_precursor_mz(7.7e-304));
    assert!(!plausible_precursor_mz(f64::NAN));
}

#[test]
fn ms3_without_precursor_mz_keeps_its_parent_scan() {
    let mut ms3 = record(9, 3);
    ms3.precursor = Some(PrecursorInfo {
        master_scan_number: Some(8),
        ..Default::default()
    });
    let spectrum = converted(&ThermoRawReader::with_file_id(0), &ms3);
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
        .convert(&ms2, &ms2.mz, ms2.intensity.clone())
        .precursors
        .is_empty());
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
        let mut groups = std::collections::BTreeMap::new();
        for spectrum in spectra.iter().filter(|spectrum| spectrum.ms_level == 2) {
            *groups.entry(spectrum.acquisition.label()).or_insert(0usize) += 1;
        }
        println!(
            "{path}: levels {observed:?}, {missing} MS2 scans without a precursor, \
             {unlinked} MSn scans above MS2 without a parent scan, MS2 groups {groups:?}"
        );
    }
}

#[test]
fn records_without_a_filter_have_an_unknown_acquisition_group() {
    use sage_core::spectrum::{AcquisitionGroup, MassAnalyzer};
    // A scan without a filter has an unknown group rather than a guess,
    // and unknown is not low accuracy.
    let unfiltered = record(5, 2);
    let spectrum = converted(&ThermoRawReader::with_file_id(0), &unfiltered);
    assert_eq!(spectrum.acquisition, AcquisitionGroup::default());
    assert_eq!(spectrum.acquisition.analyzer, MassAnalyzer::Unknown);
    assert!(!spectrum.acquisition.analyzer.is_low_accuracy());
}
