use super::*;
use mzdata::io::{
    mzmlb::{MzMLbWriter, MzMLbWriterBuilder},
    SpectrumWriter,
};
use mzdata::spectrum::{
    bindata::{ArrayType, BinaryArrayMap, BinaryDataArrayType, DataArray},
    Acquisition, IsolationWindow, IsolationWindowState, Precursor as MzDataPrecursor,
    RawSpectrum as MzDataRawSpectrum, ScanEvent, SelectedIon, SpectrumDescription,
};
use mzdata::{mzpeaks::CentroidPeak, mzpeaks::MZPeakSetType};

#[test]
fn reads_mzmlb_spectra_and_precursor_metadata() -> io::Result<()> {
    let directory = tempfile::tempdir()?;
    let path = directory.path().join("tiny.mzMLb");
    let peaks = MZPeakSetType::new(vec![
        CentroidPeak::new(100.0, 10.0, 0),
        CentroidPeak::new(200.0, 20.0, 1),
    ]);
    let mut acquisition = Acquisition::default();
    acquisition.scans.push(ScanEvent {
        start_time: 12.5,
        injection_time: 7.0,
        ..ScanEvent::default()
    });
    let mut arrays = BinaryArrayMap::from(&peaks);
    let mut charge_array =
        DataArray::from_name_and_type(&ArrayType::ChargeArray, BinaryDataArrayType::Int32);
    charge_array.extend(&[2i32, 0]).unwrap();
    arrays.add(charge_array);
    let spectrum = MzDataRawSpectrum {
        description: SpectrumDescription {
            id: "scan=42".into(),
            ms_level: 2,
            signal_continuity: SignalContinuity::Centroid,
            acquisition,
            precursor: vec![MzDataPrecursor {
                ions: vec![SelectedIon {
                    mz: 500.25,
                    intensity: 123.0,
                    charge: Some(2),
                    ..SelectedIon::default()
                }],
                isolation_window: IsolationWindow::new(
                    500.0,
                    499.5,
                    500.75,
                    IsolationWindowState::Complete,
                ),
                precursor_id: Some("scan=41".into()),
                ..MzDataPrecursor::default()
            }],
            ..SpectrumDescription::default()
        },
        arrays,
    };
    let mut writer: MzMLbWriter = MzMLbWriterBuilder::new(&path).create()?;
    writer.write(&spectrum)?;
    writer.close().map_err(io::Error::from)?;

    let spectra = MzMLbReader::with_file_id(3).parse(&path)?;
    assert_eq!(spectra.len(), 1);
    let spectrum = &spectra[0];
    assert_eq!(spectrum.file_id, 3);
    assert_eq!(spectrum.id, "scan=42");
    assert_eq!(spectrum.scan_start_time, 12.5);
    assert_eq!(spectrum.ion_injection_time, 7.0);
    assert_eq!(spectrum.mz, [100.0, 200.0]);
    assert_eq!(spectrum.intensity, [10.0, 20.0]);
    assert_eq!(spectrum.fragment_charges, Some(vec![2, 0]));
    assert_eq!(spectrum.precursors[0].mz, 500.25);
    assert_eq!(spectrum.precursors[0].charge, Some(2));
    assert_eq!(
        spectrum.precursors[0].isolation_window,
        Some(Tolerance::Da(-0.5, 0.75))
    );
    Ok(())
}

fn written_spectrum(id: &str, ms_level: u8, start_time: f64, mz: f64) -> MzDataRawSpectrum {
    let peaks = MZPeakSetType::new(vec![
        CentroidPeak::new(mz, 10.0, 0),
        CentroidPeak::new(mz + 100.0, 20.0, 1),
    ]);
    let mut acquisition = Acquisition::default();
    acquisition.scans.push(ScanEvent {
        start_time,
        injection_time: 3.0,
        ..ScanEvent::default()
    });
    let precursor = (ms_level == 2)
        .then(|| MzDataPrecursor {
            ions: vec![SelectedIon {
                mz: 600.5,
                charge: Some(2),
                ..SelectedIon::default()
            }],
            ..MzDataPrecursor::default()
        })
        .into_iter()
        .collect();
    MzDataRawSpectrum {
        description: SpectrumDescription {
            id: id.into(),
            ms_level,
            signal_continuity: SignalContinuity::Centroid,
            acquisition,
            precursor,
            ..SpectrumDescription::default()
        },
        arrays: BinaryArrayMap::from(&peaks),
    }
}

#[test]
fn ms1_arrays_are_not_read_when_the_search_reads_only_ms2() -> io::Result<()> {
    let directory = tempfile::tempdir()?;
    let path = directory.path().join("dda.mzMLb");
    let mut writer: MzMLbWriter = MzMLbWriterBuilder::new(&path).create()?;
    for spectrum in [
        written_spectrum("scan=1", 1, 1.0, 400.0),
        written_spectrum("scan=2", 2, 1.5, 150.0),
        written_spectrum("scan=3", 2, 2.0, 160.0),
        written_spectrum("scan=4", 1, 2.5, 410.0),
    ] {
        writer.write(&spectrum)?;
    }
    writer.close().map_err(io::Error::from)?;

    let full = MzMLbReader::with_file_id(2).parse(&path)?;
    let skipped = MzMLbReader::with_file_id(2)
        .set_decode_ms1(false)
        .parse(&path)?;

    assert_eq!(skipped.len(), 4);
    assert_eq!(skipped.len(), full.len());
    for (skipped, full) in skipped.iter().zip(&full) {
        assert_eq!(skipped.id, full.id);
        assert_eq!(skipped.file_id, full.file_id);
        assert_eq!(skipped.ms_level, full.ms_level);
        assert_eq!(skipped.scan_start_time, full.scan_start_time);
        assert_eq!(skipped.ion_injection_time, full.ion_injection_time);
        assert_eq!(skipped.representation, full.representation);
        assert_eq!(skipped.acquisition, full.acquisition);
        if full.ms_level == 1 {
            assert!(skipped.metadata_only, "{}", skipped.id);
            assert!(skipped.mz.is_empty() && skipped.intensity.is_empty());
        } else {
            assert!(!skipped.metadata_only);
            assert_eq!(skipped.mz, full.mz);
            assert_eq!(skipped.intensity, full.intensity);
            assert_eq!(skipped.total_ion_current, full.total_ion_current);
            assert_eq!(
                format!("{:?}", skipped.precursors),
                format!("{:?}", full.precursors)
            );
        }
    }
    assert!(full.iter().all(|s| !s.metadata_only && s.mz.len() == 2));
    Ok(())
}
