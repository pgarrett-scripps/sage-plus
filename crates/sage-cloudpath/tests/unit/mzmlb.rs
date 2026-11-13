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
