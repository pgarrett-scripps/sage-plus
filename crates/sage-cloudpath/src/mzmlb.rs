use std::io;
use std::path::Path;

use mzdata::io::MzMLbReader as MzMLbReaderImpl;
use mzdata::prelude::{IonMobilityMeasure, PrecursorSelection, SpectrumLike};
use mzdata::spectrum::{RawSpectrum as MzDataSpectrum, SignalContinuity};
use sage_core::mass::Tolerance;
use sage_core::spectrum::{Precursor, RawSpectrum, Representation};

pub struct MzMLbReader {
    file_id: usize,
}

fn normalize_fragment_charges(
    charges: &[i32],
    peak_count: usize,
    spectrum_id: &str,
) -> Option<Vec<u8>> {
    if charges.len() != peak_count {
        log::warn!(
            "ignoring fragment charge array with length different from m/z array for spectrum {}",
            spectrum_id
        );
        return None;
    }
    let converted = charges
        .iter()
        .map(|&charge| u8::try_from(charge).ok())
        .collect::<Option<Vec<_>>>();
    if converted.is_none() {
        log::warn!(
            "ignoring fragment charge array with values outside 0 through 255 for spectrum {}",
            spectrum_id
        );
    }
    converted
}

impl MzMLbReader {
    pub fn with_file_id(file_id: usize) -> Self {
        Self { file_id }
    }

    pub fn parse(&self, path: impl AsRef<Path>) -> io::Result<Vec<RawSpectrum>> {
        let reader = MzMLbReaderImpl::new(&path.as_ref().to_path_buf())?;
        Ok(reader
            .map(|spectrum| {
                let ms_level = spectrum.ms_level();
                let id = spectrum.id().to_string();
                let scan_start_time = spectrum.start_time() as f32;
                let ion_injection_time = spectrum
                    .acquisition()
                    .first_scan()
                    .map(|scan| scan.injection_time)
                    .unwrap_or_default();
                let scan_mobility = spectrum.ion_mobility().map(|value| value as f32);
                let representation = match spectrum.signal_continuity() {
                    SignalContinuity::Centroid => Representation::Centroid,
                    SignalContinuity::Profile | SignalContinuity::Unknown => {
                        Representation::Profile
                    }
                };
                let total_ion_current = spectrum.peaks().tic();
                let precursors = spectrum
                    .precursor_iter()
                    .flat_map(|precursor| {
                        let window = precursor.isolation_window();
                        let isolation_window = (!window.is_empty()).then_some({
                            Tolerance::Da(
                                window.lower_bound - window.target,
                                window.upper_bound - window.target,
                            )
                        });
                        let spectrum_ref = precursor.precursor_id().cloned();
                        let ions = precursor
                            .iter()
                            .map(|ion| Precursor {
                                mz: ion.mz as f32,
                                intensity: Some(ion.intensity),
                                charge: ion.charge.and_then(|charge| u8::try_from(charge).ok()),
                                spectrum_ref: spectrum_ref.clone(),
                                isolation_window,
                                inverse_ion_mobility: ion
                                    .ion_mobility()
                                    .map(|value| value as f32)
                                    .or(scan_mobility),
                            })
                            .collect::<Vec<_>>();
                        if ions.is_empty() && window.target > 0.0 {
                            vec![Precursor {
                                mz: window.target,
                                spectrum_ref,
                                isolation_window,
                                inverse_ion_mobility: scan_mobility,
                                ..Precursor::default()
                            }]
                        } else {
                            ions
                        }
                    })
                    .collect();
                let raw: MzDataSpectrum = spectrum.into();
                let fragment_charges =
                    raw.arrays.charges().ok().and_then(|charges| {
                        normalize_fragment_charges(&charges, raw.mzs().len(), &id)
                    });

                RawSpectrum {
                    file_id: self.file_id,
                    ms_level,
                    id,
                    precursors,
                    representation,
                    scan_start_time,
                    ion_injection_time,
                    total_ion_current,
                    mz: raw.mzs().iter().map(|value| *value as f32).collect(),
                    intensity: raw.intensities().to_vec(),
                    fragment_charges,
                    mobility: None,
                }
            })
            .collect())
    }
}

#[cfg(test)]
mod tests {
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
}
