use opentfraw::{iter_spectra, RawFileReader, SpectrumRecord};
use sage_core::{
    mass::Tolerance,
    spectrum::{Precursor, RawSpectrum, Representation},
};
use std::{fs::File, io::BufReader, path::Path};

/// Reads a local Thermo Fisher `.raw` file directly into Sage spectra.
///
/// OpenTFRaw resolves profile-mode scans to their centroid peak lists when
/// `include_profile` is false, which matches Sage's search input requirements.
pub struct ThermoRawReader {
    file_id: usize,
}

impl ThermoRawReader {
    pub fn with_file_id(file_id: usize) -> Self {
        Self { file_id }
    }

    pub fn parse(&self, path: impl AsRef<Path>) -> opentfraw::Result<Vec<RawSpectrum>> {
        let path = path.as_ref();
        let raw = RawFileReader::open_path(path)?;
        let mut source = BufReader::new(File::open(path)?);

        let expected = raw.num_scans as usize;
        let spectra: Vec<_> = iter_spectra(&raw, &mut source, false)
            .map(|record| self.convert(record))
            .collect();
        if spectra.len() != expected {
            log::warn!(
                "OpenTFRaw decoded {} of {} scans from {}",
                spectra.len(),
                expected,
                path.display()
            );
        }
        Ok(spectra)
    }

    fn convert(&self, record: SpectrumRecord) -> RawSpectrum {
        let precursor = record.precursor.and_then(|value| {
            let mz = value.selected_mz.or(value.target_mz)? as f32;
            let isolation_window = value.isolation_width.map(|width| {
                let half_width = width as f32 / 2.0;
                Tolerance::Da(-half_width, half_width)
            });

            Some(Precursor {
                mz,
                charge: value.charge.and_then(|charge| u8::try_from(charge).ok()),
                spectrum_ref: value
                    .master_scan_number
                    .map(|scan| format!("controllerType=0 controllerNumber=1 scan={scan}")),
                isolation_window,
                ..Default::default()
            })
        });

        RawSpectrum {
            file_id: self.file_id,
            ms_level: record.ms_level as u8,
            id: format!(
                "controllerType=0 controllerNumber=1 scan={}",
                record.scan_number
            ),
            precursors: precursor.into_iter().collect(),
            // `iter_spectra(..., false)` resolves every scan to its centroid
            // peak list, even when the instrument's nominal mode is profile.
            representation: Representation::Centroid,
            scan_start_time: record.retention_time_min as f32,
            ion_injection_time: record.ion_injection_time_ms.unwrap_or_default() as f32,
            total_ion_current: record.total_ion_current as f32,
            mz: record.mz.into_iter().map(|mz| mz as f32).collect(),
            intensity: record.intensity,
            mobility: None,
        }
    }
}

#[cfg(test)]
mod tests {
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
    #[ignore = "requires SAGE_THERMO_RAW_TEST_FILE to point to a real Thermo RAW file"]
    fn parses_real_raw_file() {
        let path = std::env::var("SAGE_THERMO_RAW_TEST_FILE")
            .expect("set SAGE_THERMO_RAW_TEST_FILE to a local .raw file");
        let spectra = ThermoRawReader::with_file_id(3).parse(path).unwrap();

        assert!(!spectra.is_empty());
        assert!(spectra.iter().all(|spectrum| {
            spectrum.file_id == 3 && spectrum.mz.len() == spectrum.intensity.len()
        }));
        assert!(spectra.iter().any(|spectrum| spectrum.ms_level > 1));
    }
}
