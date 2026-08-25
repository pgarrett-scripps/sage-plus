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
#[path = "../tests/unit/thermoraw.rs"]
mod tests;
