use sage_core::{
    mass::Tolerance,
    spectrum::{AcquisitionGroup, Precursor, RawSpectrum, Representation},
};
use sage_plus_raw::{try_iter_spectra, PrecursorInfo, RawFileReader, SpectrumRecord};
use std::{fs::File, io::BufReader, path::Path};

/// Reads a local Thermo Fisher `.raw` file directly into Sage spectra.
///
/// sage-plus-raw resolves profile-mode scans to their centroid peak lists
/// when `include_profile` is false, which matches Sage's search input
/// requirements.
pub struct ThermoRawReader {
    file_id: usize,
}

impl ThermoRawReader {
    pub fn with_file_id(file_id: usize) -> Self {
        Self { file_id }
    }

    pub fn parse(&self, path: impl AsRef<Path>) -> sage_plus_raw::Result<Vec<RawSpectrum>> {
        let path = path.as_ref();
        let raw = RawFileReader::open_path(path)?;
        let mut source = BufReader::new(File::open(path)?);

        let expected = raw.num_scans as usize;
        let (records, skipped) =
            decoded_scans(expected, try_iter_spectra(&raw, &mut source, false));
        if let Some((failed, error)) = skipped {
            log::warn!(
                "skipped {} of {} scans in {} that could not be decoded; first error: {}",
                failed,
                expected,
                path.display(),
                error
            );
        }
        let dropped: usize = records.iter().map(|record| record.dropped_peaks).sum();
        if dropped > 0 {
            log::warn!(
                "dropped {} peaks with an invalid m/z or intensity from {}",
                dropped,
                path.display()
            );
        }
        let unsearchable = records
            .iter()
            .filter(|record| record.ms_level == 2 && precursor_mz(record).is_none())
            .count();
        if unsearchable > 0 {
            log::warn!(
                "{} MS2 scans in {} have no plausible precursor m/z and will not be \
                 searched; convert the file to mzML with msconvert to search them",
                unsearchable,
                path.display()
            );
        }
        Ok(records
            .into_iter()
            .map(|record| self.convert(record))
            .collect())
    }

    fn convert(&self, record: SpectrumRecord) -> RawSpectrum {
        let ms_level = record.ms_level;
        let acquisition = record
            .filter
            .as_deref()
            .map(AcquisitionGroup::from_thermo_filter)
            .unwrap_or_default();
        let precursor = record.precursor.and_then(|value| {
            // MS3 reporter-ion quantification (SPS-MS3 TMT) needs only the
            // link to the MS2 scan, which the reader keeps even when the m/z
            // is unknown.
            let mz = match plausible_mz(&value) {
                Some(mz) => mz as f32,
                None if ms_level > 2 && value.master_scan_number.is_some() => 0.0,
                None => return None,
            };
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
            // `try_iter_spectra(..., false)` resolves every scan to its
            // centroid peak list, even when the nominal mode is profile.
            representation: Representation::Centroid,
            scan_start_time: record.retention_time_min as f32,
            ion_injection_time: record.ion_injection_time_ms.unwrap_or_default() as f32,
            total_ion_current: record.total_ion_current as f32,
            mz: record.mz.into_iter().map(|mz| mz as f32).collect(),
            intensity: record.intensity,
            fragment_charges: None,
            mobility: None,
            acquisition,
        }
    }
}

/// The scans that decoded, and the number that did not with the first
/// error. A scan that fails to decode is skipped rather than failing the file.
pub(crate) fn decoded_scans(
    expected: usize,
    results: impl Iterator<Item = sage_plus_raw::Result<SpectrumRecord>>,
) -> (Vec<SpectrumRecord>, Option<(usize, sage_plus_raw::Error)>) {
    let mut records = Vec::with_capacity(expected);
    let mut skipped: Option<(usize, sage_plus_raw::Error)> = None;
    for result in results {
        match result {
            Ok(record) => records.push(record),
            Err(error) => match &mut skipped {
                Some((failed, _)) => *failed += 1,
                None => skipped = Some((1, error)),
            },
        }
    }
    (records, skipped)
}

pub(crate) fn plausible_precursor_mz(mz: f64) -> bool {
    (50.0..20_000.0).contains(&mz)
}

/// The selected m/z, or else the isolation target, when it is plausible.
fn plausible_mz(precursor: &PrecursorInfo) -> Option<f64> {
    [precursor.selected_mz, precursor.target_mz]
        .into_iter()
        .flatten()
        .find(|&mz| plausible_precursor_mz(mz))
}

fn precursor_mz(record: &SpectrumRecord) -> Option<f64> {
    record.precursor.as_ref().and_then(plausible_mz)
}

#[cfg(test)]
#[path = "../tests/unit/thermoraw.rs"]
mod tests;
