use opentfraw::{iter_spectra, PrecursorInfo, RawFileReader, ScanParams, SpectrumRecord};
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
        let mut records: Vec<_> = iter_spectra(&raw, &mut source, false).collect();
        let first_scan = raw.run_header.sample_info.first_scan_number;
        let masters = (0..raw.num_scans)
            .map(|idx| {
                raw.scan_params(first_scan + idx)
                    .and_then(|params| params.master_scan_number())
                    .map(|master| master.max(0) as u32)
            })
            .collect::<Vec<_>>();
        let levels = trailer_levels(first_scan, &masters);
        let disagreements = event_level_disagreements(first_scan, &records, &levels);
        if events_misaligned(disagreements, records.len()) {
            log::warn!(
                "OpenTFRaw scan events disagree with the master scan numbers of {} of {} scans in {}; \
                 deriving MS levels and precursors from the scan trailers",
                disagreements,
                records.len(),
                path.display()
            );
            for record in &mut records {
                let idx = (record.scan_number - first_scan) as usize;
                if let Some(level) = levels[idx] {
                    let params = raw.scan_params(record.scan_number);
                    apply_trailer_level(record, level, params.as_ref());
                }
            }
            let unsearchable = records
                .iter()
                .filter(|record| {
                    record.ms_level == 2
                        && record
                            .precursor
                            .as_ref()
                            .and_then(|precursor| precursor.selected_mz.or(precursor.target_mz))
                            .is_none()
                })
                .count();
            if unsearchable > 0 {
                log::warn!(
                    "{} MS2 scans in {} have no precursor m/z in their trailer and will not be \
                     searched; convert the file to mzML with msconvert to search them",
                    unsearchable,
                    path.display()
                );
            }
        }
        let spectra: Vec<_> = records
            .into_iter()
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
        let ms_level = record.ms_level;
        let precursor = record.precursor.and_then(|value| {
            // MS3 reporter-ion quantification needs only the link to the MS2
            // scan, which the trailer keeps even when the m/z is unknown.
            let mz = match value.selected_mz.or(value.target_mz) {
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
            // `iter_spectra(..., false)` resolves every scan to its centroid
            // peak list, even when the instrument's nominal mode is profile.
            representation: Representation::Centroid,
            scan_start_time: record.retention_time_min as f32,
            ion_injection_time: record.ion_injection_time_ms.unwrap_or_default() as f32,
            total_ion_current: record.total_ion_current as f32,
            mz: record.mz.into_iter().map(|mz| mz as f32).collect(),
            intensity: record.intensity,
            fragment_charges: None,
            mobility: None,
        }
    }
}

/// MS levels implied by the trailer master scan numbers, indexed from the
/// first scan. A scan without a master is MS1, and a dependent scan is one
/// level above its master, so SPS-MS3 scans resolve to MS3. Scans without a
/// trailer master, or whose master is not an earlier scan, have no level.
pub(crate) fn trailer_levels(first_scan: u32, masters: &[Option<u32>]) -> Vec<Option<u32>> {
    let mut levels: Vec<Option<u32>> = Vec::with_capacity(masters.len());
    for (idx, master) in masters.iter().enumerate() {
        let scan = first_scan + idx as u32;
        let level = match *master {
            None => None,
            Some(0) => Some(1),
            Some(master) if master >= first_scan && master < scan => {
                levels[(master - first_scan) as usize].map(|level| level + 1)
            }
            Some(_) => None,
        };
        levels.push(level);
    }
    levels
}

/// Decoded scans whose event MS level differs from the trailer level.
pub(crate) fn event_level_disagreements(
    first_scan: u32,
    records: &[SpectrumRecord],
    levels: &[Option<u32>],
) -> usize {
    records
        .iter()
        .filter(|record| {
            levels
                .get((record.scan_number - first_scan) as usize)
                .copied()
                .flatten()
                .is_some_and(|level| level != record.ms_level)
        })
        .count()
}

/// OpenTFRaw can decode scan events out of step with the scans on some
/// Orbitrap Fusion files, so that about half of the dependent scans read as
/// MS1. A few disagreements are tolerated; beyond that the trailers win.
pub(crate) fn events_misaligned(disagreements: usize, scans: usize) -> bool {
    disagreements > 10.max(scans / 100)
}

/// Replace the event-derived MS level and precursor of `record` with values
/// from its scan trailer, ignoring the event's reaction entirely.
pub(crate) fn apply_trailer_level(
    record: &mut SpectrumRecord,
    level: u32,
    params: Option<&ScanParams<'_>>,
) {
    record.ms_level = level;
    record.is_ms1 = level == 1;
    record.precursor = match level {
        1 => None,
        _ => params.map(|params| {
            let monoisotopic = params.monoisotopic_mz().filter(|&mz| mz > 0.0);
            let target = params
                .isolation_target_mz()
                .filter(|&mz| mz > 0.0)
                .or(monoisotopic);
            PrecursorInfo {
                target_mz: target,
                selected_mz: monoisotopic.or(target),
                isolation_width: params.isolation_width_mz(),
                charge: params.charge_state().filter(|&charge| charge > 0),
                master_scan_number: params
                    .master_scan_number()
                    .filter(|&scan| scan > 0)
                    .map(|scan| scan as u32),
                ..Default::default()
            }
        }),
    };
}

#[cfg(test)]
#[path = "../tests/unit/thermoraw.rs"]
mod tests;
