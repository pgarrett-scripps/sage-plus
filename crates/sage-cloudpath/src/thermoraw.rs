use opentfraw::{iter_spectra, PrecursorInfo, RawFileReader, ScanParams, SpectrumRecord};
use sage_core::{
    mass::Tolerance,
    spectrum::{AcquisitionGroup, Precursor, RawSpectrum, Representation},
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
                    .and_then(|params| trailer_master_scan(&params))
            })
            .collect::<Vec<_>>();
        let levels = trailer_levels(first_scan, &masters);
        let plausible_events = records
            .iter()
            .map(|record| {
                raw.scan_events
                    .get((record.scan_number - first_scan) as usize)
                    .and_then(|event| event.reactions.first())
                    .is_some_and(|reaction| plausible_precursor_mz(reaction.precursor_mz))
            })
            .collect::<Vec<_>>();
        // Once any event contradicts its trailer, the events of the whole file
        // are out of step: their precursors belong to other scans and a
        // plausible-looking masterless MSn event is no evidence of DIA.
        let misaligned = records
            .iter()
            .zip(&plausible_events)
            .any(|(record, &plausible)| {
                let idx = (record.scan_number - first_scan) as usize;
                corrected_level(record.ms_level, levels[idx], plausible).is_some()
            });
        let mut corrected = 0;
        for (record, &plausible) in records.iter_mut().zip(&plausible_events) {
            let idx = (record.scan_number - first_scan) as usize;
            let level = match misaligned {
                true => levels[idx],
                false => corrected_level(record.ms_level, levels[idx], plausible),
            };
            let Some(level) = level else {
                continue;
            };
            if level != record.ms_level {
                corrected += 1;
            }
            let params = raw.scan_params(record.scan_number);
            apply_trailer_level(record, level, params.as_ref());
        }
        let unsearchable = records
            .iter()
            .filter(|record| record.ms_level == 2 && precursor_mz(record).is_none())
            .count();
        if corrected > 0 {
            log::warn!(
                "OpenTFRaw scan events contradict the scan trailers of {} of {} scans in {}; \
                 MS levels and precursors come from the trailers",
                corrected,
                records.len(),
                path.display()
            );
        }
        if unsearchable > 0 {
            log::warn!(
                "{} MS2 scans in {} have no plausible precursor m/z and will not be \
                 searched; convert the file to mzML with msconvert to search them",
                unsearchable,
                path.display()
            );
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
        let acquisition = record
            .filter
            .as_deref()
            .map(AcquisitionGroup::from_thermo_filter)
            .unwrap_or_default();
        let precursor = record.precursor.and_then(|value| {
            // MS3 reporter-ion quantification needs only the link to the MS2
            // scan, which the trailer keeps even when the m/z is unknown.
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
            acquisition,
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

/// The MS level to use when a scan's event contradicts its trailer, or
/// `None` to keep the event.
///
/// OpenTFRaw can decode scan events out of step with the scans on some
/// Orbitrap Fusion files, while the trailer master scan numbers stay correct.
/// A dependent scan (one with a master) takes its level from the master
/// chain. A scan without a master is MS1 unless its event is an MSn scan with
/// a plausible precursor, as DIA and targeted scans have no master scan. Once
/// any scan is corrected, the reader applies the trailer levels to every scan
/// instead, since the file's events can no longer be trusted.
pub(crate) fn corrected_level(
    event_level: u32,
    trailer_level: Option<u32>,
    plausible_event_precursor: bool,
) -> Option<u32> {
    match trailer_level? {
        level if level == event_level => None,
        1 if plausible_event_precursor => None,
        level => Some(level),
    }
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

/// The trailer "Master Scan Number". OpenTFRaw falls back to "Master Index",
/// which is a 0/1 flag on QE and LTQ Orbitrap files rather than a scan number.
pub(crate) fn trailer_master_scan(params: &ScanParams<'_>) -> Option<u32> {
    params
        .record()
        .get_i32("Master Scan Number:")
        .map(|master| master.max(0) as u32)
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
            // `isolation_target_mz` can be an isolation offset, so only
            // plausible m/z values are kept.
            let monoisotopic = params
                .monoisotopic_mz()
                .filter(|&mz| plausible_precursor_mz(mz));
            let target = params
                .isolation_target_mz()
                .filter(|&mz| plausible_precursor_mz(mz))
                .or(monoisotopic);
            PrecursorInfo {
                target_mz: target,
                selected_mz: monoisotopic.or(target),
                isolation_width: params.isolation_width_mz(),
                charge: params.charge_state().filter(|&charge| charge > 0),
                master_scan_number: trailer_master_scan(params).filter(|&scan| scan > 0),
                ..Default::default()
            }
        }),
    };
}

#[cfg(test)]
#[path = "../tests/unit/thermoraw.rs"]
mod tests;
