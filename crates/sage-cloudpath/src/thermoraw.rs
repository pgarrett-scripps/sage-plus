use sage_core::{
    mass::Tolerance,
    spectrum::{AcquisitionGroup, Precursor, RawSpectrum, Representation},
};
use sage_plus_raw::{PrecursorInfo, SpectrumReader, SpectrumRecord};
use std::path::Path;

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

    /// Read the file's spectra. MS1 scans are decoded only when
    /// `requires_ms1` is set: a search reads MS2 scans, and only
    /// quantification (and DIA) reads MS1, which is the same rule the Bruker
    /// reader already follows.
    pub fn parse(
        &self,
        path: impl AsRef<Path>,
        requires_ms1: bool,
    ) -> sage_plus_raw::Result<Vec<RawSpectrum>> {
        let path = path.as_ref();
        // One handle for the structure and the scans, and one m/z buffer and
        // one block buffer (inside the reader) for the whole file, instead of
        // a second open and three allocations per scan. Each scan is turned
        // into its Sage spectrum as it is read, so the f64 m/z arrays of the
        // whole file are never alive at once.
        let mut reader = SpectrumReader::open_path(path)?;
        let expected = reader.num_scans() as usize;

        let mut spectra = Vec::with_capacity(expected);
        let mut mz = Vec::new();
        let mut intensity = Vec::new();
        let mut dropped = 0usize;
        let mut unsearchable = 0usize;
        let mut skipped: Option<(usize, sage_plus_raw::Error)> = None;
        let mut undecoded_ms1 = 0usize;
        for idx in 0..reader.num_scans() {
            // The MS level and the scan's own metadata are in the scan index
            // and the scan event, so a scan the run will not search is never
            // decoded: its peak block is not read, parsed or centroided, and
            // on a 3 GB DDA file those are most of the bytes. `requires_ms1`
            // is tested first, so a run that reads MS1 does not pay for a
            // metadata lookup per scan.
            if !requires_ms1 {
                let metadata = sage_plus_raw::scan_metadata(reader.raw(), idx);
                if !decodes_scan(metadata.as_ref().map(|meta| meta.ms_level), false) {
                    // The scan still dates the file. Quality control takes a
                    // file's retention-time span, and the gradient bins over
                    // it, from the scans it is given, so the spectrum is kept
                    // with its metadata and without peaks: it says when the
                    // instrument scanned and measures nothing, and nothing
                    // counts it as an MS1 spectrum that was read.
                    if let Some(metadata) = metadata {
                        spectra.push(self.undecoded(&metadata));
                    }
                    undecoded_ms1 += 1;
                    continue;
                }
            }
            match reader.read_spectrum_into(idx, false, &mut mz, &mut intensity) {
                Ok(header) => {
                    dropped += header.dropped_peaks;
                    // The arrays stay in the buffers: `into_record` takes the
                    // metadata, the m/z buffer is read again next scan, and
                    // the intensities move into the spectrum, which holds
                    // exactly the `Vec<f32>` the reader filled.
                    let record = header.into_record(Vec::new(), Vec::new());
                    if record.ms_level == 2 && precursor_mz(&record).is_none() {
                        unsearchable += 1;
                    }
                    spectra.push(self.convert(&record, &mz, std::mem::take(&mut intensity)));
                }
                Err(error) => note_skipped(&mut skipped, error),
            }
        }

        if undecoded_ms1 > 0 {
            log::debug!(
                "{}: {} MS1 scans were not decoded; this search reads only MS2",
                path.display(),
                undecoded_ms1
            );
        }
        if let Some((failed, error)) = skipped {
            log::warn!(
                "skipped {} of {} scans in {} that could not be decoded; first error: {}",
                failed,
                expected,
                path.display(),
                error
            );
        }
        if dropped > 0 {
            log::warn!(
                "dropped {} peaks with an invalid m/z or intensity from {}",
                dropped,
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
        Ok(spectra)
    }

    /// A kept scan whose peaks were not decoded: the same fields
    /// [`Self::convert`] fills from a decoded scan's record, taken from the
    /// scan metadata that decided to skip it, and no peaks. It dates the file
    /// for quality control and measures nothing.
    fn undecoded(&self, metadata: &sage_plus_raw::ScanMetadata) -> RawSpectrum {
        RawSpectrum {
            file_id: self.file_id,
            ms_level: metadata.ms_level as u8,
            id: scan_id(metadata.scan_number),
            precursors: Vec::new(),
            representation: Representation::Centroid,
            scan_start_time: metadata.retention_time_min as f32,
            ion_injection_time: metadata.ion_injection_time_ms.unwrap_or_default() as f32,
            total_ion_current: metadata.total_ion_current as f32,
            mz: Vec::new(),
            intensity: Vec::new(),
            fragment_charges: None,
            mobility: None,
            acquisition: metadata
                .filter
                .as_deref()
                .map(AcquisitionGroup::from_thermo_filter)
                .unwrap_or_default(),
            metadata_only: true,
        }
    }

    /// Build a Sage spectrum from one decoded scan. The peak arrays are
    /// passed apart from the scan's metadata so the read loop can keep its
    /// m/z buffer across scans and hand the intensities over without a copy:
    /// `RawSpectrum` holds exactly the `Vec<f32>` the reader filled.
    fn convert(&self, record: &SpectrumRecord, mz: &[f64], intensity: Vec<f32>) -> RawSpectrum {
        let ms_level = record.ms_level;
        let acquisition = record
            .filter
            .as_deref()
            .map(AcquisitionGroup::from_thermo_filter)
            .unwrap_or_default();
        let precursor = record.precursor.as_ref().and_then(|value| {
            // MS3 reporter-ion quantification (SPS-MS3 TMT) needs only the
            // link to the MS2 scan, which the reader keeps even when the m/z
            // is unknown.
            let mz = match plausible_mz(value) {
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
            id: scan_id(record.scan_number),
            precursors: precursor.into_iter().collect(),
            // Scans are read with `include_profile` false, which resolves
            // every scan to its centroid peak list even when the nominal
            // mode is profile.
            representation: Representation::Centroid,
            scan_start_time: record.retention_time_min as f32,
            ion_injection_time: record.ion_injection_time_ms.unwrap_or_default() as f32,
            total_ion_current: record.total_ion_current as f32,
            mz: mz.iter().map(|&mz| mz as f32).collect(),
            intensity,
            fragment_charges: None,
            mobility: None,
            acquisition,
            metadata_only: false,
        }
    }
}

/// Record a scan that could not be decoded: how many failed, and the first
/// error. A scan that fails is skipped rather than failing the whole file.
/// A scan's Thermo native identifier.
fn scan_id(scan_number: u32) -> String {
    format!("controllerType=0 controllerNumber=1 scan={scan_number}")
}

/// Whether the scan at an index must be decoded, from its MS level alone.
///
/// A run that reads MS1 decodes every scan. Otherwise an MS1 scan is passed
/// over, and a scan whose metadata could not be read is decoded, so nothing is
/// dropped on a guess.
fn decodes_scan(level: Option<u32>, requires_ms1: bool) -> bool {
    requires_ms1 || level != Some(1)
}

pub(crate) fn note_skipped(
    skipped: &mut Option<(usize, sage_plus_raw::Error)>,
    error: sage_plus_raw::Error,
) {
    match skipped {
        Some((failed, _)) => *failed += 1,
        None => *skipped = Some((1, error)),
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

#[cfg(test)]
#[path = "../tests/unit/thermoraw.rs"]
mod tests;
