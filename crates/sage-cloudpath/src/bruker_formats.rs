//! Bruker TSF and ProteoScape miniTDF spectra, read through `sage-plus-tdf`.
//!
//! TSF (`analysis.tsf` + `analysis.tsf_bin`, one centroided spectrum per frame
//! and no ion mobility) is read with its calibrated m/z model; its MS/MS frames
//! take their precursor from `FrameMsMsInfo`. timsrust 0.6 gave TSF MS/MS
//! spectra no precursor, so Beta 16 searched none of them.
//!
//! miniTDF (`<name>.ms2spectrum.bin` + `.parquet`) holds centroided MS2 spectra
//! with their precursors; it is read as timsrust 0.6 did, with the m/z values
//! and isolation widths that reader reported.
use crate::tdf::TIMS_TOF;
use sage_core::{
    mass::Tolerance,
    spectrum::{Precursor, RawSpectrum, Representation},
};
use sage_plus_tdf::{MiniTdfPrecursor, MiniTdfReader, TsfReader};
use std::{
    collections::HashMap,
    path::{Path, PathBuf},
};

/// The `.d` directory holding `analysis.tsf` and `analysis.tsf_bin`, from the
/// directory itself or a file inside it.
pub fn tsf_directory(path: &Path) -> Option<PathBuf> {
    path.ancestors()
        .filter(|directory| !directory.as_os_str().is_empty())
        .find(|directory| {
            directory.join("analysis.tsf").is_file() && directory.join("analysis.tsf_bin").is_file()
        })
        .map(Path::to_path_buf)
}

/// The directory holding miniTDF `ms2spectrum` files, from the directory
/// itself or one of the files.
pub fn minitdf_directory(path: &Path) -> Result<Option<PathBuf>, crate::Error> {
    let directory = if path.is_file() {
        match path.parent() {
            Some(parent) => parent,
            None => return Ok(None),
        }
    } else {
        path
    };
    Ok(MiniTdfReader::find(directory)?.map(|_| directory.to_path_buf()))
}

/// MS2 spectra of a TSF acquisition (plus MS1 when `requires_ms1`), with
/// calibrated m/z. The id is the frame `Id`.
pub fn read_tsf(
    path: &Path,
    directory: &Path,
    file_id: usize,
    requires_ms1: bool,
) -> Result<Vec<RawSpectrum>, crate::Error> {
    let start = std::time::Instant::now();
    let reader = TsfReader::open(directory)?;
    if !reader.has_line_spectra() {
        return Err(crate::Error::Unsupported(format!(
            "{}: TSF without line (centroid) spectra; profile spectra are not read",
            path.display()
        )));
    }
    let msms: HashMap<u64, sage_plus_tdf::TsfFrameMsMsInfo> = reader
        .frame_msms_info()?
        .into_iter()
        .map(|row| (row.frame, row))
        .collect();
    let mut spectra = Vec::new();
    let mut missing_precursor = 0usize;
    for frame in reader.frames() {
        let level = frame.ms_level();
        if level == 1 && !requires_ms1 {
            continue;
        }
        let precursor = if level == 1 {
            None
        } else {
            match msms.get(&frame.id) {
                Some(row) => Some(row),
                None => {
                    missing_precursor += 1;
                    continue;
                }
            }
        };
        let model = reader.mz_model(frame.id)?;
        let line = reader.read_line_spectrum(frame.id)?;
        let mz = line
            .indices
            .iter()
            .map(|&index| model.mz(index).map(|mz| mz as f32))
            .collect::<Result<Vec<f32>, _>>()?;
        let rt = frame.retention_time_seconds;
        let spectrum = match precursor {
            None => RawSpectrum {
                file_id,
                precursors: vec![],
                representation: Representation::Centroid,
                scan_start_time: rt as f32 / 60.0,
                // As for TDF MS1 frames, which record no injection time either.
                ion_injection_time: 100.0,
                total_ion_current: line.intensities.iter().sum(),
                mz,
                ms_level: 1,
                id: frame.id.to_string(),
                intensity: line.intensities,
                fragment_charges: None,
                mobility: None,
                acquisition: TIMS_TOF,
            },
            Some(row) => {
                let half = row.isolation_width as f32 / 2.0;
                let precursor = Precursor {
                    mz: row.trigger_mass as f32,
                    charge: row
                        .precursor_charge
                        .and_then(|charge| u8::try_from(charge).ok())
                        .filter(|&charge| charge > 0 && charge <= i8::MAX as u8),
                    spectrum_ref: row.parent.map(|parent| parent.to_string()),
                    isolation_window: Some(Tolerance::Da(-half, half)),
                    ..Precursor::default()
                };
                ms2_spectrum(
                    file_id,
                    frame.id.to_string(),
                    precursor,
                    rt,
                    mz,
                    line.intensities,
                )
            }
        };
        spectra.push(spectrum);
    }
    if missing_precursor > 0 {
        log::warn!(
            "{}: skipped {missing_precursor} TSF MS/MS frames without a FrameMsMsInfo row",
            path.display()
        );
    }
    log::info!(
        "read {} TSF spectra in {:#?}",
        spectra.len(),
        start.elapsed()
    );
    Ok(spectra)
}

/// miniTDF MS2 spectra, one per precursor row. The id is the row index, as in
/// timsrust 0.6. MS1 frames (`msframe.*`) are not read.
pub fn read_minitdf(
    path: &Path,
    directory: &Path,
    file_id: usize,
    requires_ms1: bool,
) -> Result<Vec<RawSpectrum>, crate::Error> {
    if requires_ms1 {
        return Err(crate::Error::Unsupported(format!(
            "{}: miniTDF MS1 frames are not read, so LFQ is unavailable for miniTDF input",
            path.display()
        )));
    }
    let start = std::time::Instant::now();
    let reader = MiniTdfReader::open(directory)?;
    let mut spectra = Vec::with_capacity(reader.len());
    let mut dropped = 0usize;
    for (index, row) in reader.precursors().iter().enumerate() {
        let Some(mz) = row.monoisotopic_mz else {
            dropped += 1;
            continue;
        };
        let spectrum = reader.read_spectrum(index)?;
        let half = MiniTdfPrecursor::assumed_isolation_width(mz) as f32 / 2.0;
        let precursor = Precursor {
            mz: mz as f32,
            charge: row
                .charge
                .and_then(|charge| u8::try_from(charge).ok())
                .filter(|&charge| charge > 0 && charge <= i8::MAX as u8),
            intensity: row.intensity.map(|value| value as f32),
            spectrum_ref: row.parent_frame.map(|frame| frame.to_string()),
            inverse_ion_mobility: row.one_over_k0.map(|value| value as f32),
            isolation_window: Some(Tolerance::Da(-half, half)),
        };
        spectra.push(ms2_spectrum(
            file_id,
            index.to_string(),
            precursor,
            row.retention_time_seconds,
            spectrum.mz.iter().map(|&mz| mz as f32).collect(),
            spectrum.intensities,
        ));
    }
    if dropped > 0 {
        log::warn!(
            "{}: skipped {dropped} miniTDF spectra without a monoisotopic m/z",
            path.display()
        );
    }
    log::info!(
        "read {} miniTDF spectra in {:#?}",
        spectra.len(),
        start.elapsed()
    );
    Ok(spectra)
}

/// An MS2 spectrum with the conventions of the TDF reader: the retention time
/// in minutes as the scan start and in seconds as the injection time.
fn ms2_spectrum(
    file_id: usize,
    id: String,
    precursor: Precursor,
    rt_seconds: f64,
    mz: Vec<f32>,
    intensity: Vec<f32>,
) -> RawSpectrum {
    RawSpectrum {
        file_id,
        precursors: vec![precursor],
        representation: Representation::Centroid,
        scan_start_time: rt_seconds as f32 / 60.0,
        ion_injection_time: rt_seconds as f32,
        total_ion_current: 0.0,
        mz,
        ms_level: 2,
        id,
        intensity,
        fragment_charges: None,
        mobility: None,
        acquisition: TIMS_TOF,
    }
}

#[cfg(test)]
#[path = "../tests/unit/bruker_formats.rs"]
mod tests;
