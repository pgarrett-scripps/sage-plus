use crate::denoise::{BrukerDenoiseConfig, DenoiseError, Ms1Denoiser};
use crate::tdf_spectra::{self, DiaSlice, Expansion, IsolationWindow};
use crate::tims_mobility::{
    analysis_tdf, BrukerMobilityScale, LinearMobilityScale, MobilityCalibration,
};
use rayon::prelude::*;
use sage_core::{
    mass::Tolerance,
    spectrum::{
        AcquisitionGroup, Activation, MassAnalyzer, Precursor, RawSpectrum, Representation,
    },
};
use sage_plus_tdf::{AcquisitionType, Frame, LinearMzScale};
use serde::{Deserialize, Serialize};
use std::sync::atomic::{AtomicU64, Ordering as AtomicOrdering};
use std::sync::Mutex;
use std::{cmp::Ordering, collections::HashMap, path::Path};

/// timsTOF spectra are all measured by the TOF analyzer after CID.
pub(crate) const TIMS_TOF: AcquisitionGroup = AcquisitionGroup {
    analyzer: MassAnalyzer::Tof,
    activation: Activation::Cid,
};

pub struct TdfReader;

#[derive(Deserialize, Serialize, Debug, Clone, Copy, schemars::JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct BrukerSpectrumProcessingConfig {
    pub smoothing_window: u32,
    pub centroiding_window: u32,
    pub calibration_tolerance: f64,
    pub calibrate: bool,
}

impl Default for BrukerSpectrumProcessingConfig {
    fn default() -> Self {
        Self {
            smoothing_window: 1,
            centroiding_window: 1,
            calibration_tolerance: 0.1,
            calibrate: false,
        }
    }
}

#[derive(Deserialize, Serialize, Debug, Clone, Copy, schemars::JsonSchema)]
pub enum BrukerQuadWindowExpansionStrategy {
    None,
    Even(usize),
    UniformMobility((f64, f64), Option<()>),
    UniformScan((usize, usize)),
}

impl Default for BrukerQuadWindowExpansionStrategy {
    fn default() -> Self {
        Self::Even(1)
    }
}

#[derive(Deserialize, Serialize, Debug, Clone, Copy, schemars::JsonSchema)]
pub enum BrukerFrameWindowSplittingConfig {
    Quadrupole(BrukerQuadWindowExpansionStrategy),
    Window(BrukerQuadWindowExpansionStrategy),
}

impl Default for BrukerFrameWindowSplittingConfig {
    fn default() -> Self {
        Self::Quadrupole(BrukerQuadWindowExpansionStrategy::default())
    }
}

#[derive(Default, Deserialize, Serialize, Debug, Clone, Copy, schemars::JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct BrukerSpectrumConfig {
    pub spectrum_processing_params: BrukerSpectrumProcessingConfig,
    pub frame_splitting_params: BrukerFrameWindowSplittingConfig,
}

#[derive(Deserialize, Serialize, Debug, Clone, Copy, schemars::JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct BrukerMS1CentoidingConfig {
    pub mz_ppm: f32,
    pub ims_pct: f32,
}

impl Default for BrukerMS1CentoidingConfig {
    fn default() -> Self {
        BrukerMS1CentoidingConfig {
            mz_ppm: 5.0,
            ims_pct: 3.0,
        }
    }
}

#[derive(Default, Deserialize, Serialize, Debug, Clone, Copy, schemars::JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct BrukerProcessingConfig {
    #[serde(default)]
    pub ms2: BrukerSpectrumConfig,
    #[serde(default)]
    pub ms1: BrukerMS1CentoidingConfig,
    /// Scan-to-1/K0 scale for reported ion mobilities. `calibrated` applies the
    /// acquisition's TimsCalibration model; `linear` uses timsrust's interpolation
    /// between the acquisition range limits, as in Beta 6 and earlier.
    #[serde(default)]
    pub ion_mobility_scale: BrukerMobilityScale,
    /// Opt-in dnoise denoising of MS1 frames before centroiding. Off by default.
    #[serde(default)]
    pub denoise: BrukerDenoiseConfig,
}

impl BrukerSpectrumConfig {
    /// diaPASEF splitting: `(per isolation window, expansion)`.
    fn splitting(self) -> (bool, Expansion) {
        let expansion = |strategy| match strategy {
            BrukerQuadWindowExpansionStrategy::None => Expansion::None,
            BrukerQuadWindowExpansionStrategy::Even(count) => Expansion::Even(count),
            BrukerQuadWindowExpansionStrategy::UniformMobility((span, step), _) => {
                Expansion::UniformMobility { span, step }
            }
            BrukerQuadWindowExpansionStrategy::UniformScan((span, step)) => {
                Expansion::UniformScan { span, step }
            }
        };
        match self.frame_splitting_params {
            BrukerFrameWindowSplittingConfig::Quadrupole(strategy) => (true, expansion(strategy)),
            BrukerFrameWindowSplittingConfig::Window(strategy) => (false, expansion(strategy)),
        }
    }
}

/// Scan-to-1/K0 conversion for one acquisition.
enum MobilityScale {
    Calibrated {
        calibration: MobilityCalibration,
        precursors: HashMap<usize, (usize, f64)>,
    },
    /// The scale of Beta 6 and earlier.
    Linear(LinearMobilityScale),
}

impl MobilityScale {
    fn new(
        path: &Path,
        scale: BrukerMobilityScale,
    ) -> Result<Self, crate::tims_mobility::MobilityCalibrationError> {
        Ok(match scale {
            BrukerMobilityScale::Calibrated => Self::Calibrated {
                calibration: MobilityCalibration::from_path(path)?,
                precursors: MobilityCalibration::dda_precursor_scans(path)?,
            },
            BrukerMobilityScale::Linear => Self::Linear(LinearMobilityScale::from_path(path)?),
        })
    }

    /// Precursor 1/K0. DDA precursors convert the fractional average scan with
    /// their parent frame's model; DIA window centers use the run's dominant model.
    /// `scan` is the truncated scan timsrust reported.
    fn precursor(&self, index: usize, frame: usize, scan: u32) -> f32 {
        match self {
            Self::Linear(linear) => linear.one_over_k0(scan) as f32,
            Self::Calibrated {
                calibration,
                precursors,
            } => match precursors.get(&index) {
                Some(&(parent, scan)) if parent == frame => {
                    calibration.frame(parent).one_over_k0(scan) as f32
                }
                _ => calibration.dominant().one_over_k0(f64::from(scan)) as f32,
            },
        }
    }
}

/// The scan range `begin..end` of `frame`, or an error if the frame has fewer scans.
fn scan_slice(frame: &Frame, begin: usize, end: usize) -> Result<(&[u32], &[u32]), String> {
    let (Some(&start), Some(&stop)) = (frame.scan_offsets.get(begin), frame.scan_offsets.get(end))
    else {
        return Err(format!(
            "frame {} has {} scans, scans {begin}..{end} were requested",
            frame.id,
            frame.scan_count()
        ));
    };
    if start > stop {
        return Err(format!(
            "frame {}: scans {begin}..{end} are reversed",
            frame.id
        ));
    }
    Ok((
        &frame.tof_indices[start..stop],
        &frame.intensities[start..stop],
    ))
}

/// Count of MS2 spectra that could not be read, with the first reason, logged once.
#[derive(Default)]
struct Dropped {
    count: AtomicU64,
    first: Mutex<Option<String>>,
}

impl Dropped {
    fn drop_spectrum(&self, why: String) {
        self.count.fetch_add(1, AtomicOrdering::Relaxed);
        self.first.lock().unwrap().get_or_insert(why);
    }

    fn log(self, path: &Path, what: &str) {
        let count = self.count.into_inner();
        if count > 0 {
            log::warn!(
                "{}: skipped {count} {what} spectra that could not be read, first: {}",
                path.display(),
                self.first.into_inner().unwrap().unwrap_or_default()
            );
        }
    }
}

impl TdfReader {
    pub fn parse(
        &self,
        path_name: impl AsRef<Path>,
        file_id: usize,
        config: BrukerProcessingConfig,
        requires_ms1: bool,
    ) -> Result<Vec<RawSpectrum>, crate::Error> {
        let path = path_name.as_ref();
        // timsrust's detection order: TDF, then TSF, then miniTDF.
        let Some(tdf) = analysis_tdf(path) else {
            if let Some(directory) = crate::bruker_formats::tsf_directory(path) {
                return crate::bruker_formats::read_tsf(path, &directory, file_id, requires_ms1);
            }
            if let Some(directory) = crate::bruker_formats::minitdf_directory(path)? {
                return crate::bruker_formats::read_minitdf(
                    path,
                    &directory,
                    file_id,
                    requires_ms1,
                );
            }
            return Err(crate::Error::Unsupported(format!(
                "{}: not a Bruker TDF, TSF or miniTDF acquisition (a .d directory with \
                 analysis.tdf and analysis.tdf_bin, or analysis.tsf and analysis.tsf_bin, or a \
                 directory with ms2spectrum.bin and ms2spectrum.parquet); timsrust parquet \
                 spectra are not supported",
                path.display()
            )));
        };
        let directory = tdf.parent().expect("analysis.tdf lives in a .d directory");
        let reader = sage_plus_tdf::TdfReader::open(directory)?;
        let scale = MobilityScale::new(path, config.ion_mobility_scale)?;
        let mz_scale = reader.linear_mz_scale()?;
        let mut spectra = match reader.acquisition_type() {
            AcquisitionType::DdaPasef => {
                self.read_dda_spectra(path, file_id, &reader, &mz_scale, &scale, config.ms2)?
            }
            AcquisitionType::DiaPasef => {
                self.read_dia_spectra(path, file_id, &reader, &mz_scale, &scale, config.ms2)?
            }
            AcquisitionType::Unknown => {
                return Err(crate::Error::Unsupported(format!(
                    "{}: the acquisition is neither ddaPASEF nor diaPASEF",
                    path.display()
                )))
            }
        };
        if requires_ms1 {
            let ms1s = self.read_ms1_spectra(
                path,
                file_id,
                &reader,
                &mz_scale,
                config.ms1,
                config.denoise,
                &scale,
            )?;
            spectra.extend(ms1s);
        }

        Ok(spectra)
    }

    #[allow(clippy::too_many_arguments)]
    fn read_ms1_spectra(
        &self,
        path: &Path,
        file_id: usize,
        reader: &sage_plus_tdf::TdfReader,
        mz_scale: &LinearMzScale,
        config: BrukerMS1CentoidingConfig,
        denoise: BrukerDenoiseConfig,
        scale: &MobilityScale,
    ) -> Result<Vec<RawSpectrum>, crate::Error> {
        let start = std::time::Instant::now();
        // The denoiser borrows its parameters, so both live for the whole read.
        let denoise_params = denoise.enabled.then(|| denoise.params());
        let denoiser = match &denoise_params {
            Some(params) => Some(params.denoiser(&reader.directory().join("analysis.tdf"))?),
            None => None,
        };
        let raw_points = AtomicU64::new(0);
        let kept_points = AtomicU64::new(0);
        let tol_ppm = config.mz_ppm;
        let im_tol_pct = config.ims_pct;

        let frames: Vec<(u64, f64, usize)> = reader
            .frames()
            .filter(|info| info.ms_level() == Some(1))
            .map(|info| (info.id, info.retention_time_seconds, info.scan_count))
            .collect();
        let ms1_spectra: Vec<RawSpectrum> = frames
            .par_iter()
            .map_init(
                || PeakBuffer::with_capacity(2 * MAX_PEAKS),
                |buffer, &(id, rt, num_scans)| match reader.read_frame(id) {
                    Ok(frame) => {
                        buffer.clear();
                        let survivors = match &denoiser {
                            Some(denoiser) => {
                                let survivors = denoise_frame(denoiser, &frame, num_scans)?;
                                raw_points.fetch_add(
                                    frame.intensities.len() as u64,
                                    AtomicOrdering::Relaxed,
                                );
                                kept_points
                                    .fetch_add(survivors.len() as u64, AtomicOrdering::Relaxed);
                                Some(survivors)
                            }
                            None => None,
                        };
                        let survivors = survivors.as_deref();
                        // Convert every scan before centroiding so the reported
                        // mobility is an average of calibrated values.
                        match scale {
                            MobilityScale::Calibrated { calibration, .. } => {
                                let model = calibration.frame(id as usize);
                                buffer.load(
                                    &frame,
                                    survivors,
                                    |scan| model.one_over_k0(scan as f64) as f32,
                                    mz_scale,
                                )
                            }
                            MobilityScale::Linear(linear) => buffer.load(
                                &frame,
                                survivors,
                                |scan| {
                                    let scan =
                                        u32::try_from(scan).expect("scan index exceeds u32 range");
                                    linear.one_over_k0(scan) as f32
                                },
                                mz_scale,
                            ),
                        }

                        // Squash the mobility dimension
                        let (mz, (intensity, mobility)): (Vec<f32>, (Vec<f32>, Vec<f32>)) =
                            buffer.fastcentroid_frame(tol_ppm, im_tol_pct);

                        let scan_start_time = rt as f32 / 60.0;
                        let ion_injection_time = 100.0; // This is made up, in theory we can read
                                                        // if from the tdf file
                        let total_ion_current = intensity.iter().sum::<f32>();

                        let spec = RawSpectrum {
                            file_id,
                            precursors: vec![],
                            representation: Representation::Centroid,
                            scan_start_time,
                            ion_injection_time,
                            mz,
                            ms_level: 1,
                            id: id.to_string(),
                            intensity,
                            total_ion_current,
                            fragment_charges: None,
                            mobility: Some(mobility),
                            acquisition: TIMS_TOF,
                        };
                        Ok(Some(spec))
                    }
                    Err(x) => {
                        log::error!("error parsing spectrum: {x}");
                        Ok(None)
                    }
                },
            )
            .filter_map(Result::transpose)
            .collect::<Result<_, DenoiseError>>()?;
        log::info!(
            "read {} ms1 spectra in {:#?}",
            ms1_spectra.len(),
            start.elapsed()
        );
        if denoiser.is_some() {
            let raw = raw_points.into_inner();
            let kept = kept_points.into_inner();
            log::info!(
                "{}: denoising kept {kept} of {raw} MS1 points ({:.1}%)",
                path.display(),
                100.0 * kept as f64 / raw.max(1) as f64
            );
        }
        Ok(ms1_spectra)
    }

    /// ddaPASEF spectra, one per precursor, in precursor `Id` order. The id is
    /// the precursor `Id`.
    fn read_dda_spectra(
        &self,
        path: &Path,
        file_id: usize,
        reader: &sage_plus_tdf::TdfReader,
        mz_scale: &LinearMzScale,
        scale: &MobilityScale,
        ms2: BrukerSpectrumConfig,
    ) -> Result<Vec<RawSpectrum>, crate::Error> {
        let start = std::time::Instant::now();
        let rows = reader.pasef_frame_msms_info()?;
        // timsrust paired the n-th spectrum with the n-th Precursors row; pairing by
        // Id is the same when every precursor has PASEF rows, and correct otherwise.
        let precursors: HashMap<u64, sage_plus_tdf::Precursor> = reader
            .precursors()?
            .into_iter()
            .map(|precursor| (precursor.id, precursor))
            .collect();
        let processing = ms2.spectrum_processing_params;
        let dropped = Dropped::default();
        let spectra: Vec<RawSpectrum> = tdf_spectra::dda_groups(&rows)
            .par_iter()
            .filter_map(|(id, indices)| {
                let spectrum = (|| {
                    let precursor = precursors
                        .get(id)
                        .ok_or_else(|| format!("precursor {id} is not in Precursors"))?;
                    let parent = precursor
                        .parent
                        .ok_or_else(|| format!("precursor {id} has no Parent frame"))?;
                    let scan = precursor
                        .scan_number
                        .ok_or_else(|| format!("precursor {id} has no ScanNumber"))?;
                    let rt = reader
                        .frame(parent)
                        .ok_or_else(|| format!("precursor {id}: unknown Parent frame {parent}"))?
                        .retention_time_seconds;
                    let mut tofs = Vec::new();
                    let mut intensities = Vec::new();
                    let mut last = &rows[indices[0]];
                    for &index in indices {
                        let row = &rows[index];
                        last = row;
                        let frame = reader.read_frame(row.frame).map_err(|e| e.to_string())?;
                        if frame.intensities.is_empty() {
                            continue;
                        }
                        let (t, i) = scan_slice(
                            &frame,
                            row.scan_num_begin as usize,
                            row.scan_num_end as usize,
                        )?;
                        tofs.extend_from_slice(t);
                        intensities.extend(i.iter().map(|&x| u64::from(x)));
                    }
                    let (tofs, intensities) = tdf_spectra::process(
                        tofs,
                        intensities,
                        processing.smoothing_window,
                        processing.centroiding_window,
                    );
                    let scan = scan as u32;
                    let inverse_ion_mobility = scale.precursor(*id as usize, parent as usize, scan);
                    let precursor = Precursor {
                        mz: precursor.monoisotopic_mz.unwrap_or_default() as f32,
                        charge: precursor
                            .charge
                            .and_then(|charge| i8::try_from(charge).ok())
                            .map(|charge| charge as u8),
                        intensity: Some(precursor.intensity as f32),
                        spectrum_ref: Some(parent.to_string()),
                        inverse_ion_mobility: Some(inverse_ion_mobility),
                        ..Precursor::default()
                    };
                    let window = IsolationWindow::from_center(
                        last.isolation_mz,
                        last.isolation_width,
                        last.collision_energy,
                    );
                    Ok::<_, String>(Self::ms2_spectrum(
                        file_id,
                        id.to_string(),
                        precursor,
                        rt,
                        window,
                        &tofs,
                        &intensities,
                        mz_scale,
                    ))
                })();
                spectrum.map_err(|why| dropped.drop_spectrum(why)).ok()
            })
            .collect();
        dropped.log(path, "ddaPASEF");
        log::info!(
            "read {} ddaPASEF spectra in {:#?}",
            spectra.len(),
            start.elapsed()
        );
        Ok(spectra)
    }

    /// diaPASEF MS2 spectra, one per window split, with the window center as the precursor
    /// m/z and the full window as the isolation window. The id is the spectrum index.
    ///
    /// timsrust 0.6's `SpectrumReader` sent diaPASEF to its precursor-anchored centroid
    /// reader, which ignores `bruker_config.ms2`; Sage used its TDF window reader
    /// instead, reproduced here: one spectrum per diaPASEF window, split by
    /// `ms2.frame_splitting_params`.
    fn read_dia_spectra(
        &self,
        path: &Path,
        file_id: usize,
        reader: &sage_plus_tdf::TdfReader,
        mz_scale: &LinearMzScale,
        scale: &MobilityScale,
        ms2: BrukerSpectrumConfig,
    ) -> Result<Vec<RawSpectrum>, crate::Error> {
        let start = std::time::Instant::now();
        let failed =
            |why: String| crate::Error::Unsupported(format!("{}: diaPASEF {why}", path.display()));
        let groups = tdf_spectra::quad_groups(&reader.dia_frame_msms_windows()?)
            .map_err(|why| failed(format!("windows could not be read: {why}")))?;
        let (per_window, expansion) = ms2.splitting();
        // UniformMobility splits on timsrust's uncalibrated linear scale.
        let mobility = reader.linear_mobility_scale()?;
        let slices = tdf_spectra::dia_slices(
            &reader.dia_frame_msms_info()?,
            &groups,
            per_window,
            expansion,
            &mobility,
        )
        .map_err(|why| failed(format!("windows could not be read: {why}")))?;
        let processing = ms2.spectrum_processing_params;
        let dropped = Dropped::default();
        // Consecutive slices of one frame share a frame read.
        let mut by_frame: Vec<(usize, &[DiaSlice])> = Vec::new();
        let mut offset = 0;
        for chunk in slices.chunk_by(|a, b| a.frame == b.frame) {
            by_frame.push((offset, chunk));
            offset += chunk.len();
        }
        let spectra: Vec<RawSpectrum> = by_frame
            .par_iter()
            .flat_map_iter(|&(offset, chunk)| {
                let frame_id = chunk[0].frame;
                let frame = reader.read_frame(frame_id).map_err(|e| e.to_string());
                // timsrust took the retention time of the frame before (Id - 1).
                let rt = reader
                    .frame(frame_id.wrapping_sub(1))
                    .map(|info| info.retention_time_seconds)
                    .ok_or_else(|| format!("frame {frame_id} has no preceding frame"));
                let dropped = &dropped;
                chunk.iter().enumerate().filter_map(move |(k, slice)| {
                    let index = offset + k;
                    let spectrum = (|| {
                        let frame = frame.as_ref().map_err(Clone::clone)?;
                        let rt = rt.clone()?;
                        let (t, i) = scan_slice(frame, slice.scan_start, slice.scan_end)?;
                        let (tofs, intensities) = tdf_spectra::process(
                            t.to_vec(),
                            i.iter().map(|&x| u64::from(x)).collect(),
                            processing.smoothing_window,
                            processing.centroiding_window,
                        );
                        let scan = ((slice.scan_start + slice.scan_end) as f32 / 2.0) as u32;
                        let precursor = Precursor {
                            mz: slice.window.center() as f32,
                            spectrum_ref: Some(slice.frame.to_string()),
                            inverse_ion_mobility: Some(scale.precursor(
                                index,
                                slice.frame as usize,
                                scan,
                            )),
                            ..Precursor::default()
                        };
                        let window = IsolationWindow::from_center(
                            slice.window.center(),
                            slice.window.width(),
                            slice.window.collision_energy,
                        );
                        Ok::<_, String>(Self::ms2_spectrum(
                            file_id,
                            index.to_string(),
                            precursor,
                            rt,
                            window,
                            &tofs,
                            &intensities,
                            mz_scale,
                        ))
                    })();
                    spectrum.map_err(|why| dropped.drop_spectrum(why)).ok()
                })
            })
            .collect();
        dropped.log(path, "diaPASEF");
        log::info!(
            "read {} diaPASEF window spectra in {:#?}",
            spectra.len(),
            start.elapsed()
        );
        Ok(spectra)
    }

    #[allow(clippy::too_many_arguments)]
    fn ms2_spectrum(
        file_id: usize,
        id: String,
        mut precursor: Precursor,
        rt_seconds: f64,
        window: IsolationWindow,
        tof_indices: &[u32],
        intensities: &[u64],
        mz_scale: &LinearMzScale,
    ) -> RawSpectrum {
        let isolation_width = window.width();
        precursor.isolation_window = Some(Tolerance::Da(
            -isolation_width as f32 / 2.0,
            isolation_width as f32 / 2.0,
        ));
        RawSpectrum {
            file_id,
            precursors: vec![precursor],
            representation: Representation::Centroid,
            scan_start_time: rt_seconds as f32 / 60.0,
            ion_injection_time: rt_seconds as f32,
            total_ion_current: 0.0,
            mz: tof_indices
                .iter()
                .map(|&tof| mz_scale.mz(f64::from(tof)) as f32)
                .collect(),
            ms_level: 2,
            id,
            // timsrust reported f64 intensities; keep its double rounding.
            intensity: intensities
                .iter()
                .map(|&value| value as f64 as f32)
                .collect(),
            fragment_charges: None,
            mobility: None,
            acquisition: TIMS_TOF,
        }
    }
}

/// Denoise one MS1 frame, returning its surviving `(scan, tof, intensity)` points.
///
/// LZF frames can store more scans than `Frames.NumScans`; trailing empty
/// scans beyond it are dropped so the denoiser sees the frame's real extent.
fn denoise_frame(
    denoiser: &Ms1Denoiser<'_>,
    frame: &Frame,
    num_scans: usize,
) -> Result<Vec<(u32, u32, u32)>, DenoiseError> {
    if frame.intensities.is_empty() {
        return Ok(Vec::new());
    }
    let mut offsets = frame.scan_offsets.as_slice();
    while offsets.len() > num_scans + 1 && offsets[offsets.len() - 1] == offsets[offsets.len() - 2]
    {
        offsets = &offsets[..offsets.len() - 1];
    }
    denoiser.denoise(
        frame.id as usize,
        offsets,
        frame.tof_indices.clone(),
        frame.intensities.clone(),
    )
}

#[derive(Clone, Copy)]
struct ImsPeak {
    mz: f32,
    intensity: f32,
    im: f32,
}
const MAX_PEAKS: usize = 10_000;

/// Buffer that gets re-used on each thread to store the intermediates
/// of the centroiding for a single frame.
#[derive(Clone)]
struct PeakBuffer {
    peaks: Vec<ImsPeak>,
    order: Vec<usize>,
    agg_buff: Vec<ImsPeak>,
}

impl PeakBuffer {
    fn with_capacity(capacity: usize) -> Self {
        Self {
            peaks: Vec::with_capacity(capacity),
            order: Vec::with_capacity(capacity),
            agg_buff: Vec::with_capacity(MAX_PEAKS),
        }
    }

    fn with_frame(
        &mut self,
        frame: &Frame,
        scan_to_im: impl Fn(usize) -> f32,
        mz_scale: &LinearMzScale,
    ) {
        let expect_len = frame.tof_indices.len();
        self.expand_to_capacity(expect_len);

        let mz_iter = frame
            .tof_indices
            .iter()
            .map(|&index| mz_scale.mz(f64::from(index)) as f32);
        let intensities_iter = frame.intensities.iter().map(|&value| value as f32);
        let imss_iter = Self::expand_mobility_iter(&frame.scan_offsets, &scan_to_im);

        let peak_iter = mz_iter
            .zip(intensities_iter)
            .zip(imss_iter)
            .map(|((mz, intensity), im)| ImsPeak { mz, intensity, im });
        self.peaks.extend(peak_iter);
        assert_eq!(self.peaks.len(), expect_len);
        self.sort_and_order();
    }

    /// Load the frame's raw points, or only `survivors` of denoising when set.
    fn load(
        &mut self,
        frame: &Frame,
        survivors: Option<&[(u32, u32, u32)]>,
        scan_to_im: impl Fn(usize) -> f32,
        mz_scale: &LinearMzScale,
    ) {
        match survivors {
            Some(survivors) => self.with_survivors(survivors, scan_to_im, mz_scale),
            None => self.with_frame(frame, scan_to_im, mz_scale),
        }
    }

    /// Load denoised `(scan, tof, intensity)` points, which arrive in scan order.
    fn with_survivors(
        &mut self,
        survivors: &[(u32, u32, u32)],
        scan_to_im: impl Fn(usize) -> f32,
        mz_scale: &LinearMzScale,
    ) {
        self.expand_to_capacity(survivors.len());
        let mut current = None;
        let mut im = 0.0;
        for &(scan, tof, intensity) in survivors {
            if current != Some(scan) {
                current = Some(scan);
                im = scan_to_im(scan as usize);
            }
            self.peaks.push(ImsPeak {
                mz: mz_scale.mz(f64::from(tof)) as f32,
                intensity: intensity as f32,
                im,
            });
        }
        self.sort_and_order();
    }

    fn sort_and_order(&mut self) {
        // sort by mz ... bc binary searching on the mz space
        // for neighbors is the fastest way to find neighbors that I have tried.
        self.peaks.sort_by(|a, b| a.mz.partial_cmp(&b.mz).unwrap());

        // The "order" is sorted by intensity
        // This will be used later during the centroiding (for details check that implementation)
        self.order.extend(0..self.len());
        self.order.sort_unstable_by(|&a, &b| {
            self.peaks[b]
                .intensity
                .partial_cmp(&self.peaks[a].intensity)
                .unwrap_or(Ordering::Equal)
        });
    }

    fn clear(&mut self) {
        self.peaks.clear();
        self.order.clear();
        self.agg_buff.clear();
    }

    fn expand_to_capacity(&mut self, capacity: usize) {
        if capacity <= self.len() {
            return;
        }
        let diff = capacity - self.len();
        // Grow by whatever is the largest 20% of the current capacity
        // or the difference.
        let diff = diff.max(self.len() / 5);

        self.peaks.reserve(diff);
        self.order.reserve(diff);
        self.agg_buff.reserve(capacity);
    }

    fn len(&self) -> usize {
        self.peaks.len()
    }

    /// Expand the scan offset slice to mobilities.
    ///
    /// The scan offsets is in essence a run-length
    /// encoded vector of scan numbers that can be converter to the 1/k0
    /// values.
    ///
    /// Essentially ... the slice [0,4,5,5], would expand to
    /// [0,0,0,0,1]; 0 to 4 have index 0, 4 to 5 have index 1, 5 to 5 would
    /// have index 2 but its empty!
    ///
    /// Then this index can be converted to 1/K0.
    fn expand_mobility_iter<'a>(
        scan_offsets: &'a [usize],
        scan_to_im: &'a impl Fn(usize) -> f32,
    ) -> impl Iterator<Item = f32> + 'a {
        scan_offsets
            .windows(2)
            .enumerate()
            .filter_map(|(scan, w)| {
                let (lo, hi) = (w[0], w[1]);
                (hi > lo).then(|| (scan_to_im(scan), lo, hi))
            })
            .flat_map(|(im, lo, hi)| (lo..hi).map(move |_| im))
    }

    /// Centroiding of the IM-containing spectra
    ///
    /// This is a very rudimentary centroiding algorithm but... it seems to work well.
    /// It iterativelty goes over the peaks in decreasing intensity order and
    /// accumulates the intensity of the peaks surrounding the peak. (sort of
    /// like the first pass in dbscan).
    ///
    /// The preserved mobility and mz are the ones from the apex peak.
    /// A more complex version where the weighted mean is preserved is possible
    /// but I have seen only marginal gains and a lot more complexity + time.
    ///
    /// This dramatically reduces the number of peaks in the spectra
    /// which saves a ton of memory and time when doing LFQ, since we
    /// iterate over each peak.
    fn fastcentroid_frame(
        &mut self,
        mz_tol_ppm: f32,
        im_tol_pct: f32,
    ) -> (Vec<f32>, (Vec<f32>, Vec<f32>)) {
        // Make sure the array is mz sorted ... I should delete
        // this assertions once I am confident of the implementation.
        // but tbh, its not that slow and its simple.
        assert!(
            self.peaks.windows(2).all(|x| x[0].mz <= x[1].mz),
            "mz_array is not sorted"
        );
        assert!(self.agg_buff.is_empty(), "agg_buff is not empty");

        let mut global_num_included = 0;

        let utol = mz_tol_ppm / 1e6;
        let im_tol = im_tol_pct / 100.0;

        for &idx in &self.order {
            if self.peaks[idx].intensity <= 0.0 {
                continue;
            }
            if self.agg_buff.len() > MAX_PEAKS {
                let curr_loc_int = self.peaks[idx].intensity;
                if curr_loc_int > 200.0 {
                    log::debug!(
                        "Reached limit of the agg buffer at index {}/{} curr int={}",
                        idx,
                        self.len(),
                        curr_loc_int
                    );
                }
                break;
            }

            let mz = self.peaks[idx].mz;
            let im = self.peaks[idx].im;
            let da_tol = mz * utol;
            let left_e = mz - da_tol;
            let right_e = mz + da_tol;

            let ss_start = self.peaks.partition_point(|&x| x.mz < left_e);
            let ss_end = self.peaks.partition_point(|&x| x.mz <= right_e);

            let abs_im_tol = im * im_tol;
            let left_im = im - abs_im_tol;
            let right_im = im + abs_im_tol;

            let mut curr_intensity = 0.0;

            let mut num_includable = 0;
            for i in ss_start..ss_end {
                let im_i = self.peaks[i].im;
                if (self.peaks[i].intensity > 0.0) && im_i >= left_im && im_i <= right_im {
                    curr_intensity += self.peaks[i].intensity;
                    self.peaks[i].intensity = -1.0;
                    num_includable += 1;
                }
            }

            assert!(num_includable > 0, "At least 'itself' should be included");

            self.agg_buff.push(ImsPeak {
                mz,
                intensity: curr_intensity,
                im,
            });
            global_num_included += num_includable;

            if global_num_included == self.len() {
                log::debug!("All peaks were included in the centroiding");
                break;
            }
        }

        self.agg_buff
            .sort_unstable_by(|a, b| a.mz.partial_cmp(&b.mz).unwrap());
        // println!("Centroiding: Start len: {}; end len: {};", arr_len, result.len());
        // Ultra data is usually start: 40k end 10k,
        // HT2 data is usually start 400k end 40k, limiting to 10k
        // rarely leaves peaks with intensity > 200 ... ive never seen
        // it happen. -JSP 2025-Jan

        self.agg_buff
            .drain(..)
            .map(|x| (x.mz, (x.intensity, x.im)))
            .unzip()
    }
}

#[cfg(test)]
#[path = "../tests/unit/tdf.rs"]
mod tests;
