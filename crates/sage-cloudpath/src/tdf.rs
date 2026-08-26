use rayon::prelude::*;
use sage_core::{
    mass::Tolerance,
    spectrum::{Precursor, RawSpectrum, Representation},
};
use serde::{Deserialize, Serialize};
use std::{cmp::Ordering, path::Path};
use timsrust::{
    core::{Converter, Frame, Im, MSLevel, Precursor as TimsrustPrecursor, ScanIndex},
    tdf::{
        FrameWindowSplittingConfiguration, QuadWindowExpansionStrategy, SpectrumProcessingParams,
        SpectrumReaderConfig,
    },
    ImConverter, MzConverter, SpectrumReader, TimsTofPath,
};
pub struct TdfReader;

#[derive(Deserialize, Serialize, Debug, Clone, Copy)]
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

#[derive(Deserialize, Serialize, Debug, Clone, Copy)]
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

#[derive(Deserialize, Serialize, Debug, Clone, Copy)]
pub enum BrukerFrameWindowSplittingConfig {
    Quadrupole(BrukerQuadWindowExpansionStrategy),
    Window(BrukerQuadWindowExpansionStrategy),
}

impl Default for BrukerFrameWindowSplittingConfig {
    fn default() -> Self {
        Self::Quadrupole(BrukerQuadWindowExpansionStrategy::default())
    }
}

#[derive(Default, Deserialize, Serialize, Debug, Clone, Copy)]
pub struct BrukerSpectrumConfig {
    pub spectrum_processing_params: BrukerSpectrumProcessingConfig,
    pub frame_splitting_params: BrukerFrameWindowSplittingConfig,
}

impl BrukerSpectrumConfig {
    fn into_timsrust(self) -> SpectrumReaderConfig<ImConverter> {
        let processing = self.spectrum_processing_params;
        let expansion = |strategy| match strategy {
            BrukerQuadWindowExpansionStrategy::None => QuadWindowExpansionStrategy::None,
            BrukerQuadWindowExpansionStrategy::Even(count) => {
                QuadWindowExpansionStrategy::Even(count)
            }
            BrukerQuadWindowExpansionStrategy::UniformMobility(span_step, _) => {
                QuadWindowExpansionStrategy::UniformMobility(span_step, None)
            }
            BrukerQuadWindowExpansionStrategy::UniformScan(span_step) => {
                QuadWindowExpansionStrategy::UniformScan(span_step)
            }
        };
        let frame_splitting_params = match self.frame_splitting_params {
            BrukerFrameWindowSplittingConfig::Quadrupole(strategy) => {
                FrameWindowSplittingConfiguration::Quadrupole(expansion(strategy))
            }
            BrukerFrameWindowSplittingConfig::Window(strategy) => {
                FrameWindowSplittingConfiguration::Window(expansion(strategy))
            }
        };
        SpectrumReaderConfig {
            spectrum_processing_params: SpectrumProcessingParams {
                smoothing_window: processing.smoothing_window,
                centroiding_window: processing.centroiding_window,
                calibration_tolerance: processing.calibration_tolerance,
                calibrate: processing.calibrate,
            },
            frame_splitting_params,
        }
    }
}

#[derive(Deserialize, Serialize, Debug, Clone, Copy)]
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

#[derive(Default, Deserialize, Serialize, Debug, Clone, Copy)]
pub struct BrukerProcessingConfig {
    pub ms2: BrukerSpectrumConfig,
    pub ms1: BrukerMS1CentoidingConfig,
}

impl TdfReader {
    pub fn parse(
        &self,
        path_name: impl AsRef<Path>,
        file_id: usize,
        config: BrukerProcessingConfig,
        requires_ms1: bool,
    ) -> Result<Vec<RawSpectrum>, timsrust::TimsRustError> {
        let path = TimsTofPath::new(path_name.as_ref().to_string_lossy())?;
        let spectrum_reader = SpectrumReader::build()
            .with_path(&path)
            .with_config(config.ms2.into_timsrust())
            .finalize()?;
        let mut spectra = self.read_msn_spectra(file_id, &spectrum_reader)?;
        if requires_ms1 {
            let ms1s = self.read_ms1_spectra(&path_name, file_id, config.ms1)?;
            spectra.extend(ms1s);
        }

        Ok(spectra)
    }

    fn read_ms1_spectra(
        &self,
        path_name: impl AsRef<Path>,
        file_id: usize,
        config: BrukerMS1CentoidingConfig,
    ) -> Result<Vec<RawSpectrum>, timsrust::TimsRustError> {
        let start = std::time::Instant::now();
        let path = TimsTofPath::new(path_name.as_ref().to_string_lossy())?;
        let frame_reader = path.frame_reader().map_err(|error| match error {
            timsrust::TimsTofFrameReaderError::FrameReaderError(error) => error,
            timsrust::TimsTofFrameReaderError::NotSupported => {
                unreachable!("TimsTofPath recognized a non-TDF path as a Bruker directory")
            }
        })?;
        let mz_converter = path
            .mz_converter()
            .expect("TDF paths always provide an m/z converter");
        let ims_converter = path
            .im_converter()
            .expect("TDF paths always provide an ion mobility converter");
        let tol_ppm = config.mz_ppm;
        let im_tol_pct = config.ims_pct;

        let ms1_spectra: Vec<RawSpectrum> = frame_reader
            .parallel_filter(|frame| frame.info().ms_level() == MSLevel::MS1)
            .map_init(
                || PeakBuffer::with_capacity(2 * MAX_PEAKS),
                |buffer, frame| match frame {
                    Ok(frame) => {
                        buffer.clear();
                        buffer.with_frame(&frame, &ims_converter, &mz_converter);

                        // Squash the mobility dimension
                        let (mz, (intensity, mobility)): (Vec<f32>, (Vec<f32>, Vec<f32>)) =
                            buffer.fastcentroid_frame(tol_ppm, im_tol_pct);

                        let scan_start_time = frame.info().rt_in_seconds() as f32 / 60.0;
                        let ion_injection_time = 100.0; // This is made up, in theory we can read
                                                        // if from the tdf file
                        let total_ion_current = intensity.iter().sum::<f32>();
                        let id = frame.index().to_string();

                        let spec = RawSpectrum {
                            file_id,
                            precursors: vec![],
                            representation: Representation::Centroid,
                            scan_start_time,
                            ion_injection_time,
                            mz,
                            ms_level: 1,
                            id,
                            intensity,
                            total_ion_current,
                            fragment_charges: None,
                            mobility: Some(mobility),
                        };
                        Some(spec)
                    }
                    Err(x) => {
                        log::error!("error parsing spectrum: {:?}", x);
                        None
                    }
                },
            )
            .flatten()
            .collect();
        log::info!(
            "read {} ms1 spectra in {:#?}",
            ms1_spectra.len(),
            start.elapsed()
        );
        Ok(ms1_spectra)
    }

    fn read_msn_spectra(
        &self,
        file_id: usize,
        spectrum_reader: &SpectrumReader,
    ) -> Result<Vec<RawSpectrum>, timsrust::TimsRustError> {
        let spectra: Vec<RawSpectrum> = spectrum_reader
            .par_iter()
            .filter_map(|result| match result {
                Ok(dda_spectrum) => match dda_spectrum.precursor() {
                    Some(dda_precursor) => {
                        let mut precursor = Self::parse_precursor(dda_precursor);
                        let isolation_width = f64::from(dda_spectrum.isolation_window().width());
                        precursor.isolation_window = Option::from(Tolerance::Da(
                            -isolation_width as f32 / 2.0,
                            isolation_width as f32 / 2.0,
                        ));
                        let spectrum: RawSpectrum = RawSpectrum {
                            file_id,
                            precursors: vec![precursor],
                            representation: Representation::Centroid,
                            scan_start_time: f64::from(dda_precursor.rt()) as f32 / 60.0,
                            ion_injection_time: f64::from(dda_precursor.rt()) as f32,
                            total_ion_current: 0.0,
                            mz: dda_spectrum
                                .mz_values()
                                .iter()
                                .map(|&value| f64::from(value) as f32)
                                .collect(),
                            ms_level: 2,
                            id: dda_spectrum.index().to_string(),
                            intensity: dda_spectrum
                                .intensities()
                                .iter()
                                .map(|&value| value as f32)
                                .collect(),
                            fragment_charges: None,
                            mobility: None,
                        };
                        Some(spectrum)
                    }
                    None => None,
                },
                Err(error) => {
                    log::warn!("error parsing Bruker MS2 spectrum: {error}");
                    None
                }
            })
            .collect();
        Ok(spectra)
    }

    fn parse_precursor(dda_precursor: &TimsrustPrecursor) -> Precursor {
        Precursor {
            mz: f64::from(dda_precursor.mz()) as f32,
            charge: dda_precursor.charge().map(|charge| i8::from(charge) as u8),
            intensity: dda_precursor.intensity().map(|value| value as f32),
            spectrum_ref: Option::from(dda_precursor.frame_index().to_string()),
            inverse_ion_mobility: Option::from(f64::from(dda_precursor.im()) as f32),
            ..Precursor::default()
        }
    }
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
        ims_converter: &ImConverter,
        mz_converter: &MzConverter,
    ) {
        let expect_len = frame.ions().tof_indices().len();
        self.expand_to_capacity(expect_len);

        let mz_iter = frame
            .ions()
            .tof_indices()
            .iter()
            .map(|&index| f64::from(mz_converter.convert(index)) as f32);
        let intensities_iter = frame
            .ions()
            .intensities()
            .iter()
            .map(|&value| u32::from(value) as f32);
        let imss_iter = Self::expand_mobility_iter(frame.ions().scan_offsets(), ims_converter);

        let peak_iter = mz_iter
            .zip(intensities_iter)
            .zip(imss_iter)
            .map(|((mz, intensity), im)| ImsPeak { mz, intensity, im });
        self.peaks.extend(peak_iter);
        assert_eq!(self.peaks.len(), expect_len);

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
    /// Then this index can be converted using the Scan2ImConverter.convert
    ///
    /// ... This should problably be implemented and exposed in timsrust.
    fn expand_mobility_iter<'a, C>(
        scan_offsets: &'a [usize],
        ims_converter: &'a C,
    ) -> impl Iterator<Item = f32> + 'a
    where
        C: Converter<ScanIndex, Im> + 'a,
    {
        let ims_iter = scan_offsets
            .windows(2)
            .enumerate()
            .filter_map(|(i, w)| {
                let num = w[1] - w[0];
                if num == 0 {
                    return None;
                }
                let lo = w[0];
                let hi = w[1];

                let scan_index = ScanIndex::try_from(i).expect("scan index exceeds u32 range");
                let im = f64::from(ims_converter.convert(scan_index)) as f32;
                Some((im, lo, hi))
            })
            .flat_map(|(im, lo, hi)| (lo..hi).map(move |_| im));
        ims_iter
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
