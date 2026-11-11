//! The retention-time by m/z grid the MS1 traces are built on, and its peaks.
//!
//! Moved out of `lfq.rs` unchanged; every item keeps its `lfq::` path.

use super::*;
use crate::ml::matrix::Matrix;

pub struct Grid {
    pub(super) rt_min: f32,
    pub(super) rt_step: f32,
    files: usize,
    /// File with the most confident PSM
    pub(super) reference_file_id: usize,
    /// Relative theoretical isotopic abundances
    pub distribution: [f32; N_ISOTOPES],
    /// Matrix of summed intensities for each isotopic trace in each file, divided
    /// among equally spaced retention time bins. This is a [N_FILES * N_ISOTOPES, GRID_SIZE]
    /// sized matrix.
    ///
    /// Isotopic summed intensities are arranged in consequtive rows, ordered by file. The first
    /// N_ISOTOPE rows correspond to file 0, then the following N_ISOTOPE rows correspond to file 1, etc.
    /// Indexing into the matrix is done by [(n_file * N_ISOTOPES + isotope, rt_window)]
    pub matrix: Matrix,
}

#[derive(Clone)]
pub struct Traces {
    /// Matrix of dot(MS1 ions, Grid.distribution). This collapses our N_FILES * N_ISOTOPES rows
    /// down to just N_FILES
    pub dot_product: Matrix,
    /// Matrix of spectral angles at each retention time for each file
    pub spectral_angle: Matrix,
    /// File with the most confident PSM
    pub(super) reference_file_id: usize,
    /// Consensus RT of bin 0 and the bin width.
    pub geometry: GridGeometry,
}

/// Position of a trace grid on the consensus RT axis.
#[derive(Copy, Clone, Debug, Default)]
pub struct GridGeometry {
    pub rt_min: f32,
    pub rt_step: f32,
}

impl GridGeometry {
    fn rt(&self, bin: f64) -> f32 {
        self.rt_min + (bin as f32) * self.rt_step
    }
}

#[derive(Clone, Debug, Default)]
pub struct Peak {
    /// Discretized retention time
    pub rt: usize,
    /// Intensity weighted normalized spectral angle
    pub spectral_angle: f64,
    /// Peak score
    pub score: f64,

    pub q_value: f32,
}

impl Traces {
    /// Calculate and apply time warping factors
    fn warp(&mut self) -> Vec<isize> {
        let time_warps = self.find_time_warps(&self.dot_product, WARP_SLACK);
        Self::apply_time_warps(&mut self.spectral_angle, &time_warps);
        Self::apply_time_warps(&mut self.dot_product, &time_warps);
        time_warps
    }

    /// Choose, align and integrate the peak: around the identification RTs'
    /// apexes when `seeds` are given (apex recentering), otherwise the
    /// best-scoring bin of the traced window.
    #[allow(clippy::type_complexity)]
    pub fn pick(
        &mut self,
        settings: &LfqSettings,
        seeds: Option<&[Option<usize>]>,
    ) -> Option<(
        Peak,
        Vec<Option<f64>>,
        Vec<Option<FileEvidence>>,
        PeakWindow,
    )> {
        // Precursors without an identification inside the grid (e.g. label
        // channels seeded by another channel) keep the window search.
        match seeds.filter(|seeds| seeds.iter().any(Option::is_some)) {
            Some(seeds) => self.integrate_apex(settings, seeds),
            None => self.integrate_window(settings),
        }
    }

    /// Isotope-consistent trace: signal where the isotope pattern matches.
    fn consistent(&self, file: usize, settings: &LfqSettings) -> Vec<f64> {
        self.dot_product
            .row_slice(file)
            .iter()
            .zip(self.spectral_angle.row_slice(file))
            .map(|(&dot, &angle)| {
                if angle >= settings.spectral_angle {
                    dot
                } else {
                    0.0
                }
            })
            .collect()
    }

    /// Integrate the peak whose apex is nearest the identification RTs.
    ///
    /// * Each seeded (identified) row climbs from its identification bin to
    ///   the apex of the isotope-consistent trace it sits on.
    /// * The consensus apex is the median of those apexes; each seeded row is
    ///   shifted so that its own apex lands on it.
    /// * Unseeded (transfer) rows get the usual warp search towards the sum of
    ///   the aligned seeded traces.
    /// * Integration bounds follow the summed seeded trace down to a valley
    ///   or [`APEX_BOUND_FRACTION`] of its apex, without fixed caps.
    #[allow(clippy::type_complexity)]
    pub fn integrate_apex(
        &mut self,
        settings: &LfqSettings,
        seeds: &[Option<usize>],
    ) -> Option<(
        Peak,
        Vec<Option<f64>>,
        Vec<Option<FileEvidence>>,
        PeakWindow,
    )> {
        let rows = self.dot_product.rows;
        let cols = self.dot_product.cols;
        let apexes = (0..rows)
            .map(|file| {
                let seed = (*seeds.get(file)?)?;
                climb_to_apex(&self.consistent(file, settings), seed)
            })
            .collect::<Vec<_>>();
        let mut found = apexes.iter().flatten().copied().collect::<Vec<_>>();
        if found.is_empty() {
            return None;
        }
        found.sort_unstable();
        let apex = found[(found.len() - 1) / 2];

        let mut shifts = vec![0isize; rows];
        let mut reference = vec![0.0; cols];
        let mut bounds_trace = vec![0.0; cols];
        for (file, file_apex) in apexes.iter().enumerate() {
            if let Some(file_apex) = file_apex {
                shifts[file] = *file_apex as isize - apex as isize;
                let consistent = self.consistent(file, settings);
                for (i, (r, b)) in reference.iter_mut().zip(&mut bounds_trace).enumerate() {
                    let j = i as isize + shifts[file];
                    if j >= 0 && (j as usize) < cols {
                        *r += self.dot_product[(file, j as usize)];
                        *b += consistent[j as usize];
                    }
                }
            }
        }
        for (file, file_apex) in apexes.iter().enumerate() {
            if file_apex.is_none() {
                shifts[file] =
                    Self::best_shift(&reference, self.dot_product.row_slice(file), WARP_SLACK);
            }
        }
        Self::apply_time_warps(&mut self.spectral_angle, &shifts);
        Self::apply_time_warps(&mut self.dot_product, &shifts);

        let (left, right) = peak_bounds(&bounds_trace, apex, APEX_BOUND_FRACTION);

        // Same form as the hybrid score, with retention-time proximity
        // measured from the identification RTs rather than the window center.
        let mut max = 0.0f64;
        let mut at_apex = (0.0, 1.0);
        for col in 0..cols {
            let mut summed_int = 1.0;
            let mut weighted = 0.0;
            for (sa, dotp) in self.spectral_angle.col(col).zip(self.dot_product.col(col)) {
                weighted += sa * dotp;
                summed_int += dotp;
            }
            max = max.max(summed_int);
            if col == apex {
                at_apex = (weighted / summed_int, summed_int);
            }
        }
        let mut seeded = seeds.iter().flatten().copied().collect::<Vec<_>>();
        seeded.sort_unstable();
        let center = seeded[(seeded.len() - 1) / 2];
        let half = (cols / 2).max(1) as f64;
        let proximity = (1.0 - (apex as f64 - center as f64).abs() / half).max(0.0);
        let (spectral_angle, intensity) = at_apex;
        let score =
            spectral_angle.max(0.0).powi(3) * proximity.powf(0.33) * (intensity / max).sqrt();
        if score <= 0.0 || !score.is_finite() {
            return None;
        }
        let peak = Peak {
            rt: apex,
            spectral_angle,
            score,
            q_value: 0.0,
        };
        let window = PeakWindow {
            rt: apex,
            left,
            right,
            shifts,
            reference,
        };
        let (areas, evidence) = self.window_evidence(&window, settings);
        Some((peak, areas, evidence, window))
    }

    /// Find time warping offsets for each file that maximize the dot product
    /// with the most intense run
    ///
    /// * Use the LC-MS run with the most confident PSM for a peptide as the reference run
    /// * For each LC-MS run, find the time warping shif that maximizes dot product
    ///   with the `reference` run
    pub fn find_time_warps(&self, matrix: &Matrix, slack: isize) -> Vec<isize> {
        Self::time_warps_to(matrix.row_slice(self.reference_file_id), matrix, slack)
    }

    /// Time warping offsets for each row of `matrix` that maximize the dot
    /// product with `reference`
    fn time_warps_to(reference: &[f64], matrix: &Matrix, slack: isize) -> Vec<isize> {
        (0..matrix.rows)
            .map(|row| Self::best_shift(reference, matrix.row_slice(row), slack))
            .collect()
    }

    /// Offset of `run` that maximizes its dot product with `reference`.
    fn best_shift(reference: &[f64], run: &[f64], slack: isize) -> isize {
        let mut best_offset = (0isize, 0.0);
        for offset in -slack..=slack {
            let mut dot = 0.0;
            for (i, ref_int) in reference.iter().enumerate() {
                let j = i as isize + offset;
                if j >= 0 && j < run.len() as isize {
                    dot += ref_int * run[j as usize];
                }
            }

            // Without evidence for a shift, preserve the trace. For tied
            // correlations, prefer the smallest displacement.
            if dot > best_offset.1 || (dot == best_offset.1 && offset.abs() < best_offset.0.abs()) {
                best_offset = (offset, dot);
            }
        }
        best_offset.0
    }

    /// Perform local Correlation Optimization Warping
    pub(super) fn apply_time_warps(matrix: &mut Matrix, time_warps: &[isize]) {
        for (row, warp) in time_warps.iter().enumerate() {
            let run = matrix.row_slice_mut(row);
            let mut shifted = vec![0.0; run.len()];
            for (i, val) in shifted.iter_mut().enumerate() {
                let j = i as isize + warp;
                if j >= 0 && j < run.len() as isize {
                    *val = run[j as usize];
                }
            }
            run.copy_from_slice(&shifted);
        }
    }

    pub fn scores(&self, strategy: PeakScoringStrategy) -> (Vec<f64>, Vec<f64>) {
        let mut spectral = Vec::with_capacity(self.spectral_angle.cols);
        let mut intensity = Vec::with_capacity(self.spectral_angle.cols);
        let mut max = 0.0f64;
        for col in 0..self.spectral_angle.cols {
            let mut summed_int = 1.0;
            let mut weighted = 0.0;
            for (sa, dotp) in self.spectral_angle.col(col).zip(self.dot_product.col(col)) {
                weighted += sa * dotp;
                summed_int += dotp;
            }
            spectral.push(weighted / summed_int);
            intensity.push(summed_int);
            max = max.max(summed_int);
        }

        let center = self.spectral_angle.cols as isize / 2;
        let scores = spectral
            .iter()
            .zip(intensity.iter())
            .enumerate()
            .map(|(rt, (s, i))| match strategy {
                PeakScoringStrategy::RetentionTime => {
                    (1.0 - ((rt as isize - center).abs() as f64 / center as f64)).powf(0.33)
                }
                PeakScoringStrategy::SpectralAngle => *s,
                PeakScoringStrategy::Intensity => (*i / max).sqrt(),
                PeakScoringStrategy::Hybrid => {
                    let rt = 1.0 - ((rt as isize - center).abs() as f64 / center as f64);
                    s.powi(3) * rt.powf(0.33) * (*i / max).sqrt()
                    // s.powi(3) * (*i / max).sqrt()
                }
            })
            .collect();
        (scores, spectral)
    }

    /// Align and integrate MS1 traces across files
    ///
    /// * Calculate time warping factors for each file, aligning them so that
    ///   the correlation between them is maximized (we just do a dot product)
    /// * Locate the retention time window corresponding to the maximum average
    ///   angle observed across all of the files
    /// * Integrate all of the MS1 traces within said window, returning a vector
    ///   of length `n_files` containing the summed MS1 intensities
    pub fn integrate(&mut self, settings: &LfqSettings) -> Option<(Peak, Vec<Option<f64>>)> {
        self.integrate_with_evidence(settings)
            .map(|(peak, areas, _)| (peak, areas))
    }

    #[allow(clippy::type_complexity)]
    pub fn integrate_with_evidence(
        &mut self,
        settings: &LfqSettings,
    ) -> Option<(Peak, Vec<Option<f64>>, Vec<Option<FileEvidence>>)> {
        self.integrate_window(settings)
            .map(|(peak, areas, evidence, _)| (peak, areas, evidence))
    }

    /// Integrate the best peak and also return its window, so that a paired
    /// decoy can be evaluated at exactly the same place (see [`Traces::paired_evidence`]).
    #[allow(clippy::type_complexity)]
    pub fn integrate_window(
        &mut self,
        settings: &LfqSettings,
    ) -> Option<(
        Peak,
        Vec<Option<f64>>,
        Vec<Option<FileEvidence>>,
        PeakWindow,
    )> {
        let shifts = self.warp();

        let (scores, spectral) = self.scores(settings.peak_scoring);
        let mut best = Peak::default();
        for (rt, s) in scores.iter().enumerate() {
            if *s > best.score && spectral[rt] >= settings.spectral_angle {
                best.score = *s;
                best.rt = rt;
            }
        }

        if best.score == 0.0 {
            return None;
        }

        // Inclusive boundaries contain only bins that pass both thresholds.
        let mut left = best.rt;
        let mut right = best.rt;

        let threshold = best.score * 0.50;

        // Don't let peaks extend more than one fifth of the grid to either side.
        let max_width = scores.len() / 5;
        while left > best.rt.saturating_sub(max_width)
            && scores[left - 1] >= threshold
            && spectral[left - 1] >= settings.spectral_angle
        {
            left -= 1;
        }

        while right < (scores.len() - 1).min(best.rt + max_width)
            && scores[right + 1] >= threshold
            && spectral[right + 1] >= settings.spectral_angle
        {
            right += 1;
        }

        let mut summed_int = 1.0;
        let mut weighted = 0.0;
        for (sa, dotp) in self
            .spectral_angle
            .col(best.rt)
            .zip(self.dot_product.col(best.rt))
        {
            weighted += sa * dotp;
            summed_int += dotp;
        }
        best.spectral_angle = weighted / summed_int;
        let window = PeakWindow {
            rt: best.rt,
            left,
            right,
            shifts,
            reference: self.dot_product.row_slice(self.reference_file_id).to_vec(),
        };
        let (areas, evidence) = self.window_evidence(&window, settings);
        Some((best, areas, evidence, window))
    }

    /// Evaluate these (decoy) traces at a peak chosen from other traces: the
    /// same apex, integration window and reference trace. Each file gets the
    /// same local warp search the target files got (towards that reference),
    /// so a file row answers the same question as the target row it is paired
    /// with: can an isotope envelope be aligned to the expected elution?
    ///
    /// With apex recentering (`seeds`), a seeded row instead climbs from the
    /// same identification bin as the target row and is shifted so that its
    /// own apex lands on the target's apex, exactly as the target row was. A
    /// row whose climb finds no signal keeps the warp search, as a target row
    /// does in [`Traces::integrate_apex`].
    pub fn paired_evidence(
        &mut self,
        window: &PeakWindow,
        settings: &LfqSettings,
        seeds: Option<&[Option<usize>]>,
    ) -> Vec<Option<FileEvidence>> {
        let mut shifts = Self::time_warps_to(&window.reference, &self.dot_product, WARP_SLACK);
        for (file, seed) in seeds.unwrap_or_default().iter().enumerate() {
            if let Some(apex) =
                seed.and_then(|seed| climb_to_apex(&self.consistent(file, settings), seed))
            {
                shifts[file] = apex as isize - window.rt as isize;
            }
        }
        Self::apply_time_warps(&mut self.spectral_angle, &shifts);
        Self::apply_time_warps(&mut self.dot_product, &shifts);
        let window = PeakWindow {
            shifts,
            reference: window.reference.clone(),
            ..*window
        };
        self.window_evidence(&window, settings).1
    }

    fn window_evidence(
        &self,
        window: &PeakWindow,
        settings: &LfqSettings,
    ) -> (Vec<Option<f64>>, Vec<Option<FileEvidence>>) {
        let mut areas = Vec::with_capacity(self.dot_product.rows);
        for file in 0..self.dot_product.rows {
            let area = match settings.integration {
                IntegrationStrategy::Sum => self.dot_product.row_slice(file)
                    [window.left..=window.right]
                    .iter()
                    .sum::<f64>(),
                IntegrationStrategy::Apex => self.dot_product.row_slice(file)[window.rt],
            };

            areas.push((area.is_finite() && area > 0.0).then_some(area));
        }

        let reference = &window.reference;
        let reference_norm = reference.iter().map(|x| x * x).sum::<f64>().sqrt();
        let evidence = areas
            .iter()
            .enumerate()
            .map(|(file, area)| {
                area.map(|_| {
                    let trace = self.dot_product.row_slice(file);
                    let norm = trace.iter().map(|x| x * x).sum::<f64>().sqrt();
                    let dot = trace.iter().zip(reference).map(|(x, y)| x * y).sum::<f64>();
                    let trace_cosine = if norm > 0.0 && reference_norm > 0.0 {
                        (dot / (norm * reference_norm)).clamp(0.0, 1.0)
                    } else {
                        0.0
                    };
                    let spectral_angle = self.spectral_angle[(file, window.rt)].clamp(0.0, 1.0);
                    let proximity = (1.0
                        - window.shifts[file].unsigned_abs() as f64 / self.dot_product.cols as f64)
                        .max(0.0);
                    // Positions in the unwarped grid of this file, as consensus RT.
                    let shift = window.shifts[file] as f64;
                    let fwhm = half_maximum_width(trace, window.rt)
                        .map(|bins| bins as f32 * self.geometry.rt_step);
                    FileEvidence {
                        transfer_candidate: false,
                        spectral_angle,
                        trace_cosine,
                        rt_shift_bins: window.shifts[file] as i32,
                        score: spectral_angle.powi(3) * trace_cosine * proximity,
                        extraction_q_value: None,
                        apex_rt: self.geometry.rt(window.rt as f64 + shift),
                        peak_start_rt: self.geometry.rt(window.left as f64 + shift),
                        peak_end_rt: self.geometry.rt(window.right as f64 + shift),
                        fwhm,
                        id_apex_offset: None,
                    }
                })
            })
            .collect();
        (areas, evidence)
    }
}

/// Where a precursor's cross-run peak was found: local warps, apex and
/// integration bounds, and the warped reference trace.
#[derive(Clone, Debug)]
pub struct PeakWindow {
    pub rt: usize,
    pub left: usize,
    pub right: usize,
    pub shifts: Vec<isize>,
    pub reference: Vec<f64>,
}

impl Grid {
    pub fn new(
        entry: &PrecursorRange,
        rt_tol: f32,
        distribution: [f32; N_ISOTOPES],
        files: usize,
        grid_size: usize,
    ) -> Grid {
        let matrix = Matrix::new(
            vec![0.0; grid_size * files * N_ISOTOPES],
            files * N_ISOTOPES,
            grid_size,
        );
        let rt_step = (rt_tol * 2.0) / (grid_size) as f32;

        Grid {
            rt_min: entry.rt - rt_tol,
            rt_step,
            distribution,
            matrix,
            files,
            reference_file_id: entry.file_id,
        }
    }

    /// Add a data point to the integration grid
    pub fn add_entry(&mut self, spectrum_rt: f32, isotope: usize, file_id: usize, intensity: f32) {
        let bin_lo = ((spectrum_rt - self.rt_min) / self.rt_step).floor() as usize;
        let bin_lo = bin_lo.min(self.matrix.cols - 1);
        let bin_hi = (bin_lo + 1).min(self.matrix.cols - 1);

        let bin_lo_rt = bin_lo as f32 * self.rt_step + self.rt_min;
        // what fraction [0.0, 1.0] of the way are we to the higher bin?
        let interp = (spectrum_rt - bin_lo_rt) / self.rt_step;

        self.matrix[(file_id * N_ISOTOPES + isotope, bin_lo)] +=
            ((1.0 - interp) * intensity) as f64;
        self.matrix[(file_id * N_ISOTOPES + isotope, bin_hi)] += (interp * intensity) as f64;
    }

    /// Combine individual isotopic traces across files into aligned, summed
    /// MS1 traces for each file
    ///
    /// * Perform gaussian smoothing on summed intensities
    /// * Calculate normalized spectral angle for observed isotopic distribution
    ///   relative to theoretical distribution
    pub fn summarize_traces(&mut self) -> Traces {
        let k = gaussian_kernel(0.5, K_WIDTH);

        let mut spectral_angle = Matrix::new(
            vec![0.0; self.files * self.matrix.cols],
            self.files,
            self.matrix.cols,
        );

        let mut dot_product = spectral_angle.clone();

        // square root of the summed squared relative abundances of theoretical
        // isotopic distribution
        let ss_dist = self
            .distribution
            .iter()
            .map(|&x| (x as f64).powi(2))
            .sum::<f64>()
            .sqrt();

        for file in 0..self.files {
            let mut summed_squared_intensities = vec![0.0; self.matrix.cols];
            for isotope in 0..N_ISOTOPES {
                let convolved = convolve(self.matrix.row_slice(file * N_ISOTOPES + isotope), &k);
                for (col, intensity) in convolved.iter().enumerate() {
                    spectral_angle[(file, col)] += intensity * self.distribution[isotope] as f64;
                    summed_squared_intensities[col] += intensity.powi(2);
                }
                self.matrix
                    .row_slice_mut(file * N_ISOTOPES + isotope)
                    .copy_from_slice(&convolved);
            }

            for (col, ss) in summed_squared_intensities.iter().enumerate() {
                let dot = spectral_angle[(file, col)];
                let similarity = if *ss > 0.0 {
                    dot / (ss.sqrt() * ss_dist)
                } else {
                    0.0
                };

                // Calculate the normalized spectral angle
                spectral_angle[(file, col)] =
                    1.0 - 2.0 * similarity.clamp(-1.0, 1.0).acos() / std::f64::consts::PI;
                dot_product[(file, col)] = dot;
            }
        }

        Traces {
            dot_product,
            spectral_angle,
            reference_file_id: self.reference_file_id,
            geometry: GridGeometry {
                rt_min: self.rt_min,
                rt_step: self.rt_step,
            },
        }
    }
}
