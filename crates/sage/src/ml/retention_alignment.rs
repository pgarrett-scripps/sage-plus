//! Perform robust global retention time alignment using shared peptide landmarks.
//!
//! The alignment is deliberately split into two stages:
//! 1. A RANSAC-style least-median affine fit rejects grossly inconsistent landmarks.
//! 2. Median landmarks are fit with isotonic regression and linearly interpolated.
//!
//! The resulting warp can represent local gradient distortion while remaining monotone.
//! Runs with too few landmarks fall back to the affine model. If LFQ is enabled, MS1
//! apex times are transformed through the same warp as PSM retention times.

use std::collections::HashMap;
use std::hash::BuildHasherDefault;
use std::sync::atomic::AtomicU32;

use super::matrix::Matrix;
use crate::database::PeptideIx;
use crate::scoring::Feature;
use dashmap::DashMap;
use fnv::FnvHasher;
use rayon::prelude::*;
use serde::{Deserialize, Serialize};

type FnvDashMap<K, V> = DashMap<K, V, BuildHasherDefault<FnvHasher>>;

const MIN_NONLINEAR_LANDMARKS: usize = 16;
const MIN_NONLINEAR_SPAN: f64 = 0.25;
const MAX_RANSAC_POINTS: usize = 2_048;
const RANSAC_TRIALS: usize = 256;
const MIN_X_SEPARATION: f64 = 0.05;

#[derive(Copy, Clone, Debug, Default, Deserialize, Eq, PartialEq, Serialize)]
#[serde(rename_all = "snake_case")]
pub enum AlignmentMethod {
    /// Preserve Sage's existing ordinary least-squares alignment.
    #[default]
    Linear,
    /// Robust outlier filtering followed by a monotone piecewise-linear warp.
    Nonlinear,
}

fn median(values: &mut [f64]) -> f64 {
    let mid = values.len() / 2;
    values.select_nth_unstable_by(mid, f64::total_cmp);
    if values.len() % 2 == 1 {
        values[mid]
    } else {
        let lower = values[..mid]
            .iter()
            .copied()
            .max_by(f64::total_cmp)
            .unwrap_or(values[mid]);
        (lower + values[mid]) / 2.0
    }
}

fn max_rt_by_file(features: &[Feature], n_files: usize) -> Vec<f64> {
    let max_rt = (0..n_files).map(|_| AtomicU32::new(0)).collect::<Vec<_>>();

    features.par_iter().for_each(|feat| {
        max_rt[feat.file_id].fetch_max(feat.rt.ceil() as u32, std::sync::atomic::Ordering::SeqCst);
    });

    max_rt
        .into_iter()
        .map(|v| v.load(std::sync::atomic::Ordering::Acquire).max(1) as f64)
        .collect()
}

/// Return a map from peptide to the earliest high-confidence RT in each file.
fn mean_rt_by_file(features: &[Feature]) -> FnvDashMap<PeptideIx, HashMap<usize, f64>> {
    let rts: FnvDashMap<PeptideIx, HashMap<usize, f64>> = DashMap::default();
    features
        .par_iter()
        .filter(|feat| feat.label == 1 && feat.spectrum_q <= 0.01)
        .for_each(|feat| {
            rts.entry(feat.peptide_idx)
                .or_default()
                .entry(feat.file_id)
                .and_modify(|f| *f = f.min(feat.rt as f64))
                .or_insert(feat.rt as f64);
        });
    rts
}

fn rt_matrix(features: &[Feature], max_rt: &[f64], shared_only: bool) -> Matrix {
    let mean_rt = mean_rt_by_file(features);

    let mat: Vec<Vec<f64>> = mean_rt
        .par_iter()
        .filter(|entry| !shared_only || entry.value().len() >= 2)
        .map(|entry| {
            let mut row = vec![f64::NAN; max_rt.len()];
            for (&file_id, &rt) in entry.value() {
                row[file_id] = rt / max_rt[file_id];
            }
            row
        })
        .collect();
    let n = mat.len();
    let mat: Vec<f64> = mat.into_par_iter().flatten().collect();

    Matrix::new(mat, n, max_rt.len())
}

#[derive(Clone, Debug)]
pub struct Alignment {
    pub file_id: usize,
    pub max_rt: f32,
    pub slope: f32,
    pub intercept: f32,
    /// Monotone `(observed_rt, consensus_rt)` knots in normalized RT units.
    /// Empty for an affine-only alignment.
    pub knots: Vec<(f32, f32)>,
}

impl Alignment {
    /// Transform a retention time in the original run's time units into consensus RT.
    pub fn transform(&self, rt: f32) -> f32 {
        let x = rt / self.max_rt;
        if self.knots.len() < 2 {
            return x * self.slope + self.intercept;
        }

        let upper = self.knots.partition_point(|&(knot_x, _)| knot_x <= x);
        let (left, right) = match upper {
            0 => (self.knots[0], self.knots[1]),
            n if n == self.knots.len() => (self.knots[n - 2], self.knots[n - 1]),
            n => (self.knots[n - 1], self.knots[n]),
        };

        let width = right.0 - left.0;
        if width <= f32::EPSILON {
            left.1
        } else {
            left.1 + (x - left.0) * (right.1 - left.1) / width
        }
    }
}

fn ordinary_least_squares(points: &[(f64, f64)]) -> (f64, f64) {
    if points.len() < 2 {
        return (1.0, 0.0);
    }

    let n = points.len() as f64;
    let (sum_x, sum_y) = points
        .iter()
        .fold((0.0, 0.0), |(sx, sy), &(x, y)| (sx + x, sy + y));
    let mean_x = sum_x / n;
    let mean_y = sum_y / n;
    let (ss_xy, ss_xx) = points.iter().fold((0.0, 1e-8), |(xy, xx), &(x, y)| {
        (xy + (x - mean_x) * (y - mean_y), xx + (x - mean_x).powi(2))
    });
    let slope = ss_xy / ss_xx;
    let intercept = mean_y - slope * mean_x;

    if slope.is_finite() && slope > 0.0 && intercept.is_finite() {
        (slope, intercept)
    } else {
        (1.0, 0.0)
    }
}

fn legacy_linear_alignment(file_id: usize, max_rt: f64, points: &[(f64, f64)]) -> Alignment {
    if points.len() < 2 {
        return Alignment {
            file_id,
            max_rt: max_rt as f32,
            slope: 1.0,
            intercept: 0.0,
            knots: Vec::new(),
        };
    }

    let n = points.len() as f64;
    let (sum_x, sum_y, dot) = points.iter().fold((0.0, 0.0, 0.0), |acc, &(x, y)| {
        (acc.0 + x, acc.1 + y, acc.2 + x * y)
    });
    let mean_x = sum_x / n;
    let mean_y = sum_y / n;
    let ss_xy = dot - n * mean_x * mean_y;
    let ss_xx = points
        .iter()
        .fold(1e-8, |sum, &(x, _)| sum + (x - mean_x).powi(2));
    let mut slope = ss_xy / ss_xx;
    let mut intercept = mean_y - slope * mean_x;
    if !slope.is_finite() {
        slope = 1.0;
    }
    if !intercept.is_finite() {
        intercept = 0.0;
    }

    Alignment {
        file_id,
        max_rt: max_rt as f32,
        slope: slope as f32,
        intercept: intercept as f32,
        knots: Vec::new(),
    }
}

/// Find an affine trend using a deterministic least-median-of-squares variant of RANSAC,
/// then return that trend and the indices of landmarks consistent with it.
fn robust_affine_inliers(points: &[(f64, f64)]) -> ((f64, f64), Vec<usize>) {
    if points.len() < 3 {
        return (ordinary_least_squares(points), (0..points.len()).collect());
    }

    let sample: Vec<_> = if points.len() <= MAX_RANSAC_POINTS {
        points.to_vec()
    } else {
        (0..MAX_RANSAC_POINTS)
            .map(|i| points[i * (points.len() - 1) / (MAX_RANSAC_POINTS - 1)])
            .collect()
    };

    let mut state = 0x9e37_79b9_7f4a_7c15_u64 ^ points.len() as u64;
    let mut next_index = || {
        state = state
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1_442_695_040_888_963_407);
        (state as usize) % sample.len()
    };

    let mut best = ordinary_least_squares(&sample);
    let mut residuals: Vec<f64> = sample
        .iter()
        .map(|&(x, y)| (y - (best.0 * x + best.1)).abs())
        .collect();
    let mut best_median = median(&mut residuals);

    for _ in 0..RANSAC_TRIALS {
        let a = sample[next_index()];
        let b = sample[next_index()];
        if (b.0 - a.0).abs() < MIN_X_SEPARATION {
            continue;
        }
        let slope = (b.1 - a.1) / (b.0 - a.0);
        if !slope.is_finite() || slope <= 0.0 {
            continue;
        }
        let candidate = (slope, a.1 - slope * a.0);
        let mut residuals: Vec<f64> = sample
            .iter()
            .map(|&(x, y)| (y - (candidate.0 * x + candidate.1)).abs())
            .collect();
        let candidate_median = median(&mut residuals);
        if candidate_median < best_median {
            best = candidate;
            best_median = candidate_median;
        }
    }

    let residuals: Vec<f64> = points
        .iter()
        .map(|&(x, y)| (y - (best.0 * x + best.1)).abs())
        .collect();
    let mut residual_distribution = residuals.clone();
    let residual_median = median(&mut residual_distribution);
    let mut deviations: Vec<f64> = residuals
        .iter()
        .map(|residual| (residual - residual_median).abs())
        .collect();
    let mad = median(&mut deviations);
    // A liberal adaptive cutoff preserves real curvature for the nonlinear stage.
    let cutoff = 0.02_f64
        .max(residual_median * 4.0)
        .max(residual_median + 4.0 * 1.4826 * mad);
    let inliers = residuals
        .iter()
        .enumerate()
        .filter_map(|(idx, &residual)| (residual <= cutoff).then_some(idx))
        .collect();

    (best, inliers)
}

/// Pool-adjacent-violators isotonic regression for monotone knot targets.
fn isotonic(values: &[f64], weights: &[usize]) -> Vec<f64> {
    #[derive(Clone)]
    struct Block {
        start: usize,
        end: usize,
        weighted_sum: f64,
        weight: usize,
    }

    let mut blocks: Vec<Block> = Vec::with_capacity(values.len());
    for (idx, (&value, &weight)) in values.iter().zip(weights).enumerate() {
        blocks.push(Block {
            start: idx,
            end: idx + 1,
            weighted_sum: value * weight as f64,
            weight,
        });
        while blocks.len() >= 2 {
            let n = blocks.len();
            let left_mean = blocks[n - 2].weighted_sum / blocks[n - 2].weight as f64;
            let right_mean = blocks[n - 1].weighted_sum / blocks[n - 1].weight as f64;
            if left_mean <= right_mean {
                break;
            }
            let right = blocks.pop().unwrap();
            let left = blocks.last_mut().unwrap();
            left.end = right.end;
            left.weighted_sum += right.weighted_sum;
            left.weight += right.weight;
        }
    }

    let mut fitted = vec![0.0; values.len()];
    for block in blocks {
        let mean = block.weighted_sum / block.weight as f64;
        fitted[block.start..block.end].fill(mean);
    }
    fitted
}

fn nonlinear_knots(points: &[(f64, f64)], inliers: &[usize]) -> Vec<(f32, f32)> {
    if inliers.len() < MIN_NONLINEAR_LANDMARKS {
        return Vec::new();
    }

    let mut clean: Vec<_> = inliers.iter().map(|&idx| points[idx]).collect();
    clean.sort_unstable_by(|a, b| a.0.total_cmp(&b.0));
    if clean.last().unwrap().0 - clean.first().unwrap().0 < MIN_NONLINEAR_SPAN {
        return Vec::new();
    }
    let bins = (clean.len() / 8).clamp(4, 24);
    let mut xs = Vec::with_capacity(bins);
    let mut ys = Vec::with_capacity(bins);
    let mut weights = Vec::with_capacity(bins);

    for bin in 0..bins {
        let start = bin * clean.len() / bins;
        let end = (bin + 1) * clean.len() / bins;
        let mut x: Vec<_> = clean[start..end].iter().map(|point| point.0).collect();
        let mut y: Vec<_> = clean[start..end].iter().map(|point| point.1).collect();
        xs.push(median(&mut x));
        ys.push(median(&mut y));
        weights.push(end - start);
    }

    let ys = isotonic(&ys, &weights);
    xs.into_iter()
        .zip(ys)
        .filter(|(x, y)| x.is_finite() && y.is_finite())
        .map(|(x, y)| (x as f32, y as f32))
        .collect()
}

fn fit_alignment(file_id: usize, max_rt: f64, mut points: Vec<(f64, f64)>) -> Alignment {
    points.retain(|(x, y)| x.is_finite() && y.is_finite());
    points.sort_unstable_by(|a, b| a.0.total_cmp(&b.0));
    let (robust, inliers) = robust_affine_inliers(&points);
    let affine_points: Vec<_> = inliers.iter().map(|&idx| points[idx]).collect();
    let affine = if affine_points.len() >= 2 {
        ordinary_least_squares(&affine_points)
    } else {
        robust
    };
    let knots = nonlinear_knots(&points, &inliers);

    Alignment {
        file_id,
        max_rt: max_rt as f32,
        slope: affine.0 as f32,
        intercept: affine.1 as f32,
        knots,
    }
}

/// Align runs with Sage's existing linear method.
///
/// Kept as the default entry point for API compatibility. Use
/// [`global_alignment_with_method`] to request nonlinear alignment.
pub fn global_alignment(features: &mut [Feature], n_files: usize) -> Vec<Alignment> {
    global_alignment_with_method(features, n_files, AlignmentMethod::Linear)
}

pub fn global_alignment_with_method(
    features: &mut [Feature],
    n_files: usize,
    method: AlignmentMethod,
) -> Vec<Alignment> {
    let max_rt = max_rt_by_file(features, n_files);
    let rt = rt_matrix(features, &max_rt, method == AlignmentMethod::Nonlinear);

    // Nonlinear mode uses a robust median consensus; linear mode retains the legacy mean.
    let mean_rts: Vec<f64> = (0..rt.rows)
        .into_par_iter()
        .map(|row| {
            let mut values: Vec<_> = rt.row(row).filter(|rt| rt.is_finite()).collect();
            match method {
                AlignmentMethod::Linear => values.iter().sum::<f64>() / values.len() as f64,
                AlignmentMethod::Nonlinear => median(&mut values),
            }
        })
        .collect();

    let alignments = (0..n_files)
        .into_par_iter()
        .map(|file_id| {
            let points: Vec<_> = rt
                .col(file_id)
                .zip(mean_rts.iter().copied())
                .filter(|(x, y)| x.is_finite() && y.is_finite())
                .collect();
            let alignment = match method {
                AlignmentMethod::Linear => {
                    legacy_linear_alignment(file_id, max_rt[file_id], &points)
                }
                AlignmentMethod::Nonlinear => fit_alignment(file_id, max_rt[file_id], points),
            };
            match method {
                AlignmentMethod::Linear => log::info!(
                    "aligning file #{file}: y = {m:.4}x + {b:.4}",
                    file = file_id,
                    m = alignment.slope,
                    b = alignment.intercept,
                ),
                AlignmentMethod::Nonlinear => log::info!(
                    "aligning file #{file}: robust affine y = {m:.4}x + {b:.4}, nonlinear knots = {knots}",
                    file = file_id,
                    m = alignment.slope,
                    b = alignment.intercept,
                    knots = alignment.knots.len(),
                ),
            }
            alignment
        })
        .collect::<Vec<_>>();

    log::info!("aligned retention times across {} files", n_files);

    features.par_iter_mut().for_each(|feature| {
        feature.aligned_rt = alignments[feature.file_id].transform(feature.rt);
    });

    alignments
}

#[cfg(test)]
mod tests {
    use super::*;

    fn transform_normalized(alignment: &Alignment, x: f64) -> f64 {
        alignment.transform((x * alignment.max_rt as f64) as f32) as f64
    }

    #[test]
    fn nonlinear_warp_beats_affine_fit_with_outliers() {
        let points: Vec<_> = (0..160)
            .map(|i| {
                let x = i as f64 / 159.0;
                let expected = 0.04 + 0.9 * x + 0.055 * (std::f64::consts::PI * x).sin();
                let y = if i % 29 == 0 {
                    1.0 - expected
                } else {
                    expected
                };
                (x, y)
            })
            .collect();
        let alignment = fit_alignment(0, 100.0, points);
        assert!(alignment.knots.len() >= 4);

        let (nonlinear_error, affine_error) = (0..100).fold((0.0, 0.0), |(nl, linear), i| {
            let x = i as f64 / 99.0;
            let expected = 0.04 + 0.9 * x + 0.055 * (std::f64::consts::PI * x).sin();
            (
                nl + (transform_normalized(&alignment, x) - expected).abs(),
                linear + (alignment.slope as f64 * x + alignment.intercept as f64 - expected).abs(),
            )
        });
        assert!(nonlinear_error < affine_error * 0.45);
    }

    #[test]
    fn fitted_warp_is_monotone() {
        let points: Vec<_> = (0..80)
            .map(|i| {
                let x = i as f64 / 79.0;
                let noise = if i % 7 == 0 { -0.04 } else { 0.0 };
                (x, x.powf(1.35) + noise)
            })
            .collect();
        let alignment = fit_alignment(0, 60.0, points);
        let transformed: Vec<_> = (0..200)
            .map(|i| transform_normalized(&alignment, i as f64 / 199.0))
            .collect();
        assert!(transformed.windows(2).all(|pair| pair[0] <= pair[1]));
    }

    #[test]
    fn sparse_landmarks_use_affine_fallback() {
        let points = vec![(0.1, 0.25), (0.4, 0.46), (0.8, 0.74), (0.95, 0.845)];
        let alignment = fit_alignment(2, 80.0, points);
        assert!(alignment.knots.is_empty());
        assert!((transform_normalized(&alignment, 0.6) - 0.6).abs() < 0.01);
    }

    #[test]
    fn robust_fit_handles_large_shift_and_bad_landmarks() {
        let points: Vec<_> = (0..100)
            .map(|i| {
                let x = i as f64 / 99.0;
                let expected = 0.2 + 0.65 * x;
                (
                    x,
                    if i % 17 == 0 {
                        1.0 - expected
                    } else {
                        expected
                    },
                )
            })
            .collect();
        let alignment = fit_alignment(0, 120.0, points);
        assert!((transform_normalized(&alignment, 0.5) - 0.525).abs() < 0.01);
    }

    #[test]
    fn alignment_method_is_explicit_and_defaults_to_linear() {
        assert_eq!(AlignmentMethod::default(), AlignmentMethod::Linear);
    }
}
