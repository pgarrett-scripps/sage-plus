//! Robust models for correcting systematic mass error.
//!
//! The sign convention used here is `observed - theoretical`, in ppm. A model
//! prediction is therefore removed from an observed mass before searching.

use serde::{Deserialize, Serialize};

/// A mass-error observation from a confidently identified PSM.
#[derive(Copy, Clone, Debug, PartialEq)]
pub struct CalibrationPoint {
    pub rt_minutes: f32,
    pub error_ppm: f32,
}

/// Controls when an RT-dependent model is allowed to replace a static offset.
#[derive(Copy, Clone, Debug, PartialEq)]
pub struct FitOptions {
    /// Minimum finite observations required for any model.
    pub min_points: usize,
    /// Minimum RT range represented by the inliers.
    pub min_rt_span: f32,
    /// Required fractional reduction in median absolute residual.
    pub min_linear_improvement: f32,
    /// Robust outlier cutoff, expressed as scaled median absolute deviations.
    pub outlier_mads: f32,
}

impl Default for FitOptions {
    fn default() -> Self {
        Self {
            min_points: 50,
            min_rt_span: 10.0,
            min_linear_improvement: 0.15,
            outlier_mads: 4.0,
        }
    }
}

#[derive(Copy, Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum CalibrationKind {
    Static,
    RetentionTimeLinear,
}

/// `error_ppm = intercept_ppm + slope_ppm_per_min * (rt - rt_center)`.
#[derive(Copy, Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct CalibrationModel {
    pub kind: CalibrationKind,
    pub intercept_ppm: f32,
    pub slope_ppm_per_min: f32,
    pub rt_center: f32,
}

impl CalibrationModel {
    #[inline]
    pub fn predict_ppm(&self, rt_minutes: f32) -> f32 {
        self.intercept_ppm + self.slope_ppm_per_min * (rt_minutes - self.rt_center)
    }

    /// Remove the predicted systematic error from an observed mass.
    #[inline]
    pub fn correct_observed(&self, observed: f32, rt_minutes: f32) -> f32 {
        observed / (1.0 + self.predict_ppm(rt_minutes) * 1e-6)
    }
}

#[derive(Copy, Clone, Debug, PartialEq)]
pub struct CalibrationFit {
    pub model: CalibrationModel,
    pub points: usize,
    pub inliers: usize,
    pub static_median_abs_residual: f32,
    pub model_median_abs_residual: f32,
}

/// Signed ppm error using the module-wide `observed - theoretical` convention.
#[inline]
pub fn ppm_error(observed: f32, theoretical: f32) -> f32 {
    (observed - theoretical) * 1e6 / theoretical
}

/// Translate an intensity-weighted absolute fragment error by a predicted
/// signed center without retaining every matched ion error. This is exact when
/// matched errors share a sign and preserves the original within-PSM spread.
#[inline]
pub fn align_fragment_error(average_abs_ppm: f32, signed_mean_ppm: f32, predicted_ppm: f32) -> f32 {
    (average_abs_ppm + (signed_mean_ppm - predicted_ppm).abs() - signed_mean_ppm.abs()).max(0.0)
}

/// Fit a robust static offset, promoting it to an RT-linear model only when
/// the latter has adequate RT coverage and materially improves residual error.
pub fn fit(points: &[CalibrationPoint], options: FitOptions) -> Option<CalibrationFit> {
    let finite = points
        .iter()
        .copied()
        .filter(|p| p.rt_minutes.is_finite() && p.error_ppm.is_finite())
        .collect::<Vec<_>>();
    if finite.len() < options.min_points {
        return None;
    }

    let errors = finite.iter().map(|p| p.error_ppm).collect::<Vec<_>>();
    let static_offset = median(&errors);
    let deviations = errors
        .iter()
        .map(|error| (error - static_offset).abs())
        .collect::<Vec<_>>();
    let mad = median(&deviations);
    // A small floor keeps a near-perfect run from rejecting harmless rounding
    // noise while still removing grossly incorrect PSMs.
    let cutoff = (options.outlier_mads * 1.4826 * mad).max(0.25);
    let inliers = finite
        .iter()
        .copied()
        .filter(|p| (p.error_ppm - static_offset).abs() <= cutoff)
        .collect::<Vec<_>>();
    if inliers.len() < options.min_points {
        return None;
    }

    let rts = inliers.iter().map(|p| p.rt_minutes).collect::<Vec<_>>();
    let rt_center = median(&rts);
    let static_residual = inliers
        .iter()
        .map(|p| (p.error_ppm - static_offset).abs())
        .collect::<Vec<_>>();
    let static_mar = median(&static_residual);

    let mean_x = inliers
        .iter()
        .map(|p| p.rt_minutes - rt_center)
        .sum::<f32>()
        / inliers.len() as f32;
    let mean_y = inliers.iter().map(|p| p.error_ppm).sum::<f32>() / inliers.len() as f32;
    let (covariance, variance) = inliers.iter().fold((0.0, 0.0), |(cov, var), p| {
        let x = p.rt_minutes - rt_center - mean_x;
        (cov + x * (p.error_ppm - mean_y), var + x * x)
    });
    let slope = if variance > f32::EPSILON {
        covariance / variance
    } else {
        0.0
    };
    let intercept = mean_y - slope * mean_x;
    let linear_residual = inliers
        .iter()
        .map(|p| (p.error_ppm - intercept - slope * (p.rt_minutes - rt_center)).abs())
        .collect::<Vec<_>>();
    let linear_mar = median(&linear_residual);
    let rt_span = rts.iter().copied().fold(f32::NEG_INFINITY, f32::max)
        - rts.iter().copied().fold(f32::INFINITY, f32::min);
    let improvement = if static_mar > f32::EPSILON {
        1.0 - linear_mar / static_mar
    } else {
        0.0
    };

    let use_linear = rt_span >= options.min_rt_span
        && improvement >= options.min_linear_improvement
        && slope.is_finite();
    let model = if use_linear {
        CalibrationModel {
            kind: CalibrationKind::RetentionTimeLinear,
            intercept_ppm: intercept,
            slope_ppm_per_min: slope,
            rt_center,
        }
    } else {
        CalibrationModel {
            kind: CalibrationKind::Static,
            intercept_ppm: static_offset,
            slope_ppm_per_min: 0.0,
            rt_center,
        }
    };

    Some(CalibrationFit {
        model,
        points: finite.len(),
        inliers: inliers.len(),
        static_median_abs_residual: static_mar,
        model_median_abs_residual: if use_linear { linear_mar } else { static_mar },
    })
}

fn median(values: &[f32]) -> f32 {
    debug_assert!(!values.is_empty());
    let mut values = values.to_vec();
    values.sort_unstable_by(f32::total_cmp);
    let middle = values.len() / 2;
    if values.len() % 2 == 0 {
        (values[middle - 1] + values[middle]) / 2.0
    } else {
        values[middle]
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn options() -> FitOptions {
        FitOptions {
            min_points: 10,
            min_rt_span: 5.0,
            min_linear_improvement: 0.15,
            outlier_mads: 4.0,
        }
    }

    #[test]
    fn fits_static_offset() {
        let points = (0..40)
            .map(|i| CalibrationPoint {
                rt_minutes: i as f32,
                error_ppm: 3.0 + (i % 3) as f32 * 0.05,
            })
            .collect::<Vec<_>>();
        let fit = fit(&points, options()).unwrap();
        assert_eq!(fit.model.kind, CalibrationKind::Static);
        assert!((fit.model.intercept_ppm - 3.05).abs() < 0.06);
    }

    #[test]
    fn selects_material_rt_drift() {
        let points = (0..60)
            .map(|i| {
                let rt = i as f32;
                CalibrationPoint {
                    rt_minutes: rt,
                    error_ppm: -2.0 + 0.12 * (rt - 30.0) + (i % 2) as f32 * 0.02,
                }
            })
            .collect::<Vec<_>>();
        let fit = fit(&points, options()).unwrap();
        assert_eq!(fit.model.kind, CalibrationKind::RetentionTimeLinear);
        assert!((fit.model.slope_ppm_per_min - 0.12).abs() < 0.002);
        assert!(fit.model_median_abs_residual < fit.static_median_abs_residual);
    }

    #[test]
    fn ignores_gross_outliers() {
        let mut points = (0..40)
            .map(|i| CalibrationPoint {
                rt_minutes: i as f32,
                error_ppm: 2.0 + (i % 2) as f32 * 0.1,
            })
            .collect::<Vec<_>>();
        points.extend((0..5).map(|i| CalibrationPoint {
            rt_minutes: i as f32,
            error_ppm: 100.0,
        }));
        let fit = fit(&points, options()).unwrap();
        assert_eq!(fit.inliers, 40);
        assert!((fit.model.intercept_ppm - 2.05).abs() < 0.06);
    }

    #[test]
    fn correction_removes_predicted_error() {
        let model = CalibrationModel {
            kind: CalibrationKind::RetentionTimeLinear,
            intercept_ppm: 2.0,
            slope_ppm_per_min: 0.1,
            rt_center: 20.0,
        };
        let theoretical = 1000.0;
        let observed = theoretical * (1.0 + model.predict_ppm(30.0) * 1e-6);
        assert!((model.correct_observed(observed, 30.0) - theoretical).abs() < 1e-4);
    }

    #[test]
    fn fragment_alignment_preserves_spread() {
        assert_eq!(align_fragment_error(5.0, 5.0, 5.0), 0.0);
        assert_eq!(align_fragment_error(3.0, 2.0, 2.0), 1.0);
        assert_eq!(align_fragment_error(3.0, 2.0, 0.0), 3.0);
    }
}
