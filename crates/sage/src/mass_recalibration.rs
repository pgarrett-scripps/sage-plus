//! Validation-gated models of systematic precursor and fragment mass error
//! used to recalibrate spectra before a second search pass.
//!
//! Errors follow the `observed - theoretical` convention in ppm, as in
//! [`crate::mass_calibration`]. A model predicts that error from retention
//! time and/or m/z. Its prediction is removed from observed m/z values with
//! `corrected = observed / (1 + ppm * 1e-6)`.
//!
//! The model family is deliberately small and additive:
//!
//! 1. `none`: no correction;
//! 2. `static`: one offset;
//! 3. `linear`: an offset plus a linear term in RT, m/z, or both;
//! 4. `smooth`: an offset plus a penalized piecewise-linear curve in RT, m/z,
//!    or both, with at most [`RecalibrationOptions::max_intervals`] intervals
//!    per axis. There is no RT x m/z interaction surface.
//!
//! PSMs, not individual fragment ions, are split deterministically into fit
//! and validation sets. A more complex model replaces a simpler one only when
//! its median absolute validation residual improves by a set fraction and no
//! RT or m/z tercile of the validation set gets worse.

use crate::mass::PROTON;
use crate::ml::regression::cholesky_solve;
use crate::spectrum::ProcessedSpectrum;
use serde::{Deserialize, Serialize};

/// How much search-time mass recalibration a run may apply.
#[derive(
    Copy, Clone, Debug, Default, Deserialize, Eq, PartialEq, Serialize, schemars::JsonSchema,
)]
#[serde(rename_all = "snake_case")]
pub enum MassRecalibrationMode {
    /// Search once with observed masses (post-hoc alignment still feeds the
    /// discriminant model).
    #[default]
    Off,
    /// Correct each file by at most a validated constant offset.
    Static,
    /// Allow validated offsets and linear RT/m/z terms.
    Linear,
    /// Allow validated offsets, linear terms, and smooth additive RT/m/z
    /// curves.
    Auto,
}

impl MassRecalibrationMode {
    pub fn max_kind(self) -> MassModelKind {
        match self {
            Self::Off => MassModelKind::None,
            Self::Static => MassModelKind::Static,
            Self::Linear => MassModelKind::Linear,
            Self::Auto => MassModelKind::Smooth,
        }
    }

    pub fn as_str(self) -> &'static str {
        match self {
            Self::Off => "off",
            Self::Static => "static",
            Self::Linear => "linear",
            Self::Auto => "auto",
        }
    }
}

/// Model complexity, ordered from simplest to most flexible.
#[derive(Copy, Clone, Debug, Default, PartialEq, Eq, PartialOrd, Ord, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum MassModelKind {
    #[default]
    None,
    Static,
    Linear,
    Smooth,
}

/// Axes a linear or smooth model depends on.
#[derive(Copy, Clone, Debug, Default, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum MassModelAxes {
    #[default]
    None,
    Rt,
    Mz,
    RtMz,
}

impl MassModelAxes {
    fn uses_rt(self) -> bool {
        matches!(self, Self::Rt | Self::RtMz)
    }

    fn uses_mz(self) -> bool {
        matches!(self, Self::Mz | Self::RtMz)
    }
}

/// One-dimensional term of an additive model.
///
/// Linear terms contribute `slope * (clamp(x) - center)`. Smooth terms
/// interpolate `values` at `knots`. Inputs are clamped to `[min, max]`, the
/// range covered by the fitting data, so the model never extrapolates.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct AxisTerm {
    pub min: f32,
    pub max: f32,
    pub center: f32,
    /// ppm per axis unit (minute or m/z); zero for smooth terms.
    pub slope: f32,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub knots: Vec<f32>,
    /// ppm at each knot; empty for linear terms.
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub values: Vec<f32>,
}

impl AxisTerm {
    #[inline]
    fn predict(&self, x: f32) -> f32 {
        let x = x.clamp(self.min, self.max);
        if self.knots.len() < 2 {
            return self.slope * (x - self.center);
        }
        let (j, w) = hat_position(&self.knots, x);
        self.values[j] * (1.0 - w) + self.values[j + 1] * w
    }
}

/// An additive model of mass error in ppm.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct MassErrorModel {
    pub kind: MassModelKind,
    pub axes: MassModelAxes,
    pub intercept_ppm: f32,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub rt: Option<AxisTerm>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub mz: Option<AxisTerm>,
    /// Predictions are clamped to `[-max_abs_ppm, max_abs_ppm]`.
    pub max_abs_ppm: f32,
}

impl MassErrorModel {
    pub fn none() -> Self {
        Self {
            kind: MassModelKind::None,
            axes: MassModelAxes::None,
            intercept_ppm: 0.0,
            rt: None,
            mz: None,
            max_abs_ppm: f32::INFINITY,
        }
    }

    /// Predicted `observed - theoretical` error, in ppm.
    #[inline]
    pub fn predict_ppm(&self, rt_minutes: f32, mz: f32) -> f32 {
        let mut ppm = self.intercept_ppm;
        if let Some(term) = &self.rt {
            ppm += term.predict(rt_minutes);
        }
        if let Some(term) = &self.mz {
            ppm += term.predict(mz);
        }
        if ppm.is_finite() {
            ppm.clamp(-self.max_abs_ppm, self.max_abs_ppm)
        } else {
            0.0
        }
    }

    /// Remove the predicted error from an observed m/z.
    #[inline]
    pub fn correct_mz(&self, observed_mz: f32, rt_minutes: f32) -> f32 {
        observed_mz / (1.0 + self.predict_ppm(rt_minutes, observed_mz) * 1e-6)
    }

    pub fn is_identity(&self) -> bool {
        self.kind == MassModelKind::None
    }

    /// Number of free parameters, used for reporting.
    pub fn parameters(&self) -> usize {
        let axis = |term: &Option<AxisTerm>| match term {
            Some(term) if term.knots.len() >= 2 => term.knots.len(),
            Some(_) => 1,
            None => 0,
        };
        match self.kind {
            MassModelKind::None => 0,
            _ => 1 + axis(&self.rt) + axis(&self.mz),
        }
    }
}

/// Per-file corrections applied during the second search pass.
#[derive(Clone, Debug, Default, PartialEq, Serialize, Deserialize)]
pub struct FileMassCorrection {
    /// Correction applied to precursor m/z.
    pub precursor: Option<MassErrorModel>,
    /// Correction applied to fragment m/z, as a function of fragment m/z.
    pub fragment: Option<MassErrorModel>,
}

impl FileMassCorrection {
    pub fn is_identity(&self) -> bool {
        self.precursor
            .as_ref()
            .is_none_or(MassErrorModel::is_identity)
            && self
                .fragment
                .as_ref()
                .is_none_or(MassErrorModel::is_identity)
    }

    #[inline]
    pub fn precursor_ppm(&self, rt_minutes: f32, mz: f32) -> f32 {
        self.precursor
            .as_ref()
            .map_or(0.0, |model| model.predict_ppm(rt_minutes, mz))
    }

    #[inline]
    pub fn fragment_ppm(&self, rt_minutes: f32, mz: f32) -> f32 {
        self.fragment
            .as_ref()
            .map_or(0.0, |model| model.predict_ppm(rt_minutes, mz))
    }
}

/// Search-time corrections for every input file, indexed by `file_id`.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct MassRecalibration {
    pub files: Vec<FileMassCorrection>,
}

impl MassRecalibration {
    /// The correction for a file, or `None` when nothing would change.
    #[inline]
    pub fn file(&self, file_id: usize) -> Option<&FileMassCorrection> {
        self.files
            .get(file_id)
            .filter(|correction| !correction.is_identity())
    }

    /// A copy of `spectrum` with corrected precursor and fragment m/z, or
    /// `None` when its file has no correction. Peak order, intensities and
    /// charges are unchanged; a fragment correction is far smaller than peak
    /// spacing, so ascending order is preserved up to ties.
    pub fn recalibrate(&self, spectrum: &ProcessedSpectrum) -> Option<ProcessedSpectrum> {
        let correction = self.file(spectrum.file_id)?;
        let rt = spectrum.scan_start_time;
        let mut out = spectrum.clone();
        if let Some(model) = correction.precursor.as_ref().filter(|m| !m.is_identity()) {
            for precursor in &mut out.precursors {
                precursor.mz = model.correct_mz(precursor.mz, rt);
            }
        }
        if let Some(model) = correction.fragment.as_ref().filter(|m| !m.is_identity()) {
            for (mass, &charge) in out.masses.iter_mut().zip(&spectrum.charges) {
                let z = charge.max(1) as f32;
                let mz = *mass / z + PROTON;
                *mass = (model.correct_mz(mz, rt) - PROTON) * z;
            }
            // Corrections that vary with m/z could, in principle, swap two
            // nearly coincident peaks; restore ascending order if so.
            if out.masses.windows(2).any(|w| w[0] > w[1]) {
                sort_peaks(&mut out);
            }
        }
        Some(out)
    }
}

fn sort_peaks(spectrum: &mut ProcessedSpectrum) {
    let mut order = (0..spectrum.masses.len()).collect::<Vec<_>>();
    order.sort_by(|&a, &b| spectrum.masses[a].total_cmp(&spectrum.masses[b]));
    fn permute<T: Copy>(values: &mut Vec<T>, order: &[usize]) {
        if values.len() == order.len() {
            *values = order.iter().map(|&i| values[i]).collect();
        }
    }
    permute(&mut spectrum.masses, &order);
    permute(&mut spectrum.intensities, &order);
    permute(&mut spectrum.charges, &order);
    permute(&mut spectrum.charge_is_known, &order);
    permute(&mut spectrum.mobilities, &order);
}

/// A mass-error observation. `group` identifies the PSM that produced it, so
/// all fragment ions of a PSM land on the same side of the validation split.
#[derive(Copy, Clone, Debug, PartialEq)]
pub struct MassErrorPoint {
    pub rt_minutes: f32,
    pub mz: f32,
    pub error_ppm: f32,
    pub group: u64,
}

#[derive(Copy, Clone, Debug, PartialEq)]
pub struct RecalibrationOptions {
    /// Most flexible model family that may be selected.
    pub max_kind: MassModelKind,
    /// Largest correction the model may apply, normally the search tolerance.
    pub max_abs_ppm: f32,
    pub min_static_psms: usize,
    pub min_linear_psms: usize,
    pub min_smooth_psms: usize,
    /// Minimum fitting PSMs per smooth interval.
    pub min_psms_per_interval: usize,
    /// Maximum smooth intervals per axis (knots = intervals + 1).
    pub max_intervals: usize,
    pub min_rt_span: f32,
    pub min_mz_span: f32,
    /// Validation PSMs are those with `hash(group) % 10 < validation_tenths`.
    pub validation_tenths: u64,
    /// Fractional validation improvement needed to accept a static offset.
    pub min_static_improvement: f32,
    /// Fractional validation improvement needed to step up to linear or smooth.
    pub min_step_improvement: f32,
    /// A tercile may not get worse than this fraction.
    pub max_bin_regression: f32,
    pub outlier_mads: f32,
}

impl Default for RecalibrationOptions {
    fn default() -> Self {
        Self {
            max_kind: MassModelKind::Smooth,
            max_abs_ppm: 20.0,
            min_static_psms: 50,
            min_linear_psms: 200,
            min_smooth_psms: 1000,
            min_psms_per_interval: 50,
            max_intervals: 5,
            min_rt_span: 10.0,
            min_mz_span: 100.0,
            validation_tenths: 3,
            min_static_improvement: 0.02,
            min_step_improvement: 0.05,
            max_bin_regression: 0.02,
            outlier_mads: 4.0,
        }
    }
}

/// Validation result for one candidate model.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct CandidateScore {
    pub kind: MassModelKind,
    pub axes: MassModelAxes,
    pub validation_median_abs_ppm: f32,
    pub accepted: bool,
}

/// Residual summary of the validation set within one RT or m/z bin.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct ResidualBin {
    /// `rt` or `mz`.
    pub axis: String,
    pub lower: f32,
    pub upper: f32,
    pub points: usize,
    pub median_ppm_before: f32,
    pub mad_ppm_before: f32,
    pub median_ppm_after: f32,
    pub mad_ppm_after: f32,
}

/// The chosen model together with the evidence that selected it.
#[derive(Clone, Debug, Default, PartialEq, Serialize, Deserialize)]
pub struct ModelSelection {
    pub model: Option<MassErrorModel>,
    pub psms: usize,
    pub points: usize,
    pub fit_psms: usize,
    pub validation_psms: usize,
    pub inlier_points: usize,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub skipped: Option<String>,
    pub candidates: Vec<CandidateScore>,
    /// Validation residuals before (raw) and after (fit-set model).
    pub bins: Vec<ResidualBin>,
}

impl ModelSelection {
    fn skipped(psms: usize, points: usize, reason: &str) -> Self {
        Self {
            psms,
            points,
            skipped: Some(reason.to_string()),
            ..Default::default()
        }
    }

    pub fn kind(&self) -> MassModelKind {
        self.model.as_ref().map_or(MassModelKind::None, |m| m.kind)
    }
}

/// FNV-1a hash of a string, used to split PSMs deterministically.
pub fn stable_hash(value: &str) -> u64 {
    let mut hash = 0xcbf29ce484222325u64;
    for byte in value.bytes() {
        hash ^= byte as u64;
        hash = hash.wrapping_mul(0x100000001b3);
    }
    hash
}

#[inline]
fn is_validation(group: u64, tenths: u64) -> bool {
    // Mix once more so sequential integer groups are split evenly.
    let mut x = group ^ (group >> 33);
    x = x.wrapping_mul(0xff51afd7ed558ccd);
    x ^= x >> 33;
    x % 10 < tenths
}

/// Select and fit a model following the none -> static -> linear -> smooth
/// hierarchy described in the module documentation.
pub fn select_model(points: &[MassErrorPoint], options: RecalibrationOptions) -> ModelSelection {
    let finite = points
        .iter()
        .copied()
        .filter(|p| p.rt_minutes.is_finite() && p.mz.is_finite() && p.error_ppm.is_finite())
        .collect::<Vec<_>>();
    let psms = count_groups(&finite);
    if options.max_kind == MassModelKind::None {
        return ModelSelection::skipped(psms, finite.len(), "disabled");
    }
    if psms < options.min_static_psms {
        return ModelSelection::skipped(psms, finite.len(), "too_few_psms");
    }

    let (fit, validation): (Vec<_>, Vec<_>) = finite
        .iter()
        .copied()
        .partition(|p| !is_validation(p.group, options.validation_tenths));
    if fit.is_empty() || validation.is_empty() {
        return ModelSelection::skipped(psms, finite.len(), "empty_split");
    }

    // Remove grossly wrong observations around the fit-set median.
    let center = median(&fit.iter().map(|p| p.error_ppm).collect::<Vec<_>>());
    let mad = median(
        &fit.iter()
            .map(|p| (p.error_ppm - center).abs())
            .collect::<Vec<_>>(),
    );
    let cutoff = (options.outlier_mads * 1.4826 * mad).max(0.25);
    let keep = |p: &MassErrorPoint| (p.error_ppm - center).abs() <= cutoff;
    let fit = fit.into_iter().filter(keep).collect::<Vec<_>>();
    let validation = validation.into_iter().filter(keep).collect::<Vec<_>>();
    let fit_psms = count_groups(&fit);
    let validation_psms = count_groups(&validation);
    let mut selection = ModelSelection {
        model: None,
        psms,
        points: finite.len(),
        fit_psms,
        validation_psms,
        inlier_points: fit.len() + validation.len(),
        skipped: None,
        candidates: Vec::new(),
        bins: Vec::new(),
    };
    if fit_psms < options.min_static_psms || validation.is_empty() {
        selection.skipped = Some("too_few_inliers".into());
        return selection;
    }

    let rt_bins = tercile_edges(&validation.iter().map(|p| p.rt_minutes).collect::<Vec<_>>());
    let mz_bins = tercile_edges(&validation.iter().map(|p| p.mz).collect::<Vec<_>>());
    let score = |model: &MassErrorModel| evaluate(model, &validation, &rt_bins, &mz_bins);

    let mut current = MassErrorModel::none();
    let (mut current_mar, mut current_bins) = score(&current);
    selection.candidates.push(CandidateScore {
        kind: MassModelKind::None,
        axes: MassModelAxes::None,
        validation_median_abs_ppm: current_mar,
        accepted: true,
    });

    let mut consider = |candidates: Vec<MassErrorModel>,
                        min_improvement: f32,
                        selection: &mut ModelSelection,
                        current: &mut MassErrorModel| {
        let scored = candidates
            .into_iter()
            .map(|model| {
                let (mar, bins) = score(&model);
                (model, mar, bins)
            })
            .collect::<Vec<_>>();
        let best = scored
            .iter()
            .enumerate()
            .filter(|(_, (_, mar, _))| mar.is_finite())
            .min_by(|a, b| a.1 .1.total_cmp(&b.1 .1))
            .map(|(i, _)| i);
        let mut accepted = None;
        if let Some(i) = best {
            let (_, mar, bins) = &scored[i];
            let improves = *mar <= (1.0 - min_improvement) * current_mar;
            let no_bin_worse =
                bins.iter()
                    .zip(current_bins.iter())
                    .all(|(new, old)| match (new, old) {
                        (Some(new), Some(old)) => {
                            *new <= old * (1.0 + options.max_bin_regression) + 1e-6
                        }
                        _ => true,
                    });
            if improves && no_bin_worse {
                accepted = Some(i);
            }
        }
        for (i, (model, mar, _)) in scored.iter().enumerate() {
            selection.candidates.push(CandidateScore {
                kind: model.kind,
                axes: model.axes,
                validation_median_abs_ppm: *mar,
                accepted: accepted == Some(i),
            });
        }
        if let Some(i) = accepted {
            let (model, mar, bins) = scored.into_iter().nth(i).expect("index in range");
            *current = model;
            current_mar = mar;
            current_bins = bins;
        }
    };

    // Static offset.
    consider(
        vec![fit_static(&fit, options)],
        options.min_static_improvement,
        &mut selection,
        &mut current,
    );

    let rt_ok = span(fit.iter().map(|p| p.rt_minutes)) >= options.min_rt_span;
    let mz_ok = span(fit.iter().map(|p| p.mz)) >= options.min_mz_span;
    let axes = [
        (MassModelAxes::Rt, rt_ok),
        (MassModelAxes::Mz, mz_ok),
        (MassModelAxes::RtMz, rt_ok && mz_ok),
    ]
    .into_iter()
    .filter_map(|(axes, ok)| ok.then_some(axes))
    .collect::<Vec<_>>();

    if options.max_kind >= MassModelKind::Linear && fit_psms >= options.min_linear_psms {
        let candidates = axes
            .iter()
            .filter_map(|&axes| fit_additive(&fit, axes, 0, options))
            .collect::<Vec<_>>();
        consider(
            candidates,
            options.min_step_improvement,
            &mut selection,
            &mut current,
        );
    }

    let intervals = (fit_psms / options.min_psms_per_interval.max(1)).min(options.max_intervals);
    if options.max_kind >= MassModelKind::Smooth
        && fit_psms >= options.min_smooth_psms
        && intervals >= 2
    {
        let candidates = axes
            .iter()
            .filter_map(|&axes| fit_additive(&fit, axes, intervals, options))
            .collect::<Vec<_>>();
        consider(
            candidates,
            options.min_step_improvement,
            &mut selection,
            &mut current,
        );
    }

    // Report validation residuals of the fit-set model, then refit the chosen
    // form on every inlier.
    selection.bins = residual_bins(&current, &validation, &rt_bins, &mz_bins);
    if current.kind == MassModelKind::None {
        return selection;
    }
    let all = fit
        .iter()
        .chain(validation.iter())
        .copied()
        .collect::<Vec<_>>();
    let refit = match current.kind {
        MassModelKind::Static => Some(fit_static(&all, options)),
        MassModelKind::Linear => fit_additive(&all, current.axes, 0, options),
        MassModelKind::Smooth => {
            let knots = |term: &Option<AxisTerm>| term.as_ref().map(|t| t.knots.clone());
            fit_with_knots(
                &all,
                current.axes,
                knots(&current.rt),
                knots(&current.mz),
                options,
            )
        }
        MassModelKind::None => None,
    };
    selection.model = Some(refit.unwrap_or(current));
    selection
}

fn count_groups(points: &[MassErrorPoint]) -> usize {
    let mut groups = points.iter().map(|p| p.group).collect::<Vec<_>>();
    groups.sort_unstable();
    groups.dedup();
    groups.len()
}

fn span(values: impl Iterator<Item = f32>) -> f32 {
    let (lo, hi) = values.fold((f32::INFINITY, f32::NEG_INFINITY), |(lo, hi), x| {
        (lo.min(x), hi.max(x))
    });
    if hi >= lo {
        hi - lo
    } else {
        0.0
    }
}

fn fit_static(points: &[MassErrorPoint], options: RecalibrationOptions) -> MassErrorModel {
    let offset = median(&points.iter().map(|p| p.error_ppm).collect::<Vec<_>>());
    MassErrorModel {
        kind: MassModelKind::Static,
        axes: MassModelAxes::None,
        intercept_ppm: offset,
        rt: None,
        mz: None,
        max_abs_ppm: options.max_abs_ppm,
    }
}

/// Fit a linear (`intervals == 0`) or smooth additive model.
fn fit_additive(
    points: &[MassErrorPoint],
    axes: MassModelAxes,
    intervals: usize,
    options: RecalibrationOptions,
) -> Option<MassErrorModel> {
    if intervals == 0 {
        return fit_with_knots(points, axes, None, None, options);
    }
    let knots = |values: Vec<f32>| quantile_knots(values, intervals);
    let rt = axes
        .uses_rt()
        .then(|| knots(points.iter().map(|p| p.rt_minutes).collect()));
    let mz = axes
        .uses_mz()
        .then(|| knots(points.iter().map(|p| p.mz).collect()));
    // Too few distinct positions on an axis for a curve.
    if rt.as_ref().is_some_and(|k| k.len() < 3) || mz.as_ref().is_some_and(|k| k.len() < 3) {
        return None;
    }
    fit_with_knots(points, axes, rt, mz, options)
}

/// Penalized least squares for an additive model.
///
/// Linear axes (no knots) contribute one centered, scaled column. Smooth axes
/// use a piecewise-linear hat basis with a second-difference penalty on the
/// knot values, whose strong-penalty limit is a straight line. A small ridge
/// resolves the constant shared by the intercept and each hat basis.
fn fit_with_knots(
    points: &[MassErrorPoint],
    axes: MassModelAxes,
    rt_knots: Option<Vec<f32>>,
    mz_knots: Option<Vec<f32>>,
    options: RecalibrationOptions,
) -> Option<MassErrorModel> {
    if points.is_empty() || axes == MassModelAxes::None {
        return None;
    }
    let smooth = rt_knots.is_some() || mz_knots.is_some();

    struct Axis {
        min: f32,
        max: f32,
        center: f32,
        scale: f32,
        knots: Option<Vec<f32>>,
        offset: usize,
        width: usize,
    }
    let make_axis = |values: Vec<f32>, knots: Option<Vec<f32>>, offset: usize| {
        let min = values.iter().copied().fold(f32::INFINITY, f32::min);
        let max = values.iter().copied().fold(f32::NEG_INFINITY, f32::max);
        let center = median(&values);
        let scale = ((max - min) / 2.0).max(1e-3);
        let width = knots.as_ref().map_or(1, Vec::len);
        Axis {
            min,
            max,
            center,
            scale,
            knots,
            offset,
            width,
        }
    };

    let mut columns = 1;
    let rt = axes.uses_rt().then(|| {
        let axis = make_axis(
            points.iter().map(|p| p.rt_minutes).collect(),
            rt_knots.clone(),
            columns,
        );
        columns += axis.width;
        axis
    });
    let mz = axes.uses_mz().then(|| {
        let axis = make_axis(
            points.iter().map(|p| p.mz).collect(),
            mz_knots.clone(),
            columns,
        );
        columns += axis.width;
        axis
    });

    let p = columns;
    let mut a = vec![0f64; p * p];
    let mut b = vec![0f64; p];
    let mut row: Vec<(usize, f64)> = Vec::with_capacity(5);
    for point in points {
        row.clear();
        row.push((0, 1.0));
        for (axis, x) in [(&rt, point.rt_minutes), (&mz, point.mz)] {
            let Some(axis) = axis else { continue };
            let x = x.clamp(axis.min, axis.max);
            match &axis.knots {
                None => row.push((axis.offset, ((x - axis.center) / axis.scale) as f64)),
                Some(knots) => {
                    let (j, w) = hat_position(knots, x);
                    row.push((axis.offset + j, 1.0 - w as f64));
                    row.push((axis.offset + j + 1, w as f64));
                }
            }
        }
        let y = point.error_ppm as f64;
        for &(i, xi) in &row {
            b[i] += xi * y;
            for &(j, xj) in &row {
                a[i * p + j] += xi * xj;
            }
        }
    }

    let n = points.len() as f64;
    for axis in [&rt, &mz].into_iter().flatten() {
        match &axis.knots {
            None => {
                let i = axis.offset;
                a[i * p + i] += 1e-9 * n;
            }
            Some(knots) => {
                let intervals = (knots.len() - 1) as f64;
                let lambda = n / (4.0 * intervals);
                let ridge = 1e-6 * n / intervals;
                for j in 0..knots.len() {
                    let i = axis.offset + j;
                    a[i * p + i] += ridge;
                }
                // Divided second differences, normalized by mean spacing so
                // the penalty is dimensionless.
                let range = (knots[knots.len() - 1] - knots[0]).max(1e-6) as f64;
                for j in 1..knots.len() - 1 {
                    let h1 = ((knots[j] - knots[j - 1]) as f64 / range).max(1e-9);
                    let h2 = ((knots[j + 1] - knots[j]) as f64 / range).max(1e-9);
                    let hbar = (h1 + h2) / 2.0;
                    let d = [
                        (axis.offset + j - 1, hbar / h1),
                        (axis.offset + j, -hbar * (1.0 / h1 + 1.0 / h2)),
                        (axis.offset + j + 1, hbar / h2),
                    ];
                    for &(r, dr) in &d {
                        for &(c, dc) in &d {
                            a[r * p + c] += lambda * dr * dc;
                        }
                    }
                }
            }
        }
    }

    let beta = cholesky_solve(a, b)?;
    let mut intercept = beta[0];
    let mut term = |axis: &Option<Axis>| -> Option<AxisTerm> {
        let axis = axis.as_ref()?;
        Some(match &axis.knots {
            None => AxisTerm {
                min: axis.min,
                max: axis.max,
                center: axis.center,
                slope: (beta[axis.offset] / axis.scale as f64) as f32,
                knots: Vec::new(),
                values: Vec::new(),
            },
            Some(knots) => {
                let values = &beta[axis.offset..axis.offset + knots.len()];
                // Move each curve's mean into the intercept for readability.
                let mean = values.iter().sum::<f64>() / values.len() as f64;
                intercept += mean;
                AxisTerm {
                    min: axis.min,
                    max: axis.max,
                    center: axis.center,
                    slope: 0.0,
                    knots: knots.clone(),
                    values: values.iter().map(|v| (v - mean) as f32).collect(),
                }
            }
        })
    };
    let rt_term = term(&rt);
    let mz_term = term(&mz);
    let model = MassErrorModel {
        kind: if smooth {
            MassModelKind::Smooth
        } else {
            MassModelKind::Linear
        },
        axes,
        intercept_ppm: intercept as f32,
        rt: rt_term,
        mz: mz_term,
        max_abs_ppm: options.max_abs_ppm,
    };
    model.intercept_ppm.is_finite().then_some(model)
}

/// Knots at evenly spaced quantiles, deduplicated.
fn quantile_knots(mut values: Vec<f32>, intervals: usize) -> Vec<f32> {
    values.sort_unstable_by(f32::total_cmp);
    let last = values.len() - 1;
    let mut knots = (0..=intervals)
        .map(|i| values[(i * last + intervals / 2) / intervals.max(1)])
        .collect::<Vec<_>>();
    knots.dedup_by(|a, b| (*a - *b).abs() <= f32::EPSILON * b.abs().max(1.0));
    knots
}

/// Interval index and interpolation weight of `x` within sorted `knots`.
#[inline]
fn hat_position(knots: &[f32], x: f32) -> (usize, f32) {
    let last = knots.len() - 1;
    let j = knots
        .partition_point(|&k| k <= x)
        .saturating_sub(1)
        .min(last - 1);
    let width = knots[j + 1] - knots[j];
    let w = if width > 0.0 {
        ((x - knots[j]) / width).clamp(0.0, 1.0)
    } else {
        0.0
    };
    (j, w)
}

fn tercile_edges(values: &[f32]) -> [f32; 2] {
    if values.is_empty() {
        return [f32::NAN; 2];
    }
    let mut values = values.to_vec();
    values.sort_unstable_by(f32::total_cmp);
    let n = values.len();
    [values[n / 3], values[(2 * n) / 3]]
}

#[inline]
fn tercile(edges: &[f32; 2], x: f32) -> usize {
    if x < edges[0] {
        0
    } else if x < edges[1] {
        1
    } else {
        2
    }
}

/// Overall and per-tercile (3 RT, then 3 m/z) median absolute residuals.
fn evaluate(
    model: &MassErrorModel,
    points: &[MassErrorPoint],
    rt_bins: &[f32; 2],
    mz_bins: &[f32; 2],
) -> (f32, Vec<Option<f32>>) {
    let mut all = Vec::with_capacity(points.len());
    let mut bins = vec![Vec::new(); 6];
    for p in points {
        let r = (p.error_ppm - model.predict_ppm(p.rt_minutes, p.mz)).abs();
        all.push(r);
        bins[tercile(rt_bins, p.rt_minutes)].push(r);
        bins[3 + tercile(mz_bins, p.mz)].push(r);
    }
    let overall = if all.is_empty() {
        f32::NAN
    } else {
        median(&all)
    };
    let bins = bins
        .into_iter()
        .map(|values| (values.len() >= 20).then(|| median(&values)))
        .collect();
    (overall, bins)
}

fn residual_bins(
    model: &MassErrorModel,
    points: &[MassErrorPoint],
    rt_bins: &[f32; 2],
    mz_bins: &[f32; 2],
) -> Vec<ResidualBin> {
    let mut out = Vec::with_capacity(6);
    for (axis, edges) in [("rt", rt_bins), ("mz", mz_bins)] {
        for bin in 0..3 {
            let members = points
                .iter()
                .filter(|p| {
                    let x = if axis == "rt" { p.rt_minutes } else { p.mz };
                    tercile(edges, x) == bin
                })
                .collect::<Vec<_>>();
            if members.is_empty() {
                continue;
            }
            let x = |p: &MassErrorPoint| if axis == "rt" { p.rt_minutes } else { p.mz };
            let before = members.iter().map(|p| p.error_ppm).collect::<Vec<_>>();
            let after = members
                .iter()
                .map(|p| p.error_ppm - model.predict_ppm(p.rt_minutes, p.mz))
                .collect::<Vec<_>>();
            let (median_before, mad_before) = median_mad(&before);
            let (median_after, mad_after) = median_mad(&after);
            out.push(ResidualBin {
                axis: axis.to_string(),
                lower: members.iter().map(|p| x(p)).fold(f32::INFINITY, f32::min),
                upper: members
                    .iter()
                    .map(|p| x(p))
                    .fold(f32::NEG_INFINITY, f32::max),
                points: members.len(),
                median_ppm_before: median_before,
                mad_ppm_before: mad_before,
                median_ppm_after: median_after,
                mad_ppm_after: mad_after,
            });
        }
    }
    out
}

/// Median and unscaled median absolute deviation.
pub fn median_mad(values: &[f32]) -> (f32, f32) {
    if values.is_empty() {
        return (f32::NAN, f32::NAN);
    }
    let center = median(values);
    let deviations = values
        .iter()
        .map(|v| (v - center).abs())
        .collect::<Vec<_>>();
    (center, median(&deviations))
}

fn median(values: &[f32]) -> f32 {
    if values.is_empty() {
        return f32::NAN;
    }
    let mut values = values.to_vec();
    let middle = values.len() / 2;
    let (_, upper, _) = values.select_nth_unstable_by(middle, f32::total_cmp);
    let upper = *upper;
    if values.len().is_multiple_of(2) {
        let lower = values[..middle]
            .iter()
            .copied()
            .fold(f32::NEG_INFINITY, f32::max);
        (lower + upper) / 2.0
    } else {
        upper
    }
}

#[cfg(test)]
#[path = "../tests/unit/mass_recalibration.rs"]
mod tests;
