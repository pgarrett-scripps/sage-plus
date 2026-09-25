//! Linear Discriminant Analysis for FDR refinement
//!
//! "What I cannot create, I do not understand" - Richard Feynman
//!
//! One of the major reasons for the creation of Sage is to develop a search
//! engine from first principles - And when I mean first principles, I mean
//! first principles - we are going to implement a basic linear algebra system
//! (complete with Gauss-Jordan elimination and eigenvector calculation) from scratch
//! to enable LDA.

use super::gauss::Gauss;
use super::matrix::Matrix;
use rayon::prelude::*;

use crate::mass::Tolerance;
use crate::scoring::Feature;

// Declare, so that we have compile time checking of matrix dimensions
const FEATURES: usize = 20;
const FEATURE_NAMES: [&str; FEATURES] = [
    "rank",
    "charge",
    "ln1p(hyperscore)",
    "ln1p(delta_next)",
    "ln1p(delta_best)",
    "delta_mass_model",
    "isotope_error",
    "aligned_fragment_ppm",
    "ln1p(-poisson)",
    "ln1p(matched_intensity_pct)",
    "ln1p(matched_peaks)",
    "ln1p(longest_b)",
    "ln1p(longest_y)",
    "longest_y_pct",
    "ln1p(peptide_len)",
    "missed_cleavages",
    "rt",
    "ims",
    "sqrt(delta_rt_model)",
    "sqrt(delta_ims_model)",
];

struct Features<'a>(&'a [f64]);

impl std::fmt::Debug for Features<'_> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_map()
            .entries(FEATURE_NAMES.iter().zip(self.0))
            .finish()
    }
}

/// Fewest target and fewest decoy PSMs [`score_psms`] fits a model on: one
/// per feature. With fewer, the within-class scatter is rank deficient and
/// the fit is driven by the solver's ridge rather than the data.
pub const MIN_CLASS_PSMS: usize = FEATURES;

/// Largest accepted residual of the scatter solve, relative to the largest
/// class-mean difference. A larger residual means the scatter matrix was
/// singular in a direction the classes differ in, and the solver's ridge,
/// not the data, set the coefficients.
const MAX_RELATIVE_RESIDUAL: f64 = 1e-3;

/// A feature column whose within-class standard deviation is below this
/// fraction of `1 + |mean|` counts as constant: it only differs by rounding.
const CONSTANT_TOLERANCE: f64 = 1e-9;

/// Why a linear discriminant model could not be used.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum LdaFailure {
    /// Too few target or decoy PSMs to estimate class means and scatter.
    TooFewPsms { targets: usize, decoys: usize },
    /// A feature row, the scatter matrix, or the fitted coefficients held
    /// NaN or infinite values.
    NonFinite,
    /// The within-class scatter matrix could not be solved.
    SingularScatter,
    /// The model does not separate the classes: identical class means, or
    /// the same score for every PSM.
    Degenerate,
}

impl std::fmt::Display for LdaFailure {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            LdaFailure::TooFewPsms { targets, decoys } => write!(
                f,
                "too few PSMs to fit ({targets} targets, {decoys} decoys; need {MIN_CLASS_PSMS} of each)"
            ),
            LdaFailure::NonFinite => write!(f, "features or coefficients are not finite"),
            LdaFailure::SingularScatter => write!(f, "within-class scatter matrix is singular"),
            LdaFailure::Degenerate => write!(f, "model does not separate targets from decoys"),
        }
    }
}

impl std::error::Error for LdaFailure {}

pub struct LinearDiscriminantAnalysis {
    coef: Vec<f64>,
}

impl LinearDiscriminantAnalysis {
    /// Fit LDA over `items`, where `feat_fn(item)` produces a feature row.
    /// `decoy[i] == true` means item i is a decoy.
    ///
    /// Two streaming passes -- class means then within-class scatter -- so
    /// the `n x D` feature matrix is never materialized.
    pub fn train<T, const D: usize>(
        items: &[T],
        decoy: &[bool],
        feat_fn: impl Fn(&T) -> [f64; D],
    ) -> Result<LinearDiscriminantAnalysis, LdaFailure> {
        Self::train_regularized(items, decoy, feat_fn, 0.0)
    }

    /// Fit LDA with diagonal regularization for correlated or sparse features.
    pub fn train_regularized<T, const D: usize>(
        items: &[T],
        decoy: &[bool],
        feat_fn: impl Fn(&T) -> [f64; D],
        regularization: f64,
    ) -> Result<LinearDiscriminantAnalysis, LdaFailure> {
        assert_eq!(items.len(), decoy.len());

        // Pass 1: per-class sums -> per-class means.
        // Index 0 = decoy, 1 = target.
        let mut class_sum = [[0.0f64; D]; 2];
        let mut class_count = [0usize; 2];
        for (item, &is_decoy) in items.iter().zip(decoy) {
            let row = feat_fn(item);
            if row.iter().any(|x| !x.is_finite()) {
                return Err(LdaFailure::NonFinite);
            }
            let cls = if is_decoy { 0 } else { 1 };
            for j in 0..D {
                class_sum[cls][j] += row[j];
            }
            class_count[cls] += 1;
        }
        if class_count[0] == 0 || class_count[1] == 0 {
            return Err(LdaFailure::TooFewPsms {
                targets: class_count[1],
                decoys: class_count[0],
            });
        }
        let mut class_mean = [[0.0f64; D]; 2];
        for c in 0..2 {
            let nf = class_count[c] as f64;
            for j in 0..D {
                class_mean[c][j] = class_sum[c][j] / nf;
            }
        }

        // Pass 2: per-class within-class scatter sum_i (x_i - mu_c)(x_i - mu_c)^T.
        let mut scatter_per_class: [Matrix; 2] = [Matrix::zeros(D, D), Matrix::zeros(D, D)];
        for (item, &is_decoy) in items.iter().zip(decoy) {
            let row = feat_fn(item);
            let cls = if is_decoy { 0 } else { 1 };
            let mu = &class_mean[cls];
            let centered: [f64; D] = std::array::from_fn(|j| row[j] - mu[j]);
            for j in 0..D {
                for k in 0..D {
                    scatter_per_class[cls][(j, k)] += centered[j] * centered[k];
                }
            }
        }
        // scatter_within = sum_c (X_c - mu_c)^T (X_c - mu_c) / n_c
        let mut scatter_within = Matrix::zeros(D, D);
        for c in 0..2 {
            scatter_within += scatter_per_class[c].clone() / class_count[c] as f64;
        }
        let within_variance: [f64; D] = std::array::from_fn(|j| scatter_within[(j, j)]);
        if regularization > 0.0 {
            let scale = (0..D)
                .map(|index| scatter_within[(index, index)])
                .sum::<f64>()
                / D as f64;
            let ridge = scale.max(1.0) * regularization;
            for index in 0..D {
                scatter_within[(index, index)] += ridge;
            }
        }

        // For two-class LDA, Sb is rank-1 in the direction of (mu_t - mu_d), so
        // the dominant eigenvector of Sw^-1 Sb is parallel to Sw^-1 (mu_t - mu_d).
        // Solve Sw * w = (mu_t - mu_d) directly. Target projects higher than
        // decoy by construction since Sw^-1 is positive-definite.
        let mu_diff: Vec<f64> = (0..D)
            .map(|j| class_mean[1][j] - class_mean[0][j])
            .collect();
        if mu_diff.iter().any(|x| !x.is_finite())
            || (0..D).any(|j| !scatter_within[(j, j)].is_finite())
        {
            return Err(LdaFailure::NonFinite);
        }

        // A column that is constant across all PSMs (ion mobility on Orbitrap
        // data, rank when only rank 1 is reported, model deltas without a
        // model) has no within-class variance and no class-mean difference.
        // It carries no information but makes the scatter matrix singular, so
        // it gets zero weight and the model is solved over the other columns.
        // A column constant within each class but differing between them
        // separates the classes perfectly; no finite model exists.
        let total = (class_count[0] + class_count[1]) as f64;
        let mut active = Vec::with_capacity(D);
        for j in 0..D {
            let mean = (class_sum[0][j] + class_sum[1][j]) / total;
            let scale = 1.0 + mean.abs();
            if within_variance[j] > (CONSTANT_TOLERANCE * scale).powi(2) {
                active.push(j);
            } else if mu_diff[j].abs() > CONSTANT_TOLERANCE * scale {
                return Err(LdaFailure::SingularScatter);
            }
        }
        if active.len() < D {
            log::debug!(
                "linear model ignores {} constant feature column(s): {:?}",
                D - active.len(),
                (0..D).filter(|j| !active.contains(j)).collect::<Vec<_>>()
            );
        }

        let n = active.len();
        let mut scatter = Matrix::zeros(n, n);
        for (a, &j) in active.iter().enumerate() {
            for (b, &k) in active.iter().enumerate() {
                scatter[(a, b)] = scatter_within[(j, k)];
            }
        }
        let target: Vec<f64> = active.iter().map(|&j| mu_diff[j]).collect();
        let max_diff = target.iter().fold(0.0f64, |acc, x| acc.max(x.abs()));
        if max_diff == 0.0 {
            return Err(LdaFailure::Degenerate);
        }
        let solved = Gauss::solve(scatter.clone(), Matrix::col_vector(target.clone()))
            .ok_or(LdaFailure::SingularScatter)?
            .take();
        if solved.iter().any(|w| !w.is_finite()) {
            return Err(LdaFailure::NonFinite);
        }
        // The solver adds a growing ridge until elimination succeeds, so it
        // "solves" singular systems too. Accept the fit only if it solves the
        // unridged system.
        let residual = scatter
            .dotv(&solved)
            .iter()
            .zip(&target)
            .fold(0.0f64, |acc, (fitted, target)| {
                acc.max((fitted - target).abs())
            });
        let mut coef = vec![0.0; D];
        for (&j, &w) in active.iter().zip(&solved) {
            coef[j] = w;
        }
        if residual.is_nan() || residual > MAX_RELATIVE_RESIDUAL * max_diff {
            log::debug!(
                "linear model residual {:e} exceeds {:e} of the class-mean difference {:e}",
                residual,
                MAX_RELATIVE_RESIDUAL,
                max_diff
            );
            return Err(LdaFailure::SingularScatter);
        }

        log::trace!("- linear model fit with {:?}", Features(&coef));

        Ok(LinearDiscriminantAnalysis { coef })
    }

    /// Project a single feature row.
    pub fn score(&self, row: &[f64]) -> f64 {
        debug_assert_eq!(row.len(), self.coef.len());
        self.coef.iter().zip(row).map(|(w, x)| w * x).sum()
    }
}

/// Fit the linear discriminant on `scores` and set each PSM's
/// `discriminant_score` and `posterior_error`.
///
/// On failure `scores` are left untouched and the reason is returned; callers
/// then rank PSMs with [`score_psms_fallback`].
pub fn score_psms(scores: &mut [Feature], precursor_tol: Tolerance) -> Result<(), LdaFailure> {
    log::trace!("fitting linear discriminant model...");
    let decoys = scores
        .par_iter()
        .map(|sc| sc.label == -1)
        .collect::<Vec<_>>();
    let decoy_count = decoys.iter().filter(|&&decoy| decoy).count();
    let target_count = decoys.len() - decoy_count;
    if target_count < MIN_CLASS_PSMS || decoy_count < MIN_CLASS_PSMS {
        return Err(LdaFailure::TooFewPsms {
            targets: target_count,
            decoys: decoy_count,
        });
    }

    let mass_error = match precursor_tol {
        Tolerance::Ppm(_, _) => |feat: &Feature| feat.aligned_delta_mass as f64,
        Tolerance::Pct(_, _) => unreachable!("Pct tolerance should never be used on mz"),
        Tolerance::Da(_, _) => |feat: &Feature| (feat.expmass - feat.calcmass) as f64,
    };

    let (bw_adjust, bin_size) = match precursor_tol {
        Tolerance::Ppm(lo, hi) => (2.0f64, (hi - lo).max(100.0)),
        Tolerance::Pct(_, _) => unreachable!("Pct tolerance should never be used on mz"),
        Tolerance::Da(lo, hi) => (0.1f64, (hi - lo).max(1000.0)),
    };

    let delta_mass = scores.par_iter().map(mass_error).collect::<Vec<_>>();

    let mass_model = super::kde::Builder::default()
        .monotonic(false)
        .bw_adjust(move |x| x * bw_adjust)
        .bins(bin_size.ceil().abs() as usize)
        .build(&delta_mass, &decoys);

    // Compute a feature row on demand. Called twice in `train` (means + scatter)
    // and once during scoring -- no `n_psms x FEATURES` matrix is materialized.
    let compute_features = |perc: &Feature| -> [f64; FEATURES] {
        let poisson = match (-perc.poisson).ln_1p() {
            x if x.is_finite() => x,
            _ => 3.5,
        };

        // Transform features - LDA requires that each feature is normally
        // distributed. This is not true for all of our inputs, so we log
        // transform many of them to get them closer to a gaussian distr.
        [
            (perc.rank as f64),
            (perc.charge as f64),
            (perc.hyperscore).ln_1p(),
            (perc.delta_next).ln_1p(),
            (perc.delta_best).ln_1p(),
            mass_model.posterior_error(mass_error(perc)),
            (perc.isotope_error as f64),
            (perc.aligned_average_ppm as f64),
            (poisson),
            (perc.matched_intensity_pct as f64).ln_1p(),
            (perc.matched_peaks as f64),
            (perc.longest_b as f64).ln_1p(),
            (perc.longest_y as f64).ln_1p(),
            (perc.longest_y as f64 / perc.peptide_len as f64),
            (perc.peptide_len as f64).ln_1p(),
            (perc.missed_cleavages as f64),
            (perc.aligned_rt as f64),
            (perc.ims as f64),
            (perc.delta_rt_model as f64).clamp(0.001, 0.999).sqrt(),
            (perc.delta_ims_model as f64).clamp(0.001, 0.999).sqrt(),
        ]
    };

    let lda = LinearDiscriminantAnalysis::train::<_, FEATURES>(scores, &decoys, &compute_features)
        .inspect_err(|failure| {
            if *failure == LdaFailure::NonFinite {
                if let Some(row) = scores
                    .iter()
                    .map(&compute_features)
                    .find(|row| row.iter().any(|f| !f.is_finite()))
                {
                    log::warn!("example feature vector with NaN: {:?}", Features(&row));
                }
            }
        })?;
    let discriminants: Vec<f64> = scores
        .par_iter()
        .map(|perc| lda.score(&compute_features(perc)))
        .collect();
    if discriminants.iter().any(|score| !score.is_finite()) {
        return Err(LdaFailure::NonFinite);
    }
    let (lo, hi) = discriminants
        .iter()
        .fold((f64::INFINITY, f64::NEG_INFINITY), |(lo, hi), &x| {
            (lo.min(x), hi.max(x))
        });
    if lo >= hi {
        return Err(LdaFailure::Degenerate);
    }

    log::trace!("- fitting non-parametric model for posterior error probabilities");
    let kde = super::kde::Builder::default().build(&discriminants, &decoys);

    scores
        .par_iter_mut()
        .zip(&discriminants)
        .for_each(|(perc, score)| {
            perc.discriminant_score = *score as f32;
            perc.posterior_error = kde.posterior_error(*score).log10() as f32;
            if perc.posterior_error.is_infinite() {
                // This is approximately the log10 of the smallest positive
                // non-zero f64
                perc.posterior_error = -324.0;
            }
        });

    Ok(())
}

/// Heuristic discriminant used when [`score_psms`] cannot fit a model:
/// `ln(1 - poisson) + longest_y_pct / 3`. A non-finite Poisson term (an
/// underflowed match probability) counts as [`FALLBACK_POISSON_CAP`] so the
/// strongest matches rank first, and NaN counts as zero.
pub fn fallback_discriminant(feature: &Feature) -> f32 {
    let poisson = (-feature.poisson as f32).ln_1p();
    let poisson = if poisson.is_nan() {
        0.0
    } else {
        poisson.min(FALLBACK_POISSON_CAP)
    };
    let longest_y = if feature.longest_y_pct.is_finite() {
        feature.longest_y_pct
    } else {
        0.0
    };
    poisson + longest_y / 3.0
}

/// Above `ln(1 + x)` for any finite f64 Poisson log-probability (at most
/// about 745 in magnitude).
pub const FALLBACK_POISSON_CAP: f32 = 8.0;

/// Rank PSMs with [`fallback_discriminant`] and estimate posterior error
/// probabilities from it. Where no estimate is possible (no decoys or no
/// targets, or a non-finite density) `posterior_error` is 0, a PEP of 1.
pub fn score_psms_fallback(scores: &mut [Feature]) {
    scores.par_iter_mut().for_each(|feat| {
        feat.discriminant_score = fallback_discriminant(feat);
        feat.posterior_error = 0.0;
    });
    let decoys = scores
        .iter()
        .map(|feat| feat.label == -1)
        .collect::<Vec<_>>();
    if decoys.iter().all(|&decoy| decoy) || !decoys.iter().any(|&decoy| decoy) {
        return;
    }
    let discriminants = scores
        .iter()
        .map(|feat| feat.discriminant_score as f64)
        .collect::<Vec<_>>();
    let kde = super::kde::Builder::default().build(&discriminants, &decoys);
    scores
        .par_iter_mut()
        .zip(&discriminants)
        .for_each(|(feat, score)| {
            let pep = kde.posterior_error(*score).log10() as f32;
            feat.posterior_error = if pep.is_nan() {
                0.0
            } else if pep.is_infinite() {
                -324.0
            } else {
                pep
            };
        });
}

#[cfg(test)]
#[path = "../../tests/unit/ml/linear_discriminant.rs"]
mod test;
