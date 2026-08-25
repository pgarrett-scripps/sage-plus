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
    ) -> Option<LinearDiscriminantAnalysis> {
        Self::train_regularized(items, decoy, feat_fn, 0.0)
    }

    /// Fit LDA with diagonal regularization for correlated or sparse features.
    pub fn train_regularized<T, const D: usize>(
        items: &[T],
        decoy: &[bool],
        feat_fn: impl Fn(&T) -> [f64; D],
        regularization: f64,
    ) -> Option<LinearDiscriminantAnalysis> {
        assert_eq!(items.len(), decoy.len());

        // Pass 1: per-class sums -> per-class means.
        // Index 0 = decoy, 1 = target.
        let mut class_sum = [[0.0f64; D]; 2];
        let mut class_count = [0usize; 2];
        for (item, &is_decoy) in items.iter().zip(decoy) {
            let row = feat_fn(item);
            let cls = if is_decoy { 0 } else { 1 };
            for j in 0..D {
                class_sum[cls][j] += row[j];
            }
            class_count[cls] += 1;
        }
        if class_count[0] == 0 || class_count[1] == 0 {
            return None;
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
        let coef = Gauss::solve(scatter_within, Matrix::col_vector(mu_diff))?.take();

        log::trace!("- linear model fit with {:?}", Features(&coef));

        Some(LinearDiscriminantAnalysis { coef })
    }

    /// Project a single feature row.
    pub fn score(&self, row: &[f64]) -> f64 {
        debug_assert_eq!(row.len(), self.coef.len());
        self.coef.iter().zip(row).map(|(w, x)| w * x).sum()
    }
}

const LIBRARY_FEATURES: usize = 12;

fn library_feature_row(feature: &Feature) -> [f64; LIBRARY_FEATURES] {
    let finite = |value: f32| {
        if value.is_finite() {
            f64::from(value)
        } else {
            0.0
        }
    };
    let positive = |value: f32| finite(value.max(0.0));
    let rt_available = feature.predicted_rt.is_finite() && feature.predicted_rt > 0.0;
    let mobility_available = feature.ims.is_finite()
        && feature.ims > 0.0
        && feature.predicted_ims.is_finite()
        && feature.predicted_ims > 0.0;
    [
        positive(feature.spectral_angle),
        positive(feature.delta_next as f32),
        positive(feature.explained_library_intensity),
        positive(feature.explained_query_intensity),
        positive(feature.matched_peaks as f32).ln_1p(),
        finite(feature.aligned_delta_mass.abs()).ln_1p(),
        finite(feature.aligned_average_ppm.abs()).ln_1p(),
        finite(feature.isotope_error.abs()).ln_1p(),
        if rt_available {
            positive(feature.delta_rt_model)
        } else {
            0.0
        },
        if rt_available { 1.0 } else { 0.0 },
        if mobility_available {
            positive(feature.delta_ims_model)
        } else {
            0.0
        },
        if mobility_available { 1.0 } else { 0.0 },
    ]
}

/// Rescore library PSMs using library-specific target-decoy evidence.
pub fn score_library_psms(scores: &mut [Feature]) -> Option<()> {
    let decoys = scores
        .iter()
        .map(|score| score.label == -1)
        .collect::<Vec<_>>();
    let decoy_count = decoys.iter().filter(|decoy| **decoy).count();
    let target_count = decoys.len() - decoy_count;
    if target_count < 20 || decoy_count < 20 {
        return None;
    }

    let lda = LinearDiscriminantAnalysis::train_regularized::<_, LIBRARY_FEATURES>(
        scores,
        &decoys,
        library_feature_row,
        1e-3,
    )?;
    if !lda.coef.iter().all(|coefficient| coefficient.is_finite()) {
        return None;
    }
    log::trace!("library linear model coefficients: {:?}", lda.coef);
    let discriminants = scores
        .par_iter()
        .map(|score| lda.score(&library_feature_row(score)))
        .collect::<Vec<_>>();
    let kde = super::kde::Builder::default().build(&discriminants, &decoys);
    scores
        .par_iter_mut()
        .zip(discriminants)
        .for_each(|(feature, score)| {
            feature.discriminant_score = score as f32;
            feature.posterior_error = kde.posterior_error(score).log10() as f32;
            if feature.posterior_error.is_infinite() {
                feature.posterior_error = -324.0;
            }
        });
    Some(())
}

pub fn score_psms(scores: &mut [Feature], precursor_tol: Tolerance) -> Option<()> {
    log::trace!("fitting linear discriminant model...");
    let decoys = scores
        .par_iter()
        .map(|sc| sc.label == -1)
        .collect::<Vec<_>>();

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

    let lda = LinearDiscriminantAnalysis::train::<_, FEATURES>(scores, &decoys, &compute_features)?;
    if !lda.coef.iter().all(|f| f.is_finite()) {
        log::error!(
            "linear model coefficients include NaN: this likely indicates a bug, please report!"
        );
        for perc in scores.iter() {
            let row = compute_features(perc);
            if row.iter().any(|f| !f.is_finite()) {
                log::error!("example feature vector with NaN: {:?}", row);
                break;
            }
        }
        return None;
    }
    let discriminants: Vec<f64> = scores
        .par_iter()
        .map(|perc| lda.score(&compute_features(perc)))
        .collect();

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

    Some(())
}

#[cfg(test)]
#[path = "../../tests/unit/ml/linear_discriminant.rs"]
mod test;
