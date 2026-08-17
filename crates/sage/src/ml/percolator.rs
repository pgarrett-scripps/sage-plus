//! Percolator-style semi-supervised PSM rescoring with out-of-fold predictions.

use crate::scoring::Feature;
use crate::{mass::Tolerance, ml::kde::Estimator};
use fnv::FnvHasher;
use perpetual::{objective::Objective, Matrix, PerpetualBooster};
use serde::{Deserialize, Serialize};
use std::collections::HashMap;
use std::fmt::{Display, Formatter};
use std::hash::{Hash, Hasher};

const FEATURES: usize = 20;

#[derive(Clone, Copy, Debug, Default, Serialize, Deserialize, PartialEq, Eq)]
#[serde(rename_all = "snake_case")]
pub enum PercolatorModel {
    #[default]
    Svm,
    Perpetual,
}

#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(default)]
pub struct PercolatorSettings {
    /// Enable Percolator-style rescoring. LDA remains the default.
    pub enabled: bool,
    /// Rescoring model: linear SVM or Perpetual boosted trees.
    pub model: PercolatorModel,
    /// Perpetual model complexity budget.
    pub budget: f32,
    /// Linear SVM hinge-loss penalty.
    pub svm_c: f64,
    /// Number of batch-gradient epochs for each SVM fit.
    pub svm_epochs: usize,
    /// Number of cross-validation folds.
    pub folds: usize,
    /// Maximum number of semi-supervised positive-set updates.
    pub iterations: usize,
    /// Maximum q-value for provisional positive training examples.
    pub train_fdr: f32,
    /// Minimum provisional positives required in every training split.
    pub min_positive_psms: usize,
    /// Minimum decoys required in every training split.
    pub min_decoy_psms: usize,
    /// Seed used for deterministic spectrum-level fold assignment.
    pub seed: u64,
}

impl Default for PercolatorSettings {
    fn default() -> Self {
        Self {
            enabled: false,
            model: PercolatorModel::Svm,
            budget: 0.3,
            svm_c: 1.0,
            svm_epochs: 100,
            folds: 3,
            iterations: 3,
            train_fdr: 0.01,
            min_positive_psms: 50,
            min_decoy_psms: 50,
            seed: 42,
        }
    }
}

#[derive(Debug)]
pub struct RescoreError(String);

impl Display for RescoreError {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        f.write_str(&self.0)
    }
}

impl std::error::Error for RescoreError {}

fn err(message: impl Into<String>) -> RescoreError {
    RescoreError(message.into())
}

fn validate(settings: &PercolatorSettings) -> Result<(), RescoreError> {
    if !(2..=10).contains(&settings.folds) {
        return Err(err("percolator.folds must be between 2 and 10"));
    }
    if settings.iterations == 0 {
        return Err(err("percolator.iterations must be at least 1"));
    }
    if settings.model == PercolatorModel::Perpetual
        && (!settings.budget.is_finite() || settings.budget <= 0.0)
    {
        return Err(err("percolator.budget must be finite and greater than 0"));
    }
    if settings.model == PercolatorModel::Svm
        && (!settings.svm_c.is_finite() || settings.svm_c <= 0.0)
    {
        return Err(err("percolator.svm_c must be finite and greater than 0"));
    }
    if settings.model == PercolatorModel::Svm && settings.svm_epochs == 0 {
        return Err(err("percolator.svm_epochs must be at least 1"));
    }
    if !settings.train_fdr.is_finite() || settings.train_fdr <= 0.0 || settings.train_fdr >= 0.5 {
        return Err(err("percolator.train_fdr must be between 0 and 0.5"));
    }
    Ok(())
}

/// Label-independent features shared by both rescorers. `mass_error_pep` is
/// fitted using only the current fold's training PSMs.
fn feature_row(psm: &Feature, mass_error_pep: f64) -> [f64; FEATURES] {
    let poisson = match (-psm.poisson).ln_1p() {
        value if value.is_finite() => value,
        _ => 3.5,
    };
    let delta_mass = psm.aligned_delta_mass as f64;
    let fragment_ppm = psm.aligned_average_ppm as f64;

    [
        psm.charge as f64,
        (psm.hyperscore).ln_1p(),
        (psm.delta_next).ln_1p(),
        (psm.delta_best).ln_1p(),
        delta_mass,
        delta_mass.abs(),
        mass_error_pep,
        psm.isotope_error as f64,
        fragment_ppm,
        fragment_ppm.abs(),
        poisson,
        (psm.matched_intensity_pct as f64).ln_1p(),
        psm.matched_peaks as f64,
        (psm.longest_b as f64).ln_1p(),
        (psm.longest_y as f64).ln_1p(),
        psm.longest_y_pct as f64,
        (psm.peptide_len as f64).ln_1p(),
        psm.missed_cleavages as f64,
        (psm.delta_rt_model as f64).clamp(0.001, 0.999).sqrt(),
        (psm.delta_ims_model as f64).clamp(0.001, 0.999).sqrt(),
    ]
}

fn fold(psm: &Feature, settings: &PercolatorSettings) -> usize {
    let mut hasher = FnvHasher::default();
    settings.seed.hash(&mut hasher);
    psm.file_id.hash(&mut hasher);
    psm.spec_id.hash(&mut hasher);
    hasher.finish() as usize % settings.folds
}

/// Retain the best rescored candidate per spectrum for target-decoy confidence
/// estimation. All candidates remain available for model training and output.
fn competition(features: &[Feature], scores: &[f64], indices: &[usize]) -> Vec<usize> {
    let mut winners: HashMap<(usize, &str), usize> = HashMap::new();
    for &idx in indices {
        let key = (features[idx].file_id, features[idx].spec_id.as_str());
        winners
            .entry(key)
            .and_modify(|winner| {
                if scores[idx].total_cmp(&scores[*winner]).is_gt() {
                    *winner = idx;
                }
            })
            .or_insert(idx);
    }
    winners.into_values().collect()
}

/// Compute target-decoy q-values for the supplied, already competed indices.
fn q_values(features: &[Feature], scores: &[f64], indices: &[usize]) -> Vec<(usize, f32)> {
    let mut order = indices.to_vec();
    order.sort_unstable_by(|&a, &b| scores[b].total_cmp(&scores[a]));

    let mut decoys = 1usize;
    let mut targets = 0usize;
    let mut values = Vec::with_capacity(order.len());
    for &idx in &order {
        if features[idx].label == -1 {
            decoys += 1;
        } else {
            targets += 1;
        }
        values.push((idx, decoys as f32 / targets as f32));
    }

    let mut minimum = 1.0f32;
    for (_, q) in values.iter_mut().rev() {
        minimum = minimum.min(*q);
        *q = minimum;
    }
    values
}

fn positives(
    features: &[Feature],
    scores: &[f64],
    indices: &[usize],
    threshold: f32,
) -> Vec<usize> {
    let winners = competition(features, scores, indices);
    q_values(features, scores, &winners)
        .into_iter()
        .filter_map(|(idx, q)| (features[idx].label != -1 && q <= threshold).then_some(idx))
        .collect()
}

fn mass_error(psm: &Feature, precursor_tol: Tolerance) -> f64 {
    match precursor_tol {
        Tolerance::Ppm(_, _) => psm.aligned_delta_mass as f64,
        Tolerance::Pct(_, _) => unreachable!("Pct tolerance should never be used on mz"),
        Tolerance::Da(_, _) => (psm.expmass - psm.calcmass) as f64,
    }
}

fn fold_rows(
    features: &[Feature],
    training: &[usize],
    precursor_tol: Tolerance,
) -> Result<Vec<[f64; FEATURES]>, RescoreError> {
    let errors = training
        .iter()
        .map(|&idx| mass_error(&features[idx], precursor_tol))
        .collect::<Vec<_>>();
    let decoys = training
        .iter()
        .map(|&idx| features[idx].label == -1)
        .collect::<Vec<_>>();
    let (bw_adjust, bins) = match precursor_tol {
        Tolerance::Ppm(lo, hi) => (2.0, (hi - lo).max(100.0).ceil().abs() as usize),
        Tolerance::Pct(_, _) => unreachable!("Pct tolerance should never be used on mz"),
        Tolerance::Da(lo, hi) => (0.1, (hi - lo).max(1000.0).ceil().abs() as usize),
    };
    let model: Estimator = super::kde::Builder::default()
        .monotonic(false)
        .bw_adjust(move |value| value * bw_adjust)
        .bins(bins)
        .build(&errors, &decoys);
    let rows = features
        .iter()
        .map(|psm| {
            let pep = model.posterior_error(mass_error(psm, precursor_tol));
            feature_row(psm, pep)
        })
        .collect::<Vec<_>>();
    if rows.iter().flatten().any(|value| !value.is_finite()) {
        return Err(err(
            "fold-specific mass-error KDE produced a non-finite feature",
        ));
    }
    Ok(rows)
}

/// Perpetual expects a contiguous column-major matrix.
fn matrix(rows: &[[f64; FEATURES]], indices: &[usize]) -> Vec<f64> {
    let mut columns = Vec::with_capacity(indices.len() * FEATURES);
    for column in 0..FEATURES {
        columns.extend(indices.iter().map(|&idx| rows[idx][column]));
    }
    columns
}

#[derive(Clone)]
struct LinearSvm {
    weights: [f64; FEATURES],
    bias: f64,
    means: [f64; FEATURES],
    scales: [f64; FEATURES],
}

impl LinearSvm {
    fn fit(
        rows: &[[f64; FEATURES]],
        features: &[Feature],
        selected: &[usize],
        positive_count: usize,
        decoy_count: usize,
        settings: &PercolatorSettings,
    ) -> Self {
        let mut means = [0.0; FEATURES];
        for &idx in selected {
            for (mean, value) in means.iter_mut().zip(rows[idx]) {
                *mean += value;
            }
        }
        for mean in &mut means {
            *mean /= selected.len() as f64;
        }
        let mut scales = [0.0; FEATURES];
        for &idx in selected {
            for column in 0..FEATURES {
                scales[column] += (rows[idx][column] - means[column]).powi(2);
            }
        }
        for scale in &mut scales {
            *scale = (*scale / selected.len() as f64).sqrt().max(1e-8);
        }

        let positive_weight = selected.len() as f64 / (2.0 * positive_count as f64);
        let decoy_weight = selected.len() as f64 / (2.0 * decoy_count as f64);
        let mut weights = [0.0; FEATURES];
        let mut bias = 0.0;
        for epoch in 0..settings.svm_epochs {
            let mut gradient = weights;
            let mut bias_gradient = 0.0;
            for &idx in selected {
                let label = if features[idx].label == -1 { -1.0 } else { 1.0 };
                let sample_weight = if label > 0.0 {
                    positive_weight
                } else {
                    decoy_weight
                };
                let margin = bias
                    + (0..FEATURES)
                        .map(|column| {
                            weights[column] * (rows[idx][column] - means[column]) / scales[column]
                        })
                        .sum::<f64>();
                if label * margin < 1.0 {
                    let factor = settings.svm_c * sample_weight * label / selected.len() as f64;
                    for column in 0..FEATURES {
                        gradient[column] -=
                            factor * (rows[idx][column] - means[column]) / scales[column];
                    }
                    bias_gradient -= factor;
                }
            }
            let learning_rate = 0.2 / (1.0 + epoch as f64 * 0.05);
            for column in 0..FEATURES {
                weights[column] -= learning_rate * gradient[column];
            }
            bias -= learning_rate * bias_gradient;
        }
        Self {
            weights,
            bias,
            means,
            scales,
        }
    }

    fn predict(&self, row: &[f64; FEATURES]) -> f64 {
        self.bias
            + (0..FEATURES)
                .map(|column| {
                    self.weights[column] * (row[column] - self.means[column]) / self.scales[column]
                })
                .sum::<f64>()
    }
}

enum TrainedModel {
    Svm(LinearSvm),
    Perpetual(PerpetualBooster),
}

impl TrainedModel {
    fn predict(&self, rows: &[[f64; FEATURES]], indices: &[usize]) -> Vec<f64> {
        match self {
            Self::Svm(model) => indices
                .iter()
                .map(|&idx| model.predict(&rows[idx]))
                .collect(),
            Self::Perpetual(model) => {
                let data = matrix(rows, indices);
                let matrix = Matrix::new(&data, indices.len(), FEATURES);
                model.predict(&matrix, true)
            }
        }
    }
}

fn train(
    features: &[Feature],
    rows: &[[f64; FEATURES]],
    scores: &mut [f64],
    training: &[usize],
    settings: &PercolatorSettings,
    model_seed: u64,
) -> Result<TrainedModel, RescoreError> {
    let decoys = training
        .iter()
        .copied()
        .filter(|&idx| features[idx].label == -1)
        .collect::<Vec<_>>();
    if decoys.len() < settings.min_decoy_psms {
        return Err(err(format!(
            "only {} decoys available; {} required",
            decoys.len(),
            settings.min_decoy_psms
        )));
    }

    let mut last_positives = Vec::new();
    let mut final_model = None;
    for _ in 0..settings.iterations {
        let positive = positives(features, scores, training, settings.train_fdr);
        if positive.len() < settings.min_positive_psms {
            return Err(err(format!(
                "only {} positive PSMs at {:.2}% FDR; {} required",
                positive.len(),
                settings.train_fdr * 100.0,
                settings.min_positive_psms
            )));
        }

        let mut selected = positive.clone();
        selected.extend_from_slice(&decoys);
        let model = match settings.model {
            PercolatorModel::Svm => TrainedModel::Svm(LinearSvm::fit(
                rows,
                features,
                &selected,
                positive.len(),
                decoys.len(),
                settings,
            )),
            PercolatorModel::Perpetual => {
                let data = matrix(rows, &selected);
                let selected_matrix = Matrix::new(&data, selected.len(), FEATURES);
                let labels = selected
                    .iter()
                    .map(|&idx| if features[idx].label == -1 { 0.0 } else { 1.0 })
                    .collect::<Vec<_>>();
                let positive_weight = selected.len() as f64 / (2.0 * positive.len() as f64);
                let decoy_weight = selected.len() as f64 / (2.0 * decoys.len() as f64);
                let weights = labels
                    .iter()
                    .map(|&label| {
                        if label == 1.0 {
                            positive_weight
                        } else {
                            decoy_weight
                        }
                    })
                    .collect::<Vec<_>>();
                let mut model = PerpetualBooster::default()
                    .set_objective(Objective::LogLoss)
                    .set_budget(settings.budget)
                    .set_seed(model_seed);
                model
                    .fit(&selected_matrix, &labels, Some(&weights), None)
                    .map_err(|error| err(format!("Perpetual training failed: {error}")))?;
                TrainedModel::Perpetual(model)
            }
        };

        for (&idx, prediction) in training.iter().zip(model.predict(rows, training)) {
            scores[idx] = prediction;
        }
        final_model = Some(model);

        let mut sorted_positive = positive;
        sorted_positive.sort_unstable();
        if sorted_positive == last_positives {
            break;
        }
        last_positives = sorted_positive;
    }

    final_model.ok_or_else(|| err("Perpetual did not produce a model"))
}

fn normalize_folds(
    features: &[Feature],
    scores: &mut [f64],
    assignments: &[usize],
    settings: &PercolatorSettings,
) -> Result<(), RescoreError> {
    for held_out in 0..settings.folds {
        let indices = assignments
            .iter()
            .enumerate()
            .filter_map(|(idx, &fold)| (fold == held_out).then_some(idx))
            .collect::<Vec<_>>();
        let winners = competition(features, scores, &indices);
        let q = q_values(features, scores, &winners);
        let cutoff = q
            .iter()
            .filter_map(|&(idx, q)| {
                (features[idx].label != -1 && q <= settings.train_fdr).then_some(scores[idx])
            })
            .min_by(f64::total_cmp)
            .ok_or_else(|| {
                err(format!(
                    "fold {held_out} has no target PSMs at {:.2}% FDR",
                    settings.train_fdr * 100.0
                ))
            })?;

        let mut decoys = winners
            .iter()
            .filter_map(|&idx| (features[idx].label == -1).then_some(scores[idx]))
            .collect::<Vec<_>>();
        decoys.sort_unstable_by(f64::total_cmp);
        let median_decoy = decoys
            .get(decoys.len() / 2)
            .copied()
            .ok_or_else(|| err(format!("fold {held_out} has no held-out decoys")))?;
        let scale = cutoff - median_decoy;
        if !scale.is_finite() || scale <= f64::EPSILON {
            return Err(err(format!("fold {held_out} has an invalid score scale")));
        }
        for &idx in &indices {
            scores[idx] = (scores[idx] - cutoff) / scale;
        }
    }
    Ok(())
}

/// Score every PSM with a model that did not train on that spectrum.
///
/// Fold scores are aligned using Percolator's convention: the score at the
/// training FDR is mapped to 0 and the median held-out decoy is mapped to -1.
pub fn score_psms(
    features: &mut [Feature],
    settings: &PercolatorSettings,
    precursor_tol: Tolerance,
) -> Result<(), RescoreError> {
    validate(settings)?;
    if features.is_empty() {
        return Err(err("no PSMs available for Percolator-style rescoring"));
    }

    let assignments = features
        .iter()
        .map(|psm| fold(psm, settings))
        .collect::<Vec<_>>();
    // Poisson probability is the stable, label-independent initial direction.
    let initial_scores = features.iter().map(|psm| -psm.poisson).collect::<Vec<_>>();
    let mut scores = vec![f64::NAN; features.len()];

    for held_out in 0..settings.folds {
        let training = assignments
            .iter()
            .enumerate()
            .filter_map(|(idx, &fold)| (fold != held_out).then_some(idx))
            .collect::<Vec<_>>();
        let testing = assignments
            .iter()
            .enumerate()
            .filter_map(|(idx, &fold)| (fold == held_out).then_some(idx))
            .collect::<Vec<_>>();
        if testing.is_empty() {
            return Err(err(format!("fold {held_out} is empty")));
        }
        let rows = fold_rows(features, &training, precursor_tol)?;

        // Each fold's semi-supervised iterations are independent. Only the
        // held-out predictions are copied into the final score vector.
        let mut training_scores = initial_scores.clone();
        let model = train(
            features,
            &rows,
            &mut training_scores,
            &training,
            settings,
            settings.seed.wrapping_add(held_out as u64),
        )?;
        for (&idx, prediction) in testing.iter().zip(model.predict(&rows, &testing)) {
            scores[idx] = prediction;
        }
    }

    normalize_folds(features, &mut scores, &assignments, settings)?;
    if scores.iter().any(|score| !score.is_finite()) {
        return Err(err("Perpetual produced a non-finite score"));
    }
    if settings.model == PercolatorModel::Perpetual {
        // Tree ensembles produce discrete leaf scores and fold normalization
        // can create extreme tails. Compress them monotonically for stable KDE
        // PEPs, then use Poisson only to break exact leaf-score ties.
        for (score, initial) in scores.iter_mut().zip(&initial_scores) {
            *score = score.asinh() + 1e-6 * initial.atan();
        }
    }
    let decoys = features
        .iter()
        .map(|psm| psm.label == -1)
        .collect::<Vec<_>>();
    let pep = super::kde::Builder::default().build(&scores, &decoys);
    for (feature, score) in features.iter_mut().zip(scores) {
        feature.discriminant_score = score as f32;
        feature.posterior_error = pep.posterior_error(score).log10() as f32;
        if feature.posterior_error.is_infinite() {
            feature.posterior_error = -324.0;
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn q_values_include_plus_one_correction() {
        let mut features = vec![Feature::default(); 4];
        features[0].label = 1;
        features[1].label = 1;
        features[2].label = -1;
        features[3].label = 1;
        let scores = vec![4.0, 3.0, 2.0, 1.0];
        let q = q_values(&features, &scores, &[0, 1, 2, 3]);
        assert_eq!(q[0], (0, 0.5));
        assert_eq!(q[1], (1, 0.5));
        assert!((q[2].1 - 2.0 / 3.0).abs() < 1e-6);
    }

    #[test]
    fn spectrum_candidates_share_a_fold() {
        let settings = PercolatorSettings::default();
        let mut first = Feature::default();
        first.file_id = 7;
        first.spec_id = "scan=10".into();
        let mut second = first.clone();
        second.peptide_idx.0 = second.peptide_idx.0.wrapping_add(1);
        assert_eq!(fold(&first, &settings), fold(&second, &settings));
    }

    #[test]
    fn feature_rows_use_aligned_mass_errors() {
        let feature = Feature {
            delta_mass: 100.0,
            average_ppm: 200.0,
            aligned_delta_mass: 1.5,
            aligned_average_ppm: 2.5,
            ..Feature::default()
        };
        let row = feature_row(&feature, 0.0);
        assert_eq!(row[4], 1.5);
        assert_eq!(row[8], 2.5);
    }

    #[test]
    fn held_out_labels_do_not_change_mass_error_features() {
        let mut features = Vec::new();
        for idx in 0..120 {
            let mut feature = Feature::default();
            feature.label = if idx % 3 == 0 { -1 } else { 1 };
            feature.delta_mass = if feature.label == -1 {
                4.0 + (idx % 7) as f32 * 0.2
            } else {
                (idx % 7) as f32 * 0.02
            };
            feature.aligned_delta_mass = feature.delta_mass;
            features.push(feature);
        }
        let training = (0..100).collect::<Vec<_>>();
        let before = fold_rows(&features, &training, Tolerance::Ppm(-10.0, 10.0)).unwrap();
        features[110].label *= -1;
        let after = fold_rows(&features, &training, Tolerance::Ppm(-10.0, 10.0)).unwrap();
        assert_eq!(before[110], after[110]);
    }

    fn assert_learns_separable_psms(model: PercolatorModel) {
        let mut features = Vec::new();
        for idx in 0..900 {
            let target = idx % 3 != 0;
            let mut feature = Feature::default();
            feature.spec_id = format!("scan={idx}");
            feature.label = if target { 1 } else { -1 };
            feature.poisson = if target { -20.0 } else { -0.1 };
            feature.hyperscore = if target { 100.0 } else { 1.0 };
            feature.matched_peaks = if target { 20 } else { 2 };
            feature.delta_mass = if target {
                0.1 + (idx % 5) as f32 * 0.01
            } else {
                5.0 + (idx % 5) as f32 * 0.1
            };
            feature.aligned_delta_mass = feature.delta_mass;
            features.push(feature);
        }
        let settings = PercolatorSettings {
            model,
            min_positive_psms: 20,
            min_decoy_psms: 20,
            ..PercolatorSettings::default()
        };
        score_psms(&mut features, &settings, Tolerance::Ppm(-10.0, 10.0)).unwrap();
        let target_mean = features
            .iter()
            .filter(|psm| psm.label == 1)
            .map(|psm| psm.discriminant_score as f64)
            .sum::<f64>()
            / 600.0;
        let decoy_mean = features
            .iter()
            .filter(|psm| psm.label == -1)
            .map(|psm| psm.discriminant_score as f64)
            .sum::<f64>()
            / 300.0;
        assert!(target_mean > decoy_mean);
    }

    #[test]
    fn svm_learns_separable_psms_out_of_fold() {
        assert_learns_separable_psms(PercolatorModel::Svm);
    }

    #[test]
    fn perpetual_learns_separable_psms_out_of_fold() {
        assert_learns_separable_psms(PercolatorModel::Perpetual);
    }
}
