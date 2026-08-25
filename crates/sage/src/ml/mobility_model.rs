//! Peptide ion-mobility prediction using cross-fitted linear regression.

use super::regression::LinearRegression;
use super::retention_model::{peptide_fold, variable_mod_count};
use crate::database::IndexedDatabase;
use crate::mass::VALID_AA;
use crate::modification::ModificationSpecificity;
use crate::peptide::Peptide;
use crate::scoring::Feature;
use rayon::prelude::*;
use serde::{Deserialize, Serialize};

#[derive(Clone, Copy, Debug, Default, Serialize, Deserialize, PartialEq, Eq)]
#[serde(rename_all = "snake_case")]
pub enum IonMobilityFeatureSet {
    #[default]
    Basic,
    Physicochemical,
    AdditivePtm,
}

#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(default)]
pub struct IonMobilitySettings {
    /// Disable ion-mobility prediction even when mobility observations are present.
    pub enabled: bool,
    /// Sequence features used by the linear regression model.
    pub features: IonMobilityFeatureSet,
    /// Number of peptide-grouped cross-validation folds.
    pub folds: usize,
    /// Seed used for deterministic peptide-level fold assignment.
    pub seed: u64,
    /// Ridge penalty for additive variable-PTM mobility offsets.
    pub ptm_regularization: f64,
    /// Minimum number of high-confidence mobility observations required for fitting.
    pub min_training_psms: usize,
}

impl Default for IonMobilitySettings {
    fn default() -> Self {
        Self {
            enabled: true,
            features: IonMobilityFeatureSet::Basic,
            folds: 3,
            seed: 42,
            ptm_regularization: 25.0,
            min_training_psms: 200,
        }
    }
}

fn valid_mobility(feature: &Feature) -> bool {
    feature.ims.is_finite() && feature.ims > 0.0
}

/// Fit ion mobility only when the input contains enough usable observations.
pub fn predict(
    db: &IndexedDatabase,
    features: &mut [Feature],
    settings: &IonMobilitySettings,
) -> Option<()> {
    if !settings.enabled {
        return None;
    }
    if !(2..=10).contains(&settings.folds) {
        log::warn!("ion-mobility folds must be between 2 and 10; skipping prediction");
        return None;
    }

    let training_count = features
        .iter()
        .filter(|feature| {
            feature.label == 1 && feature.spectrum_q <= 0.01 && valid_mobility(feature)
        })
        .count();
    if training_count == 0 {
        log::debug!("no ion-mobility observations present; skipping mobility prediction");
        return None;
    }
    if training_count < settings.min_training_psms {
        log::warn!(
            "only {training_count} high-confidence ion-mobility observations present (minimum {}); skipping prediction",
            settings.min_training_psms
        );
        return None;
    }

    let requested = settings.features;
    let (fit, fitted_features) = match fit_feature_set(db, features, settings, requested) {
        Ok(fit) => (fit, requested),
        Err(error) if requested != IonMobilityFeatureSet::Basic => {
            log::warn!(
                "{requested:?} ion-mobility model failed ({error}); falling back to basic features"
            );
            match fit_feature_set(db, features, settings, IonMobilityFeatureSet::Basic) {
                Ok(fit) => (fit, IonMobilityFeatureSet::Basic),
                Err(error) => {
                    log::warn!("basic ion-mobility fallback failed: {error}");
                    return None;
                }
            }
        }
        Err(error) => {
            log::warn!("ion-mobility model failed: {error}");
            return None;
        }
    };

    let (mut predictions, base_r2, base_mae) = fit;
    log::info!(
        "- fit cross-validated {fitted_features:?} ion-mobility model, rsq = {base_r2:.4}, mae = {base_mae:.4}"
    );

    if requested == IonMobilityFeatureSet::AdditivePtm {
        let sequence_only = predictions.clone();
        match apply_ptm_offsets(db, features, settings, &mut predictions) {
            Ok((r2, mae)) if r2 >= base_r2 && mae <= base_mae => {
                log::info!(
                    "- applied cross-validated additive variable-PTM mobility offsets, rsq = {r2:.4}, mae = {mae:.4}"
                )
            }
            Ok((r2, mae)) => {
                predictions = sequence_only;
                log::info!(
                    "- rejected additive variable-PTM mobility offsets (rsq = {r2:.4}, mae = {mae:.4}); cross-validated metrics did not both improve"
                );
            }
            Err(error) => log::warn!(
                "additive PTM mobility offsets failed ({error}); retaining the sequence-only predictions"
            ),
        }
    }

    features
        .par_iter_mut()
        .zip(predictions.into_par_iter())
        .for_each(|(feature, prediction)| {
            if valid_mobility(feature) && prediction.is_finite() {
                let bounded = prediction.clamp(0.0, 2.0) as f32;
                feature.predicted_ims = bounded;
                feature.delta_ims_model = (feature.ims - bounded).abs();
            }
        });
    Some(())
}

type FitResult = (Vec<f64>, f64, f64);

fn fit_feature_set(
    db: &IndexedDatabase,
    features: &[Feature],
    settings: &IonMobilitySettings,
    feature_set: IonMobilityFeatureSet,
) -> Result<FitResult, String> {
    match feature_set {
        IonMobilityFeatureSet::Basic => {
            cross_fit::<BASIC_FEATURES>(db, features, settings, basic_embed)
        }
        IonMobilityFeatureSet::Physicochemical | IonMobilityFeatureSet::AdditivePtm => {
            cross_fit::<ENRICHED_FEATURES>(db, features, settings, enriched_embed)
        }
    }
}

fn cross_fit<const D: usize>(
    db: &IndexedDatabase,
    features: &[Feature],
    settings: &IonMobilitySettings,
    embed: impl Fn(&Peptide, u8, &[usize; 26]) -> [f64; D] + Copy + Sync,
) -> Result<FitResult, String> {
    let map = amino_acid_map();
    let training = features
        .iter()
        .enumerate()
        .filter_map(|(idx, feature)| {
            (feature.label == 1 && feature.spectrum_q <= 0.01 && valid_mobility(feature))
                .then_some(idx)
        })
        .collect::<Vec<_>>();
    let assignments = features
        .iter()
        .map(|feature| {
            peptide_fold(
                &db[feature.peptide_idx].sequence,
                settings.folds,
                settings.seed,
            )
        })
        .collect::<Vec<_>>();
    let mut predictions = vec![f64::NAN; features.len()];

    for held_out in 0..settings.folds {
        let train = training
            .iter()
            .copied()
            .filter(|&idx| assignments[idx] != held_out)
            .map(|idx| &features[idx])
            .collect::<Vec<_>>();
        if train.len() <= D {
            return Err(format!(
                "fold {held_out} has {} training observations for {D} coefficients",
                train.len()
            ));
        }
        let model = LinearRegression::fit::<_, D>(
            &train,
            |_| true,
            |feature| embed(&db[feature.peptide_idx], feature.charge, &map),
            |feature| feature.ims as f64,
        )
        .ok_or_else(|| format!("linear fit failed for fold {held_out}"))?;

        for (idx, feature) in features.iter().enumerate() {
            if assignments[idx] == held_out && valid_mobility(feature) {
                let row = embed(&db[feature.peptide_idx], feature.charge, &map);
                predictions[idx] = row
                    .iter()
                    .zip(&model.beta)
                    .map(|(value, weight)| value * weight)
                    .sum();
            }
        }
    }

    if training.iter().any(|&idx| !predictions[idx].is_finite()) {
        return Err("cross-validation did not predict every training observation".into());
    }
    let (r2, mae) = prediction_metrics(features, &predictions, &training);
    if !r2.is_finite() || !mae.is_finite() {
        return Err("cross-validated prediction metrics were not finite".into());
    }
    Ok((predictions, r2, mae))
}

fn prediction_metrics(features: &[Feature], predictions: &[f64], training: &[usize]) -> (f64, f64) {
    let mean = training
        .iter()
        .map(|&idx| features[idx].ims as f64)
        .sum::<f64>()
        / training.len() as f64;
    let residual = training
        .iter()
        .map(|&idx| (features[idx].ims as f64 - predictions[idx]).powi(2))
        .sum::<f64>();
    let total = training
        .iter()
        .map(|&idx| (features[idx].ims as f64 - mean).powi(2))
        .sum::<f64>();
    let mae = training
        .iter()
        .map(|&idx| (features[idx].ims as f64 - predictions[idx]).abs())
        .sum::<f64>()
        / training.len() as f64;
    (1.0 - residual / total, mae)
}

fn amino_acid_map() -> [usize; 26] {
    let mut map = [0; 26];
    for (idx, aa) in VALID_AA.iter().enumerate() {
        map[(aa - b'A') as usize] = idx;
    }
    map
}

const BULKY_AA_IDXS: [usize; 6] = [
    b'L' as usize - b'A' as usize,
    b'V' as usize - b'A' as usize,
    b'I' as usize - b'A' as usize,
    b'F' as usize - b'A' as usize,
    b'W' as usize - b'A' as usize,
    b'Y' as usize - b'A' as usize,
];
const UNCHARGED_POLAR_AA_IDXS: [usize; 4] = [
    b'S' as usize - b'A' as usize,
    b'T' as usize - b'A' as usize,
    b'N' as usize - b'A' as usize,
    b'Q' as usize - b'A' as usize,
];
const POSITIVE_AA_IDXS: [usize; 3] = [
    b'R' as usize - b'A' as usize,
    b'K' as usize - b'A' as usize,
    b'H' as usize - b'A' as usize,
];
const NEGATIVE_AA_IDXS: [usize; 2] = [b'D' as usize - b'A' as usize, b'E' as usize - b'A' as usize];
const TINY_AA_IDXS: [usize; 3] = [
    b'G' as usize - b'A' as usize,
    0,
    b'S' as usize - b'A' as usize,
];
const BRANCHED_AA_IDXS: [usize; 3] = [
    b'L' as usize - b'A' as usize,
    b'I' as usize - b'A' as usize,
    b'V' as usize - b'A' as usize,
];

const BASIC_FEATURES: usize = VALID_AA.len() * 4 + 12;
const BASIC_PCT_START: usize = VALID_AA.len();
const BASIC_N_TERMINAL: usize = VALID_AA.len() * 2;
const BASIC_C_TERMINAL: usize = VALID_AA.len() * 3;
const BASIC_NUM_BRANCHED: usize = BASIC_FEATURES - 12;
const BASIC_NUM_TINY: usize = BASIC_FEATURES - 11;
const BASIC_NUM_UC_POLAR: usize = BASIC_FEATURES - 10;
const BASIC_NUM_BULKY: usize = BASIC_FEATURES - 9;
const BASIC_NUM_POSITIVE: usize = BASIC_FEATURES - 8;
const BASIC_NUM_NEGATIVE: usize = BASIC_FEATURES - 7;
const BASIC_INV_CHARGE: usize = BASIC_FEATURES - 6;
const BASIC_CHARGE: usize = BASIC_FEATURES - 5;
const BASIC_MZ: usize = BASIC_FEATURES - 4;
const BASIC_LEN: usize = BASIC_FEATURES - 3;
const BASIC_MASS: usize = BASIC_FEATURES - 2;
const BASIC_INTERCEPT: usize = BASIC_FEATURES - 1;

fn basic_embed(peptide: &Peptide, charge: u8, map: &[usize; 26]) -> [f64; BASIC_FEATURES] {
    let mut embedding = [0.0; BASIC_FEATURES];
    let cterm = peptide.sequence.len().saturating_sub(3);
    let length = peptide.sequence.len().max(1) as f64;
    for (aa_idx, residue) in peptide.sequence.iter().enumerate() {
        let idx = map[(residue - b'A') as usize];
        embedding[idx] += 1.0;
        match aa_idx {
            0 | 1 => embedding[BASIC_N_TERMINAL + idx] += 1.0,
            x if x > cterm => embedding[BASIC_C_TERMINAL + idx] += 1.0,
            _ => {}
        }
        embedding[BASIC_NUM_BULKY] += usize::from(BULKY_AA_IDXS.contains(&idx)) as f64;
        embedding[BASIC_NUM_UC_POLAR] += usize::from(UNCHARGED_POLAR_AA_IDXS.contains(&idx)) as f64;
        embedding[BASIC_NUM_POSITIVE] += usize::from(POSITIVE_AA_IDXS.contains(&idx)) as f64;
        embedding[BASIC_NUM_NEGATIVE] += usize::from(NEGATIVE_AA_IDXS.contains(&idx)) as f64;
        embedding[BASIC_NUM_TINY] += usize::from(TINY_AA_IDXS.contains(&idx)) as f64;
        embedding[BASIC_NUM_BRANCHED] += usize::from(BRANCHED_AA_IDXS.contains(&idx)) as f64;
    }
    for idx in 0..VALID_AA.len() {
        embedding[BASIC_PCT_START + idx] = embedding[idx] / length;
    }
    let charge = charge.max(1) as f64;
    embedding[BASIC_CHARGE] = charge;
    embedding[BASIC_INV_CHARGE] = 1.0 / charge;
    embedding[BASIC_LEN] = length;
    embedding[BASIC_MASS] = peptide.monoisotopic as f64 / 1000.0;
    embedding[BASIC_MZ] = peptide.monoisotopic as f64 / charge / 1000.0;
    embedding[BASIC_INTERCEPT] = 1.0;
    embedding
}

const ENRICHED_AA_FEATURES: usize = VALID_AA.len() - 1;
const ENRICHED_TERMINAL_FEATURES: usize = ENRICHED_AA_FEATURES * 2;
const ENRICHED_HYDROPHOBIC_FEATURES: usize = 4;
const ENRICHED_PROPERTY_FEATURES: usize = 6;
const ENRICHED_GLOBAL_FEATURES: usize = 6;
const ENRICHED_FEATURES: usize = ENRICHED_AA_FEATURES
    + ENRICHED_TERMINAL_FEATURES
    + ENRICHED_HYDROPHOBIC_FEATURES
    + ENRICHED_PROPERTY_FEATURES
    + ENRICHED_GLOBAL_FEATURES;

fn hydrophobicity(residue: u8) -> f64 {
    match residue {
        b'I' => 4.5,
        b'V' => 4.2,
        b'L' => 3.8,
        b'F' => 2.8,
        b'C' => 2.5,
        b'M' => 1.9,
        b'A' => 1.8,
        b'G' => -0.4,
        b'T' => -0.7,
        b'S' => -0.8,
        b'W' => -0.9,
        b'Y' => -1.3,
        b'P' => -1.6,
        b'H' => -3.2,
        b'E' | b'Q' | b'D' | b'N' => -3.5,
        b'K' => -3.9,
        b'R' => -4.5,
        _ => 0.0,
    }
}

fn enriched_embed(peptide: &Peptide, charge: u8, map: &[usize; 26]) -> [f64; ENRICHED_FEATURES] {
    let mut embedding = [0.0; ENRICHED_FEATURES];
    let mut output = 0;
    let length = peptide.sequence.len().max(1);

    let mut counts = [0.0; ENRICHED_AA_FEATURES];
    let mut properties = [0.0; ENRICHED_PROPERTY_FEATURES];
    let mut hydro_bins = [0.0; ENRICHED_HYDROPHOBIC_FEATURES];
    let mut hydro_counts = [0usize; ENRICHED_HYDROPHOBIC_FEATURES];
    for (idx, &residue) in peptide.sequence.iter().enumerate() {
        let aa = map[(residue - b'A') as usize];
        if aa < ENRICHED_AA_FEATURES {
            counts[aa] += 1.0;
        }
        properties[0] += usize::from(BULKY_AA_IDXS.contains(&aa)) as f64;
        properties[1] += usize::from(UNCHARGED_POLAR_AA_IDXS.contains(&aa)) as f64;
        properties[2] += usize::from(POSITIVE_AA_IDXS.contains(&aa)) as f64;
        properties[3] += usize::from(NEGATIVE_AA_IDXS.contains(&aa)) as f64;
        properties[4] += usize::from(TINY_AA_IDXS.contains(&aa)) as f64;
        properties[5] += usize::from(BRANCHED_AA_IDXS.contains(&aa)) as f64;
        let bin =
            (idx * ENRICHED_HYDROPHOBIC_FEATURES / length).min(ENRICHED_HYDROPHOBIC_FEATURES - 1);
        hydro_bins[bin] += hydrophobicity(residue);
        hydro_counts[bin] += 1;
    }
    for value in counts {
        embedding[output] = value;
        output += 1;
    }
    for sequence_idx in [0, length.saturating_sub(1)] {
        if let Some(&residue) = peptide.sequence.get(sequence_idx) {
            let aa = map[(residue - b'A') as usize];
            if aa < ENRICHED_AA_FEATURES {
                embedding[output + aa] = 1.0;
            }
        }
        output += ENRICHED_AA_FEATURES;
    }
    for (sum, count) in hydro_bins.into_iter().zip(hydro_counts) {
        embedding[output] = sum / count.max(1) as f64;
        output += 1;
    }
    for value in properties {
        embedding[output] = value / length as f64;
        output += 1;
    }
    let charge = charge.max(1) as f64;
    embedding[output] = charge;
    embedding[output + 1] = 1.0 / charge;
    embedding[output + 2] = peptide.monoisotopic as f64 / 1000.0;
    embedding[output + 3] = peptide.monoisotopic as f64 / charge / 1000.0;
    embedding[output + 4] = length as f64;
    embedding[output + 5] = 1.0;
    debug_assert_eq!(output + ENRICHED_GLOBAL_FEATURES, embedding.len());
    embedding
}

#[derive(Clone)]
struct MobilityPtmOffsetModel {
    keys: Vec<(ModificationSpecificity, f32)>,
    charges: Vec<u8>,
    offsets: Vec<f64>,
}

impl MobilityPtmOffsetModel {
    fn fit(
        db: &IndexedDatabase,
        features: &[Feature],
        predictions: &[f64],
        indices: &[usize],
        regularization: f64,
    ) -> Option<Self> {
        use super::{gauss::Gauss, matrix::Matrix};

        let mut keys = db.model_mods.clone();
        keys.sort_unstable_by(|(a_spec, a_mass), (b_spec, b_mass)| {
            a_spec.cmp(b_spec).then_with(|| a_mass.total_cmp(b_mass))
        });
        keys.dedup_by(|(a_spec, a_mass), (b_spec, b_mass)| {
            a_spec == b_spec && a_mass.to_bits() == b_mass.to_bits()
        });
        if keys.is_empty() {
            return None;
        }
        let mut charges = indices
            .iter()
            .map(|&idx| features[idx].charge)
            .collect::<Vec<_>>();
        charges.sort_unstable();
        charges.dedup();

        let dimensions = keys.len() * (charges.len() + 1);
        let mut covariance = vec![0.0; dimensions * dimensions];
        let mut response = vec![0.0; dimensions];
        for &idx in indices {
            let row = Self::row(
                &db[features[idx].peptide_idx],
                features[idx].charge,
                &keys,
                &charges,
            );
            let residual = features[idx].ims as f64 - predictions[idx];
            for column in 0..dimensions {
                response[column] += row[column] * residual;
                for other in 0..dimensions {
                    covariance[column * dimensions + other] += row[column] * row[other];
                }
            }
        }
        for diagonal in 0..dimensions {
            let charge_specific = diagonal >= keys.len();
            covariance[diagonal * dimensions + diagonal] +=
                regularization * if charge_specific { 4.0 } else { 1.0 };
        }
        let offsets = Gauss::solve(
            Matrix::new(covariance, dimensions, dimensions),
            Matrix::col_vector(response),
        )?
        .take();
        Some(Self {
            keys,
            charges,
            offsets,
        })
    }

    fn row(
        peptide: &Peptide,
        charge: u8,
        keys: &[(ModificationSpecificity, f32)],
        charges: &[u8],
    ) -> Vec<f64> {
        let counts = keys
            .iter()
            .map(|&(specificity, mass)| variable_mod_count(peptide, specificity, mass))
            .collect::<Vec<_>>();
        let mut row = Vec::with_capacity(keys.len() * (charges.len() + 1));
        row.extend(&counts);
        for &model_charge in charges {
            row.extend(
                counts
                    .iter()
                    .map(|&count| if charge == model_charge { count } else { 0.0 }),
            );
        }
        row
    }

    fn predict(&self, peptide: &Peptide, charge: u8) -> f64 {
        Self::row(peptide, charge, &self.keys, &self.charges)
            .iter()
            .zip(&self.offsets)
            .map(|(value, offset)| value * offset)
            .sum()
    }
}

fn apply_ptm_offsets(
    db: &IndexedDatabase,
    features: &[Feature],
    settings: &IonMobilitySettings,
    predictions: &mut [f64],
) -> Result<(f64, f64), String> {
    if !settings.ptm_regularization.is_finite() || settings.ptm_regularization <= 0.0 {
        return Err("ptm_regularization must be finite and greater than zero".into());
    }
    if db.model_mods.is_empty() {
        return Err("no variable modifications are configured".into());
    }
    let training = features
        .iter()
        .enumerate()
        .filter_map(|(idx, feature)| {
            (feature.label == 1 && feature.spectrum_q <= 0.01 && valid_mobility(feature))
                .then_some(idx)
        })
        .collect::<Vec<_>>();
    let assignments = features
        .iter()
        .map(|feature| {
            peptide_fold(
                &db[feature.peptide_idx].sequence,
                settings.folds,
                settings.seed,
            )
        })
        .collect::<Vec<_>>();
    let mut corrections = vec![0.0; features.len()];

    for held_out in 0..settings.folds {
        let train_indices = training
            .iter()
            .copied()
            .filter(|&idx| assignments[idx] != held_out)
            .collect::<Vec<_>>();
        let model = MobilityPtmOffsetModel::fit(
            db,
            features,
            predictions,
            &train_indices,
            settings.ptm_regularization,
        )
        .ok_or_else(|| format!("offset fit failed for fold {held_out}"))?;
        for (idx, feature) in features.iter().enumerate() {
            if assignments[idx] == held_out && valid_mobility(feature) {
                corrections[idx] = model.predict(&db[feature.peptide_idx], feature.charge);
            }
        }
    }
    for (idx, correction) in corrections.into_iter().enumerate() {
        if predictions[idx].is_finite() {
            predictions[idx] += correction;
        }
    }
    Ok(prediction_metrics(features, predictions, &training))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn ion_mobility_settings_have_safe_defaults() {
        let settings = IonMobilitySettings::default();
        assert!(settings.enabled);
        assert_eq!(settings.features, IonMobilityFeatureSet::Basic);
        assert_eq!(settings.folds, 3);
        assert_eq!(settings.ptm_regularization, 25.0);
        assert_eq!(settings.min_training_psms, 200);
    }

    #[test]
    fn enriched_embedding_is_finite_and_charge_aware() {
        let map = amino_acid_map();
        let peptide = Peptide {
            sequence: b"ACDEFGHIK".to_vec().into(),
            modifications: vec![0.0; 9],
            monoisotopic: 1000.0,
            ..Peptide::default()
        };
        let charge_two = enriched_embed(&peptide, 2, &map);
        let charge_three = enriched_embed(&peptide, 3, &map);
        assert!(charge_two.iter().all(|value| value.is_finite()));
        assert_ne!(charge_two, charge_three);
        assert_eq!(charge_two.len(), ENRICHED_FEATURES);
    }

    #[test]
    fn zero_mobility_is_not_a_training_observation() {
        let mut feature = Feature::default();
        assert!(!valid_mobility(&feature));
        feature.ims = 1.1;
        assert!(valid_mobility(&feature));
        feature.ims = f32::NAN;
        assert!(!valid_mobility(&feature));
    }

    #[test]
    fn ptm_row_has_global_and_charge_specific_effects() {
        let peptide = Peptide {
            sequence: b"AMPEPTIDEK".to_vec().into(),
            modifications: vec![0.0, 15.994_915, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
            ..Peptide::default()
        };
        let keys = [(ModificationSpecificity::Residue(b'M'), 15.994_915)];
        assert_eq!(
            MobilityPtmOffsetModel::row(&peptide, 2, &keys, &[2, 3]),
            vec![1.0, 1.0, 0.0]
        );
    }
}
