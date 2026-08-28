//! Retention time prediction using linear regression
//!
//! See Klammer et al., Anal. Chem. 2007, 79, 16, 6111–6118
//! <https://doi.org/10.1021/ac070262k>

use super::regression::LinearRegression;
use crate::database::IndexedDatabase;
use crate::enzyme::Position;
use crate::mass::VALID_AA;
use crate::modification::ModificationSpecificity;
use crate::peptide::Peptide;
use crate::scoring::Feature;
use rayon::prelude::*;
use serde::{Deserialize, Serialize};

#[derive(
    Clone, Copy, Debug, Default, Serialize, Deserialize, schemars::JsonSchema, PartialEq, Eq,
)]
#[serde(rename_all = "snake_case")]
pub enum RetentionTimeFeatureSet {
    #[default]
    Basic,
    Physicochemical,
    AdditivePtm,
}

#[derive(Clone, Debug, Serialize, Deserialize, schemars::JsonSchema)]
#[serde(default, deny_unknown_fields)]
pub struct RetentionTimeSettings {
    /// Sequence features used by the linear regression model.
    pub features: RetentionTimeFeatureSet,
    /// Number of peptide-grouped cross-validation folds.
    #[schemars(range(min = 2, max = 10))]
    pub folds: usize,
    /// Seed used for deterministic peptide-level fold assignment.
    pub seed: u64,
    /// Ridge penalty for additive variable-PTM retention-time offsets.
    #[schemars(range(min = 0.0))]
    pub ptm_regularization: f64,
}

impl Default for RetentionTimeSettings {
    fn default() -> Self {
        Self {
            features: RetentionTimeFeatureSet::Basic,
            folds: 3,
            seed: 42,
            ptm_regularization: 25.0,
        }
    }
}

/// Try to fit a retention time prediction model
pub fn predict(
    db: &IndexedDatabase,
    features: &mut [Feature],
    settings: &RetentionTimeSettings,
) -> Option<()> {
    predict_linear(db, features, settings)
}

fn predict_linear(
    db: &IndexedDatabase,
    features: &mut [Feature],
    settings: &RetentionTimeSettings,
) -> Option<()> {
    // Training LR might fail - not enough values, or r-squared is < 0.7
    let lr = RetentionModel::fit(db, features, settings.features)?;
    features.par_iter_mut().for_each(|feat| {
        // LR can sometimes predict crazy values - clamp predicted RT
        let rt = lr.predict_peptide(db, feat);
        let bounded = rt.clamp(0.0, 1.0) as f32;
        feat.predicted_rt = bounded;
        feat.delta_rt_model = (feat.aligned_rt - bounded).abs();
    });
    if settings.features == RetentionTimeFeatureSet::AdditivePtm {
        apply_ptm_offsets(db, features, settings);
    }
    Some(())
}
pub struct RetentionModel {
    beta: Vec<f64>,
    map: [usize; 26],
    feature_set: RetentionTimeFeatureSet,
    pub r2: f64,
}

const BASIC_FEATURES: usize = VALID_AA.len() * 3 + 3;
const N_TERMINAL: usize = VALID_AA.len();
const C_TERMINAL: usize = VALID_AA.len() * 2;
const PEPTIDE_LEN: usize = BASIC_FEATURES - 3;
const PEPTIDE_MASS: usize = BASIC_FEATURES - 2;
const INTERCEPT: usize = BASIC_FEATURES - 1;

const TERMINAL_POSITIONS: usize = VALID_AA.len() * 4;
const HYDROPHOBIC_FEATURES: usize = 9;
const PROPERTY_FEATURES: usize = 8;
const MODIFICATION_FEATURES: usize = 4;
const PHYSICOCHEMICAL_SOURCE_FEATURES: usize = BASIC_FEATURES
    + TERMINAL_POSITIONS
    + HYDROPHOBIC_FEATURES
    + PROPERTY_FEATURES
    + MODIFICATION_FEATURES;
const POSITIONAL_START: usize = BASIC_FEATURES;
const HYDROPHOBIC_START: usize = POSITIONAL_START + TERMINAL_POSITIONS;
const PROPERTY_START: usize = HYDROPHOBIC_START + HYDROPHOBIC_FEATURES;
const MODIFICATION_START: usize = PROPERTY_START + PROPERTY_FEATURES;

// Reference-code one amino acid in each group to avoid exact dummy-variable
// collinearity in ordinary least squares. The omitted O residue is recoverable
// from peptide length or the one-hot group's intercept.
const LINEAR_AA_FEATURES: usize = VALID_AA.len() - 1;
const LINEAR_POSITIONAL_FEATURES: usize = LINEAR_AA_FEATURES * 4;
const LINEAR_HYDROPHOBIC_FEATURES: usize = 4;
const LINEAR_PROPERTY_FEATURES: usize = 7;
const LINEAR_MODIFICATION_FEATURES: usize = 4;
const LINEAR_GLOBAL_FEATURES: usize = 3;
const LINEAR_PHYSICOCHEMICAL_FEATURES: usize = LINEAR_AA_FEATURES
    + LINEAR_POSITIONAL_FEATURES
    + LINEAR_HYDROPHOBIC_FEATURES
    + LINEAR_PROPERTY_FEATURES
    + LINEAR_MODIFICATION_FEATURES
    + LINEAR_GLOBAL_FEATURES;
const LINEAR_ADDITIVE_PTM_FEATURES: usize = LINEAR_AA_FEATURES
    + LINEAR_AA_FEATURES * 2
    + LINEAR_HYDROPHOBIC_FEATURES
    + LINEAR_PROPERTY_FEATURES
    + LINEAR_GLOBAL_FEATURES;

impl RetentionModel {
    /// One-hot encoding of peptide sequences into feature vector
    /// Note that this currently does not take into account any modifications
    fn embed(peptide: &Peptide, map: &[usize; 26]) -> [f64; BASIC_FEATURES] {
        let mut embedding = [0.0; BASIC_FEATURES];
        let cterm = peptide.sequence.len().saturating_sub(3);
        for (aa_idx, residue) in peptide.sequence.iter().enumerate() {
            let idx = map[(residue - b'A') as usize];
            embedding[idx] += 1.0;
            // Embed N- and C-terminal AA's (2 on each end, excluding K/R)
            match aa_idx {
                0 | 1 => embedding[N_TERMINAL + idx] += 1.0,
                x if x == cterm || x == cterm + 1 => embedding[C_TERMINAL + idx] += 1.0,
                _ => {}
            }
        }
        embedding[PEPTIDE_LEN] = peptide.sequence.len() as f64;
        embedding[PEPTIDE_MASS] = (peptide.monoisotopic as f64).ln_1p();
        embedding[INTERCEPT] = 1.0;
        embedding
    }

    /// Attempt to fit a linear regression model: peptide sequence ~ retention time
    pub fn fit(
        db: &IndexedDatabase,
        training_set: &[Feature],
        feature_set: RetentionTimeFeatureSet,
    ) -> Option<Self> {
        // Create a mapping from amino acid character to vector embedding
        let mut map = [0; 26];
        for (idx, aa) in VALID_AA.iter().enumerate() {
            map[(aa - b'A') as usize] = idx;
        }

        let lr = match feature_set {
            RetentionTimeFeatureSet::Basic => LinearRegression::fit::<_, BASIC_FEATURES>(
                training_set,
                |feat| feat.label == 1 && feat.spectrum_q <= 0.01,
                |psm| Self::embed(&db[psm.peptide_idx], &map),
                |psm| psm.aligned_rt as f64,
            )?,
            RetentionTimeFeatureSet::Physicochemical => {
                LinearRegression::fit::<_, LINEAR_PHYSICOCHEMICAL_FEATURES>(
                    training_set,
                    |feat| feat.label == 1 && feat.spectrum_q <= 0.01,
                    |psm| linear_physicochemical_embed(&db[psm.peptide_idx], &map),
                    |psm| psm.aligned_rt as f64,
                )?
            }
            RetentionTimeFeatureSet::AdditivePtm => {
                LinearRegression::fit::<_, LINEAR_ADDITIVE_PTM_FEATURES>(
                    training_set,
                    |feat| feat.label == 1 && feat.spectrum_q <= 0.01,
                    |psm| linear_additive_ptm_embed(&db[psm.peptide_idx], &map),
                    |psm| psm.aligned_rt as f64,
                )?
            }
        };

        log::info!(
            "- fit {:?} linear retention time model, rsq = {}",
            feature_set,
            lr.r2
        );
        Some(Self {
            beta: lr.beta,
            map,
            feature_set,
            r2: lr.r2,
        })
    }

    /// Predict retention times for a collection of PSMs
    pub fn predict_peptide(&self, db: &IndexedDatabase, psm: &Feature) -> f64 {
        match self.feature_set {
            RetentionTimeFeatureSet::Basic => Self::embed(&db[psm.peptide_idx], &self.map)
                .into_iter()
                .zip(&self.beta)
                .fold(0.0f64, |sum, (x, y)| sum + x * y),
            RetentionTimeFeatureSet::Physicochemical => {
                linear_physicochemical_embed(&db[psm.peptide_idx], &self.map)
                    .into_iter()
                    .zip(&self.beta)
                    .fold(0.0f64, |sum, (x, y)| sum + x * y)
            }
            RetentionTimeFeatureSet::AdditivePtm => {
                linear_additive_ptm_embed(&db[psm.peptide_idx], &self.map)
                    .into_iter()
                    .zip(&self.beta)
                    .fold(0.0f64, |sum, (x, y)| sum + x * y)
            }
        }
    }
}

fn hydrophobicity(residue: u8) -> f64 {
    // Kyte-Doolittle hydropathy scale. Non-canonical U/O are neutral.
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

fn physicochemical_source_embed(
    peptide: &Peptide,
    map: &[usize; 26],
    feature_set: RetentionTimeFeatureSet,
) -> [f64; PHYSICOCHEMICAL_SOURCE_FEATURES] {
    let mut embedding = [0.0; PHYSICOCHEMICAL_SOURCE_FEATURES];
    embedding[..BASIC_FEATURES].copy_from_slice(&RetentionModel::embed(peptide, map));
    if feature_set == RetentionTimeFeatureSet::Basic {
        return embedding;
    }

    let len = peptide.sequence.len();
    let positions = [0, 1, len.saturating_sub(2), len.saturating_sub(1)];
    for (position, sequence_idx) in positions.into_iter().enumerate() {
        if let Some(&residue) = peptide.sequence.get(sequence_idx) {
            let aa = map[(residue - b'A') as usize];
            embedding[POSITIONAL_START + position * VALID_AA.len() + aa] = 1.0;
        }
    }

    let mut hydro_sum = 0.0;
    let mut hydro_sq_sum = 0.0;
    let mut hydro_min = f64::INFINITY;
    let mut hydro_max = f64::NEG_INFINITY;
    let mut bins = [0.0; 4];
    let mut bin_counts = [0usize; 4];
    let mut positive = 0usize;
    let mut negative = 0usize;
    let mut polar = 0usize;
    let mut aromatic = 0usize;
    let mut aliphatic = 0usize;
    let mut proline = 0usize;
    let mut glycine = 0usize;

    for (idx, &residue) in peptide.sequence.iter().enumerate() {
        let hydro = hydrophobicity(residue);
        hydro_sum += hydro;
        hydro_sq_sum += hydro * hydro;
        hydro_min = hydro_min.min(hydro);
        hydro_max = hydro_max.max(hydro);
        let bin = (idx * 4 / len.max(1)).min(3);
        bins[bin] += hydro;
        bin_counts[bin] += 1;
        positive += usize::from(matches!(residue, b'K' | b'R' | b'H'));
        negative += usize::from(matches!(residue, b'D' | b'E'));
        polar += usize::from(matches!(residue, b'S' | b'T' | b'N' | b'Q' | b'C' | b'Y'));
        aromatic += usize::from(matches!(residue, b'F' | b'W' | b'Y'));
        aliphatic += usize::from(matches!(residue, b'I' | b'L' | b'V'));
        proline += usize::from(residue == b'P');
        glycine += usize::from(residue == b'G');
    }

    let len_f = len.max(1) as f64;
    let hydro_mean = hydro_sum / len_f;
    embedding[HYDROPHOBIC_START] = hydro_sum;
    embedding[HYDROPHOBIC_START + 1] = hydro_mean;
    embedding[HYDROPHOBIC_START + 2] = hydro_min;
    embedding[HYDROPHOBIC_START + 3] = hydro_max;
    embedding[HYDROPHOBIC_START + 4] = (hydro_sq_sum / len_f - hydro_mean * hydro_mean)
        .max(0.0)
        .sqrt();
    for bin in 0..4 {
        embedding[HYDROPHOBIC_START + 5 + bin] = bins[bin] / bin_counts[bin].max(1) as f64;
    }

    embedding[PROPERTY_START] = positive as f64 / len_f;
    embedding[PROPERTY_START + 1] = negative as f64 / len_f;
    embedding[PROPERTY_START + 2] = (positive as f64 - negative as f64) / len_f;
    embedding[PROPERTY_START + 3] = polar as f64 / len_f;
    embedding[PROPERTY_START + 4] = aromatic as f64 / len_f;
    embedding[PROPERTY_START + 5] = aliphatic as f64 / len_f;
    embedding[PROPERTY_START + 6] = proline as f64 / len_f;
    embedding[PROPERTY_START + 7] = glycine as f64 / len_f;

    let modification_count = peptide
        .sequence
        .iter()
        .enumerate()
        .filter(|(index, _)| peptide.modification_at(*index) != 0.0)
        .count()
        + usize::from(peptide.nterm.is_some())
        + usize::from(peptide.cterm.is_some());
    let modification_mass = (0..peptide.sequence.len())
        .map(|index| peptide.modification_at(index))
        .sum::<f32>()
        + peptide.nterm.unwrap_or_default()
        + peptide.cterm.unwrap_or_default();
    embedding[MODIFICATION_START] = modification_count as f64 / len_f;
    embedding[MODIFICATION_START + 1] = modification_mass as f64 / 100.0;
    embedding[MODIFICATION_START + 2] = peptide.nterm.unwrap_or_default() as f64 / 100.0;
    embedding[MODIFICATION_START + 3] = peptide.cterm.unwrap_or_default() as f64 / 100.0;
    embedding
}

fn linear_physicochemical_embed(
    peptide: &Peptide,
    map: &[usize; 26],
) -> [f64; LINEAR_PHYSICOCHEMICAL_FEATURES] {
    let source =
        physicochemical_source_embed(peptide, map, RetentionTimeFeatureSet::Physicochemical);
    let mut embedding = [0.0; LINEAR_PHYSICOCHEMICAL_FEATURES];
    let mut output = 0;

    // Amino-acid counts with O as the reference category.
    for value in source.iter().take(LINEAR_AA_FEATURES) {
        embedding[output] = *value;
        output += 1;
    }
    // Exact N1/N2/C2/C1 one-hot values, also reference-coded.
    for position in 0..4 {
        let start = POSITIONAL_START + position * VALID_AA.len();
        for value in source.iter().skip(start).take(LINEAR_AA_FEATURES) {
            embedding[output] = *value;
            output += 1;
        }
    }
    // Four sequence-order-aware hydrophobicity-bin means. Total and mean
    // hydrophobicity are omitted because residue counts already encode them.
    for value in source
        .iter()
        .skip(HYDROPHOBIC_START + 5)
        .take(LINEAR_HYDROPHOBIC_FEATURES)
    {
        embedding[output] = *value;
        output += 1;
    }
    // Positive, negative, polar, aromatic, aliphatic, proline and glycine
    // fractions. Net charge is omitted because it is positive minus negative.
    for property in [0, 1, 3, 4, 5, 6, 7] {
        embedding[output] = source[PROPERTY_START + property];
        output += 1;
    }
    for value in source
        .iter()
        .skip(MODIFICATION_START)
        .take(LINEAR_MODIFICATION_FEATURES)
    {
        embedding[output] = *value;
        output += 1;
    }
    embedding[output] = source[PEPTIDE_LEN];
    embedding[output + 1] = source[PEPTIDE_MASS];
    embedding[output + 2] = 1.0;
    debug_assert_eq!(output + LINEAR_GLOBAL_FEATURES, embedding.len());
    embedding
}

fn linear_additive_ptm_embed(
    peptide: &Peptide,
    map: &[usize; 26],
) -> [f64; LINEAR_ADDITIVE_PTM_FEATURES] {
    let source =
        physicochemical_source_embed(peptide, map, RetentionTimeFeatureSet::Physicochemical);
    let mut embedding = [0.0; LINEAR_ADDITIVE_PTM_FEATURES];
    let mut output = 0;

    for value in source.iter().take(LINEAR_AA_FEATURES) {
        embedding[output] = *value;
        output += 1;
    }
    // Only the first and last residues are position encoded.
    for position in [0, 3] {
        let start = POSITIONAL_START + position * VALID_AA.len();
        for value in source.iter().skip(start).take(LINEAR_AA_FEATURES) {
            embedding[output] = *value;
            output += 1;
        }
    }
    for value in source
        .iter()
        .skip(HYDROPHOBIC_START + 5)
        .take(LINEAR_HYDROPHOBIC_FEATURES)
    {
        embedding[output] = *value;
        output += 1;
    }
    for property in [0, 1, 3, 4, 5, 6, 7] {
        embedding[output] = source[PROPERTY_START + property];
        output += 1;
    }
    embedding[output] = source[PEPTIDE_LEN];
    embedding[output + 1] = source[PEPTIDE_MASS];
    embedding[output + 2] = 1.0;
    debug_assert_eq!(output + LINEAR_GLOBAL_FEATURES, embedding.len());
    embedding
}

fn mass_matches(observed: Option<f32>, expected: f32) -> bool {
    observed
        .map(|mass| (mass - expected).abs() <= 1e-3)
        .unwrap_or(false)
}

pub(crate) fn variable_mod_count(
    peptide: &Peptide,
    specificity: ModificationSpecificity,
    mass: f32,
) -> f64 {
    let first = peptide.sequence.first().copied();
    let last = peptide.sequence.last().copied();
    let first_mass = nonzero_modification(peptide.modification_at(0));
    let last_mass =
        nonzero_modification(peptide.modification_at(peptide.sequence.len().saturating_sub(1)));
    match specificity {
        ModificationSpecificity::PeptideN(None) => {
            usize::from(mass_matches(peptide.nterm, mass)) as f64
        }
        ModificationSpecificity::PeptideC(None) => {
            usize::from(mass_matches(peptide.cterm, mass)) as f64
        }
        ModificationSpecificity::ProteinN(None)
            if matches!(peptide.position, Position::Nterm | Position::Full) =>
        {
            usize::from(mass_matches(peptide.nterm, mass)) as f64
        }
        ModificationSpecificity::ProteinC(None)
            if matches!(peptide.position, Position::Cterm | Position::Full) =>
        {
            usize::from(mass_matches(peptide.cterm, mass)) as f64
        }
        ModificationSpecificity::PeptideN(Some(residue)) if first == Some(residue) => {
            usize::from(mass_matches(first_mass, mass)) as f64
        }
        ModificationSpecificity::PeptideC(Some(residue)) if last == Some(residue) => {
            usize::from(mass_matches(last_mass, mass)) as f64
        }
        ModificationSpecificity::ProteinN(Some(residue))
            if first == Some(residue)
                && matches!(peptide.position, Position::Nterm | Position::Full) =>
        {
            usize::from(mass_matches(first_mass, mass)) as f64
        }
        ModificationSpecificity::ProteinC(Some(residue))
            if last == Some(residue)
                && matches!(peptide.position, Position::Cterm | Position::Full) =>
        {
            usize::from(mass_matches(last_mass, mass)) as f64
        }
        ModificationSpecificity::Residue(residue) => peptide
            .sequence
            .iter()
            .enumerate()
            .filter(|(index, aa)| {
                **aa == residue && (peptide.modification_at(*index) - mass).abs() <= 1e-3
            })
            .count() as f64,
        _ => 0.0,
    }
}

fn nonzero_modification(mass: f32) -> Option<f32> {
    (mass != 0.0).then_some(mass)
}

#[derive(Clone)]
struct PtmOffsetModel {
    keys: Vec<(ModificationSpecificity, f32)>,
    offsets: Vec<f64>,
}

impl PtmOffsetModel {
    fn fit(
        db: &IndexedDatabase,
        features: &[Feature],
        indices: &[usize],
        regularization: f64,
    ) -> Option<Self> {
        use super::{gauss::Gauss, matrix::Matrix as SageMatrix};

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

        let dimensions = keys.len();
        let mut covariance = vec![0.0; dimensions * dimensions];
        let mut response = vec![0.0; dimensions];
        for &idx in indices {
            let peptide = &db[features[idx].peptide_idx];
            let row = keys
                .iter()
                .map(|&(specificity, mass)| variable_mod_count(peptide, specificity, mass))
                .collect::<Vec<_>>();
            let residual = features[idx].aligned_rt as f64 - features[idx].predicted_rt as f64;
            for column in 0..dimensions {
                response[column] += row[column] * residual;
                for other in 0..dimensions {
                    covariance[column * dimensions + other] += row[column] * row[other];
                }
            }
        }
        for diagonal in 0..dimensions {
            covariance[diagonal * dimensions + diagonal] += regularization;
        }
        let offsets = Gauss::solve(
            SageMatrix::new(covariance, dimensions, dimensions),
            SageMatrix::col_vector(response),
        )?
        .take();
        Some(Self { keys, offsets })
    }

    fn predict(&self, peptide: &Peptide) -> f64 {
        self.keys
            .iter()
            .zip(&self.offsets)
            .map(|(&(specificity, mass), offset)| {
                variable_mod_count(peptide, specificity, mass) * offset
            })
            .sum()
    }
}

fn apply_ptm_offsets(
    db: &IndexedDatabase,
    features: &mut [Feature],
    settings: &RetentionTimeSettings,
) {
    if !(2..=10).contains(&settings.folds)
        || !settings.ptm_regularization.is_finite()
        || settings.ptm_regularization <= 0.0
    {
        log::warn!("invalid additive PTM settings; retaining the sequence-only RT model");
        return;
    }
    let training = features
        .iter()
        .enumerate()
        .filter_map(|(idx, feature)| {
            (feature.label == 1 && feature.spectrum_q <= 0.01).then_some(idx)
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
    let mut fitted = 0usize;
    for held_out in 0..settings.folds {
        let train_indices = training
            .iter()
            .copied()
            .filter(|&idx| assignments[idx] != held_out)
            .collect::<Vec<_>>();
        let Some(model) =
            PtmOffsetModel::fit(db, features, &train_indices, settings.ptm_regularization)
        else {
            continue;
        };
        fitted += 1;
        for (idx, &fold) in assignments.iter().enumerate() {
            if fold == held_out {
                corrections[idx] = model.predict(&db[features[idx].peptide_idx]);
            }
        }
    }
    if fitted != settings.folds {
        log::warn!("additive PTM offset fitting failed for one or more folds; retaining zero offsets there");
    }
    features
        .par_iter_mut()
        .zip(corrections.into_par_iter())
        .for_each(|(feature, correction)| {
            let predicted = (feature.predicted_rt as f64 + correction).clamp(0.0, 1.0) as f32;
            feature.predicted_rt = predicted;
            feature.delta_rt_model = (feature.aligned_rt - predicted).abs();
        });
    log::info!(
        "- fit cross-validated additive variable-PTM offsets with ridge penalty {}",
        settings.ptm_regularization
    );
}

pub(crate) fn peptide_fold(sequence: &[u8], folds: usize, seed: u64) -> usize {
    // FNV-1a keeps identical sequences together, including peptide records that
    // differ only by a modification the current embedding does not represent.
    let mut hash = 0xcbf2_9ce4_8422_2325u64 ^ seed;
    for residue in sequence {
        hash ^= *residue as u64;
        hash = hash.wrapping_mul(0x0000_0100_0000_01b3);
    }
    (hash % folds as u64) as usize
}

#[cfg(test)]
#[path = "../../tests/unit/ml/retention_model.rs"]
mod tests;
