//! Candidate scoring: the fragment index search, the score itself and the reported features.

use crate::database::{
    same_peptidoform, IndexedDatabase, MassOffset, MassOffsetAssignment, PeptideIx,
};
use crate::ion_series::{IonGroupSeries, Kind};
use crate::mass::{Tolerance, PROTON};
use crate::mass_recalibration::{FileMassCorrection, MassRecalibration};
use crate::peptide::{Peptide, Site};
use crate::spectrum::ProcessedSpectrum;
use serde::{Deserialize, Serialize};
use std::borrow::Cow;
use std::cell::RefCell;
use std::collections::HashMap;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::Arc;

#[derive(Copy, Clone, Debug, Serialize, Deserialize, schemars::JsonSchema)]
pub enum ScoreType {
    SageHyperScore,
    OpenMSHyperScore,
}

/// Structure to hold temporary scores
#[derive(Copy, Clone, Default, Debug, PartialEq)]
struct Score {
    peptide: PeptideIx,
    matched_b: u16,
    matched_y: u16,
    summed_b: f32,
    summed_y: f32,
    longest_b: usize,
    longest_y: usize,
    hyperscore: f64,
    /// Fragment errors after search-time recalibration (equal to the raw
    /// errors when no correction applies).
    ppm_difference: f32,
    signed_ppm_difference: f32,
    /// Fragment errors of the observed, uncorrected peaks.
    raw_ppm_difference: f32,
    raw_signed_ppm_difference: f32,
    /// Intensity-weighted sum (then mean) of squared raw signed errors.
    raw_squared_ppm_difference: f32,
    /// Cleavages whose generic fragment-loss form matched a peak, when loss
    /// ions are scored as separate features.
    matched_loss: u16,
    summed_loss: f32,
    precursor_charge: u8,
    isotope_error: i8,
}

impl Eq for Score {}

impl PartialOrd for Score {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for Score {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.hyperscore.total_cmp(&other.hyperscore)
    }
}

/// Preliminary score - # of matched peaks for each candidate peptide
#[derive(Copy, Clone, Default, Debug, PartialEq, Eq, PartialOrd, Ord)]
struct PreScore {
    matched: u16,
    peptide: PeptideIx,
    precursor_charge: u8,
    isotope_error: i8,
    /// Zero for the unshifted search, otherwise one plus the index of the
    /// tested [`crate::database::MassOffset`].
    offset: u8,
}

/// A fully scored candidate peptidoform.
struct Candidate<'db> {
    score: Score,
    fragments: Option<Fragments>,
    coverage: Coverage,
    peptide: Cow<'db, Peptide>,
    mass_offset: Option<MassOffsetAssignment>,
}

/// A spectrum-specific mass-offset hypothesis for
/// [`Scorer::score_offset_hypotheses`]: the observed-scale precursor mass,
/// the offset carried by the searched peptide, and the fragment mass shifts
/// used for preliminary matching (for example the intact offset and the
/// stubs left when a cleavable offset fragments).
#[derive(Clone, Debug)]
pub struct OffsetHypothesis {
    pub precursor_mass: f32,
    pub precursor_charge: u8,
    pub offset: MassOffset,
    pub preliminary_shifts: Vec<f32>,
}

/// One fully scored placement for an [`OffsetHypothesis`].
#[derive(Clone, Debug)]
pub struct OffsetMatch {
    pub peptide_idx: PeptideIx,
    /// The indexed peptide with the offset placed at `site`.
    pub peptide: Peptide,
    pub site: Site,
    pub preliminary_matched: u16,
    pub hyperscore: f64,
    pub matched_b: u16,
    pub matched_y: u16,
    pub summed_b: f32,
    pub summed_y: f32,
    pub longest_b: u32,
    pub longest_y: u32,
    /// Intensity-weighted mean absolute fragment error, ppm.
    pub fragment_ppm: f32,
    pub fragments: Fragments,
}

static PSM_COUNTER: AtomicUsize = AtomicUsize::new(1);

thread_local! {
    /// Reuse dense candidate-counting storage on each Rayon worker.
    static PRE_SCORE_SCRATCH: RefCell<Vec<PreScore>> = const { RefCell::new(Vec::new()) };
}

fn increment_psm_counter() -> usize {
    PSM_COUNTER.fetch_add(1, Ordering::Relaxed)
}

/// Stirling's approximation for log factorial
fn lnfact(n: u16) -> f64 {
    if n == 0 {
        1.0
    } else {
        let n = n as f64;
        n * n.ln() - n + 0.5 * n.ln() + 0.5 * (std::f64::consts::PI * 2.0 * n).ln()
    }
}

impl ScoreType {
    pub fn score(&self, matched_b: u16, matched_y: u16, summed_b: f32, summed_y: f32) -> f64 {
        let score = match self {
            // Calculate the X!Tandem hyperscore
            Self::SageHyperScore => {
                let i = (summed_b + 1.0) as f64 * (summed_y + 1.0) as f64;

                i.ln() + lnfact(matched_b) + lnfact(matched_y)
            }
            // Calculate the OpenMS flavour hyperscore
            Self::OpenMSHyperScore => {
                let summed_intensity = summed_b + summed_y;

                summed_intensity.ln_1p() as f64 + lnfact(matched_b) + lnfact(matched_y)
            }
        };
        if score.is_finite() {
            score
        } else {
            255.0
        }
    }
}

impl Score {
    /// Calculate the hyperscore for a given PSM choosing between implementations based on `score_type`
    fn hyperscore(&self, score_type: ScoreType) -> f64 {
        score_type.score(self.matched_b, self.matched_y, self.summed_b, self.summed_y)
    }
}

#[derive(Clone)]
pub struct Scorer<'db> {
    pub db: &'db IndexedDatabase,
    pub precursor_tol: Tolerance,
    pub fragment_tol: Tolerance,
    /// What is the minimum number of matched b and y ion peaks to report PSMs for?
    pub min_matched_peaks: u16,
    /// Precursor isotope error lower bounds (e.g. -1)
    pub min_isotope_err: i8,
    /// Precursor isotope error upper bounds (e.g. 3)
    pub max_isotope_err: i8,
    pub min_precursor_charge: u8,
    pub max_precursor_charge: u8,
    pub override_precursor_charge: bool,
    pub max_fragment_charge: Option<u8>,
    pub chimera: bool,
    pub report_psms: usize,

    // Rather than use a fixed precursor tolerance, dynamically alter
    // the precursor tolerance window based on MS2 isolation window and charge
    pub wide_window: bool,
    pub annotate_matches: bool,
    /// A precursor delta mass (`expmass - calcmass`) within this many ppm of the
    /// calculated mass is treated as no shift for sequence-ambiguity annotation.
    pub mass_shift_ppm: f32,
    pub score_type: ScoreType,
    /// Per-file search-time mass corrections. Spectra are corrected before
    /// matching; reported `expmass`, `delta_mass`, `average_ppm`,
    /// `signed_fragment_ppm` and fragment `mz_experimental` stay raw.
    pub mass_recalibration: Option<Arc<MassRecalibration>>,
}

#[inline(always)]
/// Calculate upper bound (excluded) of the charge state range to use for
/// searching fragment ions (1..N)
/// If user has configured max_fragment_charge, potentially override precursor
/// charge
pub(crate) fn max_fragment_charge(max_fragment_charge: Option<u8>, precursor_charge: u8) -> u8 {
    precursor_charge
        .min(
            max_fragment_charge
                .map(|c| c + 1)
                .unwrap_or(precursor_charge),
        )
        .max(2)
}

/// Translate an observed precursor query by a mass offset `delta`. The
/// tolerance is evaluated at the observed (modified) precursor mass and then
/// shifted, so the searched base-peptide window keeps the configured width.
pub(crate) fn offset_query(mass: f32, tolerance: Tolerance, delta: f32) -> (f32, Tolerance) {
    let (lo, hi) = tolerance.bounds(mass);
    (mass - delta, Tolerance::Da(lo - mass, hi - mass))
}

mod chimera;
mod feature;
mod matching;
mod standard;

pub use feature::*;
pub use matching::*;

impl<'db> Scorer<'db> {
    /// This scorer's settings against another database, e.g. a module's own
    /// candidate set. Override single settings with struct update syntax:
    /// `Scorer { chimera: false, ..scorer.with_db(&db) }`.
    pub fn with_db<'other>(&self, db: &'other IndexedDatabase) -> Scorer<'other> {
        // Listed field by field (the lifetime changes), so a new setting
        // fails to compile here until it is carried over.
        Scorer {
            db,
            precursor_tol: self.precursor_tol,
            fragment_tol: self.fragment_tol,
            min_matched_peaks: self.min_matched_peaks,
            min_isotope_err: self.min_isotope_err,
            max_isotope_err: self.max_isotope_err,
            min_precursor_charge: self.min_precursor_charge,
            max_precursor_charge: self.max_precursor_charge,
            override_precursor_charge: self.override_precursor_charge,
            max_fragment_charge: self.max_fragment_charge,
            chimera: self.chimera,
            report_psms: self.report_psms,
            wide_window: self.wide_window,
            annotate_matches: self.annotate_matches,
            mass_shift_ppm: self.mass_shift_ppm,
            score_type: self.score_type,
            mass_recalibration: self.mass_recalibration.clone(),
        }
    }
    /// Offset hypotheses searched for every precursor: zero is the ordinary
    /// search, followed by one hypothesis per configured mass offset.
    fn offsets(&self) -> std::ops::RangeInclusive<u8> {
        0..=self.db.mass_offsets.len() as u8
    }

    /// Translate a precursor query for an offset hypothesis. The tolerance is
    /// evaluated at the observed (modified) precursor mass and then shifted,
    /// so the searched base-peptide window keeps the configured width.
    fn offset_query(&self, mass: f32, tolerance: Tolerance, offset: u8) -> (f32, Tolerance) {
        if offset == 0 {
            return (mass, tolerance);
        }
        offset_query(
            mass,
            tolerance,
            self.db.mass_offsets[offset as usize - 1].mass(),
        )
    }

    fn fragment_shift(&self, offset: u8) -> Option<f32> {
        (offset > 0).then(|| self.db.mass_offsets[offset as usize - 1].fragment_shift())
    }
    /// The correction applied to spectra from `file_id`, if any.
    #[inline]
    fn correction(&self, file_id: usize) -> Option<&FileMassCorrection> {
        self.mass_recalibration
            .as_deref()
            .and_then(|recalibration| recalibration.file(file_id))
    }

    /// The spectrum as it is matched: corrected when its file has a
    /// search-time mass correction, otherwise borrowed unchanged.
    pub fn recalibrated<'a>(&self, query: &'a ProcessedSpectrum) -> Cow<'a, ProcessedSpectrum> {
        match self
            .mass_recalibration
            .as_deref()
            .and_then(|recalibration| recalibration.recalibrate(query))
        {
            Some(corrected) => Cow::Owned(corrected),
            None => Cow::Borrowed(query),
        }
    }

    pub fn score(&self, query: &ProcessedSpectrum) -> Vec<Feature> {
        assert_eq!(
            query.level, 2,
            "internal bug, trying to score a non-MS2 scan!"
        );
        let query = self.recalibrated(query);
        let query = query.as_ref();
        match self.chimera {
            true => self.score_chimera_fast(query),
            false => self.score_standard(query),
        }
    }
    /// Whether generic fragment-loss ions are scored, as separate rescoring
    /// features (never inside the hyperscore; see
    /// `benchmarks/FRAGMENT_LOSSES.md`).
    fn loss_features(&self) -> bool {
        self.db.fragment_losses.is_some()
    }

    /// Remove peaks matching a PSM from a query spectrum
    fn remove_matched_peaks(&self, query: &mut ProcessedSpectrum, psm: &Feature) {
        let peptide = self.db.resolve_peptide(psm);
        let fragments = self
            .db
            .ion_kinds
            .iter()
            .flat_map(|kind| {
                IonGroupSeries::with_fragment_losses(
                    &peptide,
                    *kind,
                    self.db.fragment_losses.as_deref(),
                )
            })
            // Loss peaks belong to this peptide too, so they are removed.
            .flat_map(|group| group.variants);

        let max_fragment_charge = max_fragment_charge(self.max_fragment_charge, psm.charge);
        let fragment_index = FragmentMatchIndex::new(query, max_fragment_charge);

        // Remove MS2 peaks matched by previous match
        let mut to_remove = vec![false; query.masses.len()];
        for frag in fragments {
            for charge in 1..max_fragment_charge {
                if let Some(peak_idx) = fragment_index.select_peak(
                    query,
                    frag.monoisotopic_mass,
                    charge,
                    self.fragment_tol,
                ) {
                    to_remove[peak_idx] = true;
                }
            }
        }

        let mut masses = Vec::with_capacity(query.masses.len());
        let mut intensities = Vec::with_capacity(query.intensities.len());
        let mut charges = Vec::with_capacity(query.charges.len());
        let mut charge_is_known = Vec::with_capacity(query.charge_is_known.len());
        let mut mobilities = Vec::with_capacity(query.mobilities.len());

        for (idx, removed) in to_remove.iter().enumerate() {
            if !removed {
                masses.push(query.masses[idx]);
                intensities.push(query.intensities[idx]);
                charges.push(query.charges[idx]);
                if !query.charge_is_known.is_empty() {
                    charge_is_known.push(query.charge_is_known[idx]);
                }
                if !query.mobilities.is_empty() {
                    mobilities.push(query.mobilities[idx]);
                }
            }
        }

        query.masses = masses;
        query.intensities = intensities;
        query.charges = charges;
        query.charge_is_known = charge_is_known;
        query.mobilities = mobilities;
        query.total_ion_current = query.intensities.iter().sum::<f32>();
    }
    /// Calculate full hyperscore for a resolved peptidoform
    fn score_peptide(
        &self,
        query: &ProcessedSpectrum,
        peptide: &Peptide,
        pre_score: &PreScore,
        collect_fragments: bool,
        fragment_index: &FragmentMatchIndex,
    ) -> (Score, Option<Fragments>, Coverage) {
        let mut score = Score {
            peptide: pre_score.peptide,
            precursor_charge: pre_score.precursor_charge,
            isotope_error: pre_score.isotope_error,
            ..Default::default()
        };
        let max_fragment_charge =
            max_fragment_charge(self.max_fragment_charge, score.precursor_charge);
        let fragment_correction = self
            .correction(query.file_id)
            .and_then(|correction| correction.fragment_model(query.acquisition));

        // Regenerate theoretical ions - initial database search might be
        // using only a subset of all possible ions (e.g. no b1/b2/y1/y2)
        // so we need to completely re-score this candidate
        let fragment_groups = self.db.ion_kinds.iter().flat_map(|kind| {
            IonGroupSeries::with_fragment_losses(peptide, *kind, self.db.fragment_losses.as_deref())
        });
        // Generic loss forms are scored separately as features, not as
        // alternatives of their cleavage in the hyperscore.
        let score_losses = self.loss_features();

        let mut b_run = Run::default();
        let mut y_run = Run::default();

        let mut fragments_details = Fragments::default();

        let n = peptide.sequence.len();
        let mut coverage = Coverage {
            forward: vec![0u16; n],
            reverse: vec![0u16; n],
        };

        for group in fragment_groups {
            for charge in 1..max_fragment_charge {
                // Neutral-loss forms are alternatives for the same cleavage
                // and charge. Select at most one, so extra configured variants
                // cannot inflate matched-ion counts or hyperscore factorials.
                let best = self.select_variant(
                    query,
                    fragment_index,
                    group
                        .variants
                        .iter()
                        .filter(|variant| !variant.fragment_loss),
                    charge,
                );
                if score_losses {
                    if let Some((_, peak_idx)) = self.select_variant(
                        query,
                        fragment_index,
                        group
                            .variants
                            .iter()
                            .filter(|variant| variant.fragment_loss),
                        charge,
                    ) {
                        score.matched_loss += 1;
                        score.summed_loss += query.intensities[peak_idx];
                    }
                }

                if let Some((frag, peak_idx)) = best {
                    let peak_mass = query.masses[peak_idx];
                    let peak_intensity = query.intensities[peak_idx];
                    let expected_mass = if query.has_known_charge(peak_idx) {
                        frag.monoisotopic_mass
                    } else {
                        frag.monoisotopic_mass / charge as f32
                    };

                    score.ppm_difference +=
                        peak_intensity * (expected_mass - peak_mass).abs() * 2E6
                            / (expected_mass + peak_mass);
                    score.signed_ppm_difference +=
                        peak_intensity * (peak_mass - expected_mass) * 2E6
                            / (expected_mass + peak_mass);

                    // Undo any search-time correction to report the error
                    // and m/z of the observed peak.
                    let mut exp_mz = query.peak_mz(peak_idx);
                    let raw_peak_mass = match fragment_correction {
                        Some(model) => {
                            let raw_mz = exp_mz
                                * (1.0 + model.predict_ppm(query.scan_start_time, exp_mz) * 1e-6);
                            let raw_mass =
                                peak_mass + (raw_mz - exp_mz) * query.charges[peak_idx] as f32;
                            exp_mz = raw_mz;
                            raw_mass
                        }
                        None => peak_mass,
                    };
                    score.raw_ppm_difference +=
                        peak_intensity * (expected_mass - raw_peak_mass).abs() * 2E6
                            / (expected_mass + raw_peak_mass);
                    score.raw_signed_ppm_difference +=
                        peak_intensity * (raw_peak_mass - expected_mass) * 2E6
                            / (expected_mass + raw_peak_mass);
                    let raw_signed_ppm =
                        (raw_peak_mass - expected_mass) * 2E6 / (expected_mass + raw_peak_mass);
                    score.raw_squared_ppm_difference +=
                        peak_intensity * raw_signed_ppm * raw_signed_ppm;
                    let calc_mz = frag.monoisotopic_mass / charge as f32 + PROTON;

                    match frag.kind {
                        Kind::A | Kind::B | Kind::C => {
                            score.matched_b += 1;
                            score.summed_b += peak_intensity;
                            b_run.matched(group.series_index);
                            coverage.forward[group.series_index] += 1;
                        }
                        Kind::X | Kind::Y | Kind::Z | Kind::ZDot => {
                            score.matched_y += 1;
                            score.summed_y += peak_intensity;
                            y_run.matched(group.series_index);
                            coverage.reverse[group.series_index + 1] += 1;
                        }
                    }

                    if collect_fragments {
                        let idx = match frag.kind {
                            Kind::A | Kind::B | Kind::C => group.series_index as i32 + 1,
                            Kind::X | Kind::Y | Kind::Z | Kind::ZDot => {
                                peptide.sequence.len().saturating_sub(1) as i32
                                    - group.series_index as i32
                            }
                        };
                        fragments_details.kinds.push(frag.kind);
                        fragments_details.charges.push(charge as i32);
                        fragments_details.mz_experimental.push(exp_mz);
                        fragments_details.mz_calculated.push(calc_mz);
                        fragments_details.fragment_ordinals.push(idx);
                        fragments_details.intensities.push(peak_intensity);
                        fragments_details
                            .neutral_losses
                            .push(frag.neutral_loss.unwrap_or(0.0));
                        fragments_details.peak_indices.push(peak_idx as u32);
                    }
                }
            }
        }

        score.hyperscore = score.hyperscore(self.score_type);
        score.longest_b = b_run.longest;
        score.longest_y = y_run.longest;
        score.ppm_difference /= score.summed_b + score.summed_y;
        score.signed_ppm_difference /= score.summed_b + score.summed_y;
        score.raw_ppm_difference /= score.summed_b + score.summed_y;
        score.raw_signed_ppm_difference /= score.summed_b + score.summed_y;
        score.raw_squared_ppm_difference /= score.summed_b + score.summed_y;

        if collect_fragments {
            (score, Some(fragments_details), coverage)
        } else {
            // drop(fragments_details);
            (score, None, coverage)
        }
    }
}

/// Remove repeated peptidoforms for the same precursor charge, keeping the
/// first (best-ranked) occurrence. Overlapping offsets, or an offset that
/// reproduces an indexed peptidoform, must not produce artificial ties.
fn dedup_peptidoforms(candidates: &mut Vec<Candidate<'_>>) {
    let mut seen: HashMap<(u8, bool, &[u8]), Vec<usize>> = HashMap::new();
    let mut keep = vec![true; candidates.len()];
    for (index, candidate) in candidates.iter().enumerate() {
        let key = (
            candidate.score.precursor_charge,
            candidate.peptide.decoy,
            candidate.peptide.sequence.as_ref(),
        );
        let previous = seen.entry(key).or_default();
        if previous
            .iter()
            .any(|&kept| same_peptidoform(&candidates[kept].peptide, &candidate.peptide))
        {
            keep[index] = false;
        } else {
            previous.push(index);
        }
    }
    let mut keep = keep.into_iter();
    candidates.retain(|_| keep.next().unwrap_or(true));
}

/// Maintain information about the longest continous ion ladder for a series
#[derive(Default)]
struct Run {
    start: usize,
    length: usize,
    last: usize,
    pub longest: usize,
}

impl Run {
    pub fn matched(&mut self, index: usize) {
        if self.last == index {
            return;
        } else if self.start + self.length == index {
            self.length += 1;
            self.longest = self.longest.max(self.length);
        } else {
            self.start = index;
            self.length = 1;
            self.longest = self.longest.max(self.length);
        }
        self.last = index;
    }
}

#[cfg(test)]
#[path = "../tests/unit/scoring.rs"]
mod tests;
