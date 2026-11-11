//! The reported result types and the code that fills them.
//!
//! Moved out of `scoring.rs` unchanged; every item keeps its `scoring::` path.

use crate::database::{MassOffsetAssignment, PeptideIx};
use crate::ion_series::{IonVariant, Kind};
use crate::mass::{NEUTRON, PROTON};
use crate::peptide::Peptide;
use crate::spectrum::{Precursor, ProcessedSpectrum};
use serde::Serialize;

use super::*;

#[derive(Serialize, Clone, Debug, Default)]
/// Features of a candidate peptide spectrum match
pub struct Feature {
    #[serde(skip_serializing)]
    pub peptide_idx: PeptideIx,
    // psm_id help to match with matched fragments table.
    pub psm_id: usize,
    pub peptide_len: usize,
    /// Spectrum id
    pub spec_id: String,
    /// File identifier
    pub file_id: usize,
    /// PSM rank
    pub rank: u32,
    /// Target/Decoy label, -1 is decoy, 1 is target
    pub label: i32,
    /// Experimental mass
    pub expmass: f32,
    /// Calculated mass
    pub calcmass: f32,
    /// Reported precursor charge
    pub charge: u8,
    /// Retention time
    pub rt: f32,
    /// Globally aligned retention time
    pub aligned_rt: f32,
    /// Predicted RT, if enabled
    pub predicted_rt: f32,
    /// Difference between predicted & observed RT
    pub delta_rt_model: f32,
    /// Ion mobility
    pub ims: f32,
    /// Predicted ion mobility, if enabled
    pub predicted_ims: f32,
    /// Difference between predicted & observed ion mobility
    pub delta_ims_model: f32,
    /// Difference between expmass and calcmass
    pub delta_mass: f32,
    /// C13 isotope error
    pub isotope_error: f32,
    /// Average ppm delta mass for matched fragments
    pub average_ppm: f32,
    /// Signed, intensity-weighted fragment error (observed - theoretical).
    #[serde(skip_serializing)]
    pub signed_fragment_ppm: f32,
    /// Intensity-weighted standard deviation of the raw signed fragment
    /// errors within this PSM, around `signed_fragment_ppm`.
    #[serde(skip_serializing)]
    pub fragment_ppm_sd: f32,
    /// Precursor mass error after per-file retention-time alignment.
    #[serde(skip_serializing)]
    pub aligned_delta_mass: f32,
    /// Fragment mass error after per-file retention-time alignment.
    #[serde(skip_serializing)]
    pub aligned_average_ppm: f32,
    /// Signed fragment error after search-time recalibration; equal to
    /// `signed_fragment_ppm` when no fragment correction applies.
    #[serde(skip_serializing)]
    pub aligned_signed_fragment_ppm: f32,
    /// X!Tandem hyperscore
    pub hyperscore: f64,
    /// Difference between hyperscore of this candidate, and the next best candidate
    pub delta_next: f64,
    /// Difference between hyperscore of this candidate, and the best candidate
    pub delta_best: f64,
    /// Number of matched theoretical fragment ions
    pub matched_peaks: u32,
    /// Longest b-ion series
    pub longest_b: u32,
    /// Longest y-ion series
    pub longest_y: u32,
    /// Longest y-ion series, divided by peptide length
    pub longest_y_pct: f32,
    /// Number of missed cleavages
    pub missed_cleavages: u8,
    /// Fraction of matched MS2 intensity
    pub matched_intensity_pct: f32,
    /// Number of scored candidates for this spectrum
    pub scored_candidates: u32,
    /// Probability of matching exactly N peaks across all candidates Pr(x=k)
    pub poisson: f64,
    /// Combined score from linear discriminant analysis, used for FDR calc
    pub discriminant_score: f32,
    /// Posterior error probability for this PSM / local FDR
    pub posterior_error: f32,
    /// Assigned q_value
    pub spectrum_q: f32,
    pub peptide_q: f32,
    pub protein_q: f32,
    pub protein_group_q: f32,

    pub ms2_intensity: f32,

    pub protein_groups: Option<String>,
    pub num_protein_groups: u32,

    pub fragments: Option<Fragments>,

    /// Generic fragment-loss evidence, when `database.fragment_losses` is
    /// configured.
    #[serde(skip_serializing)]
    pub fragment_loss: Option<FragmentLossFeatures>,

    /// Sequence-ambiguity annotation: the peptide string with residues lacking
    /// flanking fragment-ion evidence wrapped in `(?...)`, plus any residual
    /// mass-shift placement.
    pub ambiguity_sequence: String,
    /// Residual precursor mass shift (`expmass - calcmass - isotope_error`)
    /// placed during ambiguity annotation; 0.0 when within the closed-search
    /// tolerance.
    pub mass_shift: f32,

    /// Per-modification PTM site localization, if localization is enabled
    #[serde(skip_serializing)]
    pub localization: Option<crate::ptm::Localization>,

    /// Search-time mass offset placed on `peptide_idx`. Before
    /// [`IndexedDatabase::materialize_mass_offsets`] runs, `peptide_idx` is
    /// the unmodified indexed peptide; afterwards it addresses the placed
    /// peptidoform.
    #[serde(skip_serializing)]
    pub mass_offset: Option<MassOffsetAssignment>,

    /// Immonium-ion evidence; `None` unless `immonium` is enabled.
    #[serde(skip_serializing)]
    pub immonium: Option<crate::immonium::ImmoniumEvidence>,
}

/// Matched generic fragment-loss ions (`database.fragment_losses`) of a PSM.
/// They do not enter the hyperscore or `matched_peaks`; they are separate
/// rescoring features.
#[derive(Serialize, Clone, Copy, Debug, Default, PartialEq)]
pub struct FragmentLossFeatures {
    /// Cleavage and charge pairs whose loss form matched a peak.
    pub matched_peaks: u32,
    /// Percent of MS2 intensity matched by loss ions.
    pub intensity_pct: f32,
}

/// Matching Fragment details
#[derive(Serialize, Default, Clone, Debug, PartialEq)]
pub struct Fragments {
    /// Observed fragment charge state.
    #[serde(skip_serializing)]
    pub charges: Vec<i32>,
    pub kinds: Vec<Kind>,
    pub fragment_ordinals: Vec<i32>,
    pub intensities: Vec<f32>,
    pub mz_calculated: Vec<f32>,
    pub mz_experimental: Vec<f32>,
    /// Neutral loss applied to the matched theoretical fragment; zero is the
    /// retained (no-loss) variant.
    pub neutral_losses: Vec<f32>,
    /// Index of the matched peak in the scored spectrum, so callers combining
    /// several matches (for example two crosslinked chains) can count shared
    /// peaks once.
    #[serde(skip_serializing)]
    pub peak_indices: Vec<u32>,
}

/// Per-residue fragment-ion coverage, used to compute sequence-ambiguity
/// annotation. `forward[i]` counts matched a/b/c ions mapping to residue `i`
/// (ion ordinal `i + 1`); `reverse[i]` counts matched x/y/z ions mapping to
/// residue `i` (ordinal `n - i`).
#[derive(Default, Clone, Debug)]
pub struct Coverage {
    pub forward: Vec<u16>,
    pub reverse: Vec<u16>,
}

impl<'db> Scorer<'db> {
    /// Given a set of [`InitialHits`] against a query spectrum, prepare N=`report_psms`
    /// best PSMs ([`Feature`])
    pub(super) fn build_features(
        &self,
        query: &ProcessedSpectrum,
        precursor: &Precursor,
        hits: &InitialHits,
        report_psms: usize,
        features: &mut Vec<Feature>,
    ) {
        let max_charge = hits
            .preliminary
            .iter()
            .map(|pre| pre.precursor_charge)
            .max()
            .unwrap_or(self.min_precursor_charge);
        let fragment_index = FragmentMatchIndex::new(
            query,
            max_fragment_charge(self.max_fragment_charge, max_charge),
        );
        let mut score_vector = Vec::with_capacity(hits.preliminary.len());
        for pre in hits
            .preliminary
            .iter()
            .filter(|score| score.peptide != PeptideIx::default())
        {
            self.score_hypothesis(
                query,
                pre,
                self.annotate_matches,
                &fragment_index,
                &mut score_vector,
            );
        }
        score_vector.retain(|s| (s.score.matched_b + s.score.matched_y) >= self.min_matched_peaks);

        // Hyperscore is primary. Among exact ties the smaller precursor isotope
        // error wins: an isotope step must not turn one peptidoform into a
        // lighter one that scores identically (for example a hydrolyzed DSSO
        // monolink, 0.984 Da heavier than the amidated one, read as amidated
        // at isotope +1). Peptidoform identity, charge, and isotope then make
        // exact ties deterministic across database layouts and chunk filtering.
        score_vector.sort_unstable_by(|a, b| {
            b.score
                .hyperscore
                .total_cmp(&a.score.hyperscore)
                .then_with(|| {
                    a.score
                        .isotope_error
                        .unsigned_abs()
                        .cmp(&b.score.isotope_error.unsigned_abs())
                })
                .then_with(|| {
                    a.peptide
                        .monoisotopic
                        .total_cmp(&b.peptide.monoisotopic)
                        .then_with(|| a.peptide.initial_sort(&b.peptide))
                })
                .then_with(|| a.score.precursor_charge.cmp(&b.score.precursor_charge))
                .then_with(|| a.score.isotope_error.cmp(&b.score.isotope_error))
                .then_with(|| a.mass_offset.is_some().cmp(&b.mass_offset.is_some()))
                .then_with(|| a.score.peptide.cmp(&b.score.peptide))
        });
        if !self.db.mass_offsets.is_empty() {
            dedup_peptidoforms(&mut score_vector);
        }

        // Expected value for poisson distribution
        // (average # of matches peaks/peptide candidate)
        let lambda = hits.matched_peaks as f64 / hits.scored_candidates as f64;

        // Sage operates on masses without protons; [M] instead of [MH+]
        let mz = precursor.mz - PROTON;
        // `precursor` may be recalibrated; errors are also reported against
        // the observed precursor m/z.
        let raw_mz = match self.correction(query.file_id) {
            Some(correction) => {
                let ppm = correction.precursor_ppm(query.scan_start_time, precursor.mz);
                precursor.mz * (1.0 + ppm * 1e-6) - PROTON
            }
            None => mz,
        };

        for idx in 0..report_psms.min(score_vector.len()) {
            let score = score_vector[idx].score;
            let fragments: Option<Fragments> = score_vector[idx].fragments.take();
            let coverage = std::mem::take(&mut score_vector[idx].coverage);
            let mass_offset = score_vector[idx].mass_offset;
            let psm_id = increment_psm_counter();

            let peptide: &Peptide = &score_vector[idx].peptide;
            let precursor_mass = mz * score.precursor_charge as f32;
            let raw_precursor_mass = raw_mz * score.precursor_charge as f32;

            let next = score_vector
                .get(idx + 1)
                .map(|candidate| candidate.score.hyperscore)
                .unwrap_or_default();

            let best = score_vector
                .first()
                .map(|candidate| candidate.score.hyperscore)
                .expect("we know that index 0 is valid");

            // Poisson distribution log10 probability mass function
            // Computed directly in log space to avoid overflow from lambda.powi(k)
            // log10(PMF) = (k*ln(lambda) - lambda - lnfact(k)) / ln(10)
            let k = score.matched_b + score.matched_y;
            let log10_poisson =
                (k as f64 * lambda.ln() - lambda - lnfact(k)) / std::f64::consts::LN_10;

            let isotope_error = score.isotope_error as f32 * NEUTRON;
            let ppm_delta = |mass: f32| {
                (mass - peptide.monoisotopic - isotope_error) * 2E6
                    / (mass - isotope_error + peptide.monoisotopic)
            };
            let delta_mass = ppm_delta(raw_precursor_mass);
            let aligned_delta_mass = ppm_delta(precursor_mass);

            // Sequence-ambiguity annotation. A residual precursor mass shift is
            // only placed when it exceeds the closed-search tolerance (a small
            // fixed ppm threshold, independent of the precursor search window so
            // that wide/open searches still surface real shifts). The matched
            // isotope error explains part of the precursor delta and is not a
            // residual modification mass.
            let raw_mass_shift = precursor_mass - isotope_error - peptide.monoisotopic;
            let mass_shift =
                if (raw_mass_shift / peptide.monoisotopic * 1e6).abs() <= self.mass_shift_ppm {
                    None
                } else {
                    Some(raw_mass_shift)
                };
            let ambiguity = crate::ambiguity::annotate(
                peptide,
                &coverage.forward,
                &coverage.reverse,
                mass_shift,
            );

            // let (num_proteins, proteins) = self.db.assign_proteins(peptide);

            features.push(Feature {
                // Identifiers
                psm_id,
                peptide_idx: score.peptide,
                spec_id: query.id.clone(),
                file_id: query.file_id,
                rank: idx as u32 + 1,
                label: peptide.label(),
                expmass: raw_precursor_mass,
                calcmass: peptide.monoisotopic,
                // Features
                charge: score.precursor_charge,
                rt: query.scan_start_time,
                ims: precursor.inverse_ion_mobility.unwrap_or(0.0),
                delta_mass,
                aligned_delta_mass,
                isotope_error,
                average_ppm: score.raw_ppm_difference,
                signed_fragment_ppm: score.raw_signed_ppm_difference,
                fragment_ppm_sd: (score.raw_squared_ppm_difference
                    - score.raw_signed_ppm_difference.powi(2))
                .max(0.0)
                .sqrt(),
                aligned_average_ppm: score.ppm_difference,
                aligned_signed_fragment_ppm: score.signed_ppm_difference,
                hyperscore: score.hyperscore,
                delta_next: score.hyperscore - next,
                delta_best: best - score.hyperscore,
                matched_peaks: k as u32,
                matched_intensity_pct: 100.0 * (score.summed_b + score.summed_y)
                    / query.total_ion_current,
                poisson: if log10_poisson.is_finite() {
                    log10_poisson
                } else {
                    f64::NEG_INFINITY
                },
                longest_b: score.longest_b as u32,
                longest_y: score.longest_y as u32,
                longest_y_pct: score.longest_y as f32 / (peptide.sequence.len() as f32),
                peptide_len: peptide.sequence.len(),
                scored_candidates: hits.scored_candidates as u32,
                missed_cleavages: peptide.missed_cleavages,

                // Outputs
                discriminant_score: 0.0,
                posterior_error: 1.0,
                spectrum_q: 1.0,
                protein_q: 1.0,
                peptide_q: 1.0,
                predicted_rt: 0.0,
                predicted_ims: 0.0,
                aligned_rt: query.scan_start_time,
                delta_rt_model: 0.999,
                delta_ims_model: 0.999,
                ms2_intensity: score.summed_b + score.summed_y,
                fragment_loss: self.loss_features().then(|| FragmentLossFeatures {
                    matched_peaks: score.matched_loss as u32,
                    intensity_pct: 100.0 * score.summed_loss / query.total_ion_current,
                }),

                //Fragments
                protein_groups: None,
                num_protein_groups: 0,
                fragments,
                ambiguity_sequence: ambiguity.sequence,
                mass_shift: ambiguity.mass_shift,
                protein_group_q: 1.0,
                localization: None,
                immonium: None,
                mass_offset,
            })
        }
    }

    /// The most intense peak matched by any of `variants` at `charge`.
    pub(super) fn select_variant<'v>(
        &self,
        query: &ProcessedSpectrum,
        fragment_index: &FragmentMatchIndex,
        variants: impl Iterator<Item = &'v IonVariant>,
        charge: u8,
    ) -> Option<(&'v IonVariant, usize)> {
        variants
            .filter_map(|variant| {
                fragment_index
                    .select_peak(query, variant.monoisotopic_mass, charge, self.fragment_tol)
                    .map(|peak_idx| (variant, peak_idx))
            })
            .max_by(|a, b| query.intensities[a.1].total_cmp(&query.intensities[b.1]))
    }
    /// Reconstruct detailed fragment matches for an already-scored PSM.
    /// This intentionally reuses the same candidate-matching path as the
    /// primary search so neutral-loss selection, charge assignment, and mass
    /// calculations cannot drift between scoring and deferred annotation.
    pub fn annotate_candidate(&self, query: &ProcessedSpectrum, feature: &Feature) -> Fragments {
        self.annotate_recalibrated(self.recalibrated(query).as_ref(), feature)
    }

    /// [`Self::annotate_candidate`] for a spectrum already passed through
    /// [`Self::recalibrated`].
    fn annotate_recalibrated(&self, query: &ProcessedSpectrum, feature: &Feature) -> Fragments {
        let pre_score = PreScore {
            peptide: feature.peptide_idx,
            precursor_charge: feature.charge,
            ..Default::default()
        };
        let fragment_index = FragmentMatchIndex::new(
            query,
            max_fragment_charge(self.max_fragment_charge, feature.charge),
        );
        let peptide = self.db.resolve_peptide(feature);
        self.score_peptide(query, &peptide, &pre_score, true, &fragment_index)
            .1
            .expect("fragment collection was explicitly requested")
    }

    /// Reconstruct fragment details for the reported PSMs from one spectrum.
    /// In chimera mode every preceding rank is replayed, even when it is not
    /// selected for output, because later ranks were scored after its peaks
    /// had been removed.
    pub fn annotate_ranked_candidates(
        &self,
        query: &ProcessedSpectrum,
        features: &[&Feature],
        selected: &[bool],
    ) -> Vec<Option<Fragments>> {
        assert_eq!(features.len(), selected.len());
        let query = self.recalibrated(query);
        let query = query.as_ref();
        if !self.chimera {
            return features
                .iter()
                .zip(selected)
                .map(|(feature, selected)| {
                    selected.then(|| self.annotate_recalibrated(query, feature))
                })
                .collect();
        }

        let mut residual = query.clone();
        features
            .iter()
            .zip(selected)
            .map(|(feature, selected)| {
                let fragments = selected.then(|| self.annotate_recalibrated(&residual, feature));
                self.remove_matched_peaks(&mut residual, feature);
                fragments
            })
            .collect()
    }
}
