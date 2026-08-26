use crate::database::{IndexedDatabase, PeptideIx};
use crate::heap::bounded_min_heapify;
use crate::ion_series::{IonGroupSeries, Kind};
use crate::mass::{Tolerance, NEUTRON, PROTON};
use crate::spectrum::{Precursor, ProcessedSpectrum};
use serde::{Deserialize, Serialize};
use std::cell::RefCell;
use std::ops::AddAssign;
use std::sync::atomic::{AtomicBool, AtomicUsize, Ordering};

#[derive(Copy, Clone, Debug, Serialize, Deserialize)]
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
    ppm_difference: f32,
    signed_ppm_difference: f32,
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
        self.hyperscore
            .partial_cmp(&other.hyperscore)
            .unwrap_or(std::cmp::Ordering::Less)
    }
}

/// Preliminary score - # of matched peaks for each candidate peptide
#[derive(Copy, Clone, Default, Debug, PartialEq, Eq, PartialOrd, Ord)]
struct PreScore {
    matched: u16,
    peptide: PeptideIx,
    precursor_charge: u8,
    isotope_error: i8,
}

#[derive(Copy, Clone)]
struct FragmentMatchPeak {
    neutral_mass: f32,
    query_index: usize,
    charge: u8,
    charge_is_known: bool,
}

struct FragmentMatchIndex {
    peaks: Vec<FragmentMatchPeak>,
}

impl FragmentMatchIndex {
    fn new(query: &ProcessedSpectrum, max_charge: u8) -> Self {
        let capacity = query.masses.len() * usize::from(max_charge.saturating_sub(1));
        let mut peaks = Vec::with_capacity(capacity);
        for (query_index, &mass) in query.masses.iter().enumerate() {
            if query.has_known_charge(query_index) {
                let charge = query.charges[query_index];
                if charge < max_charge {
                    peaks.push(FragmentMatchPeak {
                        neutral_mass: mass,
                        query_index,
                        charge,
                        charge_is_known: true,
                    });
                }
            } else {
                for charge in 1..max_charge {
                    peaks.push(FragmentMatchPeak {
                        neutral_mass: mass * charge as f32,
                        query_index,
                        charge,
                        charge_is_known: false,
                    });
                }
            }
        }
        peaks.sort_unstable_by(|left, right| left.neutral_mass.total_cmp(&right.neutral_mass));
        Self { peaks }
    }

    fn select_peak(
        &self,
        query: &ProcessedSpectrum,
        neutral_mass: f32,
        charge: u8,
        tolerance: Tolerance,
    ) -> Option<usize> {
        let scaled_tolerance = match tolerance {
            Tolerance::Da(lo, hi) => Tolerance::Da(lo * charge as f32, hi * charge as f32),
            _ => tolerance,
        };
        let (lo, hi) = scaled_tolerance.bounds(neutral_mass);
        let start = self
            .peaks
            .partition_point(|peak| peak.neutral_mass.total_cmp(&lo).is_lt())
            .saturating_sub(1);
        let mut best_peak = None;
        let mut max_intensity = 0.0;
        for peak in &self.peaks[start..] {
            if peak.neutral_mass.total_cmp(&hi).is_gt() {
                break;
            }
            if peak.charge != charge {
                continue;
            }
            let compatible = if peak.charge_is_known {
                tolerance.contains(neutral_mass, peak.neutral_mass)
            } else {
                scaled_tolerance.contains(neutral_mass, peak.neutral_mass)
            };
            if !compatible {
                continue;
            }
            let intensity = query.intensities[peak.query_index];
            let later_tie =
                intensity == max_intensity && best_peak.is_none_or(|best| peak.query_index >= best);
            if intensity > max_intensity || later_tie {
                max_intensity = intensity;
                best_peak = Some(peak.query_index);
            }
        }
        best_peak
    }
}

/// Store preliminary scores & stats for first pass search for a query spectrum
#[derive(Clone, Default)]
struct InitialHits {
    matched_peaks: usize,
    // Number of peptide candidates with > 0 matched peaks
    scored_candidates: usize,
    preliminary: Vec<PreScore>,
}

impl AddAssign<InitialHits> for InitialHits {
    fn add_assign(&mut self, rhs: InitialHits) {
        self.matched_peaks += rhs.matched_peaks;
        self.scored_candidates += rhs.scored_candidates;

        self.preliminary.extend(rhs.preliminary);
    }
}

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
    /// Precursor mass error after per-file retention-time alignment.
    #[serde(skip_serializing)]
    pub aligned_delta_mass: f32,
    /// Fragment mass error after per-file retention-time alignment.
    #[serde(skip_serializing)]
    pub aligned_average_ppm: f32,
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
    /// Spectral angle for empirical library matching (zero for database search).
    pub spectral_angle: f32,
    /// Fraction of empirical library intensity explained by matched peaks.
    pub explained_library_intensity: f32,
    /// Fraction of query intensity explained by matched library peaks.
    pub explained_query_intensity: f32,
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

    /// Sequence-ambiguity annotation: the peptide string with residues lacking
    /// flanking fragment-ion evidence wrapped in `(?...)`, plus any residual
    /// mass-shift placement.
    pub ambiguity_sequence: String,
    /// Residual precursor mass shift (`expmass - calcmass`) placed during
    /// ambiguity annotation; 0.0 when within the closed-search tolerance.
    pub mass_shift: f32,

    /// Per-modification PTM site localization, if localization is enabled
    #[serde(skip_serializing)]
    pub localization: Option<crate::ptm::Localization>,
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
}

#[inline(always)]
/// Calculate upper bound (excluded) of the charge state range to use for
/// searching fragment ions (1..N)
/// If user has configured max_fragment_charge, potentially override precursor
/// charge
fn max_fragment_charge(max_fragment_charge: Option<u8>, precursor_charge: u8) -> u8 {
    precursor_charge
        .min(
            max_fragment_charge
                .map(|c| c + 1)
                .unwrap_or(precursor_charge),
        )
        .max(2)
}

impl<'db> Scorer<'db> {
    /// Perform a quick first-pass scoring, where we consider a peptide "identified"
    /// if it meets the following criterion:
    ///  * prefilter_low_memory = true: in the top `report_psms` hits for a spectrum
    ///  * prefilter_low_memory = false: has at least `min_matched_peaks` fragment ion matches
    /// * `keep`: A vector of atomic bools is used to maintain an identification list across scans
    pub fn quick_score(
        &self,
        query: &ProcessedSpectrum,
        prefilter_low_memory: bool,
        keep: &[AtomicBool],
    ) {
        assert_eq!(
            query.level, 2,
            "internal bug, trying to score a non-MS2 scan!"
        );
        let precursor = query.precursors.first().unwrap_or_else(|| {
            panic!("missing MS1 precursor for {}", query.id);
        });
        let hits = self.initial_hits(query, precursor);
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

        if prefilter_low_memory {
            let mut score_vector = hits
                .preliminary
                .iter()
                .filter_map(|pre| {
                    if pre.peptide == PeptideIx::default() {
                        return None;
                    }
                    let (score, _, _) =
                        self.score_candidate_with_index(query, pre, false, &fragment_index);
                    if (score.matched_b + score.matched_y) < self.min_matched_peaks {
                        return None;
                    }
                    Some(score)
                })
                .collect::<Vec<_>>();

            let k = self.report_psms.min(score_vector.len());
            bounded_min_heapify(&mut score_vector, k);
            for score in &score_vector[..k] {
                keep[score.peptide.0 as usize].store(true, Ordering::Relaxed);
            }
        } else {
            for pre in &hits.preliminary {
                if pre.peptide != PeptideIx::default() {
                    keep[pre.peptide.0 as usize].store(true, Ordering::Relaxed);
                }
            }
        }
    }

    pub fn score(&self, query: &ProcessedSpectrum) -> Vec<Feature> {
        assert_eq!(
            query.level, 2,
            "internal bug, trying to score a non-MS2 scan!"
        );
        match self.chimera {
            true => self.score_chimera_fast(query),
            false => self.score_standard(query),
        }
    }

    /// Perform a k-select and truncation of an [`InitialHits`] list.
    ///
    /// Determine how many candidates to actually calculate hyperscore for.
    /// Hyperscore is relatively computationally expensive, so we don't want
    /// to calculate it for every possible candidate (100s - 10,000s depending on search)
    /// when we are only going to report a few PSMs. But we also want to calculate
    /// it for enough candidates that we don't accidentally miss the best hit!
    ///
    /// Given that hyperscore is dominated by the number of matched peaks, it seems
    /// reasonable to assume that the highest hyperscore will belong to one of the
    /// top 50 candidates sorted by # of matched peaks.
    fn trim_hits(&self, hits: &mut InitialHits) {
        let k = 50.clamp(
            (self.report_psms * 2).min(hits.preliminary.len()),
            hits.preliminary.len(),
        );
        bounded_min_heapify(&mut hits.preliminary, k);
        hits.preliminary.truncate(k);
    }

    /// Preliminary Score, return # of matched peaks per candidate
    /// Returned hits are guaranteed to be the top-K hits (see above comment)
    /// from among all potential candidates, but the returned vector is not
    /// in sorted order.
    fn matched_peaks_with_isotope(
        &self,
        fragment_index: &FragmentMatchIndex,
        precursor_mass: f32,
        precursor_charge: u8,
        precursor_tol: Tolerance,
        isotope_error: i8,
    ) -> InitialHits {
        let candidates = self.db.query(
            precursor_mass - isotope_error as f32 * NEUTRON,
            precursor_tol,
            self.fragment_tol,
        );

        let potential = candidates.pre_idx_hi - candidates.pre_idx_lo + 1;

        PRE_SCORE_SCRATCH.with(|scratch| {
            let mut preliminary = scratch.borrow_mut();
            preliminary.resize(potential, PreScore::default());
            preliminary.fill(PreScore::default());
            let mut matched_peaks = 0;
            let mut scored_candidates = 0;

            for peak in &fragment_index.peaks {
                for frag in candidates.page_search(peak.neutral_mass) {
                    let idx = frag.peptide_index.0 as usize - candidates.pre_idx_lo;
                    let sc = &mut preliminary[idx];
                    if sc.matched == 0 {
                        scored_candidates += 1;
                        sc.precursor_charge = precursor_charge;
                        sc.peptide = frag.peptide_index;
                        sc.isotope_error = isotope_error;
                    }

                    sc.matched += 1;
                    matched_peaks += 1;
                }
            }

            if matched_peaks == 0 {
                return InitialHits::default();
            }

            let k = 50.clamp(
                (self.report_psms * 2).min(preliminary.len()),
                preliminary.len(),
            );
            bounded_min_heapify(&mut preliminary, k);
            InitialHits {
                matched_peaks,
                scored_candidates,
                preliminary: preliminary[..k].to_vec(),
            }
        })
    }

    fn matched_peaks(
        &self,
        query: &ProcessedSpectrum,
        precursor_mass: f32,
        precursor_charge: u8,
        precursor_tol: Tolerance,
    ) -> InitialHits {
        let max_fragment_charge = max_fragment_charge(self.max_fragment_charge, precursor_charge);
        let fragment_index = FragmentMatchIndex::new(query, max_fragment_charge);
        if self.min_isotope_err != self.max_isotope_err {
            let mut hits = (self.min_isotope_err..=self.max_isotope_err).fold(
                InitialHits::default(),
                |mut hits, isotope| {
                    hits += self.matched_peaks_with_isotope(
                        &fragment_index,
                        precursor_mass,
                        precursor_charge,
                        precursor_tol,
                        isotope,
                    );
                    hits
                },
            );
            self.trim_hits(&mut hits);
            hits
        } else {
            self.matched_peaks_with_isotope(
                &fragment_index,
                precursor_mass,
                precursor_charge,
                precursor_tol,
                self.min_isotope_err,
            )
        }
    }

    fn initial_hits(&self, query: &ProcessedSpectrum, precursor: &Precursor) -> InitialHits {
        // Sage operates on masses without protons; [M] instead of [MH+]
        let mz = precursor.mz - PROTON;

        // Search in wide-window/DIA mode
        let mut hits = if self.wide_window {
            (self.min_precursor_charge..=self.max_precursor_charge).fold(
                InitialHits::default(),
                |mut hits, precursor_charge| {
                    let precursor_mass = mz * precursor_charge as f32;
                    let precursor_tol = precursor
                        .isolation_window
                        .unwrap_or(Tolerance::Da(-2.4, 2.4))
                        * precursor_charge as f32;
                    hits +=
                        self.matched_peaks(query, precursor_mass, precursor_charge, precursor_tol);
                    hits
                },
            )
        } else if precursor.charge.is_some() && !self.override_precursor_charge {
            let charge = precursor.charge.unwrap();
            // Charge state is already annotated for this precusor, only search once
            let precursor_mass = mz * charge as f32;
            self.matched_peaks(query, precursor_mass, charge, self.precursor_tol)
        } else {
            // Not all selected ion precursors have charge states annotated (or user has set
            // `override_precursor_charge`)
            // assume it could be z=2, z=3, z=4 and search all three
            (self.min_precursor_charge..=self.max_precursor_charge).fold(
                InitialHits::default(),
                |mut hits, precursor_charge| {
                    let precursor_mass = mz * precursor_charge as f32;
                    hits += self.matched_peaks(
                        query,
                        precursor_mass,
                        precursor_charge,
                        self.precursor_tol,
                    );
                    hits
                },
            )
        };
        self.trim_hits(&mut hits);
        hits
    }

    /// Score a single [`ProcessedSpectrum`] against the database
    pub fn score_standard(&self, query: &ProcessedSpectrum) -> Vec<Feature> {
        let precursor = query.precursors.first().unwrap_or_else(|| {
            panic!("missing MS1 precursor for {}", query.id);
        });

        let hits = self.initial_hits(query, precursor);
        let mut features = Vec::with_capacity(self.report_psms);
        self.build_features(query, precursor, &hits, self.report_psms, &mut features);
        features
    }

    /// Given a set of [`InitialHits`] against a query spectrum, prepare N=`report_psms`
    /// best PSMs ([`Feature`])
    fn build_features(
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
        let mut score_vector = hits
            .preliminary
            .iter()
            .filter(|score| score.peptide != PeptideIx::default())
            .map(|pre| {
                self.score_candidate_with_index(query, pre, self.annotate_matches, &fragment_index)
            })
            .filter(|s| (s.0.matched_b + s.0.matched_y) >= self.min_matched_peaks)
            .collect::<Vec<_>>();

        // Hyperscore is our primary score function for PSMs
        score_vector.sort_by(|a, b| b.0.hyperscore.total_cmp(&a.0.hyperscore));

        // Expected value for poisson distribution
        // (average # of matches peaks/peptide candidate)
        let lambda = hits.matched_peaks as f64 / hits.scored_candidates as f64;

        // Sage operates on masses without protons; [M] instead of [MH+]
        let mz = precursor.mz - PROTON;

        for idx in 0..report_psms.min(score_vector.len()) {
            let score = score_vector[idx].0;
            let fragments: Option<Fragments> = score_vector[idx].1.take();
            let coverage = std::mem::take(&mut score_vector[idx].2);
            let psm_id = increment_psm_counter();

            let peptide = &self.db[score.peptide];
            let precursor_mass = mz * score.precursor_charge as f32;

            let next = score_vector
                .get(idx + 1)
                .map(|score| score.0.hyperscore)
                .unwrap_or_default();

            let best = score_vector
                .first()
                .map(|score| score.0.hyperscore)
                .expect("we know that index 0 is valid");

            // Poisson distribution log10 probability mass function
            // Computed directly in log space to avoid overflow from lambda.powi(k)
            // log10(PMF) = (k*ln(lambda) - lambda - lnfact(k)) / ln(10)
            let k = score.matched_b + score.matched_y;
            let log10_poisson =
                (k as f64 * lambda.ln() - lambda - lnfact(k)) / std::f64::consts::LN_10;

            let isotope_error = score.isotope_error as f32 * NEUTRON;
            let delta_mass = (precursor_mass - peptide.monoisotopic - isotope_error) * 2E6
                / (precursor_mass - isotope_error + peptide.monoisotopic);

            // Sequence-ambiguity annotation. A residual precursor mass shift is
            // only placed when it exceeds the closed-search tolerance (a small
            // fixed ppm threshold, independent of the precursor search window so
            // that wide/open searches still surface real shifts).
            let raw_mass_shift = precursor_mass - peptide.monoisotopic;
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
                expmass: precursor_mass,
                calcmass: peptide.monoisotopic,
                // Features
                charge: score.precursor_charge,
                rt: query.scan_start_time,
                ims: query
                    .precursors
                    .first()
                    .unwrap()
                    .inverse_ion_mobility
                    .unwrap_or(0.0),
                delta_mass,
                aligned_delta_mass: delta_mass,
                isotope_error,
                average_ppm: score.ppm_difference,
                signed_fragment_ppm: score.signed_ppm_difference,
                aligned_average_ppm: score.ppm_difference,
                hyperscore: score.hyperscore,
                delta_next: score.hyperscore - next,
                delta_best: best - score.hyperscore,
                matched_peaks: k as u32,
                matched_intensity_pct: 100.0 * (score.summed_b + score.summed_y)
                    / query.total_ion_current,
                spectral_angle: 0.0,
                explained_library_intensity: 0.0,
                explained_query_intensity: 0.0,
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

                //Fragments
                protein_groups: None,
                num_protein_groups: 0,
                fragments,
                ambiguity_sequence: ambiguity.sequence,
                mass_shift: ambiguity.mass_shift,
                protein_group_q: 1.0,
                localization: None,
            })
        }
    }

    /// Remove peaks matching a PSM from a query spectrum
    fn remove_matched_peaks(&self, query: &mut ProcessedSpectrum, psm: &Feature) {
        let peptide = &self.db[psm.peptide_idx];
        let fragments = self
            .db
            .ion_kinds
            .iter()
            .flat_map(|kind| IonGroupSeries::new(peptide, *kind))
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

    /// Return multiple PSMs for each spectra - first is the best match, second PSM is the best match
    /// after all theoretical peaks assigned to the best match are removed, etc
    pub fn score_chimera_fast(&self, query: &ProcessedSpectrum) -> Vec<Feature> {
        let precursor = query.precursors.first().unwrap_or_else(|| {
            panic!("missing MS1 precursor for {}", query.id);
        });

        let mut query = query.clone();
        let hits = self.initial_hits(&query, precursor);

        let mut candidates: Vec<Feature> = Vec::with_capacity(self.report_psms);

        let mut prev = 0;
        while candidates.len() < self.report_psms {
            self.build_features(&query, precursor, &hits, 1, &mut candidates);
            if candidates.len() > prev {
                if let Some(feat) = candidates.get_mut(prev) {
                    self.remove_matched_peaks(&mut query, feat);
                    feat.rank = prev as u32 + 1;
                }
                prev = candidates.len()
            } else {
                break;
            }
        }
        candidates
    }

    /// Calculate full hyperscore for a given PSM
    fn score_candidate(
        &self,
        query: &ProcessedSpectrum,
        pre_score: &PreScore,
        collect_fragments: bool,
    ) -> (Score, Option<Fragments>, Coverage) {
        let max_fragment_charge =
            max_fragment_charge(self.max_fragment_charge, pre_score.precursor_charge);
        let fragment_index = FragmentMatchIndex::new(query, max_fragment_charge);
        self.score_candidate_with_index(query, pre_score, collect_fragments, &fragment_index)
    }

    fn score_candidate_with_index(
        &self,
        query: &ProcessedSpectrum,
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
        let peptide = &self.db[score.peptide];
        let max_fragment_charge =
            max_fragment_charge(self.max_fragment_charge, score.precursor_charge);

        // Regenerate theoretical ions - initial database search might be
        // using only a subset of all possible ions (e.g. no b1/b2/y1/y2)
        // so we need to completely re-score this candidate
        let fragment_groups = self
            .db
            .ion_kinds
            .iter()
            .flat_map(|kind| IonGroupSeries::new(peptide, *kind));

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
                let best = group
                    .variants
                    .iter()
                    .filter_map(|variant| {
                        fragment_index
                            .select_peak(
                                query,
                                variant.monoisotopic_mass,
                                charge,
                                self.fragment_tol,
                            )
                            .map(|peak_idx| (variant, peak_idx))
                    })
                    .max_by(|a, b| query.intensities[a.1].total_cmp(&query.intensities[b.1]));

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

                    let exp_mz = query.peak_mz(peak_idx);
                    let calc_mz = frag.monoisotopic_mass / charge as f32 + PROTON;

                    match frag.kind {
                        Kind::A | Kind::B | Kind::C => {
                            score.matched_b += 1;
                            score.summed_b += peak_intensity;
                            b_run.matched(group.series_index);
                            coverage.forward[group.series_index] += 1;
                        }
                        Kind::X | Kind::Y | Kind::Z => {
                            score.matched_y += 1;
                            score.summed_y += peak_intensity;
                            y_run.matched(group.series_index);
                            coverage.reverse[group.series_index + 1] += 1;
                        }
                    }

                    if collect_fragments {
                        let idx = match frag.kind {
                            Kind::A | Kind::B | Kind::C => group.series_index as i32 + 1,
                            Kind::X | Kind::Y | Kind::Z => {
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
                    }
                }
            }
        }

        score.hyperscore = score.hyperscore(self.score_type);
        score.longest_b = b_run.longest;
        score.longest_y = y_run.longest;
        score.ppm_difference /= score.summed_b + score.summed_y;
        score.signed_ppm_difference /= score.summed_b + score.summed_y;

        if collect_fragments {
            (score, Some(fragments_details), coverage)
        } else {
            // drop(fragments_details);
            (score, None, coverage)
        }
    }

    /// Reconstruct detailed fragment matches for an already-scored PSM.
    /// This intentionally reuses the same candidate-matching path as the
    /// primary search so neutral-loss selection, charge assignment, and mass
    /// calculations cannot drift between scoring and deferred annotation.
    pub fn annotate_candidate(&self, query: &ProcessedSpectrum, feature: &Feature) -> Fragments {
        let pre_score = PreScore {
            peptide: feature.peptide_idx,
            precursor_charge: feature.charge,
            ..Default::default()
        };
        self.score_candidate(query, &pre_score, true)
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
        if !self.chimera {
            return features
                .iter()
                .zip(selected)
                .map(|(feature, selected)| {
                    selected.then(|| self.annotate_candidate(query, feature))
                })
                .collect();
        }

        let mut residual = query.clone();
        features
            .iter()
            .zip(selected)
            .map(|(feature, selected)| {
                let fragments = selected.then(|| self.annotate_candidate(&residual, feature));
                self.remove_matched_peaks(&mut residual, feature);
                fragments
            })
            .collect()
    }
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
