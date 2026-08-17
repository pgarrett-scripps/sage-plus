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
#[derive(Serialize, Default, Clone, Debug)]
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

    /// Use the bitmap-based preliminary search instead of the bucketed binary search.
    /// Requires deisotoped spectra (each peak must carry a resolved neutral monoisotopic mass).
    pub use_bitmap: bool,
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

        if prefilter_low_memory {
            let mut score_vector = hits
                .preliminary
                .iter()
                .filter_map(|pre| {
                    if pre.peptide == PeptideIx::default() {
                        return None;
                    }
                    let (score, _, _) = self.score_candidate(query, pre, false);
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
        query: &ProcessedSpectrum,
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

        let max_fragment_charge = max_fragment_charge(self.max_fragment_charge, precursor_charge);
        let potential = candidates.pre_idx_hi - candidates.pre_idx_lo + 1;

        PRE_SCORE_SCRATCH.with(|scratch| {
            let mut preliminary = scratch.borrow_mut();
            preliminary.resize(potential, PreScore::default());
            preliminary.fill(PreScore::default());
            let mut matched_peaks = 0;
            let mut scored_candidates = 0;

            for peak_mass in query.masses.iter() {
                for charge in 1..max_fragment_charge {
                    let mass = peak_mass * charge as f32;
                    for frag in candidates.page_search(mass) {
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
        if self.min_isotope_err != self.max_isotope_err {
            let mut hits = (self.min_isotope_err..=self.max_isotope_err).fold(
                InitialHits::default(),
                |mut hits, isotope| {
                    hits += self.matched_peaks_with_isotope(
                        query,
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
                query,
                precursor_mass,
                precursor_charge,
                precursor_tol,
                self.min_isotope_err,
            )
        }
    }

    // ── Bitmap-based preliminary search ──────────────────────────────────────

    /// Inner bitmap scoring for a single (precursor_mass, precursor_charge, isotope_error) triple.
    fn matched_peaks_bitmap_with_isotope(
        &self,
        exp_bitmap: &[u64],
        precursor_mass: f32,
        precursor_charge: u8,
        precursor_tol: Tolerance,
        isotope_error: i8,
    ) -> InitialHits {
        use crate::mass::NEUTRON;
        let search_mass = precursor_mass - isotope_error as f32 * NEUTRON;
        let (lo, hi) = precursor_tol.bounds(search_mass);

        let bm = &self.db.bitmap_index;
        let start = bm.precursor_masses.partition_point(|&m| m < lo);
        let end = bm.precursor_masses.partition_point(|&m| m <= hi);

        if start >= end {
            return InitialHits::default();
        }

        let potential = end - start;

        PRE_SCORE_SCRATCH.with(|scratch| {
            let mut preliminary = scratch.borrow_mut();
            preliminary.resize(potential, PreScore::default());
            preliminary.fill(PreScore::default());
            let mut matched_peaks = 0;
            let mut scored_candidates = 0;

            for i in start..end {
                let (fwd, rev) = bm.score_peptide(exp_bitmap, i);
                let total = fwd + rev;
                if total > 0 {
                    let idx = i - start;
                    let sc = &mut preliminary[idx];
                    sc.matched = total;
                    sc.peptide = bm.peptide_indices[i];
                    sc.precursor_charge = precursor_charge;
                    sc.isotope_error = isotope_error;
                    scored_candidates += 1;
                    matched_peaks += total as usize;
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

    /// Bitmap preliminary search for one (precursor_mass, precursor_charge),
    /// iterating over isotope errors when configured.
    fn matched_peaks_bitmap(
        &self,
        exp_bitmap: &[u64],
        precursor_mass: f32,
        precursor_charge: u8,
        precursor_tol: Tolerance,
    ) -> InitialHits {
        if self.min_isotope_err != self.max_isotope_err {
            let mut hits = (self.min_isotope_err..=self.max_isotope_err).fold(
                InitialHits::default(),
                |mut hits, isotope| {
                    hits += self.matched_peaks_bitmap_with_isotope(
                        exp_bitmap,
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
            self.matched_peaks_bitmap_with_isotope(
                exp_bitmap,
                precursor_mass,
                precursor_charge,
                precursor_tol,
                self.min_isotope_err,
            )
        }
    }

    /// Bitmap-based `initial_hits` — mirrors the structure of the bucketed path
    /// but uses `BitmapIndex` for scoring.
    ///
    /// **Requires deisotoped peaks**: each `peak.mass` must be the neutral
    /// monoisotopic mass M = (mz − H) × z as resolved by deisotoping.
    fn initial_hits_bitmap(&self, query: &ProcessedSpectrum, precursor: &Precursor) -> InitialHits {
        assert!(
            self.db.peptides.is_empty() || !self.db.bitmap_index.precursor_masses.is_empty(),
            "bitmap scoring requested for a database built without the bitmap index; set Parameters::use_bitmap before building"
        );
        let exp_bitmap = self
            .db
            .bitmap_index
            .experimental_bitmap(&query.masses, self.fragment_tol);

        let mz = precursor.mz - PROTON;

        let mut hits = if self.wide_window {
            (self.min_precursor_charge..=self.max_precursor_charge).fold(
                InitialHits::default(),
                |mut hits, precursor_charge| {
                    let precursor_mass = mz * precursor_charge as f32;
                    let precursor_tol = precursor
                        .isolation_window
                        .unwrap_or(Tolerance::Da(-2.4, 2.4))
                        * precursor_charge as f32;
                    hits += self.matched_peaks_bitmap(
                        &exp_bitmap,
                        precursor_mass,
                        precursor_charge,
                        precursor_tol,
                    );
                    hits
                },
            )
        } else if precursor.charge.is_some() && !self.override_precursor_charge {
            let charge = precursor.charge.unwrap();
            let precursor_mass = mz * charge as f32;
            self.matched_peaks_bitmap(&exp_bitmap, precursor_mass, charge, self.precursor_tol)
        } else {
            (self.min_precursor_charge..=self.max_precursor_charge).fold(
                InitialHits::default(),
                |mut hits, precursor_charge| {
                    let precursor_mass = mz * precursor_charge as f32;
                    hits += self.matched_peaks_bitmap(
                        &exp_bitmap,
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

    // ── Bucketed binary-search preliminary search (existing) ─────────────────

    fn initial_hits(&self, query: &ProcessedSpectrum, precursor: &Precursor) -> InitialHits {
        if self.use_bitmap {
            return self.initial_hits_bitmap(query, precursor);
        }

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
        let mut score_vector = hits
            .preliminary
            .iter()
            .filter(|score| score.peptide != PeptideIx::default())
            .map(|pre| self.score_candidate(query, pre, self.annotate_matches))
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

        // Remove MS2 peaks matched by previous match
        let mut to_remove = vec![false; query.masses.len()];
        for frag in fragments {
            for charge in 1..max_fragment_charge {
                // Experimental peaks are multipled by charge, therefore theoretical are divided
                if let Some(peak_idx) = crate::spectrum::select_most_intense_peak(
                    &query.masses,
                    &query.intensities,
                    frag.monoisotopic_mass / charge as f32,
                    self.fragment_tol,
                    None,
                ) {
                    to_remove[peak_idx] = true;
                }
            }
        }

        let mut masses = Vec::with_capacity(query.masses.len());
        let mut intensities = Vec::with_capacity(query.intensities.len());
        let mut charges = Vec::with_capacity(query.charges.len());
        let mut mobilities = Vec::with_capacity(query.mobilities.len());

        for (idx, removed) in to_remove.iter().enumerate() {
            if !removed {
                masses.push(query.masses[idx]);
                intensities.push(query.intensities[idx]);
                charges.push(query.charges[idx]);
                if !query.mobilities.is_empty() {
                    mobilities.push(query.mobilities[idx]);
                }
            }
        }

        query.masses = masses;
        query.intensities = intensities;
        query.charges = charges;
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
                        let mz = variant.monoisotopic_mass / charge as f32;
                        crate::spectrum::select_most_intense_peak(
                            &query.masses,
                            &query.intensities,
                            mz,
                            self.fragment_tol,
                            None,
                        )
                        .map(|peak_idx| (variant, mz, peak_idx))
                    })
                    .max_by(|a, b| query.intensities[a.2].total_cmp(&query.intensities[b.2]));

                if let Some((frag, mz, peak_idx)) = best {
                    let peak_mass = query.masses[peak_idx];
                    let peak_intensity = query.intensities[peak_idx];
                    let fragment_charge = query.charges[peak_idx].max(charge);

                    score.ppm_difference +=
                        peak_intensity * (mz - peak_mass).abs() * 2E6 / (mz + peak_mass);
                    score.signed_ppm_difference +=
                        peak_intensity * (peak_mass - mz) * 2E6 / (mz + peak_mass);

                    let exp_mz = query.peak_mz(peak_idx);
                    let calc_mz = frag.monoisotopic_mass / fragment_charge as f32 + PROTON;

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
                        fragments_details.charges.push(fragment_charge as i32);
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
mod tests {
    use super::*;
    use crate::database::Builder;
    use crate::enzyme::Digest;
    use crate::ion_series::IonSeries;
    use crate::modification::{ModificationDefinition, ModificationSpecificity, NeutralLossMode};
    use crate::peptide::Peptide;
    use std::{collections::HashMap, sync::Arc};

    #[test]
    fn score_ordering_uses_hyperscore() {
        let database_first = Score {
            peptide: PeptideIx(0),
            hyperscore: 1.0,
            ..Default::default()
        };
        let database_last = Score {
            peptide: PeptideIx(10_000),
            hyperscore: 100.0,
            ..Default::default()
        };

        assert!(database_last > database_first);
        assert_eq!(
            database_last.partial_cmp(&database_first),
            Some(database_last.cmp(&database_first))
        );

        let mut candidates = vec![database_first, database_last];
        bounded_min_heapify(&mut candidates, 1);
        assert_eq!(candidates[0].peptide, database_last.peptide);
    }

    #[test]
    fn longest_series() {
        let mut run = Run::default();

        run.matched(1);
        run.matched(2);
        run.matched(3);
        run.matched(3);
        run.matched(3);

        assert_eq!(run.length, 3);
        assert_eq!(run.longest, 3);

        run.matched(5);
        run.matched(5);
        assert_eq!(run.length, 1);
        assert_eq!(run.longest, 3);
        run.matched(6);
        assert_eq!(run.length, 2);
    }

    #[test]
    fn test_max_fragment_charge() {
        assert_eq!(max_fragment_charge(None, 1), 2);
        assert_eq!(max_fragment_charge(None, 2), 2);
        assert_eq!(max_fragment_charge(None, 3), 3);
        assert_eq!(max_fragment_charge(None, 4), 4);
        assert_eq!(max_fragment_charge(Some(1), 2), 2);
        assert_eq!(max_fragment_charge(Some(1), 3), 2);
        assert_eq!(max_fragment_charge(Some(2), 4), 3);
        assert_eq!(max_fragment_charge(Some(4), 1), 2);
    }

    #[test]
    fn equal_nonzero_isotope_bounds_are_honored() {
        let peptide = crate::peptide::Peptide::try_from(Digest {
            sequence: "PEPTIDER".into(),
            protein: Arc::from("protein"),
            ..Digest::default()
        })
        .unwrap();
        let fragment_masses = [Kind::B, Kind::Y]
            .into_iter()
            .flat_map(|kind| IonSeries::new(&peptide, kind))
            .map(|ion| ion.monoisotopic_mass)
            .collect::<Vec<_>>();

        for use_bitmap in [false, true] {
            let mut parameters = Builder::default().make_parameters();
            parameters.use_bitmap = use_bitmap;
            let database = parameters.build_from_peptides(vec![peptide.clone()]);
            let precursor_charge = 2;
            let precursor = Precursor {
                mz: (peptide.monoisotopic + NEUTRON) / precursor_charge as f32 + PROTON,
                charge: Some(precursor_charge),
                ..Precursor::default()
            };
            let mut query = ProcessedSpectrum {
                level: 2,
                id: "isotope-test".into(),
                precursors: vec![precursor],
                masses: fragment_masses.clone(),
                intensities: vec![1.0; fragment_masses.len()],
                charges: vec![1; fragment_masses.len()],
                total_ion_current: fragment_masses.len() as f32,
                ..ProcessedSpectrum::default()
            };
            query.masses.sort_by(f32::total_cmp);

            let scorer = Scorer {
                db: &database,
                precursor_tol: Tolerance::Da(-0.01, 0.01),
                fragment_tol: Tolerance::Da(-0.01, 0.01),
                min_matched_peaks: 1,
                min_isotope_err: 1,
                max_isotope_err: 1,
                min_precursor_charge: 2,
                max_precursor_charge: 2,
                override_precursor_charge: false,
                max_fragment_charge: Some(1),
                chimera: false,
                report_psms: 1,
                wide_window: false,
                annotate_matches: false,
                mass_shift_ppm: crate::ambiguity::DEFAULT_MASS_SHIFT_PPM,
                score_type: ScoreType::SageHyperScore,
                use_bitmap,
            };

            let features = scorer.score(&query);
            assert_eq!(features.len(), 1, "use_bitmap={use_bitmap}");
            assert_eq!(features[0].isotope_error, NEUTRON);
        }
    }

    #[test]
    fn neutral_loss_alternatives_count_once_per_cleavage_and_charge() {
        let modification = Arc::new(ModificationDefinition {
            mass: 20.0,
            name: Some(Arc::from("TestMod")),
            neutral_losses: Arc::from([10.0]),
            neutral_loss_mode: NeutralLossMode::Optional,
        });
        let peptide = Peptide::try_from(Digest {
            sequence: "AMK".into(),
            ..Default::default()
        })
        .unwrap()
        .apply(
            &[(
                ModificationSpecificity::Residue(b'M'),
                modification,
                Some(1),
            )],
            &HashMap::default(),
            1,
            None,
        )
        .into_iter()
        .find(|peptide| peptide.to_string().contains("TestMod"))
        .unwrap();

        let group = IonGroupSeries::new(&peptide, Kind::B).nth(1).unwrap();
        assert_eq!(group.variants.len(), 2);
        let mut variants = group.variants;
        variants.sort_by(|a, b| a.monoisotopic_mass.total_cmp(&b.monoisotopic_mass));

        let db = IndexedDatabase {
            peptides: vec![peptide],
            ion_kinds: vec![Kind::B],
            ..Default::default()
        };
        let scorer = Scorer {
            db: &db,
            precursor_tol: Tolerance::Da(-0.01, 0.01),
            fragment_tol: Tolerance::Da(-0.01, 0.01),
            min_matched_peaks: 1,
            min_isotope_err: 0,
            max_isotope_err: 0,
            min_precursor_charge: 2,
            max_precursor_charge: 2,
            override_precursor_charge: false,
            max_fragment_charge: Some(1),
            chimera: false,
            report_psms: 1,
            wide_window: false,
            annotate_matches: true,
            mass_shift_ppm: crate::ambiguity::DEFAULT_MASS_SHIFT_PPM,
            score_type: ScoreType::SageHyperScore,
            use_bitmap: false,
        };
        let query = ProcessedSpectrum {
            masses: variants
                .iter()
                .map(|variant| variant.monoisotopic_mass)
                .collect(),
            intensities: vec![100.0, 10.0],
            charges: vec![1, 1],
            total_ion_current: 110.0,
            ..Default::default()
        };
        let pre_score = PreScore {
            peptide: PeptideIx(0),
            precursor_charge: 2,
            ..Default::default()
        };

        let (score, fragments, _) = scorer.score_candidate(&query, &pre_score, true);
        assert_eq!(score.matched_b, 1);
        assert_eq!(score.summed_b, 100.0);
        let fragments = fragments.unwrap();
        assert_eq!(fragments.fragment_ordinals.len(), 1);
        assert_eq!(fragments.neutral_losses, vec![10.0]);

        let deferred = scorer.annotate_candidate(
            &query,
            &Feature {
                peptide_idx: PeptideIx(0),
                charge: 2,
                ..Default::default()
            },
        );
        assert_eq!(deferred.kinds, fragments.kinds);
        assert_eq!(deferred.charges, fragments.charges);
        assert_eq!(deferred.fragment_ordinals, fragments.fragment_ordinals);
        assert_eq!(deferred.intensities, fragments.intensities);
        assert_eq!(deferred.mz_calculated, fragments.mz_calculated);
        assert_eq!(deferred.mz_experimental, fragments.mz_experimental);
        assert_eq!(deferred.neutral_losses, fragments.neutral_losses);
    }

    #[test]
    fn deferred_chimera_annotation_replays_filtered_preceding_ranks() {
        let peptide = Peptide::try_from(Digest {
            sequence: "PEPTIDER".into(),
            ..Default::default()
        })
        .unwrap();
        let mut peaks = [Kind::B, Kind::Y]
            .into_iter()
            .flat_map(|kind| IonSeries::new(&peptide, kind))
            .map(|ion| (ion.monoisotopic_mass, 100.0, 1))
            .collect::<Vec<_>>();
        peaks.sort_unstable_by(|a, b| a.0.total_cmp(&b.0));
        let query = ProcessedSpectrum {
            level: 2,
            masses: peaks.iter().map(|peak| peak.0).collect(),
            intensities: peaks.iter().map(|peak| peak.1).collect(),
            charges: peaks.iter().map(|peak| peak.2).collect(),
            total_ion_current: peaks.iter().map(|peak| peak.1).sum(),
            ..Default::default()
        };
        let database = IndexedDatabase {
            peptides: vec![peptide],
            ion_kinds: vec![Kind::B, Kind::Y],
            ..Default::default()
        };
        let scorer = |chimera| Scorer {
            db: &database,
            precursor_tol: Tolerance::Da(-0.01, 0.01),
            fragment_tol: Tolerance::Da(-0.01, 0.01),
            min_matched_peaks: 1,
            min_isotope_err: 0,
            max_isotope_err: 0,
            min_precursor_charge: 2,
            max_precursor_charge: 2,
            override_precursor_charge: false,
            max_fragment_charge: Some(1),
            chimera,
            report_psms: 2,
            wide_window: false,
            annotate_matches: false,
            mass_shift_ppm: crate::ambiguity::DEFAULT_MASS_SHIFT_PPM,
            score_type: ScoreType::SageHyperScore,
            use_bitmap: false,
        };
        let rank_one = Feature {
            peptide_idx: PeptideIx(0),
            charge: 2,
            rank: 1,
            ..Default::default()
        };
        let rank_two = Feature {
            rank: 2,
            ..rank_one.clone()
        };
        let features = [&rank_one, &rank_two];

        let replayed = scorer(true).annotate_ranked_candidates(&query, &features, &[false, true]);
        assert!(replayed[0].is_none());
        assert_eq!(
            replayed[1].as_ref().unwrap().fragment_ordinals.len(),
            0,
            "rank one must remove its peaks even when it is filtered from output"
        );

        let independent =
            scorer(false).annotate_ranked_candidates(&query, &features, &[false, true]);
        assert!(!independent[1]
            .as_ref()
            .unwrap()
            .fragment_ordinals
            .is_empty());
    }
}
