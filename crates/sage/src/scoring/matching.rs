//! Fragment peak matching: the prefilter, the initial hits and their trimming.
//!
//! Moved out of `scoring.rs` unchanged; every item keeps its `scoring::` path.

use crate::database::{FragmentWindowPlan, Theoretical};
use crate::heap::bounded_min_heapify;
use crate::mass::{Tolerance, NEUTRON, PROTON};
use crate::spectrum::{Precursor, ProcessedSpectrum};
use std::ops::AddAssign;
use std::sync::atomic::{AtomicU64, Ordering};

use super::*;

/// A compact thread-safe set used while spectra are searched in parallel.
pub struct AtomicBitSet {
    words: Vec<AtomicU64>,
    len: usize,
}

impl AtomicBitSet {
    pub fn new(len: usize) -> Self {
        let words = (0..len.div_ceil(u64::BITS as usize))
            .map(|_| AtomicU64::new(0))
            .collect();
        Self { words, len }
    }

    #[inline]
    pub fn insert(&self, index: usize) {
        assert!(index < self.len);
        let word = index / u64::BITS as usize;
        let bit = index % u64::BITS as usize;
        self.words[word].fetch_or(1u64 << bit, Ordering::Relaxed);
    }

    #[inline]
    pub fn contains(&self, index: usize) -> bool {
        assert!(index < self.len);
        let word = index / u64::BITS as usize;
        let bit = index % u64::BITS as usize;
        self.words[word].load(Ordering::Relaxed) & (1u64 << bit) != 0
    }

    pub fn len(&self) -> usize {
        self.len
    }

    pub fn is_empty(&self) -> bool {
        self.len == 0
    }
}

#[derive(Copy, Clone)]
pub(crate) struct FragmentMatchPeak {
    pub(crate) neutral_mass: f32,
    pub(crate) query_index: usize,
    charge: u8,
    charge_is_known: bool,
}

pub(crate) struct FragmentMatchIndex {
    pub(crate) peaks: Vec<FragmentMatchPeak>,
}

impl FragmentMatchIndex {
    pub(crate) fn new(query: &ProcessedSpectrum, max_charge: u8) -> Self {
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

    pub(super) fn select_peak(
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
pub(super) struct InitialHits {
    pub(super) matched_peaks: usize,
    // Number of peptide candidates with > 0 matched peaks
    pub(super) scored_candidates: usize,
    pub(super) preliminary: Vec<PreScore>,
}

impl AddAssign<InitialHits> for InitialHits {
    fn add_assign(&mut self, rhs: InitialHits) {
        self.matched_peaks += rhs.matched_peaks;
        self.scored_candidates += rhs.scored_candidates;

        self.preliminary.extend(rhs.preliminary);
    }
}

impl<'db> Scorer<'db> {
    /// Mark every peptide that contributes at least one preliminary fragment
    /// match for this spectrum. Keeping all contributors preserves the normal
    /// top-K candidate set and its candidate statistics after chunk merging.
    pub fn exact_prefilter(&self, query: &ProcessedSpectrum, keep: &AtomicBitSet) {
        assert_eq!(
            query.level, 2,
            "internal bug, trying to score a non-MS2 scan!"
        );
        assert_eq!(keep.len(), self.db.peptides.len());
        let query = self.recalibrated(query);
        let query = query.as_ref();
        // Readers keep MS2 scans with no selected precursor on purpose (see
        // `spectrum.rs`); in-tree callers filter them with `is_searchable`, but a
        // library caller can hand one straight here. Such a scan has no
        // precursor window, so it marks no candidate rather than panicking.
        let Some(precursor) = query.precursors.first() else {
            return;
        };
        let mz = precursor.mz - PROTON;

        let mark = |precursor_mass: f32, precursor_charge: u8, tolerance: Tolerance| {
            // A non-finite precursor has no finite candidate window.
            if !precursor_mass.is_finite() {
                return;
            }
            let fragment_index = FragmentMatchIndex::new(
                query,
                max_fragment_charge(self.max_fragment_charge, precursor_charge),
            );
            let mut plan = FragmentWindowPlan::default();
            for offset in self.offsets() {
                // The peaks, the fragment tolerance and the shift are the same
                // for every isotope error, so the windows and their buckets are
                // planned once per offset.
                let shift = self.fragment_shift(offset);
                self.plan_windows(&fragment_index, shift.as_slice(), &mut plan);
                for isotope_error in self.min_isotope_err..=self.max_isotope_err {
                    let (mass, tolerance) = self.offset_query(
                        precursor_mass - isotope_error as f32 * NEUTRON,
                        tolerance,
                        offset,
                    );
                    let candidates = self.db.query(mass, tolerance, self.fragment_tol);
                    candidates.page_search_planned(&mut plan, |fragment| {
                        keep.insert(fragment.peptide_index.0 as usize);
                    });
                }
            }
        };

        if self.wide_window {
            for charge in self.min_precursor_charge..=self.max_precursor_charge {
                let tolerance = precursor
                    .isolation_window
                    .unwrap_or(Tolerance::Da(-2.4, 2.4))
                    * charge as f32;
                mark(mz * charge as f32, charge, tolerance);
            }
        } else if let Some(charge) = precursor.charge.filter(|_| !self.override_precursor_charge) {
            mark(mz * charge as f32, charge, self.precursor_tol);
        } else {
            for charge in self.min_precursor_charge..=self.max_precursor_charge {
                mark(mz * charge as f32, charge, self.precursor_tol);
            }
        }
    }
    /// Plan the fragment windows of one spectrum's peaks, including the
    /// shifted lookups, so several precursor windows can reuse them.
    pub(super) fn plan_windows(
        &self,
        fragment_index: &FragmentMatchIndex,
        shifts: &[f32],
        plan: &mut FragmentWindowPlan,
    ) {
        plan.build(
            self.db,
            self.fragment_tol,
            fragment_index.peaks.iter().map(|peak| peak.neutral_mass),
            shifts,
        );
    }

    /// Keep the top-K preliminary candidates of each offset hypothesis. Shifted
    /// hypotheses combine two fragment lookups, so their preliminary counts are
    /// not comparable with unshifted counts until full scoring.
    fn trim_hits_by_offset(&self, hits: &mut InitialHits) {
        if self.db.mass_offsets.is_empty() {
            self.trim_hits(hits);
            return;
        }
        let mut preliminary = std::mem::take(&mut hits.preliminary);
        preliminary.sort_unstable_by_key(|score| score.offset);
        let mut trimmed = Vec::with_capacity(preliminary.len());
        for group in preliminary.chunk_by_mut(|left, right| left.offset == right.offset) {
            let k = 50.clamp((self.report_psms * 2).min(group.len()), group.len());
            bounded_min_heapify(group, k);
            trimmed.extend_from_slice(&group[..k]);
        }
        hits.preliminary = trimmed;
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
        plan: &mut FragmentWindowPlan,
        precursor_mass: f32,
        precursor_charge: u8,
        precursor_tol: Tolerance,
        isotope_error: i8,
        offset: u8,
    ) -> InitialHits {
        let (mass, tolerance) = self.offset_query(
            precursor_mass - isotope_error as f32 * NEUTRON,
            precursor_tol,
            offset,
        );
        self.preliminary_hits(
            plan,
            mass,
            tolerance,
            PreScore {
                precursor_charge,
                isotope_error,
                offset,
                ..Default::default()
            },
        )
    }

    /// Count preliminary fragment matches for every indexed peptide whose
    /// mass lies in `tolerance` around `base_mass`. Fragments not carrying an
    /// offset match the indexed masses; fragments carrying it match at the
    /// indexed mass plus any of `shifts`. `template` supplies the charge,
    /// isotope and offset recorded on each hit.
    pub(super) fn preliminary_hits(
        &self,
        plan: &mut FragmentWindowPlan,
        base_mass: f32,
        tolerance: Tolerance,
        template: PreScore,
    ) -> InitialHits {
        let PreScore {
            precursor_charge,
            isotope_error,
            offset,
            ..
        } = template;
        let candidates = self.db.query(base_mass, tolerance, self.fragment_tol);

        let potential = candidates.pre_idx_hi - candidates.pre_idx_lo + 1;

        PRE_SCORE_SCRATCH.with(|scratch| {
            let mut preliminary = scratch.borrow_mut();
            preliminary.resize(potential, PreScore::default());
            preliminary.fill(PreScore::default());
            let mut matched_peaks = 0;
            let mut scored_candidates = 0;

            let mut count = |frag: Theoretical| {
                let idx = frag.peptide_index.0 as usize - candidates.pre_idx_lo;
                let sc = &mut preliminary[idx];
                if sc.matched == 0 {
                    scored_candidates += 1;
                    sc.precursor_charge = precursor_charge;
                    sc.peptide = frag.peptide_index;
                    sc.isotope_error = isotope_error;
                    sc.offset = offset;
                }

                sc.matched += 1;
                matched_peaks += 1;
            };
            // Fragments retaining the offset appear at the indexed base
            // fragment mass plus each shift; `plan` holds both lookups.
            candidates.page_search_planned(plan, &mut count);

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

    pub(super) fn matched_peaks(
        &self,
        query: &ProcessedSpectrum,
        precursor_mass: f32,
        precursor_charge: u8,
        precursor_tol: Tolerance,
    ) -> InitialHits {
        let max_fragment_charge = max_fragment_charge(self.max_fragment_charge, precursor_charge);
        let fragment_index = FragmentMatchIndex::new(query, max_fragment_charge);
        let mut plan = FragmentWindowPlan::default();
        self.offsets()
            .fold(InitialHits::default(), |mut hits, offset| {
                hits += self.matched_peaks_with_offset(
                    &fragment_index,
                    &mut plan,
                    precursor_mass,
                    precursor_charge,
                    precursor_tol,
                    offset,
                );
                hits
            })
    }

    fn matched_peaks_with_offset(
        &self,
        fragment_index: &FragmentMatchIndex,
        plan: &mut FragmentWindowPlan,
        precursor_mass: f32,
        precursor_charge: u8,
        precursor_tol: Tolerance,
        offset: u8,
    ) -> InitialHits {
        // One plan for every isotope error of this offset: they search the
        // same peaks and differ only in the precursor window.
        let shift = self.fragment_shift(offset);
        self.plan_windows(fragment_index, shift.as_slice(), plan);
        if self.min_isotope_err != self.max_isotope_err {
            let mut hits = (self.min_isotope_err..=self.max_isotope_err).fold(
                InitialHits::default(),
                |mut hits, isotope| {
                    hits += self.matched_peaks_with_isotope(
                        plan,
                        precursor_mass,
                        precursor_charge,
                        precursor_tol,
                        isotope,
                        offset,
                    );
                    hits
                },
            );
            self.trim_hits(&mut hits);
            hits
        } else {
            self.matched_peaks_with_isotope(
                plan,
                precursor_mass,
                precursor_charge,
                precursor_tol,
                self.min_isotope_err,
                offset,
            )
        }
    }

    pub(super) fn initial_hits(
        &self,
        query: &ProcessedSpectrum,
        precursor: &Precursor,
    ) -> InitialHits {
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
        self.trim_hits_by_offset(&mut hits);
        hits
    }
}
