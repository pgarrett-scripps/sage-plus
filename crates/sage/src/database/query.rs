//! Fragment index queries and the per-thread search scratch space.
//!
//! Moved out of `database.rs` unchanged; every item keeps its `database::` path.

use crate::mass::Tolerance;
use std::cmp::Ordering;

use super::*;

/// One spectrum's fragment windows and the buckets they touch, each distinct
/// bucket listed once.
///
/// A plan depends on the peak list, the fragment tolerance and the shifts, and
/// not on the precursor window, so every isotope-error and charge hypothesis
/// of one spectrum can share it: only the per-bucket precursor ranges have to
/// be resolved again. Its allocations are reused by every `build`.
#[derive(Default)]
pub(crate) struct FragmentWindowPlan {
    /// Per window: fragment bounds and its span in `entries`.
    windows: Vec<(f32, f32, u32, u32)>,
    /// One bucket per (window, bucket) pair, in emission order.
    entries: Vec<u32>,
    /// `(bucket, entry)` pairs sorted by bucket for de-duplication.
    order: Vec<(u32, u32)>,
    /// Distinct buckets and their resolved ranges.
    distinct: Vec<u32>,
    resolved: Vec<(usize, usize)>,
    /// Resolved range per entry.
    ranges: Vec<(usize, usize)>,
    /// Peptide-index range per precursor window, for the batched resolve.
    peptide_windows: Vec<(u32, u32)>,
}

impl FragmentWindowPlan {
    /// Lay out one window per mass and one per shifted mass, and list the
    /// distinct buckets they touch.
    pub(crate) fn build(
        &mut self,
        db: &IndexedDatabase,
        fragment_tol: Tolerance,
        masses: impl IntoIterator<Item = f32>,
        shifts: &[f32],
    ) {
        let Self {
            windows,
            entries,
            order,
            distinct,
            ..
        } = self;
        windows.clear();
        entries.clear();
        let mut push_window = |fragment_lo: f32, fragment_hi: f32| {
            let (left, right) = fragment_bucket_slice(&db.min_value, fragment_lo, fragment_hi);
            let first = entries.len() as u32;
            entries.extend(left as u32..right as u32);
            windows.push((fragment_lo, fragment_hi, first, entries.len() as u32));
        };
        for mass in masses {
            let (fragment_lo, fragment_hi) = fragment_tol.bounds(mass);
            push_window(fragment_lo, fragment_hi);
            for &shift in shifts {
                push_window(fragment_lo - shift, fragment_hi - shift);
            }
        }

        // Each distinct bucket is resolved once for the whole spectrum.
        order.clear();
        order.extend(
            entries
                .iter()
                .enumerate()
                .map(|(entry, &bucket)| (bucket, entry as u32)),
        );
        order.sort_unstable();
        distinct.clear();
        distinct.extend(order.iter().map(|&(bucket, _)| bucket));
        distinct.dedup();
    }
}

thread_local! {
    /// Plan for [`IndexedQuery::page_search_batch`], which plans and searches
    /// in one call. Callers that search the same peaks against several
    /// precursor windows keep their own plan and call
    /// [`IndexedQuery::page_search_planned`]; nothing may call
    /// `page_search_batch` while this plan is borrowed.
    static BATCH_SEARCH_PLAN: std::cell::RefCell<FragmentWindowPlan> =
        std::cell::RefCell::new(FragmentWindowPlan::default());
}

pub(super) struct FragmentIter<'a> {
    pub(super) fragments: &'a [PackedFragment],
    pub(super) mass_prefix: u32,
    pub(super) next: usize,
}

impl Iterator for FragmentIter<'_> {
    type Item = Theoretical;

    #[inline(always)]
    fn next(&mut self) -> Option<Self::Item> {
        let fragment = self.fragments.get(self.next).copied()?;
        self.next += 1;
        Some(fragment.decode(self.mass_prefix))
    }
}

pub struct IndexedQuery<'d> {
    pub(super) db: &'d IndexedDatabase,
    pub(super) precursor_mass: f32,
    pub(super) precursor_tol: Tolerance,
    pub(super) fragment_tol: Tolerance,
    pub pre_idx_lo: usize,
    pub pre_idx_hi: usize,
}

impl IndexedQuery<'_> {
    /// Search for a specified `fragment_mz` within the database
    pub fn page_search(&self, mass: f32) -> impl Iterator<Item = Theoretical> + '_ {
        self.page_search_shifted(mass, 0.0)
    }

    /// Search for indexed fragments at `mass - shift`. The tolerance window is
    /// evaluated at the observed mass and then translated, so a shifted lookup
    /// has the same width as matching the modified theoretical fragment.
    pub fn page_search_shifted(
        &self,
        mass: f32,
        shift: f32,
    ) -> impl Iterator<Item = Theoretical> + '_ {
        let (fragment_lo, fragment_hi) = self.fragment_tol.bounds(mass);
        let (fragment_lo, fragment_hi) = (fragment_lo - shift, fragment_hi - shift);
        let (precursor_lo, precursor_hi) = self.precursor_tol.bounds(self.precursor_mass);

        // Locate the mass-prefix buckets that can contain matching fragments.
        let (left_idx, right_idx) =
            fragment_bucket_slice(&self.db.min_value, fragment_lo, fragment_hi);

        let peptide_lo = self.pre_idx_lo.min(u32::MAX as usize) as u32;
        let peptide_hi = self.pre_idx_hi.min(u32::MAX as usize) as u32;

        // It is absolutely critical that we do not cross page boundaries.
        // Otherwise we can no longer rely on peptide-index ordering.
        (left_idx..right_idx).flat_map(move |page| {
            self.db
                .fragments
                .bucket_search(page, peptide_lo, peptide_hi)
                .filter(move |frag| {
                    self.accepts(frag, fragment_lo, fragment_hi, precursor_lo, precursor_hi)
                })
        })
    }

    /// Final per-fragment check shared by the single and batched searches.
    #[inline(always)]
    fn accepts(
        &self,
        frag: &Theoretical,
        fragment_lo: f32,
        fragment_hi: f32,
        precursor_lo: f32,
        precursor_hi: f32,
    ) -> bool {
        // This looks somewhat complicated, but it's a consequence of
        // how the `binary_search_slice` function works - it will return
        // the set of indices that maximally cover the desired range - the exact
        // `left` and `right` indices may be valid, or just outside of the range.
        // Anything interior of `left` and `right` is guaranteed to be within the
        // precursor tolerance, so we just need to check the edge cases
        //
        // Previously, a direct lookup to check the mass of the current fragment was
        // performed, but the pointer indirection + float comparison can slow down
        // open searches by as much as 2x!!
        // e.g. used to be `self.db[frag.peptide_index].monoisotopic >= precursor_lo`
        (frag.peptide_index.0 > self.pre_idx_lo as u32
            || (frag.peptide_index.0 == self.pre_idx_lo as u32
                && self.db[frag.peptide_index].monoisotopic >= precursor_lo))
            && (frag.peptide_index.0 < self.pre_idx_hi as u32
                || (frag.peptide_index.0 == self.pre_idx_hi as u32
                    && self.db[frag.peptide_index].monoisotopic <= precursor_hi))
            && frag.fragment_mz >= fragment_lo
            && frag.fragment_mz <= fragment_hi
    }

    /// Batched equivalent of calling [`Self::page_search`] and then
    /// [`Self::page_search_shifted`] for every shift, for each mass in turn.
    ///
    /// Matches are emitted in exactly that order. Every bucket's
    /// precursor-scoped range depends only on the bucket, because the
    /// peptide range is fixed for this query, so each distinct bucket is
    /// resolved once for the whole spectrum and all of them are resolved
    /// with interleaved binary searches. Adapted from Matteo Lacki's
    /// `page_search_batch` (MatteoLacki/sage 062f7b3 and 950641f, MIT).
    pub fn page_search_batch(
        &self,
        masses: impl IntoIterator<Item = f32>,
        shifts: &[f32],
        on_match: impl FnMut(Theoretical),
    ) {
        BATCH_SEARCH_PLAN.with(|plan| {
            let mut plan = plan.borrow_mut();
            plan.build(self.db, self.fragment_tol, masses, shifts);
            self.page_search_planned(&mut plan, on_match);
        });
    }

    /// [`Self::page_search_batch`] against a plan built once for the spectrum,
    /// so a second precursor window over the same peaks pays only for
    /// resolving the buckets again.
    pub(crate) fn page_search_planned(
        &self,
        plan: &mut FragmentWindowPlan,
        mut on_match: impl FnMut(Theoretical),
    ) {
        Self::page_search_planned_multi(std::slice::from_ref(self), plan, |_, frag| on_match(frag));
    }

    /// [`Self::page_search_planned`] for several precursor windows over one
    /// plan, resolving every bucket for every window in a single pass and
    /// then emitting each window's matches with its index.
    ///
    /// `queries` must share the database, the fragment tolerance and the peaks
    /// the plan was built from, and differ only in the precursor window --
    /// which is what the isotope-error hypotheses of one spectrum are. Each
    /// window's matches are emitted in the order a single
    /// [`Self::page_search_planned`] call would emit them, windows in order.
    pub(crate) fn page_search_planned_multi(
        queries: &[Self],
        plan: &mut FragmentWindowPlan,
        mut on_match: impl FnMut(usize, Theoretical),
    ) {
        let Some(first) = queries.first() else {
            return;
        };
        let db = first.db;
        let FragmentWindowPlan {
            windows,
            entries,
            order,
            distinct,
            resolved,
            ranges,
            peptide_windows,
        } = plan;

        peptide_windows.clear();
        peptide_windows.extend(queries.iter().map(|query| {
            (
                query.pre_idx_lo.min(u32::MAX as usize) as u32,
                query.pre_idx_hi.min(u32::MAX as usize) as u32,
            )
        }));
        db.fragments
            .resolve_bucket_ranges_multi(distinct, peptide_windows, resolved);

        ranges.clear();
        ranges.resize(entries.len(), (0, 0));
        for (index, query) in queries.iter().enumerate() {
            let (precursor_lo, precursor_hi) = query.precursor_tol.bounds(query.precursor_mass);
            let resolved = &resolved[index * distinct.len()..(index + 1) * distinct.len()];
            let mut group = 0;
            for &(bucket, entry) in order.iter() {
                if distinct[group] != bucket {
                    group += 1;
                }
                ranges[entry as usize] = resolved[group];
            }

            for &(fragment_lo, fragment_hi, first, last) in windows.iter() {
                for entry in first as usize..last as usize {
                    let (lo, hi) = ranges[entry];
                    for frag in db.fragments.range_iter(entries[entry] as usize, lo, hi) {
                        if query.accepts(
                            &frag,
                            fragment_lo,
                            fragment_hi,
                            precursor_lo,
                            precursor_hi,
                        ) {
                            on_match(index, frag);
                        }
                    }
                }
            }
        }
    }
}

pub(super) fn fragment_bucket_slice(min_values: &[f32], low: f32, high: f32) -> (usize, usize) {
    let (mut left, right) =
        binary_search_slice(min_values, |min, bounds| min.total_cmp(bounds), low, high);
    if let Some(&left_value) = min_values.get(left) {
        left = min_values.partition_point(|value| value.total_cmp(&left_value).is_lt());
    }
    (left, right)
}

/// Return the widest `left` and `right` indices into a `slice` (sorted by the
/// function `key`) such that all values between `low` and `high` are
/// contained in `slice[left..right]`
///
/// # Invariants
///
/// * `slice[left] <= low || left == 0`
/// * `slice[right] > high || right == slice.len()`
/// * `0 <= left <= right <= slice.len()`
#[inline]
pub fn binary_search_slice<T, F, S>(slice: &[T], key: F, low: S, high: S) -> (usize, usize)
where
    F: Fn(&T, &S) -> Ordering,
{
    let left_idx = slice
        .partition_point(|a| key(a, &low) == Ordering::Less)
        .saturating_sub(1);

    let right_idx =
        slice[left_idx..].partition_point(|a| key(a, &high) != Ordering::Greater) + left_idx;

    (left_idx, right_idx)
}
