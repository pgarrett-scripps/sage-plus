//! The fragment index and the `IndexedDatabase` built on top of it.
//!
//! Moved out of `database.rs` unchanged; every item keeps its `database::` path.

use crate::cleavage::ValidatedCustomCleavageLibrary;
use crate::fasta::Fasta;
use crate::fragment_loss::FragmentLosses;
use crate::ion_series::{IonGroupSeries, Kind};
use crate::mass::Tolerance;
use crate::modification::{ModificationDefinition, ModificationSpecificity, SiteMode};
use crate::peptide::{Peptide, Site};
use crate::ptm_library::PtmLibrary;
use crate::scoring::Feature;
use rayon::prelude::*;
use serde::Serialize;
use std::borrow::Cow;
use std::collections::HashMap;
use std::hash::Hash;
use std::sync::Arc;

use super::*;

#[derive(Hash, Copy, Clone, PartialEq, Eq, PartialOrd, Ord, Debug, Serialize)]
#[repr(transparent)]
pub struct PeptideIx(pub u32);

// This is unsafe for use outside of this crate
impl Default for PeptideIx {
    fn default() -> Self {
        Self(u32::MAX)
    }
}

/// Extra preliminary fragment masses of a peptide, appended to the buffer.
/// See [`Parameters::build_from_peptides_with_extra_fragments`].
pub type ExtraFragments = dyn Fn(&Peptide, &mut Vec<f32>) + Sync;

/// Marks the peaks of a spectrum that full scoring must not match as b/y
/// ions of a candidate. Called with the resolved peptidoform, the precursor
/// charge, the spectrum and a mask with one entry per peak (all false).
/// See [`IndexedDatabase::peak_exclusion`].
pub type PeakExclusion =
    dyn Fn(&Peptide, u8, &crate::spectrum::ProcessedSpectrum, &mut [bool]) + Send + Sync;

#[derive(Copy, Clone, Debug, PartialEq, Serialize)]
pub struct Theoretical {
    pub peptide_index: PeptideIx,
    pub fragment_mz: f32,
}

/// Lower `f32` mass bits stored in each packed fragment. The remaining upper
/// bits are shared by every search bucket with the same mass prefix.
pub(super) const FRAGMENT_MASS_SUFFIX_BITS: u32 = 12;

#[repr(C, packed)]
#[derive(Copy, Clone)]
pub(super) struct PackedFragment {
    pub(super) peptide_index: u32,
    pub(super) mass_suffix: u16,
}

impl PackedFragment {
    #[inline(always)]
    fn peptide_index(self) -> u32 {
        self.peptide_index
    }

    #[inline(always)]
    pub(super) fn decode(self, mass_prefix: u32) -> Theoretical {
        Theoretical {
            peptide_index: PeptideIx(self.peptide_index),
            fragment_mz: f32::from_bits(mass_prefix | u32::from(self.mass_suffix)),
        }
    }
}

#[derive(Copy, Clone)]
pub(super) struct FragmentBucket {
    pub(super) mass_prefix: u32,
    pub(super) start: usize,
    pub(super) end: usize,
}

#[derive(Default)]
/// Lossless theoretical-fragment index using exact packed mass-prefix buckets.
pub struct FragmentIndex {
    pub(super) buckets: Vec<FragmentBucket>,
    pub(super) fragments: Vec<PackedFragment>,
}

#[derive(Copy, Clone)]
struct SharedFragmentWriter(*mut std::mem::MaybeUninit<PackedFragment>);

// SAFETY: each parallel peptide chunk receives precomputed, non-overlapping
// output ranges for every mass prefix. No output position has multiple writers.
unsafe impl Send for SharedFragmentWriter {}
// SAFETY: writes through the shared pointer are confined to the disjoint ranges
// described above, and the allocation remains alive until all workers finish.
unsafe impl Sync for SharedFragmentWriter {}

impl SharedFragmentWriter {
    unsafe fn write(self, index: usize, fragment: PackedFragment) {
        // SAFETY: callers use the disjoint and in-bounds positions assigned by
        // the counting pass.
        unsafe {
            self.0
                .add(index)
                .write(std::mem::MaybeUninit::new(fragment))
        };
    }
}

impl FragmentIndex {
    fn build(
        parameters: &Parameters,
        peptides: &[Peptide],
        extra: Option<&ExtraFragments>,
    ) -> (Self, Vec<f32>) {
        let suffix_bits = FRAGMENT_MASS_SUFFIX_BITS;
        let suffix_mask = (1u32 << suffix_bits) - 1;
        let bucket_size = parameters.bucket_size.max(1);
        let requested_chunks = rayon::current_num_threads().max(1) * 4;
        let chunk_size = peptides.len().div_ceil(requested_chunks).max(1);
        let ranges = (0..peptides.len())
            .step_by(chunk_size)
            .map(|start| start..(start + chunk_size).min(peptides.len()))
            .collect::<Vec<_>>();

        let counts = ranges
            .par_iter()
            .map(|range| {
                let mut counts = HashMap::<u32, u32>::new();
                let mut buffer = Vec::new();
                for peptide in &peptides[range.clone()] {
                    for mass in preliminary_fragment_masses(parameters, peptide) {
                        let prefix = mass.to_bits() >> suffix_bits;
                        *counts.entry(prefix).or_default() += 1;
                    }
                    if let Some(extra) = extra {
                        buffer.clear();
                        extra(peptide, &mut buffer);
                        for mass in &buffer {
                            let prefix = mass.to_bits() >> suffix_bits;
                            *counts.entry(prefix).or_default() += 1;
                        }
                    }
                }
                counts
            })
            .collect::<Vec<_>>();

        let mut prefixes = counts
            .iter()
            .flat_map(|counts| counts.keys().copied())
            .collect::<Vec<_>>();
        prefixes.sort_unstable();
        prefixes.dedup();

        let prefix_counts = prefixes
            .into_iter()
            .map(|prefix| {
                let count = counts
                    .iter()
                    .map(|counts| counts.get(&prefix).copied().unwrap_or_default() as usize)
                    .sum::<usize>();
                (prefix, count)
            })
            .collect::<Vec<_>>();
        let bucket_count = prefix_counts
            .iter()
            .map(|(_, count)| count.div_ceil(bucket_size))
            .sum();

        let mut positions = vec![HashMap::<u32, usize>::new(); ranges.len()];
        let mut limits = vec![HashMap::<u32, usize>::new(); ranges.len()];
        let mut buckets = Vec::with_capacity(bucket_count);
        let mut total = 0usize;
        for (prefix, prefix_count) in prefix_counts {
            let start = total;
            for (chunk, counts) in counts.iter().enumerate() {
                let count = counts.get(&prefix).copied().unwrap_or_default() as usize;
                if count > 0 {
                    positions[chunk].insert(prefix, total);
                    total += count;
                    limits[chunk].insert(prefix, total);
                }
            }
            assert_eq!(total - start, prefix_count);
            for bucket_start in (start..total).step_by(bucket_size) {
                let bucket_end = (bucket_start + bucket_size).min(total);
                buckets.push(FragmentBucket {
                    mass_prefix: prefix << suffix_bits,
                    start: bucket_start,
                    end: bucket_end,
                });
            }
        }

        let mut uninitialized = Vec::<std::mem::MaybeUninit<PackedFragment>>::with_capacity(total);
        // SAFETY: the second generation pass writes every position exactly once
        // using the counts and disjoint offsets computed above.
        unsafe { uninitialized.set_len(total) };
        let writer = SharedFragmentWriter(uninitialized.as_mut_ptr());
        ranges
            .into_par_iter()
            .zip(positions.into_par_iter().zip(limits.into_par_iter()))
            .for_each(|(range, (mut positions, limits))| {
                let mut buffer = Vec::new();
                for (peptide_index, peptide) in peptides[range.clone()].iter().enumerate() {
                    let peptide_index = range.start + peptide_index;
                    let mut write = |mass: f32| {
                        let bits = mass.to_bits();
                        let prefix = bits >> suffix_bits;
                        let position = positions
                            .get_mut(&prefix)
                            .expect("counting and generation passes disagree");
                        assert!(
                            *position < limits[&prefix],
                            "fragment generation exceeded the counted output range"
                        );
                        // SAFETY: every chunk owns a disjoint position range for
                        // this prefix, and generation order stays within that range.
                        unsafe {
                            writer.write(
                                *position,
                                PackedFragment {
                                    peptide_index: u32::try_from(peptide_index)
                                        .expect("peptide index exceeds 32 bits"),
                                    mass_suffix: (bits & suffix_mask) as u16,
                                },
                            )
                        };
                        *position += 1;
                    };
                    for mass in preliminary_fragment_masses(parameters, peptide) {
                        write(mass);
                    }
                    if let Some(extra) = extra {
                        buffer.clear();
                        extra(peptide, &mut buffer);
                        for &mass in &buffer {
                            write(mass);
                        }
                    }
                }
                assert_eq!(
                    positions, limits,
                    "fragment generation did not fill every counted output position"
                );
            });

        let pointer = uninitialized.as_mut_ptr().cast::<PackedFragment>();
        let len = uninitialized.len();
        let capacity = uninitialized.capacity();
        std::mem::forget(uninitialized);
        // SAFETY: every element was initialized once above. `MaybeUninit<T>` and
        // `T` have identical allocation layouts.
        let fragments = unsafe { Vec::from_raw_parts(pointer, len, capacity) };
        let min_value = buckets
            .iter()
            .map(|bucket| f32::from_bits(bucket.mass_prefix))
            .collect();
        (Self { buckets, fragments }, min_value)
    }

    pub fn len(&self) -> usize {
        self.fragments.len()
    }

    pub fn is_empty(&self) -> bool {
        self.fragments.is_empty()
    }

    pub fn allocated_bytes(&self) -> usize {
        self.buckets.capacity() * std::mem::size_of::<FragmentBucket>()
            + self.fragments.capacity() * std::mem::size_of::<PackedFragment>()
    }

    pub(super) fn bucket_search(
        &self,
        bucket_index: usize,
        peptide_lo: u32,
        peptide_hi: u32,
    ) -> FragmentIter<'_> {
        let bucket = self.buckets[bucket_index];
        let fragments = &self.fragments[bucket.start..bucket.end];
        let first = fragments.partition_point(|fragment| fragment.peptide_index() < peptide_lo);
        let last = fragments.partition_point(|fragment| fragment.peptide_index() <= peptide_hi);
        FragmentIter {
            fragments: &fragments[first..last],
            mass_prefix: bucket.mass_prefix,
            next: 0,
        }
    }

    pub fn bucket(&self, bucket: usize) -> impl Iterator<Item = Theoretical> + '_ {
        self.bucket_search(bucket, 0, u32::MAX)
    }

    /// Resolve the precursor-scoped fragment range of many buckets at once.
    ///
    /// Writes one absolute `[first, last)` range into `fragments` per entry of
    /// `buckets`, identical to the range [`Self::bucket_search`] selects. A
    /// binary search is a chain of dependent loads, so one search at a time
    /// leaves the core waiting on a single cache miss. Here `LANES` buckets
    /// are stepped in lockstep, each with independent lower and upper
    /// searches, so up to `2 * LANES` misses are in flight together. The
    /// design follows Matteo Lacki's interleaved page-bound resolver
    /// (MatteoLacki/sage 950641f, MIT).
    ///
    /// Several precursor `windows` are resolved in the same pass, with `out`
    /// filled window-major: `out[window * buckets.len() + bucket]`. A bucket's
    /// search path is then read once from memory and from cache for the
    /// remaining windows. The isotope-error hypotheses of one spectrum are
    /// exactly this shape: the same buckets, a different precursor window
    /// each.
    pub(super) fn resolve_bucket_ranges_multi(
        &self,
        buckets: &[u32],
        windows: &[(u32, u32)],
        out: &mut Vec<(usize, usize)>,
    ) {
        /// Buckets resolved together. Two searches each, sized against the
        /// core's outstanding-miss capacity rather than a vector width.
        const LANES: usize = 16;
        const SLOTS: usize = 2 * LANES;

        out.clear();
        out.resize(buckets.len() * windows.len(), (0, 0));
        if self.fragments.is_empty() {
            for slot in out.iter_mut().enumerate() {
                let (index, range) = slot;
                let start = self.buckets[buckets[index % buckets.len()] as usize].start;
                *range = (start, start);
            }
            return;
        }
        let fragments = self.fragments.as_slice();

        for (chunk_index, chunk) in buckets.chunks(LANES).enumerate() {
            for (window, &(peptide_lo, peptide_hi)) in windows.iter().enumerate() {
                // Even slots search `id < peptide_lo`, odd slots `id <= peptide_hi`,
                // exactly the two `partition_point` predicates of `bucket_search`.
                let mut base = [0usize; SLOTS];
                let mut size = [0usize; SLOTS];
                let mut longest = 0usize;
                for (lane, &bucket) in chunk.iter().enumerate() {
                    let bucket = self.buckets[bucket as usize];
                    let (start, len) = (bucket.start, bucket.end - bucket.start);
                    // An empty bucket probes index 0 (always valid here) without
                    // moving; its result is taken from `start` below.
                    let start = if len == 0 { 0 } else { start };
                    base[2 * lane] = start;
                    base[2 * lane + 1] = start;
                    size[2 * lane] = len;
                    size[2 * lane + 1] = len;
                    longest = longest.max(len);
                }

                // Branchless halving with a fixed trip count: a converged or
                // unused slot has `half == 0`, re-reads its own `base` and stays.
                let steps = usize::BITS - longest.saturating_sub(1).leading_zeros();
                for _ in 0..steps {
                    for slot in 0..SLOTS {
                        let half = size[slot] / 2;
                        let mid = base[slot] + half;
                        // SAFETY: `mid < base + size <= bucket.end <= fragments.len()`
                        // when `size > 1`; otherwise `half == 0` and `mid == base`,
                        // which is in bounds for a non-empty bucket and is index 0
                        // for an empty or unused slot.
                        let id = unsafe { fragments.get_unchecked(mid) }.peptide_index();
                        let go_right = if slot % 2 == 0 {
                            id < peptide_lo
                        } else {
                            id <= peptide_hi
                        };
                        base[slot] = if go_right { mid } else { base[slot] };
                        size[slot] -= half;
                    }
                }

                for (lane, &bucket) in chunk.iter().enumerate() {
                    let slot = window * buckets.len() + chunk_index * LANES + lane;
                    let bucket = self.buckets[bucket as usize];
                    if bucket.start == bucket.end {
                        out[slot] = (bucket.start, bucket.start);
                        continue;
                    }
                    let (lo, hi) = (base[2 * lane], base[2 * lane + 1]);
                    let first = lo + usize::from(fragments[lo].peptide_index() < peptide_lo);
                    let last = hi + usize::from(fragments[hi].peptide_index() <= peptide_hi);
                    out[slot] = (first, last);
                }
            }
        }
    }

    #[inline(always)]
    pub(super) fn range_iter(&self, bucket: usize, first: usize, last: usize) -> FragmentIter<'_> {
        FragmentIter {
            fragments: &self.fragments[first..last],
            mass_prefix: self.buckets[bucket].mass_prefix,
            next: 0,
        }
    }
}

pub fn preliminary_fragment_masses<'a>(
    parameters: &'a Parameters,
    peptide: &'a Peptide,
) -> impl Iterator<Item = f32> + 'a {
    parameters
        .ion_kinds
        .iter()
        .flat_map(|kind| IonGroupSeries::new(peptide, *kind))
        .filter(|group| match group.kind {
            Kind::A | Kind::B | Kind::C => (group.series_index + 1) > parameters.min_ion_index,
            Kind::X | Kind::Y | Kind::Z | Kind::ZDot => {
                peptide.sequence.len().saturating_sub(1) - group.series_index
                    > parameters.min_ion_index
            }
        })
        .filter_map(|group| {
            group
                .variants
                .into_iter()
                .next()
                .map(|ion| ion.monoisotopic_mass)
        })
}

#[derive(Default)]
pub struct IndexedDatabase {
    /// Indexed peptides, sorted by monoisotopic mass.
    pub peptides: Vec<Peptide>,
    /// Offset-placed peptidoforms reported by the search. They are addressed
    /// by [`PeptideIx`] values following `peptides` and are not indexed.
    pub offset_peptides: Vec<Peptide>,
    /// Search-time mass offsets. Their masses are never part of the index.
    pub mass_offsets: Vec<MassOffset>,
    /// Site library restricting `library` mode offsets, when configured.
    pub offset_library: Option<Arc<PtmLibrary>>,
    pub fragments: FragmentIndex,
    pub ion_kinds: Vec<Kind>,
    /// Generic fragment losses matched during full candidate scoring. Never
    /// part of the fragment index.
    pub fragment_losses: Option<Arc<FragmentLosses>>,
    pub min_value: Vec<f32>,
    /// Variable modification candidates used by PTM localization.
    pub potential_mods: Vec<(ModificationSpecificity, f32)>,
    /// Full definitions keep equal-mass modifications distinct during localization.
    pub localization_mods: Vec<(ModificationSpecificity, Arc<ModificationDefinition>)>,
    /// Variable and precursor-label modifications used by property models.
    /// Label modifications are intentionally excluded from `potential_mods`
    /// because they are not PTM-localization candidates.
    pub model_mods: Vec<(ModificationSpecificity, f32)>,
    pub label_reference: Option<Arc<str>>,
    pub label_channels: Vec<Arc<str>>,
    pub generate_decoys: bool,
    pub decoy_tag: String,
    /// Optional explicit target pairing for non-reversal decoy peptides.
    pub decoy_pairing: Vec<PeptideIx>,
    /// Search only the mass-offset hypotheses, not the unmodified precursor.
    /// For open searches where every precursor carries one offset (for
    /// example the innermost HexNAc of an N-glycan): the offset hypothesis
    /// already counts unshifted fragments, so the unmodified pass is
    /// redundant work. Ignored when no mass offsets are configured.
    pub offsets_only: bool,
    /// Peaks each candidate may not claim as b/y ions in full scoring, for
    /// example the Y and oxonium ions of a glycopeptide. Preliminary
    /// retrieval is unaffected. `None` (the default) excludes nothing.
    pub peak_exclusion: Option<Arc<PeakExclusion>>,
}

impl IndexedDatabase {
    /// Find the paired target or decoy for a peptide. Explicit library pairing
    /// takes precedence. Generated FASTA decoys are located by their canonical
    /// reversed peptidoform identity.
    pub fn paired_peptide_index(&self, peptide_index: PeptideIx) -> Option<PeptideIx> {
        if let Some(&paired) = self.decoy_pairing.get(peptide_index.0 as usize) {
            if paired != PeptideIx::default() {
                return Some(paired);
            }
        }
        if !self.generate_decoys {
            return None;
        }

        let peptide = self.peptides.get(peptide_index.0 as usize)?;
        let paired = peptide.reverse();
        let index = self
            .peptides
            .binary_search_by(|candidate| {
                candidate
                    .monoisotopic
                    .total_cmp(&paired.monoisotopic)
                    .then_with(|| candidate.initial_sort(&paired))
            })
            .ok()?;
        (self.peptides[index].decoy != peptide.decoy).then_some(PeptideIx(index as u32))
    }

    /// Create a new [`IndexedQuery`] for a specific
    /// [`ProcessedSpectrum`](crate::spectrum::ProcessedSpectrum)
    ///
    /// All matches returned by the query will be within the specified tolerance
    /// parameters
    pub fn query(
        &self,
        precursor_mass: f32,
        precursor_tol: Tolerance,
        fragment_tol: Tolerance,
    ) -> IndexedQuery<'_> {
        let (precursor_lo, precursor_hi) = precursor_tol.bounds(precursor_mass);

        let (pre_idx_lo, pre_idx_hi) = binary_search_slice(
            &self.peptides,
            |p, bounds| p.monoisotopic.total_cmp(bounds),
            precursor_lo,
            precursor_hi,
        );

        IndexedQuery {
            db: self,
            precursor_mass,
            precursor_tol,
            fragment_tol,
            pre_idx_lo,
            pre_idx_hi,
        }
    }

    pub fn size(&self) -> usize {
        self.fragments.len()
    }

    pub fn buckets(&self) -> &[f32] {
        &self.min_value
    }
}

impl IndexedDatabase {
    /// Candidate placements for one offset on an indexed peptide. Sites that
    /// already carry a static or variable modification are excluded, matching
    /// the database expansion and PTM localization rules.
    pub fn mass_offset_sites(&self, peptide: &Peptide, offset: &MassOffset) -> Vec<Site> {
        let mut sites = Vec::new();
        for specificity in &offset.specificities {
            peptide.compatible_sites(*specificity, &mut sites);
        }
        sites.retain(|site| match site {
            Site::Nterm => peptide.nterm.is_none(),
            Site::Cterm => peptide.cterm.is_none(),
            Site::Sequence(index) => peptide.modification_at(*index as usize) == 0.0,
        });
        if offset.site_mode == SiteMode::Library {
            let supported = self.library_offset_positions(peptide, offset);
            sites.retain(|site| supported.contains(site));
        }
        sites.sort_unstable();
        sites.dedup();
        sites
    }

    /// Peptide-local positions supported by the site library. Decoys use the
    /// mirrored target coordinates of their reversed sequence, so target and
    /// decoy hypotheses receive the same number of library placements.
    fn library_offset_positions(&self, peptide: &Peptide, offset: &MassOffset) -> Vec<Site> {
        let (Some(library), Some(name)) = (&self.offset_library, offset.definition.name.as_deref())
        else {
            return Vec::new();
        };
        let length = peptide.sequence.len() as u32;
        let last = length.saturating_sub(1);
        let target_sequence = if peptide.decoy {
            Cow::Owned(peptide.sequence.reversed_internal())
        } else {
            Cow::Borrowed(&peptide.sequence)
        };
        let mut positions = Vec::new();
        for occurrence in peptide.protein_sites.iter() {
            let Some(start) = occurrence.start else {
                continue;
            };
            for site in library.sites_for(&occurrence.protein) {
                if site.modification.as_ref() != name
                    || !(start..start.saturating_add(length)).contains(&site.position)
                {
                    continue;
                }
                let position = site.position - start;
                if target_sequence.get(position as usize) != Some(&site.residue) {
                    continue;
                }
                let index = if peptide.decoy
                    && (1..last).contains(&position)
                    && site.attachment == crate::ptm_library::Attachment::Residue
                {
                    last - position
                } else {
                    position
                };
                if let Some(site) =
                    site.attachment
                        .site(index, peptide.sequence.len(), peptide.position)
                {
                    positions.push(site)
                }
            }
        }
        positions
    }

    /// Resolve the peptidoform scored for a feature, including an offset
    /// placement that has not yet been materialized.
    pub fn resolve_peptide(&self, feature: &Feature) -> Cow<'_, Peptide> {
        match feature.mass_offset {
            Some(assignment) if (feature.peptide_idx.0 as usize) < self.peptides.len() => {
                Cow::Owned(
                    self.peptides[feature.peptide_idx.0 as usize].with_mass_offset(
                        assignment.site,
                        &self.mass_offsets[assignment.offset as usize].definition,
                    ),
                )
            }
            _ => Cow::Borrowed(&self[feature.peptide_idx]),
        }
    }

    /// Give every offset-placed PSM a stable peptide identity. Placements that
    /// reproduce an indexed peptidoform reuse that peptide; other placements
    /// are stored once in `offset_peptides`, so FDR, LFQ, and site reports see
    /// ordinary peptides. Returns the number of materialized peptidoforms.
    pub fn materialize_mass_offsets(&mut self, features: &mut [Feature]) -> usize {
        let base = self.peptides.len();
        let mut identities: HashMap<(String, bool), PeptideIx> = self
            .offset_peptides
            .iter()
            .enumerate()
            .map(|(index, peptide)| {
                (
                    (peptide.to_string(), peptide.decoy),
                    PeptideIx((base + index) as u32),
                )
            })
            .collect();
        for feature in features.iter_mut() {
            let Some(assignment) = feature.mass_offset else {
                continue;
            };
            if feature.peptide_idx.0 as usize >= base {
                continue;
            }
            let peptide = self.peptides[feature.peptide_idx.0 as usize].with_mass_offset(
                assignment.site,
                &self.mass_offsets[assignment.offset as usize].definition,
            );
            if let Some(indexed) = self.find_indexed(&peptide) {
                feature.peptide_idx = indexed;
                feature.mass_offset = None;
                continue;
            }
            let key = (peptide.to_string(), peptide.decoy);
            feature.peptide_idx = *identities.entry(key).or_insert_with(|| {
                self.offset_peptides.push(peptide);
                PeptideIx((base + self.offset_peptides.len() - 1) as u32)
            });
        }
        self.offset_peptides.len()
    }

    fn find_indexed(&self, peptide: &Peptide) -> Option<PeptideIx> {
        const MASS_EPSILON: f32 = 1e-3;
        let start = self.peptides.partition_point(|candidate| {
            candidate.monoisotopic < peptide.monoisotopic - MASS_EPSILON
        });
        self.peptides[start..]
            .iter()
            .take_while(|candidate| candidate.monoisotopic <= peptide.monoisotopic + MASS_EPSILON)
            .position(|candidate| same_peptidoform(candidate, peptide))
            .map(|offset| PeptideIx((start + offset) as u32))
    }
}

impl std::ops::Index<PeptideIx> for IndexedDatabase {
    type Output = Peptide;

    fn index(&self, index: PeptideIx) -> &Self::Output {
        let index = index.0 as usize;
        match self.peptides.get(index) {
            Some(peptide) => peptide,
            None => &self.offset_peptides[index - self.peptides.len()],
        }
    }
}

impl Parameters {
    pub fn build(self, fasta: Fasta) -> IndexedDatabase {
        self.build_with_custom_cleavages(fasta, None)
    }

    pub fn build_with_custom_cleavages(
        self,
        fasta: Fasta,
        custom_cleavages: Option<&ValidatedCustomCleavageLibrary>,
    ) -> IndexedDatabase {
        let target_decoys = self.digest_with_custom_cleavages(&fasta, custom_cleavages);
        self.build_from_peptides(target_decoys)
    }

    pub fn build_from_peptides(self, target_decoys: Vec<Peptide>) -> IndexedDatabase {
        log::trace!("generating fragments");
        let (compressed_fragments, min_value) = FragmentIndex::build(&self, &target_decoys, None);
        self.assemble_database(target_decoys, compressed_fragments, min_value)
    }

    /// Like [`Parameters::build_from_peptides`], but `extra` adds neutral
    /// fragment masses of each peptide to the preliminary fragment index,
    /// after its regular ions. They only count toward preliminary candidate
    /// retrieval (and, for mass-offset hypotheses, are also looked up with
    /// the offset's fragment shift); full scoring is unchanged. `extra` must
    /// be deterministic: it is called twice per peptide.
    pub fn build_from_peptides_with_extra_fragments(
        self,
        target_decoys: Vec<Peptide>,
        extra: &ExtraFragments,
    ) -> IndexedDatabase {
        log::trace!("generating fragments");
        let (compressed_fragments, min_value) =
            FragmentIndex::build(&self, &target_decoys, Some(extra));
        self.assemble_database(target_decoys, compressed_fragments, min_value)
    }

    /// Build the peptide table and metadata without a fragment index. Used by
    /// the spectrum-indexed prefilter for decoy-pair closure.
    pub fn build_peptide_table(self, target_decoys: Vec<Peptide>) -> IndexedDatabase {
        self.assemble_database(target_decoys, FragmentIndex::default(), Vec::new())
    }

    fn assemble_database(
        self,
        target_decoys: Vec<Peptide>,
        compressed_fragments: FragmentIndex,
        min_value: Vec<f32>,
    ) -> IndexedDatabase {
        // Preserve names for mass-only consumers. Localization also retains
        // full definitions so equal-mass modifications remain distinct.
        for entries in self.variable_mods.values() {
            for entry in entries {
                let definition = entry.definition();
                if let Some(name) = definition.name.as_deref() {
                    if definition.mass.abs() >= 1e-5 {
                        crate::unimod::register_label(definition.mass, name);
                    }
                    for offset in definition.channel_offsets.values() {
                        let mass = definition.mass + offset;
                        if mass.abs() >= 1e-5 {
                            crate::unimod::register_label(mass, name);
                        }
                    }
                }
            }
        }
        for entry in self.static_mods.values() {
            let definition = entry.definition();
            if let Some(name) = definition.name.as_deref() {
                for offset in definition.channel_offsets.values() {
                    let mass = definition.mass + offset;
                    if mass.abs() >= 1e-5 {
                        crate::unimod::register_label(mass, name);
                    }
                }
            }
        }

        let mut localization_mods = self
            .variable_mods
            .iter()
            .flat_map(|(specificity, entries)| {
                entries.iter().flat_map(move |entry| {
                    let definition = entry.definition();
                    let mut definitions = vec![(*specificity, Arc::new(definition.clone()))];
                    definitions.extend(definition.channel_offsets.values().map(|offset| {
                        (
                            *specificity,
                            Arc::new(definition.with_mass(definition.mass + offset)),
                        )
                    }));
                    definitions
                })
            })
            .collect::<Vec<_>>();
        localization_mods.sort_by(|a, b| a.0.cmp(&b.0).then_with(|| a.1.cmp(&b.1)));
        localization_mods.dedup();

        let mut potential_mods = self
            .variable_mods
            .iter()
            .flat_map(|(specificity, entries)| {
                entries.iter().flat_map(move |entry| {
                    let definition = entry.definition();
                    let mut masses = vec![(*specificity, definition.mass)];
                    masses.extend(
                        definition
                            .channel_offsets
                            .values()
                            .map(|offset| (*specificity, definition.mass + offset)),
                    );
                    masses
                })
            })
            .collect::<Vec<(ModificationSpecificity, f32)>>();
        potential_mods.sort_unstable_by(|left, right| {
            left.0
                .cmp(&right.0)
                .then_with(|| left.1.total_cmp(&right.1))
        });
        potential_mods.dedup();
        let mut model_mods = potential_mods.clone();
        model_mods.extend(
            self.static_mods
                .iter()
                .filter(|(_, entry)| !entry.channel_offsets().is_empty())
                .flat_map(|(specificity, entry)| {
                    let definition = entry.definition();
                    definition
                        .channel_offsets
                        .values()
                        .map(|offset| (*specificity, definition.mass + offset))
                        .collect::<Vec<_>>()
                }),
        );
        model_mods.sort_unstable_by(|left, right| {
            left.0
                .cmp(&right.0)
                .then_with(|| left.1.total_cmp(&right.1))
        });
        model_mods.dedup();

        let label_reference = self.label_reference();
        let label_channels = self.label_channels();
        let mass_offsets = self.mass_offset_modifications();
        let offset_library = mass_offsets
            .iter()
            .any(|offset| offset.site_mode == SiteMode::Library)
            .then(|| self.loaded_ptm_library.clone())
            .flatten();
        IndexedDatabase {
            peptides: target_decoys,
            offset_peptides: Vec::new(),
            mass_offsets,
            offset_library,
            fragments: compressed_fragments,
            min_value,
            fragment_losses: self.fragment_losses().map(Arc::new),
            ion_kinds: self.ion_kinds,
            generate_decoys: self.generate_decoys,
            potential_mods,
            localization_mods,
            model_mods,
            label_reference,
            label_channels,
            decoy_tag: self.decoy_tag,
            decoy_pairing: Vec::new(),
            offsets_only: false,
            peak_exclusion: None,
        }
    }
}
