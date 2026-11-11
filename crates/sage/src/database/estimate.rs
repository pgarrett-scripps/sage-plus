//! Memory estimation for a database before it is built.
//!
//! Moved out of `database.rs` unchanged; every item keeps its `database::` path.

use crate::cleavage::ValidatedCustomCleavageLibrary;
use crate::enzyme::{Digest, DigestGroup, ProteinOccurrence};
use crate::fasta::Fasta;
use crate::modification::SiteMode;
use crate::peptide::{AppliedModification, Peptide, Site, VariableRule, INLINE_PROTEINS};
use fnv::FnvHashMap;
use rayon::prelude::*;
use std::collections::HashSet;
use std::sync::Arc;

use super::*;

/// Conservative peak-memory estimates for the major database-build stages.
#[derive(Clone, Copy, Debug, Default)]
pub struct DatabaseMemoryEstimate {
    pub unmodified_peptides: u64,
    pub modified_peptides: u64,
    pub fragments: u64,
    pub unmodified_peak_bytes: u64,
    pub modified_peak_bytes: u64,
    pub fragment_peak_bytes: u64,
}

impl DatabaseMemoryEstimate {
    /// Largest of the three stage peaks.
    pub fn peak_bytes(&self) -> u64 {
        self.unmodified_peak_bytes
            .max(self.modified_peak_bytes)
            .max(self.fragment_peak_bytes)
    }
}

/// Estimated bytes one unmodified digest of `sequence_len` residues holds.
fn digest_bytes(sequence_len: u64) -> u64 {
    const ALLOCATION_OVERHEAD: u64 = 16;
    (std::mem::size_of::<Digest>() as u64)
        .saturating_add(sequence_len)
        .saturating_add(ALLOCATION_OVERHEAD)
}

/// Per-protein sums for [`Parameters::estimate_memory_with_custom_cleavages`].
#[derive(Clone, Copy, Default)]
struct EstimateTotals {
    unmodified_peptides: u64,
    modified_peptides: u64,
    fragments: u64,
    digest_bytes: u64,
    peptide_bytes: u64,
}

impl EstimateTotals {
    fn add(self, other: Self) -> Self {
        Self {
            unmodified_peptides: self
                .unmodified_peptides
                .saturating_add(other.unmodified_peptides),
            modified_peptides: self
                .modified_peptides
                .saturating_add(other.modified_peptides),
            fragments: self.fragments.saturating_add(other.fragments),
            digest_bytes: self.digest_bytes.saturating_add(other.digest_bytes),
            peptide_bytes: self.peptide_bytes.saturating_add(other.peptide_bytes),
        }
    }
}

fn with_estimation_margin(bytes: u64) -> u64 {
    // Parallel collection, allocator size classes, and sorting create overhead that
    // cannot be derived exactly from item counts. Use a conservative 50% margin.
    bytes.saturating_add(bytes / 2)
}

/// Per-peptide variant count above which the preflight estimate switches from
/// exact counting to [`variant_upper_bound`].
const EXACT_VARIANT_COUNT_LIMIT: u64 = 1 << 20;

/// Count-only mirror of `ModificationEnumeration`. Candidates are
/// `(site, count_group, library_supported)` sorted by site, with at most one
/// entry per site and group.
struct VariantCounter<'a> {
    candidates: &'a [(usize, usize, bool)],
    /// First candidate index at a later site than each candidate.
    next_site: Vec<usize>,
    /// `(max_count, total limit)` per count group.
    limits: &'a [(Option<usize>, Option<usize>)],
    new_counts: Vec<usize>,
    total_counts: Vec<usize>,
    max_exhaustive: usize,
    max_total: usize,
    limit: u64,
    variants: u64,
}

impl<'a> VariantCounter<'a> {
    fn new(
        candidates: &'a [(usize, usize, bool)],
        limits: &'a [(Option<usize>, Option<usize>)],
        max_exhaustive: usize,
        max_total: usize,
        limit: u64,
    ) -> Self {
        let mut next_site = vec![candidates.len(); candidates.len()];
        for idx in (0..candidates.len().saturating_sub(1)).rev() {
            next_site[idx] = if candidates[idx + 1].0 != candidates[idx].0 {
                idx + 1
            } else {
                next_site[idx + 1]
            };
        }
        Self {
            candidates,
            next_site,
            limits,
            new_counts: vec![0; limits.len()],
            total_counts: vec![0; limits.len()],
            max_exhaustive,
            max_total,
            limit,
            // The unmodified peptide.
            variants: 1,
        }
    }

    fn count(&mut self, start: usize, total: usize, exhaustive: usize) {
        if total == self.max_total {
            return;
        }
        for idx in start..self.candidates.len() {
            if self.variants >= self.limit {
                return;
            }
            let (_, group, library) = self.candidates[idx];
            let is_new = !library;
            let (max_count, total_limit) = self.limits[group];
            if (is_new && exhaustive == self.max_exhaustive)
                || (is_new && max_count.is_some_and(|limit| self.new_counts[group] >= limit))
                || total_limit.is_some_and(|limit| self.total_counts[group] >= limit)
            {
                continue;
            }
            self.new_counts[group] += usize::from(is_new);
            self.total_counts[group] += 1;
            self.variants += 1;
            self.count(
                self.next_site[idx],
                total + 1,
                exhaustive + usize::from(is_new),
            );
            self.new_counts[group] -= usize::from(is_new);
            self.total_counts[group] -= 1;
        }
    }
}

/// Upper bound on variants that ignores per-modification caps.
fn variant_upper_bound(
    candidates: &[(usize, usize, bool)],
    max_variable_mods: usize,
    max_total_variable_mods: usize,
) -> u64 {
    let mut choices_by_site: FnvHashMap<usize, (u64, u64)> = FnvHashMap::default();
    for (site, _, library_supported) in candidates.iter().copied() {
        let choices = choices_by_site.entry(site).or_default();
        if library_supported {
            choices.0 += 1;
        } else {
            choices.1 += 1;
        }
    }

    let max_total = max_total_variable_mods.min(choices_by_site.len());
    let max_exhaustive = max_variable_mods.min(max_total);
    let mut counts = vec![vec![0u64; max_exhaustive + 1]; max_total + 1];
    counts[0][0] = 1;
    for (library_choices, exhaustive_choices) in choices_by_site.values().copied() {
        let mut next = counts.clone();
        for total in 0..max_total {
            for exhaustive in 0..=max_exhaustive {
                let current = counts[total][exhaustive];
                next[total + 1][exhaustive] = next[total + 1][exhaustive]
                    .saturating_add(current.saturating_mul(library_choices));
                if exhaustive < max_exhaustive {
                    next[total + 1][exhaustive + 1] = next[total + 1][exhaustive + 1]
                        .saturating_add(current.saturating_mul(exhaustive_choices));
                }
            }
        }
        counts = next;
    }
    counts.into_iter().flatten().fold(0u64, u64::saturating_add)
}

impl Parameters {
    /// Estimate database expansion without retaining digests or modified peptides.
    ///
    /// Counts raw enzymatic digests rather than assuming deduplication or mass filtering,
    /// making this a conservative upper bound for rejecting unsafe searches before the
    /// variable-modification expansion begins.
    pub fn estimate_memory(&self, fasta: &Fasta) -> DatabaseMemoryEstimate {
        self.estimate_memory_with_custom_cleavages(fasta, None)
    }

    pub fn estimate_memory_with_custom_cleavages(
        &self,
        fasta: &Fasta,
        custom_cleavages: Option<&ValidatedCustomCleavageLibrary>,
    ) -> DatabaseMemoryEstimate {
        const ALLOCATION_OVERHEAD: u64 = 16;

        let enzyme = self.enzyme_parameters();
        let decoy_multiplier = if self.generate_decoys { 2 } else { 1 };
        let rules = self.variable_modifications();

        // Proteins are estimated in parallel. Integer sums do not depend on
        // the reduction order.
        let totals = fasta
            .targets
            .par_iter()
            .map(|(protein, sequence)| {
                let boundaries = custom_cleavages
                    .map(|library| library.boundaries_for(protein))
                    .unwrap_or_default();
                let mut totals = EstimateTotals::default();
                for digest in
                    enzyme.digest_with_custom_cleavages(sequence, protein.clone(), boundaries)
                {
                    let sequence_len = digest.sequence.len() as u64;
                    let origin = ProteinOccurrence::of_protein_digest(&digest);
                    let variants = self
                        .variable_variant_count_with(&rules, &digest, std::slice::from_ref(&origin))
                        .saturating_mul(decoy_multiplier);
                    let fragments_per_variant = sequence_len
                        .saturating_sub(1)
                        .saturating_sub(self.min_ion_index as u64)
                        .saturating_mul(self.ion_kinds.len() as u64);

                    // Peptide clones share some Arc allocations, but charging sequence and
                    // protein-reference storage to every variant keeps the estimate safely high.
                    let bytes_per_variant = (std::mem::size_of::<Peptide>() as u64)
                        .saturating_add(sequence_len)
                        .saturating_add(
                            sequence_len.saturating_mul(std::mem::size_of::<f32>() as u64),
                        )
                        .saturating_add(std::mem::size_of::<Arc<str>>() as u64)
                        .saturating_add(ALLOCATION_OVERHEAD.saturating_mul(3));

                    totals = totals.add(EstimateTotals {
                        unmodified_peptides: 1,
                        modified_peptides: variants,
                        fragments: variants.saturating_mul(fragments_per_variant),
                        digest_bytes: digest_bytes(sequence_len),
                        peptide_bytes: variants.saturating_mul(bytes_per_variant),
                    });
                }
                totals
            })
            .reduce(EstimateTotals::default, EstimateTotals::add);

        let mut estimate = DatabaseMemoryEstimate {
            unmodified_peptides: totals.unmodified_peptides,
            modified_peptides: totals.modified_peptides,
            fragments: totals.fragments,
            ..DatabaseMemoryEstimate::default()
        };
        let digest_bytes = totals.digest_bytes;
        let peptide_bytes = totals.peptide_bytes;

        let fragment_bytes = estimate
            .fragments
            .saturating_mul(std::mem::size_of::<PackedFragment>() as u64);
        let bucket_bytes = self.estimated_bucket_bytes(estimate.fragments);

        estimate.unmodified_peak_bytes = with_estimation_margin(digest_bytes.saturating_mul(2));
        estimate.modified_peak_bytes =
            with_estimation_margin(digest_bytes.saturating_add(peptide_bytes));
        estimate.fragment_peak_bytes = with_estimation_margin(
            peptide_bytes
                .saturating_add(fragment_bytes)
                .saturating_add(bucket_bytes),
        );
        estimate
    }

    /// Re-estimate the fragment/index stage from the peptides that actually survived
    /// modification, filtering, prefiltering, and deduplication.
    pub fn estimate_index_memory(&self, peptides: &[Peptide]) -> DatabaseMemoryEstimate {
        const ALLOCATION_OVERHEAD: u64 = 16;

        let mut estimate = DatabaseMemoryEstimate {
            modified_peptides: peptides.len() as u64,
            ..DatabaseMemoryEstimate::default()
        };
        let mut peptide_bytes = 0u64;
        for peptide in peptides {
            let sequence_len = peptide.sequence.len() as u64;
            let fragments = sequence_len
                .saturating_sub(1)
                .saturating_sub(self.min_ion_index as u64)
                .saturating_mul(self.ion_kinds.len() as u64);
            estimate.fragments = estimate.fragments.saturating_add(fragments);
            let protein_bytes = if peptide.proteins.spilled() {
                (peptide.proteins.len() as u64)
                    .saturating_mul(std::mem::size_of::<Arc<str>>() as u64)
                    .saturating_add(ALLOCATION_OVERHEAD)
            } else {
                0
            };
            peptide_bytes =
                peptide_bytes.saturating_add(
                    (std::mem::size_of::<Peptide>() as u64)
                        .saturating_add(sequence_len)
                        .saturating_add(peptide.modifications.heap_bytes() as u64)
                        .saturating_add(protein_bytes)
                        .saturating_add(
                            (peptide.protein_sites.len() as u64)
                                .saturating_mul(std::mem::size_of::<ProteinOccurrence>() as u64),
                        )
                        .saturating_add(ALLOCATION_OVERHEAD.saturating_mul(3)),
                );
        }

        let fragment_bytes = estimate
            .fragments
            .saturating_mul(std::mem::size_of::<PackedFragment>() as u64);
        estimate.modified_peak_bytes = with_estimation_margin(peptide_bytes);
        estimate.fragment_peak_bytes = with_estimation_margin(
            peptide_bytes
                .saturating_add(fragment_bytes)
                .saturating_add(self.estimated_bucket_bytes(estimate.fragments)),
        );
        estimate
    }

    /// Estimate modification expansion from the deduplicated, unmodified digest.
    pub fn estimate_modified_memory(&self, digests: &[DigestGroup]) -> DatabaseMemoryEstimate {
        const ALLOCATION_OVERHEAD: u64 = 16;

        let decoy_multiplier = if self.generate_decoys { 2 } else { 1 };
        let mut estimate = DatabaseMemoryEstimate {
            unmodified_peptides: digests.len() as u64,
            ..DatabaseMemoryEstimate::default()
        };
        let mut peptide_bytes = 0u64;
        for digest in digests {
            let sequence_len = digest.reference.sequence.len() as u64;
            let variants = if self.loaded_ptm_library.is_some() {
                digest
                    .origins
                    .iter()
                    .map(|origin| {
                        let mut reference = digest.reference.clone();
                        reference.protein = origin.protein.clone();
                        reference.protein_start = origin.start;
                        reference.prev_aa = origin.prev_aa;
                        reference.next_aa = origin.next_aa;
                        self.variable_variant_count(&reference, std::slice::from_ref(origin))
                    })
                    .fold(0u64, u64::saturating_add)
            } else {
                self.variable_variant_count(&digest.reference, &digest.origins)
            }
            .saturating_mul(decoy_multiplier);
            estimate.modified_peptides = estimate.modified_peptides.saturating_add(variants);
            estimate.fragments = estimate.fragments.saturating_add(
                variants.saturating_mul(
                    sequence_len
                        .saturating_sub(1)
                        .saturating_sub(self.min_ion_index as u64)
                        .saturating_mul(self.ion_kinds.len() as u64),
                ),
            );

            let bytes_per_variant = (std::mem::size_of::<Peptide>() as u64)
                .saturating_add(sequence_len)
                .saturating_add(sequence_len.saturating_mul(std::mem::size_of::<f32>() as u64))
                .saturating_add(
                    sequence_len.saturating_mul(std::mem::size_of::<AppliedModification>() as u64),
                )
                .saturating_add(if digest.origins.len() > INLINE_PROTEINS {
                    (digest.origins.len() as u64)
                        .saturating_mul(std::mem::size_of::<Arc<str>>() as u64)
                        .saturating_add(ALLOCATION_OVERHEAD)
                } else {
                    0
                })
                .saturating_add(
                    (digest.origins.len() as u64)
                        .saturating_mul(std::mem::size_of::<ProteinOccurrence>() as u64),
                )
                .saturating_add(ALLOCATION_OVERHEAD.saturating_mul(4));
            peptide_bytes =
                peptide_bytes.saturating_add(variants.saturating_mul(bytes_per_variant));
        }
        estimate.modified_peak_bytes = with_estimation_margin(peptide_bytes);
        estimate
    }

    fn estimated_bucket_bytes(&self, fragments: u64) -> u64 {
        let maximum_prefixes = 1u64 << (u32::BITS - FRAGMENT_MASS_SUFFIX_BITS);
        let split_buckets = fragments.div_ceil(self.bucket_size.max(1) as u64);
        let maximum_buckets = fragments.min(maximum_prefixes.saturating_add(split_buckets));
        maximum_buckets.saturating_mul(
            (std::mem::size_of::<FragmentBucket>() + std::mem::size_of::<f32>()) as u64,
        )
    }

    pub(super) fn variable_variant_count(
        &self,
        digest: &Digest,
        origins: &[ProteinOccurrence],
    ) -> u64 {
        self.variable_variant_count_with(&self.variable_modifications(), digest, origins)
    }

    fn variable_variant_count_with(
        &self,
        rules: &[VariableRule],
        digest: &Digest,
        origins: &[ProteinOccurrence],
    ) -> u64 {
        let sequence = digest.sequence.as_bytes();
        let library_sites = self
            .loaded_ptm_library
            .as_deref()
            .map(|library| {
                let start = digest.protein_start.unwrap_or_default();
                let end = start.saturating_add(sequence.len() as u32);
                library
                    .sites_for(&digest.protein)
                    .iter()
                    .filter(|site| (start..end).contains(&site.position))
                    .filter(|site| {
                        sequence.get((site.position - start) as usize) == Some(&site.residue)
                    })
                    .filter_map(|site| {
                        site.attachment
                            .site(site.position - start, sequence.len(), digest.position)
                            .map(|attachment| (attachment, site.modification.clone()))
                    })
                    .collect::<HashSet<_>>()
            })
            .unwrap_or_default();
        let mut candidates: FnvHashMap<(usize, usize), bool> = FnvHashMap::default();
        let nterm = sequence.len();
        let cterm = sequence.len().saturating_add(1);

        for rule in rules {
            let mut add_site = |site: usize| {
                let library_site = if site == nterm {
                    Site::Nterm
                } else if site == cterm {
                    Site::Cterm
                } else {
                    Site::Sequence(site as u32)
                };
                let supported = rule.site_mode != SiteMode::Exhaustive
                    && rule.modification.name.as_ref().is_some_and(|name| {
                        !sequence.is_empty()
                            && library_sites.contains(&(library_site, name.clone()))
                    });
                if rule.site_mode != SiteMode::Library || supported {
                    candidates
                        .entry((site, rule.count_group))
                        .and_modify(|library| *library |= supported)
                        .or_insert(supported);
                }
            };

            for site in rule.specificity.sites_for_occurrences(
                sequence,
                digest.position,
                digest.decoy,
                origins,
            ) {
                add_site(match site {
                    Site::Nterm => nterm,
                    Site::Cterm => cterm,
                    Site::Sequence(i) => i as usize,
                });
            }
        }

        let group_count = rules
            .iter()
            .map(|rule| rule.count_group + 1)
            .max()
            .unwrap_or(0);
        let mut limits = vec![(None, None); group_count];
        for rule in rules {
            limits[rule.count_group] = (rule.max_count, rule.total_limit());
        }
        let mut candidates = candidates
            .into_iter()
            .map(|((site, group), library)| (site, group, library))
            .collect::<Vec<_>>();
        candidates.sort_unstable();

        // Count exactly what `ModificationEnumeration` emits. Stop early at
        // `max_combinations`, or fall back to a cap-free upper bound when a
        // single peptide has too many variants to count cheaply.
        let limit = self
            .max_combinations
            .map_or(EXACT_VARIANT_COUNT_LIMIT, |cap| {
                (cap as u64).min(EXACT_VARIANT_COUNT_LIMIT)
            });
        let mut counter = VariantCounter::new(
            &candidates,
            &limits,
            self.max_variable_mods,
            self.max_total_variable_mods,
            limit,
        );
        counter.count(0, 0, 0);
        let variable_variants = if counter.variants < limit
            || self
                .max_combinations
                .is_some_and(|cap| cap as u64 <= EXACT_VARIANT_COUNT_LIMIT)
        {
            counter.variants
        } else {
            let bound = variant_upper_bound(
                &candidates,
                self.max_variable_mods,
                self.max_total_variable_mods,
            );
            self.max_combinations
                .map_or(bound, |cap| bound.min(cap as u64))
                .max(limit)
        };
        variable_variants.saturating_mul(self.label_channels().len().max(1) as u64)
    }
}
