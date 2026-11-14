//! Digestion, modification expansion and decoy generation.
//!
//! Moved out of `database.rs` unchanged; every item keeps its `database::` path.

use crate::ambiguous_residues;
use crate::cleavage::ValidatedCustomCleavageLibrary;
use crate::enzyme::{group_protein_digests, Digest, DigestGroup, Position, ProteinOccurrence};
use crate::fasta::Fasta;
use crate::peptide::{
    LabelModificationCache, LibrarySite, ModificationKind, ModificationLookup, ModificationPlan,
    Peptide,
};
use crate::ptm_library::PtmLibrary;
use crate::sequence::{PeptideSequence, ProteinSequence};
use dashmap::DashSet;
use fnv::FnvBuildHasher;
use rayon::prelude::*;
use std::collections::HashSet;
use std::sync::Arc;

use super::*;

/// Hashes of a peptide sequence and of its internal reversal, the sequence
/// of the decoy generated from it. Equal sequences hash equally; unequal
/// sequences rarely do.
pub fn sequence_hashes(sequence: &[u8]) -> (u64, u64) {
    canonical_sequence_hashes(sequence, false)
}

/// [`sequence_hashes`], with I, L and J hashed as one residue when
/// `merge_isoleucine_leucine` is set, so peptides that
/// [`Parameters::merge_isoleucine_leucine`] merges share their hashes.
pub fn canonical_sequence_hashes(sequence: &[u8], merge_isoleucine_leucine: bool) -> (u64, u64) {
    const OFFSET: u64 = 0xcbf2_9ce4_8422_2325;
    const PRIME: u64 = 0x0000_0100_0000_01b3;
    let step = |hash: u64, byte: &u8| {
        let byte = match merge_isoleucine_leucine {
            true => crate::ambiguous_residues::isoleucine_leucine_canonical(*byte),
            false => *byte,
        };
        (hash ^ u64::from(byte)).wrapping_mul(PRIME)
    };
    let forward = sequence.iter().fold(OFFSET, step);
    let reversed = match sequence {
        [first, middle @ .., last] if !middle.is_empty() => {
            let hash = step(OFFSET, first);
            let hash = middle.iter().rev().fold(hash, step);
            step(hash, last)
        }
        _ => forward,
    };
    // FNV's low bits mix poorly; finish with a SplitMix64 round.
    let mix = |mut hash: u64| {
        hash = (hash ^ (hash >> 30)).wrapping_mul(0xbf58_476d_1ce4_e5b9);
        hash = (hash ^ (hash >> 27)).wrapping_mul(0x94d0_49bb_1331_11eb);
        hash ^ (hash >> 31)
    };
    (mix(forward), mix(reversed))
}

/// Digest groups expanded per parallel pass in
/// [`Parameters::modify_digests_with_target_sequences`].
const MODIFY_DIGEST_CHUNK_GROUPS: usize = 1 << 16;

/// Expands unmodified digest groups into modified target and decoy peptides.
/// See [`Parameters::with_digest_expander`].
pub struct DigestExpander<'a> {
    parameters: &'a Parameters,
    plan: ModificationPlan<'a>,
    label_channels: &'a [Arc<str>],
    label_modifications: &'a LabelModificationCache,
    label_reference: Option<Arc<str>>,
    library: Option<&'a PtmLibrary>,
}

impl DigestExpander<'_> {
    /// Database parameters the expander was built from.
    pub fn parameters(&self) -> &Parameters {
        self.parameters
    }

    /// Sort and deduplicate expanded peptides, as
    /// [`Parameters::modify_digests`] does.
    pub fn reorder(&self, target_decoys: &mut Vec<Peptide>) {
        Parameters::reorder_peptides_with_reference(
            target_decoys,
            self.label_reference.as_deref(),
            self.parameters.merge_isoleucine_leucine,
        );
    }

    /// Modified peptides and generated decoys of one digest group, unsorted.
    /// Decoys whose sequence is in `target_sequences` are dropped; `None`
    /// skips that check for groups known not to collide with any target.
    pub fn expand(
        &self,
        group: DigestGroup,
        target_sequences: Option<&HashSet<PeptideSequence>>,
    ) -> Vec<Peptide> {
        let parameters = self.parameters;
        let expand = |peptide: Peptide, library_sites: &[LibrarySite]| {
            let decoy_sequence = parameters
                .generate_decoys
                .then(|| peptide.sequence.reversed_internal());
            peptide
                .apply_rules(&self.plan, library_sites)
                .into_iter()
                .flat_map(|peptide| {
                    if self.label_channels.is_empty() {
                        vec![peptide]
                    } else {
                        self.label_channels
                            .iter()
                            .map(|channel| {
                                peptide
                                    .clone()
                                    .apply_label_channel(channel.clone(), self.label_modifications)
                            })
                            .collect()
                    }
                })
                .filter(|peptide| {
                    peptide.monoisotopic >= parameters.peptide_min_mass
                        && peptide.monoisotopic <= parameters.peptide_max_mass
                })
                .flat_map(|peptide| {
                    if let Some(sequence) = &decoy_sequence {
                        vec![peptide.reverse_with_sequence(sequence.clone()), peptide].into_iter()
                    } else {
                        vec![peptide].into_iter()
                    }
                })
                .filter(|peptide| {
                    !peptide.decoy
                        || !target_sequences.is_some_and(|targets| {
                            targets.contains(&peptide.sequence)
                                || (parameters.merge_isoleucine_leucine
                                    && targets.contains(
                                        &parameters.decoy_collision_key(&peptide.sequence),
                                    ))
                        })
                })
                .collect::<Vec<_>>()
        };

        match self.library {
            None => Peptide::try_from(group)
                .map(|peptide| expand(peptide, &[]))
                .unwrap_or_default(),
            Some(library) => {
                let reference = group.reference;
                group
                    .origins
                    .into_iter()
                    .flat_map(|origin| {
                        let mut digest = reference.clone();
                        digest.protein = origin.protein.clone();
                        digest.protein_start = origin.start;
                        digest.prev_aa = origin.prev_aa;
                        digest.next_aa = origin.next_aa;
                        let Ok(mut peptide) = Peptide::try_from(digest) else {
                            return Vec::new();
                        };
                        peptide.proteins = smallvec::smallvec![origin.protein.clone()];
                        // The reference sequence views another protein;
                        // keep this origin's own source for motif rules.
                        peptide.protein_sites = Arc::from([origin.clone()]);
                        let start = origin.start.unwrap_or_default();
                        let end = start.saturating_add(peptide.sequence.len() as u32);
                        let library_sites = library
                            .sites_for(&origin.protein)
                            .iter()
                            .filter(|site| (start..end).contains(&site.position))
                            .filter_map(|site| {
                                let position = site.position - start;
                                (peptide.sequence.get(position as usize) == Some(&site.residue))
                                    .then(|| LibrarySite {
                                        attachment: site.attachment,
                                        position,
                                        modification: site.modification.clone(),
                                    })
                            })
                            .collect::<Vec<_>>();
                        expand(peptide, &library_sites)
                    })
                    .collect()
            }
        }
    }
}

fn channel_zero_provenance_eq(left: &Peptide, right: &Peptide) -> bool {
    let retained = |peptide: &Peptide| {
        peptide
            .applied_modifications()
            .filter(|applied| {
                !(applied.kind == ModificationKind::Label && applied.modification.mass == 0.0)
            })
            .map(|applied| (applied.site, applied.modification.clone(), applied.kind))
            .collect::<Vec<_>>()
    };
    retained(left) == retained(right)
}

fn chemical_modifications_eq(left: &Peptide, right: &Peptide) -> bool {
    left.nterm.unwrap_or_default() == right.nterm.unwrap_or_default()
        && left.cterm.unwrap_or_default() == right.cterm.unwrap_or_default()
        && (0..left.sequence.len())
            .all(|index| left.modification_at(index) == right.modification_at(index))
}

fn preferred_channel(
    left: Option<&str>,
    right: Option<&str>,
    reference: Option<&str>,
) -> Option<Arc<str>> {
    if let Some(reference) = reference {
        if left == Some(reference) || right == Some(reference) {
            return Some(Arc::from(reference));
        }
    }
    match (left, right) {
        (Some(left), Some(right)) => Some(Arc::from(left.min(right))),
        (Some(channel), None) | (None, Some(channel)) => Some(Arc::from(channel)),
        (None, None) => None,
    }
}

/// Chemical peptidoform identity used to merge equivalent offset and indexed
/// hypotheses.
pub fn same_peptidoform(left: &Peptide, right: &Peptide) -> bool {
    left.decoy == right.decoy
        && left.sequence == right.sequence
        && left.label_channel() == right.label_channel()
        && chemical_modifications_eq(left, right)
}

impl Parameters {
    /// Digest and group proteins without applying variable modifications.
    pub fn digest_unmodified(&self, fasta: &Fasta) -> Vec<DigestGroup> {
        self.digest_unmodified_with_custom_cleavages(fasta, None)
    }

    pub fn digest_unmodified_with_custom_cleavages(
        &self,
        fasta: &Fasta,
        custom_cleavages: Option<&ValidatedCustomCleavageLibrary>,
    ) -> Vec<DigestGroup> {
        log::trace!("digesting fasta");
        let enzyme = self.enzyme_parameters();
        let started = std::time::Instant::now();
        let digests = fasta.digest_with_custom_cleavages(&enzyme, custom_cleavages);
        log::debug!(
            "database build: cleavage {:#?} ({} digests)",
            started.elapsed(),
            digests.len()
        );

        log::trace!("grouping digests");
        let start_num = digests.len();
        let started = std::time::Instant::now();
        let digests = group_protein_digests(digests);
        log::debug!(
            "database build: grouping {:#?} ({} groups)",
            started.elapsed(),
            digests.len()
        );
        log::trace!(
            "grouped {} digests into {} groups",
            start_num,
            digests.len()
        );
        digests
    }

    /// Key under which target sequences are recorded and generated decoys
    /// looked up, so a decoy equal to a target is dropped. With
    /// `merge_isoleucine_leucine`, I, L and J count as one residue, as a
    /// decoy that is a target's I/L twin would merge with it.
    pub fn decoy_collision_key(&self, sequence: &PeptideSequence) -> PeptideSequence {
        match self.merge_isoleucine_leucine {
            true => sequence
                .as_bytes()
                .iter()
                .map(|&residue| crate::ambiguous_residues::isoleucine_leucine_canonical(residue))
                .collect::<Vec<_>>()
                .into(),
            false => sequence.clone(),
        }
    }

    /// Expand variable modifications and generate decoys from an unmodified digest.
    pub fn modify_digests(&self, digests: Vec<DigestGroup>) -> Vec<Peptide> {
        let started = std::time::Instant::now();
        // One collision key per target group, over a million of them in a
        // human FASTA: hashing them on one thread cost as much as the whole
        // modification expansion.
        let target_sequences = digests
            .par_iter()
            .filter(|digest| !digest.reference.decoy)
            .map(|digest| self.decoy_collision_key(&digest.reference.sequence))
            .collect::<HashSet<_>>();
        log::debug!(
            "database build: target sequences {:#?} ({} sequences)",
            started.elapsed(),
            target_sequences.len()
        );
        self.modify_digests_with_target_sequences(digests, &target_sequences)
    }

    /// Expand a digest chunk while checking generated decoys against every
    /// target sequence in the complete database. Record the targets by
    /// [`Self::decoy_collision_key`] so decoys that are I/L twins of a
    /// target are dropped too.
    pub fn modify_digests_with_target_sequences(
        &self,
        digests: Vec<DigestGroup>,
        target_sequences: &HashSet<PeptideSequence>,
    ) -> Vec<Peptide> {
        self.modify_digest_chunks(digests, target_sequences, MODIFY_DIGEST_CHUNK_GROUPS)
    }

    pub(super) fn modify_digest_chunks(
        &self,
        digests: Vec<DigestGroup>,
        target_sequences: &HashSet<PeptideSequence>,
        chunk_groups: usize,
    ) -> Vec<Peptide> {
        log::trace!("modifying peptides");
        // Expand in fixed-size chunks of digest groups. A single parallel
        // `collect` over every group builds per-thread pieces the size of the
        // whole peptide list; once they are concatenated and freed, the
        // allocator keeps most of those pages, so they sit under the fragment
        // index build and raise peak memory. Chunking bounds the pieces to one
        // chunk, whose pages the next chunk reuses. Chunks are appended in
        // order, so the output order is unchanged. Adapted from
        // theGreatHerrLebert/sage ccce5da (chunked peptide materialisation).
        let mut target_decoys = Vec::new();
        let started = std::time::Instant::now();
        self.with_digest_expander(|expander| {
            let mut digests = digests.into_iter();
            loop {
                let chunk = digests
                    .by_ref()
                    .take(chunk_groups.max(1))
                    .collect::<Vec<_>>();
                if chunk.is_empty() {
                    break;
                }
                target_decoys.par_extend(
                    chunk
                        .into_par_iter()
                        .flat_map_iter(|group| expander.expand(group, Some(target_sequences))),
                );
            }
        });
        log::debug!(
            "database build: modification expansion {:#?} ({} peptides)",
            started.elapsed(),
            target_decoys.len()
        );
        let started = std::time::Instant::now();
        self.reorder_peptides_with_labels(&mut target_decoys);
        log::debug!(
            "database build: sort and deduplicate {:#?} ({} peptides)",
            started.elapsed(),
            target_decoys.len()
        );
        target_decoys
    }

    /// Sort and deduplicate peptides as [`Self::modify_digests`] does,
    /// merging label channels with this database's reference channel.
    pub fn reorder_peptides_with_labels(&self, target_decoys: &mut Vec<Peptide>) {
        Self::reorder_peptides_with_reference(
            target_decoys,
            self.label_reference().as_deref(),
            self.merge_isoleucine_leucine,
        );
    }

    /// Run `f` with a [`DigestExpander`] for this database's modification
    /// rules, so the rules are prepared once for many digest groups.
    pub fn with_digest_expander<R>(&self, f: impl FnOnce(&DigestExpander<'_>) -> R) -> R {
        let mods = self.variable_modifications();
        let static_mods = self.static_modifications();
        let label_channels = self.label_channels();
        let label_modifications = LabelModificationCache::new(
            mods.iter()
                .map(|rule| &rule.modification)
                .chain(static_mods.values()),
            &label_channels,
        );
        let modification_lookup = ModificationLookup::for_rules(
            &mods,
            &static_mods,
            &label_channels,
            &label_modifications,
        )
        .unwrap_or_else(|error| panic!("{error}"));
        let modification_plan = ModificationPlan::new(
            &mods,
            &static_mods,
            modification_lookup,
            self.max_variable_mods,
            self.max_total_variable_mods,
            self.max_combinations,
        );
        f(&DigestExpander {
            parameters: self,
            plan: modification_plan,
            label_channels: &label_channels,
            label_modifications: &label_modifications,
            label_reference: self.label_reference(),
            library: self.loaded_ptm_library.as_deref(),
        })
    }
    pub fn digest(&self, fasta: &Fasta) -> Vec<Peptide> {
        self.digest_with_custom_cleavages(fasta, None)
    }

    pub fn digest_with_custom_cleavages(
        &self,
        fasta: &Fasta,
        custom_cleavages: Option<&ValidatedCustomCleavageLibrary>,
    ) -> Vec<Peptide> {
        self.modify_digests(self.digest_unmodified_with_custom_cleavages(fasta, custom_cleavages))
    }

    pub fn reorder_peptides(target_decoys: &mut Vec<Peptide>) {
        Self::reorder_peptides_with_reference(target_decoys, None, false);
    }

    /// [`Self::reorder_peptides`], merging I/L/J twins when this database's
    /// `merge_isoleucine_leucine` is set.
    pub fn reorder_merged_peptides(&self, target_decoys: &mut Vec<Peptide>) {
        Self::reorder_peptides_with_reference(target_decoys, None, self.merge_isoleucine_leucine);
    }

    /// Add reversed decoys to an already filtered target set.
    ///
    /// This is used by the low-memory prefilter so decoys are only generated
    /// for targets that survive the preliminary search.
    pub fn add_reversed_decoys(&self, targets: Vec<Peptide>) -> Vec<Peptide> {
        let target_sequences: DashSet<PeptideSequence, FnvBuildHasher> = DashSet::default();
        targets
            .iter()
            .filter(|peptide| !peptide.decoy)
            .for_each(|peptide| {
                target_sequences.insert(self.decoy_collision_key(&peptide.sequence));
            });

        let mut target_decoys = targets
            .into_par_iter()
            .flat_map_iter(|peptide| {
                if peptide.decoy {
                    return vec![peptide];
                }
                let decoy = peptide.reverse();
                if target_sequences.contains(&self.decoy_collision_key(&decoy.sequence)) {
                    vec![peptide]
                } else {
                    vec![decoy, peptide]
                }
            })
            .collect::<Vec<_>>();
        Self::reorder_peptides_with_reference(
            &mut target_decoys,
            self.label_reference().as_deref(),
            self.merge_isoleucine_leucine,
        );
        target_decoys
    }

    fn reorder_peptides_with_reference(
        target_decoys: &mut Vec<Peptide>,
        reference: Option<&str>,
        merge_isoleucine_leucine: bool,
    ) {
        use crate::ambiguous_residues::{isoleucine_leucine_canonical, isoleucine_leucine_eq};
        log::trace!("sorting and deduplicating peptides");

        let init_size = target_decoys.len();
        // This is equivalent to a stable sort. The same peptide can come from
        // digests with different enzymatic state, e.g. a protein N-terminal
        // peptide that is semi-enzymatic in another protein. The
        // kept copy is then the most enzymatic one, whatever the input order.
        // A FASTA copy is kept over a peptide TSV copy, whose placeholder
        // state (fully enzymatic, whole protein) says nothing about where the
        // peptide sits in a protein.
        let from_tsv = |peptide: &Peptide| {
            !peptide
                .protein_sites
                .iter()
                .any(|occurrence| occurrence.start.is_some())
        };
        let enzymatic_state = |a: &Peptide, b: &Peptide| {
            from_tsv(a)
                .cmp(&from_tsv(b))
                .then(a.semi_enzymatic.cmp(&b.semi_enzymatic))
                .then(a.missed_cleavages.cmp(&b.missed_cleavages))
                .then(a.position.cmp(&b.position))
        };
        match merge_isoleucine_leucine {
            // I/L/J twins sort next to each other: by sequence with I, L
            // and J as one residue, modifications, decoy flag, then the
            // most enzymatic copy first. An exact copy of a target in a
            // FASTA decoy protein sorts right after the target's twins.
            true => target_decoys.par_sort_unstable_by(|a, b| {
                fn canonical(peptide: &Peptide) -> impl Iterator<Item = u8> + '_ {
                    peptide
                        .sequence
                        .iter()
                        .map(|&residue| isoleucine_leucine_canonical(residue))
                }
                a.monoisotopic
                    .total_cmp(&b.monoisotopic)
                    .then_with(|| canonical(a).cmp(canonical(b)))
                    .then_with(|| a.modification_sort(b))
                    .then(a.decoy.cmp(&b.decoy))
                    .then_with(|| enzymatic_state(a, b))
                    .then_with(|| a.sequence.cmp(&b.sequence))
            }),
            false => target_decoys.par_sort_unstable_by(|a, b| {
                a.monoisotopic
                    .total_cmp(&b.monoisotopic)
                    .then_with(|| a.initial_sort(b))
                    .then_with(|| enzymatic_state(a, b))
            }),
        }
        // Protein occurrence that names a merged twin's displayed sequence.
        let first_occurrence = |peptide: &Peptide| {
            (
                // Peptide TSV rows, even expanded ones, come last.
                peptide
                    .protein_sites
                    .iter()
                    .all(|occurrence| occurrence.start.is_none()),
                peptide.protein_sites.iter().min().cloned(),
                peptide.proteins.iter().min().cloned(),
            )
        };
        target_decoys.dedup_by(|remove, keep| {
            let twins = remove.sequence != keep.sequence
                && merge_isoleucine_leucine
                && remove.decoy == keep.decoy
                && isoleucine_leucine_eq(&remove.sequence, &keep.sequence);
            if remove.monoisotopic == keep.monoisotopic
                && (remove.sequence == keep.sequence || twins)
                && chemical_modifications_eq(remove, keep)
                && (remove.modifications == keep.modifications
                    || channel_zero_provenance_eq(remove, keep))
            {
                if remove.label_channel() != keep.label_channel() {
                    let has_channel_site = keep
                        .applied_modifications()
                        .chain(remove.applied_modifications())
                        .any(|applied| applied.kind == ModificationKind::Label);
                    let channel = has_channel_site
                        .then(|| {
                            preferred_channel(
                                keep.label_channel_str(),
                                remove.label_channel_str(),
                                reference,
                            )
                        })
                        .flatten();
                    keep.set_label_channel(channel);
                }
                // Twins show the sequence of the one found first in
                // protein order; the rest of `keep` is the most enzymatic.
                if twins && first_occurrence(remove) < first_occurrence(keep) {
                    keep.sequence = remove.sequence.clone();
                }
                keep.proteins.extend(remove.proteins.iter().cloned());
                if !remove.protein_sites.is_empty() {
                    let mut sites = keep.protein_sites.to_vec();
                    sites.extend(remove.protein_sites.iter().cloned());
                    sites.sort_unstable();
                    sites.dedup();
                    keep.protein_sites = sites.into();
                }
                // When merging peptides from different Fastas,
                // decoys in one fasta might be targets in another
                keep.decoy &= remove.decoy;
                true
            } else {
                false
            }
        });

        if merge_isoleucine_leucine {
            // Restore the database order, which lookups such as
            // `paired_peptide_index` binary-search: by mass, then sequence as
            // written. Only peptides of equal mass move. A FASTA decoy that
            // copies a target but not its displayed twin sorts after it.
            target_decoys.par_sort_unstable_by(|a, b| {
                a.monoisotopic
                    .total_cmp(&b.monoisotopic)
                    .then_with(|| a.initial_sort(b))
                    .then(a.decoy.cmp(&b.decoy))
            });
        }
        target_decoys.par_iter_mut().for_each(|peptide| {
            // Merged I/L/J twins from one protein list it once.
            peptide.proteins.sort_unstable();
            peptide.proteins.dedup();
        });

        let num_dropped = init_size - target_decoys.len();
        log::trace!(
            "dropped {} t/d pairs, remaining {}",
            num_dropped,
            target_decoys.len(),
        );
    }

    /// Build a `Vec<Peptide>` from a pre-digested TSV file.
    ///
    /// The TSV must have a header row. The `sequence` column is required;
    /// `protein` and `decoy` are optional. Configured static, variable, and
    /// channel-aware modifications are applied.
    /// Decoys are generated by reversal when `self.generate_decoys` is true,
    /// subject to the same deduplication as the normal FASTA digest path.
    /// Rows with B, X or Z are expanded like FASTA digests when
    /// `expand_ambiguous_residues` is on, and skipped with a warning otherwise.
    pub fn peptides_from_tsv(&self, content: &str) -> Vec<Peptide> {
        let mut lines = content.lines().filter(|l| !l.trim().is_empty());

        let header = match lines.next() {
            Some(h) => h,
            None => {
                log::warn!("peptide TSV file is empty");
                return vec![];
            }
        };

        let cols: Vec<&str> = header.split('\t').collect();
        let seq_col = match cols.iter().position(|&c| c == "sequence") {
            Some(i) => i,
            None => {
                log::warn!("peptide TSV is missing required `sequence` column");
                return vec![];
            }
        };
        let protein_col = cols.iter().position(|&c| c == "protein");
        let decoy_col = cols.iter().position(|&c| c == "decoy");

        // Parse all rows into Peptide structs.
        let raw: Vec<Peptide> = lines
            .flat_map(|line| {
                let fields: Vec<&str> = line.split('\t').collect();
                let Some(seq) = fields.get(seq_col).map(|seq| seq.trim().to_string()) else {
                    return Vec::new();
                };
                if seq.is_empty() {
                    return Vec::new();
                }
                let protein: Arc<str> = protein_col
                    .and_then(|i| fields.get(i).map(|s| s.trim()))
                    .filter(|s| !s.is_empty())
                    .unwrap_or(seq.as_str())
                    .into();
                let is_decoy = decoy_col
                    .and_then(|i| fields.get(i))
                    .map(|s| s.trim().eq_ignore_ascii_case("true"))
                    .unwrap_or(false);

                let digest = |sequence: String| Digest {
                    decoy: is_decoy,
                    semi_enzymatic: false,
                    sequence: sequence.into(),
                    protein: protein.clone(),
                    protein_start: None,
                    prev_aa: None,
                    next_aa: None,
                    missed_cleavages: 0,
                    position: Position::Full,
                    expanded_from: None,
                };
                if !(self.expand_ambiguous_residues
                    && ambiguous_residues::is_ambiguous(seq.as_bytes()))
                {
                    return match Peptide::try_from(digest(seq)) {
                        Ok(p) => vec![p],
                        Err(e) => {
                            log::warn!("skipping peptide: {e}");
                            Vec::new()
                        }
                    };
                }
                // Expand B, X and Z like a FASTA digest. Each variant keeps
                // the row as written, without coordinates, for the
                // substitutions column.
                let variants = ambiguous_residues::variant_count(seq.as_bytes());
                if variants > self.max_ambiguous_variants {
                    log::warn!(
                        "skipping peptide {seq}: {variants} variant(s) exceed database.max_ambiguous_variants ({})",
                        self.max_ambiguous_variants
                    );
                    return Vec::new();
                }
                let written: ProteinSequence = seq.as_str().into();
                ambiguous_residues::expand(seq.as_bytes())
                    .into_iter()
                    .filter_map(|variant| {
                        let variant = String::from_utf8(variant).ok()?;
                        match Peptide::try_from(digest(variant)) {
                            Ok(mut peptide) => {
                                peptide.protein_sites = Arc::from([ProteinOccurrence {
                                    protein: protein.clone(),
                                    start: None,
                                    prev_aa: None,
                                    next_aa: None,
                                    source: Some(written.clone()),
                                    met_clipped: false,
                                }]);
                                Some(peptide)
                            }
                            Err(e) => {
                                log::warn!("skipping peptide: {e}");
                                None
                            }
                        }
                    })
                    .collect()
            })
            .collect();

        let variable_mods = self.variable_modifications();
        let static_mods = self.static_modifications();
        let label_channels = self.label_channels();
        let label_reference = self.label_reference();
        let label_modifications = LabelModificationCache::new(
            variable_mods
                .iter()
                .map(|rule| &rule.modification)
                .chain(static_mods.values()),
            &label_channels,
        );
        let modification_lookup = ModificationLookup::for_rules(
            &variable_mods,
            &static_mods,
            &label_channels,
            &label_modifications,
        )
        .unwrap_or_else(|error| panic!("{error}"));
        let modification_plan = ModificationPlan::new(
            &variable_mods,
            &static_mods,
            modification_lookup,
            self.max_variable_mods,
            self.max_total_variable_mods,
            self.max_combinations,
        );
        let raw = raw
            .into_iter()
            .flat_map(|peptide| peptide.apply_rules(&modification_plan, &[]))
            .flat_map(|peptide| {
                if label_channels.is_empty() {
                    vec![peptide]
                } else {
                    label_channels
                        .iter()
                        .map(|channel| {
                            peptide
                                .clone()
                                .apply_label_channel(channel.clone(), &label_modifications)
                        })
                        .collect()
                }
            })
            .filter(|peptide| {
                peptide.monoisotopic >= self.peptide_min_mass
                    && peptide.monoisotopic <= self.peptide_max_mass
            })
            .collect::<Vec<_>>();

        // Build target sequence set for decoy deduplication.
        let targets: DashSet<PeptideSequence, FnvBuildHasher> = DashSet::default();
        raw.iter().filter(|p| !p.decoy).for_each(|p| {
            targets.insert(self.decoy_collision_key(&p.sequence));
        });

        // Emit targets (+ generated decoys) into the final list.
        let mut result: Vec<Peptide> = raw
            .into_iter()
            .flat_map(|peptide| {
                if self.generate_decoys && !peptide.decoy {
                    let rev = peptide.reverse();
                    if !targets.contains(&self.decoy_collision_key(&rev.sequence)) {
                        vec![rev, peptide]
                    } else {
                        vec![peptide]
                    }
                } else {
                    vec![peptide]
                }
            })
            .collect();

        Self::reorder_peptides_with_reference(
            &mut result,
            label_reference.as_deref(),
            self.merge_isoleucine_leucine,
        );
        result
    }
}
