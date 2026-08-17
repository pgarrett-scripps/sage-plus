use crate::bitmap::BitmapIndex;
use crate::cleavage::ValidatedCustomCleavageLibrary;
use crate::enzyme::{group_digests, Digest, DigestGroup, Enzyme, EnzymeParameters, Position};
use crate::fasta::Fasta;
use crate::ion_series::{IonGroupSeries, Kind};
use crate::mass::Tolerance;
use crate::modification::{
    validate_mods, validate_var_mods, ModificationDefinition, ModificationSpecificity,
    StaticModEntry, VarModEntry,
};
use crate::peptide::{AppliedModification, Peptide};
use dashmap::DashSet;
use fnv::FnvBuildHasher;
use rayon::prelude::*;
use serde::{Deserialize, Serialize};
use std::cmp::Ordering;
use std::collections::HashMap;
use std::hash::Hash;
use std::sync::Arc;

#[derive(Deserialize, Serialize, Clone, Debug)]
pub struct EnzymeBuilder {
    /// How many missed cleavages to use
    pub missed_cleavages: Option<u8>,
    /// Minimum peptide length that will be fragmented
    pub min_len: Option<usize>,
    /// Maximum peptide length that will be fragmented
    pub max_len: Option<usize>,
    pub cleave_at: Option<String>,
    pub restrict: Option<String>,
    pub c_terminal: Option<bool>,
    pub semi_enzymatic: Option<bool>,
}

impl Default for EnzymeBuilder {
    fn default() -> Self {
        Self {
            missed_cleavages: Some(0),
            min_len: Some(5),
            max_len: Some(50),
            cleave_at: Some("KR".into()),
            restrict: Some("P".into()),
            c_terminal: Some(true),
            semi_enzymatic: Some(false),
        }
    }
}

impl From<EnzymeBuilder> for EnzymeParameters {
    fn from(en: EnzymeBuilder) -> EnzymeParameters {
        EnzymeParameters {
            missed_cleavages: en.missed_cleavages.unwrap_or(1),
            min_len: en.min_len.unwrap_or(5),
            max_len: en.max_len.unwrap_or(50),
            enzyme: Enzyme::new(
                &en.cleave_at.unwrap_or_else(|| "KR".into()),
                &en.restrict.unwrap_or_else(|| "".into()),
                en.c_terminal.unwrap_or(true),
                en.semi_enzymatic.unwrap_or(false),
            ),
        }
    }
}

#[derive(Deserialize, Default)]
/// Parameters used for generating the fragment database
pub struct Builder {
    /// This parameter allows tuning of the internal search structure
    pub bucket_size: Option<usize>,
    /// Number of u64 words per peptide in the bitmap index (default 30 → 1920 bins).
    /// Only used when `use_bitmap` is enabled on the scorer.
    pub bitmap_size: Option<usize>,

    pub enzyme: Option<EnzymeBuilder>,
    /// Minimum peptide monoisotopic mass that will be fragmented
    pub peptide_min_mass: Option<f32>,
    /// Maximum peptide monoisotopic mass that will be fragmented
    pub peptide_max_mass: Option<f32>,
    /// Which kind of fragment ions to generate (a, b, c, x, y, z)
    pub ion_kinds: Option<Vec<Kind>>,
    /// Minimum ion index to be generated: 1 will remove b1/y1 ions
    /// 2 will remove b1/b2/y1/y2 ions, etc
    pub min_ion_index: Option<usize>,
    /// Static modifications to add to matching amino acids. Entries may use
    /// the existing bare mass or a structured modification object.
    pub static_mods: Option<HashMap<String, StaticModEntry>>,
    /// Variable modifications to add to matching amino acids.
    /// Each entry is either a bare mass (`15.9949`) or a structured object with
    /// mass, limits, display name, and optional neutral-loss behavior.
    pub variable_mods: Option<HashMap<String, Vec<VarModEntry>>>,
    /// Limit number of variable modifications on a peptide
    pub max_variable_mods: Option<usize>,
    /// Hard cap on the total peptide variants generated per input peptide,
    /// including its unmodified form. Values below 1 are normalized to 1.
    /// Variants with fewer PTMs are preferred (generated first).
    pub max_combinations: Option<usize>,
    /// Use this prefix for decoy proteins
    pub decoy_tag: Option<String>,

    pub generate_decoys: Option<bool>,
    /// Path to fasta database
    pub fasta: Option<String>,
    /// Path to a pre-digested peptide TSV file (additive with `fasta`).
    /// Required column: `sequence`. Optional columns: `protein`, `decoy`.
    /// No variable/static mods are applied; sequences are used as-is.
    pub peptides: Option<String>,
    /// Path to a protein-specific custom cleavage-site TSV or Parquet file.
    /// Required columns: `protein`, `position`; optional column: `context`.
    pub custom_cleavage_sites: Option<String>,
    /// Number of sequences to handle simultaneously when pre-filtering the db
    pub prefilter_chunk_size: Option<usize>,
    /// Pre-filter the database to minimize memory usage
    pub prefilter: Option<bool>,
    /// Pre-filter the database with a minimal amount of memory at the cost of speed
    pub prefilter_low_memory: Option<bool>,
}

impl Builder {
    pub fn make_parameters(self) -> Parameters {
        let bucket_size = self.bucket_size.unwrap_or(8192).next_power_of_two();
        Parameters {
            bucket_size,
            bitmap_size: self.bitmap_size.unwrap_or(30),
            peptide_min_mass: self.peptide_min_mass.unwrap_or(500.0),
            peptide_max_mass: self.peptide_max_mass.unwrap_or(5000.0),
            ion_kinds: self.ion_kinds.unwrap_or(vec![Kind::B, Kind::Y]),
            min_ion_index: self.min_ion_index.unwrap_or(2),
            decoy_tag: self.decoy_tag.unwrap_or_else(|| "rev_".into()),
            enzyme: self.enzyme.unwrap_or_default(),
            static_mods: validate_mods(self.static_mods),
            variable_mods: validate_var_mods(self.variable_mods),
            max_variable_mods: self.max_variable_mods.map(|x| x.max(1)).unwrap_or(2),
            max_combinations: self.max_combinations.map(|x| x.max(1)),
            generate_decoys: self.generate_decoys.unwrap_or(true),
            fasta: self.fasta.unwrap_or_default(),
            peptides: self.peptides,
            custom_cleavage_sites: self.custom_cleavage_sites,
            prefilter_chunk_size: self.prefilter_chunk_size.unwrap_or(0),
            prefilter: self.prefilter.unwrap_or(false),
            prefilter_low_memory: self.prefilter_low_memory.unwrap_or(true),
            use_bitmap: false,
        }
    }

    pub fn update_fasta(&mut self, fasta: String) {
        self.fasta = Some(fasta)
    }
}

#[derive(Serialize, Clone, Debug)]
pub struct Parameters {
    pub bucket_size: usize,
    pub bitmap_size: usize,
    pub enzyme: EnzymeBuilder,
    pub peptide_min_mass: f32,
    pub peptide_max_mass: f32,
    pub ion_kinds: Vec<Kind>,
    pub min_ion_index: usize,
    pub static_mods: HashMap<ModificationSpecificity, StaticModEntry>,
    pub variable_mods: HashMap<ModificationSpecificity, Vec<VarModEntry>>,
    pub max_variable_mods: usize,
    pub max_combinations: Option<usize>,
    pub decoy_tag: String,
    pub generate_decoys: bool,
    pub fasta: String,
    pub peptides: Option<String>,
    pub custom_cleavage_sites: Option<String>,
    pub prefilter_chunk_size: usize,
    pub prefilter: bool,
    pub prefilter_low_memory: bool,
    /// Select the bitmap preliminary-search index instead of the fragment index.
    /// The CLI copies its top-level `use_bitmap` option here before construction.
    #[serde(skip)]
    pub use_bitmap: bool,
}

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

impl Parameters {
    /// Flatten variable modifications into a stable order. This matters when
    /// `max_combinations` truncates variants: equivalent configurations must
    /// retain the same variants regardless of randomized `HashMap` iteration.
    fn variable_modifications(
        &self,
    ) -> Vec<(
        ModificationSpecificity,
        Arc<ModificationDefinition>,
        Option<usize>,
    )> {
        let mut mods = self
            .variable_mods
            .iter()
            .flat_map(|(specificity, entries)| {
                entries.iter().enumerate().map(|(entry_order, entry)| {
                    (
                        *specificity,
                        entry_order,
                        Arc::new(entry.definition()),
                        entry.max_count(),
                    )
                })
            })
            .collect::<Vec<_>>();
        mods.sort_unstable_by(|a, b| a.0.cmp(&b.0).then_with(|| a.1.cmp(&b.1)));
        mods.into_iter()
            .map(|(specificity, _, modification, max_count)| (specificity, modification, max_count))
            .collect()
    }

    fn static_modifications(
        &self,
    ) -> HashMap<ModificationSpecificity, Arc<ModificationDefinition>> {
        self.static_mods
            .iter()
            .map(|(specificity, entry)| (*specificity, Arc::new(entry.definition())))
            .collect()
    }

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

        let enzyme: EnzymeParameters = self.enzyme.clone().into();
        let decoy_multiplier = if self.generate_decoys { 2 } else { 1 };
        let mut estimate = DatabaseMemoryEstimate::default();
        let mut digest_bytes = 0u64;
        let mut peptide_bytes = 0u64;

        for (protein, sequence) in &fasta.targets {
            let boundaries = custom_cleavages
                .map(|library| library.boundaries_for(protein))
                .unwrap_or_default();
            for digest in enzyme.digest_with_custom_cleavages(sequence, protein.clone(), boundaries)
            {
                let sequence_len = digest.sequence.len() as u64;
                let variants = self
                    .variable_variant_count(&digest)
                    .saturating_mul(decoy_multiplier);
                let fragments_per_variant = sequence_len
                    .saturating_sub(1)
                    .saturating_sub(self.min_ion_index as u64)
                    .saturating_mul(self.ion_kinds.len() as u64);

                estimate.unmodified_peptides = estimate.unmodified_peptides.saturating_add(1);
                estimate.modified_peptides = estimate.modified_peptides.saturating_add(variants);
                estimate.fragments = estimate
                    .fragments
                    .saturating_add(variants.saturating_mul(fragments_per_variant));

                digest_bytes = digest_bytes.saturating_add(
                    (std::mem::size_of::<Digest>() as u64)
                        .saturating_add(sequence_len)
                        .saturating_add(ALLOCATION_OVERHEAD),
                );

                // Peptide clones share some Arc allocations, but charging sequence and
                // protein-reference storage to every variant keeps the estimate safely high.
                let bytes_per_variant = (std::mem::size_of::<Peptide>() as u64)
                    .saturating_add(sequence_len)
                    .saturating_add(sequence_len.saturating_mul(std::mem::size_of::<f32>() as u64))
                    .saturating_add(std::mem::size_of::<Arc<str>>() as u64)
                    .saturating_add(ALLOCATION_OVERHEAD.saturating_mul(3));
                peptide_bytes =
                    peptide_bytes.saturating_add(variants.saturating_mul(bytes_per_variant));
            }
        }

        let fragment_bytes = if self.use_bitmap {
            0
        } else {
            estimate
                .fragments
                .saturating_mul(std::mem::size_of::<Theoretical>() as u64)
        };
        let bucket_bytes = self.estimated_bucket_bytes(estimate.fragments);
        let bitmap_bytes = self.estimated_bitmap_peak_bytes(estimate.modified_peptides);

        estimate.unmodified_peak_bytes = with_estimation_margin(digest_bytes.saturating_mul(2));
        estimate.modified_peak_bytes =
            with_estimation_margin(digest_bytes.saturating_add(peptide_bytes));
        estimate.fragment_peak_bytes = with_estimation_margin(
            peptide_bytes
                .saturating_add(fragment_bytes)
                .saturating_add(bucket_bytes)
                .saturating_add(bitmap_bytes),
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
            peptide_bytes = peptide_bytes.saturating_add(
                (std::mem::size_of::<Peptide>() as u64)
                    .saturating_add(sequence_len)
                    .saturating_add(
                        (peptide.modifications.len() as u64)
                            .saturating_mul(std::mem::size_of::<f32>() as u64),
                    )
                    .saturating_add(
                        (peptide.proteins.len() as u64)
                            .saturating_mul(std::mem::size_of::<Arc<str>>() as u64),
                    )
                    .saturating_add(ALLOCATION_OVERHEAD.saturating_mul(3)),
            );
        }

        let fragment_bytes = if self.use_bitmap {
            0
        } else {
            estimate
                .fragments
                .saturating_mul(std::mem::size_of::<Theoretical>() as u64)
        };
        estimate.modified_peak_bytes = with_estimation_margin(peptide_bytes);
        estimate.fragment_peak_bytes = with_estimation_margin(
            peptide_bytes
                .saturating_add(fragment_bytes)
                .saturating_add(self.estimated_bucket_bytes(estimate.fragments))
                .saturating_add(self.estimated_bitmap_peak_bytes(estimate.modified_peptides)),
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
            let variants = self
                .variable_variant_count(&digest.reference)
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
                .saturating_add(
                    (digest.proteins.len() as u64)
                        .saturating_mul(std::mem::size_of::<Arc<str>>() as u64),
                )
                .saturating_add(ALLOCATION_OVERHEAD.saturating_mul(4));
            peptide_bytes =
                peptide_bytes.saturating_add(variants.saturating_mul(bytes_per_variant));
        }
        estimate.modified_peak_bytes = with_estimation_margin(peptide_bytes);
        estimate
    }

    fn estimated_bucket_bytes(&self, fragments: u64) -> u64 {
        if self.use_bitmap {
            return 0;
        }
        let bucket_size = self.bucket_size.max(1) as u64;
        (fragments.saturating_add(bucket_size - 1) / bucket_size)
            .saturating_mul(std::mem::size_of::<f32>() as u64)
    }

    fn estimated_bitmap_peak_bytes(&self, peptides: u64) -> u64 {
        if !self.use_bitmap {
            return 0;
        }
        // BitmapIndex::build currently creates per-peptide temporary Vecs before
        // flattening them into the retained arrays. Include both to bound peak RSS.
        let bitmap_words = (self.bitmap_size as u64).saturating_mul(2);
        let retained = bitmap_words
            .saturating_mul(std::mem::size_of::<u64>() as u64)
            .saturating_add(std::mem::size_of::<f32>() as u64)
            .saturating_add(std::mem::size_of::<PeptideIx>() as u64);
        let temporary = bitmap_words
            .saturating_mul(std::mem::size_of::<u64>() as u64)
            .saturating_add(std::mem::size_of::<(Vec<u64>, Vec<u64>)>() as u64);
        peptides.saturating_mul(retained.saturating_add(temporary))
    }

    fn variable_variant_count(&self, digest: &Digest) -> u64 {
        let sequence = digest.sequence.as_bytes();
        let mut choices_by_site: HashMap<usize, u64> = HashMap::new();
        let nterm = sequence.len();
        let cterm = sequence.len().saturating_add(1);

        for (specificity, masses) in &self.variable_mods {
            let choices = masses.len() as u64;
            if choices == 0 {
                continue;
            }
            let mut add_site = |site: usize| {
                let entry = choices_by_site.entry(site).or_default();
                *entry = entry.saturating_add(choices);
            };

            match (*specificity, digest.position) {
                (ModificationSpecificity::PeptideN(None), _) => add_site(nterm),
                (ModificationSpecificity::PeptideN(Some(residue)), _)
                    if sequence.first() == Some(&residue) =>
                {
                    add_site(0)
                }
                (ModificationSpecificity::PeptideC(None), _) => add_site(cterm),
                (ModificationSpecificity::PeptideC(Some(residue)), _)
                    if sequence.last() == Some(&residue) =>
                {
                    add_site(sequence.len().saturating_sub(1))
                }
                (ModificationSpecificity::ProteinN(None), Position::Nterm | Position::Full) => {
                    add_site(nterm)
                }
                (
                    ModificationSpecificity::ProteinN(Some(residue)),
                    Position::Nterm | Position::Full,
                ) if sequence.first() == Some(&residue) => add_site(0),
                (ModificationSpecificity::ProteinC(None), Position::Cterm | Position::Full) => {
                    add_site(cterm)
                }
                (
                    ModificationSpecificity::ProteinC(Some(residue)),
                    Position::Cterm | Position::Full,
                ) if sequence.last() == Some(&residue) => {
                    add_site(sequence.len().saturating_sub(1))
                }
                (ModificationSpecificity::Residue(residue), _) => {
                    for (index, candidate) in sequence.iter().enumerate() {
                        if *candidate == residue {
                            add_site(index);
                        }
                    }
                }
                _ => {}
            }
        }

        let max_mods = self.max_variable_mods.min(choices_by_site.len());
        let mut counts = vec![0u64; max_mods.saturating_add(1)];
        counts[0] = 1;
        for choices in choices_by_site.values().copied() {
            for count in (1..=max_mods).rev() {
                counts[count] =
                    counts[count].saturating_add(counts[count - 1].saturating_mul(choices));
            }
        }
        let variants = counts.into_iter().fold(0u64, u64::saturating_add);
        self.max_combinations
            .map_or(variants, |cap| variants.min(cap as u64))
    }

    pub fn auto_calculate_prefilter_chunk_size(
        &mut self,
        fasta: &Fasta,
        estimated_modified_peptides: u64,
    ) {
        const MAX_PEPS_PER_CHUNK: usize = 2usize.pow(23);
        self.prefilter_chunk_size = match self.prefilter_chunk_size {
            0 => {
                let chunk_count = estimated_modified_peptides
                    .saturating_add(MAX_PEPS_PER_CHUNK as u64 - 1)
                    / MAX_PEPS_PER_CHUNK as u64;
                let chunk_count = chunk_count.max(1);
                ((fasta.targets.len() as u64).saturating_add(chunk_count - 1) / chunk_count).max(1)
                    as usize
            }
            x => x,
        };
    }

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
        let enzyme = self.enzyme.clone().into();
        let digests = fasta.digest_with_custom_cleavages(&enzyme, custom_cleavages);

        log::trace!("grouping digests");
        let start_num = digests.len();
        let digests = group_digests(digests);
        log::trace!(
            "grouped {} digests into {} groups",
            start_num,
            digests.len()
        );
        digests
    }

    /// Expand variable modifications and generate decoys from an unmodified digest.
    pub fn modify_digests(&self, digests: Vec<DigestGroup>) -> Vec<Peptide> {
        let mods = self.variable_modifications();
        let static_mods = self.static_modifications();

        let targets: DashSet<_, FnvBuildHasher> = DashSet::default();
        digests
            .par_iter()
            .filter(|digest| !digest.reference.decoy)
            .for_each(|digest| {
                targets.insert(digest.reference.sequence.clone().into_bytes());
            });

        log::trace!("modifying peptides");
        let mut target_decoys = digests
            .into_par_iter()
            .map(Peptide::try_from)
            .filter_map(Result::ok)
            .flat_map_iter(|peptide| {
                peptide
                    .apply(
                        &mods,
                        &static_mods,
                        self.max_variable_mods,
                        self.max_combinations,
                    )
                    .into_iter()
                    .filter(|peptide| {
                        peptide.monoisotopic >= self.peptide_min_mass
                            && peptide.monoisotopic <= self.peptide_max_mass
                    })
                    .flat_map(|peptide| {
                        if self.generate_decoys {
                            vec![peptide.reverse(), peptide].into_iter()
                        } else {
                            vec![peptide].into_iter()
                        }
                    })
                    .filter(|peptide| !peptide.decoy || !targets.contains(&(peptide.sequence[..])))
            })
            .collect::<Vec<_>>();

        Self::reorder_peptides(&mut target_decoys);

        target_decoys
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
        log::trace!("sorting and deduplicating peptides");

        let init_size = target_decoys.len();
        // This is equivalent to a stable sort
        target_decoys.par_sort_unstable_by(|a, b| {
            a.monoisotopic
                .total_cmp(&b.monoisotopic)
                .then_with(|| a.initial_sort(b))
        });
        target_decoys.dedup_by(|remove, keep| {
            if remove.monoisotopic == keep.monoisotopic
                && remove.sequence == keep.sequence
                && remove.modifications == keep.modifications
                && remove.nterm == keep.nterm
                && remove.cterm == keep.cterm
                && remove.applied_modifications == keep.applied_modifications
            {
                keep.proteins.extend(remove.proteins.iter().cloned());
                // When merging peptides from different Fastas,
                // decoys in one fasta might be targets in another
                keep.decoy &= remove.decoy;
                true
            } else {
                false
            }
        });

        target_decoys
            .par_iter_mut()
            .for_each(|peptide| peptide.proteins.sort_unstable());

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
    /// `protein` and `decoy` are optional. No modifications are applied —
    /// sequences are used verbatim. Decoys are generated (by reversal) if
    /// `self.generate_decoys` is true, subject to the same deduplication as
    /// the normal FASTA digest path.
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
            .filter_map(|line| {
                let fields: Vec<&str> = line.split('\t').collect();
                let seq = fields.get(seq_col)?.trim().to_string();
                if seq.is_empty() {
                    return None;
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

                let digest = Digest {
                    decoy: is_decoy,
                    semi_enzymatic: false,
                    sequence: seq,
                    protein,
                    missed_cleavages: 0,
                    position: Position::Full,
                };
                match Peptide::try_from(digest) {
                    Ok(p)
                        if p.monoisotopic >= self.peptide_min_mass
                            && p.monoisotopic <= self.peptide_max_mass =>
                    {
                        Some(p)
                    }
                    Ok(_) => None,
                    Err(e) => {
                        log::warn!("skipping invalid peptide sequence: {e:?}");
                        None
                    }
                }
            })
            .collect();

        // Build target sequence set for decoy deduplication.
        let targets: DashSet<Vec<u8>, FnvBuildHasher> = DashSet::default();
        raw.iter().filter(|p| !p.decoy).for_each(|p| {
            targets.insert(p.sequence.to_vec());
        });

        // Emit targets (+ generated decoys) into the final list.
        let mut result: Vec<Peptide> = raw
            .into_iter()
            .flat_map(|peptide| {
                if self.generate_decoys && !peptide.decoy {
                    let rev = peptide.reverse();
                    if !targets.contains(&rev.sequence[..]) {
                        vec![rev, peptide]
                    } else {
                        vec![peptide]
                    }
                } else {
                    vec![peptide]
                }
            })
            .collect();

        Self::reorder_peptides(&mut result);
        result
    }

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

        // Finally, perform in silico digest for our target sequences
        // Note that multiple charge states are actually handled by
        // [`SpectrumProcessor`] or during scoring - all theoretical
        // fragments are monoisotopic/uncharged
        let mut fragments = if self.use_bitmap {
            Vec::new()
        } else {
            target_decoys
                .par_iter()
                .enumerate()
                .flat_map_iter(|(idx, peptide)| {
                    // Generate both B and Y ions, then filter down to make sure that
                    // theoretical fragments are within the search space
                    self.ion_kinds
                        .iter()
                        .flat_map(|kind| IonGroupSeries::new(peptide, *kind))
                        .filter(|group| {
                            // Don't store b1, b2, y1, y2 ions for preliminary scoring

                            match group.kind {
                                Kind::A | Kind::B | Kind::C => {
                                    (group.series_index + 1) > self.min_ion_index
                                }
                                Kind::X | Kind::Y | Kind::Z => {
                                    peptide.sequence.len().saturating_sub(1) - group.series_index
                                        > self.min_ion_index
                                }
                            }
                        })
                        // Keep one canonical form per cleavage in the preliminary
                        // index. Full rescoring evaluates every neutral-loss
                        // alternative as a group; indexing all alternatives here
                        // would bias candidate selection toward configured mods.
                        .filter_map(move |group| {
                            group.variants.into_iter().next().map(|ion| Theoretical {
                                peptide_index: PeptideIx(idx as u32),
                                fragment_mz: ion.monoisotopic_mass,
                            })
                        })
                })
                .collect::<Vec<_>>()
        };
        log::trace!("finalizing index");

        // Sort all of our theoretical fragments by m/z, from low to high
        fragments.par_sort_unstable_by(|a, b| a.fragment_mz.total_cmp(&b.fragment_mz));

        // Now, we bucket all of our theoretical fragments, and within each bucket
        // sort by precursor m/z - and save the minimum *fragment* m/z in a separate
        // vector so that we can perform an efficient binary search to reduce
        // the number of in silico fragments we evaluate
        //
        // Imagine our theoretical fragments look like this
        //
        // Fragment        A      B       C       D       E       F       G       H
        // Fragment m/z [ 1.0    1.2     1.3     2.5     2.5     2.6     3.5     4.0 ]
        // Parent m/z   [ 500    439     291     800     142     515     517     232 ]
        //
        // If we apply a bucket size of 4 we will end up with the following:
        //
        // Fragment        C      B       A       D       E       H       F       G
        // Fragment m/z [ 1.3    1.2     1.0     2.5     2.5     4.0     3.5     2.6 ]
        // Parent m/z   [ 291    439     500     800     142     232     515     517 ]
        //              |___________________________|   |____________________________|
        //               Bucket 1: min m/z 1.0          Bucket 2: min m/z 2.5
        //
        // * Example query: Fragment m/z 1.3 - 1.9 & Precursor m/z: 450 - 900
        // 1) Perform a binary search to narrow down our window to Bucket 1 only
        //      * Bucket 2 has a min m/z outside of our query range - nothing here can match
        //
        // Fragment        C      B       A       D
        // Fragment m/z [ 1.3    1.2     1.0     2.5
        // Parent m/z   [ 291    439     500     800
        //                            |_____________|
        //                                    ^
        //                                    |
        // Window with matching precursors ___|

        // and within Bucket 1, we can perform another binary search to find fragments
        // matching our desired precursor m/z tolerance

        let min_value = fragments
            .par_chunks_mut(self.bucket_size)
            .map(|chunk| {
                // There should always be at least one item in the chunk!
                //  we know the chunk is already sorted by fragment_mz too, so this is minimum value
                let min = chunk[0].fragment_mz;
                chunk.par_sort_unstable_by(|a, b| a.peptide_index.cmp(&b.peptide_index));
                min
            })
            .collect::<Vec<_>>();

        // PTM localization works from the compact mass/specificity list below,
        // while ambiguity and site reports resolve display labels by mass.
        // Preserve names from Sage Plus's structured modification definitions.
        for entries in self.variable_mods.values() {
            for entry in entries {
                let definition = entry.definition();
                if let Some(name) = definition.name.as_deref() {
                    crate::unimod::register_label(definition.mass, name);
                }
            }
        }

        let potential_mods = self
            .variable_mods
            .iter()
            .flat_map(|(specificity, entries)| {
                entries.iter().map(|entry| (*specificity, entry.mass()))
            })
            .collect::<Vec<(ModificationSpecificity, f32)>>();

        let bitmap_index = if self.use_bitmap {
            BitmapIndex::build(
                &target_decoys,
                &self.ion_kinds,
                self.bitmap_size,
                self.peptide_min_mass,
                self.peptide_max_mass,
            )
        } else {
            BitmapIndex::default()
        };

        IndexedDatabase {
            peptides: target_decoys,
            fragments,
            min_value,
            bucket_size: self.bucket_size,
            ion_kinds: self.ion_kinds,
            generate_decoys: self.generate_decoys,
            potential_mods,
            decoy_tag: self.decoy_tag,
            bitmap_index,
        }
    }
}

fn with_estimation_margin(bytes: u64) -> u64 {
    // Parallel collection, allocator size classes, and sorting create overhead that
    // cannot be derived exactly from item counts. Use a conservative 50% margin.
    bytes.saturating_add(bytes / 2)
}

#[derive(Hash, Copy, Clone, PartialEq, Eq, PartialOrd, Ord, Debug, Serialize)]
#[repr(transparent)]
pub struct PeptideIx(pub u32);

// This is unsafe for use outside of this crate
impl Default for PeptideIx {
    fn default() -> Self {
        Self(u32::MAX)
    }
}

#[derive(Copy, Clone, Debug, PartialEq, Serialize)]
pub struct Theoretical {
    pub peptide_index: PeptideIx,
    pub fragment_mz: f32,
}

#[derive(Default)]
pub struct IndexedDatabase {
    pub peptides: Vec<Peptide>,
    pub fragments: Vec<Theoretical>,
    pub ion_kinds: Vec<Kind>,
    pub min_value: Vec<f32>,
    /// Keep a list of potential (AA, mass) modifications for RT prediction
    pub potential_mods: Vec<(ModificationSpecificity, f32)>,
    pub bucket_size: usize,
    pub generate_decoys: bool,
    pub decoy_tag: String,
    /// Bitmap-based preliminary search index (forward/reverse ion bitsets per peptide).
    pub bitmap_index: BitmapIndex,
}

impl IndexedDatabase {
    /// Create a new [`IndexedQuery`] for a specific [`ProcessedSpectrum`]
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

    pub fn serialize(&self) {
        use std::io::Write;
        let mut wtr = std::io::BufWriter::new(std::fs::File::create("fragments.bin").unwrap());
        for fragment in &self.fragments {
            let _ = wtr.write(&fragment.fragment_mz.to_le_bytes()).unwrap();
            let _ = wtr.write(&fragment.peptide_index.0.to_le_bytes()).unwrap();
        }
        wtr.flush().unwrap();

        let mut wtr = std::io::BufWriter::new(std::fs::File::create("peptides.csv").unwrap());
        writeln!(wtr, "peptide,proteins,monoisotopic,decoy").unwrap();
        for fragment in &self.peptides {
            writeln!(
                wtr,
                "{},{},{},{}",
                fragment,
                fragment.proteins(&self.decoy_tag, self.generate_decoys),
                fragment.monoisotopic,
                fragment.decoy
            )
            .unwrap();
        }
        wtr.flush().unwrap();
    }
}

impl std::ops::Index<PeptideIx> for IndexedDatabase {
    type Output = Peptide;

    fn index(&self, index: PeptideIx) -> &Self::Output {
        &self.peptides[index.0 as usize]
    }
}

pub struct IndexedQuery<'d> {
    db: &'d IndexedDatabase,
    precursor_mass: f32,
    precursor_tol: Tolerance,
    fragment_tol: Tolerance,
    pub pre_idx_lo: usize,
    pub pre_idx_hi: usize,
}

impl IndexedQuery<'_> {
    /// Search for a specified `fragment_mz` within the database
    pub fn page_search(&self, mass: f32) -> impl Iterator<Item = &Theoretical> {
        let (fragment_lo, fragment_hi) = self.fragment_tol.bounds(mass);
        let (precursor_lo, precursor_hi) = self.precursor_tol.bounds(self.precursor_mass);

        // Locate the left and right page indices that contain matching fragments
        // Note that we need to multiply by `bucket_size` to transform these into
        // indices that can be used with `self.db.fragments`
        let (left_idx, right_idx) = binary_search_slice(
            &self.db.min_value,
            |min, bounds| min.total_cmp(bounds),
            fragment_lo,
            fragment_hi,
        );

        // It is absolutely critical that we do not cross page boundaries!
        // If we do, we can no longer rely on total ordering of peptide_index (precursor m/z)
        (left_idx..right_idx).flat_map(move |page| {
            let left_idx = page * self.db.bucket_size;
            // Last chunk not guaranted to be modulo bucket size, make sure we don't
            // accidentally go out of bounds!
            let right_idx = ((page + 1) * self.db.bucket_size).min(self.db.fragments.len());

            // Narrow down into our region of interest, then perform another binary
            // search to further refine down to the slice of matching precursor mzs
            let slice = &&self.db.fragments[left_idx..right_idx];

            let (inner_left, inner_right) = binary_search_slice(
                slice,
                |frag, bounds| (frag.peptide_index.0 as usize).cmp(bounds),
                self.pre_idx_lo,
                self.pre_idx_hi,
            );

            // Finally, filter down our slice into exact matches only
            slice[inner_left..inner_right].iter().filter(move |frag| {
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
            })
        })
    }
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

#[cfg(test)]
mod test {
    use std::sync::Arc;

    use super::*;
    use crate::cleavage::CustomCleavageLibrary;

    #[test]
    fn binary_search_slice_smoke() {
        // Make sure that our query returns the maximal set of indices
        let data = [1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0];
        let bounds = binary_search_slice(&data, |a: &f64, b| a.total_cmp(b), 1.75, 3.5);
        assert_eq!(bounds, (1, 6));
        assert!(data[bounds.0] <= 1.75);
        assert_eq!(&data[bounds.0..bounds.1], &[1.5, 2.0, 2.5, 3.0, 3.5]);

        let bounds = binary_search_slice(&data, |a: &f64, b| a.total_cmp(b), 0.0, 5.0);
        assert_eq!(bounds, (0, data.len()));
    }

    #[test]
    fn binary_search_slice_run() {
        // Make sure that our query returns the maximal set of indices
        let data = [1.0, 1.5, 1.5, 1.5, 1.5, 2.0, 2.5, 3.0, 3.0, 3.5, 4.0];
        let (left, right) = binary_search_slice(&data, |a: &f64, b| a.total_cmp(b), 1.5, 3.25);
        assert!(data[left] <= 1.5);
        assert!(data[right] > 3.25);
        assert_eq!(
            &data[left..right],
            &[1.0, 1.5, 1.5, 1.5, 1.5, 2.0, 2.5, 3.0, 3.0]
        );
    }

    #[test]
    fn structured_variable_mod_config_round_trips() {
        let builder: Builder = serde_json::from_value(serde_json::json!({
            "fasta": "none",
            "static_mods": {
                "C": {
                    "mass": 57.0215,
                    "name": "Carbamidomethyl"
                }
            },
            "variable_mods": {
                "M": [15.9949],
                "K": [
                    {
                        "mass": 42.0106,
                        "max_count": 1,
                        "name": "Acetyl",
                        "neutral_losses": [17.0265],
                        "neutral_loss_mode": "required"
                    },
                    {"mass": 14.0157}
                ]
            },
            "max_variable_mods": 2,
            "max_combinations": 0
        }))
        .unwrap();

        let params = builder.make_parameters();
        assert_eq!(params.max_variable_mods, 2);
        assert_eq!(params.max_combinations, Some(1));

        let mods = params.variable_modifications();
        assert_eq!(mods.len(), 3);
        assert_eq!(mods[0].0, ModificationSpecificity::Residue(b'K'));
        assert!((mods[0].1.mass - 42.0106).abs() < 1e-4);
        assert_eq!(mods[0].2, Some(1));
        assert_eq!(mods[0].1.name.as_deref(), Some("Acetyl"));
        assert_eq!(&*mods[0].1.neutral_losses, &[17.0265]);
        assert_eq!(mods[1].0, ModificationSpecificity::Residue(b'K'));
        assert!((mods[1].1.mass - 14.0157).abs() < 1e-4);
        assert_eq!(mods[1].2, None);
        assert_eq!(mods[2].0, ModificationSpecificity::Residue(b'M'));
        assert!((mods[2].1.mass - 15.9949).abs() < 1e-4);
        assert_eq!(mods[2].2, None);

        let serialized = serde_json::to_value(params).unwrap();
        let k_entries = &serialized["variable_mods"]["K"];
        assert!(k_entries[0].is_object());
        assert_eq!(k_entries[0]["max_count"], 1);
        assert_eq!(k_entries[0]["name"], "Acetyl");
        assert_eq!(k_entries[0]["neutral_loss_mode"], "required");
        assert!(k_entries[1].is_object());
        assert!(k_entries[1].get("max_count").is_none());
        assert!(serialized["variable_mods"]["M"][0].is_number());
        assert_eq!(serialized["static_mods"]["C"]["name"], "Carbamidomethyl");
    }

    #[test]
    fn digestion() {
        let fasta = r#"
        >sp|AAAAA
        MEWKLEQSMREQALLKAQLTQLK
        >sp|BBBBB
        RMEWKLEQSMREQALLKAQLTQLK
        "#;

        let fasta = Fasta::parse(fasta.into(), "rev_", false);

        // Make sure that FASTA parsed OK
        assert_eq!(
            fasta.targets,
            vec![
                (
                    Arc::from("sp|AAAAA".to_string()),
                    "MEWKLEQSMREQALLKAQLTQLK".into()
                ),
                (
                    Arc::from("sp|BBBBB".to_string()),
                    "RMEWKLEQSMREQALLKAQLTQLK".into()
                ),
            ]
        );

        let params = Parameters {
            bucket_size: 128,
            bitmap_size: 30,
            enzyme: EnzymeBuilder {
                missed_cleavages: Some(1),
                min_len: Some(6),
                max_len: Some(10),
                ..Default::default()
            },
            peptide_min_mass: 150.0,
            peptide_max_mass: 5000.0,
            ion_kinds: vec![Kind::B, Kind::Y],
            min_ion_index: 2,
            static_mods: HashMap::default(),
            variable_mods: [(
                ModificationSpecificity::ProteinN(None),
                vec![VarModEntry::Mass(42.0)],
            )]
            .into_iter()
            .collect(),
            max_variable_mods: 2,
            max_combinations: None,
            decoy_tag: "rev_".into(),
            generate_decoys: false,
            fasta: "none".into(),
            peptides: None,
            custom_cleavage_sites: None,
            prefilter: false,
            prefilter_chunk_size: 0,
            prefilter_low_memory: true,
            use_bitmap: false,
        };

        let peptides = params.digest(&fasta);

        let expected = [
            "EQALLK",
            "LEQSMR",
            "AQLTQLK",
            "MEWKLEQSMR",
            "[+42]-MEWKLEQSMR",
        ]
        .into_iter()
        .map(String::from)
        .collect::<Vec<_>>();

        let sequences = peptides.iter().map(|p| p.to_string()).collect::<Vec<_>>();
        assert_eq!(expected, sequences);

        // All peptides are shared except for the protein N-term mod
        for peptide in &peptides[..4] {
            assert_eq!(peptide.proteins.len(), 2, "{:?}", peptide);
        }
        // Ensure that this mod is uniquely called as the first protein
        assert_eq!(
            peptides.last().unwrap().proteins,
            vec!["sp|AAAAA".to_string().into()]
        );
    }

    #[test]
    fn custom_cleavages_flow_through_modification_and_memory_paths() {
        let fasta = Fasta::parse(">P1\nAAKAPEPTIDERQQQK\n".into(), "rev_", true);
        let library =
            CustomCleavageLibrary::from_tsv("protein\tposition\tcontext\nP1\t7\tAPEPT|IDER\n")
                .unwrap()
                .validate(&fasta)
                .unwrap();
        let mut builder = Builder::default();
        builder.enzyme = Some(EnzymeBuilder {
            missed_cleavages: Some(1),
            min_len: Some(3),
            max_len: Some(50),
            ..Default::default()
        });
        builder.generate_decoys = Some(false);
        let parameters = builder.make_parameters();

        let ordinary = parameters.digest(&fasta);
        let custom = parameters.digest_with_custom_cleavages(&fasta, Some(&library));
        let custom_sequences = custom
            .iter()
            .map(|peptide| std::str::from_utf8(&peptide.sequence).unwrap())
            .collect::<Vec<_>>();
        assert!(custom.len() > ordinary.len());
        assert!(custom_sequences.contains(&"APEPT"));
        assert!(custom_sequences.contains(&"IDER"));

        let estimate = parameters.estimate_memory_with_custom_cleavages(&fasta, Some(&library));
        assert!(estimate.modified_peptides as usize >= custom.len());
        assert!(estimate.modified_peptides > parameters.estimate_memory(&fasta).modified_peptides);
    }

    #[test]
    fn builds_only_selected_search_index() {
        let peptide = Peptide::try_from(Digest {
            sequence: "PEPTIDER".into(),
            protein: Arc::from("protein"),
            ..Digest::default()
        })
        .unwrap();

        let parameters = Builder::default().make_parameters();
        let fragment_database = parameters
            .clone()
            .build_from_peptides(vec![peptide.clone()]);
        assert!(!fragment_database.fragments.is_empty());
        assert!(fragment_database.bitmap_index.forward_bitmaps.is_empty());
        assert!(fragment_database.bitmap_index.reverse_bitmaps.is_empty());

        let mut parameters = parameters;
        parameters.use_bitmap = true;
        let bitmap_database = parameters.build_from_peptides(vec![peptide]);
        assert!(bitmap_database.fragments.is_empty());
        assert!(bitmap_database.min_value.is_empty());
        assert!(!bitmap_database.bitmap_index.forward_bitmaps.is_empty());
        assert!(!bitmap_database.bitmap_index.reverse_bitmaps.is_empty());
    }

    #[test]
    fn estimates_variable_modification_expansion_before_allocation() {
        let mut builder = Builder::default();
        builder.enzyme = Some(EnzymeBuilder {
            cleave_at: Some("$".into()),
            min_len: Some(1),
            max_len: Some(50),
            ..Default::default()
        });
        builder.variable_mods = Some(
            [(
                "S".to_string(),
                vec![VarModEntry::Mass(79.9663), VarModEntry::Mass(80.0)],
            )]
            .into_iter()
            .collect(),
        );
        builder.max_variable_mods = Some(3);
        builder.generate_decoys = Some(false);
        let parameters = builder.make_parameters();
        let fasta = Fasta::parse(">protein\nSSSSSSSSSS\n".into(), "rev_", false);

        let estimate = parameters.estimate_memory(&fasta);

        // 1 + C(10,1)*2 + C(10,2)*2^2 + C(10,3)*2^3
        assert_eq!(estimate.unmodified_peptides, 1);
        assert_eq!(estimate.modified_peptides, 1_161);
        assert_eq!(estimate.fragments, 1_161 * 14);
        assert!(estimate.unmodified_peak_bytes > 0);
        assert!(estimate.modified_peak_bytes > estimate.unmodified_peak_bytes);
        assert!(estimate.fragment_peak_bytes > estimate.modified_peak_bytes);

        let digests = parameters.digest_unmodified(&fasta);
        let modification_estimate = parameters.estimate_modified_memory(&digests);
        assert_eq!(modification_estimate.modified_peptides, 1_161);
        assert_eq!(parameters.modify_digests(digests).len(), 1_161);
    }
}
