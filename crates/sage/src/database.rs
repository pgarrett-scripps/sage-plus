use crate::cleavage::ValidatedCustomCleavageLibrary;
use crate::enzyme::{
    group_digests, Digest, DigestGroup, Enzyme, EnzymeParameters, Position, ProteinOccurrence,
};
use crate::fasta::Fasta;
use crate::ion_series::{IonGroupSeries, Kind};
use crate::mass::Tolerance;
use crate::modification::{
    validate_mods, validate_var_mods, ModificationDefinition, ModificationSpecificity, SiteMode,
    StaticModEntry, VarModEntry,
};
use crate::peptide::{
    AppliedModification, LabelModificationCache, LibrarySite, ModificationKind, ModificationLookup,
    ModificationPlan, Peptide, VariableRule, INLINE_PROTEINS,
};
use crate::ptm_library::PtmLibrary;
use crate::sequence::PeptideSequence;
use dashmap::DashSet;
use fnv::FnvBuildHasher;
use rayon::prelude::*;
use serde::{Deserialize, Serialize};
use std::cmp::Ordering;
use std::collections::{HashMap, HashSet};
use std::hash::Hash;
use std::sync::Arc;

#[derive(Deserialize, Serialize, Clone, Debug, schemars::JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct EnzymeBuilder {
    /// How many missed cleavages to use
    pub missed_cleavages: Option<u8>,
    /// Minimum peptide length that will be fragmented
    #[schemars(range(min = 1))]
    pub min_len: Option<usize>,
    /// Maximum peptide length that will be fragmented
    #[schemars(range(min = 1))]
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

#[derive(Deserialize, Default, Clone, schemars::JsonSchema)]
#[serde(deny_unknown_fields)]
/// Parameters used for generating the fragment database
pub struct Builder {
    /// This parameter allows tuning of the internal search structure
    pub bucket_size: Option<usize>,
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
    /// Maximum number of variable modifications after exhaustive and
    /// library-supported placements are combined.
    pub max_total_variable_mods: Option<usize>,
    /// Optional site library. Modification definitions remain in `variable_mods`.
    pub ptm_library: Option<PtmLibrarySettings>,
    /// Use this prefix for decoy proteins
    pub decoy_tag: Option<String>,

    pub generate_decoys: Option<bool>,
    /// Path to fasta database
    pub fasta: Option<String>,
    /// Path to a pre-digested peptide TSV file (additive with `fasta`).
    /// Required column: `sequence`. Optional columns: `protein`, `decoy`.
    /// Configured static, variable, and channel-aware modifications are applied.
    pub peptides: Option<String>,
    /// Path to a protein-specific custom cleavage-site TSV or Parquet file.
    /// Required columns: `protein`, `position`; optional column: `context`.
    pub custom_cleavage_sites: Option<String>,
    /// Number of sequences to handle simultaneously when pre-filtering the db
    pub prefilter_chunk_size: Option<usize>,
    /// Pre-filter the database to minimize memory usage
    pub prefilter: Option<bool>,
    /// Deprecated compatibility option. Exact prefiltering always uses compact
    /// survivor tracking and ignores this value.
    pub prefilter_low_memory: Option<bool>,
}

impl Builder {
    pub fn make_parameters(self) -> Parameters {
        if self.prefilter_low_memory.is_some() {
            log::warn!("database.prefilter_low_memory is deprecated and ignored");
        }
        let bucket_size = self.bucket_size.unwrap_or(8192).next_power_of_two();
        let max_variable_mods = self.max_variable_mods.map(|x| x.max(1)).unwrap_or(2);
        let max_total_variable_mods = self
            .max_total_variable_mods
            .map(|x| x.max(1))
            .unwrap_or(max_variable_mods)
            .max(max_variable_mods);
        Parameters {
            bucket_size,
            peptide_min_mass: self.peptide_min_mass.unwrap_or(500.0),
            peptide_max_mass: self.peptide_max_mass.unwrap_or(5000.0),
            ion_kinds: self.ion_kinds.unwrap_or(vec![Kind::B, Kind::Y]),
            min_ion_index: self.min_ion_index.unwrap_or(2),
            decoy_tag: self.decoy_tag.unwrap_or_else(|| "rev_".into()),
            enzyme: self.enzyme.unwrap_or_default(),
            static_mods: validate_mods(self.static_mods),
            variable_mods: validate_var_mods(self.variable_mods),
            max_variable_mods,
            max_combinations: self.max_combinations.map(|x| x.max(1)),
            max_total_variable_mods,
            ptm_library: self.ptm_library,
            generate_decoys: self.generate_decoys.unwrap_or(true),
            fasta: self.fasta.unwrap_or_default(),
            peptides: self.peptides,
            custom_cleavage_sites: self.custom_cleavage_sites,
            prefilter_chunk_size: self.prefilter_chunk_size.unwrap_or(0),
            prefilter: self.prefilter.unwrap_or(false),
            loaded_ptm_library: None,
        }
    }

    pub fn update_fasta(&mut self, fasta: String) {
        self.fasta = Some(fasta)
    }
}

#[derive(Serialize, Clone, Debug)]
pub struct Parameters {
    pub bucket_size: usize,
    pub enzyme: EnzymeBuilder,
    pub peptide_min_mass: f32,
    pub peptide_max_mass: f32,
    pub ion_kinds: Vec<Kind>,
    pub min_ion_index: usize,
    pub static_mods: HashMap<ModificationSpecificity, StaticModEntry>,
    pub variable_mods: HashMap<ModificationSpecificity, Vec<VarModEntry>>,
    pub max_variable_mods: usize,
    pub max_combinations: Option<usize>,
    pub max_total_variable_mods: usize,
    pub ptm_library: Option<PtmLibrarySettings>,
    pub decoy_tag: String,
    pub generate_decoys: bool,
    pub fasta: String,
    pub peptides: Option<String>,
    pub custom_cleavage_sites: Option<String>,
    pub prefilter_chunk_size: usize,
    pub prefilter: bool,
    #[serde(skip)]
    pub loaded_ptm_library: Option<Arc<PtmLibrary>>,
}

#[derive(Deserialize, Serialize, Clone, Debug, schemars::JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct PtmLibrarySettings {
    pub path: String,
    #[serde(default = "default_true")]
    pub strict: bool,
}

fn default_true() -> bool {
    true
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
    pub fn validate_compact_modifications(&self) -> Result<(), String> {
        let max_len = self.enzyme.max_len.unwrap_or(50);
        if max_len > u8::MAX as usize {
            return Err(format!(
                "database.enzyme.max_len must not exceed {} residues for compact modification encoding, but is {max_len}",
                u8::MAX
            ));
        }

        let variable_mods = self.variable_modifications();
        let static_mods = self.static_modifications();
        let channels = self.label_channels();
        let labels = LabelModificationCache::new(
            variable_mods
                .iter()
                .map(|rule| &rule.modification)
                .chain(static_mods.values()),
            &channels,
        );
        ModificationLookup::for_rules(&variable_mods, &static_mods, &channels, &labels)
            .map(|_| ())
            .map_err(|error| error.to_string())
    }

    pub fn validate_channels(&self) -> Result<(), String> {
        let definitions = self.channel_definitions();
        let Some(first) = definitions.first() else {
            return Ok(());
        };
        let expected = first.channel_offsets.keys().cloned().collect::<Vec<_>>();
        if expected.len() < 2 {
            return Err(
                "channel_offsets must define at least two channels on every channel-aware modification"
                    .into(),
            );
        }
        for definition in definitions.iter().skip(1) {
            if definition.channel_offsets.keys().ne(expected.iter()) {
                return Err(
                    "all channel_offsets dictionaries must use exactly the same channel names"
                        .into(),
                );
            }
        }
        if definitions.iter().all(|definition| {
            definition
                .channel_offsets
                .values()
                .all(|offset| *offset == 0.0)
        }) {
            return Err("channel_offsets must contain at least one non-zero offset".into());
        }
        let mut signatures = HashSet::new();
        for channel in &expected {
            let signature = definitions
                .iter()
                .map(|definition| definition.channel_offsets[channel].to_bits())
                .collect::<Vec<_>>();
            if !signatures.insert(signature) {
                return Err(format!(
                    "channel `{channel}` is chemically identical to another configured channel"
                ));
            }
        }
        Ok(())
    }

    fn channel_definitions(&self) -> Vec<ModificationDefinition> {
        self.static_mods
            .values()
            .map(StaticModEntry::definition)
            .chain(
                self.variable_mods
                    .values()
                    .flatten()
                    .map(VarModEntry::definition),
            )
            .filter(|definition| !definition.channel_offsets.is_empty())
            .collect()
    }

    fn label_channels(&self) -> Vec<Arc<str>> {
        self.channel_definitions()
            .first()
            .map(|definition| definition.channel_offsets.keys().cloned().collect())
            .unwrap_or_default()
    }

    fn label_reference(&self) -> Option<Arc<str>> {
        let definitions = self.channel_definitions();
        self.label_channels().into_iter().find(|channel| {
            definitions
                .iter()
                .all(|definition| definition.channel_offsets[channel] == 0.0)
        })
    }

    /// Flatten variable modifications into a stable order. This matters when
    /// `max_combinations` truncates variants: equivalent configurations must
    /// retain the same variants regardless of randomized `HashMap` iteration.
    fn variable_modifications(&self) -> Vec<VariableRule> {
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
                        entry.site_mode(),
                    )
                })
            })
            .collect::<Vec<_>>();
        mods.sort_unstable_by(|a, b| a.0.cmp(&b.0).then_with(|| a.1.cmp(&b.1)));
        let mut named_groups: HashMap<Arc<str>, usize> = HashMap::new();
        let mut next_group = 0usize;
        mods.into_iter()
            .map(|(specificity, _, modification, max_count, site_mode)| {
                let count_group = if let Some(name) = modification.name.clone() {
                    *named_groups.entry(name).or_insert_with(|| {
                        let group = next_group;
                        next_group += 1;
                        group
                    })
                } else {
                    let group = next_group;
                    next_group += 1;
                    group
                };
                VariableRule {
                    specificity,
                    modification,
                    max_count,
                    site_mode,
                    count_group,
                }
            })
            .collect()
    }

    pub fn validate_ptm_library(&self, library: &PtmLibrary) -> Result<(), String> {
        if self.max_total_variable_mods < self.max_variable_mods {
            return Err(
                "database.max_total_variable_mods must be at least database.max_variable_mods"
                    .into(),
            );
        }

        let rules = self.variable_modifications();
        let mut definitions: HashMap<&str, (&ModificationDefinition, Option<usize>, SiteMode)> =
            HashMap::new();
        for rule in &rules {
            if self.ptm_library.is_some() && rule.max_count.is_none() {
                return Err(
                    "all variable modifications require `max_count` when database.ptm_library is configured"
                        .into(),
                );
            }
            if rule.site_mode != SiteMode::Exhaustive
                && (rule.modification.name.is_none() || rule.max_count.is_none())
            {
                return Err(
                    "variable modifications using `library` or `both` require `name` and `max_count`"
                        .into(),
                );
            }
            if let Some(name) = rule.modification.name.as_deref() {
                if let Some((definition, max_count, site_mode)) = definitions.get(name) {
                    if *definition != rule.modification.as_ref()
                        || *max_count != rule.max_count
                        || *site_mode != rule.site_mode
                    {
                        return Err(format!(
                            "variable modification `{name}` has inconsistent definitions across specificities"
                        ));
                    }
                } else {
                    definitions.insert(name, (&rule.modification, rule.max_count, rule.site_mode));
                }
            }
        }

        for site in library.iter() {
            match definitions.get(site.modification.as_ref()) {
                None => {
                    return Err(format!(
                        "PTM library references undefined modification `{}`",
                        site.modification
                    ))
                }
                Some((_, _, SiteMode::Exhaustive)) => {
                    return Err(format!(
                        "PTM library modification `{}` must use site_mode `library` or `both`",
                        site.modification
                    ))
                }
                _ => {}
            }
        }
        Ok(())
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

        let fragment_bytes = estimate
            .fragments
            .saturating_mul(std::mem::size_of::<Theoretical>() as u64);
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
            .saturating_mul(std::mem::size_of::<Theoretical>() as u64);
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
                        self.variable_variant_count(&reference)
                    })
                    .fold(0u64, u64::saturating_add)
            } else {
                self.variable_variant_count(&digest.reference)
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
        let bucket_size = self.bucket_size.max(1) as u64;
        (fragments.saturating_add(bucket_size - 1) / bucket_size)
            .saturating_mul(std::mem::size_of::<f32>() as u64)
    }

    fn variable_variant_count(&self, digest: &Digest) -> u64 {
        let sequence = digest.sequence.as_bytes();
        let rules = self.variable_modifications();
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
                    .map(|site| ((site.position - start) as usize, site.modification.clone()))
                    .collect::<HashSet<_>>()
            })
            .unwrap_or_default();
        let mut candidates: HashMap<(usize, usize), bool> = HashMap::new();
        let nterm = sequence.len();
        let cterm = sequence.len().saturating_add(1);

        for rule in &rules {
            let mut add_site = |site: usize| {
                let library_position = if site == nterm {
                    0
                } else if site == cterm {
                    sequence.len().saturating_sub(1)
                } else {
                    site
                };
                let supported = rule.site_mode != SiteMode::Exhaustive
                    && rule.modification.name.as_ref().is_some_and(|name| {
                        !sequence.is_empty()
                            && library_sites.contains(&(library_position, name.clone()))
                    });
                if rule.site_mode != SiteMode::Library || supported {
                    candidates
                        .entry((site, rule.count_group))
                        .and_modify(|library| *library |= supported)
                        .or_insert(supported);
                }
            };

            match (rule.specificity, digest.position) {
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

        let mut choices_by_site: HashMap<usize, (u64, u64)> = HashMap::new();
        for ((site, _), library_supported) in candidates {
            let choices = choices_by_site.entry(site).or_default();
            if library_supported {
                choices.0 += 1;
            } else {
                choices.1 += 1;
            }
        }

        let max_total = self.max_total_variable_mods.min(choices_by_site.len());
        let max_exhaustive = self.max_variable_mods.min(max_total);
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
        let variants = counts.into_iter().flatten().fold(0u64, u64::saturating_add);
        let variable_variants = self
            .max_combinations
            .map_or(variants, |cap| variants.min(cap as u64));
        variable_variants.saturating_mul(self.label_channels().len().max(1) as u64)
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
        let target_sequences = digests
            .iter()
            .filter(|digest| !digest.reference.decoy)
            .map(|digest| digest.reference.sequence.clone())
            .collect::<HashSet<_>>();
        self.modify_digests_with_target_sequences(digests, &target_sequences)
    }

    /// Expand a digest chunk while checking generated decoys against every
    /// target sequence in the complete database.
    pub fn modify_digests_with_target_sequences(
        &self,
        digests: Vec<DigestGroup>,
        target_sequences: &HashSet<PeptideSequence>,
    ) -> Vec<Peptide> {
        let mods = self.variable_modifications();
        let static_mods = self.static_modifications();
        let label_channels = self.label_channels();
        let label_reference = self.label_reference();
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
        let library = self.loaded_ptm_library.as_deref();

        log::trace!("modifying peptides");
        let mut target_decoys = digests
            .into_par_iter()
            .flat_map_iter(|group| {
                let expand = |peptide: Peptide, library_sites: &[LibrarySite]| {
                    let decoy_sequence = self
                        .generate_decoys
                        .then(|| peptide.sequence.reversed_internal());
                    peptide
                        .apply_rules(&modification_plan, library_sites)
                        .into_iter()
                        .flat_map(|peptide| {
                            if label_channels.is_empty() {
                                vec![peptide]
                            } else {
                                label_channels
                                    .iter()
                                    .map(|channel| {
                                        peptide.clone().apply_label_channel(
                                            channel.clone(),
                                            &label_modifications,
                                        )
                                    })
                                    .collect()
                            }
                        })
                        .filter(|peptide| {
                            peptide.monoisotopic >= self.peptide_min_mass
                                && peptide.monoisotopic <= self.peptide_max_mass
                        })
                        .flat_map(|peptide| {
                            if let Some(sequence) = &decoy_sequence {
                                vec![peptide.reverse_with_sequence(sequence.clone()), peptide]
                                    .into_iter()
                            } else {
                                vec![peptide].into_iter()
                            }
                        })
                        .filter(|peptide| {
                            !peptide.decoy || !target_sequences.contains(&(peptide.sequence[..]))
                        })
                        .collect::<Vec<_>>()
                };

                match library {
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
                                let start = origin.start.unwrap_or_default();
                                let end = start.saturating_add(peptide.sequence.len() as u32);
                                let library_sites = library
                                    .sites_for(&origin.protein)
                                    .iter()
                                    .filter(|site| (start..end).contains(&site.position))
                                    .filter_map(|site| {
                                        let position = site.position - start;
                                        (peptide.sequence.get(position as usize)
                                            == Some(&site.residue))
                                        .then(|| LibrarySite {
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
            })
            .collect::<Vec<_>>();

        Self::reorder_peptides_with_reference(&mut target_decoys, label_reference.as_deref());

        target_decoys
    }

    /// Partition digest groups without splitting equal raw sequences. This
    /// keeps terminal variants and protein occurrences together so a modified
    /// peptidoform is generated in exactly one prefilter chunk.
    pub fn partition_digests_by_sequence(
        mut digests: Vec<DigestGroup>,
        chunk_size: usize,
    ) -> Vec<Vec<DigestGroup>> {
        let chunk_size = chunk_size.max(1);
        digests.sort_unstable_by(|left, right| {
            left.reference
                .sequence
                .cmp(&right.reference.sequence)
                .then_with(|| left.reference.position.cmp(&right.reference.position))
                .then_with(|| left.reference.decoy.cmp(&right.reference.decoy))
        });

        let mut chunks = Vec::new();
        let mut chunk: Vec<DigestGroup> = Vec::new();
        for digest in digests {
            let sequence_changed = chunk
                .last()
                .is_some_and(|previous| previous.reference.sequence != digest.reference.sequence);
            if sequence_changed && chunk.len() >= chunk_size {
                chunks.push(std::mem::take(&mut chunk));
            }
            chunk.push(digest);
        }
        if !chunk.is_empty() {
            chunks.push(chunk);
        }
        chunks
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
        Self::reorder_peptides_with_reference(target_decoys, None);
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
                target_sequences.insert(peptide.sequence.clone());
            });

        let mut target_decoys = targets
            .into_par_iter()
            .flat_map_iter(|peptide| {
                if peptide.decoy {
                    return vec![peptide];
                }
                let decoy = peptide.reverse();
                if target_sequences.contains(&decoy.sequence[..]) {
                    vec![peptide]
                } else {
                    vec![decoy, peptide]
                }
            })
            .collect::<Vec<_>>();
        Self::reorder_peptides_with_reference(
            &mut target_decoys,
            self.label_reference().as_deref(),
        );
        target_decoys
    }

    fn reorder_peptides_with_reference(target_decoys: &mut Vec<Peptide>, reference: Option<&str>) {
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
                && chemical_modifications_eq(remove, keep)
                && (remove.modifications == keep.modifications
                    || channel_zero_provenance_eq(remove, keep))
            {
                if remove.label_channel != keep.label_channel {
                    let has_channel_site = keep
                        .applied_modifications()
                        .chain(remove.applied_modifications())
                        .any(|applied| applied.kind == ModificationKind::Label);
                    keep.label_channel = has_channel_site
                        .then(|| {
                            preferred_channel(
                                keep.label_channel.as_deref(),
                                remove.label_channel.as_deref(),
                                reference,
                            )
                        })
                        .flatten();
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
    /// `protein` and `decoy` are optional. Configured static, variable, and
    /// channel-aware modifications are applied.
    /// Decoys are generated by reversal when `self.generate_decoys` is true,
    /// subject to the same deduplication as the normal FASTA digest path.
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
                    sequence: seq.into(),
                    protein,
                    protein_start: None,
                    prev_aa: None,
                    next_aa: None,
                    missed_cleavages: 0,
                    position: Position::Full,
                };
                match Peptide::try_from(digest) {
                    Ok(p) => Some(p),
                    Err(e) => {
                        log::warn!("skipping peptide: {e}");
                        None
                    }
                }
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
            targets.insert(p.sequence.clone());
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

        Self::reorder_peptides_with_reference(&mut result, label_reference.as_deref());
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
        let mut fragments = target_decoys
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
            .collect::<Vec<_>>();
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
        IndexedDatabase {
            peptides: target_decoys,
            fragments,
            min_value,
            bucket_size: self.bucket_size,
            ion_kinds: self.ion_kinds,
            generate_decoys: self.generate_decoys,
            potential_mods,
            model_mods,
            label_reference,
            label_channels,
            decoy_tag: self.decoy_tag,
            decoy_pairing: Vec::new(),
        }
    }
}

fn with_estimation_margin(bytes: u64) -> u64 {
    // Parallel collection, allocator size classes, and sorting create overhead that
    // cannot be derived exactly from item counts. Use a conservative 50% margin.
    bytes.saturating_add(bytes / 2)
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
    /// Variable modification candidates used by PTM localization.
    pub potential_mods: Vec<(ModificationSpecificity, f32)>,
    /// Variable and precursor-label modifications used by property models.
    /// Label modifications are intentionally excluded from `potential_mods`
    /// because they are not PTM-localization candidates.
    pub model_mods: Vec<(ModificationSpecificity, f32)>,
    pub label_reference: Option<Arc<str>>,
    pub label_channels: Vec<Arc<str>>,
    pub bucket_size: usize,
    pub generate_decoys: bool,
    pub decoy_tag: String,
    /// Optional explicit target pairing for non-reversal decoy peptides.
    pub decoy_pairing: Vec<PeptideIx>,
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
#[path = "../tests/unit/database.rs"]
mod test;
