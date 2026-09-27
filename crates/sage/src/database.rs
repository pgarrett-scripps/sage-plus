use crate::cleavage::ValidatedCustomCleavageLibrary;
use crate::enzyme::{
    group_protein_digests, Digest, DigestGroup, Enzyme, EnzymeParameters, Position,
    ProteinOccurrence,
};
use crate::fasta::Fasta;
use crate::ion_series::{IonGroupSeries, Kind};
use crate::mass::Tolerance;
use crate::modification::{
    validate_mods, validate_var_mods, ModificationDefinition, ModificationSpecificity, SearchMode,
    SiteMode, StaticModEntry, VarModEntry,
};
use crate::peptide::{
    AppliedModification, LabelModificationCache, LibrarySite, ModificationKind, ModificationLookup,
    ModificationPlan, Peptide, Site, VariableRule, INLINE_PROTEINS,
};
use crate::ptm_library::PtmLibrary;
use crate::scoring::Feature;
use crate::sequence::PeptideSequence;
use dashmap::DashSet;
use fnv::FnvBuildHasher;
use rayon::prelude::*;
use serde::{Deserialize, Serialize};
use std::borrow::Cow;
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
            clip_n_term_met: false,
            missed_cleavages: en.missed_cleavages.unwrap_or(1),
            min_len: en.min_len.unwrap_or(5),
            max_len: en.max_len.unwrap_or(50),
            enzyme: Enzyme::new(
                &en.cleave_at.unwrap_or_else(|| "KR".into()),
                &en.restrict.unwrap_or_else(|| "".into()),
                en.c_terminal.unwrap_or(true),
                en.semi_enzymatic.unwrap_or(false),
            ),
            ambiguous_variants: None,
        }
    }
}

#[derive(Deserialize, Default, Clone, schemars::JsonSchema)]
#[serde(deny_unknown_fields)]
/// Parameters used for generating the fragment database
pub struct Builder {
    /// Maximum number of theoretical fragments stored in one search bucket.
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
    /// Named static definitions with mass and explicit sites.
    /// Upstream Sage symbol-keyed masses (`{"C": 57.021464}`) remain readable.
    #[serde(default, deserialize_with = "crate::modification::deserialize_mod_map")]
    #[schemars(with = "Option<crate::modification::StaticModConfig>")]
    pub static_mods: Option<HashMap<String, StaticModEntry>>,
    /// Named variable definitions with mass, explicit sites, and optional
    /// per-modification limits, library policy, and fragment behavior.
    /// Upstream Sage symbol-keyed masses (`{"M": [15.9949]}`) remain readable.
    #[serde(default, deserialize_with = "crate::modification::deserialize_mod_map")]
    #[schemars(with = "Option<crate::modification::VariableModConfig>")]
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
    /// Also search each protein with its initiator methionine removed when
    /// the second residue is G, A, S, T, C, P, or V (MetAP clipping). The
    /// clipped peptides are protein N-terminal. Enzymatic digests only.
    /// Default: true.
    pub clip_n_term_met: Option<bool>,
    /// Path to fasta database
    pub fasta: Option<String>,
    /// Expand ambiguous FASTA residues into every residue they may stand
    /// for: B to D or N, Z to E or Q, X to each of the 20 standard residues
    /// (default false: peptides containing B, X or Z are not searched).
    /// J (Ile or Leu) is always searched with the I/L mass.
    pub expand_ambiguous_residues: Option<bool>,
    /// Peptides whose ambiguous residues expand into more sequences than
    /// this are dropped (default 20, one X per peptide). Values below 1 are
    /// normalized to 1.
    pub max_ambiguous_variants: Option<usize>,
    /// Path to a pre-digested peptide TSV file (additive with `fasta`).
    /// Required column: `sequence`. Optional columns: `protein`, `decoy`.
    /// Configured static, variable, and channel-aware modifications are applied.
    pub peptides: Option<String>,
    /// Path to a protein-specific custom cleavage-site TSV or Parquet file.
    /// Required columns: `protein`, `position`; optional column: `context`.
    pub custom_cleavage_sites: Option<String>,
    /// Deprecated and ignored. The prefilter streams proteins one at a time.
    pub prefilter_chunk_size: Option<usize>,
    /// Pre-filter the database to minimize memory usage
    pub prefilter: Option<bool>,
    /// Preliminary fragment matches one precursor hypothesis of a spectrum
    /// needs for the prefilter to keep a peptide (default 4).
    pub prefilter_min_matched_peaks: Option<u16>,
    /// Only each spectrum's most intense peaks are used by the prefilter
    /// (default: every processed peak).
    pub prefilter_max_peaks: Option<usize>,
    /// Deprecated compatibility option. Exact prefiltering always uses compact
    /// survivor tracking and ignores this value.
    pub prefilter_low_memory: Option<bool>,
}

impl Builder {
    /// Reject enzyme residues that would otherwise abort database building.
    pub fn validate_enzyme(&self) -> Result<(), String> {
        let Some(enzyme) = &self.enzyme else {
            return Ok(());
        };
        Enzyme::validate_residues(
            enzyme.cleave_at.as_deref().unwrap_or("KR"),
            enzyme.restrict.as_deref().unwrap_or(""),
        )
    }

    pub fn validate_modification_keys(&self) -> Result<(), String> {
        for key in self
            .static_mods
            .iter()
            .flat_map(|mods| mods.keys())
            .chain(self.variable_mods.iter().flat_map(|mods| mods.keys()))
        {
            key.parse::<ModificationSpecificity>()
                .map_err(|_| format!("invalid modification key `{key}`"))?;
        }
        Ok(())
    }

    pub fn make_parameters(self) -> Parameters {
        if self.prefilter_low_memory.is_some() {
            log::warn!("database.prefilter_low_memory is deprecated and ignored");
        }
        if self.prefilter_chunk_size.is_some() {
            log::warn!("database.prefilter_chunk_size is deprecated and ignored");
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
            clip_n_term_met: self.clip_n_term_met.unwrap_or(true),
            fasta: self.fasta.unwrap_or_default(),
            expand_ambiguous_residues: self.expand_ambiguous_residues.unwrap_or(false),
            max_ambiguous_variants: self.max_ambiguous_variants.unwrap_or(20).max(1),
            peptides: self.peptides,
            custom_cleavage_sites: self.custom_cleavage_sites,
            prefilter: self.prefilter.unwrap_or(false),
            prefilter_min_matched_peaks: self.prefilter_min_matched_peaks.unwrap_or(4).max(1),
            prefilter_max_peaks: self.prefilter_max_peaks.filter(|&n| n > 0),
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
    #[serde(serialize_with = "crate::modification::serialize_static_mods")]
    pub static_mods: HashMap<ModificationSpecificity, StaticModEntry>,
    #[serde(serialize_with = "crate::modification::serialize_variable_mods")]
    pub variable_mods: HashMap<ModificationSpecificity, Vec<VarModEntry>>,
    pub max_variable_mods: usize,
    pub max_combinations: Option<usize>,
    pub max_total_variable_mods: usize,
    pub ptm_library: Option<PtmLibrarySettings>,
    pub decoy_tag: String,
    pub generate_decoys: bool,
    pub clip_n_term_met: bool,
    pub fasta: String,
    pub expand_ambiguous_residues: bool,
    pub max_ambiguous_variants: usize,
    pub peptides: Option<String>,
    pub custom_cleavage_sites: Option<String>,
    pub prefilter: bool,
    pub prefilter_min_matched_peaks: u16,
    pub prefilter_max_peaks: Option<usize>,
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

/// Hashes of a peptide sequence and of its internal reversal, the sequence
/// of the decoy generated from it. Equal sequences hash equally; unequal
/// sequences rarely do.
pub fn sequence_hashes(sequence: &[u8]) -> (u64, u64) {
    const OFFSET: u64 = 0xcbf2_9ce4_8422_2325;
    const PRIME: u64 = 0x0000_0100_0000_01b3;
    let step = |hash: u64, byte: &u8| (hash ^ u64::from(*byte)).wrapping_mul(PRIME);
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

impl Parameters {
    /// Digest settings, including initiator methionine clipping and
    /// ambiguous-residue expansion.
    pub fn enzyme_parameters(&self) -> EnzymeParameters {
        let mut enzyme: EnzymeParameters = self.enzyme.clone().into();
        enzyme.clip_n_term_met = self.clip_n_term_met;
        enzyme.ambiguous_variants = self
            .expand_ambiguous_residues
            .then_some(self.max_ambiguous_variants);
        enzyme
    }

    /// Log how many digests with ambiguous residues were expanded or
    /// dropped. Only proteins with such residues are digested again.
    pub fn log_ambiguous_expansion(
        &self,
        fasta: &Fasta,
        custom_cleavages: Option<&ValidatedCustomCleavageLibrary>,
    ) {
        if let Some(summary) =
            fasta.ambiguous_expansion_summary(&self.enzyme_parameters(), custom_cleavages)
        {
            log::info!(
                "expanded {} digest(s) with ambiguous residues (B, X, Z) into {} variant(s); dropped {} with more than {} variant(s) (database.max_ambiguous_variants)",
                summary.expanded,
                summary.variants,
                summary.dropped,
                self.max_ambiguous_variants,
            );
        }
    }

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
                entries
                    .iter()
                    .enumerate()
                    .filter(|(_, entry)| entry.search_mode() == SearchMode::Database)
                    .map(|(entry_order, entry)| {
                        (
                            *specificity,
                            entry_order,
                            Arc::new(entry.definition()),
                            entry.max_count(),
                            entry.max_total_count(),
                            entry.site_mode(),
                        )
                    })
            })
            .collect::<Vec<_>>();
        mods.sort_unstable_by(|a, b| a.0.cmp(&b.0).then_with(|| a.1.cmp(&b.1)));
        let mut named_groups: HashMap<Arc<str>, usize> = HashMap::new();
        let mut next_group = 0usize;
        mods.into_iter()
            .map(
                |(specificity, _, modification, max_count, max_total_count, site_mode)| {
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
                        max_total_count,
                        site_mode,
                        count_group,
                    }
                },
            )
            .collect()
    }

    /// Group search-time mass offsets by chemical definition. Each group keeps
    /// every configured specificity so placement and localization agree.
    pub fn mass_offset_modifications(&self) -> Vec<MassOffset> {
        let mut groups: Vec<MassOffset> = Vec::new();
        let mut entries = self
            .variable_mods
            .iter()
            .flat_map(|(specificity, entries)| {
                entries
                    .iter()
                    .filter(|entry| entry.search_mode() == SearchMode::MassOffset)
                    .map(move |entry| (*specificity, entry.definition(), entry.site_mode()))
            })
            .collect::<Vec<_>>();
        entries.sort_unstable_by(|left, right| {
            left.1.cmp(&right.1).then_with(|| left.0.cmp(&right.0))
        });
        for (specificity, definition, site_mode) in entries {
            match groups
                .iter_mut()
                .find(|group| *group.definition == definition)
            {
                Some(group) => {
                    if !group.specificities.contains(&specificity) {
                        group.specificities.push(specificity);
                    }
                }
                None => groups.push(MassOffset {
                    definition: Arc::new(definition),
                    specificities: vec![specificity],
                    site_mode,
                }),
            }
        }
        groups
    }

    /// Resolve localization permissions using the same typed library evidence as search.
    pub fn localization_rules(
        &self,
        peptide: &Peptide,
    ) -> Vec<crate::ptm::ResolvedLocalizationRule> {
        self.variable_mods
            .iter()
            .flat_map(|(specificity, entries)| {
                entries.iter().map(move |entry| {
                    let definition = Arc::new(entry.definition());
                    let mut sites = peptide.rule_sites(*specificity);
                    if entry.site_mode() == SiteMode::Library {
                        sites.retain(|candidate| {
                            self.loaded_ptm_library.as_ref().is_some_and(|library| {
                                peptide.protein_sites.iter().any(|occurrence| {
                                    let Some(start) = occurrence.start else {
                                        return false;
                                    };
                                    library.sites_for(&occurrence.protein).iter().any(|record| {
                                        let Some(index) = record.position.checked_sub(start) else {
                                            return false;
                                        };
                                        definition.name.as_deref()
                                            == Some(record.modification.as_ref())
                                            && peptide.sequence.get(index as usize)
                                                == Some(&record.residue)
                                            && record.attachment.site(
                                                index,
                                                peptide.sequence.len(),
                                                peptide.position,
                                            ) == Some(*candidate)
                                    })
                                })
                            })
                        });
                    }
                    crate::ptm::ResolvedLocalizationRule {
                        specificity: *specificity,
                        definition,
                        sites,
                    }
                })
            })
            .collect()
    }

    pub fn validate_ptm_library(&self, library: &PtmLibrary) -> Result<(), String> {
        let static_rules = self.static_mods.iter().collect::<Vec<_>>();
        for (index, (left, left_entry)) in static_rules.iter().enumerate() {
            let left_definition = left_entry.definition();
            for (right, right_entry) in &static_rules[index + 1..] {
                let right_definition = right_entry.definition();
                if left_definition != right_definition
                    && (left.overlaps(**right)
                        || (left_definition.name.is_some()
                            && left_definition.name == right_definition.name))
                {
                    return Err(format!("conflicting static modifications at `{}` and `{}`. A physical attachment can have only one fixed definition", left.explicit_name(), right.explicit_name()));
                }
            }
        }
        let static_names = self
            .static_mods
            .values()
            .filter_map(|entry| entry.definition().name)
            .collect::<HashSet<_>>();
        if let Some(name) = self
            .variable_mods
            .values()
            .flatten()
            .filter_map(|entry| entry.definition().name)
            .find(|name| static_names.contains(name))
        {
            return Err(format!("modification `{name}` is defined in both static_mods and variable_mods. Use distinct IDs for different application policies"));
        }
        if self.max_total_variable_mods < self.max_variable_mods {
            return Err(
                "database.max_total_variable_mods must be at least database.max_variable_mods"
                    .into(),
            );
        }

        let mut offset_names: HashMap<Arc<str>, (ModificationDefinition, SiteMode)> =
            HashMap::new();
        let mut offset_count = 0usize;
        for entries in self.variable_mods.values() {
            for entry in entries
                .iter()
                .filter(|entry| entry.search_mode() == SearchMode::MassOffset)
            {
                offset_count += 1;
                let definition = entry.definition();
                if entry.site_mode() != SiteMode::Exhaustive && definition.name.is_none() {
                    return Err(
                        "mass_offset modifications using `library` or `both` require `name`".into(),
                    );
                }
                if let Some(name) = definition.name.clone() {
                    match offset_names.get(&name) {
                        Some((existing, site_mode))
                            if *existing != definition || *site_mode != entry.site_mode() =>
                        {
                            return Err(format!(
                                "variable modification `{name}` has inconsistent definitions across specificities"
                            ));
                        }
                        Some(_) => {}
                        None => {
                            offset_names.insert(name, (definition, entry.site_mode()));
                        }
                    }
                }
            }
        }
        if self.mass_offset_modifications().len() > MAX_MASS_OFFSETS {
            return Err(format!(
                "at most {MAX_MASS_OFFSETS} distinct mass_offset modifications are supported"
            ));
        }
        if offset_count > 0 && self.label_channels().len() > 1 {
            return Err(
                "mass_offset modifications cannot be combined with channel-aware labels".into(),
            );
        }

        let rules = self.variable_modifications();
        type Definition<'a> = (
            &'a ModificationDefinition,
            Option<usize>,
            Option<usize>,
            SiteMode,
        );
        let mut definitions: HashMap<&str, Definition> = HashMap::new();
        for rule in &rules {
            if self.ptm_library.is_some() && rule.total_limit().is_none() {
                return Err(
                    "all variable modifications require `max_count` or `max_total_count` when database.ptm_library is configured"
                        .into(),
                );
            }
            if rule.site_mode != SiteMode::Exhaustive
                && (rule.modification.name.is_none() || rule.total_limit().is_none())
            {
                return Err(
                    "variable modifications using `library` or `both` require `name` and `max_count` or `max_total_count`"
                        .into(),
                );
            }
            if let Some(name) = rule.modification.name.as_deref() {
                if offset_names.contains_key(name) {
                    return Err(format!(
                        "variable modification `{name}` cannot use both `database` and `mass_offset` search modes"
                    ));
                }
                if let Some((definition, max_count, max_total_count, site_mode)) =
                    definitions.get(name)
                {
                    if *definition != rule.modification.as_ref()
                        || *max_count != rule.max_count
                        || *max_total_count != rule.max_total_count
                        || *site_mode != rule.site_mode
                    {
                        return Err(format!(
                            "variable modification `{name}` has inconsistent definitions across specificities"
                        ));
                    }
                } else {
                    definitions.insert(
                        name,
                        (
                            &rule.modification,
                            rule.max_count,
                            rule.max_total_count,
                            rule.site_mode,
                        ),
                    );
                }
            }
        }

        for site in library.iter() {
            if let Some((_, site_mode)) = offset_names.get(site.modification.as_ref()) {
                if *site_mode == SiteMode::Exhaustive {
                    return Err(format!(
                        "PTM library modification `{}` must use site_mode `library` or `both`",
                        site.modification
                    ));
                }
                continue;
            }
            match definitions.get(site.modification.as_ref()) {
                None => {
                    return Err(format!(
                        "PTM library references undefined modification `{}`",
                        site.modification
                    ))
                }
                Some((_, _, _, SiteMode::Exhaustive)) => {
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

    fn variable_variant_count(&self, digest: &Digest, origins: &[ProteinOccurrence]) -> u64 {
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
        let mut candidates: HashMap<(usize, usize), bool> = HashMap::new();
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
        let digests = fasta.digest_with_custom_cleavages(&enzyme, custom_cleavages);

        log::trace!("grouping digests");
        let start_num = digests.len();
        let digests = group_protein_digests(digests);
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
        self.modify_digest_chunks(digests, target_sequences, MODIFY_DIGEST_CHUNK_GROUPS)
    }

    fn modify_digest_chunks(
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
        self.reorder_peptides_with_labels(&mut target_decoys);
        target_decoys
    }

    /// Sort and deduplicate peptides as [`Self::modify_digests`] does,
    /// merging label channels with this database's reference channel.
    pub fn reorder_peptides_with_labels(&self, target_decoys: &mut Vec<Peptide>) {
        Self::reorder_peptides_with_reference(target_decoys, self.label_reference().as_deref());
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
        target_decoys.par_sort_unstable_by(|a, b| {
            a.monoisotopic
                .total_cmp(&b.monoisotopic)
                .then_with(|| a.initial_sort(b))
                .then_with(|| from_tsv(a).cmp(&from_tsv(b)))
                .then(a.semi_enzymatic.cmp(&b.semi_enzymatic))
                .then(a.missed_cleavages.cmp(&b.missed_cleavages))
                .then(a.position.cmp(&b.position))
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
                    expanded_from: None,
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
        let (compressed_fragments, min_value) = FragmentIndex::build(&self, &target_decoys);
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
            ion_kinds: self.ion_kinds,
            generate_decoys: self.generate_decoys,
            potential_mods,
            localization_mods,
            model_mods,
            label_reference,
            label_channels,
            decoy_tag: self.decoy_tag,
            decoy_pairing: Vec::new(),
        }
    }
}

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
        Parameters::reorder_peptides_with_reference(target_decoys, self.label_reference.as_deref());
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
                        || !target_sequences
                            .is_some_and(|targets| targets.contains(&(peptide.sequence[..])))
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

/// Upper bound on distinct search-time offsets; candidate hypotheses encode
/// the offset in one byte, with zero reserved for "no offset".
pub const MAX_MASS_OFFSETS: usize = u8::MAX as usize - 1;

/// A modification tested at search time instead of being expanded into the
/// fragment index. At most one offset copy is placed on a peptide.
#[derive(Clone, Debug)]
pub struct MassOffset {
    pub definition: Arc<ModificationDefinition>,
    /// Every configured site rule for this definition.
    pub specificities: Vec<ModificationSpecificity>,
    pub site_mode: SiteMode,
}

impl MassOffset {
    pub fn mass(&self) -> f32 {
        self.definition.mass
    }

    /// Mass difference of the first generated fragment form containing the
    /// modification. This mirrors the preliminary fragment index, which keeps
    /// the first variant of every ion group; ion series sort losses ascending,
    /// so a required loss contributes its smallest configured mass.
    pub fn fragment_shift(&self) -> f32 {
        match self.definition.neutral_loss_mode {
            crate::modification::NeutralLossMode::Required => {
                let smallest = self
                    .definition
                    .neutral_losses
                    .iter()
                    .copied()
                    .min_by(f32::total_cmp)
                    .unwrap_or_default();
                self.definition.mass - smallest
            }
            crate::modification::NeutralLossMode::Optional => self.definition.mass,
        }
    }
}

/// The offset placement selected for a candidate. `offset` indexes
/// [`IndexedDatabase::mass_offsets`].
#[derive(Copy, Clone, Debug, PartialEq, Eq, Hash)]
pub struct MassOffsetAssignment {
    pub offset: u16,
    pub site: Site,
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

/// Lower `f32` mass bits stored in each packed fragment. The remaining upper
/// bits are shared by every search bucket with the same mass prefix.
const FRAGMENT_MASS_SUFFIX_BITS: u32 = 12;

#[repr(C, packed)]
#[derive(Copy, Clone)]
struct PackedFragment {
    peptide_index: u32,
    mass_suffix: u16,
}

impl PackedFragment {
    #[inline(always)]
    fn peptide_index(self) -> u32 {
        self.peptide_index
    }

    #[inline(always)]
    fn decode(self, mass_prefix: u32) -> Theoretical {
        Theoretical {
            peptide_index: PeptideIx(self.peptide_index),
            fragment_mz: f32::from_bits(mass_prefix | u32::from(self.mass_suffix)),
        }
    }
}

#[derive(Copy, Clone)]
struct FragmentBucket {
    mass_prefix: u32,
    start: u32,
    end: u32,
}

#[derive(Default)]
/// Lossless theoretical-fragment index using exact packed mass-prefix buckets.
pub struct FragmentIndex {
    buckets: Vec<FragmentBucket>,
    fragments: Vec<PackedFragment>,
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
    fn build(parameters: &Parameters, peptides: &[Peptide]) -> (Self, Vec<f32>) {
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
                for peptide in &peptides[range.clone()] {
                    for mass in preliminary_fragment_masses(parameters, peptide) {
                        let prefix = mass.to_bits() >> suffix_bits;
                        *counts.entry(prefix).or_default() += 1;
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
                    start: u32::try_from(bucket_start)
                        .expect("fragment index exceeds 32-bit offsets"),
                    end: u32::try_from(bucket_end).expect("fragment index exceeds 32-bit offsets"),
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
                for (peptide_index, peptide) in peptides[range.clone()].iter().enumerate() {
                    let peptide_index = range.start + peptide_index;
                    for mass in preliminary_fragment_masses(parameters, peptide) {
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

    fn bucket_search(
        &self,
        bucket_index: usize,
        peptide_lo: u32,
        peptide_hi: u32,
    ) -> FragmentIter<'_> {
        let bucket = self.buckets[bucket_index];
        let fragments = &self.fragments[bucket.start as usize..bucket.end as usize];
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
    fn resolve_bucket_ranges(
        &self,
        buckets: &[u32],
        peptide_lo: u32,
        peptide_hi: u32,
        out: &mut Vec<(u32, u32)>,
    ) {
        /// Buckets resolved together. Two searches each, sized against the
        /// core's outstanding-miss capacity rather than a vector width.
        const LANES: usize = 16;
        const SLOTS: usize = 2 * LANES;

        out.clear();
        out.reserve(buckets.len());
        if self.fragments.is_empty() {
            out.extend(buckets.iter().map(|&bucket| {
                let start = self.buckets[bucket as usize].start;
                (start, start)
            }));
            return;
        }
        let fragments = self.fragments.as_slice();

        for chunk in buckets.chunks(LANES) {
            // Even slots search `id < peptide_lo`, odd slots `id <= peptide_hi`,
            // exactly the two `partition_point` predicates of `bucket_search`.
            let mut base = [0usize; SLOTS];
            let mut size = [0usize; SLOTS];
            let mut longest = 0usize;
            for (lane, &bucket) in chunk.iter().enumerate() {
                let bucket = self.buckets[bucket as usize];
                let (start, len) = (bucket.start as usize, (bucket.end - bucket.start) as usize);
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
                let bucket = self.buckets[bucket as usize];
                if bucket.start == bucket.end {
                    out.push((bucket.start, bucket.start));
                    continue;
                }
                let (lo, hi) = (base[2 * lane], base[2 * lane + 1]);
                let first = lo + usize::from(fragments[lo].peptide_index() < peptide_lo);
                let last = hi + usize::from(fragments[hi].peptide_index() <= peptide_hi);
                out.push((first as u32, last as u32));
            }
        }
    }

    #[inline(always)]
    fn range_iter(&self, bucket: usize, first: u32, last: u32) -> FragmentIter<'_> {
        FragmentIter {
            fragments: &self.fragments[first as usize..last as usize],
            mass_prefix: self.buckets[bucket].mass_prefix,
            next: 0,
        }
    }
}

/// Per-thread scratch for [`IndexedQuery::page_search_batch`], reused across
/// spectra so the batched path allocates nothing in steady state.
#[derive(Default)]
struct BatchSearchScratch {
    /// Per window: fragment bounds and its span in `entries`.
    windows: Vec<(f32, f32, u32, u32)>,
    /// One bucket per (window, bucket) pair, in emission order.
    entries: Vec<u32>,
    /// `(bucket, entry)` pairs sorted by bucket for de-duplication.
    order: Vec<(u32, u32)>,
    /// Distinct buckets and their resolved ranges.
    distinct: Vec<u32>,
    resolved: Vec<(u32, u32)>,
    /// Resolved range per entry.
    ranges: Vec<(u32, u32)>,
}

thread_local! {
    static BATCH_SEARCH_SCRATCH: std::cell::RefCell<BatchSearchScratch> =
        std::cell::RefCell::new(BatchSearchScratch::default());
}

struct FragmentIter<'a> {
    fragments: &'a [PackedFragment],
    mass_prefix: u32,
    next: usize,
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

/// Chemical peptidoform identity used to merge equivalent offset and indexed
/// hypotheses.
pub fn same_peptidoform(left: &Peptide, right: &Peptide) -> bool {
    left.decoy == right.decoy
        && left.sequence == right.sequence
        && left.label_channel == right.label_channel
        && chemical_modifications_eq(left, right)
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
        mut on_match: impl FnMut(Theoretical),
    ) {
        let (precursor_lo, precursor_hi) = self.precursor_tol.bounds(self.precursor_mass);
        let peptide_lo = self.pre_idx_lo.min(u32::MAX as usize) as u32;
        let peptide_hi = self.pre_idx_hi.min(u32::MAX as usize) as u32;

        BATCH_SEARCH_SCRATCH.with(|scratch| {
            let mut scratch = scratch.borrow_mut();
            let BatchSearchScratch {
                windows,
                entries,
                order,
                distinct,
                resolved,
                ranges,
            } = &mut *scratch;
            windows.clear();
            entries.clear();

            let mut push_window = |fragment_lo: f32, fragment_hi: f32| {
                let (left, right) =
                    fragment_bucket_slice(&self.db.min_value, fragment_lo, fragment_hi);
                let first = entries.len() as u32;
                entries.extend(left as u32..right as u32);
                windows.push((fragment_lo, fragment_hi, first, entries.len() as u32));
            };
            for mass in masses {
                let (fragment_lo, fragment_hi) = self.fragment_tol.bounds(mass);
                push_window(fragment_lo, fragment_hi);
                for &shift in shifts {
                    push_window(fragment_lo - shift, fragment_hi - shift);
                }
            }

            // Resolve each distinct bucket once.
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
            self.db
                .fragments
                .resolve_bucket_ranges(distinct, peptide_lo, peptide_hi, resolved);
            ranges.clear();
            ranges.resize(entries.len(), (0, 0));
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
                    for frag in self
                        .db
                        .fragments
                        .range_iter(entries[entry] as usize, lo, hi)
                    {
                        if self.accepts(&frag, fragment_lo, fragment_hi, precursor_lo, precursor_hi)
                        {
                            on_match(frag);
                        }
                    }
                }
            }
        });
    }
}

fn fragment_bucket_slice(min_values: &[f32], low: f32, high: f32) -> (usize, usize) {
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
    let mut choices_by_site: HashMap<usize, (u64, u64)> = HashMap::new();
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

#[cfg(test)]
#[path = "../tests/unit/database.rs"]
mod test;
