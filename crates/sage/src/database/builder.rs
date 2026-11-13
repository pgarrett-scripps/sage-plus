//! Search-space configuration: `Builder`, `Parameters` and their accessors.
//!
//! Moved out of `database.rs` unchanged; every item keeps its `database::` path.

use crate::cleavage::ValidatedCustomCleavageLibrary;
use crate::enzyme::{Enzyme, EnzymeParameters};
use crate::fasta::Fasta;
use crate::fragment_loss::{FragmentLossEntry, FragmentLosses};
use crate::ion_series::Kind;
use crate::modification::{
    validate_mods, validate_var_mods, ModificationDefinition, ModificationSpecificity, SearchMode,
    SiteMode, StaticModEntry, VarModEntry,
};
use crate::peptide::{LabelModificationCache, ModificationLookup, Peptide, VariableRule};
use crate::ptm_library::PtmLibrary;
use serde::{Deserialize, Serialize};
use std::collections::{BTreeMap, HashMap, HashSet};
use std::sync::Arc;

use super::*;

#[derive(Deserialize, Serialize, Clone, Debug, schemars::JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct EnzymeBuilder {
    /// How many missed cleavages to use, or `"unlimited"` to keep every
    /// peptide within the length limits. Required when a FASTA is digested.
    pub missed_cleavages: Option<Limit<u8>>,
    /// Minimum peptide length that will be fragmented
    #[schemars(range(min = 1))]
    pub min_len: Option<usize>,
    /// Maximum peptide length that will be fragmented, or `"unlimited"` for
    /// no length limit beyond the peptide mass range and the
    /// 255-residue encoding limit.
    pub max_len: Option<Limit<usize>>,
    /// Residues the enzyme cleaves at, `$` for no digestion (whole proteins),
    /// or an empty string for a non-specific digest. Required when a FASTA is
    /// digested.
    pub cleave_at: Option<String>,
    /// Residues that block cleavage when they follow a cleavage site (`""`
    /// for none). Required when a FASTA is digested.
    pub restrict: Option<String>,
    pub c_terminal: Option<bool>,
    /// Also generate semi-enzymatic peptides. Required when a FASTA is
    /// digested.
    pub semi_enzymatic: Option<bool>,
}

impl EnzymeBuilder {
    /// Enzyme fields every FASTA search must set.
    pub fn missing_required_fields(&self) -> Vec<&'static str> {
        [
            ("cleave_at", self.cleave_at.is_none()),
            ("restrict", self.restrict.is_none()),
            ("missed_cleavages", self.missed_cleavages.is_none()),
            ("semi_enzymatic", self.semi_enzymatic.is_none()),
        ]
        .into_iter()
        .filter_map(|(name, missing)| missing.then_some(name))
        .collect()
    }

    /// Effective maximum peptide length: `max_len`, or the 255-residue
    /// encoding limit when it is `"unlimited"`.
    pub fn effective_max_len(&self) -> usize {
        match self.max_len.unwrap_or(Limit::Limited(50)) {
            Limit::Limited(value) => value,
            Limit::Unlimited => crate::enzyme::MAX_PEPTIDE_LEN,
        }
    }
}

impl Default for EnzymeBuilder {
    fn default() -> Self {
        Self {
            missed_cleavages: Some(Limit::Limited(0)),
            min_len: Some(5),
            max_len: Some(Limit::Limited(50)),
            cleave_at: Some("KR".into()),
            restrict: Some("P".into()),
            c_terminal: Some(true),
            semi_enzymatic: Some(false),
        }
    }
}

impl From<EnzymeBuilder> for EnzymeParameters {
    fn from(en: EnzymeBuilder) -> EnzymeParameters {
        let max_len = en.effective_max_len();
        EnzymeParameters {
            clip_n_term_met: false,
            missed_cleavages: en.missed_cleavages.unwrap_or(Limit::Limited(1)).limited(),
            min_len: en.min_len.unwrap_or(5),
            max_len,
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
    /// Generic fragment neutral losses, keyed by a label, for example
    /// `{"Water": {"mass": 18.010565, "sites": ["S", "T", "E", "D"],
    /// "ion_kinds": ["b", "y"]}}`. Loss ions are matched only when a
    /// candidate is fully scored and enter rescoring as separate features.
    /// Off unless present. Separate from modification `neutral_losses`.
    pub fragment_losses: Option<BTreeMap<String, FragmentLossEntry>>,
    /// Most generic fragment losses stacked on one fragment (default 1).
    /// Requires `fragment_losses`.
    #[schemars(range(min = 1))]
    pub max_fragment_losses: Option<usize>,
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
    /// Path to a FASTA of common contaminant proteins, searched alongside
    /// `fasta`. Its proteins are ordinary targets, digested and decoyed the
    /// same way; what changes is that a PSM whose peptide belongs to any
    /// contaminant protein is reported as one, in `is_contaminant`. Unset
    /// (the default) leaves every output column as it was.
    pub contaminant_fasta: Option<String>,
    /// Path to a FASTA of entrapment proteins — a proteome that cannot be in
    /// the sample — searched alongside `fasta`. Its proteins are ordinary
    /// targets with their own decoys, so target-decoy competition stays fair;
    /// every discovery of one is false, which is what makes the false
    /// discovery proportion estimable. Unset (the default) leaves every
    /// output column as it was.
    pub entrapment_fasta: Option<String>,
    /// Expand ambiguous FASTA residues into every residue they may stand
    /// for: B to D or N, Z to E or Q, X to each of the 20 standard residues
    /// (default false: peptides containing B, X or Z are not searched).
    /// J (Ile or Leu) is always searched with the I/L mass.
    pub expand_ambiguous_residues: Option<bool>,
    /// Peptides whose ambiguous residues expand into more sequences than
    /// this are dropped (default 20, one X per peptide). Values below 1 are
    /// normalized to 1.
    pub max_ambiguous_variants: Option<usize>,
    /// Merge peptides that differ only in I, L and J (Ile or Leu) into one
    /// peptide listing every protein, when their modifications and decoy
    /// flag also match (default true). They have the same mass and fragments.
    pub merge_isoleucine_leucine: Option<bool>,
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

    /// Reject peptide length and mass ranges that cannot contain a peptide.
    pub fn validate_peptide_ranges(&self) -> Result<(), String> {
        let enzyme = self.enzyme.clone().unwrap_or_default();
        let (min_len, max_len) = (enzyme.min_len.unwrap_or(5), enzyme.effective_max_len());
        if min_len > max_len {
            return Err(format!(
                "`database.enzyme.min_len` ({min_len}) is greater than `max_len` ({max_len})"
            ));
        }
        let min_mass = self.peptide_min_mass.unwrap_or(500.0);
        let max_mass = self.peptide_max_mass.unwrap_or(5000.0);
        if min_mass
            .partial_cmp(&max_mass)
            .is_none_or(|order| order.is_gt())
        {
            return Err(format!(
                "`database.peptide_min_mass` ({min_mass}) must not exceed \
                 `peptide_max_mass` ({max_mass})"
            ));
        }
        Ok(())
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

    /// Validate `fragment_losses` and `max_fragment_losses`.
    pub fn validate_fragment_losses(&self) -> Result<(), String> {
        let ion_kinds = self.ion_kinds.clone().unwrap_or(vec![Kind::B, Kind::Y]);
        crate::fragment_loss::resolve(
            self.fragment_losses.as_ref(),
            self.max_fragment_losses,
            &ion_kinds,
        )
        .map(|_| ())
    }

    /// Resolve the configured database settings. Fails only on a
    /// modification key that is not a site specificity, which is a
    /// configuration error; see
    /// [`crate::modification::validate_mods`]. Every other setting has a
    /// default.
    pub fn make_parameters(self) -> Result<Parameters, String> {
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
        Ok(Parameters {
            bucket_size,
            peptide_min_mass: self.peptide_min_mass.unwrap_or(500.0),
            peptide_max_mass: self.peptide_max_mass.unwrap_or(5000.0),
            ion_kinds: self.ion_kinds.unwrap_or(vec![Kind::B, Kind::Y]),
            min_ion_index: self.min_ion_index.unwrap_or(2),
            max_fragment_losses: self.fragment_losses.as_ref().map(|_| {
                self.max_fragment_losses
                    .unwrap_or(crate::fragment_loss::DEFAULT_MAX_FRAGMENT_LOSSES)
            }),
            fragment_losses: self.fragment_losses,
            decoy_tag: self.decoy_tag.unwrap_or_else(|| "rev_".into()),
            enzyme: self.enzyme.unwrap_or_default(),
            static_mods: validate_mods(self.static_mods)?,
            variable_mods: validate_var_mods(self.variable_mods)?,
            max_variable_mods,
            max_combinations: self.max_combinations.map(|x| x.max(1)),
            max_total_variable_mods,
            ptm_library: self.ptm_library,
            generate_decoys: self.generate_decoys.unwrap_or(true),
            clip_n_term_met: self.clip_n_term_met.unwrap_or(true),
            fasta: self.fasta.unwrap_or_default(),
            contaminant_fasta: self.contaminant_fasta,
            entrapment_fasta: self.entrapment_fasta,
            expand_ambiguous_residues: self.expand_ambiguous_residues.unwrap_or(false),
            max_ambiguous_variants: self.max_ambiguous_variants.unwrap_or(20).max(1),
            merge_isoleucine_leucine: self.merge_isoleucine_leucine.unwrap_or(true),
            peptides: self.peptides,
            custom_cleavage_sites: self.custom_cleavage_sites,
            prefilter: self.prefilter.unwrap_or(false),
            prefilter_min_matched_peaks: self.prefilter_min_matched_peaks.unwrap_or(4).max(1),
            prefilter_max_peaks: self.prefilter_max_peaks.filter(|&n| n > 0),
            loaded_ptm_library: None,
        })
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
    #[serde(skip_serializing_if = "Option::is_none")]
    pub fragment_losses: Option<BTreeMap<String, FragmentLossEntry>>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub max_fragment_losses: Option<usize>,
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
    /// Path to a FASTA of contaminant proteins, searched alongside `fasta`.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub contaminant_fasta: Option<String>,
    /// Path to a FASTA of entrapment proteins, searched alongside `fasta`.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub entrapment_fasta: Option<String>,
    pub expand_ambiguous_residues: bool,
    pub max_ambiguous_variants: usize,
    pub merge_isoleucine_leucine: bool,
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

    /// How FASTA proteins with B, X or Z are searched: an info message when
    /// `expand_ambiguous_residues` expands their peptides, a warning when
    /// they are dropped, and `None` when there are no such proteins.
    pub fn ambiguous_proteins_message(&self, fasta: &Fasta) -> Option<(log::Level, String)> {
        let proteins = fasta.ambiguous_protein_count();
        if proteins == 0 {
            return None;
        }
        let prefix = format!(
            "{proteins} FASTA protein(s) contain residues without a single mass (B, X or Z)"
        );
        Some(match self.expand_ambiguous_residues {
            true => (
                log::Level::Info,
                format!(
                    "{prefix}; peptides containing them are expanded into up to {} variant(s) each (database.max_ambiguous_variants)",
                    self.max_ambiguous_variants
                ),
            ),
            false => (
                log::Level::Warn,
                format!(
                    "{prefix}; peptides containing them are not searched unless database.expand_ambiguous_residues is true"
                ),
            ),
        })
    }

    /// Log [`Self::ambiguous_proteins_message`].
    pub fn log_ambiguous_proteins(&self, fasta: &Fasta) {
        if let Some((level, message)) = self.ambiguous_proteins_message(fasta) {
            log::log!(level, "{message}");
        }
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

    /// Warning for a digest whose peptide length is bounded only by peptide
    /// mass: `max_len` is `"unlimited"` together with unlimited missed
    /// cleavages or a non-specific digest.
    pub fn unbounded_length_warning(&self) -> Option<String> {
        let enzyme = &self.enzyme;
        let unbounded_sites = enzyme.missed_cleavages == Some(Limit::Unlimited)
            || enzyme.cleave_at.as_deref() == Some("");
        (enzyme.max_len == Some(Limit::Unlimited) && unbounded_sites).then(|| {
            format!(
                "database.enzyme.max_len is \"unlimited\" with {}: only the peptide mass range \
                 and the {}-residue encoding limit bound peptide length, so every residue can \
                 start hundreds of peptides and the digest can be very large",
                match enzyme.cleave_at.as_deref() {
                    Some("") => "a non-specific digest",
                    _ => "unlimited missed cleavages",
                },
                crate::enzyme::MAX_PEPTIDE_LEN
            )
        })
    }

    pub fn validate_compact_modifications(&self) -> Result<(), String> {
        let max_len = self.enzyme.effective_max_len();
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

    pub(super) fn label_channels(&self) -> Vec<Arc<str>> {
        self.channel_definitions()
            .first()
            .map(|definition| definition.channel_offsets.keys().cloned().collect())
            .unwrap_or_default()
    }

    pub(super) fn label_reference(&self) -> Option<Arc<str>> {
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
    pub(super) fn variable_modifications(&self) -> Vec<VariableRule> {
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

    /// Validated generic fragment losses, or `None` when not configured.
    /// Panics on an invalid configuration; the CLI validates it first with
    /// [`Builder::validate_fragment_losses`].
    pub fn fragment_losses(&self) -> Option<FragmentLosses> {
        crate::fragment_loss::resolve(
            self.fragment_losses.as_ref(),
            self.max_fragment_losses,
            &self.ion_kinds,
        )
        .unwrap_or_else(|error| panic!("{error}"))
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
        // `variable_mods` is a HashMap, so sort the rules: the localizer numbers
        // modifications (ProForma `#g1`, `#g2`, site rows) in rule order, which
        // must not change from run to run.
        let mut rules = self
            .variable_mods
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
            .collect::<Vec<_>>();
        rules.sort_by(|a, b| {
            a.specificity
                .cmp(&b.specificity)
                .then_with(|| a.definition.cmp(&b.definition))
        });
        rules
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

    pub(super) fn static_modifications(
        &self,
    ) -> HashMap<ModificationSpecificity, Arc<ModificationDefinition>> {
        self.static_mods
            .iter()
            .map(|(specificity, entry)| (*specificity, Arc::new(entry.definition())))
            .collect()
    }
}
