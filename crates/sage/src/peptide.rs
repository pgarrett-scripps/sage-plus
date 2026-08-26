use std::cmp::Ordering;
use std::{
    collections::{BTreeMap, HashMap},
    fmt::Debug,
    sync::Arc,
};

use crate::modification::{ModificationDefinition, ModificationSpecificity, SiteMode};
use crate::{
    enzyme::{Digest, DigestGroup, Position, ProteinOccurrence},
    mass::{monoisotopic, H2O},
};
use fnv::{FnvHashMap, FnvHashSet};
use itertools::Itertools;
use smallvec::SmallVec;
use std::sync::OnceLock;

pub const INLINE_PROTEINS: usize = 1;

/// Most peptides map to one protein, so keep that accession inline.
/// Shared peptides transparently spill to heap storage.
pub type ProteinAccessions = SmallVec<[Arc<str>; INLINE_PROTEINS]>;

const MAX_COMPACT_DEFINITIONS: usize = u8::MAX as usize;

#[derive(Copy, Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
enum SiteClass {
    Nterm,
    Cterm,
    Residue,
}

impl SiteClass {
    fn from_specificity(specificity: ModificationSpecificity) -> Self {
        match specificity {
            ModificationSpecificity::PeptideN(None) | ModificationSpecificity::ProteinN(None) => {
                Self::Nterm
            }
            ModificationSpecificity::PeptideC(None) | ModificationSpecificity::ProteinC(None) => {
                Self::Cterm
            }
            _ => Self::Residue,
        }
    }

    fn from_site(site: Site) -> Self {
        match site {
            Site::Nterm => Self::Nterm,
            Site::Cterm => Self::Cterm,
            Site::Sequence(_) => Self::Residue,
        }
    }

    fn site(self, position: u8) -> Site {
        match self {
            Self::Nterm => Site::Nterm,
            Self::Cterm => Site::Cterm,
            Self::Residue => Site::Sequence(position as u32),
        }
    }
}

#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord)]
struct LookupKey {
    site: SiteClass,
    kind: ModificationKind,
    definition: ModificationDefinition,
}

#[derive(Clone, Debug)]
struct LookupRecord {
    site: SiteClass,
    kind: ModificationKind,
    definition: Arc<ModificationDefinition>,
}

/// Shared metadata addressed by the compact one-byte modification IDs stored
/// in each peptide.
#[derive(Clone, Debug, Default)]
pub struct ModificationLookup {
    records: Vec<LookupRecord>,
    ids: BTreeMap<LookupKey, u8>,
    pointer_ids: FnvHashMap<(usize, SiteClass, ModificationKind), u8>,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ModificationLookupError {
    pub definitions: usize,
}

impl std::fmt::Display for ModificationLookupError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "compact modification encoding supports at most {MAX_COMPACT_DEFINITIONS} distinct definition and site variants, but {} are required",
            self.definitions
        )
    }
}

impl std::error::Error for ModificationLookupError {}

impl ModificationLookup {
    fn from_records(
        records: impl IntoIterator<Item = LookupRecord>,
    ) -> Result<Arc<Self>, ModificationLookupError> {
        let mut unique = BTreeMap::<LookupKey, Arc<ModificationDefinition>>::new();
        let mut pointers = Vec::new();
        for record in records {
            let key = LookupKey {
                site: record.site,
                kind: record.kind,
                definition: record.definition.as_ref().clone(),
            };
            pointers.push((Arc::as_ptr(&record.definition) as usize, key.clone()));
            unique.entry(key).or_insert(record.definition);
        }
        if unique.len() > MAX_COMPACT_DEFINITIONS {
            return Err(ModificationLookupError {
                definitions: unique.len(),
            });
        }

        let mut lookup = Self::default();
        for (key, definition) in unique {
            let id = lookup.records.len() as u8;
            lookup.records.push(LookupRecord {
                site: key.site,
                kind: key.kind,
                definition,
            });
            lookup.ids.insert(key, id);
        }
        for (pointer, key) in pointers {
            lookup.pointer_ids.insert(
                (pointer, key.site, key.kind),
                *lookup.ids.get(&key).expect("compact lookup key is missing"),
            );
        }
        Ok(Arc::new(lookup))
    }

    pub fn from_definitions(
        definitions: impl IntoIterator<Item = (Site, Arc<ModificationDefinition>, ModificationKind)>,
    ) -> Result<Arc<Self>, ModificationLookupError> {
        Self::from_records(
            definitions
                .into_iter()
                .map(|(site, definition, kind)| LookupRecord {
                    site: SiteClass::from_site(site),
                    kind,
                    definition,
                }),
        )
    }

    fn id(
        &self,
        site: SiteClass,
        definition: &Arc<ModificationDefinition>,
        kind: ModificationKind,
    ) -> Option<u8> {
        self.pointer_ids
            .get(&(Arc::as_ptr(definition) as usize, site, kind))
            .copied()
            .or_else(|| self.id_value(site, definition, kind))
    }

    fn id_value(
        &self,
        site: SiteClass,
        definition: &ModificationDefinition,
        kind: ModificationKind,
    ) -> Option<u8> {
        self.ids
            .get(&LookupKey {
                site,
                kind,
                definition: definition.clone(),
            })
            .copied()
    }

    fn record(&self, id: u8) -> &LookupRecord {
        &self.records[id as usize]
    }
}

fn empty_modification_lookup() -> Arc<ModificationLookup> {
    static EMPTY: OnceLock<Arc<ModificationLookup>> = OnceLock::new();
    EMPTY
        .get_or_init(|| Arc::new(ModificationLookup::default()))
        .clone()
}

/// A residue position and shared modification-definition ID.
#[repr(C)]
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub struct EncodedModification {
    pub position: u8,
    pub modification_id: u8,
}

/// Compact peptide-local modification sites plus shared full metadata.
///
/// Four entries are stored inline. Larger collections spill to one allocation.
#[derive(Clone, Debug)]
pub struct CompactModifications {
    entries: SmallVec<[EncodedModification; 4]>,
    lookup: Arc<ModificationLookup>,
}

impl Default for CompactModifications {
    fn default() -> Self {
        Self {
            entries: SmallVec::new(),
            lookup: empty_modification_lookup(),
        }
    }
}

impl PartialEq for CompactModifications {
    fn eq(&self, other: &Self) -> bool {
        self.semantic_cmp(other) == Ordering::Equal
    }
}

impl CompactModifications {
    pub fn from_dense(masses: impl IntoIterator<Item = f32>) -> Self {
        Self::from_sparse(
            masses
                .into_iter()
                .enumerate()
                .filter(|(_, mass)| *mass != 0.0),
        )
    }

    pub fn from_sparse(masses: impl IntoIterator<Item = (usize, f32)>) -> Self {
        let masses = masses.into_iter().collect::<Vec<_>>();
        let records = masses
            .iter()
            .enumerate()
            .map(|(_, (position, mass))| {
                (
                    Site::Sequence(*position as u32),
                    Arc::new(ModificationDefinition::bare(*mass)),
                    ModificationKind::Ordinary,
                )
            })
            .collect::<Vec<_>>();
        let lookup = ModificationLookup::from_definitions(records.clone())
            .expect("dense modification array exceeds compact definition limit");
        let mut compact = Self {
            entries: SmallVec::new(),
            lookup,
        };
        for (site, definition, kind) in records {
            compact.push(site, &definition, kind);
        }
        compact.sort();
        compact
    }

    pub fn from_applied(
        modifications: impl IntoIterator<Item = AppliedModification>,
    ) -> Result<Self, ModificationLookupError> {
        let modifications = modifications.into_iter().collect::<Vec<_>>();
        let lookup = ModificationLookup::from_definitions(
            modifications
                .iter()
                .map(|applied| (applied.site, applied.modification.clone(), applied.kind)),
        )?;
        let mut compact = Self {
            entries: SmallVec::new(),
            lookup,
        };
        for applied in modifications {
            compact.push(applied.site, &applied.modification, applied.kind);
        }
        compact.sort();
        Ok(compact)
    }

    fn install_lookup(&mut self, lookup: Arc<ModificationLookup>) {
        if self.entries.is_empty() {
            self.lookup = lookup;
            return;
        }
        if Arc::ptr_eq(&self.lookup, &lookup) {
            return;
        }

        let old_lookup = self.lookup.clone();
        let records = old_lookup
            .records
            .iter()
            .chain(lookup.records.iter())
            .cloned()
            .collect::<Vec<_>>();
        let merged =
            ModificationLookup::from_records(records).unwrap_or_else(|error| panic!("{error}"));
        for encoded in &mut self.entries {
            let record = old_lookup.record(encoded.modification_id);
            encoded.modification_id = merged
                .id(record.site, &record.definition, record.kind)
                .expect("existing modification is missing from merged compact lookup");
        }
        self.lookup = merged;
    }

    fn push(
        &mut self,
        site: Site,
        definition: &Arc<ModificationDefinition>,
        kind: ModificationKind,
    ) {
        let class = SiteClass::from_site(site);
        let id = self.lookup.id(class, definition, kind).unwrap_or_else(|| {
            panic!("modification definition was not registered in the shared compact lookup")
        });
        let position = match site {
            Site::Sequence(position) => u8::try_from(position)
                .expect("peptide residue position exceeds compact modification encoding"),
            Site::Nterm | Site::Cterm => 0,
        };
        self.entries.push(EncodedModification {
            position,
            modification_id: id,
        });
    }

    fn applied(&self) -> impl ExactSizeIterator<Item = AppliedModificationRef<'_>> {
        self.entries.iter().map(|encoded| {
            let record = self.lookup.record(encoded.modification_id);
            AppliedModificationRef {
                site: record.site.site(encoded.position),
                modification: record.definition.as_ref(),
                kind: record.kind,
            }
        })
    }

    fn mass_at(&self, position: usize) -> f32 {
        self.entries
            .iter()
            .filter_map(|encoded| {
                let record = self.lookup.record(encoded.modification_id);
                (record.site == SiteClass::Residue && encoded.position as usize == position)
                    .then_some(record.definition.mass)
            })
            .sum()
    }

    /// Read all modifications at the next residue in a forward sequence scan.
    ///
    /// Compact entries are sorted by site, so callers such as ion generation
    /// can visit every entry once instead of searching the collection for each
    /// residue.
    pub(crate) fn mass_at_with_cursor(&self, position: usize, cursor: &mut usize) -> f32 {
        let mut mass = 0.0;
        while let Some(encoded) = self.entries.get(*cursor) {
            let record = self.lookup.record(encoded.modification_id);
            if record.site != SiteClass::Residue {
                *cursor += 1;
                continue;
            }

            match (encoded.position as usize).cmp(&position) {
                Ordering::Less => *cursor += 1,
                Ordering::Equal => {
                    mass += record.definition.mass;
                    *cursor += 1;
                }
                Ordering::Greater => break,
            }
        }
        mass
    }

    fn total_mass(&self) -> f32 {
        self.entries
            .iter()
            .map(|encoded| self.lookup.record(encoded.modification_id).definition.mass)
            .sum()
    }

    fn sort(&mut self) {
        let lookup = self.lookup.clone();
        self.entries.sort_unstable_by(|left, right| {
            let left_record = lookup.record(left.modification_id);
            let right_record = lookup.record(right.modification_id);
            left_record
                .site
                .site(left.position)
                .cmp(&right_record.site.site(right.position))
                .then_with(|| left_record.definition.cmp(&right_record.definition))
                .then_with(|| left_record.kind.cmp(&right_record.kind))
        });
    }

    fn semantic_cmp(&self, other: &Self) -> Ordering {
        self.applied()
            .map(|applied| (applied.site, applied.modification, applied.kind))
            .cmp(
                other
                    .applied()
                    .map(|applied| (applied.site, applied.modification, applied.kind)),
            )
    }

    fn reverse_internal(&mut self, last: usize) {
        let lookup = self.lookup.clone();
        for encoded in &mut self.entries {
            if lookup.record(encoded.modification_id).site != SiteClass::Residue {
                continue;
            }
            let position = encoded.position as usize;
            if (1..last).contains(&position) {
                encoded.position = (last - position) as u8;
            }
        }
        self.sort();
    }

    pub fn len(&self) -> usize {
        self.entries.len()
    }

    pub fn is_empty(&self) -> bool {
        self.entries.is_empty()
    }

    pub fn spilled(&self) -> bool {
        self.entries.spilled()
    }

    pub(crate) fn heap_bytes(&self) -> usize {
        usize::from(self.entries.spilled())
            * self.entries.capacity()
            * std::mem::size_of::<EncodedModification>()
    }

    fn relocate_mass(&mut self, mass: f32, candidates: &[usize], chosen: &[usize], epsilon: f32) {
        let selected_id = self.entries.iter().find_map(|encoded| {
            let record = self.lookup.record(encoded.modification_id);
            (record.site == SiteClass::Residue
                && candidates.contains(&(encoded.position as usize))
                && (record.definition.mass - mass).abs() < epsilon)
                .then_some(encoded.modification_id)
        });
        let Some(modification_id) = selected_id else {
            return;
        };
        let lookup = self.lookup.clone();
        self.entries.retain(|encoded| {
            let record = lookup.record(encoded.modification_id);
            !(record.site == SiteClass::Residue
                && candidates.contains(&(encoded.position as usize))
                && (record.definition.mass - mass).abs() < epsilon)
        });
        self.entries.extend(
            chosen.iter().map(|position| EncodedModification {
                position: u8::try_from(*position)
                    .expect("PTM position exceeds compact modification encoding"),
                modification_id,
            }),
        );
        self.sort();
    }
}

#[derive(Clone, PartialEq, Default)]
pub struct Peptide {
    pub decoy: bool,
    pub sequence: Arc<[u8]>,
    /// Compact modification sites with shared definition metadata.
    pub modifications: CompactModifications,
    /// Precursor-resolved labeling channel, when this peptide was generated
    /// from a coherent database label definition.
    pub label_channel: Option<Arc<str>>,
    /// Library-provided base peptidoform identity when applied modification
    /// provenance is unavailable after parsing ProForma mass deltas.
    pub label_group_override: Option<Arc<str>>,
    /// Modification on peptide C-terminus
    pub nterm: Option<f32>,
    /// Modification on peptide C-terminus
    pub cterm: Option<f32>,
    /// Monoisotopic mass, inclusive of N/C-terminal mods
    pub monoisotopic: f32,
    /// Number of missed cleavages for this sequence
    pub missed_cleavages: u8,
    /// Is this a semi-enzymatic peptide?
    pub semi_enzymatic: bool,
    /// Where is this peptide located in the protein?
    pub position: Position,

    pub proteins: ProteinAccessions,
    /// Protein-coordinate occurrences shared by all modification variants.
    pub protein_sites: Arc<[ProteinOccurrence]>,
}

pub trait ModificationSource {
    fn definition(&self) -> Arc<ModificationDefinition>;
}

#[derive(Clone)]
pub(crate) struct VariableRule {
    pub specificity: ModificationSpecificity,
    pub modification: Arc<ModificationDefinition>,
    pub max_count: Option<usize>,
    pub site_mode: SiteMode,
    /// Rules with the same named modification share one occurrence counter.
    pub count_group: usize,
}

/// Channel-resolved modification definitions shared by every peptide produced
/// during one database expansion.
pub(crate) struct LabelModificationCache {
    by_source: FnvHashMap<usize, BTreeMap<Arc<str>, Arc<ModificationDefinition>>>,
}

impl LabelModificationCache {
    pub(crate) fn new<'a>(
        definitions: impl IntoIterator<Item = &'a Arc<ModificationDefinition>>,
        channels: &[Arc<str>],
    ) -> Self {
        let mut by_source: FnvHashMap<usize, BTreeMap<Arc<str>, Arc<ModificationDefinition>>> =
            FnvHashMap::default();
        let mut interned: BTreeMap<ModificationDefinition, Arc<ModificationDefinition>> =
            BTreeMap::new();

        for definition in definitions {
            if definition.channel_offsets.is_empty() {
                continue;
            }
            let source = Arc::as_ptr(definition) as usize;
            let resolved_by_channel = by_source.entry(source).or_default();
            for channel in channels {
                let offset = definition.channel_offsets[channel];
                let resolved = definition.with_mass(definition.mass + offset);
                let shared = interned
                    .entry(resolved.clone())
                    .or_insert_with(|| Arc::new(resolved))
                    .clone();
                resolved_by_channel.insert(channel.clone(), shared);
            }
        }

        Self { by_source }
    }

    #[inline]
    fn resolve(
        &self,
        definition: &Arc<ModificationDefinition>,
        channel: &str,
    ) -> Option<Arc<ModificationDefinition>> {
        self.by_source
            .get(&(Arc::as_ptr(definition) as usize))?
            .get(channel)
            .cloned()
    }
}

impl ModificationLookup {
    pub(crate) fn for_rules(
        variable_mods: &[VariableRule],
        static_mods: &HashMap<ModificationSpecificity, Arc<ModificationDefinition>>,
        channels: &[Arc<str>],
        labels: &LabelModificationCache,
    ) -> Result<Arc<Self>, ModificationLookupError> {
        let mut records = Vec::new();
        for (specificity, definition) in variable_mods
            .iter()
            .map(|rule| (rule.specificity, &rule.modification))
            .chain(
                static_mods
                    .iter()
                    .map(|(specificity, definition)| (*specificity, definition)),
            )
        {
            let site = SiteClass::from_specificity(specificity);
            let kind = if definition.channel_offsets.is_empty() {
                ModificationKind::Ordinary
            } else {
                ModificationKind::ChannelBase
            };
            records.push(LookupRecord {
                site,
                kind,
                definition: definition.clone(),
            });
            if kind == ModificationKind::ChannelBase {
                for channel in channels {
                    if let Some(definition) = labels.resolve(definition, channel) {
                        records.push(LookupRecord {
                            site,
                            kind: ModificationKind::Label,
                            definition,
                        });
                    }
                }
            }
        }
        Self::from_records(records)
    }
}

#[derive(Clone, Debug)]
pub(crate) struct LibrarySite {
    pub position: u32,
    pub modification: Arc<str>,
}

#[derive(Copy, Clone)]
struct ModificationCandidate {
    site: Site,
    rule: usize,
    library_supported: bool,
}

impl ModificationSource for f32 {
    fn definition(&self) -> Arc<ModificationDefinition> {
        Arc::new(ModificationDefinition::bare(*self))
    }
}

impl ModificationSource for Arc<ModificationDefinition> {
    fn definition(&self) -> Arc<ModificationDefinition> {
        self.clone()
    }
}

impl Peptide {
    #[inline]
    pub fn modification_at(&self, index: usize) -> f32 {
        self.modifications.mass_at(index)
    }

    pub fn applied_modifications(
        &self,
    ) -> impl ExactSizeIterator<Item = AppliedModificationRef<'_>> {
        self.modifications.applied()
    }

    pub(crate) fn relocate_modification_mass(
        &mut self,
        mass: f32,
        candidates: &[usize],
        chosen: &[usize],
        epsilon: f32,
    ) {
        self.modifications
            .relocate_mass(mass, candidates, chosen, epsilon);
    }

    pub fn initial_sort(&self, other: &Self) -> std::cmp::Ordering {
        self.sequence
            .cmp(&other.sequence)
            .then_with(|| {
                (0..self.sequence.len())
                    .find_map(|index| {
                        let left = self.modification_at(index);
                        let right = other.modification_at(index);
                        (left != right).then(|| left.partial_cmp(&right).unwrap_or(Ordering::Equal))
                    })
                    .unwrap_or(Ordering::Equal)
            })
            .then_with(|| {
                self.nterm
                    .partial_cmp(&other.nterm)
                    .unwrap_or(Ordering::Equal)
            })
            .then_with(|| {
                self.cterm
                    .partial_cmp(&other.cterm)
                    .unwrap_or(Ordering::Equal)
            })
            .then_with(|| self.modifications.semantic_cmp(&other.modifications))
            .then_with(|| self.label_channel.cmp(&other.label_channel))
    }
}

impl Debug for Peptide {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("Peptide")
            .field("proteins", &self.proteins)
            .field("protein_sites", &self.protein_sites)
            .field("decoy", &self.decoy)
            .field(
                "sequence",
                &std::str::from_utf8(&self.sequence).unwrap_or("error"),
            )
            .field("nterm", &self.nterm)
            .field("cterm", &self.cterm)
            .field(
                "applied_modifications",
                &self.applied_modifications().collect::<Vec<_>>(),
            )
            .field("monoisotopic", &self.monoisotopic)
            .field("missed_cleavages", &self.missed_cleavages)
            .field("position", &self.position)
            .finish()
    }
}

impl Peptide {
    pub fn label(&self) -> i32 {
        match self.decoy {
            true => -1,
            false => 1,
        }
    }

    pub fn proteins(&self, decoy_tag: &str, generate_decoys: bool) -> String {
        if self.decoy {
            self.proteins
                .iter()
                .map(|s| {
                    if generate_decoys {
                        format!("{}{}", decoy_tag, s)
                    } else {
                        s.to_string()
                    }
                })
                .join(";")
        } else {
            self.proteins.iter().join(";")
        }
    }

    pub fn modification_count(&self, target: ModificationSpecificity, mass: f32) -> usize {
        if !self.modifications.is_empty() {
            return self
                .applied_modifications()
                .filter(|applied| applied.modification.mass == mass)
                .filter(|applied| match (target, applied.site, self.position) {
                    (ModificationSpecificity::PeptideN(None), Site::Nterm, _) => true,
                    (ModificationSpecificity::PeptideC(None), Site::Cterm, _) => true,
                    (
                        ModificationSpecificity::ProteinN(None),
                        Site::Nterm,
                        Position::Nterm | Position::Full,
                    ) => true,
                    (
                        ModificationSpecificity::ProteinC(None),
                        Site::Cterm,
                        Position::Cterm | Position::Full,
                    ) => true,
                    (ModificationSpecificity::PeptideN(Some(residue)), Site::Sequence(0), _)
                    | (
                        ModificationSpecificity::ProteinN(Some(residue)),
                        Site::Sequence(0),
                        Position::Nterm | Position::Full,
                    ) => self.sequence.first() == Some(&residue),
                    (
                        ModificationSpecificity::PeptideC(Some(residue)),
                        Site::Sequence(index),
                        _,
                    )
                    | (
                        ModificationSpecificity::ProteinC(Some(residue)),
                        Site::Sequence(index),
                        Position::Cterm | Position::Full,
                    ) => {
                        index as usize == self.sequence.len().saturating_sub(1)
                            && self.sequence.last() == Some(&residue)
                    }
                    (ModificationSpecificity::Residue(residue), Site::Sequence(index), _) => {
                        self.sequence.get(index as usize) == Some(&residue)
                    }
                    _ => false,
                })
                .count();
        }
        match target {
            ModificationSpecificity::PeptideN(r) | ModificationSpecificity::ProteinN(r) => {
                if r.map(|resi| resi == *self.sequence.first().unwrap_or(&0))
                    .unwrap_or(true)
                    && self.nterm.unwrap_or_default() == mass
                {
                    1
                } else {
                    0
                }
            }
            ModificationSpecificity::PeptideC(r) | ModificationSpecificity::ProteinC(r) => {
                if r.map(|resi| resi == *self.sequence.last().unwrap_or(&0))
                    .unwrap_or(true)
                    && self.cterm.unwrap_or_default() == mass
                {
                    1
                } else {
                    0
                }
            }
            ModificationSpecificity::Residue(resi) => self
                .sequence
                .iter()
                .enumerate()
                .filter(|(index, residue)| {
                    resi == **residue && mass == self.modification_at(*index)
                })
                .count(),
        }
    }

    fn modification_mass(&self) -> f32 {
        self.modifications.total_mass()
    }

    fn finalize_modifications(&mut self) {
        self.modifications.sort();
        self.monoisotopic += self.modification_mass();
    }

    /// Apply all variable mods in `sites` to self
    fn apply_site(&mut self, site: Site, modification: Arc<ModificationDefinition>) {
        let kind = if modification.channel_offsets.is_empty() {
            ModificationKind::Ordinary
        } else {
            ModificationKind::ChannelBase
        };
        self.apply_site_with_kind(site, modification, kind, false);
    }

    fn apply_site_with_kind(
        &mut self,
        site: Site,
        modification: Arc<ModificationDefinition>,
        kind: ModificationKind,
        stack: bool,
    ) {
        let mass = modification.mass;
        let mut applied = false;
        match site {
            Site::Nterm => {
                if stack || self.nterm.is_none() {
                    self.nterm = Some(self.nterm.unwrap_or_default() + mass);
                    applied = true;
                }
            }
            Site::Cterm => {
                if stack || self.cterm.is_none() {
                    self.cterm = Some(self.cterm.unwrap_or_default() + mass);
                    applied = true;
                }
            }
            Site::Sequence(index) => {
                if stack || self.modification_at(index as usize) == 0.0 {
                    applied = true;
                }
            }
        }
        if applied {
            self.modifications.push(site, &modification, kind);
        }
    }

    /// Resolve every channel-aware modification on this peptidoform to one
    /// coherent precursor channel.
    pub(crate) fn apply_label_channel(
        mut self,
        channel: Arc<str>,
        cache: &LabelModificationCache,
    ) -> Self {
        let before = self.modification_mass();
        let lookup = self.modifications.lookup.clone();
        for encoded in &mut self.modifications.entries {
            let record = lookup.record(encoded.modification_id);
            if record.kind != ModificationKind::ChannelBase {
                continue;
            }
            let offset = record.definition.channel_offsets[&channel];
            let site = record.site.site(encoded.position);
            match site {
                Site::Nterm => self.nterm = Some(self.nterm.unwrap_or_default() + offset),
                Site::Cterm => self.cterm = Some(self.cterm.unwrap_or_default() + offset),
                Site::Sequence(_) => {}
            }
            let resolved = cache
                .resolve(&record.definition, &channel)
                .unwrap_or_else(|| {
                    Arc::new(record.definition.with_mass(record.definition.mass + offset))
                });
            encoded.modification_id = lookup
                .id(record.site, &resolved, ModificationKind::Label)
                .expect("resolved label definition missing from compact lookup");
        }
        self.monoisotopic += self.modification_mass() - before;
        self.label_channel = Some(channel);
        self.modifications.sort();
        self
    }

    /// Exact peptidoform identity with precursor-label modifications removed.
    /// This is used for channel grouping, peptide-level FDR, and protein
    /// inference without collapsing ordinary biological modifications.
    pub fn label_group(&self) -> String {
        if let Some(group) = &self.label_group_override {
            return group.to_string();
        }
        if self.label_channel.is_none() {
            return self.to_string();
        }
        let mut base = self.clone();
        let channel = base.label_channel.clone().unwrap();
        let lookup = base.modifications.lookup.clone();
        for encoded in &mut base.modifications.entries {
            let record = lookup.record(encoded.modification_id);
            if record.kind != ModificationKind::Label {
                continue;
            }
            let offset = record.definition.channel_offsets[&channel];
            match record.site.site(encoded.position) {
                Site::Nterm => base.nterm = nonzero_mass(base.nterm.unwrap_or_default() - offset),
                Site::Cterm => base.cterm = nonzero_mass(base.cterm.unwrap_or_default() - offset),
                Site::Sequence(_) => {}
            }
            base.monoisotopic -= offset;
            let unresolved = record.definition.with_mass(record.definition.mass - offset);
            encoded.modification_id = lookup
                .id_value(record.site, &unresolved, ModificationKind::ChannelBase)
                .expect("base label definition missing from compact lookup");
        }
        base.label_channel = None;
        base.label_group_override = None;
        base.modifications.sort();
        base.to_string()
    }

    fn push_resi(
        &self,
        acc: &mut Vec<(Site, f32, usize)>,
        target: ModificationSpecificity,
        mass: f32,
        mod_idx: usize,
    ) {
        match (target, self.position) {
            (ModificationSpecificity::PeptideN(None), _) => acc.push((Site::Nterm, mass, mod_idx)),
            (ModificationSpecificity::PeptideN(Some(resi)), _)
                if resi == *self.sequence.first().unwrap_or(&0) =>
            {
                acc.push((Site::Sequence(0), mass, mod_idx))
            }
            (ModificationSpecificity::PeptideC(None), _) => acc.push((Site::Cterm, mass, mod_idx)),
            (ModificationSpecificity::PeptideC(Some(resi)), _)
                if resi == *self.sequence.last().unwrap_or(&0) =>
            {
                acc.push((
                    Site::Sequence(self.sequence.len().saturating_sub(1) as u32),
                    mass,
                    mod_idx,
                ))
            }
            (ModificationSpecificity::ProteinN(None), Position::Nterm | Position::Full) => {
                acc.push((Site::Nterm, mass, mod_idx))
            }
            (ModificationSpecificity::ProteinN(Some(resi)), Position::Nterm | Position::Full)
                if resi == *self.sequence.first().unwrap_or(&0) =>
            {
                acc.push((Site::Sequence(0), mass, mod_idx))
            }
            (ModificationSpecificity::ProteinC(None), Position::Cterm | Position::Full) => {
                acc.push((Site::Cterm, mass, mod_idx))
            }
            (ModificationSpecificity::ProteinC(Some(resi)), Position::Cterm | Position::Full)
                if resi == *self.sequence.last().unwrap_or(&0) =>
            {
                acc.push((
                    Site::Sequence(self.sequence.len().saturating_sub(1) as u32),
                    mass,
                    mod_idx,
                ))
            }
            (ModificationSpecificity::Residue(resi), _) => {
                acc.extend(
                    self.sequence
                        .iter()
                        .enumerate()
                        .filter_map(|(idx, residue)| {
                            if resi == *residue {
                                Some((Site::Sequence(idx as u32), mass, mod_idx))
                            } else {
                                None
                            }
                        }),
                );
            }
            _ => {}
        }
    }

    fn static_mods_all<M: ModificationSource>(
        &mut self,
        static_mods: &HashMap<ModificationSpecificity, M>,
    ) {
        let mut occupied = [false; u8::MAX as usize];
        for applied in self.applied_modifications() {
            if let Site::Sequence(position) = applied.site {
                occupied[position as usize] = true;
            }
        }
        let mut nterm_occupied = self.nterm.is_some();
        let mut cterm_occupied = self.cterm.is_some();

        for (target, source) in static_mods {
            let modification = source.definition();
            let kind = if modification.channel_offsets.is_empty() {
                ModificationKind::Ordinary
            } else {
                ModificationKind::ChannelBase
            };
            let mut sites = SmallVec::<[Site; 4]>::new();
            match (*target, self.position) {
                (ModificationSpecificity::PeptideN(None), _)
                | (ModificationSpecificity::ProteinN(None), Position::Nterm | Position::Full)
                    if !nterm_occupied =>
                {
                    nterm_occupied = true;
                    sites.push(Site::Nterm);
                }
                (ModificationSpecificity::PeptideC(None), _)
                | (ModificationSpecificity::ProteinC(None), Position::Cterm | Position::Full)
                    if !cterm_occupied =>
                {
                    cterm_occupied = true;
                    sites.push(Site::Cterm);
                }
                (ModificationSpecificity::PeptideN(Some(residue)), _)
                | (
                    ModificationSpecificity::ProteinN(Some(residue)),
                    Position::Nterm | Position::Full,
                ) if self.sequence.first() == Some(&residue) && !occupied[0] => {
                    occupied[0] = true;
                    sites.push(Site::Sequence(0));
                }
                (ModificationSpecificity::PeptideC(Some(residue)), _)
                | (
                    ModificationSpecificity::ProteinC(Some(residue)),
                    Position::Cterm | Position::Full,
                ) if self.sequence.last() == Some(&residue) => {
                    let position = self.sequence.len().saturating_sub(1);
                    if !occupied[position] {
                        occupied[position] = true;
                        sites.push(Site::Sequence(position as u32));
                    }
                }
                (ModificationSpecificity::Residue(residue), _) => {
                    for (position, observed) in self.sequence.iter().copied().enumerate() {
                        if observed == residue && !occupied[position] {
                            occupied[position] = true;
                            sites.push(Site::Sequence(position as u32));
                        }
                    }
                }
                _ => {}
            }
            for site in sites {
                self.apply_site_with_kind(site, modification.clone(), kind, true);
            }
        }
    }

    /// Apply variable modifications, then static modifications to a peptide.
    /// `variable_mods` entries are `(specificity, definition, per_mod_limit)`.
    /// `max_combinations` caps total variants (unmodified + modified); fewer-PTM
    /// variants are always generated first so they are preferred when truncating.
    pub fn apply<M: ModificationSource>(
        self,
        variable_mods: &[(ModificationSpecificity, M, Option<usize>)],
        static_mods: &HashMap<ModificationSpecificity, M>,
        combinations: usize,
        max_combinations: Option<usize>,
    ) -> Vec<Peptide> {
        let rules = variable_mods
            .iter()
            .enumerate()
            .map(
                |(count_group, (specificity, modification, max_count))| VariableRule {
                    specificity: *specificity,
                    modification: modification.definition(),
                    max_count: *max_count,
                    site_mode: SiteMode::Exhaustive,
                    count_group,
                },
            )
            .collect::<Vec<_>>();
        let static_mods = static_mods
            .iter()
            .map(|(specificity, modification)| (*specificity, modification.definition()))
            .collect::<HashMap<_, _>>();
        let labels = LabelModificationCache::new(
            rules
                .iter()
                .map(|rule| &rule.modification)
                .chain(static_mods.values()),
            &[],
        );
        let lookup = ModificationLookup::for_rules(&rules, &static_mods, &[], &labels)
            .unwrap_or_else(|error| panic!("{error}"));
        self.apply_rules(
            &rules,
            &[],
            &static_mods,
            lookup,
            combinations,
            combinations,
            max_combinations,
        )
    }

    /// Apply config-defined and library-supported variable modifications in a
    /// single enumeration. Library-supported placements do not consume the
    /// exhaustive budget, but all placements consume the total budget and the
    /// named modification's shared `max_count`.
    pub(crate) fn apply_rules(
        mut self,
        variable_mods: &[VariableRule],
        library_sites: &[LibrarySite],
        static_mods: &HashMap<ModificationSpecificity, Arc<ModificationDefinition>>,
        lookup: Arc<ModificationLookup>,
        max_exhaustive_mods: usize,
        max_total_mods: usize,
        max_combinations: Option<usize>,
    ) -> Vec<Peptide> {
        self.modifications.install_lookup(lookup);
        if variable_mods.is_empty() {
            self.static_mods_all(static_mods);
            self.finalize_modifications();
            vec![self]
        } else {
            let mut candidates: Vec<ModificationCandidate> = Vec::new();
            for (rule_idx, rule) in variable_mods.iter().enumerate() {
                let mut compatible = Vec::new();
                self.push_resi(
                    &mut compatible,
                    rule.specificity,
                    rule.modification.mass,
                    rule_idx,
                );
                for (site, _, _) in compatible {
                    let library_supported = match (site, rule.modification.name.as_deref()) {
                        (Site::Sequence(position), Some(name))
                            if rule.site_mode != SiteMode::Exhaustive =>
                        {
                            library_sites.iter().any(|library| {
                                library.position == position
                                    && library.modification.as_ref() == name
                            })
                        }
                        (Site::Nterm, Some(name)) if rule.site_mode != SiteMode::Exhaustive => {
                            library_sites.iter().any(|library| {
                                library.position == 0 && library.modification.as_ref() == name
                            })
                        }
                        (Site::Cterm, Some(name)) if rule.site_mode != SiteMode::Exhaustive => {
                            let position = self.sequence.len().saturating_sub(1) as u32;
                            library_sites.iter().any(|library| {
                                library.position == position
                                    && library.modification.as_ref() == name
                            })
                        }
                        _ => false,
                    };
                    if rule.site_mode != SiteMode::Library || library_supported {
                        candidates.push(ModificationCandidate {
                            site,
                            rule: rule_idx,
                            library_supported,
                        });
                    }
                }
            }

            // Multiple specificity entries can describe the same named
            // modification at one site. Collapse them and preserve library
            // support if either route supplies it.
            let mut mods: Vec<ModificationCandidate> = Vec::with_capacity(candidates.len());
            let mut candidate_indices: HashMap<(Site, usize), usize> = HashMap::new();
            for candidate in candidates {
                let key = (candidate.site, variable_mods[candidate.rule].count_group);
                if let Some(index) = candidate_indices.get(&key).copied() {
                    mods[index].library_supported |= candidate.library_supported;
                } else {
                    candidate_indices.insert(key, mods.len());
                    mods.push(candidate);
                }
            }

            let mut modified = Vec::new();
            modified.push(self.clone());
            let group_count = variable_mods
                .iter()
                .map(|rule| rule.count_group)
                .max()
                .map_or(0, |group| group + 1);
            let mut mod_counts = vec![0usize; group_count];
            let mut occupied = FnvHashSet::default();
            for n in 1..=max_total_mods.min(mods.len()) {
                if !enumerate_modifications(
                    &self,
                    &mods,
                    variable_mods,
                    0,
                    n,
                    0,
                    max_exhaustive_mods,
                    &mut occupied,
                    &mut Vec::with_capacity(n),
                    &mut mod_counts,
                    &mut modified,
                    max_combinations,
                ) {
                    break;
                }
            }

            // Apply static mods to all peptides
            for peptide in modified.iter_mut() {
                peptide.static_mods_all(static_mods);
                peptide.finalize_modifications();
            }

            modified
        }
    }

    pub fn reverse(&self) -> Peptide {
        let mut pep = self.clone();
        pep.decoy = !self.decoy;
        let n = pep.sequence.len().saturating_sub(1);
        if n > 1 {
            let mut s = Vec::from(pep.sequence.as_ref());
            s[1..n].reverse();
            pep.sequence = Arc::from(s.into_boxed_slice());
            pep.modifications.reverse_internal(n);
        }
        pep
    }

    pub(crate) fn modification_tag(&self, site: Site, mass: f32) -> String {
        let applied = self
            .applied_modifications()
            .filter(|applied| applied.site == site)
            .collect::<Vec<_>>();
        let represented_mass = applied
            .iter()
            .map(|applied| applied.modification.mass)
            .sum::<f32>();
        if !applied.is_empty() && (represented_mass - mass).abs() < 1e-4 {
            let mut tag = String::new();
            for applied in applied {
                if let Some(name) = applied.modification.name.as_deref() {
                    tag.push_str(&format!("[{name}]"));
                } else {
                    tag.push_str(&format!("[{:+}]", applied.modification.mass));
                }
            }
            tag
        } else {
            format!("[{mass:+}]")
        }
    }

    fn fmt_mod(&self, f: &mut std::fmt::Formatter<'_>, site: Site, mass: f32) -> std::fmt::Result {
        f.write_str(&self.modification_tag(site, mass))
    }
}

fn nonzero_mass(mass: f32) -> Option<f32> {
    (mass.abs() >= 1e-5).then_some(mass)
}

fn enumerate_modifications(
    peptide: &Peptide,
    modifications: &[ModificationCandidate],
    variable_mods: &[VariableRule],
    start: usize,
    remaining: usize,
    exhaustive_count: usize,
    max_exhaustive_mods: usize,
    occupied: &mut FnvHashSet<Site>,
    selected: &mut Vec<ModificationCandidate>,
    mod_counts: &mut [usize],
    output: &mut Vec<Peptide>,
    max_combinations: Option<usize>,
) -> bool {
    for idx in start..modifications.len() {
        let modification = modifications[idx];
        if occupied.insert(modification.site) {
            let rule = &variable_mods[modification.rule];
            let group = rule.count_group;
            let next_exhaustive = exhaustive_count + usize::from(!modification.library_supported);
            if next_exhaustive > max_exhaustive_mods
                || rule
                    .max_count
                    .is_some_and(|limit| mod_counts[group] >= limit)
            {
                occupied.remove(&modification.site);
                continue;
            }

            mod_counts[group] += 1;
            selected.push(modification);
            if remaining == 1 {
                if max_combinations.is_some_and(|cap| output.len() >= cap) {
                    selected.pop();
                    mod_counts[group] -= 1;
                    occupied.remove(&modification.site);
                    return false;
                }
                let mut modified = peptide.clone();
                for selected in selected.iter().copied() {
                    modified.apply_site(
                        selected.site,
                        variable_mods[selected.rule].modification.clone(),
                    );
                }
                output.push(modified);
            } else if !enumerate_modifications(
                peptide,
                modifications,
                variable_mods,
                idx + 1,
                remaining - 1,
                next_exhaustive,
                max_exhaustive_mods,
                occupied,
                selected,
                mod_counts,
                output,
                max_combinations,
            ) {
                selected.pop();
                mod_counts[group] -= 1;
                occupied.remove(&modification.site);
                return false;
            }
            selected.pop();
            mod_counts[group] -= 1;
            occupied.remove(&modification.site);
        }
    }

    true
}

#[derive(Copy, Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum Site {
    Nterm,
    Cterm,
    Sequence(u32),
}

#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub struct AppliedModification {
    pub site: Site,
    pub modification: Arc<ModificationDefinition>,
    pub kind: ModificationKind,
}

#[derive(Copy, Clone, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub struct AppliedModificationRef<'a> {
    pub site: Site,
    pub modification: &'a ModificationDefinition,
    pub kind: ModificationKind,
}

#[derive(Copy, Clone, Debug, Default, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum ModificationKind {
    #[default]
    Ordinary,
    ChannelBase,
    Label,
}

#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub enum PeptideError {
    InvalidSequence(String),
    SequenceTooLong { length: usize, maximum: usize },
}

impl std::fmt::Display for PeptideError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::InvalidSequence(sequence) => write!(f, "invalid peptide sequence: {sequence}"),
            Self::SequenceTooLong { length, maximum } => write!(
                f,
                "peptide has {length} residues, but compact modification encoding supports at most {maximum}"
            ),
        }
    }
}

impl std::error::Error for PeptideError {}

impl TryFrom<DigestGroup> for Peptide {
    type Error = PeptideError;

    fn try_from(value: DigestGroup) -> Result<Self, Self::Error> {
        let mut pep = Peptide::try_from(value.reference)?;
        pep.protein_sites = Arc::from(value.origins.clone());
        pep.proteins = value
            .origins
            .into_iter()
            .map(|origin| origin.protein)
            .collect();
        pep.proteins.sort_unstable();
        pep.proteins.dedup();
        Ok(pep)
    }
}

impl TryFrom<Digest> for Peptide {
    type Error = PeptideError;

    fn try_from(value: Digest) -> Result<Self, Self::Error> {
        if value.sequence.len() > u8::MAX as usize {
            return Err(PeptideError::SequenceTooLong {
                length: value.sequence.len(),
                maximum: u8::MAX as usize,
            });
        }
        let mut mass = H2O;
        let protein_sites: Arc<[ProteinOccurrence]> = value
            .protein_start
            .map(|start| {
                vec![ProteinOccurrence {
                    protein: value.protein.clone(),
                    start: Some(start),
                    prev_aa: value.prev_aa,
                    next_aa: value.next_aa,
                }]
                .into()
            })
            .unwrap_or_default();
        // This is an important invariant to enforce, that ensures safety
        // while reversing peptide sequences
        if !value.sequence.is_ascii() {
            return Err(PeptideError::InvalidSequence(value.sequence));
        }
        for c in value.sequence.as_bytes() {
            let mono = monoisotopic(*c);
            if mono == 0.0 {
                return Err(PeptideError::InvalidSequence(value.sequence));
            }
            mass += mono;
        }

        Ok(Peptide {
            decoy: value.decoy,
            position: value.position,
            modifications: CompactModifications::default(),
            label_channel: None,
            label_group_override: None,
            sequence: Arc::from(value.sequence.into_bytes().into_boxed_slice()),
            monoisotopic: mass,
            nterm: None,
            cterm: None,
            missed_cleavages: value.missed_cleavages,
            semi_enzymatic: value.semi_enzymatic,
            proteins: smallvec::smallvec![value.protein],
            protein_sites,
        })
    }
}

impl std::fmt::Display for Peptide {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        if let Some(m) = self.nterm {
            self.fmt_mod(f, Site::Nterm, m)?;
            write!(f, "-")?;
        }
        for (index, c) in self.sequence.iter().enumerate() {
            let modification = self.modification_at(index);
            if modification != 0.0 {
                write!(f, "{}", *c as char)?;
                self.fmt_mod(f, Site::Sequence(index as u32), modification)?;
            } else {
                write!(f, "{}", *c as char)?;
            }
        }
        if let Some(m) = self.cterm {
            write!(f, "-")?;
            self.fmt_mod(f, Site::Cterm, m)?;
        }
        Ok(())
    }
}

#[cfg(test)]
#[path = "../tests/unit/peptide.rs"]
mod test;
