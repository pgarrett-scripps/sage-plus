//! The modification lookup table and the compact per-peptide encoding of it.
//!
//! Moved out of `peptide.rs` unchanged; every item keeps its `peptide::` path.

use std::cmp::Ordering;
use std::{
    collections::{BTreeMap, HashMap},
    fmt::Debug,
    sync::Arc,
};

use super::*;
use crate::modification::{ModificationDefinition, ModificationSpecificity};
use fnv::FnvHashMap;
use smallvec::SmallVec;
use std::sync::OnceLock;

const MAX_COMPACT_DEFINITIONS: usize = u8::MAX as usize;

#[derive(Copy, Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub(super) enum SiteClass {
    Nterm,
    Cterm,
    Residue,
}

impl SiteClass {
    fn from_specificity(specificity: ModificationSpecificity) -> Self {
        match specificity {
            ModificationSpecificity::PeptideN(None)
            | ModificationSpecificity::ProteinN(None)
            | ModificationSpecificity::PeptideNTerm(_)
            | ModificationSpecificity::ProteinNTerm(_) => Self::Nterm,
            ModificationSpecificity::PeptideC(None)
            | ModificationSpecificity::ProteinC(None)
            | ModificationSpecificity::PeptideCTerm(_)
            | ModificationSpecificity::ProteinCTerm(_) => Self::Cterm,
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

    pub(super) fn site(self, position: u8) -> Site {
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
pub(super) struct LookupRecord {
    pub(super) site: SiteClass,
    pub(super) kind: ModificationKind,
    pub(super) definition: Arc<ModificationDefinition>,
    /// Channel-base definition that a `Label` record was resolved from. Kept by
    /// identity because `(mass + offset) - offset` is not exact in f32.
    pub(super) base: Option<Arc<ModificationDefinition>>,
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
        let mut unique = BTreeMap::<LookupKey, LookupRecord>::new();
        let mut pointers = Vec::new();
        for record in records {
            let key = LookupKey {
                site: record.site,
                kind: record.kind,
                definition: record.definition.as_ref().clone(),
            };
            pointers.push((Arc::as_ptr(&record.definition) as usize, key.clone()));
            unique.entry(key).or_insert(record);
        }
        if unique.len() > MAX_COMPACT_DEFINITIONS {
            return Err(ModificationLookupError {
                definitions: unique.len(),
            });
        }

        let mut lookup = Self::default();
        for (key, record) in unique {
            let id = lookup.records.len() as u8;
            lookup.records.push(record);
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
                    base: None,
                }),
        )
    }

    pub(super) fn id(
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

    pub(super) fn record(&self, id: u8) -> &LookupRecord {
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
    pub(super) entries: SmallVec<[EncodedModification; 4]>,
    pub(super) lookup: Arc<ModificationLookup>,
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
            .map(|(position, mass)| {
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

    pub(super) fn install_lookup(&mut self, lookup: Arc<ModificationLookup>) {
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

    pub(super) fn push(
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

    pub(super) fn applied(&self) -> impl ExactSizeIterator<Item = AppliedModificationRef<'_>> {
        self.entries.iter().map(|encoded| {
            let record = self.lookup.record(encoded.modification_id);
            AppliedModificationRef {
                site: record.site.site(encoded.position),
                modification: record.definition.as_ref(),
                kind: record.kind,
            }
        })
    }

    pub(super) fn mass_at(&self, position: usize) -> f32 {
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

    pub(super) fn total_mass(&self) -> f32 {
        self.entries
            .iter()
            .map(|encoded| self.lookup.record(encoded.modification_id).definition.mass)
            .sum()
    }

    pub(super) fn sort(&mut self) {
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

    pub(super) fn semantic_cmp(&self, other: &Self) -> Ordering {
        self.applied()
            .map(|applied| (applied.site, applied.modification, applied.kind))
            .cmp(
                other
                    .applied()
                    .map(|applied| (applied.site, applied.modification, applied.kind)),
            )
    }

    pub(super) fn reverse_internal(&mut self, last: usize) {
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
    pub(super) fn resolve(
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
                base: None,
            });
            if kind == ModificationKind::ChannelBase {
                for channel in channels {
                    if let Some(resolved) = labels.resolve(definition, channel) {
                        records.push(LookupRecord {
                            site,
                            kind: ModificationKind::Label,
                            definition: resolved,
                            base: Some(definition.clone()),
                        });
                    }
                }
            }
        }
        Self::from_records(records)
    }
}
