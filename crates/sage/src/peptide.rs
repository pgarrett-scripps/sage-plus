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

pub const INLINE_PROTEINS: usize = 1;

/// Most peptides map to one protein, so keep that accession inline.
/// Shared peptides transparently spill to heap storage.
pub type ProteinAccessions = SmallVec<[Arc<str>; INLINE_PROTEINS]>;

#[derive(Clone, PartialEq, Default)]
pub struct Peptide {
    pub decoy: bool,
    pub sequence: Arc<[u8]>,
    /// Per-residue modification masses. An empty vector represents an
    /// unmodified peptide and is expanded lazily when a modification is applied.
    pub modifications: Vec<f32>,
    /// Identity and fragmentation metadata for applied modifications. Masses
    /// remain in the fields above for compatibility and fast mass arithmetic.
    pub applied_modifications: Arc<Vec<AppliedModification>>,
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
        self.modifications.get(index).copied().unwrap_or_default()
    }

    fn ensure_dense_modifications(&mut self) {
        if self.modifications.is_empty() {
            self.modifications.resize(self.sequence.len(), 0.0);
        }
    }

    pub fn initial_sort(&self, other: &Self) -> std::cmp::Ordering {
        self.sequence
            .cmp(&other.sequence)
            .then_with(|| {
                self.modifications
                    .partial_cmp(&other.modifications)
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
            .then_with(|| self.applied_modifications.cmp(&other.applied_modifications))
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
            .field("applied_modifications", &self.applied_modifications)
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
        if !self.applied_modifications.is_empty() {
            return self
                .applied_modifications
                .iter()
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
                .zip(self.modifications.iter())
                .filter(|(&r, &m)| resi == r && mass == m)
                .count(),
        }
    }

    fn modification_mass(&self) -> f32 {
        self.modifications.iter().sum::<f32>()
            + self.nterm.unwrap_or(0.0)
            + self.cterm.unwrap_or(0.0)
    }

    fn finalize_modifications(&mut self) {
        Arc::make_mut(&mut self.applied_modifications).sort_unstable();
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
                self.ensure_dense_modifications();
                if stack || self.modifications[index as usize] == 0.0 {
                    self.modifications[index as usize] += mass;
                    applied = true;
                }
            }
        }
        if applied {
            Arc::make_mut(&mut self.applied_modifications).push(AppliedModification {
                site,
                modification,
                kind,
            });
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
        let applied = Arc::make_mut(&mut self.applied_modifications);
        for modification in applied.iter_mut() {
            if modification.kind != ModificationKind::ChannelBase {
                continue;
            }
            let offset = modification.modification.channel_offsets[&channel];
            match modification.site {
                Site::Nterm => self.nterm = Some(self.nterm.unwrap_or_default() + offset),
                Site::Cterm => self.cterm = Some(self.cterm.unwrap_or_default() + offset),
                Site::Sequence(index) => self.modifications[index as usize] += offset,
            }
            modification.modification = cache
                .resolve(&modification.modification, &channel)
                .unwrap_or_else(|| {
                    Arc::new(
                        modification
                            .modification
                            .with_mass(modification.modification.mass + offset),
                    )
                });
            modification.kind = ModificationKind::Label;
        }
        self.monoisotopic += self.modification_mass() - before;
        self.label_channel = Some(channel);
        Arc::make_mut(&mut self.applied_modifications).sort_unstable();
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
        let label_modifications = Arc::make_mut(&mut base.applied_modifications);
        for applied in label_modifications.iter_mut() {
            if applied.kind != ModificationKind::Label {
                continue;
            }
            let offset = applied.modification.channel_offsets[&channel];
            match applied.site {
                Site::Nterm => base.nterm = nonzero_mass(base.nterm.unwrap_or_default() - offset),
                Site::Cterm => base.cterm = nonzero_mass(base.cterm.unwrap_or_default() - offset),
                Site::Sequence(index) => {
                    if let Some(mass) = base.modifications.get_mut(index as usize) {
                        *mass -= offset;
                        if mass.abs() < 1e-5 {
                            *mass = 0.0;
                        }
                    }
                }
            }
            base.monoisotopic -= offset;
            applied.modification = Arc::new(
                applied
                    .modification
                    .with_mass(applied.modification.mass - offset),
            );
            applied.kind = ModificationKind::ChannelBase;
        }
        base.label_channel = None;
        base.label_group_override = None;
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

    fn static_mods<M: ModificationSource>(
        &mut self,
        target: ModificationSpecificity,
        modification: &M,
    ) {
        let modification = modification.definition();
        match (target, self.position) {
            (ModificationSpecificity::PeptideN(None), _) => {
                self.apply_site(Site::Nterm, modification)
            }
            (ModificationSpecificity::PeptideN(Some(resi)), _)
                if resi == *self.sequence.first().unwrap_or(&0) =>
            {
                self.apply_site(Site::Sequence(0), modification)
            }
            (ModificationSpecificity::PeptideC(None), _) => {
                self.apply_site(Site::Cterm, modification)
            }
            (ModificationSpecificity::PeptideC(Some(resi)), _)
                if resi == *self.sequence.last().unwrap_or(&0) =>
            {
                self.apply_site(
                    Site::Sequence(self.sequence.len().saturating_sub(1) as u32),
                    modification,
                )
            }
            (ModificationSpecificity::ProteinN(None), Position::Nterm | Position::Full) => {
                self.apply_site(Site::Nterm, modification)
            }
            (ModificationSpecificity::ProteinN(Some(resi)), Position::Nterm | Position::Full)
                if resi == *self.sequence.first().unwrap_or(&0) =>
            {
                self.apply_site(Site::Sequence(0), modification)
            }
            (ModificationSpecificity::ProteinC(None), Position::Cterm | Position::Full) => {
                self.apply_site(Site::Cterm, modification)
            }
            (ModificationSpecificity::ProteinC(Some(resi)), Position::Cterm | Position::Full)
                if resi == *self.sequence.last().unwrap_or(&0) =>
            {
                self.apply_site(
                    Site::Sequence(self.sequence.len().saturating_sub(1) as u32),
                    modification,
                )
            }
            (ModificationSpecificity::Residue(resi), _) if self.sequence.contains(&resi) => {
                self.ensure_dense_modifications();
                let sites = self
                    .sequence
                    .iter()
                    .enumerate()
                    .filter_map(|(idx, residue)| {
                        (resi == *residue && self.modifications[idx] == 0.0)
                            .then_some(Site::Sequence(idx as u32))
                    })
                    .collect::<Vec<_>>();
                for site in sites {
                    self.apply_site(site, modification.clone());
                }
            }
            _ => {}
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
            .collect();
        self.apply_rules(
            &rules,
            &[],
            &static_mods,
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
        max_exhaustive_mods: usize,
        max_total_mods: usize,
        max_combinations: Option<usize>,
    ) -> Vec<Peptide> {
        if variable_mods.is_empty() {
            for (target, modification) in static_mods {
                self.static_mods(*target, modification);
            }
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
                for (target, modification) in static_mods {
                    peptide.static_mods(*target, modification);
                }
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
            if !pep.modifications.is_empty() {
                pep.modifications[1..n].reverse();
            }
            for applied in Arc::make_mut(&mut pep.applied_modifications) {
                if let Site::Sequence(index) = &mut applied.site {
                    let original = *index as usize;
                    if (1..n).contains(&original) {
                        *index = (n - original) as u32;
                    }
                }
            }
            Arc::make_mut(&mut pep.applied_modifications).sort_unstable();
        }
        pep
    }

    pub(crate) fn modification_tag(&self, site: Site, mass: f32) -> String {
        let applied = self
            .applied_modifications
            .iter()
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

#[derive(Copy, Clone, Debug, Default, PartialEq, Eq, PartialOrd, Ord)]
pub enum ModificationKind {
    #[default]
    Ordinary,
    ChannelBase,
    Label,
}

#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub enum PeptideError {
    InvalidSequence(String),
}

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
            // An empty vector is the compact representation for an unmodified peptide.
            // It is expanded lazily only when a modification is applied.
            modifications: Vec::new(),
            applied_modifications: Arc::default(),
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
