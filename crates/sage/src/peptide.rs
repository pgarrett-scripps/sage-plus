use std::cmp::Ordering;
use std::{collections::HashMap, fmt::Debug, sync::Arc};

use crate::modification::{ModificationDefinition, ModificationSpecificity, SiteMode};
use crate::{
    enzyme::{Digest, DigestGroup, Position},
    mass::{monoisotopic, H2O},
};
use fnv::FnvHashSet;
use itertools::Itertools;

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

    pub proteins: Vec<Arc<str>>,
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
    pub(crate) fn apply_label_channel(mut self, channel: Arc<str>) -> Self {
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
            modification.modification = Arc::new(
                modification
                    .modification
                    .with_mass(modification.modification.mass + offset),
            );
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
            proteins: vec![value.protein],
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
mod test {
    use crate::enzyme::{Digest, Enzyme, EnzymeParameters};
    use crate::modification::NeutralLossMode;

    use super::*;

    #[test]
    fn unmodified_peptides_use_compact_modification_storage() {
        let peptide = Peptide::try_from(Digest {
            sequence: "PEPTIDER".into(),
            ..Digest::default()
        })
        .unwrap();
        assert!(peptide.modifications.is_empty());
        assert_eq!(peptide.modification_at(3), 0.0);

        let modified = peptide
            .apply(
                &[(ModificationSpecificity::Residue(b'P'), 10.0, None)],
                &HashMap::default(),
                1,
                None,
            )
            .into_iter()
            .find(|peptide| {
                peptide.modification_count(ModificationSpecificity::Residue(b'P'), 10.0) > 0
            })
            .unwrap();
        assert_eq!(modified.modifications.len(), modified.sequence.len());
    }

    fn detailed_mod(
        mass: f32,
        name: &str,
        neutral_losses: &[f32],
        neutral_loss_mode: NeutralLossMode,
    ) -> Arc<ModificationDefinition> {
        Arc::new(ModificationDefinition {
            mass,
            name: Some(Arc::from(name)),
            neutral_losses: Arc::from(neutral_losses),
            neutral_loss_mode,
            channel_offsets: Arc::default(),
        })
    }

    fn var_mod_sequence(
        peptide: &Peptide,
        mods: &[(ModificationSpecificity, f32)],
        combo: usize,
    ) -> Vec<String> {
        let static_mods = HashMap::default();
        let mods_with_limits: Vec<(ModificationSpecificity, f32, Option<usize>)> =
            mods.iter().map(|&(s, m)| (s, m, None)).collect();
        peptide
            .clone()
            .apply(&mods_with_limits, &static_mods, combo, None)
            .into_iter()
            .map(|p| p.to_string())
            .collect::<Vec<_>>()
    }

    #[test]
    fn full() {
        let sequence = "MPEPTIDEKMSAGEKEND";
        let tryp = EnzymeParameters {
            min_len: 0,
            max_len: 50,
            missed_cleavages: 0,
            enzyme: Enzyme::new("KR", "P", true, false),
        };

        let peptides = tryp
            .digest(sequence, Default::default())
            .into_iter()
            .map(Peptide::try_from)
            .collect::<Result<Vec<_>, _>>()
            .unwrap();

        assert_eq!(peptides.len(), 3);
        assert_eq!(peptides[0].to_string(), "MPEPTIDEK");
        assert_eq!(peptides[0].position, Position::Nterm);
        assert_eq!(peptides[1].to_string(), "MSAGEK");
        assert_eq!(peptides[1].position, Position::Internal);
        assert_eq!(peptides[2].to_string(), "END");
        assert_eq!(peptides[2].position, Position::Cterm);

        use ModificationSpecificity::*;

        let mods = [
            (ProteinN(None), 42.0),
            (ProteinC(None), 11.0),
            (PeptideN(None), 12.0),
            (PeptideC(None), 19.0),
        ];
        let a = var_mod_sequence(&peptides[0], &mods, 2);
        let b = var_mod_sequence(&peptides[1], &mods, 2);
        let c = var_mod_sequence(&peptides[2], &mods, 2);

        // Make sure no duplicates exist
        assert_eq!(
            a,
            vec![
                "MPEPTIDEK",
                "[+42]-MPEPTIDEK",
                "[+12]-MPEPTIDEK",
                "MPEPTIDEK-[+19]",
                "[+42]-MPEPTIDEK-[+19]",
                "[+12]-MPEPTIDEK-[+19]",
            ]
        );

        assert_eq!(
            b,
            vec![
                "MSAGEK",
                "[+12]-MSAGEK",
                "MSAGEK-[+19]",
                "[+12]-MSAGEK-[+19]",
            ]
        );

        assert_eq!(
            c,
            vec![
                "END",
                "END-[+11]",
                "[+12]-END",
                "END-[+19]",
                "[+12]-END-[+11]",
                "[+12]-END-[+19]",
            ]
        );
    }

    #[test]
    fn test_variable_mods() {
        use ModificationSpecificity::*;
        let variable_mods = [(Residue(b'M'), 16.0f32), (Residue(b'C'), 57.)];
        let peptide = Peptide::try_from(Digest {
            sequence: "GCMGCMG".into(),
            ..Default::default()
        })
        .unwrap();

        let expected = vec![
            "GCMGCMG",
            "GCM[+16]GCMG",
            "GCMGCM[+16]G",
            "GC[+57]MGCMG",
            "GCMGC[+57]MG",
            "GCM[+16]GCM[+16]G",
            "GC[+57]M[+16]GCMG",
            "GCM[+16]GC[+57]MG",
            "GC[+57]MGCM[+16]G",
            "GCMGC[+57]M[+16]G",
            "GC[+57]MGC[+57]MG",
        ];

        let peptides = var_mod_sequence(&peptide, &variable_mods, 2);
        assert_eq!(peptides, expected);
    }

    #[test]
    fn test_variable_mods_no_effeect() {
        use ModificationSpecificity::*;
        let variable_mods = [(Residue(b'M'), 16.0f32), (Residue(b'C'), 57.)];
        let peptide = Peptide::try_from(Digest {
            sequence: "AAAAAAAA".into(),
            ..Default::default()
        })
        .unwrap();

        let expected = vec!["AAAAAAAA"];
        let peptides = var_mod_sequence(&peptide, &variable_mods, usize::MAX);
        assert_eq!(peptides, expected);
    }

    #[test]
    fn test_variable_mods_nterm() {
        use ModificationSpecificity::*;
        let variable_mods = [(PeptideN(None), 42.), (Residue(b'M'), 16.)];
        let peptide = Peptide::try_from(Digest {
            sequence: "GCMGCMG".into(),
            ..Default::default()
        })
        .unwrap();

        let expected = vec![
            "GCMGCMG",
            "[+42]-GCMGCMG",
            "GCM[+16]GCMG",
            "GCMGCM[+16]G",
            "[+42]-GCM[+16]GCMG",
            "[+42]-GCMGCM[+16]G",
            "GCM[+16]GCM[+16]G",
            "[+42]-GCM[+16]GCM[+16]G",
        ];

        let peptides = var_mod_sequence(&peptide, &variable_mods, 3);
        assert_eq!(peptides, expected);
    }

    #[test]
    fn test_variable_mods_cterm() {
        use ModificationSpecificity::*;
        let variable_mods = [(PeptideC(None), 42.), (Residue(b'M'), 16.)];
        let peptide = Peptide::try_from(Digest {
            sequence: "GCMGCMG".into(),
            ..Default::default()
        })
        .unwrap();

        let expected = vec![
            "GCMGCMG",
            "GCMGCMG-[+42]",
            "GCM[+16]GCMG",
            "GCMGCM[+16]G",
            "GCM[+16]GCMG-[+42]",
            "GCMGCM[+16]G-[+42]",
            "GCM[+16]GCM[+16]G",
            "GCM[+16]GCM[+16]G-[+42]",
        ];

        let peptides = var_mod_sequence(&peptide, &variable_mods, 3);
        assert_eq!(peptides, expected);
    }

    #[test]
    fn test_variable_mods_multi() {
        use ModificationSpecificity::*;
        let variable_mods = [(Residue(b'S'), 79.), (Residue(b'S'), 541.)];
        let peptide = Peptide::try_from(Digest {
            sequence: "GGGSGGGS".into(),
            ..Default::default()
        })
        .unwrap();

        let expected = vec![
            "GGGSGGGS",
            "GGGS[+79]GGGS",
            "GGGSGGGS[+79]",
            "GGGS[+541]GGGS",
            "GGGSGGGS[+541]",
            "GGGS[+79]GGGS[+79]",
            "GGGS[+79]GGGS[+541]",
            "GGGS[+541]GGGS[+79]",
            "GGGS[+541]GGGS[+541]",
        ];

        let peptides = var_mod_sequence(&peptide, &variable_mods, 2);
        assert_eq!(peptides, expected);
    }

    /// Check that picked-peptide approach will match forward and reverse peptides
    #[test]
    fn test_psuedo_forward() {
        let trypsin = crate::enzyme::EnzymeParameters {
            missed_cleavages: 0,
            min_len: 3,
            max_len: 30,
            enzyme: Enzyme::new("KR", "P", true, false),
        };

        let fwd = "MADEEKLPPGWEKRMSRSSGRVYYFNHITNASQWERPSGN";
        for digest in trypsin.digest(fwd, Default::default()) {
            let fwd = Peptide::try_from(digest.clone()).unwrap();
            let rev = Peptide::try_from(digest.reverse()).unwrap();

            assert_eq!(fwd.decoy, false);
            assert_eq!(rev.decoy, true);
            assert!(
                fwd.sequence.len() < 4 || fwd.sequence != rev.sequence,
                "{} {}",
                fwd,
                rev
            );
            assert_eq!(rev.reverse().to_string(), fwd.to_string());
        }
    }

    #[test]
    fn apply_mods() {
        use ModificationSpecificity::*;
        let peptide = Peptide::try_from(Digest {
            sequence: "AACAACAA".into(),
            ..Default::default()
        })
        .unwrap();

        let expected = vec![
            "AAC[+57]AAC[+57]AA",
            "AAC[+30]AAC[+57]AA",
            "AAC[+57]AAC[+30]AA",
            "AAC[+30]AAC[+30]AA",
        ];

        let mut static_mods = HashMap::new();
        static_mods.insert(Residue(b'C'), 57.0);

        let variable_mods = [(Residue(b'C'), 30.0, None)];

        let peptides = peptide
            .apply(&variable_mods, &static_mods, 2, None)
            .into_iter()
            .map(|p| p.to_string())
            .collect::<Vec<_>>();

        assert_eq!(peptides, expected);
    }

    #[test]
    fn test_per_mod_limit() {
        use ModificationSpecificity::*;
        // GCMGCMG has two M residues; limit oxidation to max 1 per peptide
        let variable_mods = [(Residue(b'M'), 16.0f32, Some(1))];
        let peptide = Peptide::try_from(Digest {
            sequence: "GCMGCMG".into(),
            ..Default::default()
        })
        .unwrap();

        let static_mods = HashMap::default();
        let peptides: Vec<String> = peptide
            .clone()
            .apply(&variable_mods, &static_mods, 2, None)
            .into_iter()
            .map(|p| p.to_string())
            .collect();

        // Should get unmodified + each single-M variant, but NOT the double-M variant
        let expected = vec!["GCMGCMG", "GCM[+16]GCMG", "GCMGCM[+16]G"];
        assert_eq!(peptides, expected);
    }

    #[test]
    fn test_max_combinations() {
        use ModificationSpecificity::*;
        // GCMGCMG with oxidation and carbamidomethylation would normally yield many variants;
        // cap at 4 total (unmodified + 3 modified)
        let variable_mods = [(Residue(b'M'), 16.0f32, None), (Residue(b'C'), 57.0, None)];
        let peptide = Peptide::try_from(Digest {
            sequence: "GCMGCMG".into(),
            ..Default::default()
        })
        .unwrap();

        let static_mods = HashMap::default();
        let peptides: Vec<String> = peptide
            .clone()
            .apply(&variable_mods, &static_mods, 2, Some(4))
            .into_iter()
            .map(|p| p.to_string())
            .collect();

        // Cap at 4: unmodified + the first 3 single-mod variants (fewest PTMs first)
        assert_eq!(peptides.len(), 4);
        assert_eq!(peptides[0], "GCMGCMG");
    }

    #[test]
    fn modification_sites() {
        use Site::*;
        let peptide = Peptide::try_from(Digest {
            sequence: "AACAACAA".into(),
            ..Default::default()
        })
        .unwrap();

        let mut mods = vec![];
        peptide.push_resi(&mut mods, ModificationSpecificity::Residue(b'C'), 16.0, 0);
        assert_eq!(mods, vec![(Sequence(2), 16.0, 0), (Sequence(5), 16.0, 0)]);
        mods.clear();

        peptide.push_resi(&mut mods, ModificationSpecificity::PeptideC(None), 16.0, 0);
        assert_eq!(mods, vec![(Cterm, 16.0, 0)]);
        mods.clear();

        peptide.push_resi(&mut mods, ModificationSpecificity::PeptideN(None), 16.0, 0);
        assert_eq!(mods, vec![(Nterm, 16.0, 0)]);
        mods.clear();

        let mut mods = vec![];
        for (idx, (residue, mass)) in [("^", 12.0), ("$", 200.0), ("C", 57.0), ("A", 43.0)]
            .iter()
            .enumerate()
        {
            peptide.push_resi(&mut mods, residue.parse().unwrap(), *mass, idx);
        }

        assert_eq!(
            mods,
            vec![
                (Nterm, 12.0, 0),
                (Cterm, 200.0, 1),
                (Sequence(2), 57.0, 2),
                (Sequence(5), 57.0, 2),
                (Sequence(0), 43.0, 3),
                (Sequence(1), 43.0, 3),
                (Sequence(3), 43.0, 3),
                (Sequence(4), 43.0, 3),
                (Sequence(6), 43.0, 3),
                (Sequence(7), 43.0, 3),
            ]
        );
    }

    #[test]
    fn test_per_mod_limit_exactly_met() {
        use ModificationSpecificity::*;
        // Limit of 2 on a peptide with exactly 2 M residues — all combos should be allowed
        let variable_mods = [(Residue(b'M'), 16.0f32, Some(2))];
        let peptide = Peptide::try_from(Digest {
            sequence: "GCMGCMG".into(),
            ..Default::default()
        })
        .unwrap();

        let static_mods = HashMap::default();
        let peptides: Vec<String> = peptide
            .clone()
            .apply(&variable_mods, &static_mods, 2, None)
            .into_iter()
            .map(|p| p.to_string())
            .collect();

        // No restriction: unmodified + each single + double
        let expected = vec![
            "GCMGCMG",
            "GCM[+16]GCMG",
            "GCMGCM[+16]G",
            "GCM[+16]GCM[+16]G",
        ];
        assert_eq!(peptides, expected);
    }

    #[test]
    fn test_per_mod_limit_zero() {
        use ModificationSpecificity::*;
        // Limit of 0 means this mod is entirely suppressed
        let variable_mods = [(Residue(b'M'), 16.0f32, Some(0))];
        let peptide = Peptide::try_from(Digest {
            sequence: "GCMGCMG".into(),
            ..Default::default()
        })
        .unwrap();

        let static_mods = HashMap::default();
        let peptides: Vec<String> = peptide
            .clone()
            .apply(&variable_mods, &static_mods, 2, None)
            .into_iter()
            .map(|p| p.to_string())
            .collect();

        assert_eq!(peptides, vec!["GCMGCMG"]);
    }

    #[test]
    fn test_mixed_limited_and_unlimited() {
        use ModificationSpecificity::*;
        // M oxidation limited to 1; C carbamidomethylation unlimited
        // GCMGCMG has 2 M and 2 C
        let variable_mods = [
            (Residue(b'M'), 16.0f32, Some(1)),
            (Residue(b'C'), 57.0f32, None),
        ];
        let peptide = Peptide::try_from(Digest {
            sequence: "GCMGCMG".into(),
            ..Default::default()
        })
        .unwrap();

        let static_mods = HashMap::default();
        let peptides: Vec<String> = peptide
            .clone()
            .apply(&variable_mods, &static_mods, 2, None)
            .into_iter()
            .map(|p| p.to_string())
            .collect();

        // Should include all combos with ≤1 oxidized M,
        // but never both M residues oxidized simultaneously
        for p in &peptides {
            let oxid_count = p.matches("[+16]").count();
            assert!(oxid_count <= 1, "too many oxidations in: {}", p);
        }
        // Both C residues carbamidomethylated simultaneously should be present
        assert!(
            peptides.contains(&"GC[+57]MGC[+57]MG".to_string()),
            "expected double-C mod"
        );
        // Double oxidation should be absent
        assert!(
            !peptides.contains(&"GCM[+16]GCM[+16]G".to_string()),
            "double oxidation should be suppressed"
        );
    }

    #[test]
    fn test_limits_are_per_mod_not_per_residue() {
        use ModificationSpecificity::*;
        // Both modifications target M, but only oxidation is limited to one.
        let variable_mods = [
            (Residue(b'M'), 16.0f32, Some(1)),
            (Residue(b'M'), 32.0f32, None),
        ];
        let peptide = Peptide::try_from(Digest {
            sequence: "GCMGCMG".into(),
            ..Default::default()
        })
        .unwrap();

        let peptides = peptide.apply(&variable_mods, &HashMap::default(), 2, None);
        let peptides = peptides.iter().map(ToString::to_string).collect::<Vec<_>>();

        assert!(!peptides.contains(&"GCM[+16]GCM[+16]G".to_string()));
        assert!(peptides.contains(&"GCM[+32]GCM[+32]G".to_string()));
    }

    #[test]
    fn test_limits_support_more_than_64_mod_entries() {
        use ModificationSpecificity::*;
        let mut variable_mods = (1..=65)
            .map(|mass| (Residue(b'M'), mass as f32, None))
            .collect::<Vec<_>>();
        variable_mods[64].2 = Some(0);
        let peptide = Peptide::try_from(Digest {
            sequence: "GMG".into(),
            ..Default::default()
        })
        .unwrap();

        let peptides = peptide.apply(&variable_mods, &HashMap::default(), 1, None);

        assert_eq!(peptides.len(), 65); // unmodified + 64 allowed entries
        assert!(!peptides
            .iter()
            .any(|peptide| peptide.to_string().contains("[+65]")));
    }

    #[test]
    fn test_max_combinations_only_unmodified() {
        use ModificationSpecificity::*;
        // cap of 1 means only the unmodified peptide is returned
        let variable_mods = [(Residue(b'M'), 16.0f32, None)];
        let peptide = Peptide::try_from(Digest {
            sequence: "GCMGCMG".into(),
            ..Default::default()
        })
        .unwrap();

        let static_mods = HashMap::default();
        let peptides: Vec<String> = peptide
            .clone()
            .apply(&variable_mods, &static_mods, 2, Some(1))
            .into_iter()
            .map(|p| p.to_string())
            .collect();

        assert_eq!(peptides, vec!["GCMGCMG"]);
    }

    #[test]
    fn test_max_combinations_prefers_fewer_ptms() {
        use ModificationSpecificity::*;
        // GCMGCMG with oxidation (2 sites) — normally 3 variants (unmod + 2 single + 1 double)
        // cap at 3 means we get unmod + both singles but not the double
        let variable_mods = [(Residue(b'M'), 16.0f32, None)];
        let peptide = Peptide::try_from(Digest {
            sequence: "GCMGCMG".into(),
            ..Default::default()
        })
        .unwrap();

        let static_mods = HashMap::default();
        let peptides: Vec<String> = peptide
            .clone()
            .apply(&variable_mods, &static_mods, 2, Some(3))
            .into_iter()
            .map(|p| p.to_string())
            .collect();

        assert_eq!(peptides, vec!["GCMGCMG", "GCM[+16]GCMG", "GCMGCM[+16]G"]);
        // Double-mod must not appear — it would require cap > 3
        assert!(!peptides.contains(&"GCM[+16]GCM[+16]G".to_string()));
    }

    #[test]
    fn names_follow_exact_modification_identity_and_decoy_position() {
        use ModificationSpecificity::*;
        let peptide = Peptide::try_from(Digest {
            sequence: "AMMAK".into(),
            ..Default::default()
        })
        .unwrap();
        let mods = [
            (
                Residue(b'M'),
                detailed_mod(15.9949, "Oxidation", &[], NeutralLossMode::Optional),
                Some(1),
            ),
            (
                Residue(b'M'),
                detailed_mod(15.9949, "AlternateName", &[], NeutralLossMode::Optional),
                Some(1),
            ),
        ];
        let peptides = peptide.apply(&mods, &HashMap::default(), 1, None);
        let rendered = peptides.iter().map(ToString::to_string).collect::<Vec<_>>();

        assert!(rendered.contains(&"AM[Oxidation]MAK".to_string()));
        assert!(rendered.contains(&"AM[AlternateName]MAK".to_string()));
        assert_ne!(
            peptides[1].applied_modifications,
            peptides[3].applied_modifications
        );

        let named = peptides
            .iter()
            .find(|peptide| peptide.to_string() == "AM[Oxidation]MAK")
            .unwrap();
        assert!(named.reverse().to_string().contains("[Oxidation]"));
    }

    #[test]
    fn library_and_exhaustive_candidates_are_enumerated_together() {
        let peptide = Peptide::try_from(Digest {
            sequence: "MSS".into(),
            ..Default::default()
        })
        .unwrap();
        let phospho = Arc::new(ModificationDefinition {
            mass: 79.96633,
            name: Some(Arc::from("Phospho")),
            neutral_losses: Arc::from([]),
            neutral_loss_mode: NeutralLossMode::Optional,
            channel_offsets: Arc::default(),
        });
        let oxidation = Arc::new(ModificationDefinition {
            mass: 15.9949,
            name: Some(Arc::from("Oxidation")),
            neutral_losses: Arc::from([]),
            neutral_loss_mode: NeutralLossMode::Optional,
            channel_offsets: Arc::default(),
        });
        let rules = vec![
            VariableRule {
                specificity: ModificationSpecificity::Residue(b'S'),
                modification: phospho,
                max_count: Some(2),
                site_mode: SiteMode::Both,
                count_group: 0,
            },
            VariableRule {
                specificity: ModificationSpecificity::Residue(b'M'),
                modification: oxidation,
                max_count: Some(1),
                site_mode: SiteMode::Both,
                count_group: 1,
            },
        ];
        let library = vec![
            LibrarySite {
                position: 1,
                modification: Arc::from("Phospho"),
            },
            LibrarySite {
                position: 2,
                modification: Arc::from("Phospho"),
            },
        ];

        let variants = peptide.apply_rules(&rules, &library, &HashMap::new(), 1, 3, None);

        assert!(variants.iter().any(|peptide| {
            peptide.modification_at(0) != 0.0
                && peptide.modification_at(1) != 0.0
                && peptide.modification_at(2) != 0.0
        }));
        assert!(!variants.iter().any(|peptide| {
            peptide.modification_at(0) != 0.0
                && peptide.modification_at(1) == 0.0
                && peptide.modification_at(2) == 0.0
                && peptide.applied_modifications.len() > 1
        }));
    }

    #[test]
    fn named_max_count_is_shared_across_residue_rules() {
        let peptide = Peptide::try_from(Digest {
            sequence: "ST".into(),
            ..Default::default()
        })
        .unwrap();
        let phospho = Arc::new(ModificationDefinition {
            mass: 79.96633,
            name: Some(Arc::from("Phospho")),
            neutral_losses: Arc::from([]),
            neutral_loss_mode: NeutralLossMode::Optional,
            channel_offsets: Arc::default(),
        });
        let rules = (*b"ST").map(|residue| VariableRule {
            specificity: ModificationSpecificity::Residue(residue),
            modification: phospho.clone(),
            max_count: Some(1),
            site_mode: SiteMode::Library,
            count_group: 0,
        });
        let library = vec![
            LibrarySite {
                position: 0,
                modification: Arc::from("Phospho"),
            },
            LibrarySite {
                position: 1,
                modification: Arc::from("Phospho"),
            },
        ];

        let variants = peptide.apply_rules(&rules, &library, &HashMap::new(), 0, 2, None);
        assert_eq!(variants.len(), 3);
        assert!(variants
            .iter()
            .all(|peptide| peptide.applied_modifications.len() <= 1));
    }

    #[test]
    fn static_modification_names_are_rendered() {
        let peptide = Peptide::try_from(Digest {
            sequence: "ACK".into(),
            ..Default::default()
        })
        .unwrap();
        let static_mods = HashMap::from([(
            ModificationSpecificity::Residue(b'C'),
            detailed_mod(57.0215, "Carbamidomethyl", &[], NeutralLossMode::Optional),
        )]);

        let peptides = peptide.apply(&[], &static_mods, 0, None);
        assert_eq!(peptides[0].to_string(), "AC[Carbamidomethyl]K");
    }
}
