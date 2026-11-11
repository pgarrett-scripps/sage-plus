//! Peptides: the sequence, its modifications, and how a digest becomes one.

use std::cmp::Ordering;
use std::{collections::HashMap, fmt::Debug, sync::Arc};

use crate::modification::{ModificationDefinition, ModificationSpecificity, SiteMode};
use crate::{
    enzyme::{Digest, DigestGroup, Position, ProteinOccurrence},
    mass::{monoisotopic, H2O},
    sequence::PeptideSequence,
};
use itertools::Itertools;
use smallvec::SmallVec;

pub const INLINE_PROTEINS: usize = 1;

/// Most peptides map to one protein, so keep that accession inline.
/// Shared peptides transparently spill to heap storage.
pub type ProteinAccessions = SmallVec<[Arc<str>; INLINE_PROTEINS]>;

#[derive(Clone, PartialEq, Default)]
pub struct Peptide {
    pub decoy: bool,
    pub sequence: PeptideSequence,
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
        let length = self.sequence.len();
        let decode = |index| {
            if index == length {
                Site::Nterm
            } else if index == length + 1 {
                Site::Cterm
            } else {
                Site::Sequence(index as u32)
            }
        };
        let candidates = candidates
            .iter()
            .map(|&index| decode(index))
            .collect::<Vec<_>>();
        let selected = self
            .applied_modifications()
            .find(|applied| {
                candidates.contains(&applied.site)
                    && (applied.modification.mass - mass).abs() < epsilon
            })
            .map(|applied| (Arc::new(applied.modification.clone()), applied.kind));
        let Some((definition, kind)) = selected else {
            return;
        };
        let mut retained = self
            .applied_modifications()
            .filter(|applied| {
                !(candidates.contains(&applied.site) && applied.modification == definition.as_ref())
            })
            .map(|applied| AppliedModification {
                site: applied.site,
                modification: Arc::new(applied.modification.clone()),
                kind: applied.kind,
            })
            .collect::<Vec<_>>();
        retained.extend(chosen.iter().map(|&index| AppliedModification {
            site: decode(index),
            modification: definition.clone(),
            kind,
        }));
        self.nterm = None;
        self.cterm = None;
        for applied in &retained {
            match applied.site {
                Site::Nterm => {
                    self.nterm = Some(self.nterm.unwrap_or_default() + applied.modification.mass)
                }
                Site::Cterm => {
                    self.cterm = Some(self.cterm.unwrap_or_default() + applied.modification.mass)
                }
                _ => (),
            }
        }
        self.modifications = CompactModifications::from_applied(retained)
            .expect("localized modifications fit compact lookup");
    }

    /// Append every site compatible with `specificity`, ignoring occupancy.
    pub fn compatible_sites(&self, specificity: ModificationSpecificity, sites: &mut Vec<Site>) {
        sites.extend(self.rule_sites(specificity));
    }

    /// Return a copy carrying one additional modification at `site`. Used for
    /// search-time mass offsets, which are not part of the indexed expansion.
    pub fn with_mass_offset(
        &self,
        site: Site,
        definition: &Arc<ModificationDefinition>,
    ) -> Peptide {
        let mut peptide = self.clone();
        let lookup = ModificationLookup::from_definitions([(
            site,
            definition.clone(),
            ModificationKind::Ordinary,
        )])
        .unwrap_or_else(|error| panic!("{error}"));
        peptide.modifications.install_lookup(lookup);
        match site {
            Site::Nterm => {
                peptide.nterm = Some(peptide.nterm.unwrap_or_default() + definition.mass);
            }
            Site::Cterm => {
                peptide.cterm = Some(peptide.cterm.unwrap_or_default() + definition.mass);
            }
            Site::Sequence(_) => {}
        }
        peptide
            .modifications
            .push(site, definition, ModificationKind::Ordinary);
        peptide.modifications.sort();
        // Sum in the same order as database expansion so equivalent offset
        // and indexed peptidoforms have bit-identical masses and tie-breaks.
        // Decoys inherit the mass of their forward target.
        let target_order = if peptide.decoy && peptide.sequence.len() > 2 {
            peptide.sequence.reversed_internal()
        } else {
            peptide.sequence.clone()
        };
        peptide.monoisotopic = target_order
            .iter()
            .fold(H2O, |mass, residue| mass + monoisotopic(*residue))
            + peptide.modifications.total_mass();
        peptide
    }

    pub fn initial_sort(&self, other: &Self) -> std::cmp::Ordering {
        self.sequence
            .cmp(&other.sequence)
            .then_with(|| self.modification_sort(other))
    }

    /// The part of [`Self::initial_sort`] after the sequence: modifications
    /// by position, termini, then label channel.
    pub fn modification_sort(&self, other: &Self) -> std::cmp::Ordering {
        std::cmp::Ordering::Equal
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

    /// FASTA spans, as written, that this peptide was expanded from when they
    /// contain ambiguous residues (B, X or Z), joined by `;`. A generated
    /// decoy reports its target's span reversed the same way. `None` when no
    /// occurrence was expanded.
    pub fn database_peptide(&self) -> Option<String> {
        use crate::ambiguous_residues::{is_ambiguous, is_expansion_of};
        let len = self.sequence.len();
        let mut spans = self
            .protein_sites
            .iter()
            .filter_map(|site| {
                let start = site.start? as usize;
                let span = site.source.as_ref()?.as_bytes().get(start..start + len)?;
                if !is_ambiguous(span) {
                    return None;
                }
                if is_expansion_of(span, &self.sequence) {
                    return Some(span.to_vec());
                }
                let reversed = crate::sequence::reverse_interior(span);
                is_expansion_of(&reversed, &self.sequence).then_some(reversed)
            })
            .collect::<Vec<_>>();
        if spans.is_empty() {
            return None;
        }
        spans.sort_unstable();
        spans.dedup();
        Some(
            spans
                .iter()
                .map(|span| String::from_utf8_lossy(span))
                .join(";"),
        )
    }

    /// Ambiguous FASTA residues (B, Z or X) replaced to make this peptide,
    /// formatted as `X4K;B7D` (see
    /// [`crate::ambiguous_residues::format_substitutions`]); empty when there
    /// are none.
    ///
    /// Occurrences are read in protein order, the order of `proteins` and
    /// `protein_sites`, which also picks the displayed sequence of merged
    /// I/L twins. If any occurrence has the peptide's residues as written,
    /// the result is empty: a FASTA span without B, Z or X, a site without a
    /// FASTA span, or a protein listed without a site (a peptide TSV row
    /// without B, Z or X). Otherwise the first expanded
    /// occurrence gives the substitutions. A generated decoy reports its
    /// target span reversed like the decoy, at decoy positions.
    pub fn substitutions(&self) -> String {
        use crate::ambiguous_residues::{format_substitutions, is_ambiguous, oriented_span};
        let len = self.sequence.len();
        let mut expanded = None;
        for site in self.protein_sites.iter() {
            // A site without a start is an expanded peptide TSV row, whose
            // source is the row as written.
            let span = site.source.as_ref().and_then(|source| match site.start {
                Some(start) => source.as_bytes().get(start as usize..start as usize + len),
                None => Some(source.as_bytes()).filter(|row| row.len() == len),
            });
            match span {
                Some(span) if is_ambiguous(span) => {
                    if expanded.is_none() {
                        expanded = oriented_span(span, &self.sequence);
                    }
                }
                _ => return String::new(),
            }
        }
        let Some(written) = expanded else {
            return String::new();
        };
        let listed_without_site = self.proteins.iter().any(|protein| {
            !self
                .protein_sites
                .iter()
                .any(|site| site.protein == *protein)
        });
        if listed_without_site {
            return String::new();
        }
        format_substitutions(&written, &self.sequence)
    }

    pub fn modification_count(&self, target: ModificationSpecificity, mass: f32) -> usize {
        let sites = self.rule_sites(target);
        if !self.modifications.is_empty() {
            return self
                .applied_modifications()
                .filter(|applied| {
                    applied.modification.mass == mass && sites.contains(&applied.site)
                })
                .count();
        }
        sites
            .into_iter()
            .filter(|site| match site {
                Site::Nterm => self.nterm == Some(mass),
                Site::Cterm => self.cterm == Some(mass),
                Site::Sequence(i) => self.modification_at(*i as usize) == mass,
            })
            .count()
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
        match self.label_group_peptide() {
            Some(base) => base.to_string(),
            None => self.to_string(),
        }
    }

    /// The channel-collapsed peptidoform behind [`Peptide::label_group`]:
    /// every precursor-label modification resolved back to its channel base,
    /// with the channel offset removed from the terminal masses. `None` when
    /// this peptide carries no precursor label channel, in which case the
    /// peptide is already its own label group.
    ///
    /// Callers that need a structural identity rather than a rendered string
    /// (peptide-level FDR competition) use this instead of parsing
    /// `label_group`.
    pub fn label_group_peptide(&self) -> Option<Peptide> {
        let channel = self.label_channel.clone()?;
        let mut base = self.clone();
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
            let unresolved = record
                .base
                .as_ref()
                .expect("label definition has no recorded channel base");
            encoded.modification_id = lookup
                .id(record.site, unresolved, ModificationKind::ChannelBase)
                .expect("base label definition missing from compact lookup");
        }
        base.label_channel = None;
        base.label_group_override = None;
        base.modifications.sort();
        Some(base)
    }

    /// Physical sites selected by one rule. Motif rules are evaluated against
    /// the source protein of every occurrence; see
    /// [`ModificationSpecificity::sites_for_occurrences`].
    pub fn rule_sites(&self, target: ModificationSpecificity) -> Vec<Site> {
        target.sites_for_occurrences(
            &self.sequence,
            self.position,
            self.decoy,
            &self.protein_sites,
        )
    }

    fn push_resi(
        &self,
        acc: &mut Vec<(Site, f32, usize)>,
        target: ModificationSpecificity,
        mass: f32,
        mod_idx: usize,
    ) {
        acc.extend(
            self.rule_sites(target)
                .into_iter()
                .map(|site| (site, mass, mod_idx)),
        );
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
            let sites = self
                .rule_sites(*target)
                .into_iter()
                .filter(|site| {
                    let occupied = match site {
                        Site::Nterm => &mut nterm_occupied,
                        Site::Cterm => &mut cterm_occupied,
                        Site::Sequence(i) => &mut occupied[*i as usize],
                    };
                    let available = !*occupied;
                    *occupied = true;
                    available
                })
                .collect::<Vec<_>>();
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
                    max_total_count: None,
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
        let plan = ModificationPlan::new(
            &rules,
            &static_mods,
            lookup,
            combinations,
            combinations,
            max_combinations,
        );
        self.apply_rules(&plan, &[])
    }

    /// Apply config-defined and library-supported variable modifications in a
    /// single enumeration. Library-supported placements do not consume the
    /// exhaustive budget or the named modification's `max_count`; all
    /// placements consume the total budget and its `max_total_count`.
    pub(crate) fn apply_rules(
        mut self,
        plan: &ModificationPlan<'_>,
        library_sites: &[LibrarySite],
    ) -> Vec<Peptide> {
        self.modifications.install_lookup(plan.lookup.clone());
        if plan.variable_mods.is_empty() {
            self.static_mods_all(plan.static_mods);
            self.finalize_modifications();
            vec![self]
        } else {
            let mut candidates: Vec<ModificationCandidate> = Vec::new();
            for (rule_idx, rule) in plan.variable_mods.iter().enumerate() {
                let mut compatible = Vec::new();
                self.push_resi(
                    &mut compatible,
                    rule.specificity,
                    rule.modification.mass,
                    rule_idx,
                );
                for (site, _, _) in compatible {
                    let library_supported = rule.site_mode != SiteMode::Exhaustive
                        && rule.modification.name.as_deref().is_some_and(|name| {
                            library_sites.iter().any(|library| {
                                library.modification.as_ref() == name
                                    && library.attachment.site(
                                        library.position,
                                        self.sequence.len(),
                                        self.position,
                                    ) == Some(site)
                            })
                        });
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
                let key = (
                    candidate.site,
                    plan.variable_mods[candidate.rule].count_group,
                );
                if let Some(index) = candidate_indices.get(&key).copied() {
                    mods[index].library_supported |= candidate.library_supported;
                } else {
                    candidate_indices.insert(key, mods.len());
                    mods.push(candidate);
                }
            }

            let group_count = plan
                .variable_mods
                .iter()
                .map(|rule| rule.count_group)
                .max()
                .map_or(0, |group| group + 1);
            let mut enumeration = ModificationEnumeration::new(&self, &mods, plan, group_count);
            for n in 1..=plan.max_total_mods.min(mods.len()) {
                enumeration.selected.reserve(n);
                if !enumeration.enumerate(0, n, 0) {
                    break;
                }
            }
            let mut modified = enumeration.output;

            // Apply static mods to all peptides
            for peptide in modified.iter_mut() {
                peptide.static_mods_all(plan.static_mods);
                peptide.finalize_modifications();
            }

            modified
        }
    }

    pub fn reverse(&self) -> Peptide {
        self.reverse_with_sequence(self.sequence.reversed_internal())
    }

    pub(crate) fn reverse_with_sequence(&self, sequence: PeptideSequence) -> Peptide {
        let mut pep = self.clone();
        pep.decoy = !self.decoy;
        let n = pep.sequence.len().saturating_sub(1);
        if n > 1 {
            pep.sequence = sequence;
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

mod enumerate;
mod lookup;

pub use enumerate::*;
pub use lookup::*;

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
                let mut occurrence = ProteinOccurrence::of(&value);
                occurrence.start = Some(start);
                vec![occurrence].into()
            })
            .unwrap_or_default();
        // This is an important invariant to enforce, that ensures safety
        // while reversing peptide sequences
        if !value.sequence.is_ascii() {
            return Err(PeptideError::InvalidSequence(
                String::from_utf8_lossy(value.sequence.as_bytes()).into_owned(),
            ));
        }
        for c in value.sequence.iter() {
            let mono = monoisotopic(*c);
            if mono == 0.0 {
                return Err(PeptideError::InvalidSequence(
                    String::from_utf8_lossy(value.sequence.as_bytes()).into_owned(),
                ));
            }
            mass += mono;
        }

        Ok(Peptide {
            decoy: value.decoy,
            position: value.position,
            modifications: CompactModifications::default(),
            label_channel: None,
            label_group_override: None,
            sequence: value.sequence,
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
