use crate::ambiguous_residues;
use crate::cleavage::ValidatedCustomCleavageLibrary;
use crate::enzyme::{Digest, EnzymeParameters};
use crate::sequence::ProteinSequence;
use rayon::prelude::*;
use std::sync::Arc;

/// Which FASTA a protein was read from. The class is a second axis next to
/// target/decoy: a decoy is a decoy whatever file it came from, and only
/// targets are split by class.
#[derive(Copy, Clone, Debug, Default, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum ProteinClass {
    /// `database.fasta`, the sample proteins.
    #[default]
    Sample,
    /// `database.contaminant_fasta`.
    Contaminant,
    /// `database.entrapment_fasta`.
    Entrapment,
}

/// The classes of the proteins a peptide was digested from. A peptide can
/// belong to proteins of more than one class, so the flags are a set rather
/// than one class, and the rules read it:
///
/// - contaminant: **any** of its proteins is a contaminant;
/// - entrapment: **every** one of its proteins is entrapment.
///
/// One byte, so that carrying it on every [`crate::peptide::Peptide`] of a
/// database of tens of millions costs nothing: it fits the padding the
/// peptide's other one-byte fields already leave.
#[derive(Copy, Clone, Debug, Default, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct ProteinClassSet(u8);

impl ProteinClassSet {
    const SAMPLE: u8 = 1;
    const CONTAMINANT: u8 = 1 << 1;
    const ENTRAPMENT: u8 = 1 << 2;

    const fn flag(class: ProteinClass) -> u8 {
        match class {
            ProteinClass::Sample => Self::SAMPLE,
            ProteinClass::Contaminant => Self::CONTAMINANT,
            ProteinClass::Entrapment => Self::ENTRAPMENT,
        }
    }

    /// The set holding one class.
    pub const fn of(class: ProteinClass) -> Self {
        Self(Self::flag(class))
    }

    /// Add a protein's class to the set.
    pub fn insert(&mut self, class: ProteinClass) {
        self.0 |= Self::flag(class);
    }

    /// Add every class of another set, as merging two copies of one peptide
    /// does.
    pub fn union(&mut self, other: Self) {
        self.0 |= other.0;
    }

    /// Whether any protein of this peptide is a contaminant.
    pub const fn is_contaminant(self) -> bool {
        self.0 & Self::CONTAMINANT != 0
    }

    /// Whether every protein of this peptide is entrapment. An empty set is
    /// not entrapment.
    pub const fn is_entrapment(self) -> bool {
        self.0 == Self::ENTRAPMENT
    }

    /// Whether any class at all was recorded.
    pub const fn is_empty(self) -> bool {
        self.0 == 0
    }
}

/// The accessions of every protein whose class is not [`ProteinClass::Sample`],
/// so a peptide's classes can be read from its protein accessions. Keeping
/// them beside the database rather than on every peptide costs one set of
/// accessions instead of a byte on each of a database's tens of millions of
/// peptides, which would round [`crate::peptide::Peptide`] up by eight.
#[derive(Clone, Debug, Default)]
pub struct ProteinClasses {
    contaminant: fnv::FnvHashSet<Arc<str>>,
    entrapment: fnv::FnvHashSet<Arc<str>>,
}

impl ProteinClasses {
    /// Collect the classed accessions of a FASTA. `None` when every protein
    /// is a sample protein, which is the ordinary search.
    pub fn of_fasta(fasta: &Fasta) -> Option<Self> {
        let mut classes = Self::default();
        for (index, (accession, _)) in fasta.targets.iter().enumerate() {
            match fasta.class(index) {
                ProteinClass::Sample => {}
                ProteinClass::Contaminant => {
                    classes.contaminant.insert(accession.clone());
                }
                ProteinClass::Entrapment => {
                    classes.entrapment.insert(accession.clone());
                }
            }
        }
        (!classes.is_empty()).then_some(classes)
    }

    /// The classes of the proteins a peptide belongs to. Accessions that are
    /// in neither set are sample proteins.
    pub fn classes_of<'a>(
        &self,
        proteins: impl IntoIterator<Item = &'a Arc<str>>,
    ) -> ProteinClassSet {
        let mut set = ProteinClassSet::default();
        for protein in proteins {
            if self.contaminant.contains(protein) {
                set.insert(ProteinClass::Contaminant);
            } else if self.entrapment.contains(protein) {
                set.insert(ProteinClass::Entrapment);
            } else {
                set.insert(ProteinClass::Sample);
            }
        }
        set
    }

    /// Whether any protein of this peptide is a contaminant, which is what
    /// makes the peptide one.
    pub fn is_contaminant<'a>(&self, proteins: impl IntoIterator<Item = &'a Arc<str>>) -> bool {
        proteins
            .into_iter()
            .any(|protein| self.contaminant.contains(protein))
    }

    /// Contaminant proteins the FASTA contributed.
    pub fn contaminant_proteins(&self) -> usize {
        self.contaminant.len()
    }

    pub fn is_empty(&self) -> bool {
        self.contaminant.is_empty() && self.entrapment.is_empty()
    }
}

#[derive(Clone, Debug)]
pub struct Fasta {
    pub targets: Vec<(Arc<str>, ProteinSequence)>,
    /// The class of each protein of [`Self::targets`], at the same index.
    /// Parallel to `targets` so that reading the proteins is unchanged.
    classes: Vec<ProteinClass>,
    decoy_tag: String,
    // Should we ignore decoys in the fasta database
    // and generate them internally?
    generate_decoys: bool,
}

impl Fasta {
    // Parse a string into a fasta database
    pub fn parse<S: Into<String>>(
        contents: String,
        decoy_tag: S,
        generate_decoys: bool,
    ) -> Result<Fasta, FastaError> {
        let decoy_tag = decoy_tag.into();

        let mut targets: Vec<(Arc<str>, ProteinSequence)> = Vec::new();
        let mut last_id: Option<(&str, usize)> = None;
        let mut s = String::new();

        for (line_index, line) in contents.as_str().lines().enumerate() {
            let line_number = line_index + 1;
            let line = line.trim();
            if line.is_empty() {
                continue;
            }
            if let Some(id) = line.strip_prefix('>') {
                if !s.is_empty() {
                    let (last_id, header_line) =
                        last_id.ok_or(FastaError::MissingHeader { line: line_number })?;
                    let acc = accession(last_id, header_line)?;
                    let seq = std::mem::take(&mut s);
                    if !acc.contains(&decoy_tag) || !generate_decoys {
                        targets.push((acc, seq.into()));
                    }
                }
                accession(id, line_number)?;
                last_id = Some((id, line_number));
            } else {
                if last_id.is_none() {
                    return Err(FastaError::MissingHeader { line: line_number });
                }
                // Ambiguous uppercase residues (B, J, X, Z) are kept so the
                // rest of the protein remains searchable. J carries the I/L
                // mass; peptides with B, X or Z are expanded or dropped at
                // digestion, see `database.expand_ambiguous_residues`.
                if let Some(residue) = line.bytes().find(|residue| !residue.is_ascii_uppercase()) {
                    return Err(FastaError::InvalidResidue {
                        line: line_number,
                        residue: residue as char,
                    });
                }
                s.push_str(line);
            }
        }

        if !s.is_empty() {
            let (last_id, header_line) = last_id.ok_or(FastaError::MissingHeader {
                line: contents.lines().count().max(1),
            })?;
            let acc = accession(last_id, header_line)?;
            if !acc.contains(&decoy_tag) || !generate_decoys {
                targets.push((acc, s.into()));
            }
        }

        if targets.is_empty() {
            return Err(FastaError::NoSequences);
        }

        Ok(Fasta {
            classes: vec![ProteinClass::Sample; targets.len()],
            targets,
            decoy_tag,
            generate_decoys,
        })
    }

    /// Add every protein of another FASTA under `class`. The decoy tag and
    /// decoy strategy are the search's, so they stay this file's; proteins
    /// that repeat an accession or a sequence of another class are kept, and
    /// the per-peptide rules of [`ProteinClassSet`] decide what the peptide
    /// is.
    pub fn extend_with_class(&mut self, other: Fasta, class: ProteinClass) {
        self.classes
            .extend(std::iter::repeat_n(class, other.targets.len()));
        self.targets.extend(other.targets);
    }

    /// The class of protein `index` of [`Self::targets`].
    pub fn class(&self, index: usize) -> ProteinClass {
        self.classes.get(index).copied().unwrap_or_default()
    }

    /// How many proteins of `class` the database holds.
    pub fn class_count(&self, class: ProteinClass) -> usize {
        self.classes.iter().filter(|&&other| other == class).count()
    }

    pub fn digest(&self, enzyme: &EnzymeParameters) -> Vec<Digest> {
        self.digest_with_custom_cleavages(enzyme, None)
    }

    pub fn digest_with_custom_cleavages(
        &self,
        enzyme: &EnzymeParameters,
        custom_cleavages: Option<&ValidatedCustomCleavageLibrary>,
    ) -> Vec<Digest> {
        self.digest_where(enzyme, custom_cleavages, |_| true)
    }

    /// Digest every protein, keeping only digests accepted by `keep`. The
    /// rejected digests are dropped per protein, so a caller that selects a
    /// fraction of the digest never holds the rest.
    pub fn digest_where(
        &self,
        enzyme: &EnzymeParameters,
        custom_cleavages: Option<&ValidatedCustomCleavageLibrary>,
        keep: impl Fn(&Digest) -> bool + Sync + Send,
    ) -> Vec<Digest> {
        (0..self.targets.len())
            .into_par_iter()
            .flat_map_iter(|index| {
                let mut digests = self.digest_protein(index, enzyme, custom_cleavages);
                digests.retain(|digest| keep(digest));
                digests
            })
            .collect()
    }

    /// Digest protein `index` of [`Self::targets`]. Proteins carrying the
    /// decoy tag give decoy digests, or none when decoys are generated.
    pub fn digest_protein(
        &self,
        index: usize,
        enzyme: &EnzymeParameters,
        custom_cleavages: Option<&ValidatedCustomCleavageLibrary>,
    ) -> Vec<Digest> {
        let (protein, sequence) = &self.targets[index];
        let is_decoy = protein.contains(&self.decoy_tag);
        if is_decoy && self.generate_decoys {
            return Vec::new();
        }
        let boundaries = custom_cleavages
            .map(|library| library.boundaries_for(protein))
            .unwrap_or_default();
        let mut digests =
            enzyme.digest_protein_with_custom_cleavages(sequence, protein.clone(), boundaries);
        if is_decoy {
            for digest in &mut digests {
                digest.decoy = true;
            }
        }
        digests
    }

    /// Number of proteins containing B, X or Z, residues without a single
    /// mass, whose peptides are expanded or dropped at digestion.
    pub fn ambiguous_protein_count(&self) -> usize {
        self.targets
            .par_iter()
            .filter(|(_, sequence)| ambiguous_residues::is_ambiguous(sequence.as_bytes()))
            .count()
    }

    /// Counts of digests with ambiguous residues that `enzyme` expands or
    /// drops, or `None` when it does not expand them. Only proteins that
    /// contain B, X or Z are digested.
    pub fn ambiguous_expansion_summary(
        &self,
        enzyme: &EnzymeParameters,
        custom_cleavages: Option<&ValidatedCustomCleavageLibrary>,
    ) -> Option<AmbiguousExpansionSummary> {
        let max_variants = enzyme.ambiguous_variants?;
        let as_written = EnzymeParameters {
            ambiguous_variants: None,
            ..enzyme.clone()
        };
        let summary = (0..self.targets.len())
            .into_par_iter()
            .filter(|&index| ambiguous_residues::is_ambiguous(self.targets[index].1.as_bytes()))
            .map(|index| {
                let mut summary = AmbiguousExpansionSummary::default();
                for digest in self.digest_protein(index, &as_written, custom_cleavages) {
                    match ambiguous_residues::variant_count(digest.sequence.as_bytes()) {
                        1 => {}
                        count if count > max_variants => summary.dropped += 1,
                        count => {
                            summary.expanded += 1;
                            summary.variants += count;
                        }
                    }
                }
                summary
            })
            .reduce(AmbiguousExpansionSummary::default, |left, right| {
                AmbiguousExpansionSummary {
                    expanded: left.expanded + right.expanded,
                    variants: left.variants + right.variants,
                    dropped: left.dropped + right.dropped,
                }
            });
        Some(summary)
    }

    pub fn iter_chunks(&self, chunk_size: usize) -> impl Iterator<Item = Self> + '_ {
        self.targets
            .chunks(chunk_size)
            .enumerate()
            .map(move |(index, target_chunk)| {
                let start = index * chunk_size;
                Self {
                    targets: target_chunk.to_vec(),
                    classes: self.classes[start..start + target_chunk.len()].to_vec(),
                    decoy_tag: self.decoy_tag.clone(),
                    generate_decoys: self.generate_decoys,
                }
            })
    }
}

/// See [`Fasta::ambiguous_expansion_summary`]. Counts are per digest, before
/// identical sequences are merged and before generated decoys.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct AmbiguousExpansionSummary {
    /// Digests replaced by their expansions.
    pub expanded: usize,
    /// Expanded digests created in their place.
    pub variants: usize,
    /// Digests with more expansions than allowed, not searched.
    pub dropped: usize,
}

fn accession(id: &str, line: usize) -> Result<Arc<str>, FastaError> {
    id.split_ascii_whitespace()
        .next()
        .filter(|accession| !accession.is_empty())
        .map(Arc::from)
        .ok_or(FastaError::MissingIdentifier { line })
}

#[derive(Debug, PartialEq, Eq)]
pub enum FastaError {
    NoSequences,
    MissingHeader { line: usize },
    MissingIdentifier { line: usize },
    InvalidResidue { line: usize, residue: char },
}

impl std::fmt::Display for FastaError {
    fn fmt(&self, formatter: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::NoSequences => write!(formatter, "FASTA contains no usable protein sequences"),
            Self::MissingHeader { line } => {
                write!(
                    formatter,
                    "FASTA sequence at line {line} appears before its header"
                )
            }
            Self::MissingIdentifier { line } => {
                write!(formatter, "FASTA header at line {line} has no identifier")
            }
            Self::InvalidResidue { line, residue } => write!(
                formatter,
                "FASTA sequence at line {line} contains invalid residue `{residue}`"
            ),
        }
    }
}

impl std::error::Error for FastaError {}

#[cfg(test)]
#[path = "../tests/unit/fasta.rs"]
mod tests;
