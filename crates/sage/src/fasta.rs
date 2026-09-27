use crate::cleavage::ValidatedCustomCleavageLibrary;
use crate::enzyme::{Digest, EnzymeParameters};
use crate::mass::monoisotopic;
use crate::sequence::ProteinSequence;
use rayon::prelude::*;
use std::sync::Arc;

#[derive(Clone, Debug)]
pub struct Fasta {
    pub targets: Vec<(Arc<str>, ProteinSequence)>,
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

        let nonstandard = targets
            .iter()
            .filter(|(_, sequence)| {
                sequence
                    .as_str()
                    .bytes()
                    .any(|residue| monoisotopic(residue) == 0.0)
            })
            .count();
        if nonstandard > 0 {
            log::warn!(
                "{nonstandard} FASTA protein(s) contain residues without a defined mass (e.g. B, X, Z); peptides containing them will not be searched"
            );
        }

        Ok(Fasta {
            targets,
            decoy_tag,
            generate_decoys,
        })
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

    pub fn iter_chunks(&self, chunk_size: usize) -> impl Iterator<Item = Self> + '_ {
        self.targets
            .chunks(chunk_size)
            .map(move |target_chunk| Self {
                targets: target_chunk.to_vec(),
                decoy_tag: self.decoy_tag.clone(),
                generate_decoys: self.generate_decoys,
            })
    }
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
