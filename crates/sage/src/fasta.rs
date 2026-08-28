use crate::cleavage::ValidatedCustomCleavageLibrary;
use crate::enzyme::{Digest, EnzymeParameters};
use crate::mass::VALID_AA;
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

        let mut targets = Vec::new();
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
                if let Some(residue) = line.bytes().find(|residue| !VALID_AA.contains(residue)) {
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
        self.targets
            .par_iter()
            .flat_map_iter(|(protein, sequence)| {
                let boundaries = custom_cleavages
                    .map(|library| library.boundaries_for(protein))
                    .unwrap_or_default();
                enzyme
                    .digest_protein_with_custom_cleavages(sequence, protein.clone(), boundaries)
                    .into_iter()
                    .filter_map(|mut digest| {
                        if protein.contains(&self.decoy_tag) {
                            if !self.generate_decoys {
                                digest.decoy = true;
                                Some(digest)
                            } else {
                                None
                            }
                        } else {
                            Some(digest)
                        }
                    })
            })
            .collect()
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
