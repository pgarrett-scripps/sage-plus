//! Expansion of ambiguous FASTA residues into the residues they stand for.
//!
//! B (Asp or Asn), Z (Glu or Gln) and X (any residue) have no single mass.
//! When `database.expand_ambiguous_residues` is on, each enzymatic digest
//! containing them is replaced by one digest per combination of standard
//! residues. J (Ile or Leu) needs no expansion: it carries the I/L mass.

/// The 20 standard residues an X may stand for.
pub const STANDARD_RESIDUES: &[u8; 20] = b"ACDEFGHIKLMNPQRSTVWY";

/// Residues an ambiguous FASTA residue may stand for, or `None` when
/// `residue` is not expanded.
pub const fn alternatives(residue: u8) -> Option<&'static [u8]> {
    match residue {
        b'B' => Some(b"DN"),
        b'Z' => Some(b"EQ"),
        b'X' => Some(STANDARD_RESIDUES),
        _ => None,
    }
}

/// Does `sequence` contain a residue that [`expand`] replaces?
pub fn is_ambiguous(sequence: &[u8]) -> bool {
    sequence
        .iter()
        .any(|&residue| alternatives(residue).is_some())
}

/// Number of sequences [`expand`] yields for `sequence`, saturating at
/// `usize::MAX`. A sequence without ambiguous residues counts as one.
pub fn variant_count(sequence: &[u8]) -> usize {
    sequence
        .iter()
        .filter_map(|&residue| alternatives(residue))
        .fold(1usize, |count, choices| count.saturating_mul(choices.len()))
}

/// Every sequence `sequence` may stand for, with each ambiguous residue
/// replaced by one of its [`alternatives`], in lexicographic order of the
/// choices. Call [`variant_count`] first to bound the result.
pub fn expand(sequence: &[u8]) -> Vec<Vec<u8>> {
    let mut variants = vec![sequence.to_vec()];
    for (position, &residue) in sequence.iter().enumerate() {
        let Some(choices) = alternatives(residue) else {
            continue;
        };
        variants = variants
            .into_iter()
            .flat_map(|variant| {
                choices.iter().map(move |&choice| {
                    let mut variant = variant.clone();
                    variant[position] = choice;
                    variant
                })
            })
            .collect();
    }
    variants
}

/// Could `peptide` have been expanded from `database`? True when both have
/// the same length and every residue is equal or one of the database
/// residue's [`alternatives`].
pub fn is_expansion_of(database: &[u8], peptide: &[u8]) -> bool {
    database.len() == peptide.len()
        && database.iter().zip(peptide).all(|(&written, &observed)| {
            written == observed
                || alternatives(written).is_some_and(|choices| choices.contains(&observed))
        })
}

/// `residue` with I, L and J mapped to one symbol (L), so sequences that
/// differ only in Ile/Leu/(Ile or Leu) compare equal.
pub const fn isoleucine_leucine_canonical(residue: u8) -> u8 {
    match residue {
        b'I' | b'J' => b'L',
        other => other,
    }
}

/// Are `left` and `right` equal once I, L and J are one symbol?
pub fn isoleucine_leucine_eq(left: &[u8], right: &[u8]) -> bool {
    left.len() == right.len()
        && left
            .iter()
            .zip(right)
            .all(|(&a, &b)| isoleucine_leucine_canonical(a) == isoleucine_leucine_canonical(b))
}

#[cfg(test)]
#[path = "../tests/unit/ambiguous_residues.rs"]
mod test;
