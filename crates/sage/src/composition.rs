//! Amino acid composition of identified peptides, against the composition the
//! search database predicts for the same proteins.
//!
//! The observed profile counts residues in the distinct confident target
//! peptide sequences of one file. The expected profile counts residues in the
//! in-silico peptides of the proteins identified in that file, so it carries
//! the same enzyme, length limits and protein abundance bias as the
//! observations and needs no species table.
//!
//! Both profiles count INTERIOR residues only: the residue at a terminus the
//! enzyme produced is dictated by the enzyme, not by the protein, and counting
//! it would say more about trypsin than about the sample. A terminus the
//! enzyme did not produce keeps its residue, so a non-specific digest counts
//! every residue. Modifications are ignored; a residue is counted once per
//! distinct sequence, however many PSMs or peptidoforms carry it.
//!
//! A residue that goes missing from the observed profile is the point of the
//! comparison: no threshold is applied and no test is run here, because with
//! hundreds of thousands of residues any test is significant in every run.

use crate::digestion::{classify_peptide, TerminusClass};
use crate::enzyme::Enzyme;
use crate::peptide::Peptide;

/// The residues of `sequence` that are not at a terminus the enzyme produced.
///
/// `class` is the peptide's terminus classification, or `None` when it has no
/// recorded protein occurrence or the digest was non-specific; both keep every
/// residue. A sequence shorter than the residues to drop contributes nothing.
pub fn interior_residues(sequence: &[u8], class: Option<TerminusClass>) -> &[u8] {
    let (drop_n, drop_c) = match class {
        None | Some(TerminusClass::NonEnzymatic) => (0, 0),
        Some(TerminusClass::Enzymatic) => (1, 1),
        // The ragged terminus keeps its residue; the enzymatic one does not.
        Some(TerminusClass::SemiN) => (0, 1),
        Some(TerminusClass::SemiC) => (1, 0),
    };
    let start = drop_n.min(sequence.len());
    let end = sequence.len().saturating_sub(drop_c).max(start);
    &sequence[start..end]
}

/// Interior residue counts over a set of distinct peptide sequences.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct ResidueCounts {
    /// Counts of `A` to `Z`, by `residue - b'A'`.
    counts: [u64; 26],
    residues: u64,
    peptides: u64,
    other: u64,
}

impl ResidueCounts {
    /// Count the interior residues of one peptide, classifying its termini
    /// against `enzyme`. `enzyme` is `None` for a non-specific digest, which
    /// keeps every residue.
    pub fn add_peptide(
        &mut self,
        enzyme: Option<&Enzyme>,
        clip_n_term_met: bool,
        peptide: &Peptide,
    ) {
        let class = enzyme.and_then(|enzyme| classify_peptide(enzyme, peptide, clip_n_term_met));
        self.add_sequence(peptide.sequence.as_bytes(), class);
    }

    /// Count the interior residues of one sequence; see [`interior_residues`].
    pub fn add_sequence(&mut self, sequence: &[u8], class: Option<TerminusClass>) {
        self.peptides += 1;
        for &residue in interior_residues(sequence, class) {
            match residue.checked_sub(b'A') {
                Some(index) if residue.is_ascii_uppercase() => {
                    self.counts[index as usize] += 1;
                    self.residues += 1;
                }
                // A residue outside A-Z cannot be a fraction of a profile the
                // other side also has; it is counted so the totals stay honest.
                _ => self.other += 1,
            }
        }
    }

    /// Add another set of counts, for example a peptide shared by two files.
    pub fn merge(&mut self, other: &Self) {
        for (total, add) in self.counts.iter_mut().zip(&other.counts) {
            *total += add;
        }
        self.residues += other.residues;
        self.peptides += other.peptides;
        self.other += other.other;
    }

    /// Distinct sequences counted.
    pub fn peptides(&self) -> u64 {
        self.peptides
    }

    /// Interior residues counted, excluding residues outside `A` to `Z`.
    pub fn residues(&self) -> u64 {
        self.residues
    }

    /// Interior residues outside `A` to `Z`, excluded from every fraction.
    pub fn other_residues(&self) -> u64 {
        self.other
    }

    /// Count of one residue.
    pub fn count(&self, residue: u8) -> u64 {
        residue
            .to_ascii_uppercase()
            .checked_sub(b'A')
            .and_then(|index| self.counts.get(index as usize))
            .copied()
            .unwrap_or(0)
    }

    /// Fraction of the counted interior residues that are `residue`; zero
    /// when nothing was counted.
    pub fn fraction(&self, residue: u8) -> f64 {
        if self.residues == 0 {
            return 0.0;
        }
        self.count(residue) as f64 / self.residues as f64
    }
}

/// One residue's observed and expected share, and the ratio between them.
#[derive(Debug, Clone, PartialEq, serde::Serialize, serde::Deserialize)]
pub struct ResidueComposition {
    pub residue: char,
    pub observed: u64,
    pub observed_fraction: f64,
    pub expected: u64,
    pub expected_fraction: f64,
    /// `log2(observed_fraction / expected_fraction)`, or `None` when either
    /// fraction is zero and the ratio has no finite value.
    pub log2_ratio: Option<f64>,
}

/// Observed against expected interior composition for one file or one run.
#[derive(Debug, Clone, Default, PartialEq, serde::Serialize, serde::Deserialize)]
pub struct CompositionSummary {
    /// Distinct confident target peptide sequences observed.
    pub observed_peptides: u64,
    pub observed_residues: u64,
    /// Distinct in-silico peptide sequences of the identified proteins.
    pub expected_peptides: u64,
    pub expected_residues: u64,
    /// One row per residue present in either profile, in alphabetical order.
    pub residues: Vec<ResidueComposition>,
    /// Residue with the largest `|log2_ratio|`, and that ratio. `None` when no
    /// residue has a finite ratio.
    pub largest_shift: Option<(char, f64)>,
    /// Jensen-Shannon distance between the two profiles: the square root of
    /// the Jensen-Shannon divergence in bits, so 0 means identical profiles
    /// and 1 means profiles with no residue in common. `None` when either
    /// profile counted no residues, where the two are not both distributions
    /// and a distance would read as a shift that was never measured.
    pub jensen_shannon_distance: Option<f64>,
}

impl CompositionSummary {
    /// Compare an observed profile with the expected one.
    pub fn new(observed: &ResidueCounts, expected: &ResidueCounts) -> Self {
        let mut residues = Vec::new();
        let mut largest_shift: Option<(char, f64)> = None;
        let mut divergence = 0.0;
        for letter in b'A'..=b'Z' {
            let (observed_count, expected_count) = (observed.count(letter), expected.count(letter));
            if observed_count == 0 && expected_count == 0 {
                continue;
            }
            let (p, q) = (observed.fraction(letter), expected.fraction(letter));
            let log2_ratio = (p > 0.0 && q > 0.0).then(|| (p / q).log2());
            if let Some(ratio) = log2_ratio {
                if largest_shift.is_none_or(|(_, best)| ratio.abs() > best.abs()) {
                    largest_shift = Some((letter as char, ratio));
                }
            }
            let mean = 0.5 * (p + q);
            if mean > 0.0 {
                // 0 log 0 is 0, so a residue missing from one profile
                // contributes only through the other's term.
                if p > 0.0 {
                    divergence += 0.5 * p * (p / mean).log2();
                }
                if q > 0.0 {
                    divergence += 0.5 * q * (q / mean).log2();
                }
            }
            residues.push(ResidueComposition {
                residue: letter as char,
                observed: observed_count,
                observed_fraction: p,
                expected: expected_count,
                expected_fraction: q,
                log2_ratio,
            });
        }
        Self {
            observed_peptides: observed.peptides(),
            observed_residues: observed.residues(),
            expected_peptides: expected.peptides(),
            expected_residues: expected.residues(),
            residues,
            largest_shift,
            // Rounding can leave the divergence a hair outside [0, 1].
            jensen_shannon_distance: (observed.residues() > 0 && expected.residues() > 0)
                .then(|| divergence.clamp(0.0, 1.0).sqrt()),
        }
    }
}

#[cfg(test)]
#[path = "../tests/unit/composition.rs"]
mod test;
