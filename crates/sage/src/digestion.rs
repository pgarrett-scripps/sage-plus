//! Digestion summary of identified peptides: missed cleavages and termini
//! that the enzyme did not produce.
//!
//! Every rate is counted over distinct peptide sequences (modifications
//! ignored), not PSMs. Decoy peptides passing the same filter estimate how
//! many target peptides in each class are false, and that count is
//! subtracted class by class before rates are formed.

use crate::enzyme::{Enzyme, ProteinOccurrence, METAP_SECOND_RESIDUES};
use crate::peptide::Peptide;
use std::collections::HashSet;

/// How a peptide's termini sit against its protein.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum TerminusClass {
    /// Both termini are enzymatic or protein termini.
    Enzymatic,
    /// The N-terminus was not produced by the enzyme (ragged N-terminus).
    SemiN,
    /// The C-terminus was not produced by the enzyme (ragged C-terminus).
    SemiC,
    /// Neither terminus was produced by the enzyme.
    NonEnzymatic,
}

impl TerminusClass {
    fn from_termini(n_enzymatic: bool, c_enzymatic: bool) -> Self {
        match (n_enzymatic, c_enzymatic) {
            (true, true) => Self::Enzymatic,
            (false, true) => Self::SemiN,
            (true, false) => Self::SemiC,
            (false, false) => Self::NonEnzymatic,
        }
    }
}

/// Whether an occurrence starts at its protein's N-terminus, so its
/// N-terminus is not a cleavage product.
///
/// A peptide with no preceding residue starts the protein. With initiator
/// methionine clipping, a peptide at residue 2 after an `M` that MetAP
/// removes (see [`crate::enzyme::metap_clips`]) also starts the mature protein.
fn at_protein_n_terminus(
    sequence: &[u8],
    occurrence: &ProteinOccurrence,
    clip_n_term_met: bool,
) -> bool {
    let clipped = || {
        occurrence.start == Some(1)
            && occurrence.prev_aa == Some(b'M')
            && sequence
                .first()
                .is_some_and(|first| METAP_SECOND_RESIDUES.contains(first))
    };
    occurrence.prev_aa.is_none() || (clip_n_term_met && clipped())
}

/// Whether an occurrence ends at its protein's C-terminus.
fn at_protein_c_terminus(occurrence: &ProteinOccurrence) -> bool {
    occurrence.next_aa.is_none()
}

/// Classify the termini of `sequence` at one protein occurrence.
pub fn classify_occurrence(
    enzyme: &Enzyme,
    sequence: &[u8],
    occurrence: &ProteinOccurrence,
    clip_n_term_met: bool,
) -> TerminusClass {
    let (Some(&first), Some(&last)) = (sequence.first(), sequence.last()) else {
        return TerminusClass::Enzymatic;
    };
    let n_enzymatic = at_protein_n_terminus(sequence, occurrence, clip_n_term_met)
        || occurrence
            .prev_aa
            .is_some_and(|previous| enzyme.cleaves_between(previous, first));
    let c_enzymatic = at_protein_c_terminus(occurrence)
        || occurrence
            .next_aa
            .is_some_and(|next| enzyme.cleaves_between(last, next));
    TerminusClass::from_termini(n_enzymatic, c_enzymatic)
}

/// Classify a peptide at its most enzymatic protein occurrence (the one with
/// the most enzymatic termini, first on ties). `None` when the peptide has
/// no recorded protein occurrence.
pub fn classify_peptide(
    enzyme: &Enzyme,
    peptide: &Peptide,
    clip_n_term_met: bool,
) -> Option<TerminusClass> {
    let sequence = peptide.sequence.as_bytes();
    let rank = |class: TerminusClass| match class {
        TerminusClass::Enzymatic => 0,
        TerminusClass::SemiN | TerminusClass::SemiC => 1,
        TerminusClass::NonEnzymatic => 2,
    };
    peptide
        .protein_sites
        .iter()
        .map(|occurrence| classify_occurrence(enzyme, sequence, occurrence, clip_n_term_met))
        .min_by_key(|&class| rank(class))
}

/// Enzymatic sites inside `sequence`, counted from the residues themselves so
/// semi-enzymatic and custom-cleavage digests are counted alike.
pub fn missed_cleavages(enzyme: &Enzyme, sequence: &[u8]) -> usize {
    sequence
        .windows(2)
        .filter(|pair| enzyme.cleaves_between(pair[0], pair[1]))
        .count()
}

/// Distinct-peptide counts of one population (targets or decoys).
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct DigestionCounts {
    pub peptides: usize,
    /// Peptides with 0, 1, and 2 or more missed cleavages.
    pub missed_cleavages: [usize; 3],
    pub semi_n: usize,
    pub semi_c: usize,
    pub non_enzymatic: usize,
}

impl DigestionCounts {
    fn add(&mut self, missed: usize, class: Option<TerminusClass>) {
        self.peptides += 1;
        self.missed_cleavages[missed.min(2)] += 1;
        match class {
            Some(TerminusClass::SemiN) => self.semi_n += 1,
            Some(TerminusClass::SemiC) => self.semi_c += 1,
            Some(TerminusClass::NonEnzymatic) => self.non_enzymatic += 1,
            Some(TerminusClass::Enzymatic) | None => {}
        }
    }
}

/// Decoy-corrected digestion summary of one file or of the whole run.
///
/// Counts are distinct target peptides minus distinct decoy peptides in the
/// same class, floored at zero. Every rate shares the denominator `peptides`.
#[derive(Debug, Clone, Default, PartialEq, serde::Serialize, serde::Deserialize)]
pub struct DigestionSummary {
    pub target_peptides: usize,
    pub decoy_peptides: usize,
    /// Denominator: target minus decoy peptides.
    pub peptides: usize,
    pub missed_cleavages_0: usize,
    pub missed_cleavages_1: usize,
    pub missed_cleavages_2_plus: usize,
    pub semi_n: usize,
    pub semi_c: usize,
    pub non_enzymatic: usize,
    /// Percent of peptides with at least one missed cleavage.
    pub missed_cleavage_pct: f64,
    pub semi_n_pct: f64,
    pub semi_c_pct: f64,
    pub non_enzymatic_pct: f64,
}

impl DigestionSummary {
    pub fn new(targets: &DigestionCounts, decoys: &DigestionCounts) -> Self {
        let net = |target: usize, decoy: usize| target.saturating_sub(decoy);
        let peptides = net(targets.peptides, decoys.peptides);
        let pct = |count: usize| {
            if peptides == 0 {
                0.0
            } else {
                100.0 * count as f64 / peptides as f64
            }
        };
        let [mc0, mc1, mc2] = std::array::from_fn(|index| {
            net(
                targets.missed_cleavages[index],
                decoys.missed_cleavages[index],
            )
        });
        let semi_n = net(targets.semi_n, decoys.semi_n);
        let semi_c = net(targets.semi_c, decoys.semi_c);
        let non_enzymatic = net(targets.non_enzymatic, decoys.non_enzymatic);
        Self {
            target_peptides: targets.peptides,
            decoy_peptides: decoys.peptides,
            peptides,
            missed_cleavages_0: mc0,
            missed_cleavages_1: mc1,
            missed_cleavages_2_plus: mc2,
            semi_n,
            semi_c,
            non_enzymatic,
            missed_cleavage_pct: pct(mc1 + mc2),
            semi_n_pct: pct(semi_n),
            semi_c_pct: pct(semi_c),
            non_enzymatic_pct: pct(non_enzymatic),
        }
    }
}

/// Summarize identified peptides. `peptides` may repeat and may mix targets
/// and decoys; each distinct sequence is counted once per population. With
/// no enzyme (non-specific digestion) every peptide counts as enzymatic with
/// no missed cleavages.
pub fn summarize<'a, I>(
    enzyme: Option<&Enzyme>,
    clip_n_term_met: bool,
    peptides: I,
) -> DigestionSummary
where
    I: IntoIterator<Item = &'a Peptide>,
{
    let mut seen = HashSet::new();
    let mut targets = DigestionCounts::default();
    let mut decoys = DigestionCounts::default();
    for peptide in peptides {
        if !seen.insert((peptide.decoy, peptide.sequence.as_bytes())) {
            continue;
        }
        let (missed, class) = match enzyme {
            Some(enzyme) => (
                missed_cleavages(enzyme, peptide.sequence.as_bytes()),
                classify_peptide(enzyme, peptide, clip_n_term_met),
            ),
            None => (0, Some(TerminusClass::Enzymatic)),
        };
        if peptide.decoy {
            decoys.add(missed, class);
        } else {
            targets.add(missed, class);
        }
    }
    DigestionSummary::new(&targets, &decoys)
}

#[cfg(test)]
#[path = "../tests/unit/digestion.rs"]
mod test;
