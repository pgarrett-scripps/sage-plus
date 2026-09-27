//! Digestion summary of identified peptides: missed cleavages and termini
//! that the enzyme did not produce.
//!
//! Every rate is counted over distinct peptide sequences (modifications
//! ignored), not PSMs. Decoy peptides passing the same filter estimate how
//! many target peptides in each class are false, and that count is
//! subtracted class by class before rates are formed.

use crate::ambiguous_residues::{is_ambiguous, is_expansion_of};
use crate::cleavage::ValidatedCustomCleavageLibrary;
use crate::enzyme::{Enzyme, ProteinOccurrence, METAP_SECOND_RESIDUES};
use crate::peptide::Peptide;
use std::borrow::Cow;
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
    classify_occurrence_with_custom_cleavages(enzyme, sequence, occurrence, clip_n_term_met, &[])
}

/// [`classify_occurrence`] where a terminus at one of `custom_boundaries`
/// (zero-based offsets into the occurrence's protein, as in
/// [`ValidatedCustomCleavageLibrary::boundaries_for`]) is also enzymatic.
pub fn classify_occurrence_with_custom_cleavages(
    enzyme: &Enzyme,
    sequence: &[u8],
    occurrence: &ProteinOccurrence,
    clip_n_term_met: bool,
    custom_boundaries: &[usize],
) -> TerminusClass {
    let (Some(&first), Some(&last)) = (sequence.first(), sequence.last()) else {
        return TerminusClass::Enzymatic;
    };
    let custom = |offset: usize| custom_boundaries.contains(&offset);
    let start = occurrence.start.map(|start| start as usize);
    let n_enzymatic = at_protein_n_terminus(sequence, occurrence, clip_n_term_met)
        || start.is_some_and(custom)
        || occurrence
            .prev_aa
            .is_some_and(|previous| enzyme.cleaves_between(previous, first));
    let c_enzymatic = at_protein_c_terminus(occurrence)
        || start.is_some_and(|start| custom(start + sequence.len()))
        || occurrence
            .next_aa
            .is_some_and(|next| enzyme.cleaves_between(last, next));
    TerminusClass::from_termini(n_enzymatic, c_enzymatic)
}

/// Classify a peptide at its most enzymatic protein occurrence (the one with
/// the most enzymatic termini, first on ties). `None` when the peptide has
/// no recorded protein occurrence.
///
/// A peptide expanded from ambiguous FASTA residues (B, X or Z) is judged on
/// its residues as written, as the digest cleaved them: an X searched as P
/// after a K is still a cleavage product.
pub fn classify_peptide(
    enzyme: &Enzyme,
    peptide: &Peptide,
    clip_n_term_met: bool,
) -> Option<TerminusClass> {
    classify_cleavages(enzyme, peptide, clip_n_term_met, None).1
}

/// Missed cleavages and termini class of `peptide` at its most enzymatic
/// occurrence, fewest missed cleavages first on ties; see
/// [`classify_peptide`]. Without occurrences, missed cleavages are counted
/// on the sequence and the class is `None`.
fn classify_cleavages(
    enzyme: &Enzyme,
    peptide: &Peptide,
    clip_n_term_met: bool,
    custom_cleavages: Option<&ValidatedCustomCleavageLibrary>,
) -> (usize, Option<TerminusClass>) {
    let rank = |class: TerminusClass| match class {
        TerminusClass::Enzymatic => 0,
        TerminusClass::SemiN | TerminusClass::SemiC => 1,
        TerminusClass::NonEnzymatic => 2,
    };
    peptide
        .protein_sites
        .iter()
        // An expanded peptide TSV row keeps its row as written in an
        // occurrence without coordinates; like a plain row, it has no termini.
        .filter(|occurrence| occurrence.start.is_some() || occurrence.source.is_none())
        .map(|occurrence| {
            let written = as_written(peptide, occurrence);
            let custom_boundaries = custom_cleavages
                .map(|library| library.boundaries_for(&occurrence.protein))
                .unwrap_or_default();
            (
                classify_occurrence_with_custom_cleavages(
                    enzyme,
                    &written,
                    occurrence,
                    clip_n_term_met,
                    custom_boundaries,
                ),
                missed_cleavages(enzyme, &written),
            )
        })
        .min_by_key(|&(class, missed)| (rank(class), missed))
        .map(|(class, missed)| (missed, Some(class)))
        .unwrap_or_else(|| (missed_cleavages(enzyme, peptide.sequence.as_bytes()), None))
}

/// `peptide`'s residues as written in the FASTA at `occurrence`. Only a
/// peptide expanded from ambiguous residues differs from its sequence; a
/// generated decoy's span is reversed like the decoy.
fn as_written<'a>(peptide: &'a Peptide, occurrence: &'a ProteinOccurrence) -> Cow<'a, [u8]> {
    let sequence = peptide.sequence.as_bytes();
    let len = sequence.len();
    let span = occurrence
        .source
        .as_ref()
        .zip(occurrence.start)
        .and_then(|(protein, start)| {
            let start = start as usize;
            protein.as_bytes().get(start..start + len)
        });
    let Some(span) = span.filter(|span| is_ambiguous(span)) else {
        return Cow::Borrowed(sequence);
    };
    if is_expansion_of(span, sequence) {
        return Cow::Borrowed(span);
    }
    let mut reversed = span.to_vec();
    if len > 2 {
        reversed[1..len - 1].reverse();
    }
    match is_expansion_of(&reversed, sequence) {
        true => Cow::Owned(reversed),
        false => Cow::Borrowed(sequence),
    }
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
    summarize_with_custom_cleavages(enzyme, clip_n_term_met, None, peptides)
}

/// [`summarize`] where termini at the search's custom cleavage sites count
/// as enzymatic, as the digest produced them.
pub fn summarize_with_custom_cleavages<'a, I>(
    enzyme: Option<&Enzyme>,
    clip_n_term_met: bool,
    custom_cleavages: Option<&ValidatedCustomCleavageLibrary>,
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
            Some(enzyme) => classify_cleavages(enzyme, peptide, clip_n_term_met, custom_cleavages),
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
