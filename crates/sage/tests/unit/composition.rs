use super::*;
use crate::enzyme::{Digest, ProteinOccurrence};

fn trypsin() -> Enzyme {
    Enzyme::new("KR", "P", true, false).unwrap()
}

fn occurrence(prev_aa: Option<u8>, next_aa: Option<u8>) -> ProteinOccurrence {
    ProteinOccurrence {
        protein: "P1".into(),
        start: Some(if prev_aa.is_some() { 5 } else { 0 }),
        prev_aa,
        next_aa,
        source: None,
        met_clipped: false,
    }
}

fn peptide(sequence: &str, sites: Vec<ProteinOccurrence>) -> Peptide {
    let mut peptide = Peptide::try_from(Digest {
        sequence: sequence.into(),
        ..Default::default()
    })
    .unwrap();
    peptide.protein_sites = sites.into();
    peptide
}

fn counts(sequences: &[&str]) -> ResidueCounts {
    let mut counts = ResidueCounts::default();
    for sequence in sequences {
        counts.add_sequence(sequence.as_bytes(), Some(TerminusClass::NonEnzymatic));
    }
    counts
}

#[test]
fn interior_drops_only_the_termini_the_enzyme_made() {
    let sequence = b"ACDEK";
    assert_eq!(
        interior_residues(sequence, Some(TerminusClass::Enzymatic)),
        b"CDE"
    );
    // A ragged terminus is not the enzyme's choice, so it keeps its residue.
    assert_eq!(
        interior_residues(sequence, Some(TerminusClass::SemiN)),
        b"ACDE"
    );
    assert_eq!(
        interior_residues(sequence, Some(TerminusClass::SemiC)),
        b"CDEK"
    );
    assert_eq!(
        interior_residues(sequence, Some(TerminusClass::NonEnzymatic)),
        b"ACDEK"
    );
    // A non-specific digest has no classification and keeps everything.
    assert_eq!(interior_residues(sequence, None), b"ACDEK");
}

#[test]
fn short_sequences_contribute_nothing_instead_of_panicking() {
    for sequence in ["", "A", "AC"] {
        assert!(interior_residues(sequence.as_bytes(), Some(TerminusClass::Enzymatic)).is_empty());
    }
    let mut counted = ResidueCounts::default();
    counted.add_sequence(b"AC", Some(TerminusClass::Enzymatic));
    assert_eq!((counted.peptides(), counted.residues()), (1, 0));
    // An empty profile has no fractions to divide by.
    assert_eq!(counted.fraction(b'A'), 0.0);
}

#[test]
fn a_peptide_is_classified_against_the_enzyme_that_made_it() {
    let enzyme = trypsin();
    // K.ACDEK.A: both termini are tryptic, so both residues are dropped.
    let tryptic = peptide("ACDEK", vec![occurrence(Some(b'K'), Some(b'A'))]);
    let mut counted = ResidueCounts::default();
    counted.add_peptide(Some(&enzyme), false, &tryptic);
    assert_eq!(counted.residues(), 3);
    assert_eq!(counted.count(b'A'), 0);
    assert_eq!(counted.count(b'C'), 1);

    // A.ACDEA.A has no tryptic terminus and keeps all five.
    let ragged = peptide("ACDEA", vec![occurrence(Some(b'A'), Some(b'A'))]);
    let mut counted = ResidueCounts::default();
    counted.add_peptide(Some(&enzyme), false, &ragged);
    assert_eq!(counted.residues(), 5);
    assert_eq!(counted.count(b'A'), 2);

    // K.ACDEK.P is tryptic at the N-terminus only: trypsin does not cut
    // before proline, so the C-terminal K is the peptide's own residue.
    let blocked = peptide("ACDEK", vec![occurrence(Some(b'K'), Some(b'P'))]);
    let mut counted = ResidueCounts::default();
    counted.add_peptide(Some(&enzyme), false, &blocked);
    assert_eq!(counted.residues(), 4);
    assert_eq!(counted.count(b'K'), 1);

    // Without an enzyme, nothing is dropped.
    let mut counted = ResidueCounts::default();
    counted.add_peptide(None, false, &tryptic);
    assert_eq!(counted.residues(), 5);
}

#[test]
fn residues_outside_the_alphabet_are_counted_apart() {
    let mut counted = ResidueCounts::default();
    counted.add_sequence(b"AC*DE", Some(TerminusClass::NonEnzymatic));
    assert_eq!(counted.residues(), 4);
    assert_eq!(counted.other_residues(), 1);
    assert_eq!(counted.count(b'*'), 0);
    // Fractions are taken over the alphabet only.
    assert_eq!(counted.fraction(b'A'), 0.25);
}

#[test]
fn merging_adds_every_field() {
    let mut left = counts(&["AC"]);
    left.merge(&counts(&["CD"]));
    assert_eq!((left.peptides(), left.residues()), (2, 4));
    assert_eq!(left.count(b'C'), 2);
}

#[test]
fn identical_profiles_have_no_shift_and_no_distance() {
    let profile = counts(&["ACDEFG", "ACDE"]);
    let summary = CompositionSummary::new(&profile, &profile);
    assert_eq!(summary.observed_residues, summary.expected_residues);
    assert_eq!(summary.residues.len(), 6);
    for residue in &summary.residues {
        assert_eq!(residue.log2_ratio, Some(0.0));
        assert_eq!(residue.observed_fraction, residue.expected_fraction);
    }
    assert_eq!(summary.largest_shift, Some(('A', 0.0)));
    assert_eq!(summary.jensen_shannon_distance, Some(0.0));
}

#[test]
fn a_halved_residue_reads_as_minus_one() {
    // Observed: 2 of 8 are C. Expected: 4 of 8 are C.
    let observed = counts(&["CCAAAAAA"]);
    let expected = counts(&["CCCCAAAA"]);
    let summary = CompositionSummary::new(&observed, &expected);
    let carbon = summary
        .residues
        .iter()
        .find(|residue| residue.residue == 'C')
        .unwrap();
    assert_eq!(carbon.observed, 2);
    assert_eq!(carbon.expected, 4);
    assert_eq!(carbon.observed_fraction, 0.25);
    assert_eq!(carbon.expected_fraction, 0.5);
    assert_eq!(carbon.log2_ratio, Some(-1.0));
    // C is the residue that moved most; A moved the other way, by less.
    assert_eq!(summary.largest_shift, Some(('C', -1.0)));
    assert!(summary.jensen_shannon_distance.unwrap() > 0.0);
}

#[test]
fn a_residue_missing_from_one_side_has_no_finite_ratio() {
    let summary = CompositionSummary::new(&counts(&["AAAA"]), &counts(&["AACC"]));
    let carbon = summary
        .residues
        .iter()
        .find(|residue| residue.residue == 'C')
        .unwrap();
    assert_eq!((carbon.observed, carbon.expected), (0, 2));
    assert_eq!(carbon.log2_ratio, None);
    // The one finite ratio is A's, so it is the largest shift.
    assert_eq!(summary.largest_shift.map(|(residue, _)| residue), Some('A'));
    assert!(summary.jensen_shannon_distance.unwrap() > 0.0);
}

#[test]
fn profiles_with_no_residue_in_common_are_one_apart() {
    let summary = CompositionSummary::new(&counts(&["AAAA"]), &counts(&["CCCC"]));
    assert_eq!(summary.jensen_shannon_distance, Some(1.0));
    assert_eq!(summary.largest_shift, None);
    assert_eq!(summary.residues.len(), 2);
}

#[test]
fn an_empty_profile_is_reported_rather_than_dividing_by_zero() {
    let summary = CompositionSummary::new(&ResidueCounts::default(), &counts(&["ACDE"]));
    assert_eq!(summary.observed_peptides, 0);
    assert_eq!(summary.observed_residues, 0);
    assert_eq!(summary.expected_residues, 4);
    assert_eq!(summary.largest_shift, None);
    for residue in &summary.residues {
        assert_eq!(residue.observed_fraction, 0.0);
        assert_eq!(residue.log2_ratio, None);
    }
    // With nothing observed there is no distance to report.
    assert_eq!(summary.jensen_shannon_distance, None);
}
