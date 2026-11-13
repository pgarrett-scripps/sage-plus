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

/// A profile whose every peptide is fully cleaved, so weighting has one class
/// and cannot change a fraction.
fn profile(sequences: &[&str]) -> CompositionProfile {
    classed(&[(0, sequences)])
}

/// A profile built from `(missed cleavages, sequences)` groups.
fn classed(groups: &[(usize, &[&str])]) -> CompositionProfile {
    let mut profile = CompositionProfile::default();
    for (missed, sequences) in groups {
        for sequence in *sequences {
            profile.add_sequence(
                sequence.as_bytes(),
                Some(TerminusClass::NonEnzymatic),
                *missed,
            );
        }
    }
    profile
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
    let profile = profile(&["ACDEFG", "ACDE"]);
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
    let observed = profile(&["CCAAAAAA"]);
    let expected = profile(&["CCCCAAAA"]);
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
    let summary = CompositionSummary::new(&profile(&["AAAA"]), &profile(&["AACC"]));
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
    let summary = CompositionSummary::new(&profile(&["AAAA"]), &profile(&["CCCC"]));
    assert_eq!(summary.jensen_shannon_distance, Some(1.0));
    assert_eq!(summary.largest_shift, None);
    assert_eq!(summary.residues.len(), 2);
}

#[test]
fn an_empty_profile_is_reported_rather_than_dividing_by_zero() {
    let summary = CompositionSummary::new(&CompositionProfile::default(), &profile(&["ACDE"]));
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

#[test]
fn the_expected_profile_is_matched_to_the_observed_missed_cleavages() {
    // The database holds fully cleaved peptides with no K and one-missed-cleavage
    // peptides that keep a K inside, in equal numbers. What was identified is
    // fully cleaved except for one peptide in eight residues' worth.
    let database = classed(&[(0, &["AAAA", "AAAA"]), (1, &["AAAK", "AAAK"])]);
    let observed = classed(&[(0, &["AAAA", "AAAA", "AAAA"]), (1, &["AAAK"])]);
    let summary = CompositionSummary::new(&observed, &database);
    let lysine = summary
        .residues
        .iter()
        .find(|residue| residue.residue == 'K')
        .unwrap();
    // Unmatched, the database would expect 2 K in 16 residues, 0.125, against
    // the observed 1 in 16, 0.0625: a log2 ratio of -1 that is only the digest.
    assert_eq!(lysine.observed_fraction, 0.0625);
    assert_eq!(lysine.expected_fraction, 0.0625);
    assert_eq!(lysine.log2_ratio, Some(0.0));
    // The expected count stays the raw count of the matched database peptides.
    assert_eq!(lysine.expected, 2);
    assert_eq!(summary.matched_observed_residues, summary.observed_residues);
    assert_eq!(summary.jensen_shannon_distance, Some(0.0));
}

#[test]
fn matching_still_shows_a_residue_that_is_really_missing() {
    // Same digest on both sides, but the observed peptides lost their cysteine.
    let database = classed(&[(0, &["ACAA"]), (1, &["ACAK"])]);
    let observed = classed(&[(0, &["AAAA"]), (1, &["AAAK"])]);
    let summary = CompositionSummary::new(&observed, &database);
    let carbon = summary
        .residues
        .iter()
        .find(|residue| residue.residue == 'C')
        .unwrap();
    assert_eq!(carbon.observed, 0);
    assert_eq!(carbon.log2_ratio, None);
    let lysine = summary
        .residues
        .iter()
        .find(|residue| residue.residue == 'K')
        .unwrap();
    // K is matched class for class, so it does not move.
    assert_eq!(lysine.log2_ratio, Some(0.0));
}

#[test]
fn a_class_the_database_cannot_supply_is_left_out_of_the_weights() {
    // Nothing in the database has two missed cleavages, so those observed
    // residues carry no weight and are reported as unmatched.
    let database = classed(&[(0, &["AAAA"])]);
    let observed = classed(&[(0, &["AAAA"]), (2, &["AAKK"])]);
    let summary = CompositionSummary::new(&observed, &database);
    assert_eq!(summary.observed_residues, 8);
    assert_eq!(summary.matched_observed_residues, 4);
    // The expectation is class 0's profile: all A, so K has no finite ratio.
    let lysine = summary
        .residues
        .iter()
        .find(|residue| residue.residue == 'K')
        .unwrap();
    assert_eq!(lysine.expected_fraction, 0.0);
    assert_eq!(lysine.log2_ratio, None);
}

#[test]
fn with_no_class_in_common_the_plain_database_profile_is_used() {
    let database = classed(&[(1, &["AACC"])]);
    let observed = classed(&[(0, &["AAAA"])]);
    let summary = CompositionSummary::new(&observed, &database);
    assert_eq!(summary.matched_observed_residues, 0);
    let carbon = summary
        .residues
        .iter()
        .find(|residue| residue.residue == 'C')
        .unwrap();
    // Unweighted: half of the database's residues are C.
    assert_eq!(carbon.expected_fraction, 0.5);
}

#[test]
fn merging_profiles_keeps_the_classes_apart() {
    let mut left = classed(&[(0, &["AAAA"])]);
    left.merge(&classed(&[(2, &["KKKK"])]));
    assert_eq!(left.total().residues(), 8);
    assert_eq!(left.classes(), 3);
    assert_eq!(left.class(0).unwrap().count(b'A'), 4);
    assert_eq!(left.class(1).unwrap().residues(), 0);
    assert_eq!(left.class(2).unwrap().count(b'K'), 4);
}
