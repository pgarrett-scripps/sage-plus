use super::*;
use std::collections::HashSet;

#[test]
fn spans_share_protein_storage_but_compare_by_content() {
    let protein: ProteinSequence = "PEPTIDEPEPTIDE".into();
    let first = protein.peptide(0..7).unwrap();
    let second = protein.peptide(7..14).unwrap();
    let owned: PeptideSequence = "PEPTIDE".into();

    assert!(first.shares_storage_with(&second));
    assert_eq!(first.storage_len(), 14);
    assert_eq!(first, second);
    assert_eq!(first, owned);

    let sequences = HashSet::from([first, second, owned]);
    assert_eq!(sequences.len(), 1);
}

#[test]
fn protein_peptide_rejects_out_of_bounds_ranges() {
    let protein: ProteinSequence = String::from("MKTAYIAK").into();
    assert_eq!(protein.as_str(), "MKTAYIAK");
    assert_eq!(protein.as_bytes().len(), 8);
    assert_eq!(&*protein, "MKTAYIAK");
    assert_eq!(protein.as_ref(), "MKTAYIAK");

    #[allow(clippy::reversed_empty_ranges)]
    let reversed = protein.peptide(5..3);
    assert!(reversed.is_none());
    assert!(protein.peptide(0..9).is_none());
    assert_eq!(protein.peptide(8..8).unwrap(), "");
    assert_eq!(protein.peptide(0..8).unwrap(), "MKTAYIAK");

    let span = protein.peptide(2..6).unwrap();
    assert_eq!(span.as_str(), "TAYI");
    assert_eq!(&*span, b"TAYI");
    let (source, start) = span.source();
    assert_eq!(source, protein);
    assert_eq!(start, 2);
    assert!(span.starts_with("TA"));
    assert!(!span.starts_with("MK"));
}

#[test]
fn reversed_internal_keeps_termini() {
    let protein: ProteinSequence = "XXPEPTIDEKXX".into();
    let span = protein.peptide(2..10).unwrap();
    let rev = span.reversed_internal();
    assert_eq!(rev, "PEDITPEK");
    // The reversed peptide owns its storage instead of mutating the protein.
    assert!(!rev.shares_storage_with(&span));
    assert_eq!(protein.as_str(), "XXPEPTIDEKXX");

    assert_eq!(PeptideSequence::from("ABC").reversed_internal(), "ABC");
    assert_eq!(PeptideSequence::from("AB").reversed_internal(), "AB");
    assert_eq!(PeptideSequence::from("A").reversed_internal(), "A");
    assert_eq!(PeptideSequence::default().reversed_internal(), "");

    assert_eq!(reverse_interior(b"PEPTIDEK"), b"PEDITPEK");
    for short in [&b""[..], b"A", b"AB", b"ABC"] {
        assert_eq!(reverse_interior(short), short);
    }
    assert_eq!(reverse_interior(b"ABCD"), b"ACBD");
}

#[test]
fn conversions_produce_identical_keys() {
    let expected = PeptideSequence::from("PEPTIDE");
    let all = [
        PeptideSequence::from(String::from("PEPTIDE")),
        PeptideSequence::from(b"PEPTIDE".to_vec()),
        PeptideSequence::from(b"PEPTIDE".to_vec().into_boxed_slice()),
        PeptideSequence::from(Arc::<[u8]>::from(&b"PEPTIDE"[..])),
        PeptideSequence::from(&b"PEPTIDE"[..]),
        PeptideSequence::from(b"PEPTIDE"),
    ];
    for p in &all {
        assert_eq!(p, &expected);
        assert_eq!(p.as_ref(), b"PEPTIDE");
        let borrowed: &[u8] = p.borrow();
        assert_eq!(borrowed, b"PEPTIDE");
    }
    let set: HashSet<PeptideSequence> = all.into_iter().collect();
    assert!(set.contains(&b"PEPTIDE"[..]));
    assert!(!set.contains(&b"PEPTIDEK"[..]));
}

#[test]
fn default_is_empty_and_shared() {
    let a = PeptideSequence::default();
    let b = PeptideSequence::default();
    assert!(a.is_empty());
    assert_eq!(a, "");
    assert!(a.shares_storage_with(&b));
}

#[test]
// Owned `String`s exercise the `PartialEq<String>` impls.
#[allow(clippy::cmp_owned)]
fn string_comparisons_in_both_directions() {
    let p = PeptideSequence::from("SAGE");
    assert!(p == "SAGE");
    assert!(p == *"SAGE");
    assert!("SAGE" == p);
    assert!(*"SAGE" == p);
    assert!(p == String::from("SAGE"));
    assert!(String::from("SAGE") == p);
    assert!(p != "SAGA");
    assert!("SAGEK" != p);
    assert!(String::from("sage") != p);
}

#[test]
fn ordering_is_lexicographic_by_residues() {
    let protein: ProteinSequence = "ZZZAAB".into();
    // Offsets must not matter: the later "AA" span sorts before "ZZ".
    let late = protein.peptide(3..5).unwrap();
    let early = protein.peptide(0..2).unwrap();
    assert_eq!(late.cmp(&early), Ordering::Less);
    assert_eq!(early.partial_cmp(&late), Some(Ordering::Greater));
    let mut v = [
        PeptideSequence::from("B"),
        early,
        late,
        PeptideSequence::from("AAB"),
    ];
    v.sort();
    let strs: Vec<String> = v.iter().map(|p| p.to_string()).collect();
    assert_eq!(strs, vec!["AA", "AAB", "B", "ZZ"]);
}

#[test]
fn debug_and_display_formatting() {
    let p = PeptideSequence::from("PEK");
    assert_eq!(format!("{p:?}"), "PeptideSequence(\"PEK\")");
    assert_eq!(format!("{p}"), "PEK");
    let bad = PeptideSequence::from(vec![0xff, b'A']);
    assert_eq!(format!("{bad:?}"), "PeptideSequence([255, 65])");
    assert_eq!(format!("{bad}"), "\u{fffd}A");
}

#[test]
fn peptide_sequence_remains_two_machine_words() {
    assert_eq!(
        std::mem::size_of::<PeptideSequence>(),
        2 * std::mem::size_of::<usize>()
    );
}
