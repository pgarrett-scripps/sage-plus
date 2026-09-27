use super::*;
use crate::enzyme::Digest;

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

fn peptide(sequence: &str, decoy: bool, sites: Vec<ProteinOccurrence>) -> Peptide {
    let mut peptide = Peptide::try_from(Digest {
        sequence: sequence.into(),
        ..Default::default()
    })
    .unwrap();
    peptide.decoy = decoy;
    peptide.protein_sites = sites.into();
    peptide
}

#[test]
fn cleaves_between_follows_the_digest_rule() {
    let enzyme = trypsin();
    assert!(enzyme.cleaves_between(b'K', b'A'));
    assert!(enzyme.cleaves_between(b'R', b'G'));
    assert!(!enzyme.cleaves_between(b'K', b'P'));
    assert!(!enzyme.cleaves_between(b'A', b'K'));

    // Asp-N cuts before D.
    let asp_n = Enzyme::new("D", "", false, false).unwrap();
    assert!(asp_n.cleaves_between(b'A', b'D'));
    assert!(!asp_n.cleaves_between(b'D', b'A'));

    let none = Enzyme::new("$", "", true, false).unwrap();
    assert!(!none.cleaves_between(b'K', b'A'));
}

#[test]
fn termini_and_missed_cleavages_are_classified() {
    let enzyme = trypsin();
    let sequence = b"AKPEPTIDEK";
    // K before P is not a site, so this has no missed cleavage.
    assert_eq!(missed_cleavages(&enzyme, sequence), 0);
    assert_eq!(missed_cleavages(&enzyme, b"AKEPTIDERK"), 2);

    let class = |prev, next| classify_occurrence(&enzyme, sequence, &occurrence(prev, next), true);
    assert_eq!(class(Some(b'K'), Some(b'A')), TerminusClass::Enzymatic);
    assert_eq!(class(None, Some(b'A')), TerminusClass::Enzymatic);
    assert_eq!(class(Some(b'A'), Some(b'A')), TerminusClass::SemiN);
    assert_eq!(class(Some(b'R'), Some(b'P')), TerminusClass::SemiC);
    assert_eq!(class(Some(b'A'), Some(b'P')), TerminusClass::NonEnzymatic);
    // Protein C-terminus counts as enzymatic whatever the last residue.
    assert_eq!(
        classify_occurrence(&enzyme, b"SEPTIDE", &occurrence(Some(b'K'), None), true),
        TerminusClass::Enzymatic
    );
}

#[test]
fn most_enzymatic_occurrence_wins() {
    let enzyme = trypsin();
    let shared = peptide(
        "SEPTIDEK",
        false,
        vec![
            occurrence(Some(b'A'), Some(b'A')),
            occurrence(Some(b'K'), Some(b'A')),
        ],
    );
    assert_eq!(
        classify_peptide(&enzyme, &shared, true),
        Some(TerminusClass::Enzymatic)
    );
    assert_eq!(
        classify_peptide(&enzyme, &peptide("SEPTIDEK", false, vec![]), true),
        None
    );
}

#[test]
fn summary_counts_distinct_sequences_and_subtracts_decoys() {
    let enzyme = trypsin();
    let tryptic = || occurrence(Some(b'K'), Some(b'A'));
    let peptides = [
        peptide("SEPTIDEK", false, vec![tryptic()]),
        // Same sequence again (e.g. another PSM or modified form).
        peptide("SEPTIDEK", false, vec![tryptic()]),
        peptide("SEPKTIDEK", false, vec![tryptic()]),
        peptide("SEPKTIDRK", false, vec![tryptic()]),
        peptide("PEPTIDAK", false, vec![occurrence(Some(b'A'), Some(b'A'))]),
        peptide("SEMIK", false, vec![occurrence(Some(b'A'), Some(b'A'))]),
        peptide("DECOYK", true, vec![occurrence(Some(b'A'), Some(b'A'))]),
    ];
    let summary = summarize(Some(&enzyme), true, &peptides);
    assert_eq!(summary.target_peptides, 5);
    assert_eq!(summary.decoy_peptides, 1);
    assert_eq!(summary.peptides, 4);
    assert_eq!(summary.missed_cleavages_1, 1);
    assert_eq!(summary.missed_cleavages_2_plus, 1);
    // Two semi-N targets minus one semi-N decoy.
    assert_eq!(summary.semi_n, 1);
    assert_eq!(summary.semi_c, 0);
    assert_eq!(summary.missed_cleavages_0, 2);
    assert!((summary.semi_n_pct - 25.0).abs() < 1e-9);
    assert!((summary.missed_cleavage_pct - 50.0).abs() < 1e-9);

    let unspecific = summarize(None, true, &peptides);
    assert_eq!(unspecific.semi_n, 0);
    assert_eq!(unspecific.missed_cleavages_0, 4);
}

#[test]
fn clipped_initiator_methionine_peptides_start_the_protein() {
    let enzyme = trypsin();
    let clipped = ProteinOccurrence {
        start: Some(1),
        ..occurrence(Some(b'M'), Some(b'G'))
    };
    // ASPEPTIDEK after the initiator Met of MASPEPTIDEKG.
    assert_eq!(
        classify_occurrence(&enzyme, b"ASPEPTIDEK", &clipped, true),
        TerminusClass::Enzymatic
    );
    assert_eq!(
        classify_occurrence(&enzyme, b"ASPEPTIDEK", &clipped, false),
        TerminusClass::SemiN
    );
    // MetAP does not remove the Met before a large residue.
    assert_eq!(
        classify_occurrence(&enzyme, b"KSPEPTIDEK", &clipped, true),
        TerminusClass::SemiN
    );
    // Residue 2 after an internal-looking M elsewhere is not the protein start.
    let internal = ProteinOccurrence {
        start: Some(7),
        ..occurrence(Some(b'M'), Some(b'G'))
    };
    assert_eq!(
        classify_occurrence(&enzyme, b"ASPEPTIDEK", &internal, true),
        TerminusClass::SemiN
    );
}

#[test]
fn expanded_ambiguous_peptides_are_classified_as_written() {
    let enzyme = trypsin();
    let at = |protein: &str, start: usize, len: usize| ProteinOccurrence {
        protein: "P1".into(),
        start: Some(start as u32),
        prev_aa: start.checked_sub(1).map(|index| protein.as_bytes()[index]),
        next_aa: protein.as_bytes().get(start + len).copied(),
        source: Some(protein.into()),
        met_clipped: false,
    };
    // DEPXIDEK searched as DEPKIDEK: the X was never a cleavage site.
    let written_x = "AAKDEPXIDEKGG";
    let as_k = peptide("DEPKIDEK", false, vec![at(written_x, 3, 8)]);
    // XEPTIDEK after a K searched as PEPTIDEK: the K-X bond was cleaved.
    let leading_x = "AAKXEPTIDEKGG";
    let as_p = peptide("PEPTIDEK", false, vec![at(leading_x, 3, 8)]);
    assert_eq!(
        classify_peptide(&enzyme, &as_p, true),
        Some(TerminusClass::Enzymatic)
    );
    let summary = summarize(Some(&enzyme), true, [&as_k, &as_p]);
    assert_eq!(summary.missed_cleavages_0, 2);
    assert_eq!(summary.semi_n, 0);

    // A generated decoy carries its target's occurrence and is judged on
    // the reversed written span.
    let decoy = peptide("DEDIKPEK", true, vec![at(written_x, 3, 8)]);
    let summary = summarize(Some(&enzyme), true, [&decoy]);
    assert_eq!(summary.decoy_peptides, 1);
    assert_eq!(summary.semi_n + summary.semi_c + summary.non_enzymatic, 0);
}
