use super::*;

fn whole_protein_digest(length: usize) -> EnzymeParameters {
    EnzymeParameters {
        missed_cleavages: 0,
        min_len: length,
        max_len: length,
        enzyme: None,
    }
}

#[test]
fn parses_multiline_sequences_and_accessions() {
    let fasta = Fasta::parse(
        "\n>P1 description here\nPEP\n TIDE \n\n>P2\nACD\nEF\n".into(),
        "rev_",
        true,
    )
    .unwrap();

    assert_eq!(
        fasta
            .targets
            .iter()
            .map(|(accession, sequence)| (accession.as_ref(), sequence.as_str()))
            .collect::<Vec<_>>(),
        vec![("P1", "PEPTIDE"), ("P2", "ACDEF")]
    );
}

#[test]
fn generated_decoy_mode_filters_supplied_decoys() {
    let fasta = Fasta::parse(
        ">P1\nPEPTIDE\n>rev_P1\nEDITPEP\n>P2_rev_suffix\nACDE\n".into(),
        "rev_",
        true,
    )
    .unwrap();

    assert_eq!(fasta.targets.len(), 1);
    assert_eq!(fasta.targets[0].0.as_ref(), "P1");
}

#[test]
fn supplied_decoys_are_retained_and_marked_during_digest() {
    let fasta = Fasta::parse(">P1\nPEPT\n>rev_P1\nTPEP\n".into(), "rev_", false).unwrap();
    let mut digests = fasta.digest(&whole_protein_digest(4));
    digests.sort_by(|left, right| left.protein.cmp(&right.protein));

    assert_eq!(digests.len(), 2);
    assert!(!digests[0].decoy);
    assert_eq!(digests[0].protein.as_ref(), "P1");
    assert!(digests[1].decoy);
    assert_eq!(digests[1].protein.as_ref(), "rev_P1");
}

#[test]
fn empty_records_are_skipped() {
    let fasta = Fasta::parse(
        ">empty\n>populated\nPEPTIDE\n>also_empty\n".into(),
        "rev_",
        true,
    )
    .unwrap();

    assert_eq!(fasta.targets.len(), 1);
    assert_eq!(fasta.targets[0].0.as_ref(), "populated");
}

#[test]
fn chunks_preserve_order_and_configuration() {
    let fasta = Fasta::parse(
        ">P1\nAA\n>P2\nCC\n>P3\nDD\n>P4\nEE\n>P5\nFF\n".into(),
        "rev_",
        false,
    )
    .unwrap();
    let chunks = fasta.iter_chunks(2).collect::<Vec<_>>();

    assert_eq!(
        chunks
            .iter()
            .map(|chunk| chunk.targets.len())
            .collect::<Vec<_>>(),
        vec![2, 2, 1]
    );
    assert!(chunks.iter().all(|chunk| chunk.decoy_tag == "rev_"));
    assert!(chunks.iter().all(|chunk| !chunk.generate_decoys));
    assert_eq!(
        chunks
            .iter()
            .flat_map(|chunk| chunk.targets.iter())
            .map(|(accession, _)| accession.as_ref())
            .collect::<Vec<_>>(),
        vec!["P1", "P2", "P3", "P4", "P5"]
    );
}

#[test]
fn rejects_sequence_before_header() {
    assert_eq!(
        Fasta::parse("PEPTIDE\n>P1\nACDE\n".into(), "rev_", true).unwrap_err(),
        FastaError::MissingHeader { line: 1 }
    );
}

#[test]
fn rejects_files_without_sequences() {
    assert_eq!(
        Fasta::parse(">empty\n>also-empty\n".into(), "rev_", true).unwrap_err(),
        FastaError::NoSequences
    );
}

#[test]
fn rejects_empty_identifier() {
    assert_eq!(
        Fasta::parse(">\nPEPTIDE\n".into(), "rev_", true).unwrap_err(),
        FastaError::MissingIdentifier { line: 1 }
    );
}

#[test]
fn rejects_invalid_residue_with_line_number() {
    assert_eq!(
        Fasta::parse(">P1\nPEPtIDE\n".into(), "rev_", true).unwrap_err(),
        FastaError::InvalidResidue {
            line: 2,
            residue: 't'
        }
    );
}
