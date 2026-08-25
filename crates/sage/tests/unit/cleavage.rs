use super::*;

#[test]
fn parses_deduplicates_and_validates_sites() {
    let library = CustomCleavageLibrary::from_tsv(
        "protein\tposition\tcontext\nP1\t4\tPEPK|TIDE\nP1\t4\tPEPK|TIDE\nP2\t0\t\n",
    )
    .unwrap();
    let fasta = Fasta::parse(
        ">P1 description\nMPEPKTIDER\n>P2\nACDE\n".into(),
        "rev_",
        true,
    )
    .unwrap();
    let validated = library.validate(&fasta).unwrap();

    assert_eq!(validated.total_sites, 2);
    assert_eq!(validated.matched_sites, 2);
    assert_eq!(validated.unmatched_sites, 0);
    assert_eq!(validated.sites_without_context, 1);
    assert_eq!(validated.boundaries_for("P1"), &[5]);
    assert_eq!(validated.boundaries_for("P2"), &[1]);
}

#[test]
fn reports_context_mismatch() {
    let library =
        CustomCleavageLibrary::from_tsv("protein\tposition\tcontext\nP1\t4\tPEPR|TIDE\n").unwrap();
    let fasta = Fasta::parse(">P1\nMPEPKTIDER\n".into(), "rev_", true).unwrap();

    let error = library.validate(&fasta).unwrap_err().to_string();
    assert!(error.contains("does not match"));
}

#[test]
fn rejects_terminal_and_negative_positions() {
    let negative = CustomCleavageLibrary::from_tsv("protein\tposition\nP1\t-1\n")
        .unwrap_err()
        .to_string();
    assert!(negative.contains("zero-based residue index"));

    let library = CustomCleavageLibrary::from_tsv("protein\tposition\nP1\t3\n").unwrap();
    let fasta = Fasta::parse(">P1\nACDE\n".into(), "rev_", true).unwrap();
    let terminal = library.validate(&fasta).unwrap_err().to_string();
    assert!(terminal.contains("not internal"));
}

#[test]
fn allows_unmatched_library_subset_but_not_zero_matches() {
    let library =
        CustomCleavageLibrary::from_tsv("protein\tposition\nP1\t0\nMISSING\t1\n").unwrap();
    let fasta = Fasta::parse(">P1\nACDE\n".into(), "rev_", true).unwrap();
    let validated = library.validate(&fasta).unwrap();
    assert_eq!(validated.matched_sites, 1);
    assert_eq!(validated.unmatched_sites, 1);

    let fasta = Fasta::parse(">OTHER\nACDE\n".into(), "rev_", true).unwrap();
    assert!(library.validate(&fasta).is_err());
}
