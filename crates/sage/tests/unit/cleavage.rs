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

fn tsv_err(content: &str) -> String {
    CustomCleavageLibrary::from_tsv(content)
        .unwrap_err()
        .to_string()
}

#[test]
fn tsv_header_errors_name_the_problem() {
    assert_eq!(tsv_err(""), "custom cleavage-site TSV is empty");
    assert_eq!(tsv_err("\n  \n"), "custom cleavage-site TSV is empty");
    assert!(tsv_err("\nposition\tcontext\nP1\t1\tA|B\n")
        .contains("header on line 2 is missing required `protein` column"));
    assert!(tsv_err("protein\tcontext\nP1\tA|B\n").contains("missing required `position` column"));
    assert_eq!(
        tsv_err("protein\tposition\n"),
        "custom cleavage-site TSV contains no data rows"
    );
}

#[test]
fn tsv_row_errors_report_line_numbers() {
    assert!(tsv_err("protein\tposition\n\t3\n").contains("line 2 has an empty `protein`"));
    // A short row has no position column at all.
    assert!(tsv_err("protein\tposition\nP1\n").contains("line 2 has invalid `position` ``"));
    assert!(tsv_err("protein\tposition\nP1\tfour\n").contains("invalid `position` `four`"));
}

#[test]
fn tsv_accepts_reordered_columns_crlf_and_blank_lines() {
    let library = CustomCleavageLibrary::from_tsv(
        "context\tposition\tprotein\r\n\r\nPK|TI\t 4 \t P1 \r\n\t1\tP1\r\n",
    )
    .unwrap();
    let fasta = Fasta::parse(">P1\nMPEPKTIDER\n".into(), "rev_", true).unwrap();
    let v = library.validate(&fasta).unwrap();
    assert_eq!(v.total_sites, 2);
    assert_eq!(v.sites_without_context, 1);
    // Boundaries come out sorted by position.
    assert_eq!(v.boundaries_for("P1"), &[2, 5]);
    assert!(v.boundaries_for("P2").is_empty());
}

#[test]
fn context_must_have_uppercase_residues_on_both_sides_of_one_bar() {
    for bad in ["PEPK", "|TIDE", "PEPK|", "PE|PK|TIDE"] {
        let err = tsv_err(&format!("protein\tposition\tcontext\nP1\t4\t{bad}\n"));
        assert!(
            err.contains(&format!("line 2 has invalid context `{bad}`")),
            "{bad}: {err}"
        );
    }
    let err = tsv_err("protein\tposition\tcontext\nP1\t4\tpepk|TIDE\n");
    assert!(
        err.contains("invalid amino acids in context `pepk|TIDE`"),
        "{err}"
    );
}

#[test]
fn duplicate_sites_merge_but_conflicting_contexts_fail() {
    let err = tsv_err("protein\tposition\tcontext\nP1\t4\tPEPK|TIDE\nP1\t4\tPEPK|TIDA\n");
    assert!(
        err.contains("conflicting contexts for `P1` position 4"),
        "{err}"
    );

    // A context-free row followed by a row with context keeps the context,
    // and a later context-free duplicate does not erase it.
    let library = CustomCleavageLibrary::from_tsv(
        "protein\tposition\tcontext\nP1\t4\t\nP1\t4\tPEPK|TIDE\nP1\t4\t\n",
    )
    .unwrap();
    let v = library
        .validate(&Fasta::parse(">P1\nMPEPKTIDER\n".into(), "rev_", true).unwrap())
        .unwrap();
    assert_eq!(v.total_sites, 1);
    assert_eq!(v.sites_without_context, 0);

    // Proof the context was kept: it is checked against the sequence.
    let err = library
        .validate(&Fasta::parse(">P1\nMPEPRTIDER\n".into(), "rev_", true).unwrap())
        .unwrap_err()
        .to_string();
    assert!(
        err.contains("does not match FASTA protein `P1` at position 4"),
        "{err}"
    );
}

#[test]
fn from_records_trims_and_validates() {
    let library = CustomCleavageLibrary::from_records(vec![
        (" P1 ".to_string(), 4, Some(" PEPK|TIDE ".to_string())),
        ("P1".to_string(), 1, Some("   ".to_string())),
        ("P2".to_string(), 0, None),
    ])
    .unwrap();
    let fasta = Fasta::parse(">P1\nMPEPKTIDER\n>P2\nAC\n".into(), "rev_", true).unwrap();
    let v = library.validate(&fasta).unwrap();
    assert_eq!(v.total_sites, 3);
    assert_eq!(v.matched_sites, 3);
    assert_eq!(v.sites_without_context, 2);
    assert_eq!(v.boundaries_for("P1"), &[2, 5]);
    assert_eq!(v.boundaries_for("P2"), &[1]);

    let err = CustomCleavageLibrary::from_records(vec![(" ".to_string(), 1, None)])
        .unwrap_err()
        .to_string();
    assert_eq!(err, "custom cleavage-site record 1 has an empty `protein`");
    let err = CustomCleavageLibrary::from_records(vec![
        ("P1".to_string(), 1, None),
        ("P1".to_string(), 2, Some("AB".to_string())),
    ])
    .unwrap_err()
    .to_string();
    assert!(
        err.starts_with("custom cleavage-site record 2 has invalid context"),
        "{err}"
    );
    let err = CustomCleavageLibrary::from_records(Vec::new())
        .unwrap_err()
        .to_string();
    assert_eq!(err, "custom cleavage-site input contains no records");
}

#[test]
fn position_overflow_is_an_error_not_a_panic() {
    let library =
        CustomCleavageLibrary::from_records(vec![("P1".to_string(), usize::MAX, None)]).unwrap();
    let fasta = Fasta::parse(">P1\nACDE\n".into(), "rev_", true).unwrap();
    let err = library.validate(&fasta).unwrap_err().to_string();
    assert!(err.contains("overflows for protein `P1`"), "{err}");
}
