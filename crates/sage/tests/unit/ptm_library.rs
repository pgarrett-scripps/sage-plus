use super::*;

#[test]
fn parses_tsv_with_extra_columns_and_deduplicates() {
    let contents = concat!(
        "score\tmodification\tresidue\tposition\tprotein\n",
        "0.99\tPhospho\tS\t3\tP12345\n",
        "0.95\tPhospho\tS\t3\tP12345\n",
        "0.90\tOxidation\tm\t7\tP12345\n",
    );
    let library = PtmLibrary::from_tsv(contents).unwrap();
    assert_eq!(library.len(), 2);
    let sites = library.sites_for("P12345");
    assert_eq!(sites[0].position, 2);
    assert_eq!(sites[0].residue, b'S');
    assert_eq!(sites[1].position, 6);
    assert_eq!(sites[1].residue, b'M');
}

#[test]
fn rejects_zero_based_tsv_position() {
    let error =
        PtmLibrary::from_tsv("protein\tposition\tresidue\tmodification\nP12345\t0\tS\tPhospho\n")
            .unwrap_err();
    assert!(error.contains("positions are one-based"));
}

#[test]
fn detects_plain_and_compressed_tsv_paths() {
    assert!(is_tsv_path("sites.TSV"));
    assert!(is_tsv_path("s3://bucket/sites.tsv.gz"));
    assert!(!is_tsv_path("sites.parquet"));
}

#[test]
fn typed_library_retains_distinct_attachments_and_checks_boundaries() {
    use crate::enzyme::Position;
    use crate::peptide::Site;
    let library = PtmLibrary::from_tsv("protein\tposition\tresidue\tmodification\tattachment\nP1\t9\tK\tAcetyl\tresidue\nP1\t9\tK\tAcetyl\tpeptide_n_term\n").unwrap();
    assert_eq!(library.len(), 2);
    assert_eq!(
        Attachment::PeptideNTerm.site(0, 5, Position::Internal),
        Some(Site::Nterm)
    );
    assert_eq!(
        Attachment::PeptideNTerm.site(1, 5, Position::Internal),
        None
    );
    assert_eq!(
        Attachment::ProteinNTerm.site(0, 5, Position::Internal),
        None
    );
    assert_eq!(
        Attachment::ProteinNTerm.site(0, 5, Position::Nterm),
        Some(Site::Nterm)
    );
    assert!(PtmLibrary::from_tsv(
        "protein\tposition\tresidue\tmodification\tattachment\nP1\t9\tK\tAcetyl\twrong\n"
    )
    .is_err());
}

#[test]
fn attachment_names_round_trip_and_reject_unknown_values() {
    for attachment in [
        Attachment::Residue,
        Attachment::PeptideNTerm,
        Attachment::PeptideCTerm,
        Attachment::ProteinNTerm,
        Attachment::ProteinCTerm,
    ] {
        assert_eq!(Attachment::parse(attachment.as_str()), Ok(attachment));
    }
    let error = Attachment::parse("n_term").unwrap_err();
    assert!(error.contains("invalid PTM attachment `n_term`"), "{error}");
    assert_eq!(Attachment::default(), Attachment::Residue);
}

#[test]
fn attachment_sites_respect_peptide_and_protein_boundaries() {
    use crate::enzyme::Position;
    use crate::peptide::Site;
    // Any index past the peptide end never resolves, whatever the attachment.
    assert_eq!(Attachment::Residue.site(5, 5, Position::Full), None);
    assert_eq!(Attachment::PeptideCTerm.site(5, 5, Position::Full), None);
    assert_eq!(
        Attachment::Residue.site(3, 5, Position::Internal),
        Some(Site::Sequence(3))
    );
    assert_eq!(
        Attachment::PeptideCTerm.site(4, 5, Position::Internal),
        Some(Site::Cterm)
    );
    assert_eq!(
        Attachment::PeptideCTerm.site(3, 5, Position::Internal),
        None
    );
    // Protein C-terminal attachment needs the peptide to end the protein.
    assert_eq!(
        Attachment::ProteinCTerm.site(4, 5, Position::Internal),
        None
    );
    assert_eq!(Attachment::ProteinCTerm.site(4, 5, Position::Nterm), None);
    assert_eq!(
        Attachment::ProteinCTerm.site(4, 5, Position::Cterm),
        Some(Site::Cterm)
    );
    assert_eq!(
        Attachment::ProteinCTerm.site(4, 5, Position::Full),
        Some(Site::Cterm)
    );
    assert_eq!(Attachment::ProteinCTerm.site(3, 5, Position::Full), None);
    assert_eq!(
        Attachment::ProteinNTerm.site(0, 5, Position::Full),
        Some(Site::Nterm)
    );
    assert_eq!(Attachment::ProteinNTerm.site(0, 5, Position::Cterm), None);
}

#[test]
fn attachment_from_site_maps_terminal_groups_to_peptide_termini() {
    use crate::peptide::Site;
    assert_eq!(Attachment::from_site(Site::Nterm), Attachment::PeptideNTerm);
    assert_eq!(Attachment::from_site(Site::Cterm), Attachment::PeptideCTerm);
    assert_eq!(
        Attachment::from_site(Site::Sequence(7)),
        Attachment::Residue
    );
}

#[test]
fn tsv_sites_sort_by_position_then_modification_then_attachment() {
    let library = PtmLibrary::from_tsv(concat!(
        "protein\tposition\tresidue\tmodification\tattachment\n",
        "P1\t9\tK\tMethyl\tresidue\n",
        "P1\t9\tK\tAcetyl\tpeptide_n_term\n",
        "P1\t2\tS\tPhospho\tresidue\n",
        "P1\t9\tK\tAcetyl\tresidue\n",
        "P2\t1\tM\tAcetyl\tprotein_n_term\n",
    ))
    .unwrap();
    assert_eq!(library.len(), 5);
    assert!(!library.is_empty());
    assert_eq!(library.iter().count(), 5);
    let order: Vec<_> = library
        .sites_for("P1")
        .iter()
        .map(|site| (site.position, &*site.modification, site.attachment))
        .collect();
    assert_eq!(
        order,
        vec![
            (1, "Phospho", Attachment::Residue),
            (8, "Acetyl", Attachment::Residue),
            (8, "Acetyl", Attachment::PeptideNTerm),
            (8, "Methyl", Attachment::Residue),
        ]
    );
    assert_eq!(
        library.sites_for("P2")[0].attachment,
        Attachment::ProteinNTerm
    );
    assert!(library.sites_for("missing").is_empty());
}

#[test]
fn empty_library_is_empty() {
    let library = PtmLibrary::from_tsv("protein\tposition\tresidue\tmodification\n").unwrap();
    assert!(library.is_empty());
    assert_eq!(library.len(), 0);
    assert!(PtmLibrary::default().is_empty());
}

#[test]
fn tsv_header_with_byte_order_mark_is_accepted() {
    let library = PtmLibrary::from_tsv(
        "\u{feff}protein\tposition\tresidue\tmodification\nP1\t4\tT\tPhospho\n",
    )
    .unwrap();
    assert_eq!(library.sites_for("P1")[0].position, 3);
}

#[test]
fn tsv_errors_name_the_row_and_column() {
    let header = "protein\tposition\tresidue\tmodification\n";
    let error = PtmLibrary::from_tsv("protein\tposition\tresidue\nP1\t1\tS\n").unwrap_err();
    assert_eq!(
        error,
        "PTM library is missing required column `modification`"
    );

    let error = PtmLibrary::from_tsv(&format!("{header}P1\t1\tS\tPhospho\nP1\t\tS\tPhospho\n"))
        .unwrap_err();
    assert_eq!(error, "PTM library row 3 has an empty `position`");

    let error = PtmLibrary::from_tsv(&format!("{header}P1\t-2\tS\tPhospho\n")).unwrap_err();
    assert_eq!(error, "PTM library row 2 has an invalid `position`");

    for residue in ["ST", "1"] {
        let error =
            PtmLibrary::from_tsv(&format!("{header}P1\t2\t{residue}\tPhospho\n")).unwrap_err();
        assert_eq!(
            error,
            "PTM library row 2 has an invalid one-letter `residue`"
        );
    }

    let error = PtmLibrary::from_tsv(&format!("{header}\t2\tS\tPhospho\n")).unwrap_err();
    assert_eq!(error, "PTM library row 2 has an empty `protein`");

    let error = PtmLibrary::from_tsv(
        "protein\tposition\tresidue\tmodification\tattachment\nP1\t2\tS\tPhospho\t\n",
    )
    .unwrap_err();
    assert_eq!(error, "PTM library row 2 has an empty `attachment`");

    let error = PtmLibrary::from_tsv(&format!("{header}P1\t2\tS\n")).unwrap_err();
    assert!(error.starts_with("invalid PTM library row 2"), "{error}");
}
