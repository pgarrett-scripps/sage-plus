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
