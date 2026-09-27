use super::*;
use sage_core::database::Builder;

const FASTA: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/../../tests/Q99536.fasta");
const FASTA_SHA256: &str = "c49e2dac2407dd9f6ad41e8ae1e53e798e41af92264cfcbb84cf8cd6b156e76b";
const MZML: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../../tests/LQSRPAAPPAPGPGQLTLR.mzML"
);
const MZML_SHA256: &str = "b22b4253ab566878b74c0ade3afc4abd1986edf6b62fafb18545377c4327adb2";

fn parameters(generate_decoys: bool) -> Parameters {
    let mut builder: Builder = serde_json::from_value(serde_json::json!({
        "fasta": FASTA,
        "decoy_tag": "rev_",
        "generate_decoys": generate_decoys,
    }))
    .unwrap();
    builder.update_fasta(FASTA.into());
    builder.make_parameters()
}

#[test]
fn fasta_hash_is_stable_and_matches_sha256sum() {
    let parameters = parameters(true);
    let first = describe_fasta(&parameters, &IndexedDatabase::default());
    let second = describe_fasta(&parameters, &IndexedDatabase::default());
    assert_eq!(first, second);
    assert_eq!(first.sha256.as_deref(), Some(FASTA_SHA256));
    assert_eq!(
        first.size_bytes,
        Some(std::fs::metadata(FASTA).unwrap().len())
    );
    assert_eq!(first.proteins, 1);
    assert_eq!(first.target_proteins, 1);
    assert_eq!(first.decoy_proteins, 0);
    assert_eq!(first.decoys.strategy, "generated");
    assert_eq!(
        first.decoys.method.as_deref(),
        Some("reversed_peptide_keep_termini")
    );
    let uniprot = first.uniprot.unwrap();
    assert_eq!((uniprot.reviewed, uniprot.unreviewed), (1, 0));
    assert_eq!(
        uniprot.organisms,
        vec![Organism {
            name: Some("Homo sapiens".into()),
            taxonomy_id: Some(9606),
            proteins: 1,
        }]
    );
}

#[test]
fn supplied_decoys_are_counted_and_named() {
    let fasta = describe_fasta(&parameters(false), &IndexedDatabase::default());
    assert_eq!(fasta.decoys.strategy, "supplied");
    assert_eq!(fasta.decoys.method, None);
    assert_eq!(fasta.decoys.decoy_tag, "rev_");
}

#[test]
fn headers_count_decoys_sources_and_organisms() {
    let mut headers = FastaHeaders::new("rev_");
    for line in [
        ">sp|P1|A_HUMAN Alpha OS=Homo sapiens OX=9606 GN=A PE=1 SV=1",
        "PEPTIDEK",
        ">tr|Q2|B_HUMAN Beta OS=Homo sapiens OX=9606",
        ">sp|P3|C_YEAST Gamma OS=Saccharomyces cerevisiae (strain ATCC 204508 / S288c) OX=559292 GN=C",
        ">rev_sp|P1|A_HUMAN Alpha OS=Homo sapiens OX=9606",
        ">custom_protein no tags",
    ] {
        headers.visit(line);
    }
    let mut fasta = FastaProvenance::default();
    headers.finish(&mut fasta);
    assert_eq!(
        (fasta.proteins, fasta.target_proteins, fasta.decoy_proteins),
        (5, 4, 1)
    );
    let uniprot = fasta.uniprot.unwrap();
    assert_eq!((uniprot.reviewed, uniprot.unreviewed), (2, 1));
    assert_eq!(uniprot.organism_count, 2);
    assert_eq!(uniprot.organisms[0].taxonomy_id, Some(9606));
    assert_eq!(uniprot.organisms[0].proteins, 2);
    assert_eq!(
        uniprot.organisms[1].name.as_deref(),
        Some("Saccharomyces cerevisiae (strain ATCC 204508 / S288c)")
    );
}

#[test]
fn plain_headers_record_no_uniprot_fields() {
    let mut headers = FastaHeaders::new("DECOY_");
    headers.visit(">protein_1 something");
    headers.visit(">DECOY_protein_1");
    let mut fasta = FastaProvenance::default();
    headers.finish(&mut fasta);
    assert_eq!((fasta.proteins, fasta.decoy_proteins), (2, 1));
    assert!(fasta.uniprot.is_none());
}

#[test]
fn spectrum_hashing_is_off_by_default() {
    let output = std::env::temp_dir().join(format!(
        "sage-cli-provenance-{}-{}",
        std::process::id(),
        std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .unwrap()
            .as_nanos()
    ));
    let input: crate::input::Input = serde_json::from_value(serde_json::json!({
        "database": { "fasta": FASTA },
        "precursor_tol": { "ppm": [-10, 10] },
        "fragment_tol": { "ppm": [-10, 10] },
        "output_directory": output.to_string_lossy(),
        "mzml_paths": [MZML],
    }))
    .unwrap();
    let search = input.build().unwrap();
    std::fs::remove_dir_all(output).unwrap();
    assert!(!search.record_input_hashes);

    let url = sage_cloudpath::to_url(MZML).unwrap();
    let off = describe_spectrum_file(&url, search.record_input_hashes);
    assert_eq!(off.sha256, None);
    assert_eq!(
        off.sha256_skipped.as_deref(),
        Some("record_input_hashes is off")
    );
    assert_eq!(off.size_bytes, Some(std::fs::metadata(MZML).unwrap().len()));

    let on = describe_spectrum_file(&url, true);
    assert_eq!(on.sha256.as_deref(), Some(MZML_SHA256));
    assert_eq!(on.sha256_skipped, None);
}

#[test]
fn remote_spectrum_files_are_not_hashed() {
    let url = Url::parse("s3://bucket/run.mzML").unwrap();
    let file = describe_spectrum_file(&url, true);
    assert_eq!(file.sha256, None);
    assert_eq!(
        file.sha256_skipped.as_deref(),
        Some("remote spectrum file; not hashed")
    );
}

#[test]
fn directories_are_not_hashed() {
    let url = Url::from_directory_path(env!("CARGO_MANIFEST_DIR")).unwrap();
    let file = describe_spectrum_file(&url, true);
    assert_eq!(file.sha256, None);
    assert_eq!(
        file.sha256_skipped.as_deref(),
        Some("directory input; not hashed")
    );
}

#[test]
fn header_fields_stop_at_the_next_field() {
    let header = "sp|P1|A Alpha beta OS=Mus musculus OX=10090 GN=Abc PE=1 SV=2";
    assert_eq!(header_field(header, "OS"), Some("Mus musculus"));
    assert_eq!(header_field(header, "OX"), Some("10090"));
    assert_eq!(header_field(header, "SV"), Some("2"));
    assert_eq!(header_field(header, "XX"), None);
}
