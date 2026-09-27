use super::*;

#[test]
fn test_identify_format() {
    assert_eq!(FileFormat::from("foo.mzml"), FileFormat::MzML);
    assert_eq!(FileFormat::from("foo.mzML"), FileFormat::MzML);
    assert_eq!(FileFormat::from("foo.mzMLb"), FileFormat::MzMLb);
    assert_eq!(FileFormat::from("foo.mgf"), FileFormat::MGF);
    assert_eq!(FileFormat::from("foo.mgf.gz"), FileFormat::MGF);
    assert_eq!(FileFormat::from("foo.tdf"), FileFormat::TDF);
    assert_eq!(FileFormat::from("./tomato/foo.d"), FileFormat::TDF);
    assert_eq!(FileFormat::from("./tomato/foo.d/"), FileFormat::TDF);
    assert_eq!(FileFormat::from("foo.raw"), FileFormat::ThermoRaw);
    assert_eq!(FileFormat::from("foo.RAW"), FileFormat::ThermoRaw);
}

#[cfg(not(feature = "mzmlb"))]
#[test]
fn mzmlb_without_feature_returns_build_instructions() {
    let url = Url::parse("file:///tmp/example.mzMLb").unwrap();
    assert!(matches!(
        read_spectra(&url, 0, None, BrukerProcessingConfig::default(), false),
        Err(Error::Unsupported(message)) if message.contains("`mzmlb` feature")
    ));
}

#[test]
fn thermoraw_rejects_signal_to_noise_mode() {
    let url = Url::parse("file:///tmp/example.raw").unwrap();
    assert!(matches!(
        read_thermoraw(&url, 0, Some(2)),
        Err(Error::Unsupported(_))
    ));
}

#[test]
fn unidentified_spectra_format_returns_an_error() {
    let url = Url::parse("file:///tmp/example.unknown").unwrap();
    assert!(matches!(
        read_spectra(&url, 0, None, BrukerProcessingConfig::default(), false),
        Err(Error::Unsupported(message)) if message.contains("determine the spectra format")
    ));
}

#[test]
fn identify_format_edge_cases() {
    assert_eq!(FileFormat::from("FOO.MZML.GZ"), FileFormat::MzML);
    assert_eq!(FileFormat::from("foo.MGF.GZ"), FileFormat::MGF);
    assert_eq!(FileFormat::from("foo.tdf_bin"), FileFormat::TDF);
    assert_eq!(FileFormat::from("s3://bucket/run.mzMLb"), FileFormat::MzMLb);
    // Compressed mzMLb or raw are not recognised.
    assert_eq!(FileFormat::from("foo.raw.gz"), FileFormat::Unidentified);
    assert_eq!(FileFormat::from("foo.mzml.bz2"), FileFormat::Unidentified);
    assert_eq!(FileFormat::from("foo"), FileFormat::Unidentified);
    assert_eq!(FileFormat::from(""), FileFormat::Unidentified);
}

#[test]
fn only_bruker_reads_within_file_in_parallel() {
    for (format, expected) in [
        (FileFormat::MzML, false),
        (FileFormat::MzMLb, false),
        (FileFormat::MGF, false),
        (FileFormat::TDF, true),
        (FileFormat::ThermoRaw, false),
        (FileFormat::Unidentified, false),
    ] {
        assert_eq!(format.within_file_parallel(), expected, "{format:?}");
    }
}

#[test]
fn local_only_formats_reject_remote_urls() {
    let url = Url::parse("s3://bucket/run.raw").unwrap();
    assert!(matches!(
        read_thermoraw(&url, 0, None),
        Err(Error::InvalidUri)
    ));
    let url = Url::parse("s3://bucket/run.d").unwrap();
    assert!(matches!(
        read_tdf(&url, 0, BrukerProcessingConfig::default(), false),
        Err(Error::InvalidUri)
    ));
}

fn fixture(name: &str) -> String {
    std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../../tests")
        .join(name)
        .to_str()
        .unwrap()
        .to_string()
}

#[test]
fn read_spectra_dispatches_mzml_with_file_id() {
    let url = crate::to_url(&fixture("LQSRPAAPPAPGPGQLTLR.mzML")).unwrap();
    let spectra = read_spectra(&url, 7, None, BrukerProcessingConfig::default(), false).unwrap();
    assert_eq!(spectra.len(), 1);
    let s = &spectra[0];
    assert_eq!(s.file_id, 7);
    assert_eq!(s.ms_level, 2);
    assert!(s.id.ends_with("scan=30069"), "{}", s.id);
    assert_eq!(s.mz.len(), 299);
    assert_eq!(s.intensity.len(), 299);
    assert_eq!(s.precursors.len(), 1);
}

#[test]
fn read_spectra_dispatches_mgf_and_reports_parse_errors() {
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("two.mgf");
    std::fs::write(
        &path,
        "BEGIN IONS\nTITLE=a\nPEPMASS=500.25\nCHARGE=2+\nRTINSECONDS=60\n100.5 10\n200.25 20\nEND IONS\n\
         BEGIN IONS\nTITLE=b\nPEPMASS=600.5\nCHARGE=3+\nRTINSECONDS=120\n300.5 5\nEND IONS\n",
    )
    .unwrap();
    let url = Url::from_file_path(&path).unwrap();
    let spectra = read_spectra(&url, 3, None, BrukerProcessingConfig::default(), false).unwrap();
    assert_eq!(spectra.len(), 2);
    assert_eq!(spectra[0].file_id, 3);
    assert_eq!(spectra[0].mz, vec![100.5, 200.25]);
    assert_eq!(spectra[1].mz, vec![300.5]);
    assert!((spectra[1].precursors[0].mz - 600.5).abs() < 1e-3);

    let bad = dir.path().join("bad.mgf");
    std::fs::write(&bad, "BEGIN IONS\nPEPMASS=notanumber\nEND IONS\n").unwrap();
    let url = Url::from_file_path(&bad).unwrap();
    assert!(matches!(read_mgf(&url, 0), Err(Error::MGF(_))));
}

#[test]
fn read_text_bytes_and_json_round_trip() {
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("x.json");
    std::fs::write(&path, br#"{"a": [1, 2, 3]}"#).unwrap();
    let p = path.to_str().unwrap();

    assert_eq!(read_text(p).unwrap(), r#"{"a": [1, 2, 3]}"#);
    assert_eq!(read_bytes(p).unwrap(), br#"{"a": [1, 2, 3]}"#.to_vec());
    let v: std::collections::BTreeMap<String, Vec<u32>> = read_json(p).unwrap();
    assert_eq!(v["a"], vec![1, 2, 3]);

    let wrong: Result<Vec<u32>, _> = read_json(p);
    assert!(wrong.is_err());
    assert!(read_text(dir.path().join("missing.txt").to_str().unwrap()).is_err());
}

#[test]
fn read_fasta_honours_decoy_generation() {
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("db.fasta");
    std::fs::write(&path, ">sp|P1|A\nMKT\nLLK\n>rev_sp|P1|A\nKLLTKM\n").unwrap();
    let url = Url::from_file_path(&path).unwrap();

    let fasta = read_fasta(&url, "rev_", true).unwrap();
    assert_eq!(fasta.targets.len(), 1);
    assert_eq!(fasta.targets[0].0.as_ref(), "sp|P1|A");
    assert_eq!(fasta.targets[0].1.as_str(), "MKTLLK");

    let fasta = read_fasta(&url, "rev_", false).unwrap();
    assert_eq!(fasta.targets.len(), 2);
    assert_eq!(fasta.targets[1].1.as_str(), "KLLTKM");
}
