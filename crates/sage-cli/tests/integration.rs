use sage_core::database::Builder;
use sage_core::mass::Tolerance;
use sage_core::scoring::{ScoreType, Scorer};
use sage_core::spectrum::SpectrumProcessor;
use std::process::Command;
use std::time::{SystemTime, UNIX_EPOCH};

#[test]
fn integration() -> anyhow::Result<()> {
    let mut builder = Builder::default();
    builder.update_fasta("foo".into());

    let fasta = sage_cloudpath::util::read_fasta(
        &sage_cloudpath::to_url("../../tests/Q99536.fasta").expect("valid url"),
        "rev_",
        true,
    )?;
    let database = builder.make_parameters().build(fasta);
    let spectra = sage_cloudpath::util::read_mzml(
        &sage_cloudpath::to_url("../../tests/LQSRPAAPPAPGPGQLTLR.mzML").expect("valid url"),
        0,
        None,
    )?;
    assert_eq!(spectra.len(), 1);

    let sp = SpectrumProcessor::new(100, true, 0.0);
    let processed = sp.process(spectra[0].clone());
    assert!(processed.masses.len() <= 300);

    let scorer = Scorer {
        db: &database,
        precursor_tol: Tolerance::Ppm(-50.0, 50.0),
        fragment_tol: Tolerance::Ppm(-10.0, 10.0),
        min_matched_peaks: 4,
        min_isotope_err: -1,
        max_isotope_err: 3,
        min_precursor_charge: 2,
        max_precursor_charge: 4,
        override_precursor_charge: false,
        max_fragment_charge: Some(1),
        chimera: false,
        report_psms: 1,
        wide_window: false,
        annotate_matches: false,
        mass_shift_ppm: 50.0,
        score_type: ScoreType::SageHyperScore,
        use_bitmap: false,
    };

    let psm = scorer.score(&processed);
    assert_eq!(psm.len(), 1);
    assert_eq!(psm[0].matched_peaks, 21);
    assert!(psm[0].localization.is_none());

    Ok(())
}

#[test]
fn spectral_library_cli_writes_both_formats_and_summary() -> anyhow::Result<()> {
    let workspace = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("../..");
    let nonce = SystemTime::now().duration_since(UNIX_EPOCH)?.as_nanos();
    let output_directory = std::env::temp_dir().join(format!(
        "sage-plus-spectral-library-{}-{nonce}",
        std::process::id()
    ));

    let output = Command::new(env!("CARGO_BIN_EXE_sage"))
        .current_dir(&workspace)
        .arg(workspace.join("tests/config_spectral_library.json"))
        .arg("--output_directory")
        .arg(&output_directory)
        .arg("--disable-telemetry-i-dont-want-to-improve-sage")
        .output()?;
    assert!(
        output.status.success(),
        "sage failed:\nstdout:\n{}\nstderr:\n{}",
        String::from_utf8_lossy(&output.stdout),
        String::from_utf8_lossy(&output.stderr)
    );

    let parquet = output_directory.join("spectral_library.sage.parquet");
    let mzspeclib = output_directory.join("spectral_library.mzspeclib.txt");
    assert!(parquet.metadata()?.len() > 0);
    let text = std::fs::read_to_string(mzspeclib)?;
    assert!(text.contains("<Spectrum=1>"));
    assert!(text.contains("MS:1003270|proforma peptidoform ion notation="));

    let summary: serde_json::Value =
        serde_json::from_slice(&std::fs::read(output_directory.join("run-summary.json"))?)?;
    assert_eq!(summary["schema_version"], 3);
    assert_eq!(summary["spectral_library"]["enabled"], true);
    assert_eq!(summary["spectral_library"]["entries"], 1);
    assert_eq!(summary["spectral_library"]["transitions"], 20);
    assert_eq!(
        summary["spectral_library"]["formats"],
        serde_json::json!(["sage_parquet", "mzspeclib"])
    );

    std::fs::remove_dir_all(output_directory)?;
    Ok(())
}
