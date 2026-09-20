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
    };

    let psm = scorer.score(&processed);
    assert_eq!(psm.len(), 1);
    assert_eq!(psm[0].matched_peaks, 20);
    assert!(psm[0].localization.is_none());

    Ok(())
}

#[test]
fn empty_spectra_inputs_fail_instead_of_writing_successful_summaries() -> anyhow::Result<()> {
    let workspace = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("../..");
    let nonce = SystemTime::now().duration_since(UNIX_EPOCH)?.as_nanos();
    let root = std::env::temp_dir().join(format!(
        "sage-plus-empty-input-{}-{nonce}",
        std::process::id()
    ));
    std::fs::create_dir_all(&root)?;

    for extension in ["mgf", "mzML"] {
        let input = root.join(format!("empty.{extension}"));
        let output_directory = root.join(format!("output-{extension}"));
        std::fs::File::create(&input)?;
        let output = Command::new(env!("CARGO_BIN_EXE_sage"))
            .current_dir(&workspace)
            .arg(workspace.join("tests/config.json"))
            .arg("--output_directory")
            .arg(&output_directory)
            .arg("--disable-telemetry-i-dont-want-to-improve-sage")
            .arg(&input)
            .output()?;

        assert!(
            !output.status.success(),
            "empty {extension} unexpectedly passed"
        );
        assert!(!output_directory.join("run-summary.json").exists());
        let stderr = String::from_utf8_lossy(&output.stderr);
        assert!(stderr.contains("failed to read spectra file"), "{stderr}");
    }

    std::fs::remove_dir_all(root)?;
    Ok(())
}

#[test]
fn cli_batch_override_wins_over_configuration() -> anyhow::Result<()> {
    let workspace = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("../..");
    let root = std::env::temp_dir().join(format!(
        "sage-cli-batch-{}-{}",
        std::process::id(),
        SystemTime::now().duration_since(UNIX_EPOCH)?.as_nanos()
    ));
    std::fs::create_dir_all(&root)?;
    let mut config: serde_json::Value =
        serde_json::from_slice(&std::fs::read(workspace.join("tests/config.json"))?)?;
    config["batch_size"] = 3.into();
    config["mzml_paths"] = serde_json::json!(vec!["tests/LQSRPAAPPAPGPGQLTLR.mzML"; 3]);
    std::fs::write(root.join("config.json"), serde_json::to_vec(&config)?)?;
    let result = Command::new(env!("CARGO_BIN_EXE_sage"))
        .current_dir(&workspace)
        .arg(root.join("config.json"))
        .arg("--batch-size")
        .arg("1")
        .arg("--output_directory")
        .arg(root.join("output"))
        .arg("--events-jsonl")
        .arg(root.join("events.jsonl"))
        .arg("--disable-telemetry-i-dont-want-to-improve-sage")
        .output()?;
    assert!(
        result.status.success(),
        "{}",
        String::from_utf8_lossy(&result.stderr)
    );
    let events: Vec<serde_json::Value> = std::fs::read_to_string(root.join("events.jsonl"))?
        .lines()
        .map(serde_json::from_str)
        .collect::<Result<_, _>>()?;
    let progress: Vec<_> = events
        .iter()
        .filter(|event| event["event"] == "search_progress")
        .map(|event| event["files_completed"].as_u64().unwrap())
        .collect();
    assert_eq!(progress, vec![1, 2, 3]);
    std::fs::remove_dir_all(root)?;
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
    let config_path = std::env::temp_dir().join(format!(
        "sage-plus-spectral-library-config-{}-{nonce}.json",
        std::process::id()
    ));
    let mut config: serde_json::Value = serde_json::from_slice(&std::fs::read(
        workspace.join("tests/config_spectral_library.json"),
    )?)?;
    config["write_pin"] = true.into();
    config["write_report"] = true.into();
    config["spectral_library"]["strategy"] = "consensus".into();
    std::fs::write(&config_path, serde_json::to_vec_pretty(&config)?)?;

    let output = Command::new(env!("CARGO_BIN_EXE_sage"))
        .current_dir(&workspace)
        .arg(&config_path)
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
    assert!(output_directory.join("results.sage.pin").is_file());
    assert!(output_directory.join("results.sage.report.html").is_file());
    let text = std::fs::read_to_string(mzspeclib)?;
    assert!(text.contains("<Spectrum=1>"));
    assert!(text.contains("MS:1003270|proforma peptidoform ion notation="));
    assert!(text.contains("MS:1003067|consensus spectrum"));

    let summary: serde_json::Value =
        serde_json::from_slice(&std::fs::read(output_directory.join("run-summary.json"))?)?;
    assert_eq!(summary["schema_version"], 9);
    assert_eq!(summary["spectral_library"]["enabled"], true);
    assert_eq!(summary["spectral_library"]["entries"], 1);
    assert_eq!(summary["spectral_library"]["transitions"], 19);
    assert_eq!(summary["spectral_library"]["strategy"], "consensus");
    assert_eq!(
        summary["spectral_library"]["formats"],
        serde_json::json!(["sage_parquet", "mzspeclib"])
    );

    std::fs::remove_dir_all(output_directory)?;
    std::fs::remove_file(config_path)?;
    Ok(())
}

#[test]
fn modification_preview_cli_needs_no_search_inputs() -> anyhow::Result<()> {
    let root = std::env::temp_dir().join(format!("sage-preview-{}", std::process::id()));
    std::fs::create_dir_all(&root)?;
    let config = root.join("config.json");
    std::fs::write(
        &config,
        r#"{"database":{"variable_mods":{"~K":[42.0106]}}}"#,
    )?;
    let output = Command::new(env!("CARGO_BIN_EXE_sage"))
        .arg(&config)
        .args(["--preview-modifications", "KAKAK"])
        .output()?;
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    let result: serde_json::Value = serde_json::from_slice(&output.stdout)?;
    assert_eq!(
        result["rules"][0]["eligible_sites"],
        serde_json::json!([{"position":3}])
    );
    assert_eq!(result["variants"].as_array().unwrap().len(), 2);
    assert!(!root.join("run-summary.json").exists());
    std::fs::write(&config, r#"{"database":{"static_mods":{"KK":42}}}"#)?;
    let output = Command::new(env!("CARGO_BIN_EXE_sage"))
        .arg(&config)
        .args(["--preview-modifications", "KAKAK"])
        .output()?;
    assert!(!output.status.success());
    assert!(String::from_utf8_lossy(&output.stderr).contains("invalid modification key `KK`"));
    std::fs::remove_dir_all(root)?;
    Ok(())
}
