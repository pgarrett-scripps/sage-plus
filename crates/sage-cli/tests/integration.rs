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
    assert_eq!(summary["schema_version"], 6);
    assert_eq!(summary["spectral_library"]["enabled"], true);
    assert_eq!(summary["spectral_library"]["entries"], 1);
    assert_eq!(summary["spectral_library"]["transitions"], 20);
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
fn spectral_library_search_runs_without_database_config() -> anyhow::Result<()> {
    let workspace = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("../..");
    let nonce = SystemTime::now().duration_since(UNIX_EPOCH)?.as_nanos();
    let root = std::env::temp_dir().join(format!(
        "sage-plus-library-search-{}-{nonce}",
        std::process::id()
    ));
    let export_directory = root.join("export");
    let search_directory = root.join("search");

    let export = Command::new(env!("CARGO_BIN_EXE_sage"))
        .current_dir(&workspace)
        .arg(workspace.join("tests/config_spectral_library.json"))
        .arg("--output_directory")
        .arg(&export_directory)
        .arg("--disable-telemetry-i-dont-want-to-improve-sage")
        .output()?;
    assert!(
        export.status.success(),
        "library export failed:\n{}",
        String::from_utf8_lossy(&export.stderr)
    );

    let config_path = root.join("library-search.json");
    let events_path = root.join("events.jsonl");
    std::fs::write(
        &config_path,
        serde_json::to_vec_pretty(&serde_json::json!({
            "library_search": {
                "path": export_directory.join("spectral_library.sage.parquet")
            },
            "deisotope": true,
            "annotate_matches": true,
            "quant": {
                "tmt": "Tmt6",
                "tmt_settings": { "level": 2 },
                "lfq": true
            },
            "max_fragment_charge": 1,
            "report_psms": 1,
            "output_filter": { "psm_q_value": 1.0 },
            "precursor_tol": { "ppm": [-50, 50] },
            "fragment_tol": { "ppm": [-10, 10] },
            "mzml_paths": [workspace.join("tests/LQSRPAAPPAPGPGQLTLR.mzML")]
        }))?,
    )?;
    let search = Command::new(env!("CARGO_BIN_EXE_sage"))
        .current_dir(&workspace)
        .arg(&config_path)
        .arg("--output_directory")
        .arg(&search_directory)
        .arg("--events-jsonl")
        .arg(&events_path)
        .arg("--disable-telemetry-i-dont-want-to-improve-sage")
        .output()?;
    assert!(
        search.status.success(),
        "library search failed:\nstdout:\n{}\nstderr:\n{}",
        String::from_utf8_lossy(&search.stdout),
        String::from_utf8_lossy(&search.stderr)
    );
    assert!(search_directory.join("results.sage.parquet").is_file());
    assert!(search_directory
        .join("matched_fragments.sage.parquet")
        .is_file());
    assert!(search_directory.join("lfq.parquet").is_file());
    let summary: serde_json::Value =
        serde_json::from_slice(&std::fs::read(search_directory.join("run-summary.json"))?)?;
    assert_eq!(summary["library_search"]["enabled"], true);
    assert_eq!(summary["library_search"]["target_entries"], 1);
    assert_eq!(summary["library_search"]["decoy_entries"], 1);
    assert_eq!(summary["peptides_in_database"], 0);
    assert_eq!(summary["quantification"]["lfq_enabled"], true);
    assert_eq!(summary["quantification"]["tmt"], "tmt6");
    assert_eq!(summary["quantification"]["tmt_features"], 1);
    assert_eq!(
        summary["models"]["library_rescoring"],
        "spectral_angle_fallback"
    );
    let events = std::fs::read_to_string(events_path)?;
    assert!(events.lines().any(|line| {
        serde_json::from_str::<serde_json::Value>(line)
            .is_ok_and(|event| event["code"] == "library_source_overlap")
    }));
    assert!(events.lines().any(|line| {
        serde_json::from_str::<serde_json::Value>(line)
            .is_ok_and(|event| event["code"] == "library_rescoring_fallback")
    }));

    std::fs::remove_dir_all(root)?;
    Ok(())
}
