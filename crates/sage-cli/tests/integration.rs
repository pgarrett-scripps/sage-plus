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
        mass_recalibration: None,
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
    // One test spectrum is far below the PSMs needed for a recommendation.
    assert!(summary["recommended_tolerances"]["skipped"]
        .as_str()
        .unwrap()
        .contains("fewer than the 100 needed"));
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
fn mass_recalibration_reports_per_file_models() -> anyhow::Result<()> {
    let workspace = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("../..");
    let root = std::env::temp_dir().join(format!(
        "sage-plus-mass-recalibration-{}-{}",
        std::process::id(),
        SystemTime::now().duration_since(UNIX_EPOCH)?.as_nanos()
    ));
    std::fs::create_dir_all(&root)?;
    let mut config: serde_json::Value =
        serde_json::from_slice(&std::fs::read(workspace.join("tests/config.json"))?)?;
    config["mass_recalibration"] = "auto".into();
    std::fs::write(root.join("config.json"), serde_json::to_vec(&config)?)?;
    let result = Command::new(env!("CARGO_BIN_EXE_sage"))
        .current_dir(&workspace)
        .arg(root.join("config.json"))
        .arg("--output_directory")
        .arg(root.join("output"))
        .arg("--disable-telemetry-i-dont-want-to-improve-sage")
        .output()?;
    assert!(
        result.status.success(),
        "{}",
        String::from_utf8_lossy(&result.stderr)
    );
    let summary: serde_json::Value =
        serde_json::from_slice(&std::fs::read(root.join("output/run-summary.json"))?)?;
    let recalibration = &summary["models"]["mass_recalibration"];
    assert_eq!(recalibration["mode"], "auto");
    let file = &recalibration["files"][0];
    assert_eq!(file["file_id"], 0);
    // One spectrum cannot support a model; the search falls back to none.
    assert_eq!(file["precursor"]["skipped"], "too_few_psms");
    assert!(file["precursor"]["model"].is_null());
    assert!(summary["psms_at_one_percent_fdr"].as_u64().is_some());
    std::fs::remove_dir_all(root)?;
    Ok(())
}

#[test]
fn modification_preview_cli_needs_no_search_inputs() -> anyhow::Result<()> {
    let root = std::env::temp_dir().join(format!("sage-preview-{}", std::process::id()));
    std::fs::create_dir_all(&root)?;
    let config = root.join("config.json");
    std::fs::write(
        &config,
        r#"{"database":{"variable_mods":{"Acetyl":{"mass":42.0106,"sites":["internal_residue:K"]}}}}"#,
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
    std::fs::write(&config, r#"{"database":{"enzyme":{"cleave_at":"KB"}}}"#)?;
    let output = Command::new(env!("CARGO_BIN_EXE_sage"))
        .arg(&config)
        .args(["--preview-modifications", "KAKAK"])
        .output()?;
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(!output.status.success());
    assert!(!stderr.contains("panicked"), "{stderr}");
    assert!(stderr.contains("unsupported residues `B`"), "{stderr}");
    std::fs::remove_dir_all(root)?;
    Ok(())
}

fn run_sage_with_events(
    workspace: &std::path::Path,
    config: &std::path::Path,
    root: &std::path::Path,
) -> anyhow::Result<Vec<serde_json::Value>> {
    let result = Command::new(env!("CARGO_BIN_EXE_sage"))
        .current_dir(workspace)
        .arg(config)
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
    Ok(std::fs::read_to_string(root.join("events.jsonl"))?
        .lines()
        .map(serde_json::from_str)
        .collect::<Result<_, _>>()?)
}

#[test]
fn duplicate_spectrum_ids_are_annotated_against_their_own_spectrum() -> anyhow::Result<()> {
    let workspace = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("../..");
    let root = std::env::temp_dir().join(format!(
        "sage-cli-duplicate-ids-{}-{}",
        std::process::id(),
        SystemTime::now().duration_since(UNIX_EPOCH)?.as_nanos()
    ));
    std::fs::create_dir_all(&root)?;

    // Repeat the only spectrum, keeping its native ID, as an MGF with
    // repeated TITLE= lines would.
    let mzml = std::fs::read_to_string(workspace.join("tests/LQSRPAAPPAPGPGQLTLR.mzML"))?;
    let start = mzml.find("<spectrum index=\"0\"").unwrap();
    let end = mzml.find("</spectrum>").unwrap() + "</spectrum>".len();
    let duplicate = mzml[start..end].replacen("index=\"0\"", "index=\"1\"", 1);
    let mzml = format!("{}\n{}{}", &mzml[..end], duplicate, &mzml[end..]);
    std::fs::write(root.join("duplicate.mzML"), mzml)?;

    let mut config: serde_json::Value =
        serde_json::from_slice(&std::fs::read(workspace.join("tests/config.json"))?)?;
    config["mzml_paths"] = serde_json::json!([root.join("duplicate.mzML")]);
    std::fs::write(root.join("config.json"), serde_json::to_vec(&config)?)?;

    let events = run_sage_with_events(&workspace, &root.join("config.json"), &root)?;
    let annotation = events
        .iter()
        .find(|event| event["event"] == "fragment_annotation_completed")
        .expect("fragment annotation event");
    // Each copy is scored and annotated once, exactly like the single-copy
    // input (1 PSM, 22 fragments).
    assert_eq!(annotation["psms"], 2);
    assert_eq!(annotation["fragments"], 44);

    std::fs::remove_dir_all(root)?;
    Ok(())
}

/// MS2 scans without a precursor (Thermo RAW scans whose trailer has no
/// plausible m/z) are skipped by the discovery pass, the spectrum index, and
/// the search, instead of aborting with "missing MS1 precursor".
#[test]
fn ms2_spectra_without_precursors_are_not_searched() -> anyhow::Result<()> {
    let workspace = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("../..");
    let root = std::env::temp_dir().join(format!(
        "sage-cli-no-precursor-{}-{}",
        std::process::id(),
        SystemTime::now().duration_since(UNIX_EPOCH)?.as_nanos()
    ));
    std::fs::create_dir_all(&root)?;

    let mzml = std::fs::read_to_string(workspace.join("tests/LQSRPAAPPAPGPGQLTLR.mzML"))?;
    let start = mzml.find("<spectrum index=\"0\"").unwrap();
    let end = mzml.find("</spectrum>").unwrap() + "</spectrum>".len();
    let copy = mzml[start..end]
        .replacen("index=\"0\"", "index=\"1\"", 1)
        .replacen("scan=30069", "scan=30070", 1);
    let precursors = copy.find("<precursorList").unwrap();
    let precursors_end = copy.find("</precursorList>").unwrap() + "</precursorList>".len();
    let copy = format!("{}{}", &copy[..precursors], &copy[precursors_end..]);
    let mzml = format!("{}\n{}{}", &mzml[..end], copy, &mzml[end..]);
    std::fs::write(root.join("no-precursor.mzML"), mzml)?;

    let mut config: serde_json::Value =
        serde_json::from_slice(&std::fs::read(workspace.join("tests/config.json"))?)?;
    config["mzml_paths"] = serde_json::json!([root.join("no-precursor.mzML")]);
    config["mass_recalibration"] = "auto".into();
    config["database"]["prefilter"] = true.into();
    std::fs::write(root.join("config.json"), serde_json::to_vec(&config)?)?;

    let events = run_sage_with_events(&workspace, &root.join("config.json"), &root)?;
    let annotation = events
        .iter()
        .find(|event| event["event"] == "fragment_annotation_completed")
        .expect("fragment annotation event");
    assert_eq!(annotation["psms"], 1);
    let summary: serde_json::Value =
        serde_json::from_slice(&std::fs::read(root.join("output/run-summary.json"))?)?;
    assert_eq!(
        summary["models"]["mass_recalibration"]["files"][0]["discovery_spectra"],
        1
    );

    std::fs::remove_dir_all(root)?;
    Ok(())
}

#[test]
fn post_fdr_reread_does_not_repeat_file_events() -> anyhow::Result<()> {
    let workspace = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("../..");
    let root = std::env::temp_dir().join(format!(
        "sage-cli-reread-events-{}-{}",
        std::process::id(),
        SystemTime::now().duration_since(UNIX_EPOCH)?.as_nanos()
    ));
    std::fs::create_dir_all(&root)?;
    // tests/config.json enables annotate_matches, which rereads the input. The
    // prefilter reads the input once more before the search.
    let mut config: serde_json::Value = serde_json::from_str(&std::fs::read_to_string(
        workspace.join("tests/config.json"),
    )?)?;
    // Two proteins, so the prefilter streams more than one.
    let fasta = root.join("two-proteins.fasta");
    std::fs::write(
        &fasta,
        std::fs::read_to_string(workspace.join("tests/Q99536.fasta"))?
            + "\n>sp|P02768|ALBU_HUMAN\nMKWVTFISLLFLFSSAYSRGVFRRDAHKSEVAHRFKDLGEENFKALVLIAFAQYLQQCPFEDHVK\n",
    )?;
    config["database"]["fasta"] = serde_json::Value::String(fasta.display().to_string());
    config["database"]["prefilter"] = serde_json::Value::Bool(true);
    let prefilter_config = root.join("prefilter.json");
    std::fs::write(&prefilter_config, serde_json::to_string(&config)?)?;
    for (run, config) in [workspace.join("tests/config.json"), prefilter_config]
        .into_iter()
        .enumerate()
    {
        let run_root = root.join(run.to_string());
        std::fs::create_dir_all(&run_root)?;
        let events = run_sage_with_events(&workspace, &config, &run_root)?;
        assert!(events
            .iter()
            .any(|event| event["event"] == "fragment_annotation_completed"));
        for kind in ["file_started", "file_completed", "spectra_processed"] {
            let count = events.iter().filter(|event| event["event"] == kind).count();
            assert_eq!(
                count,
                1,
                "{kind} emitted {count} times for {}",
                config.display()
            );
        }
    }
    std::fs::remove_dir_all(root)?;
    Ok(())
}

/// A PTM site library must give the same PSMs with and without the prefilter.
/// The test spectrum is LQSRPAAPPAPGPGQLTLR. Its serine is replaced by alanine
/// in the FASTA, and a library-only +15.994915 site on that alanine restores
/// the exact residue mass, so the spectrum is identified only through the
/// library placement.
#[test]
fn ptm_library_sites_match_with_and_without_prefilter() -> anyhow::Result<()> {
    let workspace = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("../..");
    let root = std::env::temp_dir().join(format!(
        "sage-cli-ptm-library-prefilter-{}-{}",
        std::process::id(),
        SystemTime::now().duration_since(UNIX_EPOCH)?.as_nanos()
    ));
    std::fs::create_dir_all(&root)?;

    let fasta = std::fs::read_to_string(workspace.join("tests/Q99536.fasta"))?;
    // `lines` also strips the carriage returns of a Windows checkout.
    let mut lines = fasta.lines();
    let header = lines.next().expect("FASTA header");
    let sequence: String = lines.collect();
    // One-based position 66 is the S of LQSR.
    let start = sequence.find("LQSRPAAPPAPGPGQLTLR").expect("test peptide") + 2;
    assert_eq!(start + 1, 66);
    let mutated = format!("{}A{}", &sequence[..start], &sequence[start + 1..]);
    let fasta_path = root.join("proteins.fasta");
    std::fs::write(
        &fasta_path,
        format!(
            "{header}\n{mutated}\n>sp|P02768|ALBU_HUMAN\n\
             MKWVTFISLLFLFSSAYSRGVFRRDAHKSEVAHRFKDLGEENFKALVLIAFAQYLQQCPFEDHVK\n"
        ),
    )?;
    let library_path = root.join("sites.tsv");
    std::fs::write(
        &library_path,
        "protein\tposition\tresidue\tmodification\nsp|Q99536|VAT1_HUMAN\t66\tA\tHydroxyl\n",
    )?;

    let mut config: serde_json::Value =
        serde_json::from_slice(&std::fs::read(workspace.join("tests/config.json"))?)?;
    config["database"]["fasta"] = fasta_path.display().to_string().into();
    config["database"]["variable_mods"] = serde_json::json!({
        "Hydroxyl": {"mass": 15.994915, "sites": ["A"], "site_mode": "library", "max_count": 1}
    });
    config["database"]["ptm_library"] = serde_json::json!({"path": library_path});
    // Keep every peptide with a single matched fragment, so the prefilter
    // cannot drop a true hit.
    config["database"]["prefilter_min_matched_peaks"] = 1.into();
    config["report_psms"] = 5.into();
    config["write_pin"] = true.into();

    let mut results = Vec::new();
    for prefilter in [false, true] {
        let run_root = root.join(format!("prefilter-{prefilter}"));
        std::fs::create_dir_all(&run_root)?;
        let mut config = config.clone();
        config["database"]["prefilter"] = prefilter.into();
        let config_path = run_root.join("config.json");
        std::fs::write(&config_path, serde_json::to_vec(&config)?)?;
        run_sage_with_events(&workspace, &config_path, &run_root)?;

        let pin = std::fs::read_to_string(run_root.join("output/results.sage.pin"))?;
        let mut lines = pin.lines();
        let headers = lines
            .next()
            .expect("pin header")
            .split('\t')
            .collect::<Vec<_>>();
        let column = |name: &str| headers.iter().position(|h| *h == name).expect(name);
        let columns = [
            "SpecId",
            "Label",
            "rank",
            "Peptide",
            "Proteins",
            "ln(hyperscore)",
            "matched_peaks",
        ]
        .map(column);
        let mut psms = lines
            .map(|line| {
                let fields = line.split('\t').collect::<Vec<_>>();
                columns.map(|idx| fields[idx].to_string())
            })
            .collect::<Vec<_>>();
        psms.sort();
        results.push(psms);
    }

    let (without, with) = (&results[0], &results[1]);
    assert!(
        without
            .iter()
            .any(|psm| psm[3] == "LQA[Hydroxyl]RPAAPPAPGPGQLTLR"
                && psm[4].contains("sp|Q99536|VAT1_HUMAN")),
        "no PSM carries the library site: {without:?}"
    );
    assert_eq!(without, with, "prefilter changed the PTM library PSMs");

    std::fs::remove_dir_all(root)?;
    Ok(())
}

/// The synthetic spectrum is ASPEPTIDEAAK, which this FASTA holds only after
/// its initiator methionine: a clipped protein N-terminal peptide. It is found,
/// fully enzymatic, with the default initiator Met clipping and the same with
/// and without the prefilter, and not at all when clipping is off.
#[test]
fn clipped_initiator_methionine_peptides_match_with_and_without_prefilter() -> anyhow::Result<()> {
    let workspace = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("../..");
    let root = std::env::temp_dir().join(format!(
        "sage-cli-met-clip-{}-{}",
        std::process::id(),
        SystemTime::now().duration_since(UNIX_EPOCH)?.as_nanos()
    ));
    std::fs::create_dir_all(&root)?;
    let q99536 = std::fs::read_to_string(workspace.join("tests/Q99536.fasta"))?;
    let fasta_path = root.join("proteins.fasta");
    std::fs::write(
        &fasta_path,
        format!("{q99536}\n>sp|CLIP|CLIPPED\nMASPEPTIDEAAKGGLLR\n"),
    )?;

    let mut config: serde_json::Value = serde_json::from_slice(&std::fs::read(
        workspace.join("tests/synthetic/config.json"),
    )?)?;
    config["database"]["fasta"] = fasta_path.display().to_string().into();
    config["database"]["prefilter_min_matched_peaks"] = 1.into();
    config["write_pin"] = true.into();

    let run = |name: &str, prefilter: bool, clip: Option<bool>| -> anyhow::Result<_> {
        let run_root = root.join(name);
        std::fs::create_dir_all(&run_root)?;
        let mut config = config.clone();
        config["database"]["prefilter"] = prefilter.into();
        if let Some(clip) = clip {
            config["database"]["clip_n_term_met"] = clip.into();
        }
        let config_path = run_root.join("config.json");
        std::fs::write(&config_path, serde_json::to_vec(&config)?)?;
        run_sage_with_events(&workspace, &config_path, &run_root)?;

        let pin = std::fs::read_to_string(run_root.join("output/results.sage.pin"))?;
        let mut lines = pin.lines();
        let headers = lines
            .next()
            .expect("pin header")
            .split('\t')
            .collect::<Vec<_>>();
        let column = |name: &str| headers.iter().position(|h| *h == name).expect(name);
        let columns = [
            "SpecId",
            "Label",
            "Peptide",
            "Proteins",
            "semi_enzymatic",
            "missed_cleavages",
            "ln(hyperscore)",
            "matched_peaks",
        ]
        .map(column);
        let mut psms = lines
            .map(|line| {
                let fields = line.split('\t').collect::<Vec<_>>();
                columns.map(|idx| fields[idx].to_string())
            })
            .collect::<Vec<_>>();
        psms.sort();
        Ok(psms)
    };

    let without = run("prefilter-false", false, None)?;
    let with = run("prefilter-true", true, None)?;
    let clipped = without
        .iter()
        .find(|psm| psm[2] == "ASPEPTIDEAAK")
        .unwrap_or_else(|| panic!("the clipped peptide is not identified: {without:?}"));
    assert_eq!(clipped[1], "1");
    assert_eq!(clipped[3], "sp|CLIP|CLIPPED");
    assert_eq!(clipped[4], "0", "a clipped peptide is not semi-enzymatic");
    assert_eq!(without, with, "prefilter changed the clipped peptide PSMs");

    let unclipped = run("no-clip", false, Some(false))?;
    assert!(
        unclipped.iter().all(|psm| psm[2] != "ASPEPTIDEAAK"),
        "{unclipped:?}"
    );

    std::fs::remove_dir_all(root)?;
    Ok(())
}

/// Expanded ambiguous residues give the same PSMs with and without the
/// streamed prefilter. The X-containing digest is expanded in one protein, and
/// one of its variants is also written plainly in another, so both the
/// per-protein and the shared-sequence prefilter paths see expanded digests.
#[test]
fn expanded_ambiguous_residues_match_with_and_without_prefilter() -> anyhow::Result<()> {
    let workspace = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("../..");
    let root = std::env::temp_dir().join(format!(
        "sage-cli-ambiguous-prefilter-{}-{}",
        std::process::id(),
        SystemTime::now().duration_since(UNIX_EPOCH)?.as_nanos()
    ));
    std::fs::create_dir_all(&root)?;

    let fasta = std::fs::read_to_string(workspace.join("tests/Q99536.fasta"))?;
    let mut lines = fasta.lines();
    let header = lines.next().expect("FASTA header");
    let sequence: String = lines.collect();
    let peptide = "LQSRPAAPPAPGPGQLTLR";
    // Write the T of QLTLR as X.
    let start = sequence.find(peptide).expect("test peptide") + peptide.len() - 3;
    let mutated = format!("{}X{}", &sequence[..start], &sequence[start + 1..]);
    let fasta_path = root.join("proteins.fasta");
    std::fs::write(
        &fasta_path,
        format!("{header}\n{mutated}\n>sp|SHARED|TEST\nMKLQSRPAAPPAPGPGQLTLRGGK\n"),
    )?;

    let mut config: serde_json::Value =
        serde_json::from_slice(&std::fs::read(workspace.join("tests/config.json"))?)?;
    config["database"]["fasta"] = fasta_path.display().to_string().into();
    config["database"]["expand_ambiguous_residues"] = true.into();
    config["database"]["prefilter_min_matched_peaks"] = 1.into();
    config["report_psms"] = 5.into();
    config["write_pin"] = true.into();

    let mut results = Vec::new();
    for prefilter in [false, true] {
        let run_root = root.join(format!("prefilter-{prefilter}"));
        std::fs::create_dir_all(&run_root)?;
        let mut config = config.clone();
        config["database"]["prefilter"] = prefilter.into();
        let config_path = run_root.join("config.json");
        std::fs::write(&config_path, serde_json::to_vec(&config)?)?;
        run_sage_with_events(&workspace, &config_path, &run_root)?;

        let pin = std::fs::read_to_string(run_root.join("output/results.sage.pin"))?;
        let mut lines = pin.lines();
        let headers = lines
            .next()
            .expect("pin header")
            .split('\t')
            .collect::<Vec<_>>();
        let column = |name: &str| headers.iter().position(|h| *h == name).expect(name);
        let columns = [
            "SpecId",
            "Label",
            "rank",
            "Peptide",
            "Proteins",
            "ln(hyperscore)",
            "matched_peaks",
        ]
        .map(column);
        let mut psms = lines
            .map(|line| {
                let fields = line.split('\t').collect::<Vec<_>>();
                columns.map(|idx| fields[idx].to_string())
            })
            .collect::<Vec<_>>();
        psms.sort();
        results.push(psms);
    }

    let (without, with) = (&results[0], &results[1]);
    assert!(
        without.iter().any(|psm| psm[3] == peptide
            && psm[4].contains("sp|Q99536|VAT1_HUMAN")
            && psm[4].contains("sp|SHARED|TEST")),
        "no PSM maps the expanded peptide to both proteins: {without:?}"
    );
    assert_eq!(without, with, "prefilter changed the expanded PSMs");

    std::fs::remove_dir_all(root)?;
    Ok(())
}

/// End-to-end N-glycosylation motif search on a synthetic spectrum. The sequon
/// of the identified peptide is completed by the residue after it, so the
/// search, localization, and reusable library all need protein context.
#[test]
fn motif_site_search_localizes_and_exports_edge_sites() -> anyhow::Result<()> {
    use sage_core::enzyme::Position;
    use sage_core::ion_series::{IonSeries, Kind};
    use sage_core::peptide::Site;

    let nonce = SystemTime::now().duration_since(UNIX_EPOCH)?.as_nanos();
    let root = std::env::temp_dir().join(format!("sage-motif-{}-{nonce}", std::process::id()));
    std::fs::create_dir_all(&root)?;
    let fasta = ">GLYCO\nMRLSPEPTIDENKSGGAWLNGTEDVAPRQVNPTEFLKDEGNLTAYHR\n\
                 >OTHER\nMKAGDTLEWVNKPQYFLSAKGHNESTIMDR\n\
                 >TWIN\nMRLSPEPTLDENKSGGA\n";
    std::fs::write(root.join("proteins.fasta"), fasta)?;
    let database = serde_json::json!({
        "fasta": root.join("proteins.fasta"),
        "enzyme": {"min_len": 5},
        "static_mods": {},
        "variable_mods": {
            "HexNAc": {"mass": 203.079373, "sites": ["motif:N*-{P}-[ST]"], "max_count": 1}
        }
    });

    // Equal target and decoy search spaces for the motif.
    let parameters =
        serde_json::from_value::<sage_core::database::Builder>(database.clone())?.make_parameters();
    let parsed = sage_core::fasta::Fasta::parse(fasta.into(), "rev_", true)?;
    let peptides = parameters.modify_digests(parameters.digest_unmodified(&parsed));
    let placements = |decoy: bool| {
        peptides
            .iter()
            .filter(|peptide| peptide.decoy == decoy)
            .map(|peptide| peptide.modifications.len())
            .sum::<usize>()
    };
    assert!(placements(false) >= 4);
    assert_eq!(placements(false), placements(true));

    // Theoretical b/y spectrum of the edge-site glycopeptide.
    let truth = peptides
        .iter()
        .find(|peptide| !peptide.decoy && peptide.to_string() == "LSPEPTIDEN[HexNAc]K")
        .expect("edge sequon is expanded");
    assert_eq!(truth.position, Position::Internal);
    assert_eq!(
        truth.rule_sites("motif:N*-{P}-[ST]".parse().unwrap()),
        vec![Site::Sequence(9)]
    );
    let mut ions = [Kind::B, Kind::Y]
        .into_iter()
        .flat_map(|kind| IonSeries::new(truth, kind).map(|ion| ion.monoisotopic_mass))
        .collect::<Vec<_>>();
    ions.sort_by(f32::total_cmp);
    let proton = sage_core::mass::PROTON;
    let mut mgf = format!(
        "BEGIN IONS\nTITLE=glyco-1\nRTINSECONDS=60\nPEPMASS={:.6}\nCHARGE=2+\n",
        (truth.monoisotopic + 2.0 * proton) / 2.0
    );
    for (index, mass) in ions.iter().enumerate() {
        mgf.push_str(&format!("{:.6} {}\n", mass + proton, 1000 - index));
    }
    mgf.push_str("END IONS\n");
    std::fs::write(root.join("spectrum.mgf"), mgf)?;

    let config = serde_json::json!({
        "database": database,
        "mzml_paths": [root.join("spectrum.mgf")],
        "deisotope": false,
        "precursor_tol": {"ppm": [-10, 10]},
        "fragment_tol": {"ppm": [-10, 10]},
        "min_matched_peaks": 4,
        "output_filter": {"psm_q_value": 1.0},
        "ptm_localization": {"enabled": true, "psm_q_value": 1.0, "localization_q_value": 1.0}
    });
    let config_path = root.join("config.json");
    std::fs::write(&config_path, serde_json::to_vec_pretty(&config)?)?;
    let output = Command::new(env!("CARGO_BIN_EXE_sage"))
        .arg(&config_path)
        .arg("--output_directory")
        .arg(root.join("output"))
        .arg("--disable-telemetry-i-dont-want-to-improve-sage")
        .output()?;
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    let library = std::fs::read_to_string(root.join("output/results.sage.ptm-library.tsv"))?;
    // N10 of LSPEPTIDENK is residue 12 of GLYCO (one-based).
    assert!(
        library
            .lines()
            .any(|line| line.starts_with("GLYCO\t12\tN\tHexNAc\tresidue")),
        "{library}"
    );
    // TWIN carries the L twin of the merged peptide, displayed with I.
    assert!(
        library
            .lines()
            .any(|line| line.starts_with("TWIN\t12\tN\tHexNAc\tresidue")),
        "{library}"
    );

    // Preview sees the same edge site only when given the flanking residue.
    let preview = |extra: &[&str]| -> anyhow::Result<serde_json::Value> {
        let output = Command::new(env!("CARGO_BIN_EXE_sage"))
            .arg(&config_path)
            .args(["--preview-modifications", "LSPEPTIDENK"])
            .args(extra)
            .output()?;
        assert!(
            output.status.success(),
            "{}",
            String::from_utf8_lossy(&output.stderr)
        );
        Ok(serde_json::from_slice(&output.stdout)?)
    };
    assert_eq!(
        preview(&[])?["rules"][0]["eligible_sites"],
        serde_json::json!([])
    );
    assert_eq!(
        preview(&["--preview-before", "MR", "--preview-after", "SG"])?["rules"][0]
            ["eligible_sites"],
        serde_json::json!([{"position": 10}])
    );
    std::fs::remove_dir_all(root)?;
    Ok(())
}

#[test]
fn quality_control_outputs_are_written() -> anyhow::Result<()> {
    let workspace = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("../..");
    let root = std::env::temp_dir().join(format!(
        "sage-plus-qc-{}-{}",
        std::process::id(),
        SystemTime::now().duration_since(UNIX_EPOCH)?.as_nanos()
    ));
    std::fs::create_dir_all(&root)?;
    let mut config: serde_json::Value =
        serde_json::from_slice(&std::fs::read(workspace.join("tests/config.json"))?)?;
    // A window covering the whole spectrum matches its most intense peak.
    config["diagnostic_ions"] = serde_json::json!([
        {"name": "anything", "mz": 1000.0, "tolerance": {"da": [-1000.0, 1000.0]}},
        {"name": "HexNAc", "mz": 204.0867}
    ]);
    std::fs::write(root.join("config.json"), serde_json::to_vec(&config)?)?;
    let result = Command::new(env!("CARGO_BIN_EXE_sage"))
        .current_dir(&workspace)
        .arg(root.join("config.json"))
        .arg("--output_directory")
        .arg(root.join("output"))
        .arg("--disable-telemetry-i-dont-want-to-improve-sage")
        .args(["--threads", "2"])
        .output()?;
    assert!(
        result.status.success(),
        "{}",
        String::from_utf8_lossy(&result.stderr)
    );
    let stderr = String::from_utf8_lossy(&result.stderr);
    assert!(stderr.contains("digestion: "), "{stderr}");
    assert!(stderr.contains("peak memory (RSS): "), "{stderr}");
    assert_matches_published_schemas(&workspace, &root.join("output"))?;

    let digestion = std::fs::read_to_string(root.join("output/digestion.tsv"))?;
    let lines = digestion.lines().collect::<Vec<_>>();
    assert_eq!(lines.len(), 3, "{digestion}");
    assert!(lines[0].starts_with("file\ttarget_peptides\tdecoy_peptides\tpeptides\t"));
    assert!(lines[1].starts_with("LQSRPAAPPAPGPGQLTLR.mzML\t"));
    assert!(lines[2].starts_with("total\t"));

    let summary: serde_json::Value =
        serde_json::from_slice(&std::fs::read(root.join("output/run-summary.json"))?)?;
    assert_eq!(summary["schema_version"], 9);
    assert_eq!(summary["execution"]["rayon_threads"], 2);
    if cfg!(unix) {
        assert!(summary["peak_rss_bytes"].as_u64().unwrap() > 0);
    }
    let total = &summary["qc"]["digestion"]["total"];
    assert!(total["target_peptides"].as_u64().is_some());
    assert_eq!(
        summary["qc"]["digestion"]["files"][0]["file"],
        "LQSRPAAPPAPGPGQLTLR.mzML"
    );
    assert!(summary["output_paths"]
        .as_array()
        .unwrap()
        .iter()
        .any(|path| path.as_str().unwrap().ends_with("digestion.tsv")));
    // The test file has no MS1 spectra, so no polymer rows are reported.
    assert_eq!(summary["qc"]["polymers"], serde_json::json!([]));

    assert!(
        stderr.contains("diagnostic ions in 1 MS2 spectra: anything 100.00%, HexNAc 0.00%"),
        "{stderr}"
    );
    let diagnostic = std::fs::read_to_string(root.join("output/diagnostic_ions.tsv"))?;
    let lines = diagnostic.lines().collect::<Vec<_>>();
    assert_eq!(lines.len(), 2, "{diagnostic}");
    assert_eq!(lines[0], "file\tscannr\tion\tmz\trelative_intensity");
    let row = lines[1].split('\t').collect::<Vec<_>>();
    assert_eq!(row[0], "LQSRPAAPPAPGPGQLTLR.mzML");
    assert_eq!(row[2], "anything");
    let relative_intensity = row[4].parse::<f64>()?;
    assert!(relative_intensity > 0.0 && relative_intensity <= 1.0);
    let ions = &summary["qc"]["diagnostic_ions"];
    assert_eq!(ions["ms2_spectra"], 1);
    assert_eq!(ions["ions"][0]["spectra"], 1);
    assert_eq!(ions["ions"][1]["spectra"], 0);
    assert!(summary["output_paths"]
        .as_array()
        .unwrap()
        .iter()
        .any(|path| path.as_str().unwrap().ends_with("diagnostic_ions.tsv")));
    std::fs::remove_dir_all(root)?;
    Ok(())
}

/// With `immonium` absent or `false` the outputs are byte-identical and the
/// parquet immonium columns are null. On with `rescore: false` only fills
/// those columns and adds the PIN columns; scores and q-values are unchanged.
#[test]
fn immonium_off_is_identical_and_on_adds_outputs() -> anyhow::Result<()> {
    let workspace = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("../..");
    let root = std::env::temp_dir().join(format!(
        "sage-plus-immonium-{}-{}",
        std::process::id(),
        SystemTime::now().duration_since(UNIX_EPOCH)?.as_nanos()
    ));
    let mut config: serde_json::Value =
        serde_json::from_slice(&std::fs::read(workspace.join("tests/config.json"))?)?;
    config["write_pin"] = true.into();
    config["report_psms"] = 5.into();

    let run =
        |name: &str, immonium: Option<serde_json::Value>| -> anyhow::Result<std::path::PathBuf> {
            let run_root = root.join(name);
            std::fs::create_dir_all(&run_root)?;
            let mut config = config.clone();
            if let Some(value) = immonium {
                config["immonium"] = value;
            }
            let path = run_root.join("config.json");
            std::fs::write(&path, serde_json::to_vec(&config)?)?;
            run_sage_with_events(&workspace, &path, &run_root)?;
            Ok(run_root.join("output"))
        };
    let base = run("absent", None)?;
    let off = run("off", Some(false.into()))?;
    let report = run("report", Some(serde_json::json!({"rescore": false})))?;
    let on = run("on", Some(true.into()))?;

    for file in [
        "results.sage.pin",
        "results.sage.parquet",
        "matched_fragments.sage.parquet",
    ] {
        assert_eq!(
            std::fs::read(base.join(file))?,
            std::fs::read(off.join(file))?,
            "{file} differs with immonium off"
        );
    }
    for dir in [&base, &off, &report, &on] {
        assert!(!dir.join("immonium.tsv").exists());
    }
    for dir in [&base, &off] {
        let results: serde_json::Value =
            serde_json::from_slice(&std::fs::read(dir.join("results.json"))?)?;
        assert!(results.get("immonium").is_none());
    }

    // Parquet rows split into the immonium columns and everything else.
    let rows = |dir: &std::path::Path| -> anyhow::Result<(Vec<String>, Vec<String>)> {
        use parquet::file::reader::{FileReader, SerializedFileReader};
        let reader =
            SerializedFileReader::new(std::fs::File::open(dir.join("results.sage.parquet"))?)?;
        let (mut other, mut immonium) = (Vec::new(), Vec::new());
        for row in reader.get_row_iter(None)? {
            let row = row?;
            let (mut rest, mut evidence) = (Vec::new(), Vec::new());
            for (name, field) in row.get_column_iter() {
                if name.starts_with("immonium_") {
                    evidence.push(format!("{name}={field}"));
                } else {
                    rest.push(format!("{name}={field}"));
                }
            }
            assert_eq!(evidence.len(), 7, "{evidence:?}");
            other.push(rest.join("|"));
            immonium.push(evidence.join("|"));
        }
        Ok((other, immonium))
    };
    let (base_rows, base_immonium) = rows(&base)?;
    assert!(!base_rows.is_empty());
    assert!(base_immonium
        .iter()
        .all(|row| row.split('|').all(|value| value.ends_with("=null"))));
    let (report_rows, report_immonium) = rows(&report)?;
    assert_eq!(
        base_rows, report_rows,
        "report-only immonium changed results"
    );
    assert!(report_immonium
        .iter()
        .all(|row| row.split('|').all(|value| !value.ends_with("=null"))));

    // Report only: the PIN gains the immonium columns before `Peptide` and is
    // otherwise unchanged.
    let columns = sage_immonium_pin_columns();
    let strip = |pin: String| {
        let rows = pin
            .lines()
            .map(|line| line.split('\t').map(str::to_owned).collect::<Vec<_>>())
            .collect::<Vec<_>>();
        let drop = rows[0]
            .iter()
            .enumerate()
            .filter(|(_, name)| columns.contains(&name.as_str()))
            .map(|(index, _)| index)
            .collect::<Vec<_>>();
        let kept = rows
            .into_iter()
            .map(|row| {
                row.into_iter()
                    .enumerate()
                    .filter(|(index, _)| !drop.contains(index))
                    .map(|(_, value)| value)
                    .collect::<Vec<_>>()
            })
            .collect::<Vec<_>>();
        (drop.len(), kept)
    };
    let (none, base_pin) = strip(std::fs::read_to_string(base.join("results.sage.pin"))?);
    assert_eq!(none, 0);
    let report_pin = std::fs::read_to_string(report.join("results.sage.pin"))?;
    let header = report_pin
        .lines()
        .next()
        .unwrap()
        .split('\t')
        .collect::<Vec<_>>();
    let peptide = header.iter().position(|name| *name == "Peptide").unwrap();
    assert_eq!(header[peptide - columns.len()..peptide], columns);
    let (added, report_pin) = strip(report_pin);
    assert_eq!(added, columns.len());
    assert_eq!(base_pin, report_pin);

    // `true` rescores by default.
    for (dir, rescore) in [(&report, false), (&on, true)] {
        let results: serde_json::Value =
            serde_json::from_slice(&std::fs::read(dir.join("results.json"))?)?;
        assert_eq!(results["immonium"]["rescore"], rescore);
        assert_eq!(results["immonium"]["residue_ions"], true);
        assert_matches_published_schemas(&workspace, dir)?;
    }
    std::fs::remove_dir_all(root)?;
    Ok(())
}

fn sage_immonium_pin_columns() -> [&'static str; 5] {
    [
        "immonium_explained",
        "immonium_missing",
        "immonium_unexplained",
        "immonium_modified_explained",
        "immonium_modified_unexplained",
    ]
}

#[test]
fn parquet_footers_and_run_summary_record_provenance() -> anyhow::Result<()> {
    use parquet::file::reader::{FileReader, SerializedFileReader};

    let workspace = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("../..");
    let root = std::env::temp_dir().join(format!(
        "sage-plus-provenance-{}-{}",
        std::process::id(),
        SystemTime::now().duration_since(UNIX_EPOCH)?.as_nanos()
    ));
    std::fs::create_dir_all(&root)?;
    let mut config: serde_json::Value =
        serde_json::from_slice(&std::fs::read(workspace.join("tests/config.json"))?)?;
    config["quant"] = serde_json::json!({ "lfq": true });
    std::fs::write(root.join("config.json"), serde_json::to_vec(&config)?)?;
    let result = Command::new(env!("CARGO_BIN_EXE_sage"))
        .current_dir(&workspace)
        .arg(root.join("config.json"))
        .arg("--output_directory")
        .arg(root.join("output"))
        .arg("--disable-telemetry-i-dont-want-to-improve-sage")
        .output()?;
    assert!(
        result.status.success(),
        "{}",
        String::from_utf8_lossy(&result.stderr)
    );

    let summary: serde_json::Value =
        serde_json::from_slice(&std::fs::read(root.join("output/run-summary.json"))?)?;
    assert_eq!(summary["schema_version"], 9);
    let provenance = &summary["provenance"]["metadata"];
    assert_eq!(provenance["version"], env!("CARGO_PKG_VERSION"));
    let fasta = &provenance["fasta"];
    let fasta_sha256 = fasta["sha256"].as_str().expect("FASTA is always hashed");
    assert_eq!(fasta_sha256.len(), 64);
    assert_eq!(fasta["decoys"]["strategy"], "generated");
    assert_eq!(fasta["decoys"]["decoy_tag"], "rev_");
    assert_eq!(fasta["uniprot"]["organisms"][0]["taxonomy_id"], 9606);
    let input = &provenance["inputs"][0];
    assert_eq!(input["name"], "LQSRPAAPPAPGPGQLTLR.mzML");
    assert!(input["sha256"].is_null());
    assert_eq!(input["sha256_skipped"], "record_input_hashes is off");
    assert!(provenance["config"]["output_paths"].is_null());
    assert_eq!(provenance["config"]["record_input_hashes"], false);

    let mut files = std::fs::read_dir(root.join("output"))?
        .map(|entry| entry.map(|entry| entry.path()))
        .collect::<Result<Vec<_>, _>>()?;
    files.retain(|path| path.extension().is_some_and(|ext| ext == "parquet"));
    files.sort();
    let names = files
        .iter()
        .map(|path| path.file_name().unwrap().to_string_lossy().into_owned())
        .collect::<Vec<_>>();
    assert_eq!(
        names,
        [
            "lfq.parquet",
            "matched_fragments.sage.parquet",
            "results.sage.parquet"
        ]
    );
    for path in files {
        let reader = SerializedFileReader::new(std::fs::File::open(&path)?)?;
        let footer = reader
            .metadata()
            .file_metadata()
            .key_value_metadata()
            .cloned()
            .unwrap_or_default()
            .into_iter()
            .map(|entry| (entry.key, entry.value.unwrap_or_default()))
            .collect::<std::collections::HashMap<_, _>>();
        let name = path.display();
        assert!(footer.contains_key("sage.schema.name"), "{name}");
        assert_eq!(footer["sage.provenance.version"], "1", "{name}");
        assert_eq!(footer["sage.version"], env!("CARGO_PKG_VERSION"), "{name}");
        for key in [
            "sage.config",
            "sage.inputs",
            "sage.fasta",
            "sage.database_inputs",
            "sage.protein_inference",
        ] {
            let value: serde_json::Value = serde_json::from_str(&footer[key])?;
            assert_eq!(&value, &provenance[&key["sage.".len()..]], "{name}: {key}");
        }
        let fasta: serde_json::Value = serde_json::from_str(&footer["sage.fasta"])?;
        assert_eq!(fasta["sha256"], fasta_sha256, "{name}");
    }
    std::fs::remove_dir_all(root)?;
    Ok(())
}

/// Check a run's `run-summary.json` and QC TSV headers against the schemas
/// published in `schemas/`, so the files cannot drift from them.
fn assert_matches_published_schemas(
    workspace: &std::path::Path,
    output: &std::path::Path,
) -> anyhow::Result<()> {
    let schemas = workspace.join("schemas");
    let summary: serde_json::Value =
        serde_json::from_slice(&std::fs::read(output.join("run-summary.json"))?)?;
    let schema: serde_json::Value =
        serde_json::from_slice(&std::fs::read(schemas.join("run-summary.v9.schema.json"))?)?;
    let mut compiler = boon::Compiler::new();
    let mut compiled = boon::Schemas::new();
    compiler
        .add_resource("run-summary.v9.schema.json", schema)
        .map_err(|error| anyhow::anyhow!("{error}"))?;
    let index = compiler
        .compile("run-summary.v9.schema.json", &mut compiled)
        .map_err(|error| anyhow::anyhow!("{error}"))?;
    if let Err(error) = compiled.validate(&summary, index) {
        panic!("run-summary.json does not match its published schema: {error:#}");
    }
    // The schema is not vacuous: a wrongly typed field is rejected.
    let mut broken = summary.clone();
    broken["execution"]["rayon_threads"] = serde_json::json!("two");
    assert!(compiled.validate(&broken, index).is_err());

    for (file, schema, required) in [
        ("digestion.tsv", "digestion.v1.tsv.schema.json", true),
        (
            "diagnostic_ions.tsv",
            "diagnostic_ions.v1.tsv.schema.json",
            false,
        ),
    ] {
        if !required && !output.join(file).exists() {
            continue;
        }
        let table: serde_json::Value =
            serde_json::from_slice(&std::fs::read(schemas.join(schema))?)?;
        let columns = table["fields"]
            .as_array()
            .expect("table schema fields")
            .iter()
            .map(|field| field["name"].as_str().expect("field name"))
            .collect::<Vec<_>>();
        let contents = std::fs::read_to_string(output.join(file))?;
        let header = contents.lines().next().expect("TSV header");
        assert_eq!(header.split('\t').collect::<Vec<_>>(), columns, "{file}");
        for row in contents.lines().skip(1) {
            let values = row.split('\t').collect::<Vec<_>>();
            assert_eq!(values.len(), columns.len(), "{file}: {row}");
            for (value, field) in values.iter().zip(table["fields"].as_array().unwrap()) {
                let parses = match field["type"].as_str().unwrap() {
                    "integer" => value.parse::<u64>().is_ok(),
                    "number" => value.parse::<f64>().is_ok(),
                    _ => true,
                };
                assert!(parses, "{file}: `{value}` is not a {}", field["type"]);
            }
        }
    }
    Ok(())
}

fn approved_fragment_losses() -> serde_json::Value {
    serde_json::json!({
        "Water": {"mass": 18.010565, "sites": ["S", "T", "E", "D"], "ion_kinds": ["b", "y"]},
        "Ammonia": {"mass": 17.026549, "sites": ["R", "K", "N", "Q"], "ion_kinds": ["y"]}
    })
}

/// Generic fragment losses are separate rescoring features: with them
/// configured, every search score of the PSM (hyperscore, matched peaks,
/// intensities, ranks) is unchanged, and without them the loss evidence is
/// absent.
#[test]
fn fragment_losses_add_evidence_without_changing_search_scores() -> anyhow::Result<()> {
    let fasta = sage_cloudpath::util::read_fasta(
        &sage_cloudpath::to_url("../../tests/Q99536.fasta").expect("valid url"),
        "rev_",
        true,
    )?;
    let spectra = sage_cloudpath::util::read_mzml(
        &sage_cloudpath::to_url("../../tests/LQSRPAAPPAPGPGQLTLR.mzML").expect("valid url"),
        0,
        None,
    )?;
    let processed = SpectrumProcessor::new(100, true, 0.0).process(spectra[0].clone());

    let mut features = Vec::new();
    for losses in [None, Some(approved_fragment_losses())] {
        let mut builder = Builder::default();
        builder.update_fasta("foo".into());
        builder.fragment_losses = losses.map(serde_json::from_value).transpose()?;
        builder
            .validate_fragment_losses()
            .map_err(anyhow::Error::msg)?;
        let database = builder.make_parameters().build(fasta.clone());
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
            annotate_matches: true,
            mass_shift_ppm: 50.0,
            score_type: ScoreType::SageHyperScore,
            mass_recalibration: None,
        };
        let psms = scorer.score(&processed);
        assert_eq!(psms.len(), 1);
        features.push(psms.into_iter().next().unwrap());
    }
    let (off, on) = (&features[0], &features[1]);
    assert!(off.fragment_loss.is_none());
    let loss = on.fragment_loss.expect("loss evidence when configured");
    // LQSRPAAPPAPGPGQLTLR carries S, T, Q and R loss sites.
    assert!(loss.matched_peaks > 0, "{loss:?}");
    assert!(loss.intensity_pct > 0.0 && loss.intensity_pct < 100.0);
    // `fragment_loss` is not serialized: every other field but the
    // process-wide `psm_id` counter is identical, including the annotated
    // fragments.
    let fields = |feature| -> anyhow::Result<serde_json::Value> {
        let mut value = serde_json::to_value(feature)?;
        value["psm_id"] = serde_json::Value::Null;
        Ok(value)
    };
    assert_eq!(fields(off)?, fields(on)?);
    assert_eq!(off.matched_peaks, 20);
    assert_eq!(
        off.fragments.as_ref().map(|f| f.mz_experimental.clone()),
        on.fragments.as_ref().map(|f| f.mz_experimental.clone())
    );
    Ok(())
}

/// `results.sage.parquet` always has the published v5 columns; the two loss
/// columns are null without `database.fragment_losses` and filled with it.
/// The pin gains the loss features before `Peptide`, and before the immonium
/// columns when both options are on.
#[test]
fn fragment_loss_columns_are_filled_only_when_configured() -> anyhow::Result<()> {
    use parquet::file::reader::{FileReader, SerializedFileReader};
    let workspace = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("../..");
    let root = std::env::temp_dir().join(format!(
        "sage-plus-fragment-losses-{}-{}",
        std::process::id(),
        SystemTime::now().duration_since(UNIX_EPOCH)?.as_nanos()
    ));
    std::fs::create_dir_all(&root)?;
    let published = sage_cloudpath::parquet::build_schema()?
        .get_fields()
        .iter()
        .map(|field| field.name().to_string())
        .collect::<Vec<_>>();
    let loss_columns = ["matched_loss_peaks", "loss_intensity_pct"];

    let mut outputs = Vec::new();
    for (losses, immonium) in [(false, false), (true, false), (true, true)] {
        let mut config: serde_json::Value =
            serde_json::from_slice(&std::fs::read(workspace.join("tests/config.json"))?)?;
        config["write_pin"] = true.into();
        if losses {
            config["database"]["fragment_losses"] = approved_fragment_losses();
        }
        if immonium {
            config["immonium"] = true.into();
        }
        let run = root.join(format!("losses-{losses}-immonium-{immonium}"));
        std::fs::create_dir_all(&run)?;
        std::fs::write(run.join("config.json"), serde_json::to_vec(&config)?)?;
        let result = Command::new(env!("CARGO_BIN_EXE_sage"))
            .current_dir(&workspace)
            .arg(run.join("config.json"))
            .arg("--output_directory")
            .arg(run.join("output"))
            .arg("--disable-telemetry-i-dont-want-to-improve-sage")
            .args(["--threads", "2"])
            .output()?;
        assert!(
            result.status.success(),
            "{}",
            String::from_utf8_lossy(&result.stderr)
        );
        let pin = std::fs::read_to_string(run.join("output/results.sage.pin"))?;
        let pin_header = pin.lines().next().expect("pin header").to_string();
        let reader = SerializedFileReader::new(std::fs::File::open(
            run.join("output/results.sage.parquet"),
        )?)?;
        let columns = reader
            .metadata()
            .file_metadata()
            .schema()
            .get_fields()
            .iter()
            .map(|field| field.name().to_string())
            .collect::<Vec<_>>();
        assert_eq!(columns, published);
        // Loss and immonium values of every row, as strings.
        let mut values = Vec::new();
        for row in reader.get_row_iter(None)? {
            let row = row?;
            let fields = row
                .get_column_iter()
                .filter(|(name, _)| {
                    loss_columns.contains(&name.as_str()) || name.starts_with("immonium_")
                })
                .map(|(name, field)| (name.clone(), field.to_string()))
                .collect::<Vec<_>>();
            values.push(fields);
        }
        assert!(!values.is_empty());
        for fields in &values {
            for (name, value) in fields {
                let filled = if name.starts_with("immonium_") {
                    immonium
                } else {
                    losses
                };
                assert_eq!(value != "null", filled, "{name} = {value}");
            }
        }
        if losses {
            // The loss features are informative on the test file.
            assert!(values.iter().any(|fields| fields
                .iter()
                .any(|(name, value)| name == "matched_loss_peaks" && value != "0")));
        }
        outputs.push(pin_header);
    }

    let (off_pin, on_pin, both_pin) = (&outputs[0], &outputs[1], &outputs[2]);
    assert!(!off_pin.contains("loss"), "{off_pin}");
    assert_eq!(
        on_pin.as_str(),
        off_pin.replace(
            "posterior_error\tPeptide",
            "posterior_error\tmatched_loss_peaks\tln(loss_intensity_pct)\tPeptide"
        )
    );
    assert_eq!(
        both_pin.as_str(),
        on_pin.replace(
            "ln(loss_intensity_pct)\tPeptide",
            &format!(
                "ln(loss_intensity_pct)\t{}\tPeptide",
                sage_immonium_pin_columns().join("\t")
            )
        )
    );
    std::fs::remove_dir_all(root)?;
    Ok(())
}

#[test]
fn legacy_telemetry_flag_is_a_hidden_no_op() -> anyhow::Result<()> {
    let workspace = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("../..");
    let help = Command::new(env!("CARGO_BIN_EXE_sage"))
        .arg("--help")
        .output()?;
    assert!(help.status.success());
    assert!(!String::from_utf8(help.stdout)?.contains("telemetry"));

    let output = Command::new(env!("CARGO_BIN_EXE_sage"))
        .current_dir(&workspace)
        .arg(workspace.join("tests/config.json"))
        .arg("--validate-only")
        .arg("--disable-telemetry-i-dont-want-to-improve-sage")
        .output()?;
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    assert!(output.stdout.is_empty());
    Ok(())
}

#[test]
fn unusable_output_directory_error_names_the_path() -> anyhow::Result<()> {
    let workspace = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("../..");
    let nonce = SystemTime::now().duration_since(UNIX_EPOCH)?.as_nanos();
    let blocker = std::env::temp_dir().join(format!(
        "sage-plus-output-blocker-{}-{nonce}",
        std::process::id()
    ));
    std::fs::write(&blocker, b"")?;
    let output_directory = blocker.join("results");
    let output = Command::new(env!("CARGO_BIN_EXE_sage"))
        .current_dir(&workspace)
        .arg(workspace.join("tests/config.json"))
        .arg("--output_directory")
        .arg(&output_directory)
        .output()?;
    assert!(!output.status.success());
    let stderr = String::from_utf8(output.stderr)?;
    assert!(
        stderr.contains(&format!(
            "cannot create output directory `{}`",
            output_directory.display()
        )),
        "{stderr}"
    );
    std::fs::remove_file(blocker)?;
    Ok(())
}

#[test]
fn empty_search_spaces_fail_with_the_likely_cause() -> anyhow::Result<()> {
    let workspace = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("../..");
    let nonce = SystemTime::now().duration_since(UNIX_EPOCH)?.as_nanos();
    let root = std::env::temp_dir().join(format!(
        "sage-plus-empty-search-{}-{nonce}",
        std::process::id()
    ));
    std::fs::create_dir_all(&root)?;
    let empty_tsv = root.join("empty.tsv");
    std::fs::write(&empty_tsv, "sequence\n")?;
    let base: serde_json::Value =
        serde_json::from_slice(&std::fs::read(workspace.join("tests/config.json"))?)?;

    let mut lengths = base.clone();
    lengths["database"]["enzyme"]["min_len"] = 30.into();
    lengths["database"]["enzyme"]["max_len"] = 5.into();
    let mut masses = base.clone();
    masses["database"]["peptide_min_mass"] = 5000.0.into();
    masses["database"]["peptide_max_mass"] = 500.0.into();
    let mut peptides = base.clone();
    peptides["database"]
        .as_object_mut()
        .unwrap()
        .remove("fasta");
    peptides["database"]["peptides"] = empty_tsv.to_string_lossy().as_ref().into();

    for (name, config, expected) in [
        (
            "lengths",
            lengths,
            "`database.enzyme.min_len` (30) is greater than `max_len` (5)",
        ),
        (
            "masses",
            masses,
            "`database.peptide_min_mass` (5000) must not exceed",
        ),
        (
            "peptides",
            peptides,
            "the database contains no target peptides",
        ),
    ] {
        let path = root.join(format!("{name}.json"));
        std::fs::write(&path, serde_json::to_vec(&config)?)?;
        let output = Command::new(env!("CARGO_BIN_EXE_sage"))
            .current_dir(&workspace)
            .arg(&path)
            .arg("--output_directory")
            .arg(root.join(format!("output-{name}")))
            .output()?;
        let stderr = String::from_utf8(output.stderr)?;
        assert!(!output.status.success(), "{name}: {stderr}");
        assert!(stderr.contains(expected), "{name}: {stderr}");
    }
    std::fs::remove_dir_all(root)?;
    Ok(())
}
