use super::{
    assign_psm_ids, average_finite, finish_csv_writer, labeled_finite_values, median_finite,
    missing_decoy_warning, normalize_finite, passes_localization_filter, passes_output_filter,
    sort_features_by_discriminant, spectrum_id_occurrences, LabelGroupIndex, OutputTarget,
    RunSummary, SpectrumAccumulator, ToleranceRecommendation,
};

#[test]
fn missing_decoys_produce_an_actionable_warning() {
    let warning = missing_decoy_warning(false, [false, false]).unwrap();
    assert!(warning.contains("generate_decoys"));
    assert!(warning.contains("FDR"));
    assert!(missing_decoy_warning(false, [false, true]).is_none());
    assert!(missing_decoy_warning(true, [false, false])
        .unwrap()
        .contains("non-colliding"));
}

#[test]
fn label_group_closure_keeps_every_channel_partner() {
    let builder: sage_core::database::Builder = serde_json::from_value(serde_json::json!({
        "generate_decoys": false,
        "static_mods": {
            "SILAC-K": {
                "mass": 0.0,
                "channel_offsets": {"light": 0.0, "heavy": 8.014199},
                "sites": ["K"]
            }
        }
    }))
    .unwrap();
    let parameters = builder.make_parameters();
    let peptides = parameters.peptides_from_tsv("sequence\nPEPTIDEK\n");
    assert_eq!(peptides.len(), 2);
    let heavy = peptides
        .iter()
        .position(|peptide| peptide.label_channel.as_deref() == Some("heavy"))
        .unwrap();
    let keep = AtomicBitSet::new(peptides.len());
    keep.insert(heavy);

    LabelGroupIndex::new(&peptides).close(&keep);

    assert!((0..peptides.len()).all(|index| keep.contains(index)));
}

#[test]
fn pair_closure_is_target_decoy_symmetric() {
    let parameters = sage_core::database::Builder::default().make_parameters();
    let peptides = parameters.peptides_from_tsv("sequence\nPEPTIDER\n");
    let database = parameters.build_from_peptides(peptides);
    assert_eq!(
        database
            .peptides
            .iter()
            .filter(|peptide| !peptide.decoy)
            .count(),
        1
    );
    assert_eq!(
        database
            .peptides
            .iter()
            .filter(|peptide| peptide.decoy)
            .count(),
        1
    );
    for selected in 0..database.peptides.len() {
        let keep = AtomicBitSet::new(database.peptides.len());
        keep.insert(selected);

        super::close_prefilter_pairs(&database, &keep);

        assert!((0..database.peptides.len()).all(|index| keep.contains(index)));
    }
}
use rayon::prelude::*;
use sage_cloudpath::Url;
use sage_core::database::PeptideIx;
use sage_core::scoring::{AtomicBitSet, Feature};
use sage_core::spectrum::ProcessedSpectrum;
use std::io::Write;

fn temporary_output(name: &str) -> (std::path::PathBuf, std::path::PathBuf) {
    let unique = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .unwrap()
        .as_nanos();
    let directory =
        std::env::temp_dir().join(format!("sage-runner-test-{}-{unique}", std::process::id()));
    (directory.clone(), directory.join(name))
}

#[test]
fn tied_features_receive_repeatable_psm_ids() {
    let feature = |file_id, spec_id: &str, peptide_idx| Feature {
        file_id,
        spec_id: spec_id.into(),
        rank: 1,
        peptide_idx: PeptideIx(peptide_idx),
        discriminant_score: 5.0,
        ..Feature::default()
    };
    let mut forward = vec![feature(1, "scan=2", 2), feature(0, "scan=1", 1)];
    let mut reversed = forward.iter().cloned().rev().collect::<Vec<_>>();

    sort_features_by_discriminant(&mut forward);
    assign_psm_ids(&mut forward);
    sort_features_by_discriminant(&mut reversed);
    assign_psm_ids(&mut reversed);

    let identities = |features: &[Feature]| {
        features
            .iter()
            .map(|feature| (feature.file_id, feature.spec_id.clone(), feature.psm_id))
            .collect::<Vec<_>>()
    };
    assert_eq!(identities(&forward), identities(&reversed));
}

#[test]
fn localization_filter_requires_passing_target_psm() {
    let passing = Feature {
        label: 1,
        spectrum_q: 0.01,
        ..Default::default()
    };
    assert!(passes_localization_filter(&passing, 0.01));

    let failing = Feature {
        spectrum_q: 0.011,
        ..passing.clone()
    };
    assert!(!passes_localization_filter(&failing, 0.01));

    let decoy = Feature {
        label: -1,
        ..passing
    };
    assert!(!passes_localization_filter(&decoy, 0.01));
}

#[test]
fn output_filter_is_inclusive_and_applies_to_targets_and_decoys() {
    let target = Feature {
        label: 1,
        spectrum_q: 0.1,
        ..Default::default()
    };
    assert!(passes_output_filter(&target, 0.1));

    let decoy = Feature {
        label: -1,
        ..target.clone()
    };
    assert!(passes_output_filter(&decoy, 0.1));

    let failing = Feature {
        spectrum_q: 0.100_001,
        ..target
    };
    assert!(!passes_output_filter(&failing, 0.1));
}

#[test]
fn older_run_summaries_receive_compatible_defaults() {
    let summary: RunSummary = serde_json::from_value(serde_json::json!({
        "runtime_secs": 1,
        "files": 1,
        "peptides_in_database": 10,
        "fragments_in_database": 20,
        "psms_at_one_percent_fdr": 2,
        "peptides_at_one_percent_fdr": 1,
        "proteins_at_one_percent_fdr": 1,
        "protein_groups_at_one_percent_fdr": 1,
        "output_paths": []
    }))
    .unwrap();

    assert_eq!(summary.schema_version, 1);
    assert_eq!(
        summary.recommended_tolerances,
        ToleranceRecommendation::default()
    );
    assert!(!summary.ptm_localization.enabled);
    assert_eq!(summary.quantification.lfq_features, 0);
}

#[test]
fn local_output_target_creates_parents_and_flushes_contents() {
    let (directory, path) = temporary_output("nested/result.txt");
    let url = Url::from_file_path(&path).unwrap();
    let mut output = OutputTarget::new(&url).unwrap();
    output.write_all(b"sage output\n").unwrap();
    output.finish(&url).unwrap();

    assert_eq!(std::fs::read_to_string(&path).unwrap(), "sage output\n");
    std::fs::remove_dir_all(directory).unwrap();
}

#[test]
fn csv_output_is_complete_after_finalization() {
    let (directory, path) = temporary_output("nested/result.csv");
    let url = Url::from_file_path(&path).unwrap();
    let output = OutputTarget::new(&url).unwrap();
    let mut writer = csv::Writer::from_writer(output);
    writer.write_record(["peptide", "score"]).unwrap();
    writer.write_record(["PEPTIDE", "42"]).unwrap();
    finish_csv_writer(writer, &url).unwrap();

    assert_eq!(
        std::fs::read_to_string(&path).unwrap(),
        "peptide,score\nPEPTIDE,42\n"
    );
    std::fs::remove_dir_all(directory).unwrap();
}

#[cfg(feature = "cloud")]
#[test]
fn remote_output_target_flushes_through_cloud_writer() {
    let url = Url::parse("memory:///nested/result.txt").unwrap();
    let mut output = OutputTarget::new(&url).unwrap();
    output.write_all(b"remote sage output\n").unwrap();
    output.finish(&url).unwrap();
}

#[cfg(not(feature = "cloud"))]
#[test]
fn remote_output_target_requires_the_cloud_feature() {
    let url = Url::parse("s3://bucket/results.sage.tsv").unwrap();
    let message = OutputTarget::new(&url).err().unwrap().to_string();
    assert!(message.contains("`cloud` feature"), "{message}");
}

#[test]
fn repeated_spectrum_ids_are_numbered_per_file() {
    let spectrum = |file_id, id: &str| ProcessedSpectrum {
        file_id,
        id: id.into(),
        ..ProcessedSpectrum::default()
    };
    let spectra = vec![
        spectrum(0, "a"),
        spectrum(0, "b"),
        spectrum(0, "a"),
        spectrum(1, "a"),
        spectrum(0, "a"),
    ];
    assert_eq!(spectrum_id_occurrences(&spectra), vec![0, 0, 1, 0, 2]);
}

#[test]
fn spectrum_accumulator_separates_ms1_from_fragment_spectra() {
    let spectra = vec![
        ProcessedSpectrum {
            level: 1,
            id: "ms1-a".into(),
            ..ProcessedSpectrum::default()
        },
        ProcessedSpectrum {
            level: 2,
            id: "ms2-a".into(),
            ..ProcessedSpectrum::default()
        },
        ProcessedSpectrum {
            level: 1,
            id: "ms1-b".into(),
            ..ProcessedSpectrum::default()
        },
        ProcessedSpectrum {
            level: 3,
            id: "ms3-a".into(),
            ..ProcessedSpectrum::default()
        },
    ];

    let sequential = spectra.clone().into_iter().collect::<SpectrumAccumulator>();
    let parallel = spectra.into_par_iter().collect::<SpectrumAccumulator>();
    for accumulator in [sequential, parallel] {
        assert_eq!(accumulator.ms1.len(), 2);
        assert_eq!(accumulator.msn.len(), 2);
        assert!(accumulator.ms1.iter().all(|spectrum| spectrum.level == 1));
        assert!(accumulator.msn.iter().all(|spectrum| spectrum.level > 1));
    }
}

#[test]
fn report_statistics_ignore_nonfinite_values() {
    assert_eq!(
        median_finite([f32::NAN, 3.0, 1.0, f32::INFINITY, 2.0]),
        Some(2.0)
    );
    assert_eq!(average_finite([f32::NAN, 2.0, 4.0]), Some(3.0));
    assert_eq!(median_finite([f32::NAN]), None);
    assert_eq!(average_finite([f32::INFINITY]), None);
}

#[test]
fn constant_report_values_normalize_without_nan() {
    assert_eq!(normalize_finite(vec![5.0]), vec![0.0]);
    assert_eq!(normalize_finite(vec![5.0, 5.0]), vec![0.0, 0.0]);
    assert_eq!(normalize_finite(vec![2.0, 4.0]), vec![0.0, 1.0]);
}

#[test]
fn report_score_series_filters_invalid_values_and_labels() {
    let features = vec![
        Feature {
            label: 1,
            discriminant_score: 2.0,
            ..Feature::default()
        },
        Feature {
            label: -1,
            discriminant_score: f32::NAN,
            ..Feature::default()
        },
        Feature {
            label: 0,
            discriminant_score: 1.0,
            ..Feature::default()
        },
    ];
    let values = labeled_finite_values(&features, |feature| feature.discriminant_score as f64);
    assert_eq!(values, (vec![2.0], vec![1]));
}

#[test]
fn report_keeps_same_basename_files_separate() {
    let (directory, _) = temporary_output("report");
    let workspace = concat!(env!("CARGO_MANIFEST_DIR"), "/../..");
    let input: crate::input::Input = serde_json::from_value(serde_json::json!({
        "database": { "fasta": format!("{workspace}/tests/Q99536.fasta") },
        "precursor_tol": { "ppm": [-10, 10] },
        "fragment_tol": { "ppm": [-10, 10] },
        "mzml_paths": [format!("{workspace}/tests/LQSRPAAPPAPGPGQLTLR.mzML")],
        "output_directory": directory.to_string_lossy(),
    }))
    .unwrap();
    let runner = super::Runner::new(input.build().unwrap(), 1).unwrap();
    let feature = |file_id, charge| Feature {
        file_id,
        charge,
        label: 1,
        peptide_idx: PeptideIx(0),
        peptide_len: 10,
        ..Feature::default()
    };
    // a/run.mzML and b/run.mzML share a basename but are different files.
    let features = vec![feature(0, 2), feature(0, 2), feature(1, 3)];
    let filenames = vec!["run.mzML".to_string(), "run.mzML".to_string()];
    let path = runner.write_report(&features, None, &filenames).unwrap();
    let html = std::fs::read_to_string(path.to_file_path().unwrap()).unwrap();
    std::fs::remove_dir_all(directory).unwrap();

    let body = html.split("<tbody>").nth(1).unwrap();
    let body = body.split("</tbody>").next().unwrap();
    let rows = body
        .split("<tr>")
        .skip(1)
        .map(|row| {
            row.split("<td>")
                .skip(1)
                .map(|cell| cell.split("</td>").next().unwrap())
                .collect::<Vec<_>>()
        })
        .collect::<Vec<_>>();
    assert_eq!(rows.len(), 2);
    // PSM targets and average precursor charge are per input file.
    assert_eq!((rows[0][1], rows[0][12]), ("2", "2"));
    assert_eq!((rows[1][1], rows[1][12]), ("1", "3"));
}

#[test]
fn denoise_warns_about_inputs_it_cannot_change() {
    let (directory, _) = temporary_output("denoise");
    let workspace = concat!(env!("CARGO_MANIFEST_DIR"), "/../..");
    let runner = |denoise: bool, lfq: bool| {
        let input: crate::input::Input = serde_json::from_value(serde_json::json!({
            "database": { "fasta": format!("{workspace}/tests/Q99536.fasta") },
            "precursor_tol": { "ppm": [-10, 10] },
            "fragment_tol": { "ppm": [-10, 10] },
            "mzml_paths": [
                format!("{workspace}/tests/LQSRPAAPPAPGPGQLTLR.mzML"),
                format!("{workspace}/crates/sage-cloudpath/tests/data/bruker/example_dia.d"),
            ],
            "quant": { "lfq": lfq },
            "bruker_config": { "denoise": { "enabled": denoise } },
            "output_directory": directory.to_string_lossy(),
        }))
        .unwrap();
        super::Runner::new(input.build().unwrap(), 1).unwrap()
    };

    assert!(runner(false, false).denoise_warnings().is_empty());
    let warnings = runner(true, true).denoise_warnings();
    assert_eq!(warnings.len(), 1);
    assert!(warnings[0].contains("1 other input file"));
    let warnings = runner(true, false).denoise_warnings();
    assert_eq!(warnings.len(), 2);
    assert!(warnings[1].contains("quant.lfq"));
    std::fs::remove_dir_all(directory).ok();
}

#[test]
fn sidecar_outputs_are_registered_written_once_and_listed() {
    let (directory, _) = temporary_output("sidecar");
    let workspace = concat!(env!("CARGO_MANIFEST_DIR"), "/../..");
    let input: crate::input::Input = serde_json::from_value(serde_json::json!({
        "database": { "fasta": format!("{workspace}/tests/Q99536.fasta") },
        "precursor_tol": { "ppm": [-10, 10] },
        "fragment_tol": { "ppm": [-10, 10] },
        "mzml_paths": [format!("{workspace}/tests/LQSRPAAPPAPGPGQLTLR.mzML")],
        "output_directory": directory.to_string_lossy(),
    }))
    .unwrap();
    let mut runner = super::Runner::new(input.build().unwrap(), 1).unwrap();
    let before = runner.parameters.output_paths.len();

    let error = runner
        .write_sidecar("unregistered.parquet", b"x".to_vec())
        .unwrap_err();
    assert!(error.to_string().contains("not registered"));
    assert!(!directory.join("unregistered.parquet").exists());

    let name = crate::output::SIDECAR_OUTPUTS[0];
    let path = runner.write_sidecar(name, b"sidecar".to_vec()).unwrap();
    assert_eq!(
        std::fs::read(path.to_file_path().unwrap()).unwrap(),
        b"sidecar"
    );
    assert_eq!(runner.parameters.output_paths.len(), before + 1);
    assert_eq!(runner.parameters.output_paths.last(), Some(&path));

    assert!(runner.write_sidecar(name, b"again".to_vec()).is_err());
    assert_eq!(runner.parameters.output_paths.len(), before + 1);
    std::fs::remove_dir_all(directory).unwrap();
}

#[derive(Clone, Default)]
struct EventLog(std::sync::Arc<std::sync::Mutex<Vec<u8>>>);

impl std::io::Write for EventLog {
    fn write(&mut self, buffer: &[u8]) -> std::io::Result<usize> {
        self.0.lock().unwrap().extend_from_slice(buffer);
        Ok(buffer.len())
    }

    fn flush(&mut self) -> std::io::Result<()> {
        Ok(())
    }
}

/// Proteins streamed through batched spectrum indexes must keep exactly the
/// peptides the whole-digest prefilter keeps, and give the same search.
#[test]
fn streamed_prefilter_matches_the_whole_digest_search() -> anyhow::Result<()> {
    use super::prefilter::PrefilterBudgets;
    use crate::events::{CancellationToken, EventEmitter};

    let workspace = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("../..");
    let root = std::env::temp_dir().join(format!(
        "sage-prefilter-stream-{}-{}",
        std::process::id(),
        std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)?
            .as_nanos()
    ));
    std::fs::create_dir_all(&root)?;
    let fasta = root.join("proteins.fasta");
    let q99536 = std::fs::read_to_string(workspace.join("tests/Q99536.fasta"))?;
    let q99536_sequence = q99536
        .lines()
        .filter(|line| !line.starts_with('>'))
        .collect::<String>();
    // A copy of Q99536 shares every digest with it, and a protein carrying a
    // reversed albumin digest (DLGEENFK -> DFNEEGLK) collides with that
    // digest's generated decoy, so both the streamed and shared paths run.
    std::fs::write(
        &fasta,
        q99536
            + "\n>sp|P02768|ALBU_HUMAN\nMKWVTFISLLFLFSSAYSRGVFRRDAHKSEVAHRFKDLGEENFKALVLIAFAQYLQQCPFEDHVK\n"
            + &format!(">sp|COPY|Q99536_COPY\n{q99536_sequence}\n")
            + ">sp|REV|REVERSED_DIGEST\nMSSRDFNEEGLKAGSEYR\n",
    )?;
    // Three files in batches of two: a full batch and a partial one.
    let mzml = workspace.join("tests/LQSRPAAPPAPGPGQLTLR.mzML");
    let inputs = (0..3)
        .map(|idx| {
            let path = root.join(format!("run{idx}.mzML"));
            std::fs::copy(&mzml, &path).map(|_| path.display().to_string())
        })
        .collect::<Result<Vec<_>, _>>()?;
    let mut config: serde_json::Value =
        serde_json::from_slice(&std::fs::read(workspace.join("tests/config.json"))?)?;
    config["database"]["fasta"] = fasta.display().to_string().into();
    config["database"]["prefilter"] = true.into();
    config["mzml_paths"] = serde_json::json!(inputs);
    config["batch_size"] = 2.into();
    let tiny_index = PrefilterBudgets { index_bytes: 1 };

    let search = |name: &str, budgets: Option<PrefilterBudgets>| -> anyhow::Result<_> {
        let mut config = config.clone();
        config["output_directory"] = root.join(name).display().to_string().into();
        let input: crate::input::Input = serde_json::from_value(config)?;
        let mut search = input.build()?;
        search.prefilter_budgets = budgets;
        Ok(search)
    };

    // Survivors, compared directly against the whole-digest prefilter.
    let runner = super::Runner::new_with_control(
        search("peptides", Some(tiny_index))?,
        2,
        EventEmitter::from_writer(EventLog::default()),
        CancellationToken::default(),
    )?;
    let mut database = runner.database_parameters.clone();
    super::load_ptm_library(&mut database)?;
    let proteins = super::load_fasta(&database)?;
    let cleavages = super::load_custom_cleavages(&database, &proteins)?;
    let whole = runner.prefilter_whole_digest(2, &proteins, cleavages.as_ref())?;
    let (streamed, _) = runner.prefilter_peptides(2, proteins, cleavages)?;
    assert!(!whole.is_empty());
    assert!(
        whole.len() < 1000,
        "the prefilter should drop most peptides"
    );
    assert_eq!(streamed, whole);

    let run = |name: &str, budgets: Option<PrefilterBudgets>| -> anyhow::Result<_> {
        let output = root.join(name);
        let log = EventLog::default();
        let runner = super::Runner::new_with_control(
            search(name, budgets)?,
            2,
            EventEmitter::from_writer(log.clone()),
            CancellationToken::default(),
        )?;
        runner.run_with_summary(2)?;
        // Timestamps and job identity differ between runs, and files within
        // a batch may be read in parallel, so only sorted event kinds count.
        let mut events = String::from_utf8(log.0.lock().unwrap().clone())?
            .lines()
            .map(|line| {
                serde_json::from_str::<serde_json::Value>(line)
                    .map(|event| event["event"].as_str().unwrap_or("").to_string())
            })
            .collect::<Result<Vec<_>, _>>()?;
        events.sort();
        let results = ["results.sage.parquet", "matched_fragments.sage.parquet"]
            .map(|file| std::fs::read(output.join(file)))
            .into_iter()
            .collect::<Result<Vec<_>, _>>()?;
        Ok((events, results))
    };

    let (whole_events, whole) = run("whole", None)?;
    let (batched_events, batched) = run("batched", Some(tiny_index))?;
    assert_eq!(whole, batched);
    assert_eq!(whole_events, batched_events);
    for kind in ["file_started", "file_completed"] {
        let count = whole_events.iter().filter(|event| *event == kind).count();
        assert_eq!(count, 3, "{kind} emitted {count} times");
    }
    std::fs::remove_dir_all(root)?;
    Ok(())
}

#[test]
fn report_shows_signed_fragment_bias() {
    let (directory, _) = temporary_output("report-bias");
    let workspace = concat!(env!("CARGO_MANIFEST_DIR"), "/../..");
    let input: crate::input::Input = serde_json::from_value(serde_json::json!({
        "database": { "fasta": format!("{workspace}/tests/Q99536.fasta") },
        "precursor_tol": { "ppm": [-10, 10] },
        "fragment_tol": { "ppm": [-10, 10] },
        "mzml_paths": [format!("{workspace}/tests/LQSRPAAPPAPGPGQLTLR.mzML")],
        "output_directory": directory.to_string_lossy(),
    }))
    .unwrap();
    let runner = super::Runner::new(input.build().unwrap(), 1).unwrap();
    // A uniform -4 ppm fragment shift: the absolute error is 4 ppm, the bias -4.
    let feature = |delta_mass| Feature {
        label: 1,
        peptide_idx: PeptideIx(0),
        delta_mass,
        average_ppm: 4.0,
        signed_fragment_ppm: -4.0,
        ..Feature::default()
    };
    let features = vec![feature(-2.0), feature(-3.0), feature(-2.5)];
    let path = runner
        .write_report(&features, None, &["run.mzML".to_string()])
        .unwrap();
    let html = std::fs::read_to_string(path.to_file_path().unwrap()).unwrap();
    std::fs::remove_dir_all(directory).unwrap();

    assert!(html.contains("Median MS2 Mass Bias (ppm)"));
    assert!(html.contains("Median MS2 Absolute Error (ppm)"));
    let body = html.split("<tbody>").nth(1).unwrap();
    let cells = body
        .split("<td>")
        .skip(1)
        .map(|cell| cell.split("</td>").next().unwrap())
        .collect::<Vec<_>>();
    assert_eq!((cells[6], cells[7], cells[8]), ("-2.5", "-4", "4"));
}

#[test]
fn tolerance_recommendation_uses_confident_signed_errors() {
    let confident = |i: usize| Feature {
        rank: 1,
        label: 1,
        spectrum_q: 0.001,
        // Precursor: -2 ppm bias, errors -3/-2/-1.
        delta_mass: -2.0 + (i % 3) as f32 - 1.0,
        // Fragment: +3 ppm bias; each PSM's ions spread by 2 ppm.
        signed_fragment_ppm: 3.0,
        average_ppm: 3.0,
        fragment_ppm_sd: 2.0,
        ..Feature::default()
    };
    let mut features = (0..150).map(confident).collect::<Vec<_>>();
    // Decoys, lower ranks and non-confident PSMs are ignored.
    for excluded in [
        Feature {
            label: -1,
            delta_mass: 90.0,
            ..confident(0)
        },
        Feature {
            rank: 2,
            delta_mass: 90.0,
            ..confident(0)
        },
        Feature {
            spectrum_q: 0.5,
            delta_mass: 90.0,
            ..confident(0)
        },
    ] {
        features.extend(std::iter::repeat_n(excluded, 60));
    }

    let recommendation = ToleranceRecommendation::from_features(&features);
    assert_eq!(recommendation.psms, 150);
    assert!(recommendation.skipped.is_none());
    let precursor = recommendation.precursor.unwrap();
    assert_eq!(precursor.bias_ppm, -2.0);
    // |-2| + 4 * 1.4826 = 7.93 -> 10 ppm.
    assert_eq!(precursor.recommended_ppm, Some(10.0));
    let fragment = recommendation.fragment.unwrap();
    assert_eq!(fragment.bias_ppm, 3.0);
    // |3| + 4 * 2 = 11 -> 20 ppm.
    assert_eq!(fragment.recommended_ppm, Some(20.0));
    assert_eq!(
        recommendation.log_line(),
        "recommended tolerances: precursor ±10 ppm, fragment ±20 ppm (from 150 PSMs)"
    );

    let few = ToleranceRecommendation::from_features(&features[..99]);
    assert_eq!(few.psms, 99);
    assert!(few.precursor.is_none() && few.fragment.is_none());
    assert!(few.log_line().contains("skipped (99 confident PSMs"));
}

#[test]
fn diapasef_files_get_quality_control_scans() {
    let (directory, _) = temporary_output("tdf-dia-qc");
    let workspace = concat!(env!("CARGO_MANIFEST_DIR"), "/../..");
    let path = format!("{workspace}/crates/sage-cloudpath/tests/data/bruker/example_dia.d");
    let input: crate::input::Input = serde_json::from_value(serde_json::json!({
        "database": { "fasta": format!("{workspace}/tests/Q99536.fasta") },
        "precursor_tol": { "ppm": [-10, 10] },
        "fragment_tol": { "ppm": [-10, 10] },
        "mzml_paths": [path.clone()],
        "dia": { "mode": "pseudo" },
        "quant": { "lfq": true },
        "diagnostic_ions": true,
        "output_directory": directory.to_string_lossy(),
    }))
    .unwrap();
    let runner = super::Runner::new(input.build().unwrap(), 1).unwrap();
    let url = Url::from_file_path(std::fs::canonicalize(&path).unwrap()).unwrap();
    let (ms1, _) = runner
        .read_processed_spectra_with_ms1(&[url], 0, 1, true, false)
        .unwrap();
    assert!(!ms1.is_empty());

    // The raw MS1 frames and window MS2 spectra are scanned before only the
    // MS1 is kept next to the pseudo-spectra.
    let file_qc = runner.file_qc.lock().unwrap();
    let qc = file_qc.get(&0).expect("diaPASEF file scanned");
    assert!(qc.polymers.as_ref().is_some_and(|p| p.ms1_spectra > 0));
    assert!(qc
        .diagnostic_ions
        .as_ref()
        .is_some_and(|scan| scan.ms2_spectra > 0));
    drop(file_qc);
    std::fs::remove_dir_all(directory).ok();
}

#[test]
fn ptm_library_records_match_expanded_ambiguous_residues() {
    let (directory, root) = temporary_output("ptm-ambiguous");
    std::fs::create_dir_all(&root).unwrap();
    let workspace = concat!(env!("CARGO_MANIFEST_DIR"), "/../..");
    let fasta = root.join("proteins.fasta");
    // One-based position 6 is the B of MKPEPBIDEK.
    std::fs::write(&fasta, ">P1\nMKPEPBIDEKRAAAGGGK\n").unwrap();
    let library = root.join("sites.tsv");
    std::fs::write(
        &library,
        "protein\tposition\tresidue\tmodification\nP1\t6\tN\tDeamidated\n",
    )
    .unwrap();
    let runner = |expand: bool| {
        let input: crate::input::Input = serde_json::from_value(serde_json::json!({
            "database": {
                "fasta": fasta.to_string_lossy(),
                "expand_ambiguous_residues": expand,
                "prefilter": false,
                "variable_mods": {
                    "Deamidated": {
                        "mass": 0.984016, "sites": ["N"], "site_mode": "library", "max_count": 1
                    }
                },
                "ptm_library": { "path": library.to_string_lossy(), "strict": true },
            },
            "precursor_tol": { "ppm": [-10, 10] },
            "fragment_tol": { "ppm": [-10, 10] },
            "mzml_paths": [format!("{workspace}/tests/LQSRPAAPPAPGPGQLTLR.mzML")],
            "output_directory": directory.to_string_lossy(),
        }))
        .unwrap();
        super::Runner::new(input.build().unwrap(), 1)
    };

    // Without expansion the B cannot stand for the recorded N.
    let error = runner(false).err().expect("strict library aborts");
    assert!(
        format!("{error:#}").contains("expects residue N"),
        "{error:#}"
    );

    // With expansion the record validates and modifies the N variant only.
    let runner = runner(true).unwrap();
    let modified = runner
        .database
        .peptides
        .iter()
        .filter(|peptide| !peptide.decoy && !peptide.modifications.is_empty())
        .map(|peptide| String::from_utf8_lossy(&peptide.sequence).to_string())
        .collect::<Vec<_>>();
    assert!(!modified.is_empty());
    assert!(
        modified
            .iter()
            .all(|sequence| sequence.contains("PEPNIDEK")),
        "{modified:?}"
    );
    std::fs::remove_dir_all(directory).ok();
}

/// I/L/J twins of the test spectrum's peptide in other proteins merge the
/// same way whether the prefilter streams the digest or not.
#[test]
fn isoleucine_leucine_twins_merge_with_and_without_prefilter() -> anyhow::Result<()> {
    use super::prefilter::PrefilterBudgets;
    use crate::events::{CancellationToken, EventEmitter};

    let (directory, root) = temporary_output("il-merge");
    std::fs::create_dir_all(&root)?;
    let workspace = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("../..");
    let fasta = root.join("proteins.fasta");
    let q99536 = std::fs::read_to_string(workspace.join("tests/Q99536.fasta"))?;
    std::fs::write(
        &fasta,
        q99536
            + "\n>sp|TWIN_I|ISOLEUCINE\nMSSKLQSRPAAPPAPGPGQITLRGGWK\n"
            + ">sp|TWIN_J|AMBIGUOUS\nMWWKLQSRPAAPPAPGPGQJTLRAAYK\n",
    )?;
    let mut config: serde_json::Value =
        serde_json::from_slice(&std::fs::read(workspace.join("tests/config.json"))?)?;
    config["database"]["fasta"] = fasta.display().to_string().into();
    config["database"]["prefilter"] = true.into();
    config["database"]["prefilter_min_matched_peaks"] = 1.into();
    config["mzml_paths"] = serde_json::json!([workspace
        .join("tests/LQSRPAAPPAPGPGQLTLR.mzML")
        .display()
        .to_string()]);
    config["output_directory"] = directory.display().to_string().into();
    let input: crate::input::Input = serde_json::from_value(config)?;
    let mut search = input.build()?;
    // A one-byte index budget streams the digest protein by protein.
    search.prefilter_budgets = Some(PrefilterBudgets { index_bytes: 1 });
    let runner = super::Runner::new_with_control(
        search,
        1,
        EventEmitter::from_writer(EventLog::default()),
        CancellationToken::default(),
    )?;
    let database = runner.database_parameters.clone();
    let proteins = super::load_fasta(&database)?;
    let whole = runner.prefilter_whole_digest(1, &proteins, None)?;
    let (streamed, _) = runner.prefilter_peptides(1, proteins, None)?;
    assert_eq!(streamed, whole);

    let twins = |decoy: bool| {
        whole
            .iter()
            .filter(|peptide| {
                peptide.decoy == decoy
                    && peptide.modifications.is_empty()
                    && sage_core::ambiguous_residues::isoleucine_leucine_eq(
                        &peptide.sequence,
                        match decoy {
                            false => b"LQSRPAAPPAPGPGQLTLR",
                            true => b"LLTLQGPGPAPPAAPRSQR",
                        },
                    )
            })
            .collect::<Vec<_>>()
    };
    let (targets, decoys) = (twins(false), twins(true));
    assert_eq!((targets.len(), decoys.len()), (1, 1));
    // Q99536 sorts first, so its L form is displayed.
    assert_eq!(targets[0].to_string(), "LQSRPAAPPAPGPGQLTLR");
    assert_eq!(targets[0].proteins.len(), 3);
    assert_eq!(decoys[0].proteins, targets[0].proteins);
    std::fs::remove_dir_all(directory).ok();
    Ok(())
}
