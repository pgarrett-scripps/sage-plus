use super::{
    assign_psm_ids, average_finite, finish_csv_writer, labeled_finite_values, median_finite,
    missing_decoy_warning, normalize_finite, passes_localization_filter, passes_output_filter,
    sort_features_by_discriminant, spectrum_id_occurrences, LabelGroupIndex, OutputTarget,
    RunSummary, SpectrumAccumulator,
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
            "K": {
                "mass": 0.0,
                "name": "SILAC-K",
                "channel_offsets": {"light": 0.0, "heavy": 8.014199}
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

#[test]
fn remote_output_target_flushes_through_cloud_writer() {
    let url = Url::parse("memory:///nested/result.txt").unwrap();
    let mut output = OutputTarget::new(&url).unwrap();
    output.write_all(b"remote sage output\n").unwrap();
    output.finish(&url).unwrap();
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
    assert_eq!((rows[0][1], rows[0][11]), ("2", "2"));
    assert_eq!((rows[1][1], rows[1][11]), ("1", "3"));
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
