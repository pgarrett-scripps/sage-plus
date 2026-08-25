use super::{
    assign_psm_ids, average_finite, finish_csv_writer, labeled_finite_values, median_finite,
    normalize_finite, passes_localization_filter, passes_output_filter,
    sort_features_by_discriminant, OutputTarget, RunSummary, SpectrumAccumulator,
};
use rayon::prelude::*;
use sage_cloudpath::Url;
use sage_core::database::PeptideIx;
use sage_core::scoring::Feature;
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
    assert_eq!(summary.models.library_retention_time_alignment, None);
    assert_eq!(summary.models.library_retention_time_files_aligned, 0);
    assert_eq!(summary.models.library_ion_mobility_alignment, None);
    assert_eq!(summary.models.library_ion_mobility_files_aligned, 0);
    assert_eq!(summary.models.library_rescoring, None);
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
