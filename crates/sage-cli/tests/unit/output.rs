use super::*;

fn result(index: usize) -> SageResults {
    SageResults {
        ms1: vec![ProcessedSpectrum {
            id: format!("ms1-{index}"),
            ..Default::default()
        }],
        features: vec![Feature {
            spec_id: format!("feature-{index}"),
            ..Default::default()
        }],
        quant: vec![TmtQuant {
            spec_id: format!("quant-{index}"),
            file_id: index,
            ion_injection_time: index as f32,
            peaks: vec![index as f32],
        }],
    }
}

#[test]
fn sequential_collection_combines_every_result_vector() {
    let combined = (0..3).map(result).collect::<SageResults>();

    assert_eq!(
        combined
            .ms1
            .iter()
            .map(|spectrum| spectrum.id.as_str())
            .collect::<Vec<_>>(),
        vec!["ms1-0", "ms1-1", "ms1-2"]
    );
    assert_eq!(
        combined
            .features
            .iter()
            .map(|feature| feature.spec_id.as_str())
            .collect::<Vec<_>>(),
        vec!["feature-0", "feature-1", "feature-2"]
    );
    assert_eq!(
        combined
            .quant
            .iter()
            .map(|quant| quant.file_id)
            .collect::<Vec<_>>(),
        vec![0, 1, 2]
    );
}

#[test]
fn parallel_collection_keeps_every_partial_result() {
    let combined = (0..100)
        .into_par_iter()
        .map(result)
        .collect::<SageResults>();
    let mut file_ids = combined
        .quant
        .iter()
        .map(|quant| quant.file_id)
        .collect::<Vec<_>>();
    file_ids.sort_unstable();

    assert_eq!(combined.ms1.len(), 100);
    assert_eq!(combined.features.len(), 100);
    assert_eq!(file_ids, (0..100).collect::<Vec<_>>());
}

#[test]
fn existing_artifacts_require_explicit_overwrite_and_unrelated_files_survive() {
    let path = std::env::temp_dir().join(format!(
        "sage-overwrite-{}-{}",
        std::process::id(),
        std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .unwrap()
            .as_nanos()
    ));
    std::fs::create_dir_all(&path).unwrap();
    std::fs::write(path.join("run-summary.json"), "old completion").unwrap();
    std::fs::write(path.join("lfq.parquet"), "old optional output").unwrap();
    std::fs::write(path.join("notes.md"), "user notes").unwrap();
    let url = sage_cloudpath::Url::from_directory_path(&path).unwrap();
    assert!(super::prepare_local_directory(&url, false).is_err());
    assert!(path.join("run-summary.json").exists());
    super::prepare_local_directory(&url, true).unwrap();
    assert!(!path.join("run-summary.json").exists());
    assert!(!path.join("lfq.parquet").exists());
    assert_eq!(
        std::fs::read_to_string(path.join("notes.md")).unwrap(),
        "user notes"
    );
    std::fs::remove_dir_all(path).unwrap();
}
