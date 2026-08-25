use super::{resolve_batch_size, Input, OutputFilter, PtmLocalizationSettings};
use sage_core::{
    database::EnzymeBuilder,
    enzyme::EnzymeParameters,
    ml::retention_alignment::AlignmentMethod,
    ml::retention_model::{RetentionTimeFeatureSet, RetentionTimeSettings},
    spectral_library::{SpectralLibraryFormat, SpectralLibrarySettings, SpectralLibraryStrategy},
};

#[test]
fn deserialize_enriched_retention_time_settings() -> Result<(), serde_json::Error> {
    let settings: RetentionTimeSettings = serde_json::from_value(serde_json::json!({
        "features": "additive_ptm",
        "folds": 5,
        "seed": 7,
        "ptm_regularization": 12.5
    }))?;

    assert_eq!(settings.features, RetentionTimeFeatureSet::AdditivePtm);
    assert_eq!(settings.folds, 5);
    assert_eq!(settings.seed, 7);
    assert_eq!(settings.ptm_regularization, 12.5);
    Ok(())
}

#[test]
fn deserialize_ptm_localization_settings() -> Result<(), serde_json::Error> {
    let configured: PtmLocalizationSettings = serde_json::from_value(serde_json::json!({
        "enabled": true,
        "psm_q_value": 0.025,
        "localization_q_value": 0.05
    }))?;
    assert!(configured.enabled);
    assert_eq!(configured.psm_q_value, 0.025);
    assert_eq!(configured.localization_q_value, 0.05);

    let partial: PtmLocalizationSettings =
        serde_json::from_value(serde_json::json!({ "enabled": true }))?;
    assert!(partial.enabled);
    assert_eq!(partial.psm_q_value, 0.01);
    assert_eq!(partial.localization_q_value, 0.01);

    let legacy_name: PtmLocalizationSettings =
        serde_json::from_value(serde_json::json!({ "q_value": 0.02 }))?;
    assert_eq!(legacy_name.psm_q_value, 0.02);
    Ok(())
}

#[test]
fn deserialize_spectral_library_settings() -> Result<(), serde_json::Error> {
    let configured: SpectralLibrarySettings = serde_json::from_value(serde_json::json!({
        "enabled": true,
        "strategy": "consensus",
        "max_fragments": 12,
        "min_consensus_psms": 2,
        "min_fragment_frequency": 0.6,
        "formats": ["sage_parquet", "mzspeclib"]
    }))?;
    assert!(configured.enabled);
    assert_eq!(configured.psm_q_value, 0.01);
    assert_eq!(configured.peptide_q_value, 0.01);
    assert_eq!(configured.strategy, SpectralLibraryStrategy::Consensus);
    assert_eq!(configured.max_fragments, 12);
    assert_eq!(configured.min_consensus_psms, 2);
    assert_eq!(configured.min_fragment_frequency, 0.6);
    assert_eq!(
        configured.formats,
        vec![
            SpectralLibraryFormat::SageParquet,
            SpectralLibraryFormat::MzSpecLib
        ]
    );
    Ok(())
}

#[test]
fn spectral_library_settings_are_validated() {
    let input: Input = serde_json::from_value(serde_json::json!({
        "database": { "fasta": "test.fasta" },
        "precursor_tol": { "ppm": [-10, 10] },
        "fragment_tol": { "ppm": [-10, 10] },
        "mzml_paths": ["test.mzML"],
        "spectral_library": { "enabled": true, "max_fragments": 0 }
    }))
    .unwrap();
    let error = input.validate().unwrap_err().to_string();
    assert!(error.contains("spectral_library.max_fragments"));
}

fn base_search_space(value: serde_json::Value) -> Input {
    serde_json::from_value(serde_json::json!({
        "precursor_tol": { "ppm": [-10, 10] },
        "fragment_tol": { "ppm": [-10, 10] },
        "mzml_paths": ["test.mzML"],
        "database": value
    }))
    .unwrap()
}

#[test]
fn database_and_library_search_are_mutually_exclusive() {
    let database = base_search_space(serde_json::json!({ "fasta": "test.fasta" }));
    assert!(database.validate().is_ok());

    let library: Input = serde_json::from_value(serde_json::json!({
        "library_search": { "path": "library.mzspeclib.txt" },
        "precursor_tol": { "ppm": [-10, 10] },
        "fragment_tol": { "ppm": [-10, 10] },
        "isotope_errors": [-1, 2],
        "mzml_paths": ["tests/LQSRPAAPPAPGPGQLTLR.mzML"]
    }))
    .unwrap();
    assert!(library.validate().is_ok());

    let both: Input = serde_json::from_value(serde_json::json!({
        "database": { "fasta": "test.fasta" },
        "library_search": { "path": "library.sage.parquet" },
        "precursor_tol": { "ppm": [-10, 10] },
        "fragment_tol": { "ppm": [-10, 10] },
        "mzml_paths": ["test.mzML"]
    }))
    .unwrap();
    assert!(both
        .validate()
        .unwrap_err()
        .to_string()
        .contains("exactly one"));

    let neither: Input = serde_json::from_value(serde_json::json!({
        "precursor_tol": { "ppm": [-10, 10] },
        "fragment_tol": { "ppm": [-10, 10] },
        "mzml_paths": ["test.mzML"]
    }))
    .unwrap();
    assert!(neither
        .validate()
        .unwrap_err()
        .to_string()
        .contains("exactly one"));
}

#[test]
fn modification_channel_offsets_are_validated_before_search() {
    let valid = base_search_space(serde_json::json!({
        "fasta": "test.fasta",
        "static_mods": {
            "R": {
                "mass": 0.0,
                "channel_offsets": {"light": 0.0, "heavy": 10.008269}
            }
        }
    }));
    assert!(valid.validate().is_ok());

    let invalid = base_search_space(serde_json::json!({
        "fasta": "test.fasta",
        "static_mods": {
            "R": {
                "mass": 0.0,
                "channel_offsets": {"light": 0.0, "heavy": 10.008269}
            }
        },
        "variable_mods": {
            "K": [{
                "mass": 0.0,
                "channel_offsets": {"light": 0.0, "medium": 4.025107, "heavy": 8.014199}
            }]
        }
    }));
    assert!(invalid
        .validate()
        .unwrap_err()
        .to_string()
        .contains("same channel names"));
}

#[test]
fn library_search_build_does_not_create_database_parameters() {
    let spectra = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../../tests/LQSRPAAPPAPGPGQLTLR.mzML");
    let input: Input = serde_json::from_value(serde_json::json!({
        "library_search": { "path": "library.sage.parquet" },
        "precursor_tol": { "ppm": [-10, 10] },
        "fragment_tol": { "ppm": [-10, 10] },
        "isotope_errors": [-1, 2],
        "mzml_paths": [spectra]
    }))
    .unwrap();
    let search = input.build().unwrap();
    assert!(search.database.is_none());
    assert!(search.library_search.is_some());
    assert_eq!(search.isotope_errors, (-1, 2));
    assert!(!search.predict_rt);
}

#[test]
fn output_filter_defaults_and_deserializes() -> Result<(), serde_json::Error> {
    let default: OutputFilter = serde_json::from_value(serde_json::json!({}))?;
    assert_eq!(default.psm_q_value, 0.1);

    let configured: OutputFilter =
        serde_json::from_value(serde_json::json!({ "psm_q_value": 0.025 }))?;
    assert_eq!(configured.psm_q_value, 0.025);
    Ok(())
}

#[test]
fn output_filter_q_value_must_be_a_probability() {
    let input: Input = serde_json::from_value(serde_json::json!({
        "database": { "fasta": "test.fasta" },
        "precursor_tol": { "ppm": [-10, 10] },
        "fragment_tol": { "ppm": [-10, 10] },
        "mzml_paths": ["test.mzML"],
        "output_filter": { "psm_q_value": 1.1 }
    }))
    .unwrap();

    let error = input.validate().unwrap_err().to_string();
    assert!(error.contains("output_filter.psm_q_value"));
}

#[test]
fn lfq_numeric_settings_are_validated_before_search() {
    for (setting, value) in [
        ("ppm_tolerance", -1.0),
        ("rt_pct_tolerance", 0.0),
        ("mobility_pct_tolerance", -0.5),
        ("spectral_angle", 1.1),
        ("peptide_q_value", -0.1),
    ] {
        let input: Input = serde_json::from_value(serde_json::json!({
            "database": { "fasta": "test.fasta" },
            "precursor_tol": { "ppm": [-10, 10] },
            "fragment_tol": { "ppm": [-10, 10] },
            "mzml_paths": ["test.mzML"],
            "quant": {
                "lfq": true,
                "lfq_settings": { (setting): value }
            }
        }))
        .unwrap();

        let error = input.validate().unwrap_err().to_string();
        assert!(error.contains(&format!("lfq_settings.{setting}")));
    }
}

#[test]
fn deserialize_enzyme_builder() -> Result<(), serde_json::Error> {
    let a: EnzymeBuilder = serde_json::from_value(serde_json::json!({
        "cleave_at": "KR",
    }))?;
    let b: EnzymeBuilder = serde_json::from_value(serde_json::json!({
        "cleave_at": "KR",
        "restrict": "P",
    }))?;
    let c: EnzymeBuilder = serde_json::from_value(serde_json::json!({
        "cleave_at": "KR",
        "restrict": "",
    }))?;

    let a: EnzymeParameters = a.into();
    let b: EnzymeParameters = b.into();
    let c: EnzymeParameters = c.into();

    assert_eq!(a.enzyme.map(|e| e.skip_suffix), Some([false; 26]));
    {
        let mut expected = [false; 26];
        expected[(b'P' - b'A') as usize] = true;
        assert_eq!(b.enzyme.map(|e| e.skip_suffix), Some(expected));
    }
    assert_eq!(c.enzyme.map(|e| e.skip_suffix), Some([false; 26]));

    Ok(())
}

#[test]
fn deserialize_custom_cleavage_site_path() -> Result<(), serde_json::Error> {
    let input: Input = serde_json::from_value(serde_json::json!({
        "database": {
            "fasta": "proteome.fasta",
            "custom_cleavage_sites": "cleavage-sites.tsv"
        },
        "precursor_tol": { "ppm": [-10.0, 10.0] },
        "fragment_tol": { "ppm": [-20.0, 20.0] },
        "mzml_paths": ["input.mzML"]
    }))?;

    assert_eq!(
        input
            .database
            .as_ref()
            .and_then(|database| database.custom_cleavage_sites.as_deref()),
        Some("cleavage-sites.tsv")
    );
    assert!(input.validate().is_ok());
    Ok(())
}

#[test]
fn deserialize_runtime_memory_settings() -> Result<(), serde_json::Error> {
    let input: Input = serde_json::from_value(serde_json::json!({
        "database": {},
        "precursor_tol": { "ppm": [-10.0, 10.0] },
        "fragment_tol": { "ppm": [-20.0, 20.0] },
        "max_memory_gb": 12.5,
        "min_free_memory_gb": 2.0,
        "batch_size": 1
    }))?;

    assert_eq!(input.max_memory_gb, Some(12.5));
    assert_eq!(input.min_free_memory_gb, Some(2.0));
    assert_eq!(input.batch_size, Some(1));
    assert!(input.memory_limits().unwrap().is_enabled());
    Ok(())
}

#[test]
fn batch_size_must_be_positive() {
    assert!(resolve_batch_size(Some(0)).is_err());
    assert_eq!(resolve_batch_size(Some(3)).unwrap(), 3);
    assert!(resolve_batch_size(None).unwrap() >= 1);
}

#[test]
fn deserialize_nonlinear_retention_time_alignment() -> Result<(), serde_json::Error> {
    let method: AlignmentMethod = serde_json::from_value(serde_json::json!("nonlinear"))?;
    assert_eq!(method, AlignmentMethod::Nonlinear);
    Ok(())
}

#[test]
fn validation_returns_range_errors_instead_of_exiting() {
    let input: Input = serde_json::from_value(serde_json::json!({
        "database": { "fasta": "test.fasta" },
        "precursor_tol": { "ppm": [-10, 10] },
        "fragment_tol": { "ppm": [-10, 10] },
        "isotope_errors": [3, -1],
        "mzml_paths": ["test.mzML"]
    }))
    .unwrap();

    let error = input.validate().unwrap_err().to_string();
    assert!(error.contains("isotope errors"));
}
