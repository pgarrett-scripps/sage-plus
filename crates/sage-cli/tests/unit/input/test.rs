use super::{resolve_batch_size, resolve_threads, Input, OutputFilter, PtmLocalizationSettings};
use sage_core::diagnostic::{default_ions, DiagnosticIonsConfig};
use sage_core::{
    database::EnzymeBuilder,
    enzyme::EnzymeParameters,
    ml::retention_alignment::AlignmentMethod,
    ml::retention_model::{RetentionTimeFeatureSet, RetentionTimeSettings},
    spectral_library::{SpectralLibraryFormat, SpectralLibrarySettings, SpectralLibraryStrategy},
};

#[test]
fn invalid_numeric_settings_fail_logical_validation() {
    let fixture = serde_json::json!({
        "database": {"fasta": "missing.fasta"},
        "mzml_paths": ["missing.mzML"],
        "precursor_tol": {"ppm": [-10, 10]},
        "fragment_tol": {"ppm": [-20, 20]}
    });
    for (key, value) in [
        ("precursor_tol", serde_json::json!({"pct": [-1, 1]})),
        ("fragment_tol", serde_json::json!({"ppm": [20, -20]})),
        ("precursor_charge", serde_json::json!([0, 3])),
        ("max_fragment_charge", serde_json::json!(0)),
        ("min_peaks", serde_json::json!(151)),
        ("max_peaks", serde_json::json!(14)),
        ("protein_grouping_peptide_fdr", serde_json::json!(1.1)),
        ("retention_time_model", serde_json::json!({"folds": 1})),
        (
            "ion_mobility_model",
            serde_json::json!({"ptm_regularization": -1}),
        ),
    ] {
        let mut config = fixture.clone();
        config[key] = value;
        let input: Input = serde_json::from_value(config).unwrap();
        assert!(input.validate().is_err(), "accepted invalid {key}");
    }
    let mut input: Input = serde_json::from_value(fixture).unwrap();
    input.precursor_tol = sage_core::mass::Tolerance::Ppm(1.0, 10.0);
    assert!(
        input.validate().is_ok(),
        "valid asymmetric tolerance rejected"
    );
    input.precursor_tol = sage_core::mass::Tolerance::Ppm(f32::NAN, 10.0);
    assert!(input.validate().is_err());
}

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

#[test]
fn deisotope_boolean_uses_scored_defaults() {
    let input: Input = serde_json::from_value(serde_json::json!({
        "database": { "fasta": "test.fasta" },
        "precursor_tol": { "ppm": [-10, 10] },
        "fragment_tol": { "ppm": [-10, 10] },
        "mzml_paths": ["test.mzML"],
        "deisotope": true
    }))
    .unwrap();
    let settings = input.deisotope.unwrap().resolve();

    assert!(settings.enabled);
    assert_eq!(settings, sage_core::spectrum::DeisotopeSettings::default());
}

#[test]
fn deisotope_object_uses_scored_defaults_and_overrides() {
    let input: Input = serde_json::from_value(serde_json::json!({
        "database": { "fasta": "test.fasta" },
        "precursor_tol": { "ppm": [-10, 10] },
        "fragment_tol": { "ppm": [-10, 10] },
        "mzml_paths": ["test.mzML"],
        "deisotope": {
            "enabled": true,
            "ppm_tolerance": 7.5,
            "max_envelope_peaks": 3,
            "min_score": 0.6
        }
    }))
    .unwrap();
    let settings = input.deisotope.unwrap().resolve();

    assert_eq!(settings.ppm_tolerance, 7.5);
    assert_eq!(settings.max_envelope_peaks, 3);
    assert_eq!(settings.min_score, 0.6);
}

#[test]
fn deisotope_settings_are_validated() {
    let input: Input = serde_json::from_value(serde_json::json!({
        "database": { "fasta": "test.fasta" },
        "precursor_tol": { "ppm": [-10, 10] },
        "fragment_tol": { "ppm": [-10, 10] },
        "mzml_paths": ["test.mzML"],
        "deisotope": {
            "min_envelope_peaks": 4,
            "max_envelope_peaks": 3
        }
    }))
    .unwrap();

    let error = input.validate().unwrap_err().to_string();
    assert!(error.contains("min_envelope_peaks"));
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
fn database_is_required() {
    let database = base_search_space(serde_json::json!({ "fasta": "test.fasta" }));
    assert!(database.validate().is_ok());

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
        .contains("database"));
}

#[test]
fn modification_channel_offsets_are_validated_before_search() {
    let valid = base_search_space(serde_json::json!({
        "fasta": "test.fasta",
        "static_mods": {
            "Arg10": {
                "mass": 0.0,
                "channel_offsets": {"light": 0.0, "heavy": 10.008269},
                "sites": ["R"]
            }
        }
    }));
    assert!(valid.validate().is_ok());

    let invalid = base_search_space(serde_json::json!({
        "fasta": "test.fasta",
        "static_mods": {
            "Arg10": {
                "mass": 0.0,
                "channel_offsets": {"light": 0.0, "heavy": 10.008269},
                "sites": ["R"]
            }
        },
        "variable_mods": {
            "Lys": {
                "mass": 0.0,
                "channel_offsets": {"light": 0.0, "medium": 4.025107, "heavy": 8.014199},
                "sites": ["K"]
            }
        }
    }));
    assert!(invalid
        .validate()
        .unwrap_err()
        .to_string()
        .contains("same channel names"));
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
fn threads_config_parses_and_cli_takes_precedence() -> Result<(), serde_json::Error> {
    let base = serde_json::json!({
        "database": { "fasta": "test.fasta" },
        "precursor_tol": { "ppm": [-10.0, 10.0] },
        "fragment_tol": { "ppm": [-20.0, 20.0] },
        "mzml_paths": ["test.mzML"],
    });
    let unset: Input = serde_json::from_value(base.clone())?;
    assert_eq!(unset.threads, None);
    let mut config = base;
    config["threads"] = serde_json::json!(3);
    let input: Input = serde_json::from_value(config)?;
    assert_eq!(input.threads, Some(3));
    assert!(input.validate().is_ok());

    // `--threads` wins; the config applies without it; neither leaves the
    // choice to Rayon (RAYON_NUM_THREADS, then all cores).
    assert_eq!(resolve_threads(Some(8), Some(3)), Some(8));
    assert_eq!(resolve_threads(None, Some(3)), Some(3));
    assert_eq!(resolve_threads(Some(8), None), Some(8));
    assert_eq!(resolve_threads(None, None), None);
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

#[test]
fn predict_rt_default_matches_documentation() {
    let output = std::env::temp_dir().join(format!(
        "sage-cli-predict-rt-default-{}-{}",
        std::process::id(),
        std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .unwrap()
            .as_nanos()
    ));
    let mut input = base_search_space(serde_json::json!({ "fasta": "test.fasta" }));
    input.output_directory = Some(output.to_string_lossy().into_owned());
    input.mzml_paths = Some(vec![concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/../../tests/LQSRPAAPPAPGPGQLTLR.mzML"
    )
    .into()]);
    let search = input.build().unwrap();
    std::fs::remove_dir_all(output).unwrap();
    assert!(search.predict_rt);

    let docs = include_str!("../../../../../DOCS.md");
    let documented = docs
        .lines()
        .filter(|line| line.contains("predict_rt"))
        .collect::<Vec<_>>();
    assert!(!documented.is_empty());
    for line in documented {
        assert!(!line.contains("default: false"), "{line}");
        assert!(!line.contains("default=false"), "{line}");
    }
}

#[test]
fn bruker_denoise_settings_are_validated() {
    let input: Input = serde_json::from_value(serde_json::json!({
        "database": { "fasta": "test.fasta" },
        "precursor_tol": { "ppm": [-10, 10] },
        "fragment_tol": { "ppm": [-10, 10] },
        "mzml_paths": ["test.d"],
        "bruker_config": { "denoise": { "enabled": true, "halo_peak_fraction": 2.0 } }
    }))
    .unwrap();
    let error = input.validate().unwrap_err().to_string();
    assert!(error.contains("bruker_config.denoise"), "{error}");

    let input: Result<Input, _> = serde_json::from_value(serde_json::json!({
        "bruker_config": { "denoise": { "enabled": true, "halo_fraction": 0.2 } }
    }));
    assert!(input.is_err());
}

fn build_lqsr(dia: Option<serde_json::Value>, name: &str) -> super::Search {
    let output = std::env::temp_dir().join(format!(
        "sage-cli-dia-{name}-{}-{}",
        std::process::id(),
        std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .unwrap()
            .as_nanos()
    ));
    let mut config = serde_json::json!({
        "database": { "fasta": "test.fasta" },
        "precursor_tol": { "ppm": [-10, 10] },
        "fragment_tol": { "ppm": [-10, 10] },
        "wide_window": true,
        "chimera": true,
        "output_directory": output.to_string_lossy(),
        "mzml_paths": [concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/../../tests/LQSRPAAPPAPGPGQLTLR.mzML"
        )]
    });
    if let Some(dia) = dia {
        config["dia"] = dia;
    }
    let input: Input = serde_json::from_value(config).unwrap();
    let search = input.build().unwrap();
    std::fs::remove_dir_all(output).unwrap();
    search
}

#[test]
fn dia_defaults_off_and_is_omitted_from_results_json() {
    let search = build_lqsr(None, "off");
    assert!(search.dia.is_off());
    assert!(search.wide_window && search.chimera);
    let json = serde_json::to_value(&search).unwrap();
    assert!(json.get("dia").is_none());
    let explicit = build_lqsr(Some(serde_json::json!({ "mode": "off" })), "off-explicit");
    assert!(serde_json::to_value(&explicit)
        .unwrap()
        .get("dia")
        .is_none());
}

#[test]
fn dia_pseudo_searches_closed_and_is_recorded() {
    let search = build_lqsr(
        Some(serde_json::json!({ "mode": "pseudo", "min_corr": 0.6 })),
        "pseudo",
    );
    assert_eq!(search.dia.mode, sage_dia::DiaMode::Pseudo);
    assert!(!search.wide_window && !search.chimera);
    let json = serde_json::to_value(&search).unwrap();
    assert_eq!(json["dia"]["mode"], "pseudo");
    assert_eq!(json["dia"]["min_corr"], 0.6f32 as f64);
}

#[test]
fn dia_settings_are_validated() {
    let input: Input = serde_json::from_value(serde_json::json!({
        "database": { "fasta": "test.fasta" },
        "precursor_tol": { "ppm": [-10, 10] },
        "fragment_tol": { "ppm": [-10, 10] },
        "mzml_paths": ["test.mzML"],
        "dia": { "mode": "pseudo", "min_corr": 1.5 }
    }))
    .unwrap();
    let error = input.validate().unwrap_err().to_string();
    assert!(error.contains("dia.min_corr"));
    assert!(serde_json::from_value::<Input>(serde_json::json!({
        "precursor_tol": { "ppm": [-10, 10] },
        "fragment_tol": { "ppm": [-10, 10] },
        "dia": { "mode": "tiered" }
    }))
    .is_err());
}

#[test]
fn each_logical_constraint_reports_its_own_message() {
    let fixture = serde_json::json!({
        "database": {"fasta": "missing.fasta"},
        "mzml_paths": ["missing.mzML"],
        "precursor_tol": {"ppm": [-10, 10]},
        "fragment_tol": {"ppm": [-20, 20]}
    });
    let valid: Input = serde_json::from_value(fixture.clone()).unwrap();
    assert!(valid.validate().is_ok());
    for (key, value, message) in [
        (
            "database",
            serde_json::json!({}),
            "Either `database.fasta` or `database.peptides`",
        ),
        (
            "database",
            serde_json::json!({"peptides": "p.tsv", "custom_cleavage_sites": "c.tsv"}),
            "`database.custom_cleavage_sites` requires `database.fasta`",
        ),
        (
            "mzml_paths",
            serde_json::json!([]),
            "`mzml_paths` must contain at least one",
        ),
        ("isotope_errors", serde_json::json!([3, 1]), "[3, 1]"),
        ("precursor_charge", serde_json::json!([4, 2]), "[4, 2]"),
        ("mass_shift_ppm", serde_json::json!(-1.0), "mass_shift_ppm"),
        (
            "ion_mobility_model",
            serde_json::json!({"min_training_psms": 0}),
            "min_training_psms",
        ),
        (
            "retention_time_model",
            serde_json::json!({"folds": 11}),
            "retention_time_model.folds",
        ),
        ("report_psms", serde_json::json!(0), "report_psms"),
        (
            "ptm_localization",
            serde_json::json!({"psm_q_value": 2.0}),
            "ptm_localization.psm_q_value",
        ),
        (
            "ptm_localization",
            serde_json::json!({"localization_q_value": -0.5}),
            "ptm_localization.localization_q_value",
        ),
        ("max_memory_gb", serde_json::json!(-1.0), ""),
        ("batch_size", serde_json::json!(0), "batch_size"),
        (
            "threads",
            serde_json::json!(0),
            "`threads` must be greater than zero",
        ),
        (
            "precursor_tol",
            serde_json::json!({"pct": [-1, 1]}),
            "percentage precursor tolerances",
        ),
        (
            "fragment_tol",
            serde_json::json!({"da": [0.5, -0.5]}),
            "`fragment_tol` must contain finite ordered bounds",
        ),
    ] {
        let mut config = fixture.clone();
        config[key] = value;
        let input: Input = serde_json::from_value(config).unwrap();
        let error = input
            .validate()
            .expect_err(&format!("accepted invalid {key}"))
            .to_string();
        assert!(error.contains(message), "{key}: {error}");
    }

    // Percentage tolerances remain legal for fragments.
    let mut config = fixture.clone();
    config["fragment_tol"] = serde_json::json!({"pct": [-1, 1]});
    let input: Input = serde_json::from_value(config).unwrap();
    assert!(input.validate().is_ok());

    // A peptide list stands in for a FASTA.
    let mut config = fixture;
    config["database"] = serde_json::json!({"peptides": "p.tsv"});
    let input: Input = serde_json::from_value(config).unwrap();
    assert!(input.validate().is_ok());
}

#[test]
fn quant_options_fill_unset_fields_from_defaults() {
    use super::{QuantOptions, QuantSettings, TmtSettings};
    let options: QuantOptions = serde_json::from_value(serde_json::json!({
        "lfq": true,
        "tmt_settings": {"sn": true},
        "lfq_settings": {"ppm_tolerance": 7.5, "mbr": false}
    }))
    .unwrap();
    let settings: QuantSettings = options.into();
    assert!(settings.lfq);
    assert!(settings.tmt.is_none());
    assert_eq!(settings.tmt_settings.level, 3);
    assert!(settings.tmt_settings.sn);
    assert_eq!(settings.lfq_settings.ppm_tolerance, 7.5);
    assert!(!settings.lfq_settings.mbr);
    let lfq_default = sage_core::lfq::LfqSettings::default();
    assert_eq!(
        settings.lfq_settings.rt_pct_tolerance,
        lfq_default.rt_pct_tolerance
    );

    let empty: QuantSettings = QuantOptions::default().into();
    assert!(!empty.lfq);
    let tmt_default = TmtSettings::default();
    assert_eq!(
        (empty.tmt_settings.level, empty.tmt_settings.sn),
        (tmt_default.level, tmt_default.sn)
    );
    assert_eq!((tmt_default.level, tmt_default.sn), (3, false));

    let level_two: TmtSettings =
        serde_json::from_value::<super::TmtOptions>(serde_json::json!({"level": 2}))
            .unwrap()
            .into();
    assert_eq!((level_two.level, level_two.sn), (2, false));
}

#[test]
fn unsupported_enzyme_residues_fail_validation() {
    for (enzyme, expected) in [
        (serde_json::json!({"cleave_at": "KRB"}), "cleave_at"),
        (serde_json::json!({"restrict": "X"}), "restrict"),
    ] {
        let input: Input = serde_json::from_value(serde_json::json!({
            "database": { "fasta": "test.fasta", "enzyme": enzyme },
            "precursor_tol": { "ppm": [-10, 10] },
            "fragment_tol": { "ppm": [-10, 10] },
            "mzml_paths": ["test.mzML"]
        }))
        .unwrap();
        let error = input.validate().unwrap_err().to_string();
        assert!(error.contains(expected), "{error}");
        assert!(error.contains("unsupported residues"), "{error}");
    }
}

#[test]
fn diagnostic_ions_config() -> anyhow::Result<()> {
    let fixture = serde_json::json!({
        "database": {"fasta": "tests/Q99536.fasta"},
        "mzml_paths": ["tests/LQSRPAAPPAPGPGQLTLR.mzML"],
        "precursor_tol": {"ppm": [-10, 10]},
        "fragment_tol": {"ppm": [-20, 20]}
    });
    let parse = |value: Option<serde_json::Value>| {
        let mut config = fixture.clone();
        if let Some(value) = value {
            config["diagnostic_ions"] = value;
        }
        serde_json::from_value::<Input>(config).unwrap()
    };
    assert!(parse(None).diagnostic_ions.is_none());
    assert_eq!(
        parse(Some(serde_json::json!(true)))
            .diagnostic_ions
            .and_then(DiagnosticIonsConfig::resolve),
        Some(default_ions())
    );
    assert_eq!(
        parse(Some(serde_json::json!(false)))
            .diagnostic_ions
            .and_then(DiagnosticIonsConfig::resolve),
        None
    );
    let custom = parse(Some(serde_json::json!([
        {"name": "TMT126", "mz": 126.1277, "tolerance": {"da": [-0.005, 0.005]}},
        {"name": "HexNAc", "mz": 204.0867}
    ])));
    custom.validate()?;
    let ions = custom
        .diagnostic_ions
        .and_then(DiagnosticIonsConfig::resolve)
        .unwrap();
    assert_eq!(ions.len(), 2);
    assert_eq!(
        ions[1].tolerance(),
        sage_core::mass::Tolerance::Ppm(-20.0, 20.0)
    );

    for invalid in [
        serde_json::json!([]),
        serde_json::json!([{"name": "", "mz": 204.0867}]),
        serde_json::json!([{"name": "bad", "mz": -1.0}]),
        serde_json::json!([{"name": "bad", "mz": 204.0, "tolerance": {"ppm": [5, 10]}}]),
    ] {
        assert!(
            parse(Some(invalid.clone())).validate().is_err(),
            "{invalid}"
        );
    }
    let unknown_field = serde_json::json!([{"name": "x", "mz": 1.0, "charge": 2}]);
    let mut config = fixture.clone();
    config["diagnostic_ions"] = unknown_field;
    assert!(serde_json::from_value::<Input>(config).is_err());
    Ok(())
}

#[test]
fn immonium_config() -> anyhow::Result<()> {
    use sage_core::immonium::default_modified_ions;
    let fixture = serde_json::json!({
        "database": {"fasta": "tests/Q99536.fasta"},
        "mzml_paths": ["tests/LQSRPAAPPAPGPGQLTLR.mzML"],
        "precursor_tol": {"ppm": [-10, 10]},
        "fragment_tol": {"ppm": [-20, 20]}
    });
    let parse = |value: Option<serde_json::Value>| {
        let mut config = fixture.clone();
        if let Some(value) = value {
            config["immonium"] = value;
        }
        serde_json::from_value::<Input>(config)
    };
    let fragment_tol = sage_core::mass::Tolerance::Ppm(-20.0, 20.0);
    let resolve = |input: Input| {
        input
            .immonium
            .and_then(|config| config.resolve(fragment_tol))
    };

    assert!(parse(None)?.immonium.is_none());
    assert!(resolve(parse(Some(false.into()))?).is_none());
    let on = resolve(parse(Some(true.into()))?).unwrap();
    assert!(!on.rescore && on.residues);
    assert_eq!(on.modified, default_modified_ions());
    assert_eq!(on.tolerance, fragment_tol);

    let custom = parse(Some(serde_json::json!({
        "rescore": true,
        "modified": [{"name": "pH", "residue": "H", "modification": 79.96633, "mz": 190.0376}]
    })))?;
    custom.validate()?;
    let custom = resolve(custom).unwrap();
    assert!(custom.rescore);
    assert_eq!(custom.modified.len(), 1);

    for invalid in [
        serde_json::json!({"modified": [{"name": "", "residue": "Y", "modification": 79.97, "mz": 216.04}]}),
        serde_json::json!({"modified": [{"name": "x", "residue": "1", "modification": 79.97, "mz": 216.04}]}),
        serde_json::json!({"modified": [{"name": "x", "residue": "Y", "modification": 79.97, "mz": -1.0}]}),
        serde_json::json!({"tolerance": {"ppm": [5, 10]}}),
    ] {
        assert!(
            parse(Some(invalid.clone()))?.validate().is_err(),
            "{invalid}"
        );
    }
    // No formula strings or unknown keys.
    assert!(parse(Some(serde_json::json!({"loss": "H2O"}))).is_err());
    assert!(parse(Some(serde_json::json!({
        "modified": [{"name": "x", "residue": "Y", "mz": 216.04, "formula": "HPO3"}]
    })))
    .is_err());
    Ok(())
}
