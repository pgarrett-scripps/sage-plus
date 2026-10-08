use super::{list, preset, write, PRESETS};
use crate::input::Input;

fn scratch(label: &str) -> std::path::PathBuf {
    let root = std::env::temp_dir().join(format!(
        "sage-presets-{label}-{}-{}",
        std::process::id(),
        std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .unwrap()
            .as_nanos()
    ));
    std::fs::create_dir_all(&root).unwrap();
    root
}

fn parse(name: &str) -> serde_json::Value {
    serde_json::from_str(preset(name).unwrap()).unwrap()
}

#[test]
fn every_preset_passes_validation() {
    for (name, _, json) in PRESETS {
        let input: Input =
            serde_json::from_str(json).unwrap_or_else(|error| panic!("{name}: {error}"));
        input
            .validate()
            .unwrap_or_else(|error| panic!("{name}: {error}"));
        assert!(list().contains(name));
    }
}

/// `full` names every setting, and each value is the one Sage uses when the
/// setting is left out: the search it builds equals the `minimal` search.
#[test]
fn full_preset_lists_every_setting_at_its_default() {
    let full = parse("full");
    let schema: serde_json::Value =
        serde_json::from_str(&crate::config_schema::generate_config_schema()).unwrap();
    let deprecated = [
        "min_free_memory_gb",
        "prefilter_chunk_size",
        "prefilter_low_memory",
    ];
    let missing = |properties: &serde_json::Value, config: &serde_json::Value| {
        properties
            .as_object()
            .unwrap()
            .keys()
            .filter(|key| !deprecated.contains(&key.as_str()) && config.get(key).is_none())
            .cloned()
            .collect::<Vec<_>>()
    };
    assert_eq!(missing(&schema["properties"], &full), Vec::<String>::new());
    let database = &schema["properties"]["database"]["properties"];
    assert_eq!(missing(database, &full["database"]), Vec::<String>::new());
    let enzyme = &schema["$defs"]["EnzymeBuilder"]["properties"];
    assert_eq!(
        missing(enzyme, &full["database"]["enzyme"]),
        Vec::<String>::new()
    );

    let output = scratch("full");
    let build = |mut config: serde_json::Value| {
        // Building resolves the input paths, so they must exist.
        let workspace = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("../..");
        config["database"]["fasta"] = workspace
            .join("tests/Q99536.fasta")
            .display()
            .to_string()
            .into();
        config["mzml_paths"] =
            serde_json::json!([workspace.join("tests/LQSRPAAPPAPGPGQLTLR.mzML")]);
        config["output_directory"] = output.display().to_string().into();
        let input: Input = serde_json::from_value(config).unwrap();
        let mut search = serde_json::to_value(input.build().unwrap()).unwrap();
        search.as_object_mut().unwrap().remove("batch_size");
        // The search records the enzyme as written; compare effective values.
        let enzyme = &mut search["database"]["enzyme"];
        let builder: sage_core::database::EnzymeBuilder =
            serde_json::from_value(enzyme.clone()).unwrap();
        enzyme["min_len"] = builder.min_len.unwrap_or(5).into();
        enzyme["max_len"] = builder.effective_max_len().into();
        enzyme["c_terminal"] = builder.c_terminal.unwrap_or(true).into();
        search
    };
    assert_eq!(build(full), build(parse("minimal")));
    std::fs::remove_dir_all(output).unwrap();
}

#[test]
fn write_refuses_to_replace_a_file_without_overwrite() {
    let directory = scratch("write");
    let path = directory.join("config.json");
    let path = path.to_str().unwrap();
    write("minimal", Some(path), false).unwrap();
    assert_eq!(
        std::fs::read_to_string(path).unwrap(),
        preset("minimal").unwrap()
    );

    let error = write("phospho", Some(path), false).unwrap_err().to_string();
    assert!(error.contains("--overwrite"), "{error}");
    assert_eq!(
        std::fs::read_to_string(path).unwrap(),
        preset("minimal").unwrap()
    );

    write("phospho", Some(path), true).unwrap();
    assert_eq!(
        std::fs::read_to_string(path).unwrap(),
        preset("phospho").unwrap()
    );

    let error = write("tryptic", None, false).unwrap_err().to_string();
    assert!(error.contains("trypsin-hcd"), "{error}");
    std::fs::remove_dir_all(directory).unwrap();
}
