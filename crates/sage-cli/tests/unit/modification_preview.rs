use super::*;

#[test]
fn preview_h3k9_and_protein_terminal_exception() {
    let config = r#"{"database":{"variable_mods":{"Acetyl":{"mass":42.0106,"max_count":1,
        "sites":["first_residue:K","internal_residue:K","protein_last:K"]}}}}"#;
    let preview = preview(config, "KSTGGKAPR", "internal", 100).unwrap();
    assert_eq!(preview["variants"].as_array().unwrap().len(), 3);
    let terminal = super::preview(config, "KAK", "cterm", 100).unwrap();
    assert_eq!(terminal["variants"].as_array().unwrap().len(), 3);
    assert_eq!(terminal["truncated"], false);
}

#[test]
fn preview_is_bounded_and_validates_input() {
    let config = r#"{"database":{"variable_mods":{"Acetyl":{"mass":42.0,"sites":["internal_residue:K"]}},"max_variable_mods":4}}"#;
    let result = preview(config, "AKKKKA", "internal", 2).unwrap();
    assert_eq!(result["variants"].as_array().unwrap().len(), 2);
    assert_eq!(result["truncated"], true);
    assert!(preview(config, "BAD!", "internal", 2).is_err());
    assert!(preview(config, "AKA", "internal", 0).is_err());
    assert!(preview(
        r#"{"database":{"variable_mods":{"~":[42]}}}"#,
        "AKA",
        "internal",
        2
    )
    .is_err());
}

#[test]
fn preview_includes_search_time_offsets() {
    let result = preview(r#"{"database":{"variable_mods":{"Acetyl":{"mass":42.0,"max_count":1,"search_mode":"mass_offset","sites":["internal_residue:K"]}}}}"#,
        "KAKAK", "internal", 100).unwrap();
    assert_eq!(result["variants"].as_array().unwrap().len(), 2);
    assert_eq!(
        result["rules"][0]["eligible_sites"],
        json!([{"position":3}])
    );
}

fn error(result: anyhow::Result<Value>) -> String {
    format!("{:#}", result.unwrap_err())
}

#[test]
fn preview_rejects_invalid_requests_with_clear_errors() {
    let config = r#"{"database":{"variable_mods":{"Oxidation":{"mass":15.994915,"sites":["M"]}}}}"#;
    assert_eq!(
        error(preview_with_flanks(
            config,
            "PEPMK",
            "internal",
            10,
            None,
            ("m", "")
        )),
        "preview flanking residues must be uppercase one-letter codes"
    );
    assert_eq!(
        error(preview(config, "PEPMK", "internal", 10_001)),
        "preview limit must be between 1 and 10000"
    );
    assert!(error(preview("{not json", "PEPMK", "internal", 10))
        .starts_with("invalid JSON configuration"));
    assert_eq!(
        error(preview("{}", "PEPMK", "internal", 10)),
        "`database` must be configured"
    );
    assert_eq!(
        error(preview(config, "PEPMK", "middle", 10)),
        "unknown peptide position `middle`"
    );
    assert_eq!(
        error(preview_with_context(
            config,
            "PEPMK",
            "internal",
            10,
            Some(("P1", 0))
        )),
        "preview-start is one-based"
    );
    assert_eq!(
        error(preview_with_flanks(
            config,
            "PEPMK",
            "internal",
            10,
            Some(("P1", 2)),
            ("MR", "")
        )),
        "--preview-before is longer than the residues before --preview-start"
    );
    let library_mode = r#"{"database":{"variable_mods":{"Acetyl":{"mass":42.010565,"sites":["K"],"max_count":1,"site_mode":"both"}}}}"#;
    assert_eq!(
        error(preview(library_mode, "PEPKK", "internal", 10)),
        "preview requires a PTM library for non-exhaustive site_mode"
    );
}

#[test]
fn preview_reports_static_rules_offsets_and_effective_counts() {
    let config = r#"{"database":{
        "static_mods":{"Carbamidomethyl":{"mass":57.021464,"sites":["C"]}},
        "variable_mods":{
            "Oxidation":{"mass":15.994915,"sites":["M"],"max_count":2},
            "Phospho":{"mass":79.966331,"sites":["S"],"search_mode":"mass_offset"}
        },
        "max_combinations":50
    }}"#;
    let result = preview(config, "MCSMK", "full", 100).unwrap();
    let rules = result["rules"].as_array().unwrap();
    assert_eq!(rules.len(), 3);
    let rule = |name: &str| {
        rules
            .iter()
            .find(|rule| rule["name"] == name)
            .unwrap_or_else(|| panic!("missing {name}"))
    };
    let static_rule = rule("Carbamidomethyl");
    assert_eq!(static_rule["kind"], "static");
    assert_eq!(static_rule["key"], "C");
    assert_eq!(static_rule["eligible_sites"], json!([{"position": 2}]));
    let oxidation = rule("Oxidation");
    assert_eq!(oxidation["kind"], "variable");
    assert_eq!(oxidation["max_count"], 2);
    assert_eq!(oxidation["effective_max_count"], 2);
    assert_eq!(
        oxidation["eligible_sites"],
        json!([{"position": 1}, {"position": 4}])
    );
    let phospho = rule("Phospho");
    assert_eq!(phospho["search_mode"], "mass_offset");
    assert_eq!(phospho["max_count"], Value::Null);
    assert_eq!(phospho["effective_max_count"], 1);
    assert_eq!(phospho["eligible_sites"], json!([{"position": 3}]));

    assert_eq!(
        result["limits"],
        json!({"max_variable_mods": 2, "max_total_variable_mods": 2,
               "max_combinations": 50, "preview_limit": 100})
    );
    // Four indexed oxidation variants, each also previewed with phospho on S.
    let variants = result["variants"].as_array().unwrap();
    assert_eq!(variants.len(), 8);
    assert_eq!(result["truncated"], false);
    for variant in variants {
        assert!(variant["modifications"]
            .as_array()
            .unwrap()
            .iter()
            .any(|applied| applied["name"] == "Carbamidomethyl"
                && applied["site"] == json!({"position": 2})));
    }
    let phosphorylated = variants
        .iter()
        .filter(|variant| {
            variant["modifications"]
                .as_array()
                .unwrap()
                .iter()
                .any(|applied| applied["name"] == "Phospho")
        })
        .count();
    assert_eq!(phosphorylated, 4);
}

#[test]
fn library_preview_preserves_attachment_and_protein_coordinates() {
    let path = std::env::temp_dir().join(format!("sage-beta6-preview-{}.tsv", std::process::id()));
    std::fs::write(
        &path,
        "protein\tposition\tresidue\tmodification\tattachment\nP1\t9\tK\tAcetyl\tpeptide_n_term\n",
    )
    .unwrap();
    let config = json!({"database":{
        "variable_mods":{"Acetyl":{"mass":42.010565,"sites":["first_residue:K","peptide_n_term:K"],"max_count":1,"site_mode":"library"}},
        "ptm_library":{"path":path}
    }}).to_string();
    assert!(preview(&config, "KAAAK", "internal", 100).is_err());
    let value = preview_with_context(&config, "KAAAK", "internal", 100, Some(("P1", 9))).unwrap();
    assert_eq!(value["variants"].as_array().unwrap().len(), 2);
    let modified = value["variants"]
        .as_array()
        .unwrap()
        .iter()
        .find(|v| !v["modifications"].as_array().unwrap().is_empty())
        .unwrap();
    assert_eq!(modified["modifications"][0]["site"]["terminal_group"], "N");
    let other = preview_with_context(&config, "KAAAK", "internal", 100, Some(("P1", 10))).unwrap();
    assert_eq!(other["variants"].as_array().unwrap().len(), 1);
    std::fs::remove_file(path).unwrap();
}
