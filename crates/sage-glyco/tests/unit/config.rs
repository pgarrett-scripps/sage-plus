use super::*;

#[test]
fn defaults_and_validation() {
    let config = GlycoConfig::from_value(serde_json::json!({})).unwrap();
    assert_eq!(config, GlycoConfig::default());
    assert!(GlycoConfig::from_value(serde_json::json!({"isotope_errors": [2, 1]})).is_err());
    assert!(GlycoConfig::from_value(serde_json::json!({"glycans": ["Kdn(1)"]})).is_err());
    assert!(GlycoConfig::from_value(serde_json::json!({"unknown": 1})).is_err());
    assert!(GlycoConfig::from_value(serde_json::json!({"sequon": "N"})).is_err());
    assert!(GlycoConfig::from_value(serde_json::json!({"site_fragment_forms": ["bare"]})).is_err());
    assert!(
        GlycoConfig::from_value(serde_json::json!({"site_fragment_forms": ["hexnac", "x"]}))
            .is_err()
    );
    let forms = GlycoConfig::from_value(serde_json::json!({
        "site_fragment_forms": ["hexnac", "bare", "hexnac2"]
    }))
    .unwrap();
    assert_eq!(forms.site_losses(), vec![-HEXNAC as f32, HEXNAC as f32]);
    assert_eq!(GlycoConfig::default().site_losses(), vec![HEXNAC as f32]);
    assert!(GlycoConfig::from_value(serde_json::json!({"twin_feature": "decoy"})).is_err());
    for mode in TWIN_FEATURES {
        assert!(GlycoConfig::from_value(serde_json::json!({"twin_feature": mode})).is_ok());
    }
}

#[test]
fn library_merges_sources_and_requires_hexnac() {
    let config = GlycoConfig {
        glycans: vec!["HexNAc(4)Hex(5)".into()],
        high_mannose: true,
        ..Default::default()
    };
    let files = vec!["# comment\nN(2)H(9)\nHex(3)\nnot a glycan\n".to_string()];
    let library = config.library(&files).unwrap();
    // HexNAc(4)Hex(5), HexNAc(2)Hex(3..=20) (which includes N(2)H(9));
    // Hex(3) has no HexNAc and is dropped.
    assert_eq!(library.len(), 1 + 18);
    let default = GlycoConfig::default().library(&[]).unwrap();
    assert!(default.len() > 100);
}

#[test]
fn hexnac_offset_is_a_valid_modification() {
    let builder: sage_core::database::Builder = serde_json::from_value(serde_json::json!({
        "variable_mods": {"HexNAc": GlycoConfig::default().hexnac_offset()}
    }))
    .unwrap();
    builder.validate_modification_keys().unwrap();
    let parameters = builder.make_parameters();
    let offsets = parameters.mass_offset_modifications();
    assert_eq!(offsets.len(), 1);
    assert!((offsets[0].mass() as f64 - HEXNAC).abs() < 1e-3);
}
