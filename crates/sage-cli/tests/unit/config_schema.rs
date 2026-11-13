use super::*;

#[test]
fn committed_config_schema_matches_rust_types() {
    let committed = include_str!("../../../../schemas/config.schema.json");
    assert_eq!(committed, generate_config_schema());
}

#[test]
fn generated_config_schema_has_expected_contract() {
    let schema: serde_json::Value = serde_json::from_str(&generate_config_schema()).unwrap();
    assert_eq!(schema["title"], "Sage search configuration");
    assert_eq!(schema["additionalProperties"], false);
    assert!(schema["properties"]["quant"]["anyOf"].is_array());
    assert_eq!(
        schema["properties"]["database"]["additionalProperties"],
        false
    );
}
