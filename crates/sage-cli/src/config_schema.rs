use crate::input::Input;
use schemars::generate::SchemaSettings;

const SCHEMA_ID: &str =
    "https://raw.githubusercontent.com/pgarrett-scripps/sage-plus/main/schemas/config.schema.json";

pub fn generate_config_schema() -> String {
    let generator = SchemaSettings::draft2020_12().into_generator();
    let schema = generator.into_root_schema_for::<Input>();
    let mut value = serde_json::to_value(schema).expect("configuration schema is serializable");
    value
        .as_object_mut()
        .expect("root configuration schema is an object")
        .insert("$id".into(), SCHEMA_ID.into());
    value
        .pointer_mut("/properties/database/properties/prefilter_low_memory")
        .and_then(serde_json::Value::as_object_mut)
        .expect("deprecated compatibility field is present")
        .insert("deprecated".into(), true.into());

    let mut json =
        serde_json::to_string_pretty(&value).expect("configuration schema serializes as JSON");
    json.push('\n');
    json
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn committed_config_schema_matches_rust_types() {
        let committed = include_str!("../../../schemas/config.schema.json");
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
}
