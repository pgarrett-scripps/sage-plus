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

    let explicit = "^([ACDEFGHIKLMNPQRSTVWYUO]|(first_residue|internal_residue|last_residue|protein_first|protein_last):[ACDEFGHIKLMNPQRSTVWYUO]|(peptide_n_term|peptide_c_term|protein_n_term|protein_c_term)(:[ACDEFGHIKLMNPQRSTVWYUO])?|motif:[^\\s]+)(?![\\s\\S])";
    for name in ["NamedStaticModification", "NamedVariableModification"] {
        value["$defs"][name]["properties"]["sites"]["minItems"] = 1.into();
        value["$defs"][name]["properties"]["sites"]["items"]["pattern"] = explicit.into();
        value["$defs"][name]["additionalProperties"] = false.into();
    }

    value["$defs"]["NamedVariableModification"]["properties"]["max_count"]["minimum"] = 1.into();
    // Symbol-keyed maps accept upstream Sage syntax only: residue or terminal
    // symbol keys mapped to masses (variable mods: a mass or a mass array).
    let legacy = "^([ACDEFGHIKLMNPQRSTVWYUO]|[\\^$\\[\\]][ACDEFGHIKLMNPQRSTVWYUO]?)(?![\\s\\S])";
    let mass = serde_json::json!({"format": "float", "type": "number"});
    let legacy_values = [
        ("StaticModConfig", mass.clone()),
        (
            "VariableModConfig",
            serde_json::json!({"anyOf": [mass.clone(), {"items": mass, "type": "array"}]}),
        ),
    ];
    for (name, values) in legacy_values {
        value["$defs"][name]["anyOf"][1]["propertyNames"] = serde_json::json!({"pattern": legacy});
        value["$defs"][name]["anyOf"][1]["additionalProperties"] = values;
        value["$defs"][name]["anyOf"][0]["propertyNames"] = serde_json::json!({"minLength":1});
    }

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
