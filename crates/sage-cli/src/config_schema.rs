pub const CONFIG_SCHEMA: &str = include_str!("../../../schemas/config.schema.json");

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn committed_config_schema_is_valid_json() {
        let schema: serde_json::Value = serde_json::from_str(CONFIG_SCHEMA).unwrap();
        assert_eq!(schema["title"], "Sage search configuration");
        assert!(schema["properties"]["quant"]["$ref"].is_string());
        assert_eq!(schema["$defs"]["lfq"]["properties"]["mbr"]["default"], true);
        assert_eq!(
            schema["$defs"]["mobilityModel"]["properties"]["enabled"]["default"],
            true
        );
    }
}
