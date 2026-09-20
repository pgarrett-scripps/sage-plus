use anyhow::{bail, ensure, Context};
use sage_core::database::Builder;
use sage_core::modification::ModificationSpecificity;
use serde_json::{json, Map, Value};
use std::collections::{BTreeMap, BTreeSet, HashSet};

/// Convert declarations without changing the rest of the search configuration.
pub fn migrate(config: &str) -> anyhow::Result<Value> {
    let mut config: Value = serde_json::from_str(config)?;
    let database = config.get_mut("database").context("database is required")?;
    let builder: Builder = serde_json::from_value(database.clone())?;
    let parameters = builder.make_parameters();
    parameters
        .validate_ptm_library(&Default::default())
        .map_err(anyhow::Error::msg)?;
    let mut reserved = HashSet::new();
    for section in ["static_mods", "variable_mods"] {
        if let Some(entries) = database.get(section).and_then(Value::as_object) {
            for (key, value) in entries {
                if value.get("sites").is_some() {
                    reserved.insert(key.clone());
                }
                let items = value
                    .as_array()
                    .map(Vec::as_slice)
                    .unwrap_or_else(|| std::slice::from_ref(value));
                for item in items {
                    if let Some(name) = item.get("name").and_then(Value::as_str) {
                        reserved.insert(name.to_string());
                    }
                }
            }
        }
    }
    for section in ["static_mods", "variable_mods"] {
        let Some(entries) = database.get(section).and_then(Value::as_object) else {
            continue;
        };
        if entries.values().any(|entry| entry.get("sites").is_some()) {
            continue;
        }
        let mut definitions: BTreeMap<String, (Value, BTreeSet<String>)> = BTreeMap::new();
        let mut next = 1;
        for (key, value) in entries {
            let site = key
                .parse::<ModificationSpecificity>()
                .map_err(|_| anyhow::anyhow!("invalid legacy site `{key}`"))?
                .explicit_name();
            let items = value
                .as_array()
                .map(Vec::as_slice)
                .unwrap_or_else(|| std::slice::from_ref(value));
            for item in items {
                let mut definition = if item.is_number() {
                    json!({"mass": item})
                } else {
                    item.clone()
                };
                let definition_object = definition
                    .as_object_mut()
                    .context("modification must be a mass or object")?;
                let id = match definition_object.remove("name") {
                    Some(Value::String(id)) => id,
                    Some(Value::Null) | None => loop {
                        let id = format!("legacy_{section}_{next}");
                        next += 1;
                        if reserved.insert(id.clone()) {
                            break id;
                        }
                    },
                    _ => bail!("invalid modification name"),
                };
                if let Some((previous, sites)) = definitions.get_mut(&id) {
                    ensure!(*previous == definition, "modification `{id}` has inconsistent definitions. Resolve them before migration");
                    sites.insert(site.clone());
                } else {
                    definitions.insert(id, (definition, BTreeSet::from([site.clone()])));
                }
            }
        }
        let converted = definitions
            .into_iter()
            .map(|(id, (mut definition, sites))| {
                definition["sites"] = json!(sites);
                (id, definition)
            })
            .collect::<Map<_, _>>();
        database[section] = Value::Object(converted);
    }
    let migrated: Builder = serde_json::from_value(database.clone())?;
    migrated
        .make_parameters()
        .validate_ptm_library(&Default::default())
        .map_err(anyhow::Error::msg)?;
    Ok(config)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn migration_groups_named_sites_but_keeps_unnamed_equal_masses_separate() {
        let input = json!({"database":{"variable_mods":{
            "^K":[{"mass":42,"name":"Acetyl","max_count":2}],
            "~K":[{"mass":42,"name":"Acetyl","max_count":2}],
            "S":[80],"T":[80]
        }},"other":"preserved"});
        let converted = migrate(&input.to_string()).unwrap();
        let mods = &converted["database"]["variable_mods"];
        assert_eq!(mods.as_object().unwrap().len(), 3);
        assert_eq!(
            mods["Acetyl"]["sites"],
            json!(["first_residue:K", "internal_residue:K"])
        );
        assert_eq!(converted["other"], "preserved");
        assert_eq!(migrate(&converted.to_string()).unwrap(), converted);
    }

    #[test]
    fn migration_preserves_terminal_group_and_residue_distinction() {
        let converted =
            migrate(r#"{"database":{"static_mods":{"^":42,"^K":14,"$":1,"$R":2}}}"#).unwrap();
        let sites = converted["database"]["static_mods"]
            .as_object()
            .unwrap()
            .values()
            .map(|entry| entry["sites"][0].as_str().unwrap())
            .collect::<BTreeSet<_>>();
        assert_eq!(
            sites,
            BTreeSet::from([
                "peptide_n_term",
                "first_residue:K",
                "peptide_c_term",
                "last_residue:R"
            ])
        );
    }
}
