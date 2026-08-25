use super::*;

#[test]
fn parse_modifications() {
    use InvalidModification::*;
    use ModificationSpecificity::*;
    assert_eq!("[".parse::<ModificationSpecificity>(), Ok(ProteinN(None)));
    assert_eq!(
        "[M".parse::<ModificationSpecificity>(),
        Ok(ProteinN(Some(b'M')))
    );
    assert_eq!(
        "]M".parse::<ModificationSpecificity>(),
        Ok(ProteinC(Some(b'M')))
    );
    assert_eq!("M".parse::<ModificationSpecificity>(), Ok(Residue(b'M')));
    assert_eq!(
        "Z".parse::<ModificationSpecificity>(),
        Err(InvalidResidue('Z'))
    );
}

#[test]
fn var_mod_entry_bare_mass() {
    let entry = VarModEntry::Mass(15.9949);
    assert_eq!(entry.mass(), 15.9949);
    assert_eq!(entry.max_count(), None);
}

#[test]
fn var_mod_entry_detailed_with_limit() {
    let entry = VarModEntry::Detailed(VariableModification {
        mass: 15.9949,
        max_count: Some(1),
        name: None,
        neutral_losses: vec![],
        neutral_loss_mode: NeutralLossMode::Optional,
        site_mode: SiteMode::Exhaustive,
        channel_offsets: Default::default(),
    });
    assert_eq!(entry.mass(), 15.9949);
    assert_eq!(entry.max_count(), Some(1));
}

#[test]
fn deserialize_var_mod_entries() {
    let entries: Vec<VarModEntry> =
        serde_json::from_str(r#"[15.9949, {"mass": 42.0106, "max_count": 1}, {"mass": 14.0157}]"#)
            .unwrap();

    assert_eq!(entries.len(), 3);
    assert!((entries[0].mass() - 15.9949).abs() < 1e-4);
    assert_eq!(entries[0].max_count(), None);
    assert!((entries[1].mass() - 42.0106).abs() < 1e-4);
    assert_eq!(entries[1].max_count(), Some(1));
    assert!((entries[2].mass() - 14.0157).abs() < 1e-4);
    assert_eq!(entries[2].max_count(), None);

    let serialized = serde_json::to_value(&entries).unwrap();
    assert!(serialized[0].is_number());
    assert!(serialized[1].is_object());
    assert_eq!(serialized[1]["max_count"], 1);
    assert!(serialized[2].is_object());
    assert!(serialized[2].get("max_count").is_none());

    let round_tripped: Vec<VarModEntry> = serde_json::from_value(serialized).unwrap();
    assert_eq!(entries, round_tripped);
}

#[test]
fn deserialize_named_neutral_loss_modifications() {
    let entry: VarModEntry = serde_json::from_str(
        r#"{
                "mass": 79.9663,
                "max_count": 2,
                "name": "Phospho",
                "neutral_losses": [97.9769],
                "neutral_loss_mode": "required",
                "site_mode": "both"
            }"#,
    )
    .unwrap();

    let VarModEntry::Detailed(entry) = &entry else {
        panic!("expected structured modification")
    };
    assert_eq!(entry.name.as_deref(), Some("Phospho"));
    assert_eq!(entry.neutral_losses, vec![97.9769]);
    assert_eq!(entry.neutral_loss_mode, NeutralLossMode::Required);
    assert_eq!(entry.site_mode, SiteMode::Both);

    let round_trip: VarModEntry =
        serde_json::from_value(serde_json::to_value(entry).unwrap()).unwrap();
    assert_eq!(round_trip, VarModEntry::Detailed(entry.clone()));
}

#[test]
fn static_mods_accept_numeric_and_structured_entries() {
    let numeric: StaticModEntry = serde_json::from_str("57.0215").unwrap();
    assert!((numeric.definition().mass - 57.0215).abs() < 1e-4);

    let detailed: StaticModEntry = serde_json::from_str(
        r#"{
                "mass": 57.0215,
                "name": "Carbamidomethyl",
                "neutral_losses": [18.0106]
            }"#,
    )
    .unwrap();
    let definition = detailed.definition();
    assert_eq!(definition.name.as_deref(), Some("Carbamidomethyl"));
    assert_eq!(&*definition.neutral_losses, &[18.0106]);
    assert_eq!(definition.neutral_loss_mode, NeutralLossMode::Optional);
}

#[test]
fn channel_offsets_round_trip_on_static_and_variable_modifications() {
    let static_entry: StaticModEntry = serde_json::from_value(serde_json::json!({
        "mass": 0.0,
        "name": "SILAC-K",
        "channel_offsets": {"light": 0.0, "heavy": 8.014199}
    }))
    .unwrap();
    assert_eq!(static_entry.channel_offsets()["heavy"], 8.014199);

    let variable_entry: VarModEntry = serde_json::from_value(serde_json::json!({
        "mass": 0.0,
        "name": "Optional-Lys8",
        "max_count": 2,
        "channel_offsets": {"light": 0.0, "heavy": 8.014199}
    }))
    .unwrap();
    assert_eq!(variable_entry.channel_offsets()["light"], 0.0);
    let serialized = serde_json::to_value(&variable_entry).unwrap();
    let heavy = serialized["channel_offsets"]["heavy"].as_f64().unwrap();
    assert!((heavy - 8.014199).abs() < 1e-6);
}

#[test]
fn channel_offsets_reject_invalid_names_and_masses() {
    for json in [
        r#"{"mass": 0.0, "channel_offsets": {" heavy": 8.0, "light": 0.0}}"#,
        r#"{"mass": 0.0, "channel_offsets": {"heavy": 1e999, "light": 0.0}}"#,
    ] {
        assert!(serde_json::from_str::<VarModEntry>(json).is_err());
    }
}

#[test]
fn reject_invalid_neutral_loss_configuration() {
    for json in [
        r#"{"mass": 79.9663, "name": ""}"#,
        r#"{"mass": 79.9663, "neutral_losses": [0.0]}"#,
        r#"{"mass": 79.9663, "neutral_losses": [-18.0]}"#,
        r#"{"mass": 79.9663, "neutral_loss_mode": "required"}"#,
        r#"{"mass": 79.9663, "neutral_loss_mode": "sometimes"}"#,
    ] {
        assert!(
            serde_json::from_str::<VarModEntry>(json).is_err(),
            "accepted invalid configuration: {json}"
        );
    }
}

#[test]
fn reject_unsupported_var_mod_shapes_and_fields() {
    assert!(serde_json::from_str::<Vec<VarModEntry>>(r#"[[42.0106, 1]]"#).is_err());
    assert!(
        serde_json::from_str::<Vec<VarModEntry>>(r#"[{"mass": 42.0106, "max_counts": 1}]"#)
            .is_err()
    );
}

#[test]
fn validate_var_mods_mixed() {
    use ModificationSpecificity::*;
    // Mix bare masses and detailed entries
    let mut raw = HashMap::new();
    raw.insert(
        "M".to_string(),
        vec![
            VarModEntry::Mass(15.9949),
            VarModEntry::Detailed(VariableModification {
                mass: 15.9949,
                max_count: Some(1),
                name: None,
                neutral_losses: vec![],
                neutral_loss_mode: NeutralLossMode::Optional,
                site_mode: SiteMode::Exhaustive,
                channel_offsets: Default::default(),
            }),
        ],
    );
    raw.insert(
        "C".to_string(),
        vec![VarModEntry::Detailed(VariableModification {
            mass: 57.0215,
            max_count: Some(2),
            name: None,
            neutral_losses: vec![],
            neutral_loss_mode: NeutralLossMode::Optional,
            site_mode: SiteMode::Exhaustive,
            channel_offsets: Default::default(),
        })],
    );
    let result = validate_var_mods(Some(raw));

    let m_entries = result.get(&Residue(b'M')).unwrap();
    assert_eq!(m_entries.len(), 2);
    assert!((m_entries[0].mass() - 15.9949).abs() < 1e-4);
    assert_eq!(m_entries[0].max_count(), None);
    assert!((m_entries[1].mass() - 15.9949).abs() < 1e-4);
    assert_eq!(m_entries[1].max_count(), Some(1));

    let c_entries = result.get(&Residue(b'C')).unwrap();
    assert_eq!(c_entries.len(), 1);
    assert!((c_entries[0].mass() - 57.0215).abs() < 1e-4);
    assert_eq!(c_entries[0].max_count(), Some(2));
}

#[test]
fn validate_var_mods_invalid_residue_skipped() {
    let mut raw = HashMap::new();
    raw.insert("Z".to_string(), vec![VarModEntry::Mass(15.9949)]);
    raw.insert("M".to_string(), vec![VarModEntry::Mass(15.9949)]);
    let result = validate_var_mods(Some(raw));
    // Z is invalid — only M should survive
    assert_eq!(result.len(), 1);
    assert!(result.contains_key(&ModificationSpecificity::Residue(b'M')));
}
