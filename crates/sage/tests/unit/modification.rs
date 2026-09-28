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
        search_mode: SearchMode::Database,
        mass: 15.9949,
        max_count: Some(1),
        max_total_count: None,
        name: None,
        neutral_losses: vec![],
        neutral_loss_mode: NeutralLossMode::Optional,
        site_mode: SiteMode::Exhaustive,
        channel_offsets: Default::default(),
        immonium_ions: Default::default(),
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
                search_mode: SearchMode::Database,
                mass: 15.9949,
                max_count: Some(1),
                max_total_count: None,
                name: None,
                neutral_losses: vec![],
                neutral_loss_mode: NeutralLossMode::Optional,
                site_mode: SiteMode::Exhaustive,
                channel_offsets: Default::default(),
                immonium_ions: Default::default(),
            }),
        ],
    );
    raw.insert(
        "C".to_string(),
        vec![VarModEntry::Detailed(VariableModification {
            search_mode: SearchMode::Database,
            mass: 57.0215,
            max_count: Some(2),
            max_total_count: None,
            name: None,
            neutral_losses: vec![],
            neutral_loss_mode: NeutralLossMode::Optional,
            site_mode: SiteMode::Exhaustive,
            channel_offsets: Default::default(),
            immonium_ions: Default::default(),
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
#[should_panic(expected = "invalid modification key `Z`")]
fn validate_var_mods_invalid_residue_rejected() {
    let mut raw = HashMap::new();
    raw.insert("Z".to_string(), vec![VarModEntry::Mass(15.9949)]);
    raw.insert("M".to_string(), vec![VarModEntry::Mass(15.9949)]);
    validate_var_mods(Some(raw));
}

#[test]
fn search_mode_defaults_to_database_and_accepts_mass_offset() {
    let entry: VarModEntry =
        serde_json::from_str(r#"{"mass": 79.966331, "name": "Phospho"}"#).unwrap();
    assert_eq!(entry.search_mode(), SearchMode::Database);
    let entry: VarModEntry = serde_json::from_str(
        r#"{"mass": 79.966331, "name": "Phospho", "search_mode": "mass_offset"}"#,
    )
    .unwrap();
    assert_eq!(entry.search_mode(), SearchMode::MassOffset);
    assert_eq!(VarModEntry::Mass(15.99).search_mode(), SearchMode::Database);
    let serialized = serde_json::to_value(&entry).unwrap();
    assert_eq!(serialized["search_mode"], "mass_offset");
}

#[test]
fn mass_offset_rejects_labels_and_zero_mass() {
    for invalid in [
        r#"{"mass": 0.0, "search_mode": "mass_offset"}"#,
        r#"{"mass": 4.0, "search_mode": "mass_offset", "channel_offsets": {"L": 0.0, "H": 4.0}}"#,
        r#"{"mass": 4.0, "search_mode": "unknown"}"#,
    ] {
        assert!(
            serde_json::from_str::<VarModEntry>(invalid).is_err(),
            "{invalid}"
        );
    }
}

#[test]
fn positional_keys_round_trip_and_reject_malformed_keys() {
    for key in [
        "K", "^", "$", "[", "]", "^K", "$K", "[K", "]K", "~K", "~O", "~U",
    ] {
        let specificity: ModificationSpecificity = key.parse().unwrap();
        assert_eq!(specificity.to_string(), key);
        assert_eq!(serde_json::to_value(specificity).unwrap(), key);
    }
    for key in [
        "", "~", "~Z", "~k", "~KK", "KK", "K~", "^Z", "$?", "[!", "]1", "é", "~é", " K", "K ",
    ] {
        assert!(
            key.parse::<ModificationSpecificity>().is_err(),
            "accepted {key}"
        );
        for field in ["static_mods", "variable_mods"] {
            let value = if field == "static_mods" {
                serde_json::json!(42)
            } else {
                serde_json::json!([42])
            };
            let config = serde_json::json!({field: {key: value}});
            assert!(
                serde_json::from_value::<crate::database::Builder>(config).is_err(),
                "accepted {field}.{key}"
            );
        }
    }
    assert!(!ModificationSpecificity::is_internal(
        usize::MAX,
        usize::MAX
    ));
}

#[test]
fn explicit_sites_round_trip_and_distinguish_attachment() {
    use crate::enzyme::Position;
    use crate::peptide::Site;
    for (text, sites) in [
        (
            "K",
            vec![Site::Sequence(0), Site::Sequence(2), Site::Sequence(4)],
        ),
        ("first_residue:K", vec![Site::Sequence(0)]),
        ("last_residue:K", vec![Site::Sequence(4)]),
        ("internal_residue:K", vec![Site::Sequence(2)]),
        ("peptide_n_term:K", vec![Site::Nterm]),
        ("peptide_c_term:K", vec![Site::Cterm]),
        ("protein_n_term:K", vec![Site::Nterm]),
        ("protein_c_term:K", vec![Site::Cterm]),
        ("protein_first:K", vec![Site::Sequence(0)]),
        ("protein_last:K", vec![Site::Sequence(4)]),
    ] {
        let specificity: ModificationSpecificity = text.parse().unwrap();
        assert_eq!(specificity.explicit_name(), text);
        assert_eq!(specificity.sites(b"KAKAK", Position::Full), sites, "{text}");
        assert!(specificity.sites(b"", Position::Full).is_empty());
    }
    assert!("peptide_n_term:K"
        .parse::<ModificationSpecificity>()
        .unwrap()
        .sites(b"AK", Position::Internal)
        .is_empty());
    assert!("protein_n_term"
        .parse::<ModificationSpecificity>()
        .unwrap()
        .sites(b"KA", Position::Internal)
        .is_empty());
    for sequence in [b"K".as_slice(), b"KK".as_slice()] {
        assert!("internal_residue:K"
            .parse::<ModificationSpecificity>()
            .unwrap()
            .sites(sequence, Position::Full)
            .is_empty());
    }
}

#[test]
fn every_specificity_round_trips_through_explicit_name_and_display() {
    use ModificationSpecificity::*;
    for specificity in [
        Residue(b'K'),
        Internal(b'K'),
        PeptideN(None),
        PeptideC(None),
        ProteinN(None),
        ProteinC(None),
        PeptideN(Some(b'Q')),
        PeptideC(Some(b'K')),
        ProteinN(Some(b'M')),
        ProteinC(Some(b'R')),
        PeptideNTerm(b'K'),
        PeptideCTerm(b'K'),
        ProteinNTerm(b'M'),
        ProteinCTerm(b'R'),
    ] {
        let explicit = specificity.explicit_name();
        assert_eq!(explicit.parse(), Ok(specificity), "{explicit}");
        let display = specificity.to_string();
        assert_eq!(display.parse(), Ok(specificity), "{display}");
        assert!(!specificity.is_motif());
    }
    // Terminal groups gated on a residue have no symbol form, so they
    // display with their explicit name.
    assert_eq!(PeptideNTerm(b'K').to_string(), "peptide_n_term:K");
    assert_eq!(ProteinCTerm(b'R').to_string(), "protein_c_term:R");
    assert_eq!(PeptideN(None).to_string(), "^");
    assert_eq!(PeptideN(None).explicit_name(), "peptide_n_term");
}

#[test]
fn specificity_parse_errors_identify_the_problem() {
    use InvalidModification::*;
    let parse = |text: &str| text.parse::<ModificationSpecificity>();
    assert_eq!(parse(""), Err(Empty));
    assert_eq!(
        parse("first_residue:KK"),
        Err(TooLong("first_residue:KK".into()))
    );
    assert_eq!(parse("peptide_n_term:Z"), Err(InvalidResidue('Z')));
    // A residue-gated name without its residue is not a site.
    assert_eq!(parse("first_residue"), Err(TooLong("first_residue".into())));
    assert_eq!(parse("^Z"), Err(InvalidResidue('Z')));
    assert_eq!(parse("~"), Err(InvalidResidue('~')));
    assert!(matches!(parse("motif:"), Err(Motif(_))));
}

#[test]
fn residue_anchored_rules_follow_peptide_and_protein_position() {
    use crate::enzyme::Position;
    use crate::peptide::Site;
    use ModificationSpecificity::*;
    // A single residue is both the first and the last residue.
    assert_eq!(
        PeptideN(Some(b'K')).sites(b"K", Position::Internal),
        vec![Site::Sequence(0)]
    );
    assert_eq!(
        PeptideC(Some(b'K')).sites(b"K", Position::Internal),
        vec![Site::Sequence(0)]
    );
    // Residue-gated rules do not fire on a different residue.
    assert!(PeptideN(Some(b'K')).sites(b"AK", Position::Full).is_empty());
    assert!(PeptideC(Some(b'K')).sites(b"KA", Position::Full).is_empty());
    assert!(PeptideNTerm(b'K').sites(b"AK", Position::Full).is_empty());
    assert!(PeptideCTerm(b'K').sites(b"KA", Position::Full).is_empty());
    // Protein anchors need the matching protein end.
    for (position, n, c) in [
        (Position::Internal, false, false),
        (Position::Nterm, true, false),
        (Position::Cterm, false, true),
        (Position::Full, true, true),
    ] {
        let sequence = b"MAK";
        assert_eq!(!ProteinN(None).sites(sequence, position).is_empty(), n);
        assert_eq!(!ProteinC(None).sites(sequence, position).is_empty(), c);
        assert_eq!(
            !ProteinN(Some(b'M')).sites(sequence, position).is_empty(),
            n
        );
        assert_eq!(
            !ProteinC(Some(b'K')).sites(sequence, position).is_empty(),
            c
        );
        assert_eq!(!ProteinNTerm(b'M').sites(sequence, position).is_empty(), n);
        assert_eq!(!ProteinCTerm(b'K').sites(sequence, position).is_empty(), c);
        // Peptide-terminal rules never depend on the protein position.
        assert_eq!(PeptideN(None).sites(sequence, position), vec![Site::Nterm]);
        assert_eq!(PeptideC(None).sites(sequence, position), vec![Site::Cterm]);
    }
    assert!(ProteinN(Some(b'A'))
        .sites(b"MAK", Position::Full)
        .is_empty());
    assert!(ProteinCTerm(b'A').sites(b"MAK", Position::Full).is_empty());
}

#[test]
fn overlap_detects_rules_that_can_share_an_attachment() {
    use ModificationSpecificity::*;
    let overlaps = |left: ModificationSpecificity, right: ModificationSpecificity| {
        let forward = left.overlaps(right);
        assert_eq!(forward, right.overlaps(left), "{left} vs {right}");
        forward
    };
    // Same residue attachment.
    assert!(overlaps(Residue(b'K'), PeptideN(Some(b'K'))));
    assert!(overlaps(Residue(b'K'), Internal(b'K')));
    assert!(overlaps(Residue(b'K'), ProteinC(Some(b'K'))));
    // A one-residue peptide makes the first and last residue the same site.
    assert!(overlaps(PeptideN(Some(b'K')), PeptideC(Some(b'K'))));
    // Terminal groups.
    assert!(overlaps(PeptideN(None), ProteinN(None)));
    assert!(overlaps(PeptideN(None), PeptideNTerm(b'K')));
    assert!(overlaps(ProteinC(None), PeptideCTerm(b'R')));
    // Different residues or different attachments never overlap.
    assert!(!overlaps(Residue(b'K'), Residue(b'R')));
    assert!(!overlaps(Internal(b'K'), PeptideN(Some(b'K'))));
    assert!(!overlaps(Internal(b'K'), PeptideC(Some(b'K'))));
    assert!(!overlaps(PeptideNTerm(b'K'), PeptideNTerm(b'R')));
    assert!(!overlaps(PeptideN(None), PeptideC(None)));
    assert!(!overlaps(Residue(b'K'), PeptideNTerm(b'K')));
}

#[test]
fn named_definitions_reject_ambiguous_or_invalid_configuration() {
    use crate::database::Builder;
    for value in [
        serde_json::json!({"Acetyl":{"mass":42,"sites":[]}}),
        serde_json::json!({"Acetyl":{"mass":42,"sites":["^K"]}}),
        serde_json::json!({"Acetyl":{"mass":42,"sites":["peptide_n_term:KK"]}}),
        serde_json::json!({"Acetyl":{"mass":42,"sites":["first_residue"]}}),
        serde_json::json!({"Acetyl":{"mass":42,"sites":["K"],"name":"Other"}}),
        serde_json::json!({"Acetyl":{"mass":42,"sites":["K"],"max_count":0}}),
        serde_json::json!({"Acetyl":{"mass":42,"sites":["K"]},"M":[16]}),
        serde_json::json!({"Acetyl":{"mass":42,"sites":["K"],"unknown":true}}),
    ] {
        assert!(
            serde_json::from_value::<Builder>(serde_json::json!({"variable_mods":value})).is_err(),
            "{value}"
        );
    }
}

#[test]
fn named_definition_errors_explain_the_fix() {
    let error = |section: &str, value: serde_json::Value| {
        serde_json::from_value::<crate::database::Builder>(serde_json::json!({ section: value }))
            .err()
            .expect("configuration must be rejected")
            .to_string()
    };
    let cases = [
        (
            serde_json::json!({" Acetyl": {"mass": 42, "sites": ["K"]}}),
            "modification IDs must be nonempty and have no surrounding whitespace",
        ),
        (
            serde_json::json!({"Acetyl": {"mass": 42, "sites": ["K"]}, "Methyl": {"mass": 14}}),
            "modification `Methyl` requires `sites`",
        ),
        (
            serde_json::json!({"Acetyl": {"mass": 42, "sites": ["K"]}, "Methyl": 14}),
            "named and legacy modification declarations cannot be mixed",
        ),
        (
            serde_json::json!({"Acetyl": {"mass": 42, "sites": []}}),
            "modification `Acetyl` has no sites",
        ),
        (
            serde_json::json!({"Acetyl": {"mass": 42, "sites": ["K"], "name": "Other"}}),
            "modification `Acetyl` must not override its identity with a different name",
        ),
        (
            serde_json::json!({"Acetyl": {"mass": 42, "sites": ["K"], "max_count": 0}}),
            "max_count must be positive",
        ),
        (
            serde_json::json!({"Acetyl": {"mass": 42, "sites": ["Z"]}}),
            "invalid site `Z` for modification `Acetyl`",
        ),
        (
            serde_json::json!({"Acetyl": {"mass": 42, "sites": ["^K"]}}),
            "use explicit site `first_residue:K` instead of `^K`",
        ),
    ];
    for (value, expected) in cases {
        for section in ["static_mods", "variable_mods"] {
            if section == "static_mods" && expected == "max_count must be positive" {
                continue;
            }
            let message = error(section, value.clone());
            assert!(message.contains(expected), "{section} {value}: {message}");
        }
    }
    let message = error(
        "static_mods",
        serde_json::json!({
            "Carbamidomethyl": {"mass": 57.021464, "sites": ["C"]},
            "Other": {"mass": 58.0, "sites": ["C"]}
        }),
    );
    assert!(
        message.contains("different static modifications target the same site rule"),
        "{message}"
    );
}

#[test]
fn named_definitions_share_sites_and_ignore_repeated_sites() {
    let params = serde_json::from_value::<crate::database::Builder>(serde_json::json!({
        "variable_mods": {
            "Acetyl": {"mass": 42.010565, "sites": ["K", "K", "peptide_n_term"]},
            "Methyl": {"mass": 14.01565, "sites": ["K"]}
        },
        "static_mods": {
            "Carbamidomethyl": {"mass": 57.021464, "sites": ["C", "C"]}
        }
    }))
    .unwrap()
    .make_parameters();
    assert_eq!(params.static_mods.len(), 1);
    let lysine = &params.variable_mods[&ModificationSpecificity::Residue(b'K')];
    let mut names = lysine
        .iter()
        .map(|entry| entry.definition().name.unwrap().to_string())
        .collect::<Vec<_>>();
    names.sort();
    assert_eq!(names, vec!["Acetyl", "Methyl"]);
    assert_eq!(
        params.variable_mods[&ModificationSpecificity::PeptideN(None)].len(),
        1
    );
    assert_eq!(params.variable_mods.values().flatten().count(), 3);
}

#[test]
fn symbol_keyed_maps_accept_upstream_sage_syntax() {
    let builder: crate::database::Builder = serde_json::from_value(serde_json::json!({
        "static_mods": {"C": 57.021464, "^": 229.1629, "$K": 1.0, "[": 42.0, "]R": 0.5},
        "variable_mods": {"M": [15.9949], "^Q": [-17.026549], "S": [79.9663, 80.0]}
    }))
    .unwrap();
    let params = builder.make_parameters();
    assert_eq!(params.static_mods.len(), 5);
    assert_eq!(params.variable_mods.values().flatten().count(), 4);
}

#[test]
fn symbol_keyed_maps_reject_sage_plus_extensions() {
    let error = |value: serde_json::Value| {
        serde_json::from_value::<crate::database::Builder>(value)
            .err()
            .expect("configuration must be rejected")
            .to_string()
    };
    for (section, key, site) in [
        ("variable_mods", "~K", "internal_residue:K"),
        ("static_mods", "~K", "internal_residue:K"),
        ("variable_mods", "first_residue:K", "first_residue:K"),
        ("static_mods", "peptide_n_term", "peptide_n_term"),
    ] {
        let value = if section == "static_mods" {
            serde_json::json!(42.0)
        } else {
            serde_json::json!([42.0])
        };
        let message = error(serde_json::json!({ section: { key: value } }));
        assert!(message.contains(&format!("`{key}`")), "{message}");
        assert!(
            message.contains(&format!("\"sites\": [\"{site}\"]")),
            "{message}"
        );
        assert!(message.contains("DOCS.md"), "{message}");
    }
    for value in [
        serde_json::json!({"static_mods": {"C": {"mass": 57.021464, "name": "Carbamidomethyl"}}}),
        serde_json::json!({"variable_mods": {"M": [{"mass": 15.9949, "max_count": 1}]}}),
        serde_json::json!({"variable_mods": {"^K": [42.0106, {"mass": 42.0106}]}}),
    ] {
        let message = error(value);
        assert!(
            message.contains("must map to a mass, not an object"),
            "{message}"
        );
        assert!(message.contains("DOCS.md"), "{message}");
    }
}
