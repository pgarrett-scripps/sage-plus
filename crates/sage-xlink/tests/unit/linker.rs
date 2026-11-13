use super::*;

#[test]
fn stub_masses_are_consistent_with_linker_mass() {
    let dsso = CleavableLinker::dsso();
    // alkene + sulfenic acid = crosslink.
    assert!((dsso.light_stub + dsso.extra_stubs[0] - dsso.crosslink_mass).abs() < 1e-3);
    let dsbu = CleavableLinker::dsbu();
    assert!((dsbu.light_stub + dsbu.heavy_stub - dsbu.crosslink_mass).abs() < 1e-3);
}

#[test]
fn settings_parse_presets_and_custom_linkers() {
    let settings: CrosslinkSettings = serde_json::from_str(r#"{"linker": "dsbu"}"#).unwrap();
    assert_eq!(settings.linker.resolve().unwrap().name, "DSBU");
    assert_eq!(settings.residues, "K");

    let settings: CrosslinkSettings = serde_json::from_str(
        r#"{"linker": {"name": "X", "crosslink_mass": 100.0, "light_stub": 40.0, "heavy_stub": 60.0}}"#,
    )
    .unwrap();
    assert_eq!(settings.linker.resolve().unwrap().doublet_spacing(), 20.0);

    let unknown: CrosslinkSettings = serde_json::from_str(r#"{"linker": "BS3"}"#).unwrap();
    assert!(unknown.validate().is_err());
    assert!(serde_json::from_str::<CrosslinkSettings>(r#"{"bogus": 1}"#).is_err());
}
