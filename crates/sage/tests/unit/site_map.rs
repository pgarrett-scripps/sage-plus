use super::*;

#[test]
fn list_applies_to_every_site() {
    let map: SiteMap<f32> = serde_json::from_str("[126.0913]").unwrap();
    assert_eq!(map, SiteMap::All(vec![126.0913]));
    assert_eq!(map.for_site("K"), &[126.0913]);
    assert_eq!(map.for_site("internal_residue:K"), &[126.0913]);
    assert!(!map.is_empty());
    assert!(map.validate_sites("f", "Acetyl", &["K"]).is_ok());
}

#[test]
fn map_applies_to_listed_sites_only() {
    let map: SiteMap<f32> = serde_json::from_str(r#"{"Y": [216.042]}"#).unwrap();
    assert_eq!(map.for_site("Y"), &[216.042]);
    assert!(map.for_site("S").is_empty());
    assert!(map.for_site("T").is_empty());
    assert_eq!(map.values().copied().collect::<Vec<_>>(), vec![216.042]);
    assert!(map.validate_sites("f", "Phospho", &["S", "T", "Y"]).is_ok());
}

#[test]
fn map_keys_must_be_declared_sites() {
    let map: SiteMap<f32> = serde_json::from_str(r#"{"W": [1.0]}"#).unwrap();
    let error = map
        .validate_sites("immonium_ions", "Phospho", &["S", "T", "Y"])
        .unwrap_err();
    assert!(error.contains("`W`"), "{error}");
    assert!(error.contains("Phospho"), "{error}");
    assert!(error.contains("`Y`"), "{error}");
    // Keys use the site vocabulary exactly: `K` is not `internal_residue:K`.
    let map: SiteMap<f32> = serde_json::from_str(r#"{"K": [1.0]}"#).unwrap();
    assert!(map
        .validate_sites("f", "Acetyl", &["internal_residue:K"])
        .is_err());
}

#[test]
fn empty_forms_and_bad_shapes() {
    assert!(SiteMap::<f32>::default().is_empty());
    let map: SiteMap<f32> = serde_json::from_str(r#"{"Y": []}"#).unwrap();
    assert!(map.is_empty());
    assert!(serde_json::from_str::<SiteMap<f32>>("216.042").is_err());
    assert!(serde_json::from_str::<SiteMap<f32>>(r#"{"Y": 216.042}"#).is_err());
    assert!(serde_json::from_str::<SiteMap<f32>>(r#"["x"]"#).is_err());
}

#[test]
fn round_trips_through_json() {
    for text in ["[126.0913]", r#"{"Y":[216.042]}"#] {
        let map: SiteMap<f32> = serde_json::from_str(text).unwrap();
        assert_eq!(serde_json::to_string(&map).unwrap(), text);
    }
}
