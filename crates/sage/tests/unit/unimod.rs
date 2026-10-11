#![allow(clippy::excessive_precision)]

use super::*;

#[test]
fn known_accessions() {
    let oxidation = delta_mass(35).expect("UNIMOD:35 missing");
    assert!((oxidation - 15.994915).abs() < 1e-4, "got {oxidation}");
    let carbamidomethyl = delta_mass(4).expect("UNIMOD:4 missing");
    assert!((carbamidomethyl - 57.021464).abs() < 1e-4);
    let phospho = delta_mass(21).expect("UNIMOD:21 missing");
    assert!((phospho - 79.966331).abs() < 1e-4);
}

#[test]
fn unknown_accession() {
    assert!(delta_mass(9_999_999).is_none());
}

#[test]
fn names_resolve_case_insensitively() {
    let m = mass_by_name("Oxidation").expect("Oxidation missing");
    assert!((m - 15.994915).abs() < 1e-4);
    let m2 = mass_by_name("oxidation").expect("lowercase Oxidation missing");
    assert_eq!(m, m2);
    assert!(mass_by_name("totally-not-a-mod").is_none());
    assert_eq!(canonical_name("phospho"), Some("Phospho"));
}

#[test]
fn labels_first_write_wins() {
    // Use a unique mass so this test doesn't collide with other suites.
    let m = 1234.56789_f32;
    register_label(m, "MyMod");
    register_label(m, "OtherMod");
    assert_eq!(label_for(m).as_deref(), Some("MyMod"));
    assert!(label_for(0.000_001).is_none());
}

#[test]
fn unimod_licence_ships_with_release_artifacts() {
    let root = std::path::Path::new(env!("CARGO_MANIFEST_DIR"));
    let licence = std::fs::read_to_string(root.join("data/LICENSE-unimod.txt")).unwrap();
    assert!(licence.contains("Copyright (C) 2002-2006 Unimod"));
    assert!(licence.contains("DESIGN SCIENCE LICENSE"));
    assert!(licence.contains("END OF TERMS AND CONDITIONS"));

    let workspace = root.join("../..");
    let read = |path: &str| std::fs::read_to_string(workspace.join(path)).unwrap();
    // Both binary archives and the Docker context copy the licence.
    assert_eq!(
        read(".github/workflows/release.yml")
            .matches("crates/sage/data/LICENSE-unimod.txt")
            .count(),
        3
    );
    assert!(read("Dockerfile").contains("/app/licenses/unimod.txt"));
    assert!(read("THIRD_PARTY_NOTICES.md").contains("## Unimod"));
}
