use super::*;
use crate::database::{Builder, Parameters};
use crate::enzyme::Digest;
use crate::modification::ModificationDefinition;
use crate::peptide::{AppliedModification, CompactModifications, ModificationKind};
use std::sync::Arc;

fn peptide(seq: &str) -> Peptide {
    Peptide::try_from(Digest {
        sequence: seq.into(),
        ..Default::default()
    })
    .unwrap()
}

/// `peptide` carrying the named modification at each zero-based position.
fn modified(seq: &str, name: &str, mass: f32, positions: &[u32]) -> Peptide {
    let mut peptide = peptide(seq);
    let mut definition = ModificationDefinition::bare(mass);
    definition.name = Some(name.into());
    let definition = Arc::new(definition);
    peptide.modifications =
        CompactModifications::from_applied(positions.iter().map(|&index| AppliedModification {
            site: Site::Sequence(index),
            modification: definition.clone(),
            kind: ModificationKind::Ordinary,
        }))
        .unwrap();
    peptide
}

/// A processed spectrum holding singly charged peaks at the given m/z.
fn spectrum(mzs: &[f32]) -> ProcessedSpectrum {
    let mut masses = mzs.iter().map(|mz| mz - PROTON).collect::<Vec<_>>();
    masses.sort_by(f32::total_cmp);
    ProcessedSpectrum {
        level: 2,
        intensities: vec![1.0; masses.len()],
        charges: vec![1; masses.len()],
        masses,
        ..Default::default()
    }
}

fn parameters(database: serde_json::Value) -> Parameters {
    serde_json::from_value::<Builder>(database)
        .unwrap()
        .make_parameters()
}

fn ions(database: serde_json::Value) -> Vec<ModifiedImmoniumIon> {
    let parameters = parameters(database);
    modified_ions(&parameters.static_mods, &parameters.variable_mods)
}

/// pY 216.0420 (Steen et al. 2001) and acK 126.0913 (Trelle & Jensen 2008),
/// declared on the modifications as in DOCS.md.
fn published_ions() -> Vec<ModifiedImmoniumIon> {
    ions(serde_json::json!({
        "variable_mods": {
            "Phospho": {"mass": 79.966331, "sites": ["S", "T", "Y"], "immonium_ions": {"Y": [216.0420]}},
            "Acetyl": {"mass": 42.010565, "sites": ["K"], "immonium_ions": [126.0913]}
        }
    }))
}

fn settings() -> ImmoniumSettings {
    ImmoniumConfig::Enabled(true)
        .resolve(Tolerance::Ppm(-20.0, 20.0), published_ions())
        .unwrap()
}

fn phospho_y() -> f32 {
    216.042
}

#[test]
fn residue_ion_masses_match_published_values() {
    // Hohmann et al. 2008 and standard immonium tables.
    for (residue, expected) in [
        (b'P', 70.0651),
        (b'V', 72.0808),
        (b'L', 86.0964),
        (b'H', 110.0713),
        (b'F', 120.0808),
        (b'Y', 136.0757),
        (b'W', 159.0917),
    ] {
        let mz = residue_immonium_mz(residue);
        assert!((mz - expected).abs() < 5e-4, "{} {mz}", residue as char);
    }
    // The published pY ion is the Y ion plus HPO3.
    assert!((phospho_y() - (residue_immonium_mz(b'Y') + 79.96633)).abs() < 5e-4);
}

#[test]
fn config_resolves_defaults_and_off() {
    let tol = Tolerance::Ppm(-10.0, 10.0);
    assert_eq!(ImmoniumConfig::Enabled(false).resolve(tol, vec![]), None);
    let on = ImmoniumConfig::Enabled(true).resolve(tol, vec![]).unwrap();
    // Enabled: rescoring and the residue ions are on by default, and there
    // are no built-in modified ions.
    assert!(on.rescore);
    assert!(on.residue_ions);
    assert!(on.modified.is_empty());
    assert_eq!(on.tolerance, tol);
    assert!(on.validate().is_ok());

    let custom: ImmoniumConfig = serde_json::from_str(
        r#"{"rescore": false, "residue_ions": false, "tolerance": {"da": [-0.01, 0.01]}}"#,
    )
    .unwrap();
    let custom = custom.resolve(tol, published_ions()).unwrap();
    assert!(!custom.rescore);
    assert!(!custom.residue_ions);
    assert_eq!(custom.modified.len(), 2);
    assert_eq!(custom.tolerance, Tolerance::Da(-0.01, 0.01));
    assert!(custom.validate().is_ok());

    // The old list of modified ions and formula strings are rejected.
    assert!(serde_json::from_str::<ImmoniumConfig>(r#"{"formula": "H2O"}"#).is_err());
    assert!(serde_json::from_str::<ImmoniumConfig>(r#"{"modified": []}"#).is_err());
    assert!(serde_json::from_str::<ImmoniumConfig>(r#"{"residues": true}"#).is_err());
}

#[test]
fn ions_come_only_from_modification_declarations() {
    let published = published_ions();
    assert_eq!(published.len(), 2);
    assert_eq!(published[0].label, "Acetyl@K");
    assert_eq!(published[0].modification, "Acetyl");
    assert!((published[0].mz - 126.0913).abs() < 1e-4);
    assert_eq!(published[0].residues, b"K");
    assert_eq!(published[1].label, "Phospho@Y");
    assert_eq!(published[1].residues, b"Y");

    // No declaration, no ion: phospho alone gives no pY ion.
    assert!(ions(serde_json::json!({
        "variable_mods": {"Phospho": {"mass": 79.966331, "sites": ["S", "T", "Y"]}}
    }))
    .is_empty());

    // A list applies to every site; the label lists the residues.
    let listed = ions(serde_json::json!({
        "variable_mods": {"Phospho": {"mass": 79.966331, "sites": ["S", "T", "Y"], "immonium_ions": [216.042]}}
    }));
    assert_eq!(listed.len(), 1);
    assert_eq!(listed[0].label, "Phospho@S/T/Y");

    // Static modifications declare ions the same way.
    let fixed = ions(serde_json::json!({
        "static_mods": {"Acetyl": {"mass": 42.010565, "sites": ["K"], "immonium_ions": [126.0913]}}
    }));
    assert_eq!(fixed.len(), 1);
    assert_eq!(fixed[0].label, "Acetyl@K");
}

#[test]
fn modification_immonium_ions_reject_bad_keys_and_values() {
    let parse = |database: serde_json::Value| serde_json::from_value::<Builder>(database);
    let error = parse(serde_json::json!({
        "variable_mods": {"Phospho": {"mass": 79.966331, "sites": ["S", "T"], "immonium_ions": {"Y": [216.042]}}}
    }))
    .err()
    .unwrap()
    .to_string();
    assert!(
        error.contains("immonium_ions") && error.contains("`Y`"),
        "{error}"
    );
    // Keys use the site vocabulary exactly.
    assert!(parse(serde_json::json!({
        "variable_mods": {"Acetyl": {"mass": 42.010565, "sites": ["internal_residue:K"], "immonium_ions": {"K": [126.0913]}}}
    }))
    .is_err());
    assert!(parse(serde_json::json!({
        "variable_mods": {"Acetyl": {"mass": 42.010565, "sites": ["internal_residue:K"], "immonium_ions": {"internal_residue:K": [126.0913]}}}
    }))
    .is_ok());
    for bad in [
        serde_json::json!([-1.0]),
        serde_json::json!(126.0913),
        serde_json::json!({"K": 126.0913}),
        serde_json::json!(["HPO3"]),
    ] {
        assert!(
            parse(serde_json::json!({
                "variable_mods": {"Acetyl": {"mass": 42.010565, "sites": ["K"], "immonium_ions": bad.clone()}}
            }))
            .is_err(),
            "{bad}"
        );
    }
    // Mass-offset PSMs are checked unmodified, so the ions could never be
    // explained.
    assert!(parse(serde_json::json!({
        "variable_mods": {"Phospho": {"mass": 79.966331, "sites": ["Y"], "max_count": 1,
            "search_mode": "mass_offset", "immonium_ions": [216.042]}}
    }))
    .is_err());
}

#[test]
fn validation_rejects_bad_settings() {
    let mut bad = settings();
    bad.modified[0].mz = -1.0;
    assert!(bad.validate().unwrap_err().contains("Acetyl@K"));
    let mut bad = settings();
    bad.modified[0].label = "a,b@K".into();
    assert!(bad.validate().is_err());
    let mut bad = settings();
    bad.modified = vec![bad.modified[0].clone(); MAX_MODIFIED_IONS + 1];
    assert!(bad.validate().is_err());
    let mut bad = settings();
    bad.tolerance = Tolerance::Da(0.01, 0.02);
    assert!(bad.validate().is_err());
    // Nothing to look for.
    let mut empty = settings();
    empty.residue_ions = false;
    empty.modified.clear();
    assert!(empty.validate().unwrap_err().contains("residue_ions"));
}

#[test]
fn counts_explained_missing_and_unexplained_residue_ions() {
    let settings = settings();
    // Peptide has F, Y and L; spectrum has F, L (explained), W (unexplained),
    // and no Y (missing). P is below the lowest peak, so neither counts.
    let peptide = peptide("PEPTFYLK");
    let spectrum = spectrum(&[
        residue_immonium_mz(b'L'),
        residue_immonium_mz(b'F'),
        residue_immonium_mz(b'W'),
        500.0,
    ]);
    let evidence = settings.evaluate(&spectrum, &peptide);
    assert_eq!(evidence.explained, 2);
    assert_eq!(evidence.missing, 1);
    assert_eq!(evidence.unexplained, 1);
    assert_eq!(evidence.residue_names(), "L/I,F,W");
    assert_eq!(
        evidence.modified_explained + evidence.modified_unexplained,
        0
    );
}

#[test]
fn phosphotyrosine_ion_is_explained_only_by_phospho_on_y() {
    let settings = settings();
    let spectrum = spectrum(&[60.0, residue_immonium_mz(b'Y'), phospho_y()]);

    // The only Y carries the phosphate: the Y ion is unexplained and the pY
    // ion is explained.
    let pyp = modified("AAYAAK", "Phospho", 79.96633, &[2]);
    let evidence = settings.evaluate(&spectrum, &pyp);
    assert_eq!(evidence.explained, 0);
    assert_eq!(evidence.unexplained, 1);
    assert_eq!(evidence.modified_explained, 1);
    assert_eq!(evidence.modified_unexplained, 0);
    assert_eq!(
        settings.modified_names(evidence.modified_observed),
        "Phospho@Y"
    );

    // Phosphate on S instead: the pY ion is unexplained, Y explained.
    let psp = modified("ASYAAK", "Phospho", 79.96633, &[1]);
    let evidence = settings.evaluate(&spectrum, &psp);
    assert_eq!(evidence.explained, 1);
    assert_eq!(evidence.modified_explained, 0);
    assert_eq!(evidence.modified_unexplained, 1);

    // Another modification of the same mass on Y is not Phospho.
    let other = modified("AAYAAK", "Sulfo", 79.95682, &[2]);
    let evidence = settings.evaluate(&spectrum, &other);
    assert_eq!(evidence.modified_explained, 0);
    assert_eq!(evidence.modified_unexplained, 1);

    // An unnamed modification (a ProForma mass delta) matches by mass.
    let mut unnamed = peptide("AAYAAK");
    unnamed.modifications = CompactModifications::from_sparse([(2, 79.96633)]);
    assert_eq!(settings.evaluate(&spectrum, &unnamed).modified_explained, 1);

    // Unmodified peptide: the pY ion is unexplained.
    let evidence = settings.evaluate(&spectrum, &peptide("AAYAAK"));
    assert_eq!(evidence.explained, 1);
    assert_eq!(evidence.modified_unexplained, 1);
}

#[test]
fn acetyl_lysine_ion_needs_acetyl_on_k() {
    let settings = settings();
    let spectrum = spectrum(&[60.0, 126.0913]);
    let ack = modified("AAKAAR", "Acetyl", 42.010565, &[2]);
    let evidence = settings.evaluate(&spectrum, &ack);
    assert_eq!(evidence.modified_explained, 1);
    assert_eq!(
        settings.modified_names(evidence.modified_observed),
        "Acetyl@K"
    );
    let evidence = settings.evaluate(&spectrum, &peptide("AAKAAR"));
    assert_eq!(evidence.modified_unexplained, 1);
}

#[test]
fn tolerance_and_charge_limit_matches() {
    let settings = settings();
    let peptide = peptide("AAFAAK");
    let f = residue_immonium_mz(b'F');
    // 30 ppm off: outside 20 ppm.
    let evidence = settings.evaluate(&spectrum(&[60.0, f * (1.0 + 30e-6)]), &peptide);
    assert_eq!(evidence.explained, 0);
    assert_eq!(evidence.missing, 1);
    // A peak at the right neutral mass but assigned charge 2 is not an
    // immonium ion.
    let mut charged = spectrum(&[60.0, f]);
    charged.charges[1] = 2;
    assert_eq!(settings.evaluate(&charged, &peptide).explained, 0);
    // Spectra without a charge column are read as singly charged.
    let mut no_charges = spectrum(&[60.0, f]);
    no_charges.charges.clear();
    assert_eq!(settings.evaluate(&no_charges, &peptide).explained, 1);
}

#[test]
fn residue_ions_off_skips_residue_ions() {
    let mut settings = settings();
    settings.residue_ions = false;
    let evidence = settings.evaluate(
        &spectrum(&[60.0, residue_immonium_mz(b'F')]),
        &peptide("AAFAAK"),
    );
    assert_eq!(evidence, ImmoniumEvidence::default());
}
