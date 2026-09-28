use super::*;
use crate::enzyme::Digest;
use crate::peptide::CompactModifications;

fn peptide(seq: &str) -> Peptide {
    Peptide::try_from(Digest {
        sequence: seq.into(),
        ..Default::default()
    })
    .unwrap()
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

fn settings() -> ImmoniumSettings {
    ImmoniumConfig::Enabled(true)
        .resolve(Tolerance::Ppm(-20.0, 20.0))
        .unwrap()
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
    // The built-in pY ion is the Y ion plus HPO3.
    let py = default_modified_ions()[0].mz;
    assert!((py - (residue_immonium_mz(b'Y') + 79.96633)).abs() < 5e-4);
}

#[test]
fn config_resolves_defaults_and_off() {
    let tol = Tolerance::Ppm(-10.0, 10.0);
    assert_eq!(ImmoniumConfig::Enabled(false).resolve(tol), None);
    let on = ImmoniumConfig::Enabled(true).resolve(tol).unwrap();
    assert!(!on.rescore);
    assert!(on.residues);
    assert_eq!(on.modified, default_modified_ions());
    assert_eq!(on.tolerance, tol);

    let custom: ImmoniumConfig = serde_json::from_str(
        r#"{"rescore": true, "residues": false,
            "modified": [{"name": "pH", "residue": "H", "modification": 79.96633, "mz": 190.0376}],
            "tolerance": {"da": [-0.01, 0.01]}}"#,
    )
    .unwrap();
    let custom = custom.resolve(tol).unwrap();
    assert!(custom.rescore);
    assert!(!custom.residues);
    assert_eq!(custom.modified[0].name, "pH");
    assert_eq!(custom.tolerance, Tolerance::Da(-0.01, 0.01));
    assert!(custom.validate().is_ok());

    assert!(serde_json::from_str::<ImmoniumConfig>(r#"{"formula": "H2O"}"#).is_err());
}

#[test]
fn validation_rejects_bad_ions() {
    let mut bad = settings();
    bad.modified[0].residue = 'b';
    assert!(bad.validate().unwrap_err().contains("pY"));
    let mut bad = settings();
    bad.modified[0].mz = -1.0;
    assert!(bad.validate().is_err());
    let mut bad = settings();
    bad.modified[0].modification = 0.0;
    assert!(bad.validate().is_err());
    let mut bad = settings();
    bad.modified[0].name = "a,b".into();
    assert!(bad.validate().is_err());
    let mut bad = settings();
    bad.modified = vec![default_modified_ions()[0].clone(); MAX_MODIFIED_IONS + 1];
    assert!(bad.validate().is_err());
    let mut bad = settings();
    bad.tolerance = Tolerance::Da(0.01, 0.02);
    assert!(bad.validate().is_err());
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
fn modified_residue_is_not_an_unmodified_residue() {
    let settings = settings();
    let py = default_modified_ions()[0].mz;
    let spectrum = spectrum(&[60.0, residue_immonium_mz(b'Y'), py]);

    // The only Y carries the phosphate: the Y ion is unexplained and pY is
    // explained.
    let mut phospho = peptide("AAYAAK");
    phospho.modifications = CompactModifications::from_sparse([(2, 79.96633)]);
    let evidence = settings.evaluate(&spectrum, &phospho);
    assert_eq!(evidence.explained, 0);
    assert_eq!(evidence.unexplained, 1);
    assert_eq!(evidence.modified_explained, 1);
    assert_eq!(evidence.modified_unexplained, 0);
    assert_eq!(settings.modified_names(evidence.modified_observed), "pY");

    // Phosphate on S instead: pY is unexplained, Y explained.
    let mut phospho_s = peptide("ASYAAK");
    phospho_s.modifications = CompactModifications::from_sparse([(1, 79.96633)]);
    let evidence = settings.evaluate(&spectrum, &phospho_s);
    assert_eq!(evidence.explained, 1);
    assert_eq!(evidence.modified_explained, 0);
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
fn residues_off_skips_residue_ions() {
    let mut settings = settings();
    settings.residues = false;
    let evidence = settings.evaluate(
        &spectrum(&[60.0, residue_immonium_mz(b'F')]),
        &peptide("AAFAAK"),
    );
    assert_eq!(evidence, ImmoniumEvidence::default());
}
