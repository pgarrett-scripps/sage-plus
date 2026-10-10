//! Scoping checks for the rustyms successor crates (mzcore + mzannotate).
//! Run with:
//!   cargo test -p sage-core --features rustyms-fragments --release \
//!     --test rustyms_scope -- --nocapture --ignored
//! `RUSTYMS_PROFORMA_TSV` points at a TSV of `localized_peptide<TAB>calcmass`.
#![cfg(feature = "rustyms-fragments")]

use std::time::Instant;

use mzannotate::{
    annotation::model::{FragmentationModel, PrimaryIonSeries},
    prelude::*,
};
use mzcore::{
    chemistry::{MolecularFormula, Molecule, OutputMolecularFormula},
    ontology::{Ontology, STATIC_ONTOLOGIES},
    prelude::*,
    system::{dalton, e, isize::Charge},
};
use sage_core::{
    enzyme::Digest,
    ion_series::{IonSeries, Kind},
    mass::{monoisotopic, H2O, NEUTRON, NH3, PROTON},
    peptide::Peptide,
};

fn formula_mass(formula: &str) -> f64 {
    MolecularFormula::pro_forma::<true, false>(formula)
        .unwrap()
        .monoisotopic_mass()
        .get::<dalton>()
}

fn unimod(name: &str) -> f64 {
    STATIC_ONTOLOGIES
        .get_by_name(&[Ontology::Unimod], name)
        .unwrap_or_else(|| panic!("unimod {name}"))
        .monoisotopic_mass()
        .mass()
        .get::<dalton>()
}

fn report(label: &str, sage: f64, rustyms: f64, worst: &mut f64, fails: &mut Vec<String>) {
    let diff = sage - rustyms;
    let f32_exact = (rustyms as f32) == (sage as f32);
    println!("{label:<34} sage {sage:>14.7} rustyms {rustyms:>14.7} diff {diff:>+11.2e} f32-identical {f32_exact}");
    *worst = worst.max(diff.abs());
    if diff.abs() > 1e-6 {
        fails.push(format!("{label} {diff:+.2e}"));
    }
}

#[test]
#[ignore]
fn masses_agree() {
    let mut worst = 0.0f64;
    let mut fails = Vec::new();
    for &aa in b"ACDEFGHIKLMNPQRSTVWYUOJ" {
        let ours = monoisotopic(aa) as f64;
        let theirs = AminoAcid::try_from(aa)
            .unwrap()
            .monoisotopic_masses()
            .iter()
            .next()
            .unwrap()
            .mass()
            .get::<dalton>();
        report(&format!("residue {}", aa as char), ours, theirs, &mut worst, &mut fails);
    }
    // Every named modification mass Sage Plus ships in presets, docs and code.
    let mods: &[(&str, f64, &str)] = &[
        ("Carbamidomethyl", 57.021464, "presets/*.json"),
        ("Oxidation", 15.994915, "presets/*.json"),
        ("Phospho", 79.966331, "presets/phospho.json"),
        ("TMTpro", 304.207146, "presets/trypsin-hcd-tmt.json"),
        ("TMTpro", 304.2071_f32 as f64, "tmt.rs:38 (f32)"),
        ("TMT6plex", 229.162932_f32 as f64, "tmt.rs:35 (f32)"),
        ("Acetyl", 42.010565, "modification_preview.rs:318"),
        ("Acetyl", 42.0106, "modification.rs:1233 example"),
        ("HexNAc", 203.079373, "modification.rs:1269 example"),
    ];
    for (name, ours, origin) in mods {
        report(&format!("{name} [{origin}]"), *ours, unimod(name), &mut worst, &mut fails);
    }
    // Physical constants (f32 in mass.rs, immonium.rs).
    report("H2O mass.rs:5", H2O as f64, formula_mass("H2O"), &mut worst, &mut fails);
    report("NH3 mass.rs:8", NH3 as f64, formula_mass("NH3"), &mut worst, &mut fails);
    report("PROTON mass.rs:6", PROTON as f64, formula_mass("H1:z+1"), &mut worst, &mut fails);
    report("NEUTRON (13C-12C) mass.rs:7", NEUTRON as f64, formula_mass("[13C1]C-1"), &mut worst, &mut fails);
    report("CO immonium.rs:29", 27.994915_f32 as f64, formula_mass("CO"), &mut worst, &mut fails);
    println!("worst |diff| {worst:.2e} Da; over 1e-6 Da: {}", fails.len());
    for f in &fails {
        println!("  over 1e-6: {f}");
    }
}

#[test]
#[ignore]
fn proforma_round_trip() {
    let path = std::env::var("RUSTYMS_PROFORMA_TSV").expect("RUSTYMS_PROFORMA_TSV");
    let text = std::fs::read_to_string(path).unwrap();
    let (mut ok, mut failed, mut mass_off) = (0usize, 0usize, 0usize);
    let mut worst = 0.0f64;
    let mut shown = 0;
    let start = Instant::now();
    for line in text.lines() {
        let mut cols = line.split('\t');
        let proforma = cols.next().unwrap();
        let calcmass: f64 = cols.next().unwrap().parse().unwrap();
        match PeptidoformIonSet::pro_forma(proforma, &STATIC_ONTOLOGIES) {
            Ok((set, _)) => {
                ok += 1;
                let masses = set.monoisotopic_masses();
                let mass = masses.iter().next().unwrap().mass().get::<dalton>();
                let diff = (mass - calcmass).abs();
                worst = worst.max(diff);
                if diff > 1e-3 {
                    mass_off += 1;
                    if shown < 5 {
                        shown += 1;
                        println!("mass off: {proforma} sage {calcmass} rustyms {mass}");
                    }
                }
            }
            Err(errors) => {
                failed += 1;
                if shown < 10 {
                    shown += 1;
                    println!("parse failure: {proforma}: {}", errors[0]);
                }
            }
        }
    }
    println!(
        "parsed {ok} failed {failed} mass>1e-3 {mass_off} worst |mass-calcmass| {worst:.2e} Da in {:?}",
        start.elapsed()
    );
}

/// Time b/y fragment generation, Sage `IonSeries` vs mzannotate.
#[test]
#[ignore]
fn fragment_timing() {
    let path = std::env::var("RUSTYMS_PROFORMA_TSV").expect("RUSTYMS_PROFORMA_TSV");
    let text = std::fs::read_to_string(path).unwrap();
    let sequences: Vec<String> = text
        .lines()
        .map(|l| l.split('\t').nth(2).unwrap().to_string())
        .map(|p| p.chars().filter(|c| c.is_ascii_uppercase()).collect())
        .collect();

    let peptides: Vec<Peptide> = sequences
        .iter()
        .map(|s| {
            Peptide::try_from(Digest {
                sequence: s.clone().into(),
                ..Default::default()
            })
            .unwrap()
        })
        .collect();
    let rounds = 20;
    let start = Instant::now();
    let mut sage_count = 0usize;
    let mut checksum = 0.0f64;
    for _ in 0..rounds {
        for peptide in &peptides {
            for kind in [Kind::B, Kind::Y] {
                for ion in IonSeries::new(peptide, kind) {
                    sage_count += 1;
                    checksum += (ion.monoisotopic_mass + PROTON) as f64;
                }
            }
        }
    }
    let sage_time = start.elapsed();

    let model = FragmentationModel::none()
        .clone()
        .b(PrimaryIonSeries::default())
        .y(PrimaryIonSeries::default());
    let parsed: Vec<PeptidoformIonSet> = sequences
        .iter()
        .map(|s| PeptidoformIonSet::pro_forma(s, &STATIC_ONTOLOGIES).unwrap().0)
        .collect();
    // OutputMonoIsotopic yields no fragments in mzannotate 0.2.0; use formulas.
    let start = Instant::now();
    let mut rms_count = 0usize;
    let mut rms_checksum = 0.0f64;
    for _ in 0..rounds {
        for p in &parsed {
            let fragments = p.generate_theoretical_fragments::<OutputMolecularFormula>(
                Charge::new::<e>(1),
                &model,
            );
            for f in &fragments {
                if let Some(mz) = f.mz(mzcore::chemistry::MassMode::Monoisotopic) {
                    rms_count += 1;
                    rms_checksum += mz.value;
                }
            }
        }
    }
    let rms_time = start.elapsed();
    let per = |d: std::time::Duration, n: usize| d.as_nanos() as f64 / n as f64;
    println!(
        "{} peptides x {rounds}: sage {} ions in {:?} ({:.1} ns/ion, sum {checksum:.1}); \
         mzannotate {} fragments in {:?} ({:.1} ns/fragment, sum {rms_checksum:.1}); ratio {:.1}x per peptide",
        peptides.len(),
        sage_count,
        sage_time,
        per(sage_time, sage_count),
        rms_count,
        rms_time,
        per(rms_time, rms_count),
        rms_time.as_secs_f64() / sage_time.as_secs_f64()
    );
}
