use super::*;

fn spectrum(peaks: &[(f32, u8, bool, f32)]) -> ProcessedSpectrum {
    let mut peaks = peaks.to_vec();
    peaks.sort_by(|a, b| a.0.total_cmp(&b.0));
    ProcessedSpectrum {
        level: 2,
        masses: peaks.iter().map(|p| p.0).collect(),
        charges: peaks.iter().map(|p| p.1).collect(),
        charge_is_known: peaks.iter().map(|p| p.2).collect(),
        intensities: peaks.iter().map(|p| p.3).collect(),
        total_ion_current: peaks.iter().map(|p| p.3).sum(),
        ..Default::default()
    }
}

#[test]
fn parses_both_notations() {
    let long = GlycanComposition::parse("HexNAc(4)Hex(5)Fuc(1)NeuAc(2)").unwrap();
    let short = GlycanComposition::parse("N(4)H(5)F(1)A(2)").unwrap();
    assert_eq!(long, short);
    assert_eq!(long.to_string(), "HexNAc(4)Hex(5)Fuc(1)NeuAc(2)");
    // A2G2FS2, the common core-fucosylated disialylated biantennary glycan.
    assert!((long.mass() - 2350.8304).abs() < 1e-3, "{}", long.mass());
    assert!(GlycanComposition::parse("HexNAc(2)HexNAc(1)").is_err());
    assert!(GlycanComposition::parse("Kdn(1)").is_err());
}

#[test]
fn oxonium_ions_match_residue_masses() {
    let hexnac = Monosaccharide::HexNAc.mass() as f32 + PROTON;
    assert!((OXONIUM_IONS[0].mz - hexnac).abs() < 1e-3);
    let hexhexnac = (Monosaccharide::HexNAc.mass() + Monosaccharide::Hex.mass()) as f32 + PROTON;
    assert!((OXONIUM_IONS[5].mz - hexhexnac).abs() < 1e-3);
    let neuac = Monosaccharide::NeuAc.mass() as f32 + PROTON;
    assert!((OXONIUM_IONS[6].mz - neuac).abs() < 1e-3);
}

#[test]
fn explains_delta_with_isotope_error() {
    let library = GlycanLibrary::parse(
        "# comment\nHexNAc(2)Hex(5)\nHexNAc(4)Hex(5)Fuc(1)\nHexNAc(4)Hex(5)NeuAc(1)\n",
    )
    .unwrap();
    assert_eq!(library.len(), 3);
    let target = GlycanComposition::parse("HexNAc(4)Hex(5)Fuc(1)").unwrap();
    let delta = target.mass() + NEUTRON as f64;
    let hits = library.explain(delta, 0.02, 0..=1);
    assert_eq!(hits.len(), 1);
    assert_eq!(library.get(hits[0].0), Some(&target));
    assert_eq!(hits[0].1, 1);
}

#[test]
fn known_isomeric_and_near_isobaric_pairs() {
    // NeuAc + Hex and NeuGc + Fuc have the same elemental formula.
    let a = GlycanComposition::parse("HexNAc(4)Hex(5)NeuAc(1)").unwrap();
    let b = GlycanComposition::parse("HexNAc(4)Hex(4)Fuc(1)NeuGc(1)").unwrap();
    assert!((a.mass() - b.mass()).abs() < 1e-4);
    // Two Fuc sit 17 mDa from one NeuAc plus one neutron: only an
    // isotope-error assignment separates them.
    let c = GlycanComposition::parse("HexNAc(4)Hex(5)Fuc(2)").unwrap();
    let d = GlycanComposition::parse("HexNAc(4)Hex(5)NeuAc(1)").unwrap();
    let gap = c.mass() - (d.mass() + NEUTRON as f64);
    assert!(gap.abs() < 0.03, "{gap}");
}

#[test]
fn oxonium_gate_and_y_ladder() {
    let composition = GlycanComposition::parse("HexNAc(4)Hex(5)Fuc(1)").unwrap();
    let peptide = 1500.7;
    let y1 = peptide + Monosaccharide::HexNAc.mass() as f32;
    let query = spectrum(&[
        (OXONIUM_IONS[0].mz - PROTON, 1, false, 100.0),
        (OXONIUM_IONS[1].mz - PROTON, 1, false, 40.0),
        (OXONIUM_IONS[5].mz - PROTON, 1, false, 20.0),
        // Y0 observed as a deisotoped 2+ peak, Y1 as an undetermined 2+ m/z.
        (peptide, 2, true, 30.0),
        (y1 / 2.0 + PROTON - PROTON, 1, false, 50.0),
        (500.0, 1, false, 10.0),
    ]);
    let tol = Tolerance::Ppm(-20.0, 20.0);
    let oxonium = oxonium_evidence(&query, tol);
    assert_eq!(oxonium.count(), 3);
    assert!(oxonium.is_glyco(1));

    let ladder = n_glycan_core_y_ions(&composition);
    // Y0, 0,2X, and ten core fragments (five core ions, with and without Fuc).
    assert_eq!(ladder.len(), 12, "{ladder:?}");
    let evidence = y_ion_evidence(&query, peptide, &composition, 3, tol);
    assert_eq!(evidence.matched, 2);
    assert!(evidence.anchored);

    let blank = spectrum(&[(500.0, 1, false, 10.0)]);
    assert!(!oxonium_evidence(&blank, tol).is_glyco(1));
}

/// Sizes the composition ambiguity problem for the design doc:
/// `cargo test -p sage-glyco composition::test::ambiguity_report -- --ignored --nocapture`
#[test]
#[ignore]
fn ambiguity_report() {
    let space = n_glycan_composition_space();
    let human = space
        .iter()
        .copied()
        .filter(|c| c.count(Monosaccharide::NeuGc) == 0);
    for (label, library) in [
        ("human (no NeuGc)", GlycanLibrary::new(human)),
        ("mammalian (with NeuGc)", GlycanLibrary::new(space.clone())),
    ] {
        let (lo, hi) = library.mass_range().unwrap();
        println!(
            "{label}: {} compositions, mass range {lo:.1}..{hi:.1} Da",
            library.len()
        );
        for (tolerance, isotopes) in [(0.001, 0), (0.02, 0), (0.02, 1), (0.02, 2), (0.05, 2)] {
            let pairs = library.indistinguishable_pairs(tolerance, isotopes);
            let mut involved: Vec<usize> = pairs.iter().flat_map(|p| [p.0, p.1]).collect();
            involved.sort_unstable();
            involved.dedup();
            println!(
            "tolerance {tolerance} Da, isotope shifts 0..={isotopes}: {} pairs, {} of {} compositions ambiguous",
            pairs.len(),
            involved.len(),
            library.len()
        );
        }
        for (i, j, shift) in library.indistinguishable_pairs(0.02, 1).into_iter().take(6) {
            let (a, b) = (library.get(i).unwrap(), library.get(j).unwrap());
            println!("  {a} ~ {b} (+{shift} n, {:.4} Da)", b.mass() - a.mass());
        }
    }
}
