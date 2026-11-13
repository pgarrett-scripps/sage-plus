use super::*;

fn spectrum(peaks: &[(f32, u8, f32)]) -> ProcessedSpectrum {
    let mut peaks = peaks.to_vec();
    peaks.sort_by(|a, b| a.0.total_cmp(&b.0));
    ProcessedSpectrum {
        level: 2,
        masses: peaks.iter().map(|p| p.0).collect(),
        charges: peaks.iter().map(|p| p.1).collect(),
        charge_is_known: peaks.iter().map(|_| true).collect(),
        intensities: peaks.iter().map(|p| p.2).collect(),
        total_ion_current: peaks.iter().map(|p| p.2).sum(),
        ..Default::default()
    }
}

#[test]
fn decoy_twin_keeps_core_ions_and_moves_the_rest() {
    let composition = GlycanComposition::parse("HexNAc(4)Hex(5)Fuc(1)NeuAc(2)").unwrap();
    let ions = ion_set(&composition);
    assert!(ions.y.len() > 20);
    for ion in &ions.y {
        let shift = (ion.decoy_delta - ion.delta).abs();
        if ion.class == CORE as u8 {
            assert_eq!(shift, 0.0);
        } else {
            assert!((4.5..=12.5).contains(&shift), "{shift}");
        }
        // The full composition is the precursor, not a Y ion.
        assert!((ion.delta as f64 - composition.mass()).abs() > 1.0);
    }
    assert_eq!(ions.oxonium.len(), 4);
    // Deterministic.
    assert_eq!(ion_set(&composition).y, ions.y);
}

#[test]
fn true_y_ions_favour_the_target() {
    let composition = GlycanComposition::parse("HexNAc(4)Hex(5)Fuc(1)").unwrap();
    let ions = ion_set(&composition);
    let peptide = 1500.0f32;
    // Every target Y ion present, as neutral masses at charge 1.
    let peaks: Vec<_> = ions
        .y
        .iter()
        .map(|ion| (peptide + ion.delta, 1, 100.0))
        .collect();
    let query = spectrum(&peaks);
    let evidence = SpectrumEvidence::new(&query, Tolerance::Ppm(-10.0, 10.0), 3);
    let (target, decoy, fraction, anchored) = evidence.count(peptide, &ions);
    assert_eq!(target.y_matched(), ions.y.len() as u32);
    assert_eq!(target.y_generated(), decoy.y_generated());
    assert_eq!(decoy.matched[CORE], target.matched[CORE]);
    assert!(decoy.y_matched() < target.y_matched() / 2);
    assert!(anchored);
    assert!(fraction > 0.99);
    // Without NeuAc, NeuAc oxonium ions are generated as counter-evidence.
    assert_eq!(target.generated[NEUAC_ABSENT], 2);
    assert_eq!(target.matched[NEUAC_ABSENT], 0);
}

#[test]
fn ladder_counts_rungs_runs_and_control() {
    let peptide = 1500.0f32;
    let rung = |k: u8| {
        peptide + 2.0 * crate::config::HEXNAC as f32 + k as f32 * Monosaccharide::Hex.mass() as f32
    };
    // Rungs Hex 0, 1, 2, 4 and 5 present, 3 missing; one control peak.
    let mut peaks: Vec<_> = [0, 1, 2, 4, 5].map(|k| (rung(k), 1, 100.0)).to_vec();
    peaks.push((rung(1) + LADDER_CONTROL, 1, 50.0));
    let query = spectrum(&peaks);
    let evidence = SpectrumEvidence::new(&query, Tolerance::Ppm(-10.0, 10.0), 2);
    let ladder = evidence.ladder(peptide, 5);
    assert_eq!(
        ladder,
        Ladder {
            steps: 6,
            matched: 5,
            run: 3,
            control_matched: 1,
            control_run: 1,
        }
    );
    assert_eq!(evidence.ladder(peptide, 200).steps, LADDER_MAX_HEX + 1);
}

#[test]
fn random_rates_are_smoothed_probabilities() {
    let query = spectrum(&[(1000.0, 1, 1.0)]);
    let evidence = SpectrumEvidence::new(&query, Tolerance::Ppm(-10.0, 10.0), 2);
    let (y, ox) = evidence.random_rates(1500.0, 3000.0);
    assert!(y > 0.0 && y < 0.05);
    assert!(ox > 0.0 && ox < 0.05);
}
