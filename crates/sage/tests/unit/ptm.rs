use super::*;
use crate::enzyme::Digest;
use crate::mass::PROTON;

fn peptide(seq: &str) -> Peptide {
    let mut peptide = Peptide::try_from(Digest {
        sequence: seq.into(),
        ..Default::default()
    })
    .unwrap();
    peptide.modifications.resize(peptide.sequence.len(), 0.0);
    peptide
}

/// Build a synthetic spectrum from the b/y ions of `peptide` (charge 1),
/// stored the way sage stores experimental peaks (neutral mass, sorted
/// ascending).
fn synthetic_spectrum(peptide: &Peptide) -> ProcessedSpectrum {
    let mut peaks: Vec<(f32, f32)> = Vec::new();
    for kind in [Kind::B, Kind::Y] {
        for ion in IonSeries::new(peptide, kind) {
            peaks.push((ion.monoisotopic_mass, 1000.0));
        }
    }
    peaks.sort_by(|a, b| a.0.total_cmp(&b.0));
    let tic = peaks.iter().map(|p| p.1).sum();
    let (masses, intensities) = peaks.into_iter().unzip();
    ProcessedSpectrum {
        level: 2,
        masses,
        intensities,
        total_ion_current: tic,
        ..Default::default()
    }
}

const PHOSPHO: f32 = 79.96633;

#[test]
fn num_combinations_basic() {
    assert_eq!(num_combinations(4, 2), 6);
    assert_eq!(num_combinations(3, 1), 3);
    assert_eq!(num_combinations(5, 0), 1);
    assert_eq!(num_combinations(2, 3), 0);
}

#[test]
fn target_decoy_q_values_are_monotonic() {
    let evidence = [(100.0, false), (90.0, false), (80.0, true), (70.0, false)];
    let q = target_decoy_q_values(&evidence);
    assert_eq!(q, vec![0.5, 0.5, 2.0 / 3.0, 2.0 / 3.0]);
}

#[test]
fn site_determining_rule() {
    // 2 candidates, 1 mod: a prefix containing exactly one candidate is
    // determining; containing zero or both is not.
    assert!(is_site_determining(1, 2, 1));
    assert!(!is_site_determining(0, 2, 1));
    assert!(!is_site_determining(2, 2, 1));
    // When every candidate is modified there is no ambiguity.
    assert!(!is_site_determining(1, 2, 2));
}

#[test]
fn localizes_single_phospho_to_correct_residue() {
    // Two candidate sites (S at idx 2, T at idx 5); the true site is the S.
    let mut truth = peptide("AASAATAA");
    truth.modifications[2] = PHOSPHO;
    let spectrum = synthetic_spectrum(&truth);

    // The peptide handed to the localizer carries the phospho on the S as
    // sage would have reported it.
    let scored = truth.clone();
    let potential = [
        (ModificationSpecificity::Residue(b'S'), PHOSPHO),
        (ModificationSpecificity::Residue(b'T'), PHOSPHO),
    ];

    let loc = localize(
        &scored,
        &spectrum,
        &[Kind::B, Kind::Y],
        &potential,
        Tolerance::Ppm(-10.0, 10.0),
        None,
        2,
    );

    assert_eq!(loc.mods.len(), 1);
    let m = &loc.mods[0];
    assert_eq!(m.site_count, 1);
    assert_eq!(m.candidate_sites, 2);
    assert_eq!(m.best_sites.len(), 1);
    // Best site is the true serine at index 2, with high probability.
    assert_eq!(m.best_sites[0].position, 2);
    assert!(
        m.best_sites[0].probability > 0.9,
        "probability was {}",
        m.best_sites[0].probability
    );
    // Probabilities across all candidate sites sum to ~1 (single mod).
    let total: f32 = m.all_sites.iter().map(|s| s.probability).sum();
    assert!((total - 1.0).abs() < 1e-3, "sum was {}", total);
    // The correct localization should be favored over the alternative.
    assert!(m.delta_score > 0.0, "delta_score was {}", m.delta_score);
    assert!(!m.decoy_winner);
    assert!(m.target_decoy_score > 0.0);
}

#[test]
fn unambiguous_when_single_candidate() {
    // Only one S in the peptide: the phospho is trivially localized.
    let mut truth = peptide("AAASAAA");
    truth.modifications[3] = PHOSPHO;
    let spectrum = synthetic_spectrum(&truth);
    let potential = [(ModificationSpecificity::Residue(b'S'), PHOSPHO)];

    let loc = localize(
        &truth,
        &spectrum,
        &[Kind::B, Kind::Y],
        &potential,
        Tolerance::Ppm(-10.0, 10.0),
        None,
        2,
    );
    assert_eq!(loc.mods.len(), 1);
    let m = &loc.mods[0];
    assert_eq!(m.candidate_sites, 1);
    assert_eq!(m.best_sites[0].position, 3);
    assert!((m.best_sites[0].probability - 1.0).abs() < 1e-6);
    assert_eq!(m.delta_score, 0.0);
}

#[test]
fn no_localization_without_target_mod() {
    // Peptide carries no phospho; nothing to localize.
    let pep = peptide("AASAATAA");
    let spectrum = synthetic_spectrum(&pep);
    let potential = [(ModificationSpecificity::Residue(b'S'), PHOSPHO)];
    let loc = localize(
        &pep,
        &spectrum,
        &[Kind::B, Kind::Y],
        &potential,
        Tolerance::Ppm(-10.0, 10.0),
        None,
        2,
    );
    assert!(loc.mods.is_empty());
    assert!(!has_localizable_modification(&pep, &potential));
    assert!(has_localizable_modification(
        &truth_with_phospho(),
        &potential
    ));
}

fn truth_with_phospho() -> Peptide {
    let mut peptide = peptide("AASAATAA");
    peptide.modifications[2] = PHOSPHO;
    peptide
}

#[test]
fn label_is_populated_when_registered() {
    // Use a mass unique to this test so the process-global label registry
    // isn't polluted for other suites (cf. unimod::tests).
    let unique_mass = 3131.31313_f32;
    crate::unimod::register_label(unique_mass, "TestPTM");
    let mut truth = peptide("AAASAAA");
    truth.modifications[3] = unique_mass;
    let spectrum = synthetic_spectrum(&truth);
    let potential = [(ModificationSpecificity::Residue(b'S'), unique_mass)];
    let loc = localize(
        &truth,
        &spectrum,
        &[Kind::B, Kind::Y],
        &potential,
        Tolerance::Ppm(-10.0, 10.0),
        None,
        2,
    );
    assert_eq!(loc.mods[0].label.as_deref(), Some("TestPTM"));
}

// Ensure PROTON import is exercised so charge math stays consistent with
// sage's peak convention in future edits.
#[test]
fn proton_constant_available() {
    assert!(PROTON > 1.0 && PROTON < 1.01);
}
