#![allow(clippy::excessive_precision)]

use super::*;
use crate::enzyme::Digest;
use crate::mass::PROTON;
use crate::peptide::CompactModifications;

fn peptide(seq: &str) -> Peptide {
    Peptide::try_from(Digest {
        sequence: seq.into(),
        ..Default::default()
    })
    .unwrap()
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
    truth.modifications = CompactModifications::from_sparse([(2, PHOSPHO)]);
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
    truth.modifications = CompactModifications::from_sparse([(3, PHOSPHO)]);
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
    peptide.modifications = CompactModifications::from_sparse([(2, PHOSPHO)]);
    peptide
}

#[test]
fn label_is_populated_when_registered() {
    // Use a mass unique to this test so the process-global label registry
    // isn't polluted for other suites (cf. unimod::tests).
    let unique_mass = 3131.31313_f32;
    crate::unimod::register_label(unique_mass, "TestPTM");
    let mut truth = peptide("AAASAAA");
    truth.modifications = CompactModifications::from_sparse([(3, unique_mass)]);
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
    const { assert!(PROTON > 1.0 && PROTON < 1.01) }
}

#[test]
fn tied_target_arrangements_cannot_inherit_confident_competition_q_values() {
    let mut precursor = peptide("AASAATAA");
    precursor.modifications = CompactModifications::from_sparse([(2, PHOSPHO)]);
    let potential = [
        (ModificationSpecificity::Residue(b'S'), PHOSPHO),
        (ModificationSpecificity::Residue(b'T'), PHOSPHO),
    ];
    let mut localization = localize(
        &precursor,
        &ProcessedSpectrum::default(),
        &[Kind::B, Kind::Y],
        &potential,
        Tolerance::Ppm(-10.0, 10.0),
        None,
        2,
    );
    let modification = &mut localization.mods[0];
    assert_eq!(modification.delta_score, 0.0);
    assert_eq!(modification.candidate_sites, 2);
    modification.set_competition_q_value(0.001);
    assert_eq!(modification.localization_q_value, 1.0);

    modification.delta_score = f32::NAN;
    modification.set_competition_q_value(0.001);
    assert_eq!(modification.localization_q_value, 1.0);
}

#[test]
fn single_candidate_can_retain_competition_confidence() {
    let mut precursor = peptide("AAASAAA");
    precursor.modifications = CompactModifications::from_sparse([(3, PHOSPHO)]);
    let mut localization = localize(
        &precursor,
        &synthetic_spectrum(&precursor),
        &[Kind::B, Kind::Y],
        &[(ModificationSpecificity::Residue(b'S'), PHOSPHO)],
        Tolerance::Ppm(-10.0, 10.0),
        None,
        2,
    );
    let modification = &mut localization.mods[0];
    assert_eq!(modification.delta_score, 0.0);
    modification.set_competition_q_value(0.001);
    assert_eq!(modification.localization_q_value, 0.001);
}

#[test]
fn positional_localization_preserves_first_internal_and_protein_rules() {
    let mut truth = peptide("KAKAKAAK");
    truth.modifications = CompactModifications::from_sparse([(2, 42.0)]);
    let internal = [(ModificationSpecificity::Internal(b'K'), 42.0)];
    let group = modification_groups(&internal, &truth).remove(0);
    assert_eq!(group.candidates(&truth), vec![2, 4]);
    let loc = localize(
        &truth,
        &synthetic_spectrum(&truth),
        &[Kind::B, Kind::Y],
        &internal,
        Tolerance::Ppm(-10.0, 10.0),
        None,
        2,
    );
    assert_eq!(
        loc.mods[0]
            .all_sites
            .iter()
            .map(|site| site.position)
            .collect::<Vec<_>>(),
        vec![2, 4]
    );
    let combined = [
        internal[0],
        (ModificationSpecificity::PeptideN(Some(b'K')), 42.0),
        (ModificationSpecificity::ProteinC(Some(b'K')), 42.0),
    ];
    let group = modification_groups(&combined, &truth).remove(0);
    assert_eq!(group.candidates(&truth), vec![0, 2, 4]);
    truth.position = crate::enzyme::Position::Cterm;
    let group = modification_groups(&combined, &truth).remove(0);
    assert_eq!(group.candidates(&truth), vec![0, 2, 4, 7]);
    let terminal = [(ModificationSpecificity::PeptideN(None), 42.0)];
    assert_eq!(
        modification_groups(&terminal, &truth)[0].candidates(&truth),
        vec![truth.sequence.len()]
    );
}

#[test]
fn equal_mass_named_modifications_keep_identity_and_occupied_sites() {
    let definition = |name: &str| {
        Arc::new(ModificationDefinition {
            name: Some(name.into()),
            ..ModificationDefinition::bare(42.0)
        })
    };
    let first = definition("FirstOnly");
    let internal = definition("InternalOnly");
    let truth = peptide("KAKAKAK")
        .with_mass_offset(Site::Sequence(0), &first)
        .with_mass_offset(Site::Sequence(2), &internal);
    let rules = [
        (ModificationSpecificity::PeptideN(Some(b'K')), first),
        (ModificationSpecificity::Internal(b'K'), internal),
    ];
    let groups = modification_groups(&rules, &truth);
    assert_eq!(groups.len(), 2);
    assert_eq!(groups[0].candidates(&truth), vec![0]);
    assert_eq!(groups[1].candidates(&truth), vec![2, 4]);
    let loc = localize(
        &truth,
        &synthetic_spectrum(&truth),
        &[Kind::B, Kind::Y],
        &rules,
        Tolerance::Ppm(-10.0, 10.0),
        None,
        2,
    );
    assert_eq!(loc.mods.len(), 2);
    assert_eq!(loc.mods[0].label.as_deref(), Some("FirstOnly"));
    assert_eq!(loc.mods[1].label.as_deref(), Some("InternalOnly"));
    assert!(loc.mods.iter().all(|entry| entry.site_count == 1));
    let other = definition("FixedOther");
    let truth = truth.with_mass_offset(Site::Sequence(4), &other);
    assert_eq!(groups[1].candidates(&truth), vec![2]);
}

#[test]
fn terminal_localization_is_typed_and_boundary_ambiguity_is_not_promoted() {
    use crate::ptm_library::Attachment;
    let definition = Arc::new(ModificationDefinition {
        name: Some("Acetyl".into()),
        ..ModificationDefinition::bare(42.010565)
    });
    for (site, spelling, attachment) in [
        (Site::Nterm, "peptide_n_term:K", Attachment::PeptideNTerm),
        (Site::Cterm, "peptide_c_term:K", Attachment::PeptideCTerm),
    ] {
        let truth = peptide("KAAAAAAK").with_mass_offset(site, &definition);
        let mut rules = vec![(
            spelling.parse::<ModificationSpecificity>().unwrap(),
            definition.clone(),
        )];
        let mut loc = localize(
            &truth,
            &synthetic_spectrum(&truth),
            &[Kind::B, Kind::Y],
            &rules,
            Tolerance::Ppm(-10.0, 10.0),
            None,
            2,
        );
        assert_eq!(loc.mods.len(), 1);
        assert_eq!(loc.mods[0].best_sites[0].attachment, attachment);
        assert!(!loc.mods[0].decoy_winner);
        loc.mods[0].set_competition_q_value(0.001);
        assert_eq!(loc.mods[0].localization_q_value, 0.001);
        rules.push((
            if site == Site::Nterm {
                ModificationSpecificity::PeptideN(Some(b'K'))
            } else {
                ModificationSpecificity::PeptideC(Some(b'K'))
            },
            definition.clone(),
        ));
        let mut ambiguous = localize(
            &truth,
            &synthetic_spectrum(&truth),
            &[Kind::B, Kind::Y],
            &rules,
            Tolerance::Ppm(-10.0, 10.0),
            None,
            2,
        );
        assert_eq!(ambiguous.mods[0].candidate_sites, 2);
        ambiguous.mods[0].set_competition_q_value(0.001);
        assert_eq!(ambiguous.mods[0].localization_q_value, 1.0);
        assert!(ambiguous.mods[0]
            .all_sites
            .iter()
            .any(|s| s.attachment == Attachment::Residue));
        assert!(ambiguous.mods[0]
            .all_sites
            .iter()
            .any(|s| s.attachment == attachment));
    }
}
