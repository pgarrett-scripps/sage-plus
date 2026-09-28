#![allow(clippy::excessive_precision)]

use super::*;
use crate::enzyme::Digest;
use crate::mass::PROTON;
use crate::modification::NeutralLossMode;
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
fn probe_scores_read_the_population_q_curve() {
    let evidence = [(100.0, false), (90.0, false), (80.0, true), (70.0, false)];
    let q = target_decoy_q_values(&evidence);
    let probes = q_values_at_scores(&evidence, &q, &[120.0, 90.0, 85.0, 75.0, 70.0, 10.0]);
    assert_eq!(probes, vec![0.5, 0.5, 2.0 / 3.0, 2.0 / 3.0, 2.0 / 3.0, 1.0]);
    assert_eq!(q_values_at_scores(&[], &[], &[1.0]), vec![1.0]);
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

#[test]
fn oxidation_rows_do_not_change_phospho_q_values() {
    let phospho = ModificationType::new("Phospho", PHOSPHO);
    let oxidation = ModificationType::new("Oxidation", 15.994915);
    // Phospho: 200 target wins scored 300..101, one decoy win at 150.
    let mut targets = (0..200)
        .map(|i| (phospho.clone(), 300.0 - i as f32, false))
        .collect::<Vec<_>>();
    targets.push((phospho.clone(), 150.0, true));
    let phospho_only = targets.clone();
    let probes = vec![(phospho.clone(), 250.0), (oxidation.clone(), 250.0)];
    let (alone, alone_probes) = target_decoy_q_values_by_type(&phospho_only, &probes[..1]);

    // Oxidation: coin-flip competition whose decoy wins outscore every
    // phospho row. Pooled, they would dominate the top of the ranking.
    targets.extend((0..100).map(|i| (oxidation.clone(), 1000.0 - i as f32, i % 2 == 0)));
    let (q, probe_q) = target_decoy_q_values_by_type(&targets, &probes);

    assert_eq!(&q[..phospho_only.len()], &alone[..]);
    assert_eq!(probe_q[0], alone_probes[0]);
    // Best phospho q: (0 + 1) / 150 targets above the decoy.
    assert!((q[0] - 1.0 / 150.0).abs() < 1e-6, "{}", q[0]);
    // Oxidation gets its own, much higher, q-values.
    assert!(q[phospho_only.len()..].iter().all(|&value| value > 0.5));
    assert!(probe_q[1] > 0.5);

    // The pooled competition would have let oxidation decoys set phospho q.
    let pooled = targets
        .iter()
        .map(|(_, score, decoy)| (*score, *decoy))
        .collect::<Vec<_>>();
    let pooled_q = target_decoy_q_values(&pooled);
    assert!(pooled_q[0] > 0.1, "{}", pooled_q[0]);
}

#[test]
fn per_type_q_values_match_single_type_competition() {
    let key = ModificationType::new("Phospho", PHOSPHO);
    let evidence = [(100.0, false), (90.0, false), (80.0, true), (70.0, false)];
    let keyed = evidence
        .iter()
        .map(|&(score, decoy)| (key.clone(), score, decoy))
        .collect::<Vec<_>>();
    let (q, probes) = target_decoy_q_values_by_type(&keyed, &[(key.clone(), 85.0)]);
    assert_eq!(q, target_decoy_q_values(&evidence));
    assert_eq!(probes, vec![2.0 / 3.0]);
    // A probe of a type without a target population cannot be accepted.
    let other = ModificationType::new("Oxidation", 15.994915);
    let (_, probes) = target_decoy_q_values_by_type(&keyed, &[(other, 500.0)]);
    assert_eq!(probes, vec![1.0]);
}

#[test]
fn modification_type_identity_uses_name_and_mass() {
    assert_eq!(
        ModificationType::new("Phospho", 79.96633),
        ModificationType::new("Phospho", 79.9665)
    );
    assert_ne!(
        ModificationType::new("Phospho", PHOSPHO),
        ModificationType::new("Sulfo", 79.95682)
    );
    assert_ne!(
        ModificationType::new("+15.9949", 15.9949),
        ModificationType::new("+15.9949", 15.9990)
    );
}

const H3PO4: f32 = 97.976896;

fn phospho_definition(losses: &[f32], mode: NeutralLossMode) -> Arc<ModificationDefinition> {
    Arc::new(ModificationDefinition {
        name: Some("Phospho".into()),
        neutral_losses: Arc::from(losses),
        neutral_loss_mode: mode,
        ..ModificationDefinition::bare(PHOSPHO)
    })
}

fn phospho_rules(
    definition: &Arc<ModificationDefinition>,
) -> [(ModificationSpecificity, Arc<ModificationDefinition>); 3] {
    [
        (ModificationSpecificity::Residue(b'S'), definition.clone()),
        (ModificationSpecificity::Residue(b'T'), definition.clone()),
        (ModificationSpecificity::Residue(b'Y'), definition.clone()),
    ]
}

fn spectrum_from(mut peaks: Vec<(f32, f32)>) -> ProcessedSpectrum {
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

/// Spectrum from the fragments main scoring would predict (including loss
/// forms), at charge 1, plus deterministic off-series noise peaks.
fn spectrum_with_losses(peptide: &Peptide, noise: usize) -> ProcessedSpectrum {
    let mut peaks: Vec<(f32, f32)> = Vec::new();
    for kind in [Kind::B, Kind::Y] {
        for group in crate::ion_series::IonGroupSeries::new(peptide, kind) {
            peaks.extend(
                group
                    .variants
                    .iter()
                    .map(|variant| (variant.monoisotopic_mass, 1000.0)),
            );
        }
    }
    peaks.extend((0..noise).map(|i| (150.0 + i as f32 * 37.1234, 100.0)));
    spectrum_from(peaks)
}

/// Phospho-carrying fragments of `AGSAATAAYAKR` with the phospho on T5.
fn carries_t5(kind: Kind, idx: usize) -> bool {
    match kind {
        Kind::B => idx >= 5,
        _ => idx < 5,
    }
}

/// With roll-up off, localization must match the per-mod-FLR base exactly.
/// The expected strings were recorded from the base commit (ad1b86f); only
/// the new `loss_rollup_peaks: None` field is added.
#[test]
fn localization_without_rollup_is_unchanged_from_base() {
    let definition = phospho_definition(&[], NeutralLossMode::Optional);
    let rules = phospho_rules(&definition);
    let single = peptide("AGSAATAAYAKR").with_mass_offset(Site::Sequence(5), &definition);
    let double = single.with_mass_offset(Site::Sequence(8), &definition);
    let mut observed = Vec::new();
    for scored in [&single, &double] {
        for noise in [0, 25] {
            let spectrum = spectrum_with_losses(scored, noise);
            let loc = localize(
                scored,
                &spectrum,
                &[Kind::B, Kind::Y],
                &rules,
                Tolerance::Ppm(-20.0, 20.0),
                None,
                3,
            );
            let off = localize_with_rollup(
                scored,
                &spectrum,
                &[Kind::B, Kind::Y],
                &rules,
                Tolerance::Ppm(-20.0, 20.0),
                None,
                3,
                LossRollup::Off,
            );
            assert_eq!(loc, off);
            observed.push(format!("{:?}", loc));
        }
    }
    let mass_only = [
        (ModificationSpecificity::Residue(b'S'), PHOSPHO),
        (ModificationSpecificity::Residue(b'T'), PHOSPHO),
    ];
    let mut legacy = peptide("AASAATAAGK");
    legacy.modifications = CompactModifications::from_sparse([(2, PHOSPHO)]);
    let loc = localize(
        &legacy,
        &spectrum_with_losses(&legacy, 10),
        &[Kind::B, Kind::Y],
        &mass_only,
        Tolerance::Ppm(-20.0, 20.0),
        None,
        3,
    );
    observed.push(format!("{:?}", loc));
    assert_eq!(observed.join("\n"), BASE_NO_ROLLUP_OUTPUT);
}

const BASE_NO_ROLLUP_OUTPUT: &str = r#"Localization { mods: [ModLocalization { mass: 79.96633, label: Some("Phospho"), site_count: 1, candidate_sites: 3, site_determining_ions: 36, site_determining_matched: 18, delta_score: 189.58427, target_decoy_score: 65.12185, decoy_winner: false, competition_eligible: true, localization_q_value: 1.0, best_sites: [SiteScore { attachment: Residue, position: 5, residue: 84, probability: 1.0 }], all_sites: [SiteScore { attachment: Residue, position: 2, residue: 83, probability: 1.1004576e-19 }, SiteScore { attachment: Residue, position: 5, residue: 84, probability: 1.0 }, SiteScore { attachment: Residue, position: 8, residue: 89, probability: 1.1004576e-19 }] }], loss_rollup_peaks: None }
Localization { mods: [ModLocalization { mass: 79.96633, label: Some("Phospho"), site_count: 1, candidate_sites: 3, site_determining_ions: 36, site_determining_matched: 18, delta_score: 169.79123, target_decoy_score: 58.524, decoy_winner: false, competition_eligible: true, localization_q_value: 1.0, best_sites: [SiteScore { attachment: Residue, position: 5, residue: 84, probability: 1.0 }], all_sites: [SiteScore { attachment: Residue, position: 2, residue: 83, probability: 1.0492452e-17 }, SiteScore { attachment: Residue, position: 5, residue: 84, probability: 1.0 }, SiteScore { attachment: Residue, position: 8, residue: 89, probability: 1.0492452e-17 }] }], loss_rollup_peaks: None }
Localization { mods: [ModLocalization { mass: 79.96633, label: Some("Phospho"), site_count: 2, candidate_sites: 3, site_determining_ions: 36, site_determining_matched: 18, delta_score: 189.80591, target_decoy_score: 189.80591, decoy_winner: false, competition_eligible: true, localization_q_value: 1.0, best_sites: [SiteScore { attachment: Residue, position: 5, residue: 84, probability: 1.0 }, SiteScore { attachment: Residue, position: 8, residue: 89, probability: 1.0 }], all_sites: [SiteScore { attachment: Residue, position: 2, residue: 83, probability: 1.0457052e-19 }, SiteScore { attachment: Residue, position: 5, residue: 84, probability: 1.0 }, SiteScore { attachment: Residue, position: 8, residue: 89, probability: 1.0 }] }], loss_rollup_peaks: None }
Localization { mods: [ModLocalization { mass: 79.96633, label: Some("Phospho"), site_count: 2, candidate_sites: 3, site_determining_ions: 36, site_determining_matched: 18, delta_score: 170.01297, target_decoy_score: 170.01297, decoy_winner: false, competition_eligible: true, localization_q_value: 1.0, best_sites: [SiteScore { attachment: Residue, position: 5, residue: 84, probability: 1.0 }, SiteScore { attachment: Residue, position: 8, residue: 89, probability: 1.0 }], all_sites: [SiteScore { attachment: Residue, position: 2, residue: 83, probability: 9.970164e-18 }, SiteScore { attachment: Residue, position: 5, residue: 84, probability: 1.0 }, SiteScore { attachment: Residue, position: 8, residue: 89, probability: 1.0 }] }], loss_rollup_peaks: None }
Localization { mods: [ModLocalization { mass: 79.96633, label: None, site_count: 1, candidate_sites: 2, site_determining_ions: 24, site_determining_matched: 12, delta_score: 177.45218, target_decoy_score: 62.0985, decoy_winner: false, competition_eligible: true, localization_q_value: 1.0, best_sites: [SiteScore { attachment: Residue, position: 2, residue: 83, probability: 1.0 }], all_sites: [SiteScore { attachment: Residue, position: 2, residue: 83, probability: 1.0 }, SiteScore { attachment: Residue, position: 5, residue: 84, probability: 1.7979707e-18 }] }], loss_rollup_peaks: None }"#;

#[test]
fn observed_rollup_folds_losses_into_the_most_intense_parent() {
    let parent = 800.0_f32;
    let spectrum = spectrum_from(vec![
        (300.0, 50.0),
        (parent - H3PO4, 40.0),
        (parent - 2.0 * H2O, 5.0),
        (parent - H2O, 20.0),
        (parent - NH3, 10.0),
        (parent, 100.0),
        // A weaker parent for the collision peak `parent - H3PO4`: the peak
        // one water above it, i.e. the unmodified form of the same fragment.
        (parent - H3PO4 + H2O, 30.0),
    ]);
    let (rolled, folded) = roll_up_observed_losses(
        &spectrum,
        &[NH3, H2O, H3PO4],
        Tolerance::Ppm(-20.0, 20.0),
        2,
    );
    // -H2O, -NH3 and -2 H2O (through -H2O) fold into the parent; the
    // collision peak goes to the more intense of its two observed parents.
    assert_eq!(folded, 4);
    assert_eq!(rolled.masses.len(), 3);
    let at = |mass: f32| {
        let i = rolled
            .masses
            .iter()
            .position(|m| (m - mass).abs() < 1e-3)
            .unwrap();
        rolled.intensities[i]
    };
    assert_eq!(at(parent), 100.0 + 20.0 + 10.0 + 5.0 + 40.0);
    assert_eq!(at(parent - H3PO4 + H2O), 30.0);
    assert_eq!(at(300.0), 50.0);
}

#[test]
fn observed_rollup_localizes_and_reports_folded_peaks() {
    let definition = phospho_definition(&[H3PO4], NeutralLossMode::Optional);
    let rules = phospho_rules(&definition);
    let truth = peptide("AGSAATAAYAKR").with_mass_offset(Site::Sequence(5), &definition);
    let base = spectrum_with_losses(&truth, 25);
    // Add a weaker water loss for every peak.
    let mut peaks = base
        .masses
        .iter()
        .copied()
        .zip(base.intensities.iter().copied())
        .collect::<Vec<_>>();
    peaks.extend(base.masses.iter().map(|m| (m - H2O, 50.0)));
    let spectrum = spectrum_from(peaks);
    let loc = localize_with_rollup(
        &truth,
        &spectrum,
        &[Kind::B, Kind::Y],
        &rules,
        Tolerance::Ppm(-20.0, 20.0),
        None,
        3,
        LossRollup::Observed,
    );
    assert_eq!(loc.mods[0].best_sites[0].position, 5);
    assert!(!loc.mods[0].decoy_winner);
    assert!(loc.loss_rollup_peaks.unwrap() >= base.masses.len() as u32);
}

/// The H3PO4 loss of a phospho fragment has the mass of the water loss of the
/// same unmodified fragment. Placement-aware roll-up must not let that peak
/// credit either placement: with only intact ions and H3PO4-loss peaks in
/// the spectrum, every arrangement matches exactly the slots intact scoring
/// gives it.
#[test]
fn placement_rollup_drops_the_phospho_water_collision() {
    let definition = phospho_definition(&[H3PO4], NeutralLossMode::Optional);
    let truth = peptide("AGSAATAAYAKR").with_mass_offset(Site::Sequence(5), &definition);
    let mut peaks = Vec::new();
    for kind in [Kind::B, Kind::Y] {
        for (idx, ion) in IonSeries::new(&truth, kind).enumerate() {
            peaks.push((ion.monoisotopic_mass, 1000.0));
            if carries_t5(kind, idx) {
                peaks.push((ion.monoisotopic_mass - H3PO4, 500.0));
            }
        }
    }
    let spectrum = spectrum_from(peaks);
    let group = modification_groups(&phospho_rules(&definition), &truth).remove(0);
    let candidates = group.candidates(&truth);
    let tol = Tolerance::Ppm(-20.0, 20.0);
    for chosen in [[2usize], [5], [8]] {
        let variant = build_variant(&truth, PHOSPHO, &candidates, &chosen);
        let slots = score_arrangement_with_losses(
            &variant,
            &spectrum,
            &[Kind::B, Kind::Y],
            &candidates,
            1,
            &group,
            tol,
            3,
        );
        let (intact, trials) = score_arrangement(
            &variant,
            &spectrum,
            &[Kind::B, Kind::Y],
            &candidates,
            1,
            tol,
            3,
        );
        let matched = slots.iter().filter(|s| s.matched).count() as u32;
        assert_eq!(slots.len() as u32, trials);
        assert_eq!(matched, intact, "arrangement {chosen:?}");
        assert!(slots.iter().all(|s| !s.loss_only));
        // Phospho fragments keep intact, -NH3 and -H2O (not -H3PO4);
        // unmodified fragments keep intact and -NH3 (not -H2O).
        assert!(slots.iter().all(|s| s.forms == 2 || s.forms == 3));
        assert!(slots.iter().any(|s| s.forms == 2));
        assert!(slots.iter().any(|s| s.forms == 3));
    }
}

#[test]
fn placement_rollup_credits_a_fragment_through_its_own_loss() {
    let definition = phospho_definition(&[H3PO4], NeutralLossMode::Optional);
    let rules = phospho_rules(&definition);
    let truth = peptide("AGSAATAAYAKR").with_mass_offset(Site::Sequence(5), &definition);
    // Only water losses of the phospho-carrying fragments are observed (no
    // intact ions), plus noise.
    let mut peaks = (0..25)
        .map(|i| (150.0 + i as f32 * 37.1234, 100.0))
        .collect::<Vec<_>>();
    for kind in [Kind::B, Kind::Y] {
        for (idx, ion) in IonSeries::new(&truth, kind).enumerate() {
            if carries_t5(kind, idx) {
                peaks.push((ion.monoisotopic_mass - H2O, 800.0));
            }
        }
    }
    let spectrum = spectrum_from(peaks);
    let run = |rollup| {
        localize_with_rollup(
            &truth,
            &spectrum,
            &[Kind::B, Kind::Y],
            &rules,
            Tolerance::Ppm(-20.0, 20.0),
            None,
            3,
            rollup,
        )
    };
    let placement = run(LossRollup::Placement);
    let m = &placement.mods[0];
    assert_eq!(m.best_sites[0].position, 5);
    assert!(m.delta_score > 0.0);
    assert!(placement.loss_rollup_peaks.unwrap() > 0);
    let off = run(LossRollup::Off);
    assert_eq!(off.loss_rollup_peaks, None);
    assert_eq!(off.mods[0].site_determining_matched, 0);
}

#[test]
fn poisson_binomial_matches_binomial_for_equal_probabilities() {
    for (successes, trials, p) in [
        (0u32, 10u32, 0.1),
        (3, 10, 0.1),
        (7, 20, 0.05),
        (20, 20, 0.3),
    ] {
        let expected = binomial_tail(successes, trials, p);
        let observed = poisson_binomial_tail(successes, &vec![p; trials as usize]);
        assert!(
            (expected - observed).abs() <= 1e-9 * expected,
            "{successes}/{trials}: {expected} vs {observed}"
        );
    }
    assert_eq!(poisson_binomial_tail(1, &[]), 1.0);
}
