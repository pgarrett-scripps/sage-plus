use crate::enzyme::{Digest, Enzyme, EnzymeParameters};
use crate::modification::NeutralLossMode;

use super::*;

#[test]
fn protein_accessions_are_inline_without_growing_peptides() {
    assert_eq!(
        std::mem::size_of::<ProteinAccessions>(),
        std::mem::size_of::<Vec<Arc<str>>>()
    );
    #[cfg(target_pointer_width = "64")]
    assert_eq!(std::mem::size_of::<Peptide>(), 144);

    let mut proteins = ProteinAccessions::new();
    proteins.push(Arc::from("P1"));
    assert!(!proteins.spilled());
    proteins.push(Arc::from("P2"));
    assert!(proteins.spilled());
}

#[test]
fn unmodified_peptides_use_compact_modification_storage() {
    let peptide = Peptide::try_from(Digest {
        sequence: "PEPTIDER".into(),
        ..Digest::default()
    })
    .unwrap();
    assert!(peptide.modifications.is_empty());
    assert_eq!(peptide.modification_at(3), 0.0);

    let modified = peptide
        .apply(
            &[(ModificationSpecificity::Residue(b'P'), 10.0, None)],
            &HashMap::default(),
            1,
            None,
        )
        .into_iter()
        .find(|peptide| {
            peptide.modification_count(ModificationSpecificity::Residue(b'P'), 10.0) > 0
        })
        .unwrap();
    assert_eq!(modified.modifications.len(), modified.sequence.len());
}

fn detailed_mod(
    mass: f32,
    name: &str,
    neutral_losses: &[f32],
    neutral_loss_mode: NeutralLossMode,
) -> Arc<ModificationDefinition> {
    Arc::new(ModificationDefinition {
        mass,
        name: Some(Arc::from(name)),
        neutral_losses: Arc::from(neutral_losses),
        neutral_loss_mode,
        channel_offsets: Arc::default(),
    })
}

fn var_mod_sequence(
    peptide: &Peptide,
    mods: &[(ModificationSpecificity, f32)],
    combo: usize,
) -> Vec<String> {
    let static_mods = HashMap::default();
    let mods_with_limits: Vec<(ModificationSpecificity, f32, Option<usize>)> =
        mods.iter().map(|&(s, m)| (s, m, None)).collect();
    peptide
        .clone()
        .apply(&mods_with_limits, &static_mods, combo, None)
        .into_iter()
        .map(|p| p.to_string())
        .collect::<Vec<_>>()
}

#[test]
fn full() {
    let sequence = "MPEPTIDEKMSAGEKEND";
    let tryp = EnzymeParameters {
        min_len: 0,
        max_len: 50,
        missed_cleavages: 0,
        enzyme: Enzyme::new("KR", "P", true, false),
    };

    let peptides = tryp
        .digest(sequence, Default::default())
        .into_iter()
        .map(Peptide::try_from)
        .collect::<Result<Vec<_>, _>>()
        .unwrap();

    assert_eq!(peptides.len(), 3);
    assert_eq!(peptides[0].to_string(), "MPEPTIDEK");
    assert_eq!(peptides[0].position, Position::Nterm);
    assert_eq!(peptides[1].to_string(), "MSAGEK");
    assert_eq!(peptides[1].position, Position::Internal);
    assert_eq!(peptides[2].to_string(), "END");
    assert_eq!(peptides[2].position, Position::Cterm);

    use ModificationSpecificity::*;

    let mods = [
        (ProteinN(None), 42.0),
        (ProteinC(None), 11.0),
        (PeptideN(None), 12.0),
        (PeptideC(None), 19.0),
    ];
    let a = var_mod_sequence(&peptides[0], &mods, 2);
    let b = var_mod_sequence(&peptides[1], &mods, 2);
    let c = var_mod_sequence(&peptides[2], &mods, 2);

    // Make sure no duplicates exist
    assert_eq!(
        a,
        vec![
            "MPEPTIDEK",
            "[+42]-MPEPTIDEK",
            "[+12]-MPEPTIDEK",
            "MPEPTIDEK-[+19]",
            "[+42]-MPEPTIDEK-[+19]",
            "[+12]-MPEPTIDEK-[+19]",
        ]
    );

    assert_eq!(
        b,
        vec![
            "MSAGEK",
            "[+12]-MSAGEK",
            "MSAGEK-[+19]",
            "[+12]-MSAGEK-[+19]",
        ]
    );

    assert_eq!(
        c,
        vec![
            "END",
            "END-[+11]",
            "[+12]-END",
            "END-[+19]",
            "[+12]-END-[+11]",
            "[+12]-END-[+19]",
        ]
    );
}

#[test]
fn test_variable_mods() {
    use ModificationSpecificity::*;
    let variable_mods = [(Residue(b'M'), 16.0f32), (Residue(b'C'), 57.)];
    let peptide = Peptide::try_from(Digest {
        sequence: "GCMGCMG".into(),
        ..Default::default()
    })
    .unwrap();

    let expected = vec![
        "GCMGCMG",
        "GCM[+16]GCMG",
        "GCMGCM[+16]G",
        "GC[+57]MGCMG",
        "GCMGC[+57]MG",
        "GCM[+16]GCM[+16]G",
        "GC[+57]M[+16]GCMG",
        "GCM[+16]GC[+57]MG",
        "GC[+57]MGCM[+16]G",
        "GCMGC[+57]M[+16]G",
        "GC[+57]MGC[+57]MG",
    ];

    let peptides = var_mod_sequence(&peptide, &variable_mods, 2);
    assert_eq!(peptides, expected);
}

#[test]
fn test_variable_mods_no_effeect() {
    use ModificationSpecificity::*;
    let variable_mods = [(Residue(b'M'), 16.0f32), (Residue(b'C'), 57.)];
    let peptide = Peptide::try_from(Digest {
        sequence: "AAAAAAAA".into(),
        ..Default::default()
    })
    .unwrap();

    let expected = vec!["AAAAAAAA"];
    let peptides = var_mod_sequence(&peptide, &variable_mods, usize::MAX);
    assert_eq!(peptides, expected);
}

#[test]
fn test_variable_mods_nterm() {
    use ModificationSpecificity::*;
    let variable_mods = [(PeptideN(None), 42.), (Residue(b'M'), 16.)];
    let peptide = Peptide::try_from(Digest {
        sequence: "GCMGCMG".into(),
        ..Default::default()
    })
    .unwrap();

    let expected = vec![
        "GCMGCMG",
        "[+42]-GCMGCMG",
        "GCM[+16]GCMG",
        "GCMGCM[+16]G",
        "[+42]-GCM[+16]GCMG",
        "[+42]-GCMGCM[+16]G",
        "GCM[+16]GCM[+16]G",
        "[+42]-GCM[+16]GCM[+16]G",
    ];

    let peptides = var_mod_sequence(&peptide, &variable_mods, 3);
    assert_eq!(peptides, expected);
}

#[test]
fn test_variable_mods_cterm() {
    use ModificationSpecificity::*;
    let variable_mods = [(PeptideC(None), 42.), (Residue(b'M'), 16.)];
    let peptide = Peptide::try_from(Digest {
        sequence: "GCMGCMG".into(),
        ..Default::default()
    })
    .unwrap();

    let expected = vec![
        "GCMGCMG",
        "GCMGCMG-[+42]",
        "GCM[+16]GCMG",
        "GCMGCM[+16]G",
        "GCM[+16]GCMG-[+42]",
        "GCMGCM[+16]G-[+42]",
        "GCM[+16]GCM[+16]G",
        "GCM[+16]GCM[+16]G-[+42]",
    ];

    let peptides = var_mod_sequence(&peptide, &variable_mods, 3);
    assert_eq!(peptides, expected);
}

#[test]
fn test_variable_mods_multi() {
    use ModificationSpecificity::*;
    let variable_mods = [(Residue(b'S'), 79.), (Residue(b'S'), 541.)];
    let peptide = Peptide::try_from(Digest {
        sequence: "GGGSGGGS".into(),
        ..Default::default()
    })
    .unwrap();

    let expected = vec![
        "GGGSGGGS",
        "GGGS[+79]GGGS",
        "GGGSGGGS[+79]",
        "GGGS[+541]GGGS",
        "GGGSGGGS[+541]",
        "GGGS[+79]GGGS[+79]",
        "GGGS[+79]GGGS[+541]",
        "GGGS[+541]GGGS[+79]",
        "GGGS[+541]GGGS[+541]",
    ];

    let peptides = var_mod_sequence(&peptide, &variable_mods, 2);
    assert_eq!(peptides, expected);
}

/// Check that picked-peptide approach will match forward and reverse peptides
#[test]
fn test_psuedo_forward() {
    let trypsin = crate::enzyme::EnzymeParameters {
        missed_cleavages: 0,
        min_len: 3,
        max_len: 30,
        enzyme: Enzyme::new("KR", "P", true, false),
    };

    let fwd = "MADEEKLPPGWEKRMSRSSGRVYYFNHITNASQWERPSGN";
    for digest in trypsin.digest(fwd, Default::default()) {
        let fwd = Peptide::try_from(digest.clone()).unwrap();
        let rev = Peptide::try_from(digest.reverse()).unwrap();

        assert_eq!(fwd.decoy, false);
        assert_eq!(rev.decoy, true);
        assert!(
            fwd.sequence.len() < 4 || fwd.sequence != rev.sequence,
            "{} {}",
            fwd,
            rev
        );
        assert_eq!(rev.reverse().to_string(), fwd.to_string());
    }
}

#[test]
fn apply_mods() {
    use ModificationSpecificity::*;
    let peptide = Peptide::try_from(Digest {
        sequence: "AACAACAA".into(),
        ..Default::default()
    })
    .unwrap();

    let expected = vec![
        "AAC[+57]AAC[+57]AA",
        "AAC[+30]AAC[+57]AA",
        "AAC[+57]AAC[+30]AA",
        "AAC[+30]AAC[+30]AA",
    ];

    let mut static_mods = HashMap::new();
    static_mods.insert(Residue(b'C'), 57.0);

    let variable_mods = [(Residue(b'C'), 30.0, None)];

    let peptides = peptide
        .apply(&variable_mods, &static_mods, 2, None)
        .into_iter()
        .map(|p| p.to_string())
        .collect::<Vec<_>>();

    assert_eq!(peptides, expected);
}

#[test]
fn test_per_mod_limit() {
    use ModificationSpecificity::*;
    // GCMGCMG has two M residues; limit oxidation to max 1 per peptide
    let variable_mods = [(Residue(b'M'), 16.0f32, Some(1))];
    let peptide = Peptide::try_from(Digest {
        sequence: "GCMGCMG".into(),
        ..Default::default()
    })
    .unwrap();

    let static_mods = HashMap::default();
    let peptides: Vec<String> = peptide
        .clone()
        .apply(&variable_mods, &static_mods, 2, None)
        .into_iter()
        .map(|p| p.to_string())
        .collect();

    // Should get unmodified + each single-M variant, but NOT the double-M variant
    let expected = vec!["GCMGCMG", "GCM[+16]GCMG", "GCMGCM[+16]G"];
    assert_eq!(peptides, expected);
}

#[test]
fn test_max_combinations() {
    use ModificationSpecificity::*;
    // GCMGCMG with oxidation and carbamidomethylation would normally yield many variants;
    // cap at 4 total (unmodified + 3 modified)
    let variable_mods = [(Residue(b'M'), 16.0f32, None), (Residue(b'C'), 57.0, None)];
    let peptide = Peptide::try_from(Digest {
        sequence: "GCMGCMG".into(),
        ..Default::default()
    })
    .unwrap();

    let static_mods = HashMap::default();
    let peptides: Vec<String> = peptide
        .clone()
        .apply(&variable_mods, &static_mods, 2, Some(4))
        .into_iter()
        .map(|p| p.to_string())
        .collect();

    // Cap at 4: unmodified + the first 3 single-mod variants (fewest PTMs first)
    assert_eq!(peptides.len(), 4);
    assert_eq!(peptides[0], "GCMGCMG");
}

#[test]
fn modification_sites() {
    use Site::*;
    let peptide = Peptide::try_from(Digest {
        sequence: "AACAACAA".into(),
        ..Default::default()
    })
    .unwrap();

    let mut mods = vec![];
    peptide.push_resi(&mut mods, ModificationSpecificity::Residue(b'C'), 16.0, 0);
    assert_eq!(mods, vec![(Sequence(2), 16.0, 0), (Sequence(5), 16.0, 0)]);
    mods.clear();

    peptide.push_resi(&mut mods, ModificationSpecificity::PeptideC(None), 16.0, 0);
    assert_eq!(mods, vec![(Cterm, 16.0, 0)]);
    mods.clear();

    peptide.push_resi(&mut mods, ModificationSpecificity::PeptideN(None), 16.0, 0);
    assert_eq!(mods, vec![(Nterm, 16.0, 0)]);
    mods.clear();

    let mut mods = vec![];
    for (idx, (residue, mass)) in [("^", 12.0), ("$", 200.0), ("C", 57.0), ("A", 43.0)]
        .iter()
        .enumerate()
    {
        peptide.push_resi(&mut mods, residue.parse().unwrap(), *mass, idx);
    }

    assert_eq!(
        mods,
        vec![
            (Nterm, 12.0, 0),
            (Cterm, 200.0, 1),
            (Sequence(2), 57.0, 2),
            (Sequence(5), 57.0, 2),
            (Sequence(0), 43.0, 3),
            (Sequence(1), 43.0, 3),
            (Sequence(3), 43.0, 3),
            (Sequence(4), 43.0, 3),
            (Sequence(6), 43.0, 3),
            (Sequence(7), 43.0, 3),
        ]
    );
}

#[test]
fn test_per_mod_limit_exactly_met() {
    use ModificationSpecificity::*;
    // Limit of 2 on a peptide with exactly 2 M residues — all combos should be allowed
    let variable_mods = [(Residue(b'M'), 16.0f32, Some(2))];
    let peptide = Peptide::try_from(Digest {
        sequence: "GCMGCMG".into(),
        ..Default::default()
    })
    .unwrap();

    let static_mods = HashMap::default();
    let peptides: Vec<String> = peptide
        .clone()
        .apply(&variable_mods, &static_mods, 2, None)
        .into_iter()
        .map(|p| p.to_string())
        .collect();

    // No restriction: unmodified + each single + double
    let expected = vec![
        "GCMGCMG",
        "GCM[+16]GCMG",
        "GCMGCM[+16]G",
        "GCM[+16]GCM[+16]G",
    ];
    assert_eq!(peptides, expected);
}

#[test]
fn test_per_mod_limit_zero() {
    use ModificationSpecificity::*;
    // Limit of 0 means this mod is entirely suppressed
    let variable_mods = [(Residue(b'M'), 16.0f32, Some(0))];
    let peptide = Peptide::try_from(Digest {
        sequence: "GCMGCMG".into(),
        ..Default::default()
    })
    .unwrap();

    let static_mods = HashMap::default();
    let peptides: Vec<String> = peptide
        .clone()
        .apply(&variable_mods, &static_mods, 2, None)
        .into_iter()
        .map(|p| p.to_string())
        .collect();

    assert_eq!(peptides, vec!["GCMGCMG"]);
}

#[test]
fn test_mixed_limited_and_unlimited() {
    use ModificationSpecificity::*;
    // M oxidation limited to 1; C carbamidomethylation unlimited
    // GCMGCMG has 2 M and 2 C
    let variable_mods = [
        (Residue(b'M'), 16.0f32, Some(1)),
        (Residue(b'C'), 57.0f32, None),
    ];
    let peptide = Peptide::try_from(Digest {
        sequence: "GCMGCMG".into(),
        ..Default::default()
    })
    .unwrap();

    let static_mods = HashMap::default();
    let peptides: Vec<String> = peptide
        .clone()
        .apply(&variable_mods, &static_mods, 2, None)
        .into_iter()
        .map(|p| p.to_string())
        .collect();

    // Should include all combos with ≤1 oxidized M,
    // but never both M residues oxidized simultaneously
    for p in &peptides {
        let oxid_count = p.matches("[+16]").count();
        assert!(oxid_count <= 1, "too many oxidations in: {}", p);
    }
    // Both C residues carbamidomethylated simultaneously should be present
    assert!(
        peptides.contains(&"GC[+57]MGC[+57]MG".to_string()),
        "expected double-C mod"
    );
    // Double oxidation should be absent
    assert!(
        !peptides.contains(&"GCM[+16]GCM[+16]G".to_string()),
        "double oxidation should be suppressed"
    );
}

#[test]
fn test_limits_are_per_mod_not_per_residue() {
    use ModificationSpecificity::*;
    // Both modifications target M, but only oxidation is limited to one.
    let variable_mods = [
        (Residue(b'M'), 16.0f32, Some(1)),
        (Residue(b'M'), 32.0f32, None),
    ];
    let peptide = Peptide::try_from(Digest {
        sequence: "GCMGCMG".into(),
        ..Default::default()
    })
    .unwrap();

    let peptides = peptide.apply(&variable_mods, &HashMap::default(), 2, None);
    let peptides = peptides.iter().map(ToString::to_string).collect::<Vec<_>>();

    assert!(!peptides.contains(&"GCM[+16]GCM[+16]G".to_string()));
    assert!(peptides.contains(&"GCM[+32]GCM[+32]G".to_string()));
}

#[test]
fn test_limits_support_more_than_64_mod_entries() {
    use ModificationSpecificity::*;
    let mut variable_mods = (1..=65)
        .map(|mass| (Residue(b'M'), mass as f32, None))
        .collect::<Vec<_>>();
    variable_mods[64].2 = Some(0);
    let peptide = Peptide::try_from(Digest {
        sequence: "GMG".into(),
        ..Default::default()
    })
    .unwrap();

    let peptides = peptide.apply(&variable_mods, &HashMap::default(), 1, None);

    assert_eq!(peptides.len(), 65); // unmodified + 64 allowed entries
    assert!(!peptides
        .iter()
        .any(|peptide| peptide.to_string().contains("[+65]")));
}

#[test]
fn test_max_combinations_only_unmodified() {
    use ModificationSpecificity::*;
    // cap of 1 means only the unmodified peptide is returned
    let variable_mods = [(Residue(b'M'), 16.0f32, None)];
    let peptide = Peptide::try_from(Digest {
        sequence: "GCMGCMG".into(),
        ..Default::default()
    })
    .unwrap();

    let static_mods = HashMap::default();
    let peptides: Vec<String> = peptide
        .clone()
        .apply(&variable_mods, &static_mods, 2, Some(1))
        .into_iter()
        .map(|p| p.to_string())
        .collect();

    assert_eq!(peptides, vec!["GCMGCMG"]);
}

#[test]
fn test_max_combinations_prefers_fewer_ptms() {
    use ModificationSpecificity::*;
    // GCMGCMG with oxidation (2 sites) — normally 3 variants (unmod + 2 single + 1 double)
    // cap at 3 means we get unmod + both singles but not the double
    let variable_mods = [(Residue(b'M'), 16.0f32, None)];
    let peptide = Peptide::try_from(Digest {
        sequence: "GCMGCMG".into(),
        ..Default::default()
    })
    .unwrap();

    let static_mods = HashMap::default();
    let peptides: Vec<String> = peptide
        .clone()
        .apply(&variable_mods, &static_mods, 2, Some(3))
        .into_iter()
        .map(|p| p.to_string())
        .collect();

    assert_eq!(peptides, vec!["GCMGCMG", "GCM[+16]GCMG", "GCMGCM[+16]G"]);
    // Double-mod must not appear — it would require cap > 3
    assert!(!peptides.contains(&"GCM[+16]GCM[+16]G".to_string()));
}

#[test]
fn names_follow_exact_modification_identity_and_decoy_position() {
    use ModificationSpecificity::*;
    let peptide = Peptide::try_from(Digest {
        sequence: "AMMAK".into(),
        ..Default::default()
    })
    .unwrap();
    let mods = [
        (
            Residue(b'M'),
            detailed_mod(15.9949, "Oxidation", &[], NeutralLossMode::Optional),
            Some(1),
        ),
        (
            Residue(b'M'),
            detailed_mod(15.9949, "AlternateName", &[], NeutralLossMode::Optional),
            Some(1),
        ),
    ];
    let peptides = peptide.apply(&mods, &HashMap::default(), 1, None);
    let rendered = peptides.iter().map(ToString::to_string).collect::<Vec<_>>();

    assert!(rendered.contains(&"AM[Oxidation]MAK".to_string()));
    assert!(rendered.contains(&"AM[AlternateName]MAK".to_string()));
    assert_ne!(
        peptides[1].applied_modifications,
        peptides[3].applied_modifications
    );

    let named = peptides
        .iter()
        .find(|peptide| peptide.to_string() == "AM[Oxidation]MAK")
        .unwrap();
    assert!(named.reverse().to_string().contains("[Oxidation]"));
}

#[test]
fn library_and_exhaustive_candidates_are_enumerated_together() {
    let peptide = Peptide::try_from(Digest {
        sequence: "MSS".into(),
        ..Default::default()
    })
    .unwrap();
    let phospho = Arc::new(ModificationDefinition {
        mass: 79.96633,
        name: Some(Arc::from("Phospho")),
        neutral_losses: Arc::from([]),
        neutral_loss_mode: NeutralLossMode::Optional,
        channel_offsets: Arc::default(),
    });
    let oxidation = Arc::new(ModificationDefinition {
        mass: 15.9949,
        name: Some(Arc::from("Oxidation")),
        neutral_losses: Arc::from([]),
        neutral_loss_mode: NeutralLossMode::Optional,
        channel_offsets: Arc::default(),
    });
    let rules = vec![
        VariableRule {
            specificity: ModificationSpecificity::Residue(b'S'),
            modification: phospho,
            max_count: Some(2),
            site_mode: SiteMode::Both,
            count_group: 0,
        },
        VariableRule {
            specificity: ModificationSpecificity::Residue(b'M'),
            modification: oxidation,
            max_count: Some(1),
            site_mode: SiteMode::Both,
            count_group: 1,
        },
    ];
    let library = vec![
        LibrarySite {
            position: 1,
            modification: Arc::from("Phospho"),
        },
        LibrarySite {
            position: 2,
            modification: Arc::from("Phospho"),
        },
    ];

    let variants = peptide.apply_rules(&rules, &library, &HashMap::new(), 1, 3, None);

    assert!(variants.iter().any(|peptide| {
        peptide.modification_at(0) != 0.0
            && peptide.modification_at(1) != 0.0
            && peptide.modification_at(2) != 0.0
    }));
    assert!(!variants.iter().any(|peptide| {
        peptide.modification_at(0) != 0.0
            && peptide.modification_at(1) == 0.0
            && peptide.modification_at(2) == 0.0
            && peptide.applied_modifications.len() > 1
    }));
}

#[test]
fn named_max_count_is_shared_across_residue_rules() {
    let peptide = Peptide::try_from(Digest {
        sequence: "ST".into(),
        ..Default::default()
    })
    .unwrap();
    let phospho = Arc::new(ModificationDefinition {
        mass: 79.96633,
        name: Some(Arc::from("Phospho")),
        neutral_losses: Arc::from([]),
        neutral_loss_mode: NeutralLossMode::Optional,
        channel_offsets: Arc::default(),
    });
    let rules = (*b"ST").map(|residue| VariableRule {
        specificity: ModificationSpecificity::Residue(residue),
        modification: phospho.clone(),
        max_count: Some(1),
        site_mode: SiteMode::Library,
        count_group: 0,
    });
    let library = vec![
        LibrarySite {
            position: 0,
            modification: Arc::from("Phospho"),
        },
        LibrarySite {
            position: 1,
            modification: Arc::from("Phospho"),
        },
    ];

    let variants = peptide.apply_rules(&rules, &library, &HashMap::new(), 0, 2, None);
    assert_eq!(variants.len(), 3);
    assert!(variants
        .iter()
        .all(|peptide| peptide.applied_modifications.len() <= 1));
}

#[test]
fn static_modification_names_are_rendered() {
    let peptide = Peptide::try_from(Digest {
        sequence: "ACK".into(),
        ..Default::default()
    })
    .unwrap();
    let static_mods = HashMap::from([(
        ModificationSpecificity::Residue(b'C'),
        detailed_mod(57.0215, "Carbamidomethyl", &[], NeutralLossMode::Optional),
    )]);

    let peptides = peptide.apply(&[], &static_mods, 0, None);
    assert_eq!(peptides[0].to_string(), "AC[Carbamidomethyl]K");
}
