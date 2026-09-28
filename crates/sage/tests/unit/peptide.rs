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
    assert_eq!(modified.modifications.len(), 1);
    assert!(!modified.modifications.spilled());
}

#[test]
fn compact_modification_layout_and_inline_capacity() {
    assert_eq!(std::mem::size_of::<EncodedModification>(), 2);
    assert_eq!(std::mem::size_of::<CompactModifications>(), 32);
    assert_eq!(std::mem::size_of::<Peptide>(), 144);

    let inline = CompactModifications::from_sparse((0..4).map(|index| (index, index as f32 + 1.0)));
    assert_eq!(inline.len(), 4);
    assert!(!inline.spilled());

    let spilled =
        CompactModifications::from_sparse((0..5).map(|index| (index, index as f32 + 1.0)));
    assert_eq!(spilled.len(), 5);
    assert!(spilled.spilled());
}

#[test]
fn compact_forward_cursor_sums_stacked_residue_modifications() {
    let modifications = CompactModifications::from_applied([
        AppliedModification {
            site: Site::Nterm,
            modification: Arc::new(ModificationDefinition::bare(4.0)),
            kind: ModificationKind::Ordinary,
        },
        AppliedModification {
            site: Site::Sequence(2),
            modification: Arc::new(ModificationDefinition::bare(10.0)),
            kind: ModificationKind::Ordinary,
        },
        AppliedModification {
            site: Site::Sequence(2),
            modification: Arc::new(ModificationDefinition::bare(2.5)),
            kind: ModificationKind::Ordinary,
        },
        AppliedModification {
            site: Site::Sequence(4),
            modification: Arc::new(ModificationDefinition::bare(3.0)),
            kind: ModificationKind::Ordinary,
        },
    ])
    .unwrap();
    let mut cursor = 0;
    assert_eq!(modifications.mass_at_with_cursor(0, &mut cursor), 0.0);
    assert_eq!(modifications.mass_at_with_cursor(1, &mut cursor), 0.0);
    assert_eq!(modifications.mass_at_with_cursor(2, &mut cursor), 12.5);
    assert_eq!(modifications.mass_at_with_cursor(3, &mut cursor), 0.0);
    assert_eq!(modifications.mass_at_with_cursor(4, &mut cursor), 3.0);
}

#[test]
fn applying_rules_merges_an_existing_compact_lookup() {
    let peptide = Peptide {
        sequence: (&b"ACDE"[..]).into(),
        modifications: CompactModifications::from_sparse([(0, 3.0)]),
        ..Peptide::default()
    };
    let modified = peptide
        .apply(
            &[(ModificationSpecificity::Residue(b'D'), 5.0, None)],
            &HashMap::default(),
            1,
            None,
        )
        .into_iter()
        .find(|peptide| peptide.modification_at(2) == 5.0)
        .unwrap();
    assert_eq!(modified.modification_at(0), 3.0);
    assert_eq!(modified.modification_at(2), 5.0);
}

#[test]
fn compact_encoding_accepts_255_residues_and_rejects_256() {
    let peptide = Peptide::try_from(Digest {
        sequence: "A".repeat(255).into(),
        ..Digest::default()
    })
    .unwrap();
    let mut peptide = peptide;
    peptide.modifications = CompactModifications::from_sparse([(254, 12.5)]);
    assert_eq!(peptide.modification_at(254), 12.5);

    assert_eq!(
        Peptide::try_from(Digest {
            sequence: "A".repeat(256).into(),
            ..Digest::default()
        }),
        Err(PeptideError::SequenceTooLong {
            length: 256,
            maximum: 255,
        })
    );
}

#[test]
fn compact_lookup_rejects_more_than_255_definition_variants() {
    let definitions = |count| {
        (0..count).map(|index| {
            (
                Site::Sequence(0),
                Arc::new(ModificationDefinition::bare(index as f32 + 0.5)),
                ModificationKind::Ordinary,
            )
        })
    };
    assert!(ModificationLookup::from_definitions(definitions(255)).is_ok());
    assert_eq!(
        ModificationLookup::from_definitions(definitions(256)).unwrap_err(),
        ModificationLookupError { definitions: 256 }
    );
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
        site_losses: None,
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
        clip_n_term_met: false,
        min_len: 0,
        max_len: 50,
        missed_cleavages: 0,
        enzyme: Enzyme::new("KR", "P", true, false),
        ambiguous_variants: None,
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
        clip_n_term_met: false,
        missed_cleavages: 0,
        min_len: 3,
        max_len: 30,
        enzyme: Enzyme::new("KR", "P", true, false),
        ambiguous_variants: None,
    };

    let fwd = "MADEEKLPPGWEKRMSRSSGRVYYFNHITNASQWERPSGN";
    for digest in trypsin.digest(fwd, Default::default()) {
        let fwd = Peptide::try_from(digest.clone()).unwrap();
        let rev = Peptide::try_from(digest.reverse()).unwrap();

        assert!(!fwd.decoy);
        assert!(rev.decoy);
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
        peptides[1].applied_modifications().collect::<Vec<_>>(),
        peptides[3].applied_modifications().collect::<Vec<_>>()
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
        site_losses: None,
        neutral_loss_mode: NeutralLossMode::Optional,
        channel_offsets: Arc::default(),
    });
    let oxidation = Arc::new(ModificationDefinition {
        mass: 15.9949,
        name: Some(Arc::from("Oxidation")),
        neutral_losses: Arc::from([]),
        site_losses: None,
        neutral_loss_mode: NeutralLossMode::Optional,
        channel_offsets: Arc::default(),
    });
    let rules = vec![
        VariableRule {
            specificity: ModificationSpecificity::Residue(b'S'),
            modification: phospho,
            max_count: Some(2),
            max_total_count: None,
            site_mode: SiteMode::Both,
            count_group: 0,
        },
        VariableRule {
            specificity: ModificationSpecificity::Residue(b'M'),
            modification: oxidation,
            max_count: Some(1),
            max_total_count: None,
            site_mode: SiteMode::Both,
            count_group: 1,
        },
    ];
    let library = vec![
        LibrarySite {
            attachment: Default::default(),
            position: 1,
            modification: Arc::from("Phospho"),
        },
        LibrarySite {
            attachment: Default::default(),
            position: 2,
            modification: Arc::from("Phospho"),
        },
    ];

    let static_mods = HashMap::new();
    let labels = LabelModificationCache::new(rules.iter().map(|rule| &rule.modification), &[]);
    let lookup = ModificationLookup::for_rules(&rules, &static_mods, &[], &labels).unwrap();
    let plan = ModificationPlan::new(&rules, &static_mods, lookup, 1, 3, None);
    let variants = peptide.apply_rules(&plan, &library);

    assert!(variants.iter().any(|peptide| {
        peptide.modification_at(0) != 0.0
            && peptide.modification_at(1) != 0.0
            && peptide.modification_at(2) != 0.0
    }));
    assert!(!variants.iter().any(|peptide| {
        peptide.modification_at(0) != 0.0
            && peptide.modification_at(1) == 0.0
            && peptide.modification_at(2) == 0.0
            && peptide.applied_modifications().len() > 1
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
        site_losses: None,
        neutral_loss_mode: NeutralLossMode::Optional,
        channel_offsets: Arc::default(),
    });
    let rules = (*b"ST").map(|residue| VariableRule {
        specificity: ModificationSpecificity::Residue(residue),
        modification: phospho.clone(),
        max_count: Some(1),
        max_total_count: None,
        site_mode: SiteMode::Library,
        count_group: 0,
    });
    let library = vec![
        LibrarySite {
            attachment: Default::default(),
            position: 0,
            modification: Arc::from("Phospho"),
        },
        LibrarySite {
            attachment: Default::default(),
            position: 1,
            modification: Arc::from("Phospho"),
        },
    ];

    let static_mods = HashMap::new();
    let labels = LabelModificationCache::new(rules.iter().map(|rule| &rule.modification), &[]);
    let lookup = ModificationLookup::for_rules(&rules, &static_mods, &[], &labels).unwrap();
    let plan = ModificationPlan::new(&rules, &static_mods, lookup, 0, 2, None);
    let variants = peptide.apply_rules(&plan, &library);
    assert_eq!(variants.len(), 3);
    assert!(variants
        .iter()
        .all(|peptide| peptide.applied_modifications().len() <= 1));
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

fn plain_peptide(sequence: &str) -> Peptide {
    Peptide::try_from(Digest {
        sequence: sequence.into(),
        ..Default::default()
    })
    .unwrap()
}

#[test]
fn static_terminal_mod_yields_to_variable_terminal_mod() {
    use ModificationSpecificity::*;
    let static_mods = HashMap::from([(PeptideN(None), 42.0f32)]);
    let variants =
        plain_peptide("PEPK").apply(&[(PeptideN(None), 12.0, None)], &static_mods, 1, None);
    let rendered = variants.iter().map(ToString::to_string).collect::<Vec<_>>();
    // The static group fills the free N-terminus, but never stacks on top of
    // a variable modification that already occupies it.
    assert_eq!(rendered, vec!["[+42]-PEPK", "[+12]-PEPK"]);
    let base = plain_peptide("PEPK").monoisotopic;
    assert!((variants[0].monoisotopic - (base + 42.0)).abs() < 1e-3);
    assert!((variants[1].monoisotopic - (base + 12.0)).abs() < 1e-3);
}

#[test]
fn static_terminal_and_residue_mods_both_apply_and_add_mass() {
    use ModificationSpecificity::*;
    let static_mods = HashMap::from([(Residue(b'C'), 57.0f32), (PeptideC(None), 1.0)]);
    let peptide = plain_peptide("ACCK");
    let base = peptide.monoisotopic;
    let variants = peptide.apply(&[], &static_mods, 0, None);
    assert_eq!(variants.len(), 1);
    assert_eq!(variants[0].to_string(), "AC[+57]C[+57]K-[+1]");
    assert_eq!(variants[0].cterm, Some(1.0));
    assert!((variants[0].monoisotopic - (base + 115.0)).abs() < 1e-3);
}

#[test]
fn modified_decoys_mirror_internal_positions_and_keep_terminal_mods() {
    use ModificationSpecificity::*;
    let variants = plain_peptide("ACDEFK").apply(
        &[
            (Residue(b'C'), 57.0f32, None),
            (Residue(b'K'), 8.0, None),
            (PeptideN(None), 42.0, None),
        ],
        &HashMap::default(),
        3,
        None,
    );
    let target = variants
        .iter()
        .find(|peptide| peptide.to_string() == "[+42]-AC[+57]DEFK[+8]")
        .unwrap();
    let decoy = target.reverse();
    assert!(decoy.decoy);
    assert_eq!(decoy.to_string(), "[+42]-AFEDC[+57]K[+8]");
    assert_eq!(decoy.label(), -1);
    assert_eq!(target.label(), 1);
    assert_eq!(decoy.monoisotopic, target.monoisotopic);
    assert_eq!(decoy.nterm, Some(42.0));
    assert_eq!(decoy.modification_at(4), 57.0);
    assert_eq!(decoy.modification_at(1), 0.0);
    // Reversing twice restores the target exactly.
    let restored = decoy.reverse();
    assert!(!restored.decoy);
    assert_eq!(restored.to_string(), target.to_string());
    assert_eq!(restored.modification_at(1), 57.0);
}

#[test]
fn two_residue_decoys_keep_sequence_and_modifications() {
    use ModificationSpecificity::*;
    let variants = plain_peptide("CK").apply(
        &[(Residue(b'C'), 57.0f32, None)],
        &HashMap::default(),
        1,
        None,
    );
    let decoy = variants[1].reverse();
    assert!(decoy.decoy);
    assert_eq!(decoy.to_string(), "C[+57]K");
}

#[test]
fn decoy_protein_names_follow_generation_mode() {
    let mut peptide = plain_peptide("PEPK");
    peptide.proteins = smallvec::smallvec![Arc::from("P1"), Arc::from("P2")];
    assert_eq!(peptide.proteins("rev_", true), "P1;P2");
    let decoy = peptide.reverse();
    assert_eq!(decoy.proteins("rev_", true), "rev_P1;rev_P2");
    // FASTA-supplied decoys already carry their tag.
    assert_eq!(decoy.proteins("rev_", false), "P1;P2");
}

#[test]
fn mass_offset_variants_match_indexed_variant_mass_exactly() {
    use ModificationSpecificity::*;
    let peptide = plain_peptide("PEPMK");
    let oxidation = Arc::new(ModificationDefinition::bare(15.9949));
    let indexed = peptide
        .clone()
        .apply(
            &[(Residue(b'M'), 15.9949f32, None)],
            &HashMap::default(),
            1,
            None,
        )
        .pop()
        .unwrap();
    let offset = peptide.with_mass_offset(Site::Sequence(3), &oxidation);
    assert_eq!(offset.to_string(), "PEPM[+15.9949]K");
    assert_eq!(offset.to_string(), indexed.to_string());
    assert_eq!(offset.monoisotopic, indexed.monoisotopic);
    // The source peptide is untouched.
    assert!(peptide.modifications.is_empty());

    // A decoy takes its target's mass bit for bit.
    let decoy = peptide.reverse();
    assert_eq!(decoy.to_string(), "PMPEK");
    let decoy_offset = decoy.with_mass_offset(Site::Sequence(1), &oxidation);
    assert_eq!(decoy_offset.to_string(), "PM[+15.9949]PEK");
    assert_eq!(decoy_offset.monoisotopic, offset.monoisotopic);
}

#[test]
fn mass_offsets_on_terminal_sites_stack_with_existing_terminal_mods() {
    use ModificationSpecificity::*;
    let acetylated = plain_peptide("PEPK")
        .apply(
            &[(PeptideN(None), 12.0f32, None)],
            &HashMap::default(),
            1,
            None,
        )
        .pop()
        .unwrap();
    // Unlike indexed terminal mods, a search-time mass offset explains an
    // extra delta, so it adds to whatever the terminus already carries.
    let offset =
        acetylated.with_mass_offset(Site::Nterm, &Arc::new(ModificationDefinition::bare(42.0)));
    assert_eq!(offset.nterm, Some(54.0));
    assert_eq!(offset.to_string(), "[+12][+42]-PEPK");
    assert!((offset.monoisotopic - (acetylated.monoisotopic + 42.0)).abs() < 1e-3);

    let named = plain_peptide("PEPK").with_mass_offset(
        Site::Cterm,
        &detailed_mod(0.984, "Amidated", &[], NeutralLossMode::Optional),
    );
    assert_eq!(named.cterm, Some(0.984));
    assert_eq!(named.to_string(), "PEPK-[Amidated]");
}

#[test]
fn modification_count_matches_site_rule_and_mass() {
    use ModificationSpecificity::*;
    let variants = plain_peptide("MSMK").apply(
        &[(Residue(b'M'), 16.0f32, None), (PeptideN(None), 42.0, None)],
        &HashMap::default(),
        3,
        None,
    );
    let full = variants
        .iter()
        .find(|peptide| peptide.to_string() == "[+42]-M[+16]SM[+16]K")
        .unwrap();
    assert_eq!(full.modification_count(Residue(b'M'), 16.0), 2);
    assert_eq!(full.modification_count(PeptideN(None), 42.0), 1);
    assert_eq!(full.modification_count(Residue(b'M'), 32.0), 0);
    assert_eq!(full.modification_count(PeptideN(Some(b'M')), 16.0), 1);
    assert_eq!(full.modification_count(Residue(b'S'), 16.0), 0);
    assert_eq!(variants[0].modification_count(Residue(b'M'), 16.0), 0);

    // Peptides parsed from mass deltas only carry terminal masses.
    let legacy = Peptide {
        sequence: (&b"AK"[..]).into(),
        nterm: Some(42.0),
        cterm: Some(1.0),
        ..Peptide::default()
    };
    assert_eq!(legacy.modification_count(PeptideN(None), 42.0), 1);
    assert_eq!(legacy.modification_count(PeptideC(None), 1.0), 1);
    assert_eq!(legacy.modification_count(PeptideC(None), 42.0), 0);
}

#[test]
fn relocating_a_modification_moves_only_the_matching_mass() {
    use ModificationSpecificity::*;
    let peptide = plain_peptide("SASK")
        .apply(
            &[(Residue(b'S'), 79.97f32, None), (Residue(b'K'), 8.0, None)],
            &HashMap::default(),
            2,
            None,
        )
        .into_iter()
        .find(|peptide| peptide.to_string() == "S[+79.97]ASK[+8]")
        .unwrap();

    let mut moved = peptide.clone();
    moved.relocate_modification_mass(79.97, &[0, 2], &[2], 0.01);
    assert_eq!(moved.to_string(), "SAS[+79.97]K[+8]");
    assert_eq!(moved.applied_modifications().len(), 2);

    // Index `len` addresses the N-terminus.
    let mut to_nterm = peptide.clone();
    to_nterm.relocate_modification_mass(79.97, &[0, 4], &[4], 0.01);
    assert_eq!(to_nterm.nterm, Some(79.97));
    assert_eq!(to_nterm.to_string(), "[+79.97]-SASK[+8]");

    // No applied modification has this mass: nothing changes.
    let mut unchanged = peptide.clone();
    unchanged.relocate_modification_mass(42.0, &[0, 2], &[2], 0.01);
    assert_eq!(unchanged.to_string(), peptide.to_string());
}

#[test]
fn peptide_errors_describe_the_rejected_input() {
    let error = Peptide::try_from(Digest {
        sequence: "PÉP".into(),
        ..Default::default()
    })
    .unwrap_err();
    assert_eq!(error.to_string(), "invalid peptide sequence: PÉP");
    let error = PeptideError::SequenceTooLong {
        length: 300,
        maximum: 255,
    };
    assert_eq!(
        error.to_string(),
        "peptide has 300 residues, but compact modification encoding supports at most 255"
    );
}

/// Enumerate acetylation of KAKAK (lysines at 0, 2 and 4) with a library
/// that supports position 0, returning the modified positions of each variant.
fn library_acetyl_variants(
    max_count: Option<usize>,
    max_total_count: Option<usize>,
    site_mode: SiteMode,
    max_exhaustive_mods: usize,
    max_total_mods: usize,
) -> Vec<Vec<u32>> {
    let rules = vec![VariableRule {
        specificity: ModificationSpecificity::Residue(b'K'),
        modification: detailed_mod(42.0106, "Acetyl", &[], NeutralLossMode::Optional),
        max_count,
        max_total_count,
        site_mode,
        count_group: 0,
    }];
    let library = vec![LibrarySite {
        attachment: Default::default(),
        position: 0,
        modification: Arc::from("Acetyl"),
    }];
    let static_mods = HashMap::new();
    let labels = LabelModificationCache::new(rules.iter().map(|rule| &rule.modification), &[]);
    let lookup = ModificationLookup::for_rules(&rules, &static_mods, &[], &labels).unwrap();
    let plan = ModificationPlan::new(
        &rules,
        &static_mods,
        lookup,
        max_exhaustive_mods,
        max_total_mods,
        None,
    );
    let mut variants = plain_peptide("KAKAK")
        .apply_rules(&plan, &library)
        .iter()
        .map(|peptide| {
            peptide
                .applied_modifications()
                .map(|applied| match applied.site {
                    Site::Sequence(index) => index,
                    site => panic!("unexpected site {site:?}"),
                })
                .collect::<Vec<_>>()
        })
        .collect::<Vec<_>>();
    variants.sort();
    variants
}

#[test]
fn library_placements_skip_the_exhaustive_budget_but_not_the_total_budget() {
    use SiteMode::*;
    // One new placement allowed; the library site at 0 is free.
    assert_eq!(
        library_acetyl_variants(None, None, Both, 1, 3),
        vec![vec![], vec![0], vec![0, 2], vec![0, 4], vec![2], vec![4]]
    );
    // The total budget still counts library placements.
    assert_eq!(
        library_acetyl_variants(None, None, Both, 3, 1),
        vec![vec![], vec![0], vec![2], vec![4]]
    );
    // No exhaustive budget: only library placements remain.
    assert_eq!(
        library_acetyl_variants(None, None, Both, 0, 3),
        vec![vec![], vec![0]]
    );
}

#[test]
fn max_count_limits_new_sites_and_max_total_count_limits_all_sites() {
    use SiteMode::*;
    // max_count 1 alone also caps the total at 1 (max_total_count defaults to it).
    assert_eq!(
        library_acetyl_variants(Some(1), None, Both, 3, 3),
        vec![vec![], vec![0], vec![2], vec![4]]
    );
    // Raising the total lets one library site join one new site.
    assert_eq!(
        library_acetyl_variants(Some(1), Some(2), Both, 3, 3),
        vec![vec![], vec![0], vec![0, 2], vec![0, 4], vec![2], vec![4]]
    );
    // A total cap without max_count bounds every placement.
    assert_eq!(
        library_acetyl_variants(None, Some(2), Both, 3, 3).len(),
        1 + 3 + 3
    );
    // max_count 1 with a total of 3 still forbids two new sites.
    let variants = library_acetyl_variants(Some(1), Some(3), Both, 3, 3);
    assert!(!variants.contains(&vec![0, 2, 4]));
    assert!(!variants.contains(&vec![2, 4]));
    assert!(variants.contains(&vec![0, 4]));
}

#[test]
fn site_mode_decides_which_candidates_exist() {
    use SiteMode::*;
    // Library mode keeps only library-supported sites.
    assert_eq!(
        library_acetyl_variants(Some(1), None, Library, 3, 3),
        vec![vec![], vec![0]]
    );
    // Exhaustive mode ignores the library, so site 0 consumes max_count.
    assert_eq!(
        library_acetyl_variants(Some(1), Some(3), Exhaustive, 3, 3),
        vec![vec![], vec![0], vec![2], vec![4]]
    );
}
