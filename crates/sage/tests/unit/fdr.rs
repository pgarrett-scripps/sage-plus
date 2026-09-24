use super::*;
use crate::enzyme::Digest;
use crate::modification::{ModificationDefinition, NeutralLossMode};
use crate::peptide::{AppliedModification, CompactModifications, ModificationKind, Peptide, Site};
use std::sync::Arc;

fn peptide(sequence: &str, decoy: bool) -> Peptide {
    let mut peptide = Peptide::try_from(Digest {
        sequence: sequence.into(),
        ..Default::default()
    })
    .unwrap();
    peptide.decoy = decoy;
    peptide
}

#[test]
fn picked_peptide_assigns_one_to_orphaned_competition_twins() {
    let mut twin_a = peptide("PEPTIDE", false);
    twin_a.modifications = CompactModifications::from_applied([AppliedModification {
        site: Site::Sequence(1),
        modification: Arc::new(ModificationDefinition {
            mass: 10.0,
            name: None,
            neutral_losses: Arc::from([5.0]),
            neutral_loss_mode: NeutralLossMode::Optional,
            channel_offsets: Arc::default(),
        }),
        kind: ModificationKind::Ordinary,
    }])
    .unwrap();
    let mut twin_b = twin_a.clone();
    twin_b.modifications = CompactModifications::from_applied([AppliedModification {
        site: Site::Sequence(1),
        modification: Arc::new(ModificationDefinition {
            mass: 10.0,
            name: None,
            neutral_losses: Arc::from([6.0]),
            neutral_loss_mode: NeutralLossMode::Optional,
            channel_offsets: Arc::default(),
        }),
        kind: ModificationKind::Ordinary,
    }])
    .unwrap();

    assert_eq!(twin_a.to_string(), twin_b.to_string());
    let mut twins = vec![twin_a.clone(), twin_b.clone()];
    crate::database::Parameters::reorder_peptides(&mut twins);
    assert_eq!(twins.len(), 2);

    let db = IndexedDatabase {
        peptides: vec![
            twin_a,
            twin_b,
            peptide("AAAAA", false),
            peptide("CCCCC", true),
            peptide("GGGGG", true),
        ],
        generate_decoys: false,
        ..Default::default()
    };
    let mut features = [
        Feature {
            peptide_idx: PeptideIx(0),
            discriminant_score: 10.0,
            ..Default::default()
        },
        Feature {
            peptide_idx: PeptideIx(1),
            discriminant_score: 9.0,
            ..Default::default()
        },
        Feature {
            peptide_idx: PeptideIx(2),
            discriminant_score: 8.0,
            ..Default::default()
        },
        Feature {
            peptide_idx: PeptideIx(3),
            discriminant_score: 7.0,
            ..Default::default()
        },
        Feature {
            peptide_idx: PeptideIx(4),
            discriminant_score: 2.0,
            ..Default::default()
        },
    ];

    picked_peptide(&db, &mut features);

    assert_eq!(features[0].peptide_q, 1.0);
}

#[test]
fn sparse_decoys_use_finite_count_based_confidence() {
    let mut scores = FnvHashMap::default();
    for ix in 0..200_u32 {
        scores.insert(
            ix,
            Competition {
                forward: 10.0 + ix as f32 / 200.0,
                foward_ix: Some(ix),
                ..Default::default()
            },
        );
    }
    scores.insert(
        200,
        Competition {
            reverse: 1.0,
            reverse_ix: Some(200_u32),
            ..Default::default()
        },
    );
    assert!(Competition::fit_kde(&scores).is_none());
    let (q, passing) = Competition::assign_q_value(scores, 0.01);
    assert_eq!(passing, 200);
    assert_eq!(q[&0], 0.005);
    assert_eq!(q[&199], 0.005);
    assert_eq!(q[&200], 0.01);
}

#[test]
fn count_confidence_is_identical_for_ties_and_keeps_plus_one() {
    let mut rows = vec![
        Row {
            ix: 0,
            decoy: false,
            score: 2.0,
            q: 1.0,
        },
        Row {
            ix: 1,
            decoy: false,
            score: 1.0,
            q: 1.0,
        },
        Row {
            ix: 2,
            decoy: true,
            score: 1.0,
            q: 1.0,
        },
    ];
    assign_count_q_values(&mut rows, 0.01);
    assert!(rows.iter().all(|row| row.q == 1.0));
    rows.reverse();
    assign_count_q_values(&mut rows, 0.01);
    assert!(rows.iter().all(|row| row.q == 1.0));
}

#[test]
fn empty_and_decoy_only_count_confidence_are_conservative() {
    let mut empty: Vec<Row<u32>> = Vec::new();
    assert_eq!(assign_count_q_values(&mut empty, 0.01), 0);
    let mut rows = vec![Row {
        ix: 0,
        decoy: true,
        score: 1.0,
        q: 0.0,
    }];
    assert_eq!(assign_count_q_values(&mut rows, 0.01), 0);
    assert_eq!(rows[0].q, 1.0);
}

fn protein_peptide(sequence: &str, decoy: bool, protein: &str) -> Peptide {
    let mut peptide = peptide(sequence, decoy);
    peptide.proteins = [Arc::<str>::from(protein)].into_iter().collect();
    peptide
}

fn grouped(peptide_idx: u32, groups: &str, score: f32) -> Feature {
    Feature {
        peptide_idx: PeptideIx(peptide_idx),
        protein_groups: Some(groups.into()),
        num_protein_groups: 1,
        discriminant_score: score,
        ..Default::default()
    }
}

#[test]
fn decoy_protein_groups_compete_with_their_target_group() {
    for generate_decoys in [true, false] {
        let decoy_accession = |protein: &str| match generate_decoys {
            true => protein.to_string(),
            false => format!("rev_{protein}"),
        };
        let db = IndexedDatabase {
            peptides: vec![
                protein_peptide("PEPTIDEK", false, "P1"),
                protein_peptide("KEDITPEP", true, &decoy_accession("P1")),
                protein_peptide("AAAAAK", false, "P3"),
                protein_peptide("KAAAAA", true, &decoy_accession("P9")),
            ],
            decoy_tag: "rev_".into(),
            generate_decoys,
            ..Default::default()
        };
        // Targets carry parsimony groups; decoys keep their raw decoy accession.
        let features = [
            grouped(0, "P1/P2", 10.0),
            grouped(1, "rev_P1", 4.0),
            grouped(2, "P3", 8.0),
            grouped(3, "rev_P9", 3.0),
        ];
        let target_groups = target_protein_groups(&db, &features);
        assert_eq!(
            decoy_competition_group(&db, &features[1], &target_groups).as_deref(),
            Some("P1/P2")
        );
        // A decoy whose target protein has no reported group stays unpaired.
        assert_eq!(
            decoy_competition_group(&db, &features[3], &target_groups).as_deref(),
            Some("rev_P9")
        );

        // A second decoy accession reversed from the same target group.
        let mut db = db;
        db.peptides
            .push(protein_peptide("KEDITPEPP", true, &decoy_accession("P2")));
        let mut features = features.to_vec();
        features.push(grouped(4, "rev_P2", 5.0));
        picked_protein_group(&db, &mut features);
        assert!(features.iter().all(|feat| feat.protein_group_q <= 1.0));
        // Both decoys of the P1/P2 group share the competition's decoy q-value.
        assert_eq!(features[1].protein_group_q, features[4].protein_group_q);
    }
}

#[test]
fn decoys_pair_with_the_best_supported_group_of_their_target() {
    let db = IndexedDatabase {
        peptides: vec![
            protein_peptide("PEPTIDEK", false, "P1"),
            protein_peptide("LLLLLK", false, "P1"),
            protein_peptide("KEDITPEP", true, "P1"),
            protein_peptide("GGGGGK", false, "P5"),
        ],
        decoy_tag: "rev_".into(),
        generate_decoys: true,
        ..Default::default()
    };
    // A low-confidence peptide can leave P1 with a fallback group of its own;
    // the decoy still competes against the group with P1's real evidence.
    // P5 is not a member of any reported group string, so it is not paired.
    let features = [
        grouped(0, "P1/P2", 10.0),
        grouped(1, "P1", 0.5),
        grouped(2, "rev_P1", 4.0),
        grouped(3, "P1/P2", 9.0),
    ];
    let target_groups = target_protein_groups(&db, &features);
    assert_eq!(target_groups.get("P1"), Some(&"P1/P2"));
    assert_eq!(target_groups.get("P2"), Some(&"P1/P2"));
    assert_eq!(target_groups.get("P5"), None);
    assert_eq!(
        decoy_competition_group(&db, &features[2], &target_groups).as_deref(),
        Some("P1/P2")
    );
}

#[test]
fn fasta_decoy_tags_outside_the_prefix_still_pair() {
    let db = IndexedDatabase {
        peptides: vec![
            protein_peptide("PEPTIDEK", false, "sp|P1|X"),
            protein_peptide("KEDITPEP", true, "sp|rev_P1|X"),
        ],
        decoy_tag: "rev_".into(),
        generate_decoys: false,
        ..Default::default()
    };
    let features = [grouped(0, "sp|P1|X", 10.0), grouped(1, "sp|rev_P1|X", 4.0)];
    let target_groups = target_protein_groups(&db, &features);
    assert_eq!(
        decoy_competition_group(&db, &features[1], &target_groups).as_deref(),
        Some("sp|P1|X")
    );
}
