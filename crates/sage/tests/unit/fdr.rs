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
