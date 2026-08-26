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
