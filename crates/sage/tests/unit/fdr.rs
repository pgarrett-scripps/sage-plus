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

fn target(ix: u32, score: f32) -> Competition<u32> {
    Competition {
        forward: score,
        foward_ix: Some(ix),
        ..Default::default()
    }
}

fn pair(target: u32, forward: f32, decoy: u32, reverse: f32) -> Competition<u32> {
    Competition {
        forward,
        foward_ix: Some(target),
        reverse,
        reverse_ix: Some(decoy),
    }
}

#[test]
fn picked_q_values_count_winners_with_plus_one() {
    let mut scores = FnvHashMap::default();
    for ix in 0..200_u32 {
        scores.insert(ix, target(ix, 10.0 + ix as f32 / 200.0));
    }
    scores.insert(
        200,
        Competition {
            reverse: 1.0,
            reverse_ix: Some(200_u32),
            ..Default::default()
        },
    );
    let (q, passing) = Competition::assign_q_value(scores, 0.01);
    // (0 + 1) / 200 for every target; the decoy adds one: (1 + 1) / 200.
    assert_eq!(passing, 200);
    assert_eq!(q[&0], 0.005);
    assert_eq!(q[&199], 0.005);
    assert_eq!(q[&200], 0.01);
}

#[test]
fn picked_q_values_discard_the_losing_member() {
    // Pairs (target ix, target score, decoy ix, decoy score):
    //   A: 10 vs 1   -> target wins
    //   B:  9 vs 9.5 -> decoy wins; target 9 is discarded
    //   C:  8 vs 8   -> tie goes to the decoy
    //   D:  7 vs 6.9 -> target wins; its decoy 6.9 is discarded
    //   E:  6 alone  -> target wins
    // Ranked winners: A 10 (t), B' 9.5 (d), C' 8 (d), D 7 (t), E 6 (t)
    // raw (d + 1) / t: 1/1, 2/1 -> 1, 3/1 -> 1, 3/2 -> 1, 3/3 = 1
    // Everything is capped at 1 in this tiny example, so use a larger one below.
    let mut scores = FnvHashMap::default();
    scores.insert("A", pair(0, 10.0, 100, 1.0));
    scores.insert("B", pair(1, 9.0, 101, 9.5));
    scores.insert("C", pair(2, 8.0, 102, 8.0));
    scores.insert("D", pair(3, 7.0, 103, 6.9));
    scores.insert("E", target(4, 6.0));
    let (q, passing) = Competition::assign_q_value(scores, 0.5);
    assert_eq!(passing, 0);
    for ix in [0, 1, 2, 3, 4, 100, 101, 102, 103] {
        assert_eq!(q[&ix], 1.0, "{ix}");
    }

    // Ten targets at 20..11 beat their decoys at 0; target 10 loses to a
    // decoy at 15.5; one unpaired decoy at 10.5; target 11 at 10 is alone.
    // Ranked winners: t20..t16 (5), d15.5, t15..t11 (5), d10.5, t10.
    //   through t16: (0 + 1) / 5 = 0.2
    //   through d15.5 and t11: (1 + 1) / 10 = 0.2
    //   through d10.5 and t10: (2 + 1) / 11
    let mut scores = FnvHashMap::default();
    for i in 0..10_u32 {
        scores.insert(i, pair(i, 20.0 - i as f32, 100 + i, 0.0));
    }
    scores.insert(10, pair(10, 12.5, 110, 15.5));
    scores.insert(11, target(11, 10.0));
    scores.insert(
        12,
        Competition {
            reverse: 10.5,
            reverse_ix: Some(112),
            ..Default::default()
        },
    );
    let (q, passing) = Competition::assign_q_value(scores, 0.2);
    assert_eq!(passing, 10);
    for i in 0..10_u32 {
        assert!((q[&i] - 0.2).abs() < 1e-6, "target {i}: {}", q[&i]);
        // Losing decoys are not ranked at all.
        assert_eq!(q[&(100 + i)], 1.0);
    }
    // The winning decoy is ranked; the target it beat is discarded, even
    // though its own score (12.5) is above the 1-in-5 threshold.
    assert!((q[&110] - 0.2).abs() < 1e-6);
    assert_eq!(q[&10], 1.0);
    assert!((q[&112] - 3.0 / 11.0).abs() < 1e-6);
    assert!((q[&11] - 3.0 / 11.0).abs() < 1e-6);
}

#[test]
fn picked_q_values_keep_tied_winners_together() {
    // Two targets and a decoy winner share score 5: the tie group is counted
    // as a whole, so both targets see the decoy: (1 + 1) / 3.
    let mut scores = FnvHashMap::default();
    scores.insert(0, target(0, 9.0));
    scores.insert(1, target(1, 5.0));
    scores.insert(2, target(2, 5.0));
    scores.insert(3, pair(3, 1.0, 13, 5.0));
    let (q, _) = Competition::assign_q_value(scores, 0.01);
    for ix in [0, 1, 2, 13] {
        assert!((q[&ix] - 2.0 / 3.0).abs() < 1e-6, "{ix}: {}", q[&ix]);
    }
    assert_eq!(q[&3], 1.0);
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

fn scored(peptide_idx: u32, score: f32) -> Feature {
    Feature {
        peptide_idx: PeptideIx(peptide_idx),
        discriminant_score: score,
        protein_q: 42.0,
        ..Default::default()
    }
}

#[test]
fn picked_protein_counts_unique_proteins_and_skips_shared_peptides() {
    let mut shared = peptide("SHAREDK", false);
    shared.proteins = [Arc::<str>::from("P1"), Arc::<str>::from("P2")]
        .into_iter()
        .collect();
    let db = IndexedDatabase {
        peptides: vec![
            protein_peptide("AAAAK", false, "P1"),
            protein_peptide("CCCCK", false, "P2"),
            protein_peptide("DDDDK", false, "P3"),
            protein_peptide("EEEEK", false, "P4"),
            protein_peptide("KFFFF", true, "P5"),
            protein_peptide("KAAAA", true, "P1"),
            shared,
            // A second, weaker peptide of P1 must not change P1's score.
            protein_peptide("GGGGK", false, "P1"),
        ],
        decoy_tag: "rev_".into(),
        generate_decoys: true,
        ..Default::default()
    };
    let mut features = vec![
        scored(0, 10.0),
        scored(1, 9.0),
        scored(2, 8.0),
        scored(3, 7.0),
        scored(4, 1.0),
        scored(5, 0.5),
        scored(6, 50.0),
        scored(7, 0.1),
    ];
    // Picked protein scores: P1 10 beats rev_P1 0.5, which is discarded.
    // Ranked winners: P1 10, P2 9, P3 8, P4 7 (targets), rev_P5 1 (decoy)
    // raw q = (1 + decoys) / targets: 1, 1/2, 1/3, 1/4, 2/4
    // after the reverse cumulative minimum: 1/4 x4, 1/2.
    let passing = picked_protein(&db, &mut features);
    assert_eq!(passing, 0);
    for feat in &features[0..4] {
        assert!((feat.protein_q - 0.25).abs() < 1e-6, "{}", feat.protein_q);
    }
    assert!((features[4].protein_q - 0.5).abs() < 1e-6);
    // The losing decoy of P1 is not ranked.
    assert_eq!(features[5].protein_q, 1.0);
    // Shared peptides are excluded and keep their prior value.
    assert_eq!(features[6].protein_q, 42.0);
    // Every peptide of a protein inherits the protein's q-value.
    assert!((features[7].protein_q - 0.25).abs() < 1e-6);
}

#[test]
fn picked_precursor_assigns_count_q_values_per_charge_state() {
    let quantified = |score: f64| QuantifiedPeak {
        peak: crate::lfq::Peak {
            score,
            q_value: 42.0,
            ..Default::default()
        },
        intensities: Vec::new(),
        ms2_confirmed: Vec::new(),
        ms2_confirmed_strict: Vec::new(),
        file_evidence: Vec::new(),
    };
    let mut peaks = FnvHashMap::default();
    // 20 targets at 100..81, one decoy at 50, one target at 40.
    for i in 0..20u32 {
        peaks.insert(
            (PrecursorId::Charged((PeptideIx(i), 2)), false),
            quantified(100.0 - i as f64),
        );
    }
    peaks.insert(
        (PrecursorId::Charged((PeptideIx(0), 2)), true),
        quantified(50.0),
    );
    peaks.insert(
        (PrecursorId::Combined(PeptideIx(99)), false),
        quantified(40.0),
    );

    // Raw q at the 20th target is (1 + 0) / 20 = 0.05, which passes the 5%
    // threshold; the decoy raises it to 2/20 and the last target to 2/21.
    let passing = picked_precursor(&mut peaks);
    assert_eq!(passing, 20);
    for i in 0..20u32 {
        let q = peaks[&(PrecursorId::Charged((PeptideIx(i), 2)), false)]
            .peak
            .q_value;
        assert!((q - 0.05).abs() < 1e-6, "target {i}: {q}");
    }
    let decoy = peaks[&(PrecursorId::Charged((PeptideIx(0), 2)), true)]
        .peak
        .q_value;
    let last = peaks[&(PrecursorId::Combined(PeptideIx(99)), false)]
        .peak
        .q_value;
    assert!((decoy - 2.0 / 21.0).abs() < 1e-6, "{decoy}");
    assert!((last - 2.0 / 21.0).abs() < 1e-6, "{last}");
}

#[test]
fn picked_peptide_reports_a_target_beaten_by_its_decoy_with_q_one() {
    let mut peptides = Vec::new();
    let mut features = Vec::new();
    // Twenty confident target peptides whose decoys score lower, then one
    // target (index 40) whose own reversed decoy outscores it.
    for i in 0..21_u32 {
        let target = peptide(&format!("PEPT{}IDEK", "A".repeat(i as usize)), false);
        let decoy = target.reverse();
        assert!(decoy.decoy);
        peptides.push(target);
        peptides.push(decoy);
        let (target_score, decoy_score) = match i {
            20 => (5.0, 6.0),
            _ => (100.0 - i as f32, 1.0),
        };
        features.push(Feature {
            peptide_idx: PeptideIx(2 * i),
            discriminant_score: target_score,
            ..Default::default()
        });
        features.push(Feature {
            peptide_idx: PeptideIx(2 * i + 1),
            discriminant_score: decoy_score,
            ..Default::default()
        });
    }
    let db = IndexedDatabase {
        peptides,
        generate_decoys: true,
        ..Default::default()
    };
    // Winners: 20 targets (100..81), then the decoy at 6.
    // q = (0 + 1) / 20 = 0.05 for the targets; the winning decoy (1 + 1) / 20.
    let passing = picked_peptide(&db, &mut features);
    assert_eq!(passing, 0);
    for i in 0..20 {
        assert!((features[2 * i].peptide_q - 0.05).abs() < 1e-6);
        assert_eq!(features[2 * i + 1].peptide_q, 1.0);
    }
    assert_eq!(features[40].peptide_q, 1.0);
    assert!((features[41].peptide_q - 0.1).abs() < 1e-6);
}
