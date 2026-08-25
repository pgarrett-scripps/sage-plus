use super::*;
use std::sync::Arc;

use crate::peptide::Peptide;

fn get_data() -> (Vec<Vec<&'static str>>, Vec<bool>, Vec<f32>) {
    let proteins = vec![
        vec!["protein_7"],
        vec!["protein_4", "protein_6", "protein_9"],
        vec!["protein_1"],
        vec!["protein_1", "protein_5"],
        vec!["protein_7"],
        vec!["protein_3", "protein_6"],
        vec!["protein_1"],
        vec!["protein_1", "protein_2", "protein_5", "protein_8"],
        vec!["protein_1"],
        vec!["protein_4", "protein_9"],
    ];
    let decoys = vec![false; proteins.len()];
    let q_vals = vec![0.0; proteins.len()];
    (proteins, decoys, q_vals)
}

fn build_db_and_features(
    proteins: &[Vec<&str>],
    decoys: &[bool],
    q_vals: &[f32],
) -> (IndexedDatabase, Vec<Feature>) {
    build_db_and_features_with_scores(proteins, decoys, q_vals, None)
}

fn build_db_and_features_with_scores(
    proteins: &[Vec<&str>],
    decoys: &[bool],
    q_vals: &[f32],
    scores: Option<&[f32]>,
) -> (IndexedDatabase, Vec<Feature>) {
    let features: Vec<Feature> = (0..proteins.len())
        .map(|ix| Feature {
            peptide_idx: PeptideIx(ix as u32),
            label: if decoys[ix] { -1 } else { 1 },
            peptide_q: q_vals[ix],
            discriminant_score: scores.map(|s| s[ix]).unwrap_or(0.0),
            ..Default::default()
        })
        .collect();
    let db = IndexedDatabase {
        peptides: proteins
            .iter()
            .zip(decoys)
            .map(|(prots, &decoy)| Peptide {
                proteins: prots.iter().map(|&s| Arc::from(s)).collect(),
                decoy,
                ..Default::default()
            })
            .collect(),
        decoy_tag: "rev_".to_string(),
        generate_decoys: false,
        ..IndexedDatabase::default()
    };
    (db, features)
}

#[test]
fn test_protein_grouping_expected_groups() {
    let (proteins, decoys, q_vals) = get_data();
    let (db, mut features) = build_db_and_features(&proteins, &decoys, &q_vals);
    generate_protein_groups(&db, &mut features, true, Some(0.01));
    let expected = vec![
        "protein_7",
        "protein_4/protein_9;protein_6",
        "protein_1",
        "protein_1",
        "protein_7",
        "protein_6",
        "protein_1",
        "protein_1",
        "protein_1",
        "protein_4/protein_9",
    ];
    let actual: Vec<_> = features
        .iter()
        .map(|v| v.protein_groups.as_ref().unwrap().as_str())
        .collect();
    assert_eq!(actual, expected);
}

#[test]
fn test_bipartite_cover_unique_peptides() {
    // Three proteins each with a single unique peptide — all in cover
    let edges = vec![(0, 0), (1, 1), (2, 2)];
    let cover = BipartiteGraph::new(edges, 3, 3).into_cover();
    assert_eq!(cover, vec![true, true, true]);
}

#[test]
fn test_bipartite_cover_subset_protein() {
    // protein 0 -> peptides {0, 1, 2}  (superset)
    // protein 1 -> peptides {0, 1}     (subset)
    let edges = vec![(0, 0), (0, 1), (0, 2), (1, 0), (1, 1)];
    let cover = BipartiteGraph::new(edges, 2, 3).into_cover();
    assert!(cover[0], "superset protein should be covered");
    assert!(!cover[1], "subset protein should not be covered");
}

#[test]
fn test_bipartite_cover_shared_peptide() {
    // protein 0 -> peptides {0, 1}
    // protein 1 -> peptides {1, 2}
    // Both should be in cover because each has a unique peptide
    let edges = vec![(0, 0), (0, 1), (1, 1), (1, 2)];
    let cover = BipartiteGraph::new(edges, 2, 3).into_cover();
    assert_eq!(cover, vec![true, true]);
}

#[test]
fn test_bipartite_cover_empty() {
    let cover = BipartiteGraph::new(vec![], 0, 0).into_cover();
    assert!(cover.is_empty());
}

#[test]
fn test_bipartite_cover_single() {
    let cover = BipartiteGraph::new(vec![(0, 0)], 1, 1).into_cover();
    assert_eq!(cover, vec![true]);
}

#[test]
fn test_decoy_features_excluded_from_grouping() {
    let proteins = vec![vec!["protA"], vec!["protA"], vec!["protB"]];
    let decoys = vec![false, true, false];
    let q_vals = vec![0.0, 0.0, 0.0];
    let (db, mut features) = build_db_and_features(&proteins, &decoys, &q_vals);
    generate_protein_groups(&db, &mut features, true, Some(0.01));

    for feat in &features {
        assert!(
            feat.protein_groups.is_some(),
            "every feature should be annotated"
        );
    }
    assert_eq!(features[1].protein_groups.as_deref(), Some("protA"));
}

#[test]
fn test_decoy_features_with_generate_decoys() {
    let proteins = vec![vec!["protA"], vec!["protA"]];
    let decoys = vec![false, true];
    let q_vals = vec![0.0, 0.0];
    let (db, mut features) = build_db_and_features(&proteins, &decoys, &q_vals);

    // Override generate_decoys for this test
    let db = IndexedDatabase {
        generate_decoys: true,
        ..db
    };

    generate_protein_groups(&db, &mut features, false, None);
    assert_eq!(features[0].protein_groups.as_deref(), Some("protA"));
    assert_eq!(features[1].protein_groups.as_deref(), Some("rev_protA"));
}

#[test]
fn test_grouping_disabled_falls_back_to_annotate() {
    let proteins = vec![vec!["protA", "protB"], vec!["protC"]];
    let decoys = vec![false, false];
    let q_vals = vec![0.0, 0.0];
    let (db, mut features) = build_db_and_features(&proteins, &decoys, &q_vals);
    generate_protein_groups(&db, &mut features, false, None);

    assert_eq!(features[0].protein_groups.as_deref(), Some("protA;protB"));
    assert_eq!(features[0].num_protein_groups, 2);
    assert_eq!(features[1].protein_groups.as_deref(), Some("protC"));
    assert_eq!(features[1].num_protein_groups, 1);
}

#[test]
fn test_single_protein_single_peptide() {
    let proteins = vec![vec!["protA"]];
    let decoys = vec![false];
    let q_vals = vec![0.0];
    let (db, mut features) = build_db_and_features(&proteins, &decoys, &q_vals);
    generate_protein_groups(&db, &mut features, true, Some(0.01));

    assert_eq!(features[0].protein_groups.as_deref(), Some("protA"));
    assert_eq!(features[0].num_protein_groups, 1);
}

#[test]
fn test_all_shared_peptides() {
    // Every peptide maps to the same two proteins — they should form a single group
    let proteins = vec![
        vec!["protA", "protB"],
        vec!["protA", "protB"],
        vec!["protA", "protB"],
    ];
    let decoys = vec![false, false, false];
    let q_vals = vec![0.0, 0.0, 0.0];
    let (db, mut features) = build_db_and_features(&proteins, &decoys, &q_vals);
    generate_protein_groups(&db, &mut features, true, Some(0.01));

    let group = features[0].protein_groups.as_deref().unwrap();
    assert!(group.contains("protA"));
    assert!(group.contains("protB"));
    for feat in &features {
        assert_eq!(feat.protein_groups.as_deref(), Some(group));
        assert_eq!(feat.num_protein_groups, 1);
    }
}

#[test]
fn test_peptide_fdr_threshold_filtering() {
    let proteins = vec![vec!["protA"], vec!["protB"]];
    let decoys = vec![false, false];
    let q_vals = vec![0.001, 0.5];
    let (db, mut features) = build_db_and_features(&proteins, &decoys, &q_vals);
    generate_protein_groups(&db, &mut features, true, Some(0.01));

    assert!(features[0].protein_groups.is_some());
    assert!(features[1].protein_groups.is_some());
}

#[test]
fn test_all_decoy_features() {
    let proteins = vec![vec!["protA"], vec!["protB"]];
    let decoys = vec![true, true];
    let q_vals = vec![0.0, 0.0];
    let (db, mut features) = build_db_and_features(&proteins, &decoys, &q_vals);
    generate_protein_groups(&db, &mut features, true, Some(0.01));

    for feat in &features {
        assert!(feat.protein_groups.is_some());
    }
}

#[test]
fn test_proteins_with_identical_evidence_are_grouped() {
    let proteins = vec![
        vec!["protA", "protB"],
        vec!["protA", "protB"],
        vec!["protC"],
    ];
    let decoys = vec![false, false, false];
    let q_vals = vec![0.0, 0.0, 0.0];
    let (db, mut features) = build_db_and_features(&proteins, &decoys, &q_vals);
    generate_protein_groups(&db, &mut features, true, Some(0.01));

    let group_01 = features[0].protein_groups.as_deref().unwrap();
    assert!(group_01.contains("protA") && group_01.contains("protB"));
    assert_eq!(features[0].protein_groups, features[1].protein_groups);
    assert_eq!(features[2].protein_groups.as_deref(), Some("protC"));
}

#[test]
fn test_num_protein_groups_counts_distinct_groups() {
    // protA and protB each have unique evidence; peptide 2 is shared across groups
    let proteins = vec![vec!["protA"], vec!["protB"], vec!["protA", "protB"]];
    let decoys = vec![false, false, false];
    let q_vals = vec![0.0, 0.0, 0.0];
    let (db, mut features) = build_db_and_features(&proteins, &decoys, &q_vals);
    generate_protein_groups(&db, &mut features, true, Some(0.01));

    assert_eq!(features[0].num_protein_groups, 1);
    assert_eq!(features[1].num_protein_groups, 1);
    assert_eq!(features[2].num_protein_groups, 2);
}
