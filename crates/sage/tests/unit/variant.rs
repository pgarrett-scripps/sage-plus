use super::*;
use crate::fasta::{Fasta, ProteinClass};
use crate::peptide::Peptide;
use std::sync::Arc;

/// A database of `reference` reference peptides on P1, then `variant` variant
/// peptides on the variant protein P1_V, then one peptide of both, each
/// followed by its decoy, so peptide `2 * i` is a target and `2 * i + 1` its
/// decoy. Sequences are distinct, so every pair is its own competition.
fn database(reference: usize, variant: usize) -> (IndexedDatabase, ProteinClasses) {
    let mut fasta = Fasta::parse(">P1\nPEPTIDEK\n".into(), "rev_", true).unwrap();
    let variants = Fasta::parse(">P1_V\nPEPTIDEH\n".into(), "rev_", true).unwrap();
    fasta.extend_with_class(variants, ProteinClass::Variant);
    let classes = ProteinClasses::of_fasta(&fasta).unwrap();

    let mut db = IndexedDatabase::default();
    let mut push = |index: usize, proteins: &[&str]| {
        // A distinct sequence per peptide: residues spell out the index.
        let sequence = format!("PEP{}K", "AG".repeat(index % 7 + 1)).into_bytes();
        let mut sequence = sequence;
        sequence.extend(index.to_string().bytes().map(|digit| b'C' + (digit - b'0')));
        for decoy in [false, true] {
            db.peptides.push(Peptide {
                sequence: sequence.as_slice().into(),
                proteins: proteins.iter().map(|&p| Arc::from(p)).collect(),
                decoy,
                ..Peptide::default()
            });
        }
    };
    for index in 0..reference {
        push(index, &["P1"]);
    }
    for index in reference..reference + variant {
        push(index, &["P1_V"]);
    }
    push(reference + variant, &["P1", "P1_V"]);
    (db, classes)
}

fn psm(psm_id: usize, peptide: usize, decoy: bool, score: f32) -> Feature {
    Feature {
        psm_id,
        peptide_idx: PeptideIx((2 * peptide + usize::from(decoy)) as u32),
        label: if decoy { -1 } else { 1 },
        discriminant_score: score,
        spectrum_q: 1.0,
        peptide_q: 1.0,
        ..Feature::default()
    }
}

#[test]
fn reference_decoys_do_not_enter_the_variant_q_value() {
    // 100 variant targets and no variant decoys, below ten reference decoys
    // that outscore everything. Counted globally, those decoys would put every
    // variant PSM at q >= 11/100; within the class the q-value is 1/100.
    let (db, classes) = database(10, 100);
    let mut features = (0..10)
        .map(|index| psm(index, index, true, 1000.0))
        .collect::<Vec<_>>();
    features.extend((0..100).map(|index| psm(10 + index, 10 + index, false, 50.0 + index as f32)));

    let (q, psms, peptides) = variant_q_values(&db, &features, &classes);
    for feature in &features[10..] {
        assert_eq!(q.spectrum_q(feature.psm_id), Some(0.01));
        assert_eq!(q.peptide_q(feature.psm_id), Some(0.01));
    }
    for feature in &features[..10] {
        assert_eq!(q.spectrum_q(feature.psm_id), None);
    }
    assert_eq!(
        psms,
        VariantClassCounts {
            targets: 100,
            decoys: 0,
            passing_variant_q: 100,
            passing_global_q: 0,
        }
    );
    assert_eq!(peptides.targets, 100);
    assert_eq!(peptides.passing_variant_q, 100);
}

#[test]
fn a_class_without_decoys_cannot_pass_below_one_over_its_targets() {
    let (db, classes) = database(0, 50);
    let features = (0..50)
        .map(|index| psm(index, index, false, 10.0 + index as f32))
        .collect::<Vec<_>>();
    let (q, psms, _) = variant_q_values(&db, &features, &classes);
    for feature in &features {
        assert_eq!(q.spectrum_q(feature.psm_id), Some(0.02));
    }
    assert_eq!(psms.passing_variant_q, 0);
}

#[test]
fn a_variant_decoy_is_in_the_class_and_a_shared_peptide_is_not() {
    // Variant peptide 0 has a target PSM and a decoy PSM that outscores it;
    // peptide 1 is on both the reference and the variant protein.
    let (db, classes) = database(0, 1);
    let features = vec![
        psm(0, 0, false, 5.0),
        psm(1, 0, true, 6.0),
        psm(2, 1, false, 100.0),
    ];
    let (q, psms, peptides) = variant_q_values(&db, &features, &classes);
    assert_eq!(psms.targets, 1);
    assert_eq!(psms.decoys, 1);
    // The decoy wins the picked competition, so the target loses with q 1.
    assert_eq!(q.peptide_q(0), Some(1.0));
    assert_eq!((peptides.targets, peptides.decoys), (1, 1));
    // A reference protein also explains peptide 1: no variant q-value.
    assert_eq!(q.spectrum_q(2), None);
    assert_eq!(q.peptide_q(2), None);
}

#[test]
fn protein_inference_sees_no_variant_psm_and_the_order_is_kept() {
    let (db, classes) = database(2, 2);
    let mut features = vec![
        psm(0, 0, false, 9.0),
        psm(1, 2, false, 8.0),
        psm(2, 1, false, 7.0),
        psm(3, 3, true, 6.0),
        psm(4, 4, false, 5.0),
    ];
    for feature in &mut features {
        feature.protein_q = 0.5;
    }
    let seen = infer_proteins_without_variants(&db, &classes, &mut features, |features| {
        for feature in features.iter_mut() {
            feature.protein_q = 0.001;
        }
        features.iter().map(|f| f.psm_id).collect::<Vec<_>>()
    });
    // Peptides 2 and 3 are variant peptides; 4 is shared, so a reference one.
    assert_eq!(seen, vec![0, 2, 4]);
    assert_eq!(
        features.iter().map(|f| f.psm_id).collect::<Vec<_>>(),
        vec![0, 1, 2, 3, 4]
    );
    for feature in &features {
        match feature.psm_id {
            1 | 3 => {
                assert_eq!(feature.protein_q, 1.0);
                assert_eq!(feature.protein_group_q, 1.0);
                assert_eq!(feature.protein_groups.as_deref(), Some("P1_V"));
                assert_eq!(feature.num_protein_groups, 1);
            }
            _ => assert_eq!(feature.protein_q, 0.001),
        }
    }
}

#[test]
fn without_variant_psms_protein_inference_sees_every_psm() {
    let (db, classes) = database(2, 0);
    let mut features = vec![psm(0, 0, false, 9.0), psm(1, 1, true, 8.0)];
    let seen =
        infer_proteins_without_variants(&db, &classes, &mut features, |features| features.len());
    assert_eq!(seen, 2);
}

#[test]
fn a_variant_decoy_cannot_raise_the_reference_proteins_q_values() {
    // 200 reference proteins, each with one target PSM, and a variant protein
    // whose decoy outscores them all. Counted in the picked protein FDR, that
    // decoy makes every reference protein's q-value 2/200; kept out, 1/200.
    let mut db = IndexedDatabase::default();
    let mut features = Vec::new();
    for index in 0..200u32 {
        db.peptides.push(Peptide {
            sequence: format!("PEPTIDE{}K", "G".repeat(index as usize + 1))
                .as_bytes()
                .into(),
            proteins: vec![Arc::from(format!("P{index}"))].into_iter().collect(),
            ..Peptide::default()
        });
        features.push(Feature {
            psm_id: index as usize,
            peptide_idx: PeptideIx(index),
            label: 1,
            discriminant_score: 10.0 + index as f32,
            ..Feature::default()
        });
    }
    db.peptides.push(Peptide {
        sequence: b"VARIANTHK".as_slice().into(),
        proteins: vec![Arc::from("P0_V")].into_iter().collect(),
        decoy: true,
        ..Peptide::default()
    });
    features.push(Feature {
        psm_id: 200,
        peptide_idx: PeptideIx(200),
        label: -1,
        discriminant_score: 1000.0,
        ..Feature::default()
    });
    let mut fasta = Fasta::parse(">P0\nPEPTIDEK\n".into(), "rev_", true).unwrap();
    let variants = Fasta::parse(">P0_V\nVARIANTHK\n".into(), "rev_", true).unwrap();
    fasta.extend_with_class(variants, ProteinClass::Variant);
    let classes = ProteinClasses::of_fasta(&fasta).unwrap();

    let mut counted = features.clone();
    crate::fdr::picked_protein(&db, &mut counted);
    assert_eq!(counted[199].protein_q, 0.01);

    infer_proteins_without_variants(&db, &classes, &mut features, |features| {
        crate::fdr::picked_protein(&db, features)
    });
    assert_eq!(features[199].protein_q, 0.005);
    assert_eq!(features[200].protein_q, 1.0);
}

#[test]
fn a_reference_peptide_keeps_only_its_reference_proteins() {
    let mut fasta = Fasta::parse(">P1\nPEPTIDEK\n".into(), "rev_", true).unwrap();
    let variants = Fasta::parse(">P1_V\nPEPTIDEH\n".into(), "rev_", true).unwrap();
    fasta.extend_with_class(variants, ProteinClass::Variant);
    let entrapment = Fasta::parse(">AT1\nENTRAPK\n".into(), "rev_", true).unwrap();
    fasta.extend_with_class(entrapment, ProteinClass::Entrapment);
    let classes = ProteinClasses::of_fasta(&fasta).unwrap();

    let occurrence = |protein: &str| crate::enzyme::ProteinOccurrence {
        protein: Arc::from(protein),
        start: Some(0),
        prev_aa: None,
        next_aa: None,
        source: None,
        met_clipped: false,
    };
    let peptide = |proteins: &[&str], decoy: bool| Peptide {
        sequence: b"PEPTIDEK".as_slice().into(),
        proteins: proteins.iter().map(|&p| Arc::from(p)).collect(),
        protein_sites: proteins
            .iter()
            .map(|p| occurrence(p))
            .collect::<Vec<_>>()
            .into(),
        decoy,
        ..Peptide::default()
    };
    let mut peptides = vec![
        peptide(&["P1", "P1_V"], false),
        peptide(&["P1", "P1_V"], true),
        peptide(&["P1_V"], false),
        peptide(&["AT1", "P1_V"], false),
        peptide(&["P1"], false),
    ];
    let shared = keep_reference_proteins(&mut peptides, &classes);
    // Only the target that a reference protein explains is counted; its decoy
    // is mapped the same way.
    assert_eq!(shared, 1);
    let proteins = |peptide: &Peptide| {
        peptide
            .proteins
            .iter()
            .map(|p| p.to_string())
            .collect::<Vec<_>>()
    };
    assert_eq!(proteins(&peptides[0]), vec!["P1"]);
    assert_eq!(peptides[0].protein_sites.len(), 1);
    assert_eq!(peptides[0].protein_sites[0].protein.as_ref(), "P1");
    assert_eq!(proteins(&peptides[1]), vec!["P1"]);
    // A variant peptide stays a variant peptide.
    assert_eq!(proteins(&peptides[2]), vec!["P1_V"]);
    assert!(classes.is_variant(peptides[2].proteins.iter()));
    // Shared with an entrapment protein only: unchanged, so it stays an
    // original target and is not turned into an entrapment discovery.
    assert_eq!(proteins(&peptides[3]), vec!["AT1", "P1_V"]);
    assert_eq!(proteins(&peptides[4]), vec!["P1"]);
}

#[test]
fn without_variant_proteins_no_peptide_is_touched() {
    let fasta = Fasta::parse(">P1\nPEPTIDEK\n".into(), "rev_", true).unwrap();
    let mut entrapment_only = fasta.clone();
    entrapment_only.extend_with_class(
        Fasta::parse(">AT1\nENTRAPK\n".into(), "rev_", true).unwrap(),
        ProteinClass::Entrapment,
    );
    let classes = ProteinClasses::of_fasta(&entrapment_only).unwrap();
    let mut peptides = vec![Peptide {
        proteins: vec![Arc::from("P1"), Arc::from("AT1")]
            .into_iter()
            .collect(),
        ..Peptide::default()
    }];
    assert_eq!(keep_reference_proteins(&mut peptides, &classes), 0);
    assert_eq!(peptides[0].proteins.len(), 2);
}
