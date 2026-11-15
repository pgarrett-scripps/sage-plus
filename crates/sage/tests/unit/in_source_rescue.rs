use super::*;
use crate::database::{Builder, IndexedDatabase};
use crate::fasta::Fasta;
use crate::in_source::{PeakApexes, PeakShape};

/// A real database: static carbamidomethyl C, variable oxidised M and protein
/// N-terminal acetylation, so the candidates are built from peptides exactly as
/// the search expands them.
fn database() -> IndexedDatabase {
    let builder: Builder = serde_json::from_value(serde_json::json!({
        "enzyme": {"missed_cleavages": 0, "min_len": 5, "max_len": 30},
        "static_mods": {"C": 57.021464},
        "variable_mods": {"M": [15.994915], "[": [42.010565]},
        "generate_decoys": true
    }))
    .unwrap();
    let fasta = Fasta::parse(
        ">P1\nMSTEACMPLLGQRVDAHCMNEILTKSSWYPEAGCLLK\n".into(),
        "rev_",
        true,
    )
    .unwrap();
    builder.make_parameters().unwrap().build(fasta)
}

#[test]
fn the_whole_span_reproduces_the_parent() {
    let db = database();
    let mut checked = 0;
    for peptide in db.peptides.iter().filter(|p| !p.decoy) {
        let len = peptide.sequence.len();
        let copy = truncate(peptide, 0..len, None, None).unwrap();
        assert_eq!(copy.to_string(), peptide.to_string());
        assert!(
            (copy.monoisotopic - peptide.monoisotopic).abs() < 1e-3,
            "{peptide}: {} vs {}",
            copy.monoisotopic,
            peptide.monoisotopic
        );
        checked += 1;
    }
    assert!(checked >= 6, "only {checked} target peptides");
}

#[test]
fn pieces_keep_their_residues_modifications_and_the_right_mass() {
    let db = database();
    // An oxidised, acetylated, carbamidomethylated parent.
    let parent = db
        .peptides
        .iter()
        .find(|p| {
            !p.decoy
                && p.nterm.is_some()
                && p.applied_modifications()
                    .any(|m| (m.modification.mass - 15.994915).abs() < 1e-3)
        })
        .expect("an acetylated, oxidised target");
    let len = parent.sequence.len();
    let candidates = candidates_of(parent, PeptideIx(0), 4);
    // Every suffix and prefix down to 4 residues, and two losses.
    assert_eq!(candidates.len(), 2 * (len - 4) + 2);
    for candidate in &candidates {
        let piece = &candidate.peptide;
        let (start, end) = match candidate.kind.placement {
            Placement::Suffix => (len - piece.sequence.len(), len),
            Placement::Prefix => (0, piece.sequence.len()),
            Placement::SameSpan => (0, len),
            Placement::Internal => unreachable!("no internal pieces"),
        };
        assert_eq!(
            piece.sequence.as_bytes(),
            &parent.sequence.as_bytes()[start..end]
        );
        // Mass of the residues kept, their modifications, the termini kept,
        // and the deficit of the placement and loss.
        let mut expected = H2O;
        for (i, &residue) in parent.sequence.as_bytes()[start..end].iter().enumerate() {
            expected += monoisotopic(residue) + parent.modification_at(start + i);
        }
        if start == 0 {
            expected += parent.nterm.unwrap_or(0.0);
        }
        if end == len {
            expected += parent.cterm.unwrap_or(0.0);
        }
        expected -= candidate.kind.placement.structural_deficit();
        expected -= candidate.kind.loss.unwrap_or(0.0);
        assert!(
            (piece.monoisotopic - expected).abs() < 1e-3,
            "{piece}: {} vs {expected}",
            piece.monoisotopic
        );
        // Only a piece that keeps the N-terminus keeps the acetylation.
        assert_eq!(
            piece.nterm.is_some(),
            start == 0,
            "{piece} {:?}",
            candidate.kind
        );
        assert!(!piece.decoy);
        assert!(candidate.kind.placement == Placement::SameSpan || piece.semi_enzymatic);
        let site = &piece.protein_sites[0];
        assert_eq!(
            site.start,
            parent.protein_sites[0].start.map(|s| s + start as u32)
        );
    }
    // A decoy of a piece reverses its interior, as every database decoy does.
    let decoy = candidates[0].peptide.reverse();
    assert!(decoy.decoy);
    assert!((decoy.monoisotopic - candidates[0].peptide.monoisotopic).abs() < 1e-3);
}

#[test]
fn decoy_parents_and_short_parents_give_nothing_shorter_than_the_minimum() {
    let db = database();
    let decoy = db.peptides.iter().find(|p| p.decoy).unwrap();
    assert!(candidates_of(decoy, PeptideIx(0), 4).is_empty());
    let target = db.peptides.iter().find(|p| !p.decoy).unwrap();
    let len = target.sequence.len();
    // A minimum at the parent's own length leaves only the two losses.
    let only_losses = candidates_of(target, PeptideIx(0), len);
    assert_eq!(only_losses.len(), 2);
    assert!(only_losses
        .iter()
        .all(|c| c.kind.placement == Placement::SameSpan));
}

#[test]
fn windows_come_from_quantified_peaks_else_the_retention_time_tolerance() {
    let mut parents = FnvHashMap::default();
    parents.insert(PeptideIx(1), 10.0);
    parents.insert(PeptideIx(2), 20.0);
    let apexes = PeakApexes::from_peaks([(
        PeptideIx(1),
        0usize,
        PeakShape {
            apex_rt: 10.1,
            start_rt: 9.9,
            end_rt: 10.4,
            fwhm: Some(0.2),
        },
    )]);
    let windows = parent_windows(0, &parents, 0.25, Some(&apexes), None);
    assert_eq!(
        windows[&PeptideIx(1)],
        vec![Window {
            start: 9.9,
            end: 10.4
        }]
    );
    assert_eq!(
        windows[&PeptideIx(2)],
        vec![Window {
            start: 19.75,
            end: 20.25
        }]
    );
    // The null control moves each window by its own width.
    let shifted = parent_windows(0, &parents, 0.25, Some(&apexes), Some(4.0));
    let w = shifted[&PeptideIx(2)][0];
    assert!((w.start - 21.75).abs() < 1e-5 && (w.end - 22.25).abs() < 1e-5);
    assert!(!w.contains(20.0));
}

#[test]
fn q_values_are_target_decoy_competition_within_the_set() {
    // Scores 10 > 9 > 8 (decoy) > 7 > 6 (decoy) > 5.
    let evidence = [
        (7.0, false),
        (10.0, false),
        (6.0, true),
        (9.0, false),
        (8.0, true),
        (5.0, false),
    ];
    let q = competition_q_values(&evidence);
    // (decoys + 1) / targets at each threshold, then the running minimum
    // from the bottom: 10 -> 1/1, 9 -> 1/2, 8 -> 2/2, 7 -> 2/3, 6 -> 3/3, 5 -> 3/4.
    let expected = [2.0 / 3.0, 0.5, 0.75, 0.5, 2.0 / 3.0, 0.75];
    for (got, want) in q.iter().zip(expected) {
        assert!((got - want).abs() < 1e-6, "{q:?}");
    }
    assert!(competition_q_values(&[]).is_empty());
}
