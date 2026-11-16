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

fn trypsin() -> Option<crate::enzyme::Enzyme> {
    crate::enzyme::Enzyme::new("KR", "P", true, false)
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
    let candidates = candidates_of(parent, PeptideIx(0), 4, trypsin().as_ref());
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
        // Missed cleavages are the piece's own, not the parent's.
        assert_eq!(
            piece.missed_cleavages as usize,
            crate::digestion::missed_cleavages(&trypsin().unwrap(), piece.sequence.as_bytes())
        );
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
    assert!(candidates_of(decoy, PeptideIx(0), 4, trypsin().as_ref()).is_empty());
    let target = db.peptides.iter().find(|p| !p.decoy).unwrap();
    let len = target.sequence.len();
    // A minimum at the parent's own length leaves only the two losses.
    let only_losses = candidates_of(target, PeptideIx(0), len, trypsin().as_ref());
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

#[test]
fn a_b_type_piece_shifts_its_y_series_by_water_and_leaves_its_b_series() {
    use crate::ion_series::{IonSeries, Kind};
    let db = database();
    let parent = db
        .peptides
        .iter()
        .find(|p| !p.decoy && p.sequence.len() >= 8)
        .unwrap();
    let end = parent.sequence.len() - 2;
    let piece = truncate(parent, 0..end, None, Some(-H2O)).unwrap();
    let ordinary = truncate(parent, 0..end, None, None).unwrap();
    let ions = |peptide, kind| {
        IonSeries::new(peptide, kind)
            .map(|ion| ion.monoisotopic_mass)
            .collect::<Vec<_>>()
    };
    let (b_piece, b_ordinary) = (ions(&piece, Kind::B), ions(&ordinary, Kind::B));
    assert_eq!(b_piece.len(), b_ordinary.len());
    for (a, b) in b_piece.iter().zip(&b_ordinary) {
        assert!((a - b).abs() < 1e-3, "b {a} vs {b}");
    }
    let (y_piece, y_ordinary) = (ions(&piece, Kind::Y), ions(&ordinary, Kind::Y));
    assert_eq!(y_piece.len(), y_ordinary.len());
    for (a, b) in y_piece.iter().zip(&y_ordinary) {
        assert!((a - (b - H2O)).abs() < 1e-3, "y {a} vs {b} - H2O");
    }
}

#[test]
fn tied_scores_share_one_q_value() {
    // A target and a decoy tie at the top: both take (1 + 1) / 1, then the
    // running minimum from below.
    let evidence = [(9.0, false), (9.0, true), (8.0, false), (7.0, false)];
    let q = competition_q_values(&evidence);
    assert_eq!(q[0], q[1]);
    // Thresholds: 9 -> (1+1)/1 = 1; 8 -> 2/2 = 1; 7 -> 2/3.
    for (got, want) in q.iter().zip([2.0 / 3.0, 2.0 / 3.0, 2.0 / 3.0, 2.0 / 3.0]) {
        assert!((got - want).abs() < 1e-6, "{q:?}");
    }
}

#[test]
fn label_channel_masses_survive_truncation() {
    // A SILAC K8/R10 search: heavy and light partners are expanded from
    // channel offsets, and the label sits on the C-terminal K or R.
    let builder: Builder = serde_json::from_value(serde_json::json!({
        "enzyme": {"missed_cleavages": 1, "min_len": 5, "max_len": 30},
        "static_mods": {
            "CAM": {"mass": 57.021464, "sites": ["C"]},
            "SILAC-K": {"mass": 0.0, "sites": ["K"], "channel_offsets": {"light": 0.0, "heavy": 8.014199}},
            "SILAC-R": {"mass": 0.0, "sites": ["R"], "channel_offsets": {"light": 0.0, "heavy": 10.008269}}
        },
        "generate_decoys": true
    }))
    .unwrap();
    let fasta = Fasta::parse(
        ">P1\nMSTEACMPLLGQRVDAHCMNEILTKSSWYPEAGCLLK\n".into(),
        "rev_",
        true,
    )
    .unwrap();
    let db = builder.make_parameters().unwrap().build(fasta);
    let heavy = db
        .peptides
        .iter()
        .filter(|p| !p.decoy && p.label_channel_str() == Some("heavy"))
        .collect::<Vec<_>>();
    assert!(!heavy.is_empty(), "no heavy peptides were expanded");
    for parent in db.peptides.iter().filter(|p| !p.decoy) {
        let len = parent.sequence.len();
        let whole = truncate(parent, 0..len, None, None).unwrap();
        assert!(
            (whole.monoisotopic - parent.monoisotopic).abs() < 1e-3,
            "{parent}: {} vs {}",
            whole.monoisotopic,
            parent.monoisotopic
        );
        assert_eq!(whole.label_channel_str(), parent.label_channel_str());
        // The label group resolves through the channel base, as the writers
        // do for the label_group column; this panicked when truncation
        // rebuilt the modifications without their channel bases.
        assert_eq!(whole.label_group(), parent.label_group());
        // A suffix keeps the C-terminal K or R, and its label mass.
        let suffix = truncate(parent, 1..len, None, None).unwrap();
        let first = parent.sequence.as_bytes()[0];
        let expected = parent.monoisotopic
            - monoisotopic(first)
            - parent.modification_at(0)
            - parent.nterm.unwrap_or(0.0);
        let light_group = suffix.label_group();
        assert!(
            !light_group.contains("8.014") && !light_group.contains("10.008"),
            "{light_group}"
        );
        assert!(
            (suffix.monoisotopic - expected).abs() < 1e-3,
            "{parent} suffix: {} vs {expected}",
            suffix.monoisotopic
        );
    }
}

#[test]
fn the_precursor_window_must_stay_short_of_the_lightest_product() {
    // 20 ppm at 5,000 Da with isotope errors -1..3: 0.1 + 3 x 1.00335 Da.
    assert!(precursor_reach_problem(Tolerance::Ppm(-20.0, 20.0), (-1, 3), 5000.0).is_none());
    // An open search reaches the parent itself.
    let problem = precursor_reach_problem(Tolerance::Da(-500.0, 100.0), (0, 0), 5000.0).unwrap();
    assert!(problem.contains("500.00 Da"), "{problem}");
    // So do wide isotope errors on a narrow tolerance.
    assert!(precursor_reach_problem(Tolerance::Ppm(-10.0, 10.0), (-16, 0), 5000.0).is_some());
    assert!(precursor_reach_problem(Tolerance::Ppm(f32::NAN, 10.0), (0, 0), 5000.0).is_some());
}

#[test]
fn a_score_that_is_not_a_number_ends_the_competition_and_is_never_accepted() {
    // This looped forever when ties were found with `==`, which is false for
    // NaN.
    let evidence = [
        (f64::NAN, true),
        (9.0, false),
        (8.0, false),
        (f64::NAN, false),
    ];
    let q = competition_q_values(&evidence);
    assert_eq!(q[0], 1.0);
    assert_eq!(q[3], 1.0);
    // The NaN decoy is not counted: 9 -> (0+1)/1, 8 -> (0+1)/2.
    assert!((q[1] - 0.5).abs() < 1e-6, "{q:?}");
    assert!((q[2] - 0.5).abs() < 1e-6, "{q:?}");
    // Infinite scores are numbers and tie with each other.
    let q = competition_q_values(&[(f64::INFINITY, false), (f64::INFINITY, true)]);
    assert_eq!(q[0], q[1]);
}

#[test]
fn a_piece_keeps_the_flanks_and_protein_termini_of_the_ends_it_keeps() {
    let db = database();
    for parent in db
        .peptides
        .iter()
        .filter(|p| !p.decoy && p.sequence.len() >= 8)
    {
        let len = parent.sequence.len();
        let residues = parent.sequence.as_bytes();
        let prefix = truncate(parent, 0..len - 2, None, None).unwrap();
        let suffix = truncate(parent, 2..len, None, None).unwrap();
        for (piece, parent_site) in prefix.protein_sites.iter().zip(parent.protein_sites.iter()) {
            assert_eq!(piece.prev_aa, parent_site.prev_aa);
            assert_eq!(piece.next_aa, Some(residues[len - 2]));
        }
        for (piece, parent_site) in suffix.protein_sites.iter().zip(parent.protein_sites.iter()) {
            assert_eq!(piece.prev_aa, Some(residues[1]));
            assert_eq!(piece.next_aa, parent_site.next_aa);
        }
        let keeps_n = matches!(parent.position, Position::Nterm | Position::Full);
        let keeps_c = matches!(parent.position, Position::Cterm | Position::Full);
        assert_eq!(matches!(prefix.position, Position::Nterm), keeps_n);
        assert_eq!(matches!(suffix.position, Position::Cterm), keeps_c);
    }
}

#[test]
fn a_piece_without_a_labelled_residue_is_in_no_channel() {
    // Third review N5: heavy and light parents give the same unlabelled
    // piece, which the database digest would leave in no channel.
    let builder: Builder = serde_json::from_value(serde_json::json!({
        "enzyme": {"missed_cleavages": 0, "min_len": 5, "max_len": 30},
        "static_mods": {
            "SILAC-K": {"mass": 0.0, "sites": ["K"], "channel_offsets": {"light": 0.0, "heavy": 8.014199}},
            "SILAC-R": {"mass": 0.0, "sites": ["R"], "channel_offsets": {"light": 0.0, "heavy": 10.008269}}
        },
        "generate_decoys": true
    }))
    .unwrap();
    let fasta = Fasta::parse(">P1\nMSTEAGLLPVWNEILTK\n".into(), "rev_", true).unwrap();
    let db = builder.make_parameters().unwrap().build(fasta);
    let channels = db
        .peptides
        .iter()
        .filter(|p| !p.decoy && p.sequence.as_bytes() == b"MSTEAGLLPVWNEILTK")
        .collect::<Vec<_>>();
    assert_eq!(channels.len(), 2);
    assert!(channels
        .iter()
        .any(|p| p.label_channel_str() == Some("heavy")));
    assert!(channels
        .iter()
        .any(|p| p.label_channel_str() == Some("light")));
    let prefixes = channels
        .iter()
        .map(|parent| truncate(parent, 0..parent.sequence.len() - 1, None, None).unwrap())
        .collect::<Vec<_>>();
    for prefix in &prefixes {
        assert_eq!(prefix.label_channel_str(), None);
    }
    assert!(prefixes
        .windows(2)
        .all(|pair| pair[0].to_string() == pair[1].to_string()));
    assert!((prefixes[0].monoisotopic - prefixes[1].monoisotopic).abs() < 1e-3);
    // A piece that keeps the labelled residue keeps its channel.
    for parent in &channels {
        let suffix = truncate(parent, 1..parent.sequence.len(), None, None).unwrap();
        assert_eq!(suffix.label_channel_str(), parent.label_channel_str());
    }
}
