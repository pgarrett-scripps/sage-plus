use super::*;

// --- Ported reference doctests: construct_ambiguity_intervals ----------

#[test]
fn construct_forward() {
    assert_eq!(
        construct_ambiguity_intervals(&[0, 1, 1, 1, 0, 0, 0], false),
        vec![(0, 1), (4, 6)]
    );
    assert_eq!(
        construct_ambiguity_intervals(&[0, 1, 1, 1, 0, 0, 1], false),
        vec![(0, 1), (4, 6)]
    );
}

#[test]
fn construct_reverse() {
    assert_eq!(
        construct_ambiguity_intervals(&[0, 0, 1, 1, 1, 0, 0], true),
        vec![(0, 1), (4, 6)]
    );
}

// --- Ported reference doctests: combine_ambiguity_intervals ------------

#[test]
fn combine() {
    assert_eq!(
        combine_ambiguity_intervals(&[vec![(0, 1), (4, 6)], vec![(0, 1)]]),
        vec![(0, 1)]
    );
    assert_eq!(
        combine_ambiguity_intervals(&[vec![(0, 1), (4, 6)], vec![(0, 1), (4, 5)]]),
        vec![(0, 1), (4, 5)]
    );
    assert_eq!(
        combine_ambiguity_intervals(&[vec![(0, 1), (4, 6)], vec![(0, 4), (5, 6)]]),
        vec![(0, 1), (5, 6)]
    );
    assert_eq!(
        combine_ambiguity_intervals(&[vec![(2, 5)], vec![(3, 6)]]),
        vec![(3, 5)]
    );
    assert_eq!(
        combine_ambiguity_intervals(&[vec![(0, 1)], vec![(4, 6)]]),
        Vec::<(usize, usize)>::new()
    );
}

// --- Ported reference doctests: mass_shift_interval --------------------

#[test]
fn mass_shift() {
    assert_eq!(
        mass_shift_interval(&[1, 1, 1, 0, 0, 0, 0], &[0, 0, 0, 0, 1, 1, 1]),
        Some((3, 3))
    );
    assert_eq!(
        mass_shift_interval(&[1, 1, 1, 0, 0, 0, 0], &[0, 0, 0, 1, 1, 1, 1]),
        Some((3, 3))
    );
    assert_eq!(
        mass_shift_interval(&[1, 1, 0, 0, 0, 0, 0], &[0, 0, 0, 0, 1, 1, 1]),
        Some((2, 3))
    );
    assert_eq!(
        mass_shift_interval(&[0, 0, 0, 0, 0, 0, 0], &[0, 0, 0, 0, 1, 1, 1]),
        Some((0, 3))
    );
    assert_eq!(
        mass_shift_interval(&[1, 1, 1, 0, 0, 0, 0], &[0, 0, 0, 0, 0, 0, 0]),
        Some((3, 6))
    );
    assert_eq!(
        mass_shift_interval(&[1, 1, 1, 1, 1, 0, 0], &[0, 0, 0, 0, 1, 1, 1]),
        None
    );
}

// --- End-to-end rendering ---------------------------------------------

fn peptide(seq: &str) -> Peptide {
    Peptide {
        sequence: seq.as_bytes().into(),
        modifications: crate::peptide::CompactModifications::default(),
        ..Default::default()
    }
}

#[test]
fn fully_covered_has_no_intervals() {
    let p = peptide("PEPTIDE");
    let cov = vec![1u16; 7];
    let a = annotate(&p, &cov, &cov, None);
    assert_eq!(a.sequence, "PEPTIDE");
    assert_eq!(a.mass_shift, 0.0);
}

#[test]
fn equal_mass_labels_keep_their_site_specific_names() {
    let builder: crate::database::Builder = serde_json::from_value(serde_json::json!({
        "generate_decoys": false,
        "static_mods": {
            "K": {
                "mass": 0.0,
                "name": "Lys6",
                "channel_offsets": {"light": 0.0, "heavy": 6.020129}
            },
            "R": {
                "mass": 0.0,
                "name": "Arg6",
                "channel_offsets": {"light": 0.0, "heavy": 6.020129}
            }
        }
    }))
    .unwrap();
    let parameters = builder.make_parameters();
    let peptides = parameters.peptides_from_tsv("sequence\nPEPTIDEK\nPEPTIDER\n");
    let database = parameters.build_from_peptides(peptides);

    for expected in ["PEPTIDEK[Lys6]", "PEPTIDER[Arg6]"] {
        let peptide = database
            .peptides
            .iter()
            .find(|peptide| peptide.to_string() == expected)
            .unwrap();
        let coverage = vec![1; peptide.sequence.len()];
        assert_eq!(
            annotate(peptide, &coverage, &coverage, None).sequence,
            expected
        );
    }
}

#[test]
fn ambiguous_nterm_is_wrapped() {
    // The first two residues have neither forward (b1/b2) nor reverse
    // (y6/y5) evidence, so their order is ambiguous and they are wrapped.
    // (reverse[0] is always 0 since no y-ion maps to the N-terminal residue.)
    let p = peptide("PEPTIDE");
    let forward = vec![0, 0, 1, 1, 1, 1, 1];
    let reverse = vec![0, 0, 1, 1, 1, 1, 1];
    let a = annotate(&p, &forward, &reverse, None);
    assert_eq!(a.sequence, "(?PE)PTIDE");
}

#[test]
fn localized_mass_shift_is_bracketed() {
    // Forward stops after residue 2, reverse starts at residue 4 -> the
    // shift localizes to the single gap residue (index 3).
    let p = peptide("PEPTIDE");
    let forward = vec![1, 1, 1, 0, 0, 0, 0];
    let reverse = vec![0, 0, 0, 0, 1, 1, 1];
    let a = annotate(&p, &forward, &reverse, Some(79.96633));
    assert!(a.sequence.contains("T[+79.96633]"), "got: {}", a.sequence);
    assert_eq!(a.mass_shift, 79.96633);
}

#[test]
fn labile_mass_shift_is_prefixed() {
    // Forward and reverse coverage overlap -> shift cannot be localized.
    let p = peptide("PEPTIDE");
    let forward = vec![1, 1, 1, 1, 1, 0, 0];
    let reverse = vec![0, 0, 0, 0, 1, 1, 1];
    let a = annotate(&p, &forward, &reverse, Some(100.0));
    assert!(a.sequence.starts_with("{+100}"), "got: {}", a.sequence);
}
