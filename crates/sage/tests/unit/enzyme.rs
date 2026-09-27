use quickcheck_macros::quickcheck;
use std::collections::HashSet;

use super::*;

#[test]
fn hash_digest() {
    let mut digests = vec![
        Digest {
            decoy: false,
            semi_enzymatic: false,
            sequence: "MADEEK".into(),
            missed_cleavages: 0,
            position: Position::Nterm,
            protein: Arc::from(String::default()),
            protein_start: Some(0),
            prev_aa: None,
            next_aa: None,
            expanded_from: None,
        },
        Digest {
            decoy: false,
            semi_enzymatic: false,
            sequence: "MADEEK".into(),
            missed_cleavages: 0,
            position: Position::Nterm,
            protein: Arc::from(String::default()),
            protein_start: Some(0),
            prev_aa: None,
            next_aa: None,
            expanded_from: None,
        },
    ];

    // Make sure hashing a digest works
    let set = digests.drain(..).collect::<HashSet<_>>();
    assert_eq!(set.len(), 1);

    let mut digests = vec![
        Digest {
            decoy: false,
            semi_enzymatic: false,
            sequence: "MADEEK".into(),
            missed_cleavages: 0,
            position: Position::Nterm,
            protein: Arc::from(String::default()),
            protein_start: Some(0),
            prev_aa: None,
            next_aa: None,
            expanded_from: None,
        },
        Digest {
            decoy: false,
            semi_enzymatic: false,
            sequence: "MADEEK".into(),
            missed_cleavages: 0,
            position: Position::Internal,
            protein: Arc::from(String::default()),
            protein_start: Some(0),
            prev_aa: None,
            next_aa: None,
            expanded_from: None,
        },
    ];

    // // Make sure hashing a digest works
    let set = digests.drain(..).collect::<HashSet<_>>();
    assert_eq!(set.len(), 2);
}

#[test]
fn trypsin() {
    let sequence = "MADEEKLPPGWEKRMSRSSGRVYYFNHITNASQWERPSGN";
    let expected = vec![
        ("MADEEK".into(), Position::Nterm),
        ("LPPGWEK".into(), Position::Internal),
        ("MSR".into(), Position::Internal),
        ("SSGR".into(), Position::Internal),
        ("VYYFNHITNASQWERPSGN".into(), Position::Cterm),
    ];

    let tryp = EnzymeParameters {
        clip_n_term_met: false,
        min_len: 2,
        max_len: 50,
        missed_cleavages: 0,
        enzyme: Enzyme::new("KR", "P", true, false),
        ambiguous_variants: None,
    };

    assert_eq!(
        expected,
        tryp.digest(sequence, Arc::default())
            .into_iter()
            .map(|d| (d.sequence, d.position))
            .collect::<Vec<_>>()
    );
}

#[test]
fn trypsin_missed_cleavage() {
    let sequence = "MADEEKLPPGWEKRMSRSSGRVYYFNHITNASQWERPSGN";
    let expected = vec![
        "MADEEK",
        "LPPGWEK",
        "R",
        "MSR",
        "SSGR",
        "VYYFNHITNASQWERPSGN",
        "MADEEKLPPGWEK",
        "LPPGWEKR",
        "RMSR",
        "MSRSSGR",
        "SSGRVYYFNHITNASQWERPSGN",
    ];

    let tryp = EnzymeParameters {
        clip_n_term_met: false,
        min_len: 0,
        max_len: 50,
        missed_cleavages: 1,
        enzyme: Enzyme::new("KR", "P", true, false),
        ambiguous_variants: None,
    };

    assert_eq!(
        expected,
        tryp.digest(sequence, Arc::default())
            .into_iter()
            .map(|d| d.sequence)
            .collect::<Vec<_>>()
    );
}

#[test]
fn trypsin_missed_cleavage_2() {
    let sequence = "MADEEKLPPGWEKRMSRSSGRVYYFNHITNASQWERPSGN";
    let expected = vec![
        "MADEEK",
        "LPPGWEK",
        "R",
        "MSR",
        "SSGR",
        "VYYFNHITNASQWERPSGN",
        "MADEEKLPPGWEK",
        "LPPGWEKR",
        "RMSR",
        "MSRSSGR",
        "SSGRVYYFNHITNASQWERPSGN",
        "MADEEKLPPGWEKR",
        "LPPGWEKRMSR",
        "RMSRSSGR",
        "MSRSSGRVYYFNHITNASQWERPSGN",
    ];

    let tryp = EnzymeParameters {
        clip_n_term_met: false,
        min_len: 0,
        max_len: 50,
        missed_cleavages: 2,
        enzyme: Enzyme::new("KR", "P", true, false),
        ambiguous_variants: None,
    };

    assert_eq!(
        expected,
        tryp.digest(sequence, Arc::default())
            .into_iter()
            .map(|d| d.sequence)
            .collect::<Vec<_>>()
    );
}

#[test]
fn test_trypsin_pro() {
    let sequence = "MADEEKLPPGWEKRMSRSSGRVYYFNHITNASQWERPSGN";
    let expected = vec![
        "MADEEK",
        "LPPGWEK",
        "MSR",
        "SSGR",
        "VYYFNHITNASQWER",
        "PSGN",
    ];

    let tryp = EnzymeParameters {
        clip_n_term_met: false,
        min_len: 2,
        max_len: 50,
        missed_cleavages: 0,
        enzyme: Enzyme::new("KR", "", true, false),
        ambiguous_variants: None,
    };

    assert_eq!(
        expected,
        tryp.digest(sequence, Arc::default())
            .into_iter()
            .map(|d| d.sequence)
            .collect::<Vec<_>>()
    );
}

#[test]
fn test_asp_n() {
    let sequence = "MADEEKLPPGWEKRMSRSSGRVYYFNHITNASQWERPSGNW";
    let expected = vec!["MA", "DEEKLPPGWEKRMSRSSGRVYYFNHITNASQWERPSGNW"];

    let tryp = EnzymeParameters {
        clip_n_term_met: false,
        min_len: 1,
        max_len: 50,
        missed_cleavages: 0,
        enzyme: Enzyme::new("D", "", false, false),
        ambiguous_variants: None,
    };

    assert_eq!(
        expected,
        tryp.digest(sequence, Arc::default())
            .into_iter()
            .map(|d| d.sequence)
            .collect::<Vec<_>>()
    );
}

#[test]
fn test_chymotrypsin_pro() {
    let sequence = "MADEEKLPPGWEKRMSRSSGRVYYFNHITNASQWERPSGNW";
    let expected = vec![
        "MADEEKL",
        "PPGW",
        "EKRMSRSSGRVY",
        "Y",
        "F",
        "NHITNASQW",
        "ERPSGNW",
    ];

    let tryp = EnzymeParameters {
        clip_n_term_met: false,
        min_len: 1,
        max_len: 50,
        missed_cleavages: 0,
        enzyme: Enzyme::new("FYWL", "", true, false),
        ambiguous_variants: None,
    };

    assert_eq!(
        expected,
        tryp.digest(sequence, Arc::default())
            .into_iter()
            .map(|d| d.sequence)
            .collect::<Vec<_>>()
    );
}

#[test]
fn nonspecific_digest_5() {
    let sequence = "MADEEKLPPGWEKRMSRSSGRVYYFNHITNASQWERPSGNW";

    let expected = sequence
        .as_bytes()
        .windows(5)
        .flat_map(std::str::from_utf8)
        .collect::<Vec<_>>();

    let tryp = EnzymeParameters {
        clip_n_term_met: false,
        min_len: 5,
        max_len: 5,
        missed_cleavages: 0,
        enzyme: None,
        ambiguous_variants: None,
    };

    assert_eq!(
        expected,
        tryp.digest(sequence, Arc::default())
            .into_iter()
            .map(|d| d.sequence)
            .collect::<Vec<_>>()
    );
}

#[test]
fn nonspecific_digest_5_7() {
    let sequence = "MADEEKLPPGWEKRMSRSSGRVYYFNHITNASQWERPSGNW";

    let expected = (5..=7)
        .flat_map(|window| {
            sequence
                .as_bytes()
                .windows(window)
                .flat_map(std::str::from_utf8)
        })
        .collect::<Vec<_>>();

    let tryp = EnzymeParameters {
        clip_n_term_met: false,
        min_len: 5,
        max_len: 7,
        missed_cleavages: 0,
        enzyme: Enzyme::new("", "", true, false),
        ambiguous_variants: None,
    };

    assert_eq!(
        expected,
        tryp.digest(sequence, Arc::default())
            .into_iter()
            .map(|d| d.sequence)
            .collect::<Vec<_>>()
    );
}

#[test]
fn no_digest() {
    let sequence = "MADEEKLPPGWEKRMSRSSGRVYYFNHITNASQWERPSGNW";
    let expected = vec![sequence];

    let tryp = EnzymeParameters {
        clip_n_term_met: false,
        min_len: 0,
        max_len: usize::MAX,
        missed_cleavages: 0,
        enzyme: Enzyme::new("$", "", true, false),
        ambiguous_variants: None,
    };

    assert_eq!(
        expected,
        tryp.digest(sequence, Arc::default())
            .into_iter()
            .map(|d| d.sequence)
            .collect::<Vec<_>>()
    );
}

#[test]
fn preserve_repeated_sequence_coordinates() {
    let sequence = "KVEGAQNQGKKVEGAQNQGK";
    let expected = vec![
        ("VEGAQNQGK".to_string(), Some(1)),
        ("VEGAQNQGK".to_string(), Some(11)),
    ];

    let tryp = EnzymeParameters {
        clip_n_term_met: false,
        min_len: 2,
        max_len: usize::MAX,
        missed_cleavages: 0,
        enzyme: Enzyme::new("KR", "", true, false),
        ambiguous_variants: None,
    };

    assert_eq!(
        expected,
        tryp.digest(sequence, Arc::default())
            .into_iter()
            .map(|d| (d.sequence.to_string(), d.protein_start))
            .collect::<Vec<_>>()
    );
}

#[test]
fn mini_semi_trypsin() {
    let sequence = "MADEEK";
    let expected = vec![
        "MADEEK", "ADEEK", "MA", "DEEK", "MAD", "EEK", "MADE", "EK", "MADEE",
    ];

    let tryp = EnzymeParameters {
        clip_n_term_met: false,
        min_len: 2,
        max_len: 50,
        missed_cleavages: 0,
        enzyme: Enzyme::new("KR", "P", true, true),
        ambiguous_variants: None,
    };

    assert_eq!(
        expected,
        tryp.digest(sequence, Arc::default())
            .into_iter()
            .map(|d| d.sequence)
            .collect::<Vec<_>>()
    );
}

#[test]
fn semi_trypsin_trypsin_missed_cleavage() {
    let sequence = "MADEEKLPPGWEK";
    let expected = vec![
        "MADEEK",
        "LPPGWEK",       // normal KR
        "MADEEKLPPGWEK", // one missed cleavage
        "ADEEK",
        "DEEK",
        "MAD",
        "EEK",
        "MADE",
        "MADEE", // normal half-tryptics
        "PPGWEK",
        "PGWEK",
        "LPP",
        "GWEK",
        "LPPG",
        "WEK",
        "LPPGW",
        "LPPGWE",
        "ADEEKLPPGWEK",
        "DEEKLPPGWEK",
        "EEKLPPGWEK", // one missed cleavage half-tryptics
        "EKLPPGWEK",
        "KLPPGWEK",
        "MADEEKL",
        "MADEEKLP",
        "MADEEKLPP",
        "MADEEKLPPG",
        "MADEEKLPPGW",
        "MADEEKLPPGWE",
    ];

    let tryp = EnzymeParameters {
        clip_n_term_met: false,
        min_len: 3,
        max_len: 50,
        missed_cleavages: 1,
        enzyme: Enzyme::new("KR", "P", true, true),
        ambiguous_variants: None,
    };

    for (digest, expected) in tryp
        .digest(sequence, Arc::default())
        .into_iter()
        .zip(expected)
    {
        assert_eq!(digest.sequence, expected);
        // reverse and skip the first (C-terminal) AA, counting interior missed cleavages
        let missed_cleavages = digest
            .sequence
            .as_bytes()
            .iter()
            .rev()
            .skip(1)
            .map(|s| (*s == b'K' || *s == b'R') as u8)
            .sum::<u8>();
        assert_eq!(
            missed_cleavages, digest.missed_cleavages,
            "{}",
            digest.sequence
        );

        if digest.sequence.starts_with("MAD") && digest.sequence != sequence {
            assert_eq!(digest.position, Position::Nterm);
        }
    }
}

#[test]
fn custom_cleavages_add_both_sides_with_missed_cleavages() {
    let sequence = "AAKAPEPTIDERQQQK";
    let tryp = EnzymeParameters {
        clip_n_term_met: false,
        min_len: 3,
        max_len: 50,
        missed_cleavages: 1,
        enzyme: Enzyme::new("KR", "P", true, false),
        ambiguous_variants: None,
    };

    let digests = tryp.digest_with_custom_cleavages(sequence, Arc::default(), &[8]);
    let by_sequence = digests
        .iter()
        .map(|digest| (digest.sequence.as_str(), digest))
        .collect::<std::collections::HashMap<_, _>>();

    assert_eq!(by_sequence["APEPT"].missed_cleavages, 0);
    assert_eq!(by_sequence["AAKAPEPT"].missed_cleavages, 1);
    assert_eq!(by_sequence["IDER"].missed_cleavages, 0);
    assert_eq!(by_sequence["IDERQQQK"].missed_cleavages, 1);
    assert!(by_sequence["APEPT"].semi_enzymatic);
    assert!(by_sequence["IDER"].semi_enzymatic);
}

#[test]
fn existing_enzyme_boundary_does_not_add_duplicates() {
    let sequence = "AAKAPEPTIDER";
    let tryp = EnzymeParameters {
        clip_n_term_met: false,
        min_len: 3,
        max_len: 50,
        missed_cleavages: 1,
        enzyme: Enzyme::new("KR", "P", true, false),
        ambiguous_variants: None,
    };

    let ordinary = tryp.digest(sequence, Arc::default());
    let custom = tryp.digest_with_custom_cleavages(sequence, Arc::default(), &[3]);
    assert_eq!(ordinary, custom);
}

#[test]
fn custom_cleavages_are_additive_to_no_digest_and_redundant_for_nonspecific() {
    let sequence = "ACDEFGHIK";
    let no_digest = EnzymeParameters {
        clip_n_term_met: false,
        min_len: 3,
        max_len: 50,
        missed_cleavages: 0,
        enzyme: Enzyme::new("$", "", true, false),
        ambiguous_variants: None,
    };
    let sequences = no_digest
        .digest_with_custom_cleavages(sequence, Arc::default(), &[4])
        .into_iter()
        .map(|digest| digest.sequence)
        .collect::<HashSet<_>>();
    assert!(sequences.contains(&b"ACDE"[..]));
    assert!(sequences.contains(&b"FGHIK"[..]));
    assert!(sequences.contains(sequence.as_bytes()));

    let nonspecific = EnzymeParameters {
        clip_n_term_met: false,
        min_len: 3,
        max_len: 5,
        missed_cleavages: 0,
        enzyme: None,
        ambiguous_variants: None,
    };
    assert_eq!(
        nonspecific.digest(sequence, Arc::default()),
        nonspecific.digest_with_custom_cleavages(sequence, Arc::default(), &[4])
    );
}

#[test]
fn nonspecific_digest_spans_share_one_protein_allocation() {
    let sequence: ProteinSequence = "ACDEFGHIK".into();
    let nonspecific = EnzymeParameters {
        clip_n_term_met: false,
        min_len: 3,
        max_len: 3,
        missed_cleavages: 0,
        enzyme: None,
        ambiguous_variants: None,
    };
    let digests =
        nonspecific.digest_protein_with_custom_cleavages(&sequence, Arc::from("protein"), &[]);

    assert_eq!(digests.len(), 7);
    assert!(digests
        .windows(2)
        .all(|pair| pair[0].sequence.shares_storage_with(&pair[1].sequence)));
    assert!(digests
        .iter()
        .all(|digest| digest.sequence.storage_len() == sequence.len()));
}

#[test]
fn grouping_uses_sequence_content_instead_of_storage_identity() {
    let no_digest = EnzymeParameters {
        clip_n_term_met: false,
        min_len: 3,
        max_len: 50,
        missed_cleavages: 0,
        enzyme: Enzyme::new("$", "", true, false),
        ambiguous_variants: None,
    };
    let first: ProteinSequence = "PEPTIDE".into();
    let second: ProteinSequence = "PEPTIDE".into();
    let mut digests =
        no_digest.digest_protein_with_custom_cleavages(&first, Arc::from("first"), &[]);
    digests.extend(no_digest.digest_protein_with_custom_cleavages(
        &second,
        Arc::from("second"),
        &[],
    ));

    let groups = group_digests(digests);
    assert_eq!(groups.len(), 1);
    assert_eq!(groups[0].origins.len(), 2);
}

/// Helper struct for generation of random sequences of valid amino acids
#[derive(Clone, Debug)]
struct RandomSequence {
    sequence: String,
}

impl quickcheck::Arbitrary for RandomSequence {
    fn arbitrary(g: &mut quickcheck::Gen) -> Self {
        let bytes = (0..g.size())
            .filter_map(|_| g.choose(&VALID_AA))
            .copied()
            .collect();
        Self {
            sequence: String::from_utf8(bytes).unwrap(),
        }
    }
}

#[quickcheck]
/// Check that our strict ordering of missed cleavage generation is not
/// broken for arbitrary peptide sequences
fn quickcheck_semi_missed_cleavages(RandomSequence { sequence }: RandomSequence) {
    let tryp = EnzymeParameters {
        clip_n_term_met: false,
        min_len: 3,
        max_len: 50,
        missed_cleavages: 2,
        enzyme: Enzyme::new("KR", "", true, true),
        ambiguous_variants: None,
    };

    for digest in tryp.digest(&sequence, Arc::default()) {
        // reverse and skip the first (C-terminal) AA, counting interior missed cleavages
        let missed_cleavages = digest
            .sequence
            .as_bytes()
            .iter()
            .rev()
            .skip(1)
            .map(|s| (*s == b'K' || *s == b'R') as u8)
            .sum::<u8>();
        assert_eq!(
            missed_cleavages, digest.missed_cleavages,
            "{}",
            digest.sequence
        );

        assert!(digest.missed_cleavages <= 2);
    }
}

#[test]
fn unsupported_enzyme_residues_are_errors() {
    for (cleave, restrict, field) in [
        ("KB", "", "cleave_at"),
        ("KZ", "P", "cleave_at"),
        ("kr", "", "cleave_at"),
        ("KR", "X", "restrict"),
        ("KR", "J", "restrict"),
    ] {
        let error = Enzyme::try_new(cleave, restrict, true, false)
            .err()
            .unwrap_or_else(|| panic!("accepted {cleave}/{restrict}"));
        assert!(error.contains(field), "{error}");
    }
    assert!(Enzyme::try_new("KR", "P", true, false).unwrap().is_some());
    assert!(Enzyme::try_new("$", "", true, false).unwrap().is_some());
    assert!(Enzyme::try_new("", "", true, false).unwrap().is_none());
}

fn clipping_trypsin(clip_n_term_met: bool, semi_enzymatic: bool) -> EnzymeParameters {
    EnzymeParameters {
        clip_n_term_met,
        ambiguous_variants: None,
        min_len: 3,
        max_len: 50,
        missed_cleavages: 1,
        enzyme: Enzyme::new("KR", "P", true, semi_enzymatic),
    }
}

#[test]
fn metap_clips_only_before_small_residues() {
    for second in b"GASTCPV" {
        assert!(metap_clips(&[b'M', *second, b'K']));
    }
    for sequence in [&b"MKDER"[..], b"MLDER", b"MEDER", b"ASDER", b"M", b""] {
        assert!(!metap_clips(sequence), "{sequence:?}");
    }
}

#[test]
fn metap_clipping_adds_only_protein_n_terminal_peptides_from_residue_two() {
    let sequence = "MSDEREVAEAKLPPGWEKR";
    let unclipped = clipping_trypsin(false, false).digest(sequence, Arc::from("P1"));
    let clipped = clipping_trypsin(true, false).digest(sequence, Arc::from("P1"));

    // The unclipped digests stay, and only peptides from residue 2 are added.
    assert_eq!(clipped[..unclipped.len()], unclipped[..]);
    let added = &clipped[unclipped.len()..];
    let summary = added
        .iter()
        .map(|digest| {
            (
                String::from_utf8(digest.sequence.to_vec()).unwrap(),
                digest.missed_cleavages,
            )
        })
        .collect::<Vec<_>>();
    assert_eq!(
        summary,
        vec![("SDER".to_string(), 0), ("SDEREVAEAK".to_string(), 1)]
    );
    for digest in added {
        assert_eq!(digest.position, Position::Nterm);
        assert!(!digest.semi_enzymatic);
        assert!(!digest.decoy);
        // Coordinates stay relative to the real protein sequence.
        assert_eq!(digest.protein_start, Some(1));
        assert_eq!(digest.prev_aa, Some(b'M'));
    }

    for sequence in ["MKDEREVAEAK", "MLDEREVAEAK", "ASDEREVAEAK"] {
        assert_eq!(
            clipping_trypsin(true, false).digest(sequence, Arc::default()),
            clipping_trypsin(false, false).digest(sequence, Arc::default()),
            "{sequence}"
        );
    }
}

#[test]
fn metap_clipping_relabels_semi_enzymatic_spans_after_the_methionine() {
    let digests = clipping_trypsin(true, true).digest("MSDEREVAEAK", Arc::default());
    let from_residue_two = |peptide: &[u8]| {
        digests
            .iter()
            .filter(|digest| digest.protein_start == Some(1) && &digest.sequence[..] == peptide)
            .collect::<Vec<_>>()
    };
    // Once, as the clipped protein N-terminal peptide, not as a semi-enzymatic
    // internal one.
    let full = from_residue_two(b"SDER");
    assert_eq!(full.len(), 1);
    assert_eq!(full[0].position, Position::Nterm);
    assert!(!full[0].semi_enzymatic);
    // A C-terminally non-enzymatic end stays semi-enzymatic.
    let semi = from_residue_two(b"SDE");
    assert_eq!(semi.len(), 1);
    assert_eq!(semi[0].position, Position::Nterm);
    assert!(semi[0].semi_enzymatic);
}

#[test]
fn metap_clipping_ignores_nonspecific_and_clips_whole_proteins() {
    let nonspecific = |clip_n_term_met| EnzymeParameters {
        clip_n_term_met,
        ambiguous_variants: None,
        min_len: 3,
        max_len: 5,
        missed_cleavages: 0,
        enzyme: None,
    };
    assert_eq!(
        nonspecific(true).digest("MSPEPTIDE", Arc::default()),
        nonspecific(false).digest("MSPEPTIDE", Arc::default())
    );

    let no_digest = EnzymeParameters {
        clip_n_term_met: true,
        ambiguous_variants: None,
        min_len: 3,
        max_len: 50,
        missed_cleavages: 0,
        enzyme: Enzyme::new("$", "", true, false),
    };
    let digests = no_digest
        .digest("MSPEPTIDE", Arc::default())
        .into_iter()
        .map(|digest| (digest.sequence, digest.position))
        .collect::<Vec<_>>();
    assert_eq!(
        digests,
        vec![
            ("MSPEPTIDE".into(), Position::Full),
            ("SPEPTIDE".into(), Position::Full)
        ]
    );
}

#[test]
fn metap_clipping_does_not_count_a_cut_after_the_methionine_as_missed() {
    let after_met = EnzymeParameters {
        clip_n_term_met: true,
        ambiguous_variants: None,
        min_len: 3,
        max_len: 50,
        missed_cleavages: 1,
        enzyme: Enzyme::new("M", "", true, false),
    };
    let digests = after_met.digest("MSDEMKLL", Arc::default());
    let clipped = digests
        .iter()
        .filter(|digest| &digest.sequence[..] == b"SDEM")
        .collect::<Vec<_>>();
    assert_eq!(clipped.len(), 1);
    assert_eq!(clipped[0].position, Position::Nterm);
    assert_eq!(clipped[0].missed_cleavages, 0);
}

#[test]
fn group_reference_is_the_most_enzymatic_occurrence() {
    // "a" sorts before "b". In "a", ASEQK starts the protein but K-P is not a
    // trypsin site, so it is semi-enzymatic; in "b" it is the fully enzymatic
    // Met-clipped N-terminal peptide. Both are protein N-terminal.
    let enzyme = clipping_trypsin(true, true);
    let mut digests = enzyme.digest("ASEQKPLLR", Arc::from("a"));
    digests.extend(enzyme.digest("MASEQKGLLR", Arc::from("b")));
    let groups = group_digests(digests);
    let group = groups
        .iter()
        .find(|group| {
            &group.reference.sequence[..] == b"ASEQK" && group.reference.position == Position::Nterm
        })
        .unwrap();
    assert_eq!(group.origins.len(), 2);
    assert!(!group.reference.semi_enzymatic);
    assert_eq!(&*group.reference.protein, "b");
}
