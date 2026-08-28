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
        min_len: 2,
        max_len: 50,
        missed_cleavages: 0,
        enzyme: Enzyme::new("KR", "P", true, false),
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
        min_len: 0,
        max_len: 50,
        missed_cleavages: 1,
        enzyme: Enzyme::new("KR", "P", true, false),
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
        min_len: 0,
        max_len: 50,
        missed_cleavages: 2,
        enzyme: Enzyme::new("KR", "P", true, false),
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
        min_len: 2,
        max_len: 50,
        missed_cleavages: 0,
        enzyme: Enzyme::new("KR", "", true, false),
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
        min_len: 1,
        max_len: 50,
        missed_cleavages: 0,
        enzyme: Enzyme::new("D", "", false, false),
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
        min_len: 1,
        max_len: 50,
        missed_cleavages: 0,
        enzyme: Enzyme::new("FYWL", "", true, false),
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
        min_len: 5,
        max_len: 5,
        missed_cleavages: 0,
        enzyme: None,
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
        min_len: 5,
        max_len: 7,
        missed_cleavages: 0,
        enzyme: Enzyme::new("", "", true, false),
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
        min_len: 0,
        max_len: usize::MAX,
        missed_cleavages: 0,
        enzyme: Enzyme::new("$", "", true, false),
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
        min_len: 2,
        max_len: usize::MAX,
        missed_cleavages: 0,
        enzyme: Enzyme::new("KR", "", true, false),
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
        min_len: 2,
        max_len: 50,
        missed_cleavages: 0,
        enzyme: Enzyme::new("KR", "P", true, true),
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
        min_len: 3,
        max_len: 50,
        missed_cleavages: 1,
        enzyme: Enzyme::new("KR", "P", true, true),
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
        min_len: 3,
        max_len: 50,
        missed_cleavages: 1,
        enzyme: Enzyme::new("KR", "P", true, false),
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
        min_len: 3,
        max_len: 50,
        missed_cleavages: 1,
        enzyme: Enzyme::new("KR", "P", true, false),
    };

    let ordinary = tryp.digest(sequence, Arc::default());
    let custom = tryp.digest_with_custom_cleavages(sequence, Arc::default(), &[3]);
    assert_eq!(ordinary, custom);
}

#[test]
fn custom_cleavages_are_additive_to_no_digest_and_redundant_for_nonspecific() {
    let sequence = "ACDEFGHIK";
    let no_digest = EnzymeParameters {
        min_len: 3,
        max_len: 50,
        missed_cleavages: 0,
        enzyme: Enzyme::new("$", "", true, false),
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
        min_len: 3,
        max_len: 5,
        missed_cleavages: 0,
        enzyme: None,
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
        min_len: 3,
        max_len: 3,
        missed_cleavages: 0,
        enzyme: None,
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
        min_len: 3,
        max_len: 50,
        missed_cleavages: 0,
        enzyme: Enzyme::new("$", "", true, false),
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
        min_len: 3,
        max_len: 50,
        missed_cleavages: 2,
        enzyme: Enzyme::new("KR", "", true, true),
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
