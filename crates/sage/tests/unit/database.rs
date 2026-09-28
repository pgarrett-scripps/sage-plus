use crate::enzyme::group_digests;
use std::{
    collections::{BTreeMap, HashSet},
    sync::Arc,
};

use super::*;

#[test]
fn compact_modification_validation_rejects_long_peptides() {
    let builder: Builder = serde_json::from_value(serde_json::json!({
        "enzyme": {"max_len": 256}
    }))
    .unwrap();
    let error = builder
        .make_parameters()
        .validate_compact_modifications()
        .unwrap_err();
    assert!(error.contains("must not exceed 255 residues"));
}

#[test]
fn compact_modification_validation_rejects_too_many_definitions() {
    let modifications = (0..256)
        .map(|index| serde_json::Value::from(index as f64 + 0.5))
        .collect::<Vec<_>>();
    let builder: Builder = serde_json::from_value(serde_json::json!({
        "variable_mods": {"M": modifications}
    }))
    .unwrap();
    let error = builder
        .make_parameters()
        .validate_compact_modifications()
        .unwrap_err();
    assert!(error.contains("at most 255 distinct definition and site variants"));
    assert!(error.contains("256 are required"));
}
use crate::cleavage::CustomCleavageLibrary;

#[test]
fn binary_search_slice_smoke() {
    // Make sure that our query returns the maximal set of indices
    let data = [1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0];
    let bounds = binary_search_slice(&data, |a: &f64, b| a.total_cmp(b), 1.75, 3.5);
    assert_eq!(bounds, (1, 6));
    assert!(data[bounds.0] <= 1.75);
    assert_eq!(&data[bounds.0..bounds.1], &[1.5, 2.0, 2.5, 3.0, 3.5]);

    let bounds = binary_search_slice(&data, |a: &f64, b| a.total_cmp(b), 0.0, 5.0);
    assert_eq!(bounds, (0, data.len()));
}

#[test]
fn binary_search_slice_run() {
    // Make sure that our query returns the maximal set of indices
    let data = [1.0, 1.5, 1.5, 1.5, 1.5, 2.0, 2.5, 3.0, 3.0, 3.5, 4.0];
    let (left, right) = binary_search_slice(&data, |a: &f64, b| a.total_cmp(b), 1.5, 3.25);
    assert!(data[left] <= 1.5);
    assert!(data[right] > 3.25);
    assert_eq!(
        &data[left..right],
        &[1.0, 1.5, 1.5, 1.5, 1.5, 2.0, 2.5, 3.0, 3.0]
    );
}

#[test]
fn fragment_bucket_slice_includes_every_repeated_prefix() {
    let buckets = [99.96875, 100.0, 100.0, 100.0, 100.03125];
    assert_eq!(fragment_bucket_slice(&buckets, 100.01, 100.02), (1, 4));
}

#[test]
fn fragment_index_preserves_exact_ids_masses_and_ranges() {
    assert_eq!(std::mem::size_of::<PackedFragment>(), 6);
    let parameters = Builder {
        bucket_size: Some(2),
        generate_decoys: Some(false),
        ..Builder::default()
    }
    .make_parameters();
    let mut peptides = parameters
        .peptides_from_tsv("sequence\tprotein\nPEPTIDER\tprotein-a\nSEQUENCEK\tprotein-b\n");
    Parameters::reorder_peptides(&mut peptides);
    let peptides = peptides
        .into_iter()
        .flat_map(|peptide| [peptide.clone(), peptide.clone(), peptide])
        .collect::<Vec<_>>();

    let mut expected = BTreeMap::<u32, Vec<Theoretical>>::new();
    for (peptide_index, peptide) in peptides.iter().enumerate() {
        for mass in preliminary_fragment_masses(&parameters, peptide) {
            expected
                .entry(mass.to_bits() >> FRAGMENT_MASS_SUFFIX_BITS)
                .or_default()
                .push(Theoretical {
                    peptide_index: PeptideIx(peptide_index as u32),
                    fragment_mz: mass,
                });
        }
    }

    let bucket_size = parameters.bucket_size;
    let database = parameters.build_from_peptides(peptides);
    let expected = expected.into_values().flatten().collect::<Vec<_>>();
    let actual = (0..database.buckets().len())
        .flat_map(|bucket| database.fragments.bucket(bucket))
        .collect::<Vec<_>>();
    assert_eq!(actual, expected);
    assert!(
        database.buckets().windows(2).any(|pair| pair[0] == pair[1]),
        "expected a split prefix in {:?}",
        database.buckets()
    );
    assert!((0..database.buckets().len())
        .all(|bucket| database.fragments.bucket(bucket).count() <= bucket_size));
    assert_eq!(
        database.fragments.allocated_bytes(),
        expected.len() * 6 + database.buckets().len() * 12
    );
}

#[test]
fn filtered_targets_receive_paired_decoys_after_selection() {
    let parameters = Builder::default().make_parameters();
    let targets = parameters
        .peptides_from_tsv("sequence\tprotein\nPEPTIDER\tprotein-a\nSEQUENCEK\tprotein-b\n");
    let targets = targets
        .into_iter()
        .filter(|peptide| !peptide.decoy)
        .collect::<Vec<_>>();

    let peptides = parameters.add_reversed_decoys(targets);

    assert_eq!(peptides.iter().filter(|peptide| !peptide.decoy).count(), 2);
    assert_eq!(peptides.iter().filter(|peptide| peptide.decoy).count(), 2);
}

#[test]
fn filtered_decoy_generation_drops_reversals_that_collide_with_targets() {
    let builder = Builder {
        generate_decoys: Some(false),
        ..Builder::default()
    };
    let parameters = builder.make_parameters();
    let targets = parameters
        .peptides_from_tsv("sequence\tprotein\nPEPTIDER\tprotein-a\nPEDITPER\tprotein-b\n");

    let peptides = parameters.add_reversed_decoys(targets);

    assert_eq!(peptides.len(), 2);
    assert!(peptides.iter().all(|peptide| !peptide.decoy));
}

fn digest_group(sequence: &str, position: Position) -> DigestGroup {
    let reference = Digest {
        sequence: sequence.into(),
        position,
        protein: Arc::from(sequence),
        ..Digest::default()
    };
    DigestGroup {
        origins: vec![ProteinOccurrence {
            protein: reference.protein.clone(),
            start: reference.protein_start,
            prev_aa: reference.prev_aa,
            next_aa: reference.next_aa,
            source: None,
            met_clipped: false,
        }],
        reference,
    }
}

#[test]
fn chunk_decoys_are_checked_against_targets_in_other_chunks() {
    let parameters = Builder::default().make_parameters();
    let target_sequences = ["PEPTIDER".into(), "PEDITPER".into()]
        .into_iter()
        .collect::<HashSet<_>>();

    let peptides = parameters.modify_digests_with_target_sequences(
        vec![digest_group("PEPTIDER", Position::Full)],
        &target_sequences,
    );

    assert_eq!(peptides.len(), 1);
    assert!(!peptides[0].decoy);
}

#[test]
fn sequence_hashes_pair_targets_with_their_reversals() {
    for sequence in ["PEPTIDER", "PEDITPER", "AK", "K", "MPEPTIDEK", ""] {
        let reversed = PeptideSequence::from(sequence).reversed_internal();
        let (forward, reverse) = sequence_hashes(sequence.as_bytes());
        assert_eq!(reverse, sequence_hashes(reversed.as_bytes()).0);
        assert_eq!(forward, sequence_hashes(reversed.as_bytes()).1);
    }
    let hashes = (0..200)
        .map(|n| sequence_hashes(format!("PEPTIDE{n}K").as_bytes()).0)
        .collect::<HashSet<_>>();
    assert_eq!(hashes.len(), 200);
}

/// Identity of a generated peptide, independent of chunking.
type PeptideKey = (String, bool, bool, u8, Vec<String>, usize);

fn peptide_keys(mut peptides: Vec<Peptide>) -> Vec<PeptideKey> {
    Parameters::reorder_peptides(&mut peptides);
    peptides
        .iter()
        .map(|peptide| {
            (
                peptide.to_string(),
                peptide.decoy,
                peptide.semi_enzymatic,
                peptide.missed_cleavages,
                peptide.proteins.iter().map(|p| p.to_string()).collect(),
                peptide.protein_sites.len(),
            )
        })
        .collect()
}

#[test]
fn per_protein_expansion_of_unshared_digests_matches_the_whole_database() {
    // PEDITPER is the internal reversal of PEPTIDER, so its decoy collides
    // with a target from another protein; the buckets must keep them together.
    // With FASTA decoys, rev_c supplies decoys that must still be checked
    // against targets in the same bucket.
    let fasta = ">a\nMPEPTIDERSEQMENCEKAMPLIFIERKPEPTIDEKLLMSTK\n\
                 >b\nPEDITPERGGMKWHATEVERKMMKSAMPLEPEPTIDERK\n\
                 >rev_c\nREDITPEPKSEQUENCEKPEPTIDERK\n";
    let configs = [
        (
            true,
            serde_json::json!({
                "enzyme": {"missed_cleavages": 2, "min_len": 3},
                "variable_mods": {"M": [15.9949], "^E": [-18.010565]},
                "static_mods": {"C": 57.021464},
                "max_variable_mods": 2
            }),
        ),
        (
            true,
            serde_json::json!({
                "enzyme": {"missed_cleavages": 1, "min_len": 4, "semi_enzymatic": true},
                "variable_mods": {"M": [15.9949]}
            }),
        ),
        (
            false,
            serde_json::json!({
                "enzyme": {"missed_cleavages": 1, "min_len": 4},
                "variable_mods": {"M": [15.9949]}
            }),
        ),
    ];
    let mut shared_exercised = false;
    for (generate_decoys, config) in configs {
        let fasta = Fasta::parse(fasta.to_string(), "rev_", generate_decoys).unwrap();
        let mut builder: Builder = serde_json::from_value(config).unwrap();
        builder.generate_decoys = Some(generate_decoys);
        let parameters = builder.make_parameters();
        let whole = peptide_keys(parameters.digest(&fasta));
        assert!(whole.iter().any(|key| key.1), "no decoys generated");

        // As the streamed prefilter does: digests whose sequence and decoy
        // sequence are unique are expanded one protein at a time without a
        // target check; the rest are expanded per canonical key.
        let enzyme = parameters.enzyme_parameters();
        let digests = (0..fasta.targets.len())
            .map(|index| fasta.digest_protein(index, &enzyme, None))
            .collect::<Vec<_>>();
        // A sequence, its decoy, and a target equal to that decoy share a
        // canonical key; palindromes always collide with their own decoy.
        let key = |digest: &crate::enzyme::Digest| {
            let (forward, reverse) = sequence_hashes(digest.sequence.as_bytes());
            match generate_decoys {
                true if forward == reverse => u64::MAX,
                true => forward.min(reverse),
                false => forward,
            }
        };
        let mut keys = digests.iter().flatten().map(key).collect::<Vec<_>>();
        keys.sort_unstable();
        let shared = keys
            .windows(2)
            .filter(|pair| pair[0] == pair[1])
            .map(|pair| pair[0])
            .chain(generate_decoys.then_some(u64::MAX))
            .collect::<HashSet<_>>();
        let is_shared = |digest: &crate::enzyme::Digest| shared.contains(&key(digest));

        let mut deferred = Vec::new();
        let mut streamed = parameters.with_digest_expander(|expander| {
            let mut streamed = Vec::new();
            for protein in &digests {
                let mut peptides = Vec::new();
                for digest in protein {
                    if is_shared(digest) {
                        deferred.push(digest.clone());
                        continue;
                    }
                    let group = DigestGroup {
                        origins: vec![ProteinOccurrence::of_protein_digest(digest)],
                        reference: digest.clone(),
                    };
                    peptides.extend(expander.expand(group, None));
                }
                expander.reorder(&mut peptides);
                streamed.extend(peptides);
            }
            streamed
        });
        shared_exercised |= !deferred.is_empty();
        // Each shared key is filtered as its own small database.
        let mut units = BTreeMap::<u64, Vec<_>>::new();
        for digest in deferred {
            units.entry(key(&digest)).or_default().push(digest);
        }
        for unit in units.into_values() {
            streamed.extend(parameters.modify_digests(crate::enzyme::group_protein_digests(unit)));
        }
        parameters.reorder_peptides_with_labels(&mut streamed);
        assert_eq!(
            peptide_keys(streamed),
            whole,
            "generate_decoys {generate_decoys}"
        );
    }
    assert!(shared_exercised, "no shared digests exercised");
}

#[test]
fn chunked_modification_matches_a_single_pass() {
    let fasta = ">a\nMPEPTIDERSEQMENCEKAMPLIFIERKPEPTIDEKLLMSTK\n\
                 >b\nPEDITPERGGMKWHATEVERKMMKSAMPLEPEPTIDERK\n\
                 >c\nMSTYKQNNQMKPEPTIDERKSTSTSYMR\n";
    let fasta = Fasta::parse(fasta.to_string(), "rev_", true).unwrap();
    let builder: Builder = serde_json::from_value(serde_json::json!({
        "enzyme": {"missed_cleavages": 2, "min_len": 3},
        "variable_mods": {
            "M": [15.9949],
            "S": [79.966331],
            "T": [79.966331],
            "Y": [79.966331],
            "^E": [-18.010565]
        },
        "static_mods": {"C": 57.021464},
        "max_variable_mods": 2
    }))
    .unwrap();
    let parameters = builder.make_parameters();
    let digests = parameters.digest_unmodified(&fasta);
    assert!(digests.len() > 8, "too few digest groups to chunk");
    let targets = digests
        .iter()
        .filter(|digest| !digest.reference.decoy)
        .map(|digest| digest.reference.sequence.clone())
        .collect::<HashSet<_>>();
    let single = parameters.modify_digest_chunks(digests.clone(), &targets, usize::MAX);
    for chunk_groups in [0, 1, 2, 7] {
        let chunked = parameters.modify_digest_chunks(digests.clone(), &targets, chunk_groups);
        // Same peptides in the same order, not just the same set.
        assert!(
            chunked == single,
            "chunk size {chunk_groups} changed the output"
        );
    }
}

#[test]
fn chunked_modification_matches_a_single_pass_with_libraries() {
    use crate::modification::{NeutralLossMode, SiteMode, VariableModification};
    use crate::ptm_library::{PtmLibrary, PtmLibrarySite};

    // LLPEPTIDESK is shared by P1 and P3; only P3 lists its T as a site, so
    // library sites must attach per origin.
    let fasta = Fasta::parse(
        ">P1\nMSTKLLPEPTIDESKAAKAPEPTIDERQQQKSSYR\n\
         >P2\nMKSTSYKPEDITPERGGSKWHATEVERK\n\
         >P3\nMGGKLLPEPTIDESKLLSTYRSAMPLER\n"
            .into(),
        "rev_",
        true,
    )
    .unwrap();
    let cleavages = CustomCleavageLibrary::from_tsv(
        "protein\tposition\tcontext\nP1\t22\tAPEPT|IDER\nP2\t11\tPEDIT|PERG\n",
    )
    .unwrap()
    .validate(&fasta)
    .unwrap();
    let phospho = |site_mode| {
        vec![VarModEntry::Detailed(VariableModification {
            search_mode: SearchMode::Database,
            mass: 79.96633,
            max_count: Some(2),
            max_total_count: Some(3),
            name: Some("Phospho".into()),
            neutral_losses: vec![],
            neutral_loss_mode: NeutralLossMode::Optional,
            site_mode,
            channel_offsets: Default::default(),
        })]
    };
    let builder = Builder {
        enzyme: Some(EnzymeBuilder {
            missed_cleavages: Some(2),
            min_len: Some(3),
            max_len: Some(30),
            ..Default::default()
        }),
        peptide_min_mass: Some(0.0),
        max_variable_mods: Some(2),
        max_total_variable_mods: Some(3),
        variable_mods: Some(HashMap::from([
            ("S".into(), phospho(SiteMode::Both)),
            ("T".into(), phospho(SiteMode::Library)),
        ])),
        ..Default::default()
    };
    let mut parameters = builder.make_parameters();
    let site = |protein: &str, position, residue| PtmLibrarySite {
        attachment: Default::default(),
        protein: Arc::from(protein),
        position,
        residue,
        modification: Arc::from("Phospho"),
    };
    parameters.loaded_ptm_library = Some(Arc::new(PtmLibrary::new(vec![
        site("P1", 2, b'T'),
        site("P1", 13, b'S'),
        site("P2", 3, b'T'),
        site("P3", 9, b'T'),
        site("P3", 18, b'T'),
    ])));

    let digests = parameters.digest_unmodified_with_custom_cleavages(&fasta, Some(&cleavages));
    assert!(digests.len() > 8, "too few digest groups to chunk");
    assert!(
        digests.iter().any(|group| group.origins.len() > 1),
        "no shared digest group"
    );
    let targets = digests
        .iter()
        .filter(|digest| !digest.reference.decoy)
        .map(|digest| parameters.decoy_collision_key(&digest.reference.sequence))
        .collect::<HashSet<_>>();
    let single = parameters.modify_digest_chunks(digests.clone(), &targets, usize::MAX);
    assert!(
        single.iter().any(|peptide| !peptide.decoy
            && (0..peptide.sequence.len()).any(|ix| peptide.modification_at(ix) != 0.0)),
        "no modified target peptide"
    );
    for chunk_groups in [0, 1, 2, 7] {
        let chunked = parameters.modify_digest_chunks(digests.clone(), &targets, chunk_groups);
        // Same peptides in the same order, not just the same set.
        assert!(
            chunked == single,
            "chunk size {chunk_groups} changed the output"
        );
    }
}

#[test]
fn modification_variants_share_target_and_decoy_sequence_storage() {
    let builder: Builder = serde_json::from_value(serde_json::json!({
        "variable_mods": {"M": [15.9949]}
    }))
    .unwrap();
    let parameters = builder.make_parameters();
    let peptides = parameters.modify_digests(vec![digest_group("AMPEPTIDER", Position::Full)]);
    let targets = peptides
        .iter()
        .filter(|peptide| !peptide.decoy)
        .collect::<Vec<_>>();
    let decoys = peptides
        .iter()
        .filter(|peptide| peptide.decoy)
        .collect::<Vec<_>>();

    assert!(targets.len() >= 2);
    assert_eq!(targets.len(), decoys.len());
    assert!(targets
        .windows(2)
        .all(|pair| pair[0].sequence.shares_storage_with(&pair[1].sequence)));
    assert!(decoys
        .windows(2)
        .all(|pair| pair[0].sequence.shares_storage_with(&pair[1].sequence)));
}

#[test]
fn structured_variable_mod_config_round_trips() {
    let builder: Builder = serde_json::from_value(serde_json::json!({
        "fasta": "none",
        "static_mods": {
            "Carbamidomethyl": {"mass": 57.0215, "sites": ["C"]}
        },
        "variable_mods": {
            "Oxidation": {"mass": 15.9949, "sites": ["M"]},
            "Acetyl": {
                "mass": 42.0106,
                "max_count": 1,
                "neutral_losses": [17.0265],
                "neutral_loss_mode": "required",
                "sites": ["K"]
            },
            "Methyl": {"mass": 14.0157, "sites": ["K"]}
        },
        "max_variable_mods": 2,
        "max_combinations": 0
    }))
    .unwrap();

    let params = builder.make_parameters();
    assert_eq!(params.max_variable_mods, 2);
    assert_eq!(params.max_combinations, Some(1));

    let mods = params.variable_modifications();
    assert_eq!(mods.len(), 3);
    assert_eq!(mods[0].specificity, ModificationSpecificity::Residue(b'K'));
    assert!((mods[0].modification.mass - 42.0106).abs() < 1e-4);
    assert_eq!(mods[0].max_count, Some(1));
    assert_eq!(mods[0].modification.name.as_deref(), Some("Acetyl"));
    assert_eq!(&*mods[0].modification.neutral_losses, &[17.0265]);
    assert_eq!(mods[1].specificity, ModificationSpecificity::Residue(b'K'));
    assert!((mods[1].modification.mass - 14.0157).abs() < 1e-4);
    assert_eq!(mods[1].max_count, None);
    assert_eq!(mods[2].specificity, ModificationSpecificity::Residue(b'M'));
    assert!((mods[2].modification.mass - 15.9949).abs() < 1e-4);
    assert_eq!(mods[2].max_count, None);

    let serialized = serde_json::to_value(params).unwrap();
    let acetyl = &serialized["variable_mods"]["Acetyl"];
    assert_eq!(acetyl["max_count"], 1);
    assert_eq!(acetyl["neutral_loss_mode"], "required");
    assert_eq!(acetyl["sites"], serde_json::json!(["K"]));
    assert!(serialized["variable_mods"]["Methyl"]
        .get("max_count")
        .is_none());
    assert_eq!(
        serialized["variable_mods"]["Oxidation"]["sites"],
        serde_json::json!(["M"])
    );
    assert_eq!(
        serialized["static_mods"]["Carbamidomethyl"]["sites"],
        serde_json::json!(["C"])
    );
}

#[test]
fn channel_offsets_generate_complete_static_channels() {
    let builder: Builder = serde_json::from_value(serde_json::json!({
        "fasta": "none",
        "generate_decoys": false,
        "static_mods": {
            "SILAC-K": {
                "mass": 0.0,
                "channel_offsets": {"light": 0.0, "heavy": 8.014199},
                "sites": ["K"]
            },
            "SILAC-R": {
                "mass": 0.0,
                "neutral_losses": [17.026549],
                "neutral_loss_mode": "required",
                "channel_offsets": {"light": 0.0, "heavy": 10.008269},
                "sites": ["R"]
            }
        }
    }))
    .unwrap();
    let params = builder.make_parameters();
    params.validate_channels().unwrap();

    let peptides = params.peptides_from_tsv("sequence\nPEPKR\n");
    assert_eq!(peptides.len(), 2);
    let light = peptides
        .iter()
        .find(|peptide| peptide.label_channel.as_deref() == Some("light"))
        .unwrap();
    let heavy = peptides
        .iter()
        .find(|peptide| peptide.label_channel.as_deref() == Some("heavy"))
        .unwrap();
    assert!((heavy.monoisotopic - light.monoisotopic - 18.022468).abs() < 1e-4);
    assert_eq!(light.label_group(), heavy.label_group());
    assert_eq!(heavy.to_string(), "PEPK[SILAC-K]R[SILAC-R]");
    let arg10 = heavy
        .applied_modifications()
        .find(|applied| applied.modification.name.as_deref() == Some("SILAC-R"))
        .unwrap();
    assert_eq!(&*arg10.modification.neutral_losses, &[17.026549]);
    assert_eq!(
        arg10.modification.neutral_loss_mode,
        crate::modification::NeutralLossMode::Required
    );
}

#[test]
fn channels_deduplicate_peptides_without_channel_sites() {
    let builder: Builder = serde_json::from_value(serde_json::json!({
        "generate_decoys": false,
        "static_mods": {
            "SILAC-K": {
                "mass": 0.0,
                "channel_offsets": {"light": 0.0, "heavy": 8.014199},
                "sites": ["K"]
            }
        }
    }))
    .unwrap();
    let params = builder.make_parameters();
    params.validate_channels().unwrap();

    let peptides = params.peptides_from_tsv("sequence\nPEPTIDE\n");
    assert_eq!(peptides.len(), 1);
    assert_eq!(peptides[0].label_channel, None);
}

#[test]
fn variable_channel_offsets_preserve_site_variants_and_shared_light() {
    let builder: Builder = serde_json::from_value(serde_json::json!({
        "generate_decoys": false,
        "max_variable_mods": 2,
        "variable_mods": {"SILAC-K": {
            "mass": 0.0,
            "channel_offsets": {"light": 0.0, "heavy": 8.014199},
            "sites": ["K"]
        }}
    }))
    .unwrap();
    let params = builder.make_parameters();
    params.validate_channels().unwrap();
    let peptides = params.peptides_from_tsv("sequence\nPEPTIDEKK\n");
    assert_eq!(peptides.len(), 4);
    assert_eq!(
        peptides
            .iter()
            .filter(|peptide| peptide.label_channel.as_deref() == Some("light"))
            .count(),
        1
    );
    assert!(peptides
        .iter()
        .all(|peptide| peptide.label_group() == "PEPTIDEKK"));
}

#[test]
fn channel_offsets_add_to_the_modification_base_mass() {
    let builder: Builder = serde_json::from_value(serde_json::json!({
        "generate_decoys": false,
        "static_mods": {
            "TMT-SILAC-K": {
                "mass": 229.162932,
                "channel_offsets": {"light": 0.0, "heavy": 8.014199},
                "sites": ["K"]
            }
        }
    }))
    .unwrap();
    let params = builder.make_parameters();
    params.validate_channels().unwrap();
    let digest = Digest {
        sequence: "PEPTIDEK".into(),
        ..Digest::default()
    };
    let peptides = params.modify_digests(vec![DigestGroup {
        reference: digest.clone(),
        origins: vec![crate::enzyme::ProteinOccurrence {
            protein: digest.protein.clone(),
            start: digest.protein_start,
            prev_aa: None,
            next_aa: None,
            source: None,
            met_clipped: false,
        }],
    }]);
    let heavy = peptides
        .iter()
        .find(|peptide| peptide.label_channel.as_deref() == Some("heavy"))
        .unwrap();
    assert!(heavy.to_string().ends_with("K[TMT-SILAC-K]"));
    assert_eq!(
        heavy
            .applied_modifications()
            .filter(|applied| applied.site == crate::peptide::Site::Sequence(7))
            .count(),
        1
    );
}

#[test]
fn channel_resolved_modification_definitions_are_interned() {
    let builder: Builder = serde_json::from_value(serde_json::json!({
        "generate_decoys": false,
        "static_mods": {
            "Labeled-K": {
                "mass": 229.16293,
                "neutral_losses": [17.02655],
                "neutral_loss_mode": "required",
                "channel_offsets": {"light": 0.0, "heavy": 8.014199},
                "sites": ["K"]
            }
        }
    }))
    .unwrap();
    let params = builder.make_parameters();
    params.validate_channels().unwrap();

    let peptides = params.peptides_from_tsv("sequence\nPEPTIDEK\nAAAAAAK\n");
    let definition = |sequence: &[u8], channel: &str| {
        peptides
            .iter()
            .find(|peptide| {
                peptide.sequence.as_ref() == sequence
                    && peptide.label_channel.as_deref() == Some(channel)
            })
            .unwrap()
            .applied_modifications()
            .find(|applied| applied.modification.name.as_deref() == Some("Labeled-K"))
            .unwrap()
            .modification
    };

    let heavy_peptide = definition(b"PEPTIDEK", "heavy");
    let heavy_alanine = definition(b"AAAAAAK", "heavy");
    let light_peptide = definition(b"PEPTIDEK", "light");
    let light_alanine = definition(b"AAAAAAK", "light");

    assert!(std::ptr::eq(heavy_peptide, heavy_alanine));
    assert!(std::ptr::eq(light_peptide, light_alanine));
    assert!(!std::ptr::eq(heavy_peptide, light_peptide));
    assert_eq!(
        heavy_peptide.mass.to_bits(),
        (229.16293f32 + 8.014199).to_bits()
    );
    assert_eq!(heavy_peptide.name.as_deref(), Some("Labeled-K"));
    assert_eq!(&*heavy_peptide.neutral_losses, &[17.02655]);
    assert_eq!(
        heavy_peptide.neutral_loss_mode,
        crate::modification::NeutralLossMode::Required
    );
    assert_eq!(
        heavy_peptide.channel_offsets["light"].to_bits(),
        0.0f32.to_bits()
    );
    assert_eq!(
        heavy_peptide.channel_offsets["heavy"].to_bits(),
        8.014199f32.to_bits()
    );
}

#[test]
fn ptm_library_configuration_round_trips() {
    let builder: Builder = serde_json::from_value(serde_json::json!({
        "variable_mods": {
            "Phospho": {
                "mass": 79.96633,
                "max_count": 2,
                "site_mode": "both",
                "neutral_losses": [97.9769],
                "sites": ["S"]
            }
        },
        "max_variable_mods": 1,
        "max_total_variable_mods": 3,
        "max_combinations": 1000,
        "ptm_library": {"path": "sites.parquet", "strict": true}
    }))
    .unwrap();
    let params = builder.make_parameters();
    params.validate_ptm_library(&PtmLibrary::default()).unwrap();
    assert_eq!(params.max_variable_mods, 1);
    assert_eq!(params.max_total_variable_mods, 3);
    assert_eq!(params.variable_modifications()[0].site_mode, SiteMode::Both);

    let serialized = serde_json::to_value(params).unwrap();
    assert_eq!(serialized["ptm_library"]["path"], "sites.parquet");
    assert_eq!(serialized["variable_mods"]["Phospho"]["site_mode"], "both");
}

#[test]
fn digestion() {
    let fasta = r#"
        >sp|AAAAA
        MEWKLEQSMREQALLKAQLTQLK
        >sp|BBBBB
        RMEWKLEQSMREQALLKAQLTQLK
        "#;

    let fasta = Fasta::parse(fasta.into(), "rev_", false).unwrap();

    // Make sure that FASTA parsed OK
    assert_eq!(
        fasta.targets,
        vec![
            (
                Arc::from("sp|AAAAA".to_string()),
                "MEWKLEQSMREQALLKAQLTQLK".into()
            ),
            (
                Arc::from("sp|BBBBB".to_string()),
                "RMEWKLEQSMREQALLKAQLTQLK".into()
            ),
        ]
    );

    let params = Parameters {
        bucket_size: 128,
        enzyme: EnzymeBuilder {
            missed_cleavages: Some(1),
            min_len: Some(6),
            max_len: Some(10),
            ..Default::default()
        },
        peptide_min_mass: 150.0,
        peptide_max_mass: 5000.0,
        ion_kinds: vec![Kind::B, Kind::Y],
        min_ion_index: 2,
        fragment_losses: None,
        max_fragment_losses: None,
        static_mods: HashMap::default(),
        variable_mods: [(
            ModificationSpecificity::ProteinN(None),
            vec![VarModEntry::Mass(42.0)],
        )]
        .into_iter()
        .collect(),
        max_variable_mods: 2,
        max_total_variable_mods: 2,
        max_combinations: None,
        ptm_library: None,
        decoy_tag: "rev_".into(),
        generate_decoys: false,
        clip_n_term_met: false,
        fasta: "none".into(),
        expand_ambiguous_residues: false,
        max_ambiguous_variants: 20,
        merge_isoleucine_leucine: false,
        peptides: None,
        custom_cleavage_sites: None,
        prefilter: false,
        prefilter_min_matched_peaks: 1,
        prefilter_max_peaks: None,
        loaded_ptm_library: None,
    };

    let peptides = params.digest(&fasta);

    let expected = [
        "EQALLK",
        "LEQSMR",
        "AQLTQLK",
        "MEWKLEQSMR",
        "[+42]-MEWKLEQSMR",
    ]
    .into_iter()
    .map(String::from)
    .collect::<Vec<_>>();

    let sequences = peptides.iter().map(|p| p.to_string()).collect::<Vec<_>>();
    assert_eq!(expected, sequences);

    // All peptides are shared except for the protein N-term mod
    for peptide in &peptides[..4] {
        assert_eq!(peptide.proteins.len(), 2, "{:?}", peptide);
    }
    // Ensure that this mod is uniquely called as the first protein
    assert_eq!(
        peptides.last().unwrap().proteins.as_slice(),
        &[Arc::<str>::from("sp|AAAAA")]
    );
}

#[test]
fn custom_cleavages_flow_through_modification_and_memory_paths() {
    let fasta = Fasta::parse(">P1\nAAKAPEPTIDERQQQK\n".into(), "rev_", true).unwrap();
    let library =
        CustomCleavageLibrary::from_tsv("protein\tposition\tcontext\nP1\t7\tAPEPT|IDER\n")
            .unwrap()
            .validate(&fasta)
            .unwrap();
    let builder = Builder {
        enzyme: Some(EnzymeBuilder {
            missed_cleavages: Some(1),
            min_len: Some(3),
            max_len: Some(50),
            ..Default::default()
        }),
        generate_decoys: Some(false),
        ..Default::default()
    };
    let parameters = builder.make_parameters();

    let ordinary = parameters.digest(&fasta);
    let custom = parameters.digest_with_custom_cleavages(&fasta, Some(&library));
    let custom_sequences = custom
        .iter()
        .map(|peptide| std::str::from_utf8(&peptide.sequence).unwrap())
        .collect::<Vec<_>>();
    assert!(custom.len() > ordinary.len());
    assert!(custom_sequences.contains(&"APEPT"));
    assert!(custom_sequences.contains(&"IDER"));

    let estimate = parameters.estimate_memory_with_custom_cleavages(&fasta, Some(&library));
    assert!(estimate.modified_peptides as usize >= custom.len());
    assert!(estimate.modified_peptides > parameters.estimate_memory(&fasta).modified_peptides);
}

#[test]
fn peptides_with_massless_residues_are_skipped_at_digestion() {
    let fasta = Fasta::parse(
        ">P1\nPEPTIDEKAAXAWWWWWKGGBGGRSSZSSKTTJTTRLLLLLK\n>P2\nMSSWWHHK\n".into(),
        "rev_",
        true,
    )
    .unwrap();
    assert_eq!(fasta.targets.len(), 2);
    let library = CustomCleavageLibrary::from_tsv("protein\tposition\tcontext\nP1\t11\tAAXA|WW\n")
        .unwrap()
        .validate(&fasta)
        .unwrap();
    let builder = Builder {
        enzyme: Some(EnzymeBuilder {
            missed_cleavages: Some(1),
            min_len: Some(2),
            max_len: Some(50),
            ..Default::default()
        }),
        ..Default::default()
    };
    let parameters = builder.make_parameters();

    let custom = parameters.digest_with_custom_cleavages(&fasta, Some(&library));
    assert!(custom
        .iter()
        .any(|peptide| &peptide.sequence[..] == b"WWWWWK"));
    for peptides in [parameters.digest(&fasta), custom] {
        let sequences = peptides
            .iter()
            .map(|peptide| {
                (
                    std::str::from_utf8(&peptide.sequence).unwrap(),
                    peptide.decoy,
                )
            })
            .collect::<HashSet<_>>();
        assert!(sequences.contains(&("PEPTIDEK", false)));
        assert!(sequences.contains(&("LLLLLK", false)));
        assert!(sequences.contains(&("MSSWWHHK", false)));
        // J (Ile or Leu) has the I/L mass and is kept as written.
        assert!(sequences.contains(&("TTJTTR", false)));
        assert!(sequences.iter().any(|(_, decoy)| *decoy));
        for peptide in &peptides {
            assert!(
                peptide
                    .sequence
                    .iter()
                    .all(|residue| crate::mass::VALID_AA.contains(residue) || *residue == b'J'),
                "{peptide:?}"
            );
            assert!(peptide.monoisotopic > 0.0);
        }
    }

    let database = parameters.build(fasta);
    assert!(!database.peptides.is_empty());
    assert!(database.peptides.iter().all(|peptide| peptide
        .sequence
        .iter()
        .all(|residue| crate::mass::VALID_AA.contains(residue) || *residue == b'J')));
}

#[test]
fn estimates_variable_modification_expansion_before_allocation() {
    let builder = Builder {
        enzyme: Some(EnzymeBuilder {
            cleave_at: Some("$".into()),
            min_len: Some(1),
            max_len: Some(50),
            ..Default::default()
        }),
        variable_mods: Some(
            [(
                "S".to_string(),
                vec![VarModEntry::Mass(79.9663), VarModEntry::Mass(80.0)],
            )]
            .into_iter()
            .collect(),
        ),
        max_variable_mods: Some(3),
        generate_decoys: Some(false),
        ..Default::default()
    };
    let parameters = builder.make_parameters();
    let fasta = Fasta::parse(">protein\nSSSSSSSSSS\n".into(), "rev_", false).unwrap();

    let estimate = parameters.estimate_memory(&fasta);

    // 1 + C(10,1)*2 + C(10,2)*2^2 + C(10,3)*2^3
    assert_eq!(estimate.unmodified_peptides, 1);
    assert_eq!(estimate.modified_peptides, 1_161);
    assert_eq!(estimate.fragments, 1_161 * 14);
    assert!(estimate.unmodified_peak_bytes > 0);
    assert!(estimate.modified_peak_bytes > estimate.unmodified_peak_bytes);
    assert!(estimate.fragment_peak_bytes > estimate.modified_peak_bytes);

    let digests = parameters.digest_unmodified(&fasta);
    let modification_estimate = parameters.estimate_modified_memory(&digests);
    assert_eq!(modification_estimate.modified_peptides, 1_161);
    assert_eq!(parameters.modify_digests(digests).len(), 1_161);
}

/// Proteins are estimated in parallel; the totals equal the sum over any
/// split of the FASTA, as the per-chunk prefilter estimates rely on.
#[test]
fn memory_estimate_totals_are_additive_over_protein_chunks() {
    let builder = Builder {
        enzyme: Some(EnzymeBuilder {
            min_len: Some(4),
            max_len: Some(30),
            missed_cleavages: Some(1),
            ..Default::default()
        }),
        variable_mods: Some(
            [
                ("M".to_string(), vec![VarModEntry::Mass(15.9949)]),
                ("S".to_string(), vec![VarModEntry::Mass(79.9663)]),
            ]
            .into_iter()
            .collect(),
        ),
        max_variable_mods: Some(2),
        ..Default::default()
    };
    let parameters = builder.make_parameters();
    let fasta = (0..64)
        .map(|idx| {
            let residues = "MSKPEPTIDESRAGMSLKWVTFISLLFLFSSAYSR";
            let rotation = idx % residues.len();
            format!(
                ">sp|P{idx:05}|TEST\n{}{}\n",
                &residues[rotation..],
                &residues[..rotation]
            )
        })
        .collect::<String>();
    let fasta = Fasta::parse(fasta, "rev_", true).unwrap();

    let full = parameters.estimate_memory(&fasta);
    assert!(full.unmodified_peptides > 64);
    assert!(full.modified_peptides > full.unmodified_peptides);
    assert_eq!(
        format!("{:?}", parameters.estimate_memory(&fasta)),
        format!("{full:?}")
    );

    let chunks = fasta
        .iter_chunks(7)
        .map(|chunk| parameters.estimate_memory(&chunk))
        .collect::<Vec<_>>();
    let sum = |field: fn(&DatabaseMemoryEstimate) -> u64| chunks.iter().map(field).sum::<u64>();
    assert_eq!(sum(|e| e.unmodified_peptides), full.unmodified_peptides);
    assert_eq!(sum(|e| e.modified_peptides), full.modified_peptides);
    assert_eq!(sum(|e| e.fragments), full.fragments);
}

#[test]
fn protein_site_library_adds_targeted_combinations() {
    use crate::modification::{NeutralLossMode, SiteMode, VariableModification};
    use crate::ptm_library::{PtmLibrary, PtmLibrarySite};

    let builder = Builder {
        enzyme: Some(EnzymeBuilder {
            min_len: Some(1),
            max_len: Some(20),
            ..Default::default()
        }),
        peptide_min_mass: Some(0.0),
        generate_decoys: Some(false),
        max_variable_mods: Some(1),
        max_total_variable_mods: Some(2),
        variable_mods: Some(HashMap::from([(
            "S".into(),
            vec![VarModEntry::Detailed(VariableModification {
                search_mode: SearchMode::Database,
                mass: 79.96633,
                max_count: Some(2),
                max_total_count: None,
                name: Some("Phospho".into()),
                neutral_losses: vec![97.9769],
                neutral_loss_mode: NeutralLossMode::Optional,
                site_mode: SiteMode::Both,
                channel_offsets: Default::default(),
            })],
        )])),
        ..Default::default()
    };
    let mut parameters = builder.make_parameters();
    parameters.loaded_ptm_library = Some(Arc::new(PtmLibrary::new(vec![
        PtmLibrarySite {
            attachment: Default::default(),
            protein: Arc::from("P1"),
            position: 1,
            residue: b'S',
            modification: Arc::from("Phospho"),
        },
        PtmLibrarySite {
            attachment: Default::default(),
            protein: Arc::from("P1"),
            position: 2,
            residue: b'S',
            modification: Arc::from("Phospho"),
        },
    ])));
    let fasta = Fasta::parse(">P1\nMSSK\n".into(), "rev_", false).unwrap();

    let peptides = parameters.digest(&fasta);
    assert!(peptides
        .iter()
        .any(|peptide| { peptide.modification_at(1) != 0.0 && peptide.modification_at(2) != 0.0 }));
}

#[test]
fn library_sites_from_different_proteins_are_not_combined() {
    use crate::modification::{NeutralLossMode, SiteMode, VariableModification};
    use crate::ptm_library::{PtmLibrary, PtmLibrarySite};

    let builder = Builder {
        enzyme: Some(EnzymeBuilder {
            min_len: Some(1),
            max_len: Some(20),
            ..Default::default()
        }),
        peptide_min_mass: Some(0.0),
        generate_decoys: Some(false),
        max_variable_mods: Some(1),
        max_total_variable_mods: Some(2),
        variable_mods: Some(HashMap::from([(
            "S".into(),
            vec![VarModEntry::Detailed(VariableModification {
                search_mode: SearchMode::Database,
                mass: 79.96633,
                max_count: Some(2),
                max_total_count: None,
                name: Some("Phospho".into()),
                neutral_losses: vec![],
                neutral_loss_mode: NeutralLossMode::Optional,
                site_mode: SiteMode::Library,
                channel_offsets: Default::default(),
            })],
        )])),
        ..Default::default()
    };
    let mut parameters = builder.make_parameters();
    parameters.loaded_ptm_library = Some(Arc::new(PtmLibrary::new(vec![
        PtmLibrarySite {
            attachment: Default::default(),
            protein: Arc::from("P1"),
            position: 1,
            residue: b'S',
            modification: Arc::from("Phospho"),
        },
        PtmLibrarySite {
            attachment: Default::default(),
            protein: Arc::from("P2"),
            position: 2,
            residue: b'S',
            modification: Arc::from("Phospho"),
        },
    ])));
    let fasta = Fasta::parse(">P1\nMSSK\n>P2\nMSSWWHHK\n".into(), "rev_", false).unwrap();

    let peptides = parameters.digest(&fasta);
    assert!(!peptides
        .iter()
        .any(|peptide| { peptide.modification_at(1) != 0.0 && peptide.modification_at(2) != 0.0 }));
    let first_site = peptides
        .iter()
        .find(|peptide| peptide.modification_at(1) != 0.0)
        .unwrap();
    assert_eq!(first_site.proteins.as_slice(), &[Arc::from("P1")]);
}

#[test]
fn mass_offset_validation_rejects_ambiguous_definitions() {
    use crate::modification::{NeutralLossMode, VariableModification};
    let entry = |search_mode, site_mode, name: Option<&str>| {
        VarModEntry::Detailed(VariableModification {
            mass: 79.966_33,
            max_count: Some(1),
            max_total_count: None,
            name: name.map(str::to_string),
            neutral_losses: Vec::new(),
            neutral_loss_mode: NeutralLossMode::Optional,
            site_mode,
            search_mode,
            channel_offsets: Default::default(),
        })
    };
    let validate = |mods: Vec<(&str, VarModEntry)>, library: &PtmLibrary| {
        let mut variable_mods: HashMap<String, Vec<VarModEntry>> = HashMap::new();
        for (site, entry) in mods {
            variable_mods.entry(site.into()).or_default().push(entry);
        }
        Builder {
            variable_mods: Some(variable_mods),
            ..Default::default()
        }
        .make_parameters()
        .validate_ptm_library(library)
    };
    let empty = PtmLibrary::default();
    let offset = SearchMode::MassOffset;
    let database = SearchMode::Database;
    let exhaustive = SiteMode::Exhaustive;

    assert!(validate(
        vec![
            ("S", entry(offset, exhaustive, Some("Phospho"))),
            ("T", entry(offset, exhaustive, Some("Phospho"))),
        ],
        &empty
    )
    .is_ok());
    let mixed = validate(
        vec![
            ("S", entry(offset, exhaustive, Some("Phospho"))),
            ("T", entry(database, exhaustive, Some("Phospho"))),
        ],
        &empty,
    )
    .unwrap_err();
    assert!(mixed.contains("cannot use both"), "{mixed}");
    assert!(validate(
        vec![
            ("S", entry(offset, exhaustive, Some("Phospho"))),
            ("T", entry(offset, SiteMode::Both, Some("Phospho"))),
        ],
        &empty
    )
    .unwrap_err()
    .contains("inconsistent"));
    assert!(validate(vec![("S", entry(offset, SiteMode::Library, None))], &empty).is_err());

    let library = PtmLibrary::new(vec![crate::ptm_library::PtmLibrarySite {
        attachment: Default::default(),
        protein: "P1".into(),
        position: 0,
        residue: b'S',
        modification: "Phospho".into(),
    }]);
    assert!(validate(
        vec![("S", entry(offset, SiteMode::Library, Some("Phospho")))],
        &library
    )
    .is_ok());
    assert!(validate(
        vec![("S", entry(offset, exhaustive, Some("Phospho")))],
        &library
    )
    .unwrap_err()
    .contains("site_mode"));
}

fn positional_parameters(keys: &[&str], mode: &str, static_mod: bool) -> Parameters {
    let sites = keys
        .iter()
        .map(|key| {
            key.parse::<ModificationSpecificity>()
                .unwrap()
                .explicit_name()
        })
        .collect::<Vec<_>>();
    let modifications = if static_mod {
        serde_json::json!({"Acetyl": {"mass": 42.0106, "sites": sites}})
    } else {
        serde_json::json!({"Acetyl": {"mass": 42.0106, "max_count": 1, "search_mode": mode, "sites": sites}})
    };
    serde_json::from_value::<Builder>(serde_json::json!({
        if static_mod { "static_mods" } else { "variable_mods" }: modifications,
        "generate_decoys": false, "peptide_min_mass": 0, "peptide_max_mass": 100000,
        "max_variable_mods": 3,
        "enzyme": {"min_len": 1, "max_len": 50, "missed_cleavages": 0}
    }))
    .unwrap()
    .make_parameters()
}

fn positional_digest(sequence: &str, position: Position) -> Digest {
    Digest {
        sequence: sequence.into(),
        position,
        protein: "H3".into(),
        ..Default::default()
    }
}

#[test]
fn internal_placement_static_variable_and_model_counts() {
    for sequence in ["K", "KK", "KKK", "KAKAK"] {
        let internal = sequence
            .bytes()
            .enumerate()
            .filter(|(i, residue)| *residue == b'K' && *i > 0 && *i < sequence.len() - 1)
            .map(|(i, _)| i)
            .collect::<Vec<_>>();
        for static_mod in [false, true] {
            let parameters = positional_parameters(&["~K"], "database", static_mod);
            let digest = positional_digest(sequence, Position::Internal);
            let estimate = parameters
                .variable_variant_count(&digest, &[crate::enzyme::ProteinOccurrence::of(&digest)]);
            let variants = parameters.modify_digests(group_digests(vec![digest]));
            assert_eq!(
                variants.len(),
                if static_mod { 1 } else { 1 + internal.len() }
            );
            assert!(estimate >= variants.len() as u64);
            for peptide in &variants {
                assert_eq!(peptide.modification_at(0), 0.0);
                assert_eq!(peptide.modification_at(sequence.len() - 1), 0.0);
                let count = peptide.applied_modifications().len();
                assert_eq!(
                    peptide.modification_count(ModificationSpecificity::Internal(b'K'), 42.0106),
                    count
                );
                assert_eq!(
                    crate::ml::retention_model::variable_mod_count(
                        peptide,
                        ModificationSpecificity::Internal(b'K'),
                        42.0106
                    ),
                    count as f64
                );
            }
        }
    }
}

#[test]
fn h3k9_combined_keys_share_limit_and_protein_terminal_exception() {
    let parameters = positional_parameters(&["^K", "~K", "]K"], "database", false);
    let variants = parameters.modify_digests(group_digests(vec![positional_digest(
        "KSTGGKAPR",
        Position::Internal,
    )]));
    assert_eq!(variants.len(), 3);
    assert!(variants.iter().any(|p| p.modification_at(0) != 0.0));
    assert!(variants.iter().any(|p| p.modification_at(5) != 0.0));
    assert!(variants
        .iter()
        .all(|p| p.applied_modifications().len() <= 1));
    for (position, count) in [
        (Position::Internal, 2),
        (Position::Cterm, 3),
        (Position::Full, 3),
    ] {
        let variants =
            parameters.modify_digests(group_digests(vec![positional_digest("KAK", position)]));
        assert_eq!(variants.len(), count);
    }
    let overlap = positional_parameters(&["^K", "$K", "]K"], "database", false);
    let variants =
        overlap.modify_digests(group_digests(vec![positional_digest("K", Position::Full)]));
    assert_eq!(variants.len(), 2);
}

#[test]
fn positional_mass_offset_sites_match_indexed_variants() {
    for keys in [&["~K"][..], &["^K", "~K", "]K"][..]] {
        for position in [Position::Internal, Position::Cterm, Position::Full] {
            for sequence in ["K", "KK", "KAKAK", "KSTGGKAPR"] {
                let digest = positional_digest(sequence, position);
                let indexed = positional_parameters(keys, "database", false);
                let expected = indexed.modify_digests(group_digests(vec![digest.clone()]));
                let offset = positional_parameters(keys, "mass_offset", false);
                let peptides = offset.modify_digests(group_digests(vec![digest]));
                let database = offset.build_from_peptides(peptides);
                let base = &database.peptides[0];
                let rule = &database.mass_offsets[0];
                let mut actual = vec![base.clone()];
                actual.extend(
                    database
                        .mass_offset_sites(base, rule)
                        .into_iter()
                        .map(|site| base.with_mass_offset(site, &rule.definition)),
                );
                let strings = |peptides: &[Peptide]| {
                    peptides
                        .iter()
                        .map(ToString::to_string)
                        .collect::<std::collections::BTreeSet<_>>()
                };
                assert_eq!(
                    strings(&actual),
                    strings(&expected),
                    "{sequence}, {position:?}, {keys:?}"
                );
            }
        }
    }
}

#[test]
fn internal_library_sites_respect_position_in_both_search_modes() {
    use crate::ptm_library::PtmLibrarySite;
    for mode in ["database", "mass_offset"] {
        let mut parameters = positional_parameters(&["~K"], mode, false);
        if let VarModEntry::Detailed(entry) = &mut parameters
            .variable_mods
            .get_mut(&ModificationSpecificity::Internal(b'K'))
            .unwrap()[0]
        {
            entry.site_mode = SiteMode::Library;
        }
        parameters.loaded_ptm_library = Some(Arc::new(PtmLibrary::new(
            [0, 2, 4]
                .map(|position| PtmLibrarySite {
                    attachment: Default::default(),
                    protein: "H3".into(),
                    position,
                    residue: b'K',
                    modification: "Acetyl".into(),
                })
                .to_vec(),
        )));
        let mut digest = positional_digest("KAKAK", Position::Full);
        digest.protein_start = Some(0);
        let peptides = parameters.modify_digests(group_digests(vec![digest]));
        if mode == "database" {
            assert_eq!(peptides.len(), 2);
            assert!(peptides
                .iter()
                .all(|p| p.modification_at(0) == 0.0 && p.modification_at(4) == 0.0));
        } else {
            let database = parameters.build_from_peptides(peptides);
            assert_eq!(
                database.mass_offset_sites(&database.peptides[0], &database.mass_offsets[0]),
                vec![Site::Sequence(2)]
            );
        }
    }
}

#[test]
fn internal_label_channels_keep_terminal_residues_unmodified() {
    let parameters = serde_json::from_value::<Builder>(serde_json::json!({
        "static_mods":{"Label":{"mass":0.0,"channel_offsets":{"light":0.0,"heavy":8.0},"sites":["internal_residue:K"]}},
        "generate_decoys":false,"peptide_min_mass":0
    })).unwrap().make_parameters();
    let variants = parameters.modify_digests(group_digests(vec![positional_digest(
        "KAKAK",
        Position::Full,
    )]));
    assert_eq!(variants.len(), 2);
    assert!(variants
        .iter()
        .all(|p| p.modification_at(0) == 0.0 && p.modification_at(4) == 0.0));
    assert!(variants.iter().any(|p| p.modification_at(2) == 8.0));
}

#[test]
fn named_modifications_share_limits_and_preserve_terminal_identity() {
    let config = serde_json::json!({
        "variable_mods": {"Acetyl":{"mass":42.010565,"sites":["first_residue:K","internal_residue:K","peptide_n_term:K"],"max_count":1}},
        "generate_decoys":false,"peptide_min_mass":0,"max_variable_mods":3
    });
    let params = serde_json::from_value::<Builder>(config)
        .unwrap()
        .make_parameters();
    params.validate_ptm_library(&PtmLibrary::default()).unwrap();
    let output = serde_json::to_value(&params).unwrap();
    assert_eq!(
        output["variable_mods"]["Acetyl"]["sites"]
            .as_array()
            .unwrap()
            .len(),
        3
    );
    let round_trip = serde_json::from_value::<Builder>(output)
        .unwrap()
        .make_parameters();
    let digest = Digest {
        sequence: "KAKAK".into(),
        position: Position::Internal,
        ..Default::default()
    };
    let variants = round_trip.modify_digests(group_digests(vec![digest]));
    assert_eq!(variants.len(), 4);
    assert!(variants
        .iter()
        .all(|p| p.applied_modifications().count() <= 1));
    assert!(variants
        .iter()
        .any(|p| p.nterm.is_some() && p.modification_at(0) == 0.0));
    assert!(variants
        .iter()
        .any(|p| p.nterm.is_none() && p.modification_at(0) > 0.0));
}

#[test]
fn typed_library_separates_terminal_from_boundary_residue_in_both_modes() {
    use crate::ptm_library::Attachment;
    for mode in ["database", "mass_offset"] {
        for (attachment, expected) in [
            (Attachment::Residue, Site::Sequence(0)),
            (Attachment::PeptideNTerm, Site::Nterm),
        ] {
            let mut params = serde_json::from_value::<Builder>(serde_json::json!({
                "variable_mods":{"Acetyl":{"mass":42.010565,"sites":["first_residue:K","peptide_n_term:K"],"site_mode":"library","max_count":1,"search_mode":mode}},
                "ptm_library":{"path":"unused.tsv"},"generate_decoys":false,"peptide_min_mass":0
            })).unwrap().make_parameters();
            params.loaded_ptm_library = Some(Arc::new(PtmLibrary::new(vec![
                crate::ptm_library::PtmLibrarySite {
                    attachment,
                    protein: "P1".into(),
                    position: 8,
                    residue: b'K',
                    modification: "Acetyl".into(),
                },
            ])));
            params
                .validate_ptm_library(params.loaded_ptm_library.as_ref().unwrap())
                .unwrap();
            let digest = Digest {
                sequence: "KAAAK".into(),
                protein: "P1".into(),
                protein_start: Some(8),
                position: Position::Internal,
                ..Default::default()
            };
            let peptide = Peptide::try_from(digest.clone()).unwrap();
            let allowed = params
                .localization_rules(&peptide)
                .into_iter()
                .flat_map(|r| r.sites)
                .collect::<Vec<_>>();
            assert_eq!(allowed, vec![expected]);
            let peptides = params.modify_digests(group_digests(vec![digest]));
            if mode == "database" {
                let applied = peptides
                    .iter()
                    .flat_map(|p| p.applied_modifications().map(|m| m.site))
                    .collect::<Vec<_>>();
                assert_eq!(applied, vec![expected]);
            } else {
                let db = params.build_from_peptides(peptides);
                assert_eq!(
                    db.mass_offset_sites(&db.peptides[0], &db.mass_offsets[0]),
                    vec![expected]
                );
            }
        }
    }
}

#[test]
fn named_static_modifications_apply_terminal_and_residue_sites_separately() {
    let params = serde_json::from_value::<Builder>(serde_json::json!({
        "static_mods":{"Label":{"mass":42.0,"sites":["peptide_n_term:K","first_residue:K","K"]}},
        "generate_decoys":false,"peptide_min_mass":0
    }))
    .unwrap()
    .make_parameters();
    let digest = Digest {
        sequence: "KAK".into(),
        position: Position::Internal,
        ..Default::default()
    };
    let peptides = params.modify_digests(group_digests(vec![digest]));
    assert_eq!(peptides.len(), 1);
    assert_eq!(peptides[0].nterm, Some(42.0));
    assert_eq!(peptides[0].modification_at(0), 42.0);
    assert_eq!(peptides[0].modification_at(2), 42.0);
    assert_eq!(peptides[0].applied_modifications().count(), 3);
}

#[test]
fn label_group_recovers_base_definitions_with_nonzero_masses() {
    // `(mass + offset) - offset` is not exact in f32, so recovering the base
    // label by recomputing its mass used to panic for these definitions.
    for (residue, sequence, mass, heavy) in [
        ("K", "PEPKR", 28.0313_f32, 4.025107_f32),
        ("K", "PEPKR", 28.0313, 8.014199),
        ("peptide_n_term", "PEPKR", 28.0313, 4.025107),
        ("C", "PEPCR", 57.021464, 10.008269),
        ("M", "PEPMR", 15.9949, 6.020129),
    ] {
        let builder: Builder = serde_json::from_value(serde_json::json!({
            "generate_decoys": false,
            "static_mods": {
                "Label": {
                    "mass": mass,
                    "channel_offsets": {"light": 0.0, "heavy": heavy},
                    "sites": [residue]
                }
            }
        }))
        .unwrap();
        let parameters = builder.make_parameters();
        parameters.validate_channels().unwrap();
        let peptides = parameters.peptides_from_tsv(&format!("sequence\n{sequence}\n"));
        let db = parameters.build_from_peptides(peptides);
        let groups = db
            .peptides
            .iter()
            .map(|peptide| {
                (
                    peptide.label_channel.as_deref().unwrap().to_string(),
                    peptide.label_group(),
                )
            })
            .collect::<BTreeMap<_, _>>();
        assert_eq!(groups.len(), 2, "{residue} {mass} {heavy}");
        assert_eq!(groups["light"], groups["heavy"], "{residue} {mass} {heavy}");
    }

    let builder: Builder = serde_json::from_value(serde_json::json!({
        "generate_decoys": false,
        "variable_mods": {
            "Label": {
                "mass": 15.9949,
                "channel_offsets": {"light": 0.0, "heavy": 4.025107},
                "sites": ["M"]
            }
        }
    }))
    .unwrap();
    let parameters = builder.make_parameters();
    parameters.validate_channels().unwrap();
    let peptides = parameters.peptides_from_tsv("sequence\nPEPMR\n");
    let db = parameters.build_from_peptides(peptides);
    let modified = db
        .peptides
        .iter()
        .map(|peptide| peptide.label_group())
        .filter(|group| group.contains('['))
        .collect::<HashSet<_>>();
    assert_eq!(modified.len(), 1, "{modified:?}");
}

#[test]
fn required_neutral_loss_fragment_shift_matches_preliminary_index() {
    // Loss order in the configuration must not matter: the preliminary index
    // keeps the smallest total loss, so the offset shift must use it too.
    for losses in [[97.9769_f32, 18.0106], [18.0106, 97.9769]] {
        let definition = Arc::new(ModificationDefinition {
            mass: 79.96633,
            name: Some("Phospho".into()),
            neutral_losses: Arc::from(losses.as_slice()),
            neutral_loss_mode: crate::modification::NeutralLossMode::Required,
            channel_offsets: Arc::default(),
        });
        let offset = MassOffset {
            definition: definition.clone(),
            specificities: vec![ModificationSpecificity::Residue(b'S')],
            site_mode: SiteMode::Exhaustive,
        };
        assert!((offset.fragment_shift() - (79.96633 - 18.0106)).abs() < 1e-4);

        let parameters = Builder {
            generate_decoys: Some(false),
            ..Default::default()
        }
        .make_parameters();
        let peptide = parameters
            .peptides_from_tsv("sequence\nPEPSTIDEK\n")
            .pop()
            .unwrap();
        let modified = peptide.with_mass_offset(Site::Sequence(3), &definition);
        let base = preliminary_fragment_masses(&parameters, &peptide).collect::<Vec<_>>();
        let shifted = preliminary_fragment_masses(&parameters, &modified).collect::<Vec<_>>();
        assert_eq!(base.len(), shifted.len());
        let observed = base
            .iter()
            .zip(&shifted)
            .map(|(base, shifted)| shifted - base)
            .filter(|delta| delta.abs() > 1e-3)
            .collect::<Vec<_>>();
        assert!(!observed.is_empty());
        for delta in observed {
            assert!((delta - offset.fragment_shift()).abs() < 1e-3, "{delta}");
        }
    }
}

fn histone_parameters(
    max_count: usize,
    max_total_count: Option<usize>,
    max_new: usize,
    max_total: usize,
) -> Parameters {
    use crate::ptm_library::{PtmLibrary, PtmLibrarySite};

    let mut acetyl = serde_json::json!({
        "mass": 42.010565, "sites": ["K"], "max_count": max_count, "site_mode": "both"
    });
    if let Some(max_total_count) = max_total_count {
        acetyl["max_total_count"] = max_total_count.into();
    }
    let mut parameters = serde_json::from_value::<Builder>(serde_json::json!({
        "variable_mods": {
            "Acetyl": acetyl,
            "Methyl": {"mass": 14.01565, "sites": ["K", "R"], "max_count": 1, "site_mode": "both"}
        },
        "generate_decoys": false, "peptide_min_mass": 0, "peptide_max_mass": 100000,
        "max_variable_mods": max_new, "max_total_variable_mods": max_total,
        "enzyme": {"min_len": 1, "max_len": 50, "missed_cleavages": 0, "cleave_at": "$"}
    }))
    .unwrap()
    .make_parameters();
    let site = |position: u32, modification: &str| PtmLibrarySite {
        attachment: Default::default(),
        protein: Arc::from("H4"),
        position,
        residue: b'K',
        modification: Arc::from(modification),
    };
    // GKGGKGLGKGGAKKR: K4, K8 and K12 acetyl and K4 methyl are known.
    parameters.loaded_ptm_library = Some(Arc::new(PtmLibrary::new(vec![
        site(4, "Acetyl"),
        site(8, "Acetyl"),
        site(12, "Acetyl"),
        site(4, "Methyl"),
    ])));
    parameters
}

fn acetyl_sites(peptide: &Peptide, library: &[usize]) -> (usize, usize) {
    let sites = peptide
        .applied_modifications()
        .filter(|applied| applied.modification.name.as_deref() == Some("Acetyl"))
        .filter_map(|applied| match applied.site {
            crate::peptide::Site::Sequence(i) => Some(i as usize),
            _ => None,
        })
        .collect::<Vec<_>>();
    let known = sites.iter().filter(|site| library.contains(site)).count();
    (known, sites.len() - known)
}

#[test]
fn max_total_count_caps_library_and_max_count_caps_new_sites() {
    let fasta = Fasta::parse(">H4\nGKGGKGLGKGGAKKR\n".into(), "rev_", false).unwrap();
    // Sequence indices of the library acetyl sites (protein positions 4, 8, 12).
    let library = [4, 8, 12];

    // Default: max_total_count falls back to max_count, so a library acetyl
    // uses the only acetyl slot (Beta 9 behaviour).
    let peptides = histone_parameters(1, None, 1, 3).digest(&fasta);
    assert!(peptides.iter().all(|p| {
        let (known, new) = acetyl_sites(p, &library);
        known + new <= 1
    }));

    // Library acetyls are limited by max_total_count; new acetyls by max_count.
    let peptides = histone_parameters(1, Some(3), 1, 3).digest(&fasta);
    assert!(peptides.iter().any(|p| acetyl_sites(p, &library) == (3, 0)));
    assert!(peptides.iter().any(|p| acetyl_sites(p, &library) == (2, 1)));
    assert!(peptides.iter().all(|p| acetyl_sites(p, &library).1 <= 1));

    // Two new acetyls need both max_count and max_variable_mods of two.
    let peptides = histone_parameters(1, Some(3), 2, 3).digest(&fasta);
    assert!(!peptides.iter().any(|p| acetyl_sites(p, &library).1 == 2));
    let peptides = histone_parameters(2, Some(3), 2, 3).digest(&fasta);
    assert!(peptides.iter().any(|p| acetyl_sites(p, &library).1 == 2));
}

#[test]
fn memory_estimate_counts_per_modification_caps_exactly() {
    let fasta = Fasta::parse(">H4\nGKGGKGLGKGGAKKR\n".into(), "rev_", false).unwrap();
    for (max_count, max_total_count) in [
        (1, None),
        (1, Some(2)),
        (1, Some(3)),
        (2, Some(3)),
        (3, None),
    ] {
        for (max_new, max_total) in [(1, 1), (1, 3), (2, 3), (3, 5)] {
            let parameters = histone_parameters(max_count, max_total_count, max_new, max_total);
            let generated = parameters.digest(&fasta).len() as u64;
            let estimated = parameters.estimate_memory(&fasta).modified_peptides;
            assert_eq!(
                estimated, generated,
                "max_count {max_count}, max_total_count {max_total_count:?}, \
                 max_variable_mods {max_new}, max_total_variable_mods {max_total}"
            );
        }
    }
}

#[test]
fn max_total_count_must_cover_max_count() {
    let error = serde_json::from_value::<Builder>(serde_json::json!({
        "variable_mods": {
            "Acetyl": {"mass": 42.010565, "sites": ["K"], "max_count": 2, "max_total_count": 1}
        }
    }))
    .err()
    .expect("max_total_count below max_count must be rejected");
    assert!(error
        .to_string()
        .contains("max_total_count must be at least max_count"));
}

fn json_parameters(value: serde_json::Value) -> Parameters {
    serde_json::from_value::<Builder>(value)
        .unwrap()
        .make_parameters()
}

#[test]
fn make_parameters_normalizes_modification_limits() {
    let parameters = Builder::default().make_parameters();
    assert_eq!(parameters.max_variable_mods, 2);
    assert_eq!(parameters.max_total_variable_mods, 2);
    assert_eq!(parameters.max_combinations, None);
    assert!(parameters.generate_decoys);
    assert_eq!(parameters.decoy_tag, "rev_");
    assert_eq!(parameters.bucket_size, 8192);

    let parameters = Builder {
        max_variable_mods: Some(0),
        max_total_variable_mods: Some(0),
        max_combinations: Some(0),
        bucket_size: Some(1000),
        ..Default::default()
    }
    .make_parameters();
    assert_eq!(parameters.max_variable_mods, 1);
    assert_eq!(parameters.max_total_variable_mods, 1);
    assert_eq!(parameters.max_combinations, Some(1));
    assert_eq!(parameters.bucket_size, 1024);

    // The total budget can never fall below the per-peptide new-site budget.
    let parameters = Builder {
        max_variable_mods: Some(3),
        max_total_variable_mods: Some(1),
        ..Default::default()
    }
    .make_parameters();
    assert_eq!(parameters.max_total_variable_mods, 3);
    let parameters = Builder {
        max_variable_mods: Some(1),
        max_total_variable_mods: Some(4),
        ..Default::default()
    }
    .make_parameters();
    assert_eq!(
        (
            parameters.max_variable_mods,
            parameters.max_total_variable_mods
        ),
        (1, 4)
    );
    assert!(parameters
        .validate_ptm_library(&PtmLibrary::default())
        .is_ok());
}

#[test]
fn ptm_validation_rejects_conflicting_static_and_variable_definitions() {
    let empty = PtmLibrary::default();

    // Two different fixed definitions that can land on the same residue.
    let error = json_parameters(serde_json::json!({
        "static_mods": {
            "Carbamidomethyl": {"mass": 57.021464, "sites": ["C"]},
            "Other": {"mass": 58.0, "sites": ["first_residue:C"]}
        }
    }))
    .validate_ptm_library(&empty)
    .unwrap_err();
    assert!(
        error.starts_with("conflicting static modifications at"),
        "{error}"
    );

    // Distinct residues do not conflict.
    assert!(json_parameters(serde_json::json!({
        "static_mods": {
            "Carbamidomethyl": {"mass": 57.021464, "sites": ["C"]},
            "Other": {"mass": 58.0, "sites": ["K"]}
        }
    }))
    .validate_ptm_library(&empty)
    .is_ok());

    let error = json_parameters(serde_json::json!({
        "static_mods": {"Acetyl": {"mass": 42.010565, "sites": ["K"]}},
        "variable_mods": {"Acetyl": {"mass": 42.010565, "sites": ["S"]}}
    }))
    .validate_ptm_library(&empty)
    .unwrap_err();
    assert!(
        error.starts_with("modification `Acetyl` is defined in both static_mods and variable_mods"),
        "{error}"
    );

    let mut parameters = Builder::default().make_parameters();
    parameters.max_variable_mods = 3;
    parameters.max_total_variable_mods = 2;
    assert_eq!(
        parameters.validate_ptm_library(&empty).unwrap_err(),
        "database.max_total_variable_mods must be at least database.max_variable_mods"
    );
}

#[test]
fn ptm_validation_requires_limits_and_names_for_library_modifications() {
    use crate::ptm_library::PtmLibrarySite;
    let empty = PtmLibrary::default();
    let library = |modification: &str| {
        PtmLibrary::new(vec![PtmLibrarySite {
            attachment: Default::default(),
            protein: "P1".into(),
            position: 0,
            residue: b'K',
            modification: modification.into(),
        }])
    };

    // A configured PTM library requires every variable modification to be bounded.
    let mut parameters = json_parameters(serde_json::json!({
        "variable_mods": {"Oxidation": {"mass": 15.994915, "sites": ["M"]}}
    }));
    assert!(parameters.validate_ptm_library(&empty).is_ok());
    parameters.ptm_library = Some(PtmLibrarySettings {
        path: "sites.tsv".into(),
        strict: true,
    });
    assert!(parameters
        .validate_ptm_library(&empty)
        .unwrap_err()
        .starts_with("all variable modifications require `max_count` or `max_total_count`"));

    // Library-driven site modes need a name and a limit.
    let error = json_parameters(serde_json::json!({
        "variable_mods": {"Acetyl": {"mass": 42.010565, "sites": ["K"], "site_mode": "both"}}
    }))
    .validate_ptm_library(&empty)
    .unwrap_err();
    assert!(
        error.starts_with("variable modifications using `library` or `both` require `name`"),
        "{error}"
    );

    let bounded = |site_mode: &str| {
        json_parameters(serde_json::json!({
            "variable_mods": {
                "Acetyl": {"mass": 42.010565, "sites": ["K"], "max_count": 1, "site_mode": site_mode}
            }
        }))
    };
    assert!(bounded("library")
        .validate_ptm_library(&library("Acetyl"))
        .is_ok());
    assert!(bounded("both")
        .validate_ptm_library(&library("Acetyl"))
        .is_ok());
    assert_eq!(
        bounded("exhaustive")
            .validate_ptm_library(&library("Acetyl"))
            .unwrap_err(),
        "PTM library modification `Acetyl` must use site_mode `library` or `both`"
    );
    assert_eq!(
        bounded("both")
            .validate_ptm_library(&library("Methyl"))
            .unwrap_err(),
        "PTM library references undefined modification `Methyl`"
    );
}

#[test]
fn ptm_validation_rejects_inconsistent_named_database_modifications() {
    use crate::modification::{NeutralLossMode, VariableModification};
    let entry = |mass: f32, max_count: Option<usize>| {
        VarModEntry::Detailed(VariableModification {
            mass,
            max_count,
            max_total_count: None,
            name: Some("Acetyl".into()),
            neutral_losses: Vec::new(),
            neutral_loss_mode: NeutralLossMode::Optional,
            site_mode: SiteMode::Exhaustive,
            search_mode: SearchMode::Database,
            channel_offsets: Default::default(),
        })
    };
    let validate = |k: VarModEntry, s: VarModEntry| {
        Builder {
            variable_mods: Some(
                [("K".to_string(), vec![k]), ("S".to_string(), vec![s])]
                    .into_iter()
                    .collect(),
            ),
            ..Default::default()
        }
        .make_parameters()
        .validate_ptm_library(&PtmLibrary::default())
    };
    assert!(validate(entry(42.0, Some(1)), entry(42.0, Some(1))).is_ok());
    for (k, s) in [
        (entry(42.0, Some(1)), entry(43.0, Some(1))),
        (entry(42.0, Some(1)), entry(42.0, Some(2))),
        (entry(42.0, None), entry(42.0, Some(1))),
    ] {
        assert_eq!(
            validate(k, s).unwrap_err(),
            "variable modification `Acetyl` has inconsistent definitions across specificities"
        );
    }
}

#[test]
fn ptm_validation_rejects_mass_offsets_with_label_channels() {
    let error = json_parameters(serde_json::json!({
        "static_mods": {
            "SILAC-K": {"mass": 0.0, "channel_offsets": {"light": 0.0, "heavy": 8.014199}, "sites": ["K"]}
        },
        "variable_mods": {
            "Phospho": {"mass": 79.966331, "sites": ["S"], "search_mode": "mass_offset"}
        }
    }))
    .validate_ptm_library(&PtmLibrary::default())
    .unwrap_err();
    assert_eq!(
        error,
        "mass_offset modifications cannot be combined with channel-aware labels"
    );
}

#[test]
fn channel_validation_rejects_incomplete_or_redundant_channels() {
    let validate = |static_mods: serde_json::Value| {
        json_parameters(serde_json::json!({"static_mods": static_mods})).validate_channels()
    };
    assert!(validate(serde_json::json!({"C": 57.021464})).is_ok());
    assert!(validate(serde_json::json!({
        "SILAC-K": {"mass": 0.0, "channel_offsets": {"heavy": 8.014199}, "sites": ["K"]}
    }))
    .unwrap_err()
    .contains("at least two channels"));
    assert!(validate(serde_json::json!({
        "SILAC-K": {"mass": 0.0, "channel_offsets": {"light": 0.0, "heavy": 8.0}, "sites": ["K"]},
        "SILAC-R": {"mass": 0.0, "channel_offsets": {"light": 0.0, "medium": 6.0}, "sites": ["R"]}
    }))
    .unwrap_err()
    .contains("same channel names"));
    assert!(validate(serde_json::json!({
        "SILAC-K": {"mass": 0.0, "channel_offsets": {"light": 0.0, "heavy": 0.0}, "sites": ["K"]}
    }))
    .unwrap_err()
    .contains("at least one non-zero offset"));
    // `medium` and `heavy` shift K by the same amount, so they cannot be told apart.
    let error = validate(serde_json::json!({
        "SILAC-K": {
            "mass": 0.0,
            "channel_offsets": {"light": 0.0, "medium": 8.0, "heavy": 8.0},
            "sites": ["K"]
        }
    }))
    .unwrap_err();
    assert!(error.contains("is chemically identical"), "{error}");
}

#[test]
fn mass_offset_modifications_group_specificities_by_definition() {
    let parameters = json_parameters(serde_json::json!({
        "variable_mods": {
            "Phospho": {"mass": 79.966331, "sites": ["S", "T", "Y"], "search_mode": "mass_offset"},
            "Oxidation": {"mass": 15.994915, "sites": ["M"], "search_mode": "mass_offset"},
            "Acetyl": {"mass": 42.010565, "sites": ["K"]}
        },
        "generate_decoys": false, "peptide_min_mass": 0
    }));
    parameters
        .validate_ptm_library(&PtmLibrary::default())
        .unwrap();
    let offsets = parameters.mass_offset_modifications();
    assert_eq!(offsets.len(), 2);
    let group = |name: &str| {
        offsets
            .iter()
            .find(|offset| offset.definition.name.as_deref() == Some(name))
            .unwrap()
    };
    assert_eq!(
        group("Phospho").specificities,
        vec![
            ModificationSpecificity::Residue(b'S'),
            ModificationSpecificity::Residue(b'T'),
            ModificationSpecificity::Residue(b'Y'),
        ]
    );
    assert_eq!(
        group("Oxidation").specificities,
        vec![ModificationSpecificity::Residue(b'M')]
    );
    assert!(offsets
        .iter()
        .all(|offset| offset.site_mode == SiteMode::Exhaustive));

    // Only the database-mode modification is expanded into the index.
    let peptides = parameters.modify_digests(group_digests(vec![positional_digest(
        "MSKTY",
        Position::Internal,
    )]));
    let mut strings = peptides.iter().map(ToString::to_string).collect::<Vec<_>>();
    strings.sort();
    assert_eq!(strings, vec!["MSKTY", "MSK[Acetyl]TY"]);
}

#[test]
fn reversed_decoys_of_modified_targets_mirror_sites_and_keep_mass() {
    let parameters = json_parameters(serde_json::json!({
        "variable_mods": {"Oxidation": {"mass": 15.994915, "sites": ["M"]}},
        "generate_decoys": false, "peptide_min_mass": 0
    }));
    let targets = parameters.modify_digests(group_digests(vec![positional_digest(
        "PEPMK",
        Position::Full,
    )]));
    assert_eq!(targets.len(), 2);
    let peptides = parameters.add_reversed_decoys(targets);
    assert_eq!(peptides.len(), 4);
    let mut decoys = peptides
        .iter()
        .filter(|peptide| peptide.decoy)
        .map(ToString::to_string)
        .collect::<Vec<_>>();
    decoys.sort();
    assert_eq!(decoys, vec!["PMPEK", "PM[Oxidation]PEK"]);
    // Sorted by mass, each decoy sits beside its target with an identical mass.
    for pair in peptides.chunks(2) {
        assert_ne!(pair[0].decoy, pair[1].decoy);
        assert_eq!(pair[0].monoisotopic, pair[1].monoisotopic);
    }
    assert!(peptides
        .windows(2)
        .all(|pair| pair[0].monoisotopic <= pair[1].monoisotopic));

    // Already-decoy peptides pass through; re-adding decoys adds nothing.
    let again = parameters.add_reversed_decoys(peptides.clone());
    assert_eq!(again.len(), 4);

    // A reversal that reproduces another target's sequence is dropped, with
    // or without its modification (PEPM[Oxidation]K reverses to PM[Oxidation]PEK).
    let targets = parameters.modify_digests(group_digests(vec![
        positional_digest("PEPMK", Position::Full),
        positional_digest("PMPEK", Position::Full),
    ]));
    assert_eq!(targets.len(), 4);
    let peptides = parameters.add_reversed_decoys(targets);
    assert_eq!(peptides.len(), 4);
    assert!(peptides.iter().all(|peptide| !peptide.decoy));
}

#[test]
fn preflight_variant_counts_match_generated_peptides() {
    let fasta = Fasta::parse(">P1\nMCKMSTKPEPCMR\n".into(), "rev_", false).unwrap();
    let base = |modifications: serde_json::Value| {
        let mut value = serde_json::json!({
            "generate_decoys": false, "peptide_min_mass": 0, "peptide_max_mass": 100000,
            "enzyme": {"min_len": 1, "max_len": 50, "missed_cleavages": 1,
                       "cleave_at": "KR", "restrict": "P"}
        });
        for (key, entry) in modifications.as_object().unwrap() {
            value[key] = entry.clone();
        }
        value
    };
    let oxidation = serde_json::json!({"mass": 15.994915, "sites": ["M"]});
    let cases = [
        serde_json::json!({
            "static_mods": {"Carbamidomethyl": {"mass": 57.021464, "sites": ["C"]}},
            "variable_mods": {"Oxidation": oxidation},
            "max_variable_mods": 2
        }),
        serde_json::json!({
            "variable_mods": {
                "Oxidation": oxidation,
                "Acetyl": {"mass": 42.010565, "sites": ["protein_n_term", "peptide_n_term:K"]},
                "Phospho": {"mass": 79.966331, "sites": ["S", "T"], "max_count": 1}
            },
            "max_variable_mods": 3
        }),
        serde_json::json!({
            "variable_mods": {
                "Oxidation": oxidation,
                "Phospho": {"mass": 79.966331, "sites": ["S", "T"], "max_count": 1}
            },
            "max_variable_mods": 3,
            "max_combinations": 4
        }),
        serde_json::json!({
            "static_mods": {"TMT": {"mass": 229.162932, "sites": ["peptide_n_term", "K"]}},
            "variable_mods": {
                "Oxidation": oxidation,
                "Amidated": {"mass": -0.984016, "sites": ["peptide_c_term"]}
            },
            "max_variable_mods": 2
        }),
        serde_json::json!({
            "variable_mods": {
                "Oxidation": {"mass": 15.994915, "sites": ["M"], "max_count": 2},
                "Deamidated": {"mass": 0.984016, "sites": ["internal_residue:K"]}
            },
            "max_variable_mods": 1,
            "max_total_variable_mods": 3
        }),
    ];
    for case in cases {
        let parameters = json_parameters(base(case.clone()));
        parameters
            .validate_ptm_library(&PtmLibrary::default())
            .unwrap();
        let generated = parameters.digest(&fasta).len() as u64;
        let estimated = parameters.estimate_memory(&fasta).modified_peptides;
        assert_eq!(estimated, generated, "{case}");
    }
}

#[test]
fn paired_peptide_index_links_targets_and_generated_decoys() {
    let fasta = Fasta::parse(">P1\nMCKMSTKPEPCMR\n".into(), "rev_", false).unwrap();
    let build = |generate_decoys: bool, min_len: usize| {
        json_parameters(serde_json::json!({
            "variable_mods": {"Oxidation": {"mass": 15.994915, "sites": ["M"]}},
            "generate_decoys": generate_decoys, "peptide_min_mass": 0,
            "enzyme": {"min_len": min_len, "max_len": 50, "missed_cleavages": 0,
                       "cleave_at": "KR", "restrict": "P"}
        }))
        .build(fasta.clone())
    };

    // MCK reverses to itself, so neither of its forms receives a decoy or a pair.
    let database = build(true, 3);
    assert_eq!(database.peptides.len(), 10);
    let unpaired = (0..database.peptides.len())
        .filter(|&index| {
            database
                .paired_peptide_index(PeptideIx(index as u32))
                .is_none()
        })
        .map(|index| database.peptides[index].to_string())
        .collect::<HashSet<_>>();
    assert_eq!(
        unpaired,
        ["MCK", "M[Oxidation]CK"]
            .into_iter()
            .map(String::from)
            .collect()
    );

    // MSTKPEPCMR: four oxidation forms, each with a reversed decoy.
    let database = build(true, 4);
    let targets = database.peptides.iter().filter(|p| !p.decoy).count();
    assert_eq!(targets, 4);
    assert_eq!(database.peptides.len(), 8);
    for (index, peptide) in database.peptides.iter().enumerate() {
        let index = PeptideIx(index as u32);
        let paired = database
            .paired_peptide_index(index)
            .unwrap_or_else(|| panic!("{peptide} has no pair"));
        let partner = &database[paired];
        assert_ne!(partner.decoy, peptide.decoy);
        assert_eq!(partner.monoisotopic, peptide.monoisotopic);
        assert_eq!(partner.to_string(), peptide.reverse().to_string());
        assert_eq!(database.paired_peptide_index(paired), Some(index));
    }

    let database = build(false, 4);
    assert!(!database.peptides.is_empty());
    assert!((0..database.peptides.len()).all(|index| database
        .paired_peptide_index(PeptideIx(index as u32))
        .is_none()));
    assert_eq!(
        database.paired_peptide_index(PeptideIx(database.peptides.len() as u32 + 5)),
        None
    );
}

#[test]
fn same_peptidoform_compares_chemistry_decoy_state_and_sequence() {
    let parameters = json_parameters(serde_json::json!({
        "variable_mods": {"Oxidation": {"mass": 15.994915, "sites": ["M"]}},
        "generate_decoys": false, "peptide_min_mass": 0
    }));
    let variants = parameters.modify_digests(group_digests(vec![positional_digest(
        "PEPMK",
        Position::Full,
    )]));
    let unmodified = variants
        .iter()
        .find(|peptide| peptide.to_string() == "PEPMK")
        .unwrap();
    let oxidized = variants
        .iter()
        .find(|peptide| peptide.to_string() == "PEPM[Oxidation]K")
        .unwrap();
    assert!(same_peptidoform(oxidized, &oxidized.clone()));
    assert!(!same_peptidoform(unmodified, oxidized));
    let mut decoy = oxidized.clone();
    decoy.decoy = true;
    assert!(!same_peptidoform(oxidized, &decoy));
    // Reversing twice restores the same chemical peptidoform.
    assert!(same_peptidoform(oxidized, &oxidized.reverse().reverse()));
}

#[test]
fn clipped_initiator_methionine_peptides_take_protein_n_terminal_mods() {
    let fasta = Fasta::parse(">P1\nMASPEPTIDEAAKGGLLR\n".into(), "rev_", true).unwrap();
    let parameters = |clip: Option<bool>| {
        let mut config = serde_json::json!({
            "enzyme": {"missed_cleavages": 0, "min_len": 5},
            "peptide_min_mass": 100.0,
            "variable_mods": {"[": [42.010565]},
            "generate_decoys": true,
        });
        if let Some(clip) = clip {
            config["clip_n_term_met"] = clip.into();
        }
        serde_json::from_value::<Builder>(config)
            .unwrap()
            .make_parameters()
    };

    let defaults = parameters(None);
    assert!(defaults.clip_n_term_met);
    let peptides = defaults.digest(&fasta);
    let describe = |peptide: &Peptide| (peptide.decoy, peptide.to_string());
    let names = peptides.iter().map(describe).collect::<HashSet<_>>();
    for expected in [
        (false, "MASPEPTIDEAAK"),
        (false, "[+42.010567]-MASPEPTIDEAAK"),
        (false, "ASPEPTIDEAAK"),
        (false, "[+42.010567]-ASPEPTIDEAAK"),
        // The generated decoy of a clipped peptide is N-terminal too.
        (true, "AAAEDITPEPSK"),
        (true, "[+42.010567]-AAAEDITPEPSK"),
    ] {
        assert!(
            names.contains(&(expected.0, expected.1.to_string())),
            "missing {expected:?} in {names:?}"
        );
    }
    let clipped = peptides
        .iter()
        .find(|peptide| peptide.to_string() == "ASPEPTIDEAAK")
        .unwrap();
    assert_eq!(clipped.position, Position::Nterm);
    assert!(!clipped.semi_enzymatic);
    assert_eq!(clipped.protein_sites[0].start, Some(1));

    let unclipped = parameters(Some(false)).digest(&fasta);
    assert_eq!(unclipped.len() + 4, peptides.len());
    assert!(unclipped
        .iter()
        .all(|peptide| !peptide.sequence.starts_with("ASPEPT")));
}

#[test]
fn clipped_initiator_methionine_peptides_match_protein_n_terminal_motifs() {
    let fasta = Fasta::parse(">P1\nMASPEPTIDEAAKGGLLR\n".into(), "rev_", true).unwrap();
    let parameters = serde_json::from_value::<Builder>(serde_json::json!({
        "enzyme": {"missed_cleavages": 0, "min_len": 5},
        "peptide_min_mass": 100.0,
        "variable_mods": {"Nterm-A": {"mass": 10.0, "sites": ["motif:<A*"]}},
        "generate_decoys": false,
    }))
    .unwrap()
    .make_parameters();
    let names = parameters
        .digest(&fasta)
        .iter()
        .map(|peptide| peptide.to_string())
        .collect::<HashSet<_>>();
    assert!(names.contains("A[Nterm-A]SPEPTIDEAAK"), "{names:?}");
    // The anchor still needs the protein N-terminus: an internal A does not
    // qualify.
    assert!(names.iter().all(|name| !name.contains("A[Nterm-A]A")));
}

#[test]
fn fasta_copy_of_a_peptide_keeps_its_position_over_a_tsv_copy() {
    let fasta = Fasta::parse(">P1\nGGGKPEPTIDERGGGK\n".into(), "rev_", false).unwrap();
    let parameters = serde_json::from_value::<Builder>(serde_json::json!({
        "enzyme": {"missed_cleavages": 0, "min_len": 5},
        "generate_decoys": false,
    }))
    .unwrap()
    .make_parameters();
    let mut peptides = parameters.digest(&fasta);
    // The TSV copy is a whole-protein, fully enzymatic placeholder.
    peptides.extend(parameters.peptides_from_tsv("sequence\tprotein\nPEPTIDER\tT1\n"));
    Parameters::reorder_peptides(&mut peptides);
    let merged = peptides
        .iter()
        .filter(|peptide| peptide.to_string() == "PEPTIDER")
        .collect::<Vec<_>>();
    assert_eq!(merged.len(), 1);
    assert_eq!(merged[0].position, Position::Internal);
    assert_eq!(merged[0].protein_sites[0].start, Some(4));
    let proteins = merged[0]
        .proteins
        .iter()
        .map(|protein| protein.to_string())
        .collect::<Vec<_>>();
    assert_eq!(proteins, ["P1", "T1"]);
}

#[test]
fn generated_decoys_of_n_terminal_peptides_stay_balanced_and_n_terminal() {
    // ASEQK is semi-enzymatic at the N-terminus of "a" (K-P is not cut) and
    // the fully enzymatic clipped N-terminal peptide of "b".
    let fasta = Fasta::parse(
        ">a\nASEQKPLLRGGDDK\n>b\nMASEQKGLLRWWEEK\n".into(),
        "rev_",
        true,
    )
    .unwrap();
    let parameters = serde_json::from_value::<Builder>(serde_json::json!({
        "enzyme": {"missed_cleavages": 1, "min_len": 5, "semi_enzymatic": true, "restrict": "P"},
        "peptide_min_mass": 100.0,
        "generate_decoys": true,
    }))
    .unwrap()
    .make_parameters();
    let peptides = parameters.digest(&fasta);
    let n_terminal = |decoy: bool| {
        peptides
            .iter()
            .filter(|peptide| peptide.decoy == decoy && peptide.position == Position::Nterm)
            .map(|peptide| (peptide.semi_enzymatic, peptide.missed_cleavages))
            .collect::<Vec<_>>()
    };
    let (mut targets, mut decoys) = (n_terminal(false), n_terminal(true));
    targets.sort_unstable();
    decoys.sort_unstable();
    assert!(!targets.is_empty());
    assert_eq!(targets, decoys);
    let find = |name: &str| {
        peptides
            .iter()
            .find(|peptide| peptide.to_string() == name)
            .unwrap_or_else(|| panic!("missing {name}"))
    };
    for name in ["ASEQK", "AQESK"] {
        let peptide = find(name);
        assert_eq!(peptide.position, Position::Nterm, "{name}");
        assert!(!peptide.semi_enzymatic, "{name}");
        assert_eq!(peptide.protein_sites.len(), 2, "{name}");
    }
    assert!(find("AQESK").decoy);
    // The unclipped N-terminal peptide of "b" and its decoy.
    assert_eq!(find("MASEQK").position, Position::Nterm);
    assert_eq!(find("MQESAK").position, Position::Nterm);
    assert!(find("MQESAK").decoy);
}

fn ambiguous_parameters(expand: bool) -> Parameters {
    Builder {
        enzyme: Some(EnzymeBuilder {
            missed_cleavages: Some(0),
            min_len: Some(5),
            max_len: Some(50),
            ..Default::default()
        }),
        expand_ambiguous_residues: Some(expand),
        ..Default::default()
    }
    .make_parameters()
}

#[test]
fn expanded_ambiguous_peptides_keep_proteins_and_database_sequence() {
    let fasta = || Fasta::parse(">P1\nMRGEPXIDEK\n>P2\nMRGEPTIDEK\n".into(), "rev_", true).unwrap();

    let off = ambiguous_parameters(false);
    let on = ambiguous_parameters(true);
    assert!(
        on.estimate_memory(&fasta()).unmodified_peptides
            > off.estimate_memory(&fasta()).unmodified_peptides
    );

    let database = off.build(fasta());
    let peptide = database
        .peptides
        .iter()
        .find(|peptide| &peptide.sequence[..] == b"GEPTIDEK")
        .unwrap();
    assert_eq!(peptide.proteins("rev_", true), "P2");
    assert_eq!(peptide.database_peptide(), None);

    let database = on.build(fasta());
    let find = |sequence: &[u8]| {
        database
            .peptides
            .iter()
            .find(|peptide| &peptide.sequence[..] == sequence)
            .unwrap_or_else(|| panic!("{}", String::from_utf8_lossy(sequence)))
    };
    let target = find(b"GEPTIDEK");
    assert_eq!(target.proteins("rev_", true), "P1;P2");
    assert_eq!(target.database_peptide().as_deref(), Some("GEPXIDEK"));
    let variant = find(b"GEPWIDEK");
    assert_eq!(variant.proteins("rev_", true), "P1");
    assert_eq!(variant.database_peptide().as_deref(), Some("GEPXIDEK"));
    let decoy = find(b"GEDITPEK");
    assert!(decoy.decoy);
    assert_eq!(decoy.database_peptide().as_deref(), Some("GEDIXPEK"));
    assert!(database
        .peptides
        .iter()
        .all(|peptide| !crate::ambiguous_residues::is_ambiguous(&peptide.sequence)));
}

#[test]
fn ambiguous_proteins_are_reported_as_expanded_or_dropped() {
    let fasta = Fasta::parse(
        ">P1\nMRGEPXIDEK\n>P2\nMRGEPJIDEK\n>P3\nMRBZK\n".into(),
        "rev_",
        true,
    )
    .unwrap();
    assert_eq!(fasta.ambiguous_protein_count(), 2);
    let (level, message) = ambiguous_parameters(true)
        .ambiguous_proteins_message(&fasta)
        .unwrap();
    assert_eq!(level, log::Level::Info);
    assert!(message.starts_with("2 FASTA protein(s)"), "{message}");
    assert!(
        message.contains("expanded into up to 20 variant(s)"),
        "{message}"
    );
    let (level, message) = ambiguous_parameters(false)
        .ambiguous_proteins_message(&fasta)
        .unwrap();
    assert_eq!(level, log::Level::Warn);
    assert!(
        message.contains("not searched unless database.expand_ambiguous_residues"),
        "{message}"
    );
    // J alone carries the I/L mass and needs neither.
    let fasta = Fasta::parse(">P2\nMRGEPJIDEK\n".into(), "rev_", true).unwrap();
    assert_eq!(
        ambiguous_parameters(false).ambiguous_proteins_message(&fasta),
        None
    );
}

#[test]
fn expanded_peptides_report_substitutions() {
    let fasta = Fasta::parse(
        ">P1\nMRGEPXIDEK\n>P2\nMRGEPTIDEK\n>P3\nMRABEZAK\n>P4\nMRSAMPLEK\n".into(),
        "rev_",
        true,
    )
    .unwrap();
    let database = ambiguous_parameters(true).build(fasta);
    let find = |sequence: &[u8]| {
        database
            .peptides
            .iter()
            .find(|peptide| &peptide.sequence[..] == sequence)
            .unwrap_or_else(|| panic!("{}", String::from_utf8_lossy(sequence)))
    };
    // Written plainly in P2, so nothing was substituted.
    assert_eq!(find(b"GEPTIDEK").substitutions(), "");
    assert_eq!(find(b"GEPWIDEK").substitutions(), "X4W");
    // The generated decoy of GEPWIDEK, at decoy positions.
    let decoy = find(b"GEDIWPEK");
    assert!(decoy.decoy);
    assert_eq!(decoy.substitutions(), "X5W");
    assert_eq!(find(b"ADEQAK").substitutions(), "B2D;Z4Q");
    assert_eq!(find(b"ANEEAK").substitutions(), "B2N;Z4E");
    // No ambiguous residues.
    assert_eq!(find(b"SAMPLEK").substitutions(), "");
}

/// Merged copies from several proteins: the substitutions come from the first
/// expanded occurrence in protein order, and any occurrence with the residues
/// as written leaves them empty.
#[test]
fn substitutions_follow_the_first_expanded_occurrence() {
    let parameters = |merge: bool| {
        serde_json::from_value::<Builder>(serde_json::json!({
            "enzyme": {"missed_cleavages": 0, "min_len": 5},
            "generate_decoys": false,
            "expand_ambiguous_residues": true,
            "merge_isoleucine_leucine": merge,
        }))
        .unwrap()
        .make_parameters()
    };
    let peptide = |fasta: &str, merge: bool| {
        let fasta = Fasta::parse(fasta.into(), "rev_", false).unwrap();
        parameters(merge)
            .digest(&fasta)
            .into_iter()
            .find(|peptide| &peptide.sequence[..] == b"AEPTIDEK")
            .expect("AEPTIDEK")
    };

    // Both proteins are expanded, at different residues: the first in
    // protein order, A, gives them, whatever the FASTA order.
    for fasta in [
        ">A\nGGKAEPXIDEK\n>B\nGGKAEPTXDEK\n",
        ">B\nGGKAEPTXDEK\n>A\nGGKAEPXIDEK\n",
    ] {
        for merge in [false, true] {
            let merged = peptide(fasta, merge);
            assert_eq!(merged.proteins("rev_", false), "A;B");
            assert_eq!(merged.substitutions(), "X4T");
        }
    }

    // A real residue in any protein wins, before or after the expanded one.
    for fasta in [
        ">A\nGGKAEPXIDEK\n>B\nGGKAEPTIDEK\n",
        ">A\nGGKAEPTIDEK\n>B\nGGKAEPXIDEK\n",
    ] {
        let merged = peptide(fasta, true);
        assert_eq!(merged.proteins("rev_", false), "A;B");
        assert_eq!(merged.substitutions(), "");
    }
    // Also a merged I/L/J twin; unmerged, the twin is another peptide.
    let fasta = ">A\nGGKAEPXIDEK\n>B\nGGKAEPTJDEK\n";
    let merged = peptide(fasta, true);
    assert_eq!(merged.proteins("rev_", false), "A;B");
    assert_eq!(merged.substitutions(), "");
    let unmerged = peptide(fasta, false);
    assert_eq!(unmerged.proteins("rev_", false), "A");
    assert_eq!(unmerged.substitutions(), "X4T");

    // A peptide TSV row lists the sequence as written.
    let fasta = Fasta::parse(">A\nGGKAEPXIDEK\n".into(), "rev_", false).unwrap();
    let parameters = parameters(true);
    let mut peptides = parameters.digest(&fasta);
    peptides.extend(parameters.peptides_from_tsv("sequence\tprotein\nAEPTIDEK\tT1\n"));
    parameters.reorder_merged_peptides(&mut peptides);
    let merged = peptides
        .iter()
        .find(|peptide| &peptide.sequence[..] == b"AEPTIDEK")
        .unwrap();
    assert_eq!(merged.proteins("rev_", false), "A;T1");
    assert_eq!(merged.substitutions(), "");
}

/// Small deterministic generator so the randomized comparisons are repeatable.
struct Lcg(u64);

impl Lcg {
    fn next(&mut self) -> u64 {
        self.0 = self
            .0
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        self.0 >> 33
    }

    fn below(&mut self, n: u64) -> u64 {
        self.next() % n.max(1)
    }
}

#[test]
fn interleaved_bucket_ranges_match_bucket_search_on_random_queries() {
    let mut rng = Lcg(0x5eed);
    for _ in 0..200 {
        // Buckets of random length (including empty and single-fragment
        // ones), each sorted by peptide index with duplicates.
        let bucket_count = 1 + rng.below(40) as usize;
        let mut buckets = Vec::with_capacity(bucket_count);
        let mut fragments = Vec::new();
        let max_id = 1 + rng.below(3000) as u32;
        for bucket in 0..bucket_count {
            let len = match rng.below(4) {
                0 => rng.below(3),
                1 => rng.below(64),
                _ => rng.below(5000),
            } as usize;
            let mut ids = (0..len)
                .map(|_| rng.below(u64::from(max_id)) as u32)
                .collect::<Vec<_>>();
            ids.sort_unstable();
            let start = fragments.len() as u32;
            fragments.extend(ids.into_iter().map(|peptide_index| PackedFragment {
                peptide_index,
                mass_suffix: rng.below(1 << FRAGMENT_MASS_SUFFIX_BITS) as u16,
            }));
            buckets.push(FragmentBucket {
                mass_prefix: (bucket as u32 + 0x4000) << FRAGMENT_MASS_SUFFIX_BITS,
                start,
                end: fragments.len() as u32,
            });
        }
        let index = FragmentIndex { buckets, fragments };

        let mut out = Vec::new();
        for _ in 0..20 {
            let a = rng.below(u64::from(max_id) + 2) as u32;
            let b = rng.below(u64::from(max_id) + 2) as u32;
            let (lo, hi) = match rng.below(8) {
                0 => (0, u32::MAX),
                1 => (a, a),
                _ => (a.min(b), a.max(b)),
            };
            // Random subset with repeats, in random order.
            let queried = (0..rng.below(50))
                .map(|_| rng.below(bucket_count as u64) as u32)
                .collect::<Vec<_>>();
            index.resolve_bucket_ranges(&queried, lo, hi, &mut out);
            assert_eq!(out.len(), queried.len());
            for (&bucket, &(first, last)) in queried.iter().zip(&out) {
                let expected = index
                    .bucket_search(bucket as usize, lo, hi)
                    .collect::<Vec<_>>();
                let actual = index
                    .range_iter(bucket as usize, first, last)
                    .collect::<Vec<_>>();
                assert_eq!(actual, expected, "bucket={bucket} lo={lo} hi={hi}");
                let span = index.buckets[bucket as usize];
                assert!(span.start <= first && first <= last && last <= span.end);
            }
        }
    }
}

#[test]
fn batched_page_search_matches_per_peak_page_search() {
    const RESIDUES: &[u8] = b"ACDEFGHIKLMNPQRSTVWY";
    let mut rng = Lcg(0xfeed);
    let mut text = String::new();
    for protein in 0..40 {
        text.push_str(&format!(">P{protein}\n"));
        for _ in 0..(50 + rng.below(250)) {
            text.push(RESIDUES[rng.below(RESIDUES.len() as u64) as usize] as char);
        }
        text.push('\n');
    }
    let fasta = Fasta::parse(text, "rev_", true).unwrap();
    let mut matched = 0;
    for bucket_size in [1usize, 8, 64, 8192] {
        let parameters = Builder {
            bucket_size: Some(bucket_size),
            ..Builder::default()
        }
        .make_parameters();
        let database = parameters.clone().build(fasta.clone());
        assert!(!database.peptides.is_empty());

        for _ in 0..300 {
            let peptide = &database.peptides[rng.below(database.peptides.len() as u64) as usize];
            let width = [0.01, 0.5, 5.0, 200.0][rng.below(4) as usize];
            let precursor_tol = Tolerance::Da(-width, width);
            let fragment_tol = match rng.below(2) {
                0 => Tolerance::Ppm(-20.0, 20.0),
                _ => Tolerance::Da(-0.05, 0.05),
            };
            let query = database.query(peptide.monoisotopic, precursor_tol, fragment_tol);

            // Theoretical fragments of the chosen peptide plus noise peaks,
            // sorted by mass as the scorer supplies them.
            let mut masses = preliminary_fragment_masses(&parameters, peptide)
                .map(|mass| mass + (rng.below(200) as f32 - 100.0) * 1e-4)
                .collect::<Vec<_>>();
            masses.extend((0..rng.below(80)).map(|_| 100.0 + rng.below(300_000) as f32 * 0.01));
            masses.sort_unstable_by(f32::total_cmp);
            let shifts: &[f32] = match rng.below(3) {
                0 => &[],
                1 => &[15.994915],
                _ => &[79.96633, -18.010565],
            };

            let mut expected = Vec::new();
            for &mass in &masses {
                expected.extend(query.page_search(mass));
                for &shift in shifts {
                    expected.extend(query.page_search_shifted(mass, shift));
                }
            }
            let mut actual = Vec::new();
            query.page_search_batch(masses.iter().copied(), shifts, |frag| actual.push(frag));
            assert_eq!(actual, expected, "bucket_size={bucket_size}");
            matched += actual.len();
        }
    }
    assert!(matched > 10_000, "only {matched} matches compared");
}

fn isoleucine_leucine_parameters(merge: bool, extra: serde_json::Value) -> Parameters {
    let mut config = serde_json::json!({
        "enzyme": {"missed_cleavages": 0, "min_len": 5, "semi_enzymatic": false},
        "peptide_min_mass": 100.0,
        "generate_decoys": true,
        "merge_isoleucine_leucine": merge,
    });
    for (key, value) in extra.as_object().unwrap() {
        config[key] = value.clone();
    }
    serde_json::from_value::<Builder>(config)
        .unwrap()
        .make_parameters()
}

/// Unmodified peptides equal to APEPLDEK once I, L and J are one residue.
fn isoleucine_leucine_twins(peptides: &[Peptide], decoy: bool) -> Vec<&Peptide> {
    let canonical: &[u8] = if decoy { b"AEDLPEPK" } else { b"APEPLDEK" };
    peptides
        .iter()
        .filter(|peptide| {
            peptide.decoy == decoy
                && peptide.modifications.is_empty()
                && crate::ambiguous_residues::isoleucine_leucine_eq(&peptide.sequence, canonical)
        })
        .collect()
}

#[test]
fn isoleucine_leucine_twins_merge_into_one_peptide() {
    // APEPIDEK in "a" follows a G: with semi-enzymatic digestion it is
    // semi-enzymatic there, and fully enzymatic as APEPLDEK in "b".
    let fasta = Fasta::parse(
        ">c\nSSRAPEPJDEKWWK\n>b\nLLRAPEPLDEKGGWR\n>a\nGGGAPEPIDEKR\n".into(),
        "rev_",
        true,
    )
    .unwrap();
    let semi = serde_json::json!({
        "enzyme": {"missed_cleavages": 0, "min_len": 5, "semi_enzymatic": true}
    });
    let peptides = isoleucine_leucine_parameters(true, semi.clone()).digest(&fasta);
    let targets = isoleucine_leucine_twins(&peptides, false);
    let decoys = isoleucine_leucine_twins(&peptides, true);
    assert_eq!((targets.len(), decoys.len()), (1, 1), "{targets:?}");
    // The displayed sequence comes from the first protein, "a"; the rest
    // is the most enzymatic occurrence, in "b".
    let target = targets[0];
    assert_eq!(target.to_string(), "APEPIDEK");
    assert!(!target.semi_enzymatic);
    assert_eq!(
        target.proteins.as_slice(),
        &["a".into(), "b".into(), "c".into()]
    );
    assert_eq!(target.protein_sites.len(), 3);
    assert_eq!(decoys[0].to_string(), "AEDIPEPK");
    assert_eq!(decoys[0].proteins, target.proteins);
    // The decoy is the target's pair in the built database.
    let fasta_again = Fasta::parse(
        ">c\nSSRAPEPJDEKWWK\n>b\nLLRAPEPLDEKGGWR\n>a\nGGGAPEPIDEKR\n".into(),
        "rev_",
        true,
    )
    .unwrap();
    let database = isoleucine_leucine_parameters(true, semi.clone()).build(fasta_again);
    let index = database
        .peptides
        .iter()
        .position(|peptide| peptide.to_string() == "APEPIDEK" && !peptide.decoy)
        .unwrap();
    let pair = database
        .paired_peptide_index(PeptideIx(index as u32))
        .expect("paired decoy");
    assert_eq!(database.peptides[pair.0 as usize].to_string(), "AEDIPEPK");

    // Without merging, each twin is its own target with its own decoy.
    let peptides = isoleucine_leucine_parameters(false, semi).digest(&fasta);
    assert_eq!(isoleucine_leucine_twins(&peptides, false).len(), 3);
    assert_eq!(isoleucine_leucine_twins(&peptides, true).len(), 3);
}

#[test]
fn twins_from_one_protein_list_it_once() {
    let fasta = Fasta::parse(">a\nGGRAPEPIDEKAPEPLDEKR\n".into(), "rev_", true).unwrap();
    let peptides = isoleucine_leucine_parameters(true, serde_json::json!({})).digest(&fasta);
    let targets = isoleucine_leucine_twins(&peptides, false);
    assert_eq!(targets.len(), 1, "{targets:?}");
    assert_eq!(targets[0].proteins.as_slice(), &["a".into()]);
    assert_eq!(targets[0].protein_sites.len(), 2);
}

#[test]
fn modifications_on_isoleucine_or_leucine_keep_twins_apart() {
    let fasta = Fasta::parse(
        ">a\nGGRAPEPIDEKR\n>b\nLLRAPEPLDEKGGWR\n>c\nSSRAPEPJDEKWWK\n".into(),
        "rev_",
        true,
    )
    .unwrap();
    let names = |peptides: &[Peptide]| {
        let mut names = peptides
            .iter()
            .filter(|peptide| {
                !peptide.decoy
                    && peptide.sequence.len() == 8
                    && crate::ambiguous_residues::isoleucine_leucine_eq(
                        &peptide.sequence,
                        b"APEPLDEK",
                    )
            })
            .map(|peptide| (peptide.to_string(), peptide.proteins.len()))
            .collect::<Vec<_>>();
        names.sort();
        names
    };

    // A static modification on I: the I twin is heavier, L and J merge.
    let peptides =
        isoleucine_leucine_parameters(true, serde_json::json!({"static_mods": {"I": 10.0}}))
            .digest(&fasta);
    assert_eq!(
        names(&peptides),
        [
            ("APEPI[+10]DEK".to_string(), 1),
            ("APEPLDEK".to_string(), 2)
        ]
    );

    // A variable modification on L: the modified L twin stays apart, the
    // unmodified forms merge; J never carries it.
    let peptides =
        isoleucine_leucine_parameters(true, serde_json::json!({"variable_mods": {"L": [15.9949]}}))
            .digest(&fasta);
    assert_eq!(
        names(&peptides),
        [
            ("APEPIDEK".to_string(), 3),
            ("APEPL[+15.9949]DEK".to_string(), 1)
        ]
    );
}

#[test]
fn generated_decoys_that_are_twins_of_a_target_are_dropped() {
    // The decoy of DLGEENFK is DFNEEGLK, a twin of the target DFNEEGIK.
    let fasta = Fasta::parse(">c\nDLGEENFKR\n>d\nMSSRDFNEEGIKR\n".into(), "rev_", true).unwrap();
    for merge in [false, true] {
        let peptides = isoleucine_leucine_parameters(merge, serde_json::json!({})).digest(&fasta);
        let decoy = peptides
            .iter()
            .any(|peptide| peptide.decoy && peptide.to_string() == "DFNEEGLK");
        assert_eq!(decoy, !merge);
    }
}

/// Peptide TSV rows with B, X or Z are expanded like FASTA digests when
/// expansion is on, report their substitutions, and are otherwise skipped.
#[test]
fn ambiguous_peptide_tsv_rows_expand_with_substitutions() {
    let tsv = "sequence\tprotein\nPEPXIDEK\tT1\nPEBTIDZK\tT2\nSAMPLEK\tT3\n";
    let peptides = ambiguous_parameters(true).peptides_from_tsv(tsv);
    let find = |sequence: &[u8]| {
        peptides
            .iter()
            .find(|peptide| &peptide.sequence[..] == sequence)
            .unwrap_or_else(|| panic!("{}", String::from_utf8_lossy(sequence)))
    };
    let targets = |protein: &str| {
        peptides
            .iter()
            .filter(|peptide| !peptide.decoy && peptide.proteins("rev_", false) == protein)
            .count()
    };
    // 20 variants, with the I and L twins merged.
    assert_eq!(targets("T1"), 19);
    assert_eq!(targets("T2"), 4);
    let expanded = find(b"PEPWIDEK");
    assert_eq!(expanded.substitutions(), "X4W");
    assert_eq!(expanded.proteins("rev_", false), "T1");
    // No coordinates are invented: the row is still unplaced, like a plain one.
    let enzyme = ambiguous_parameters(true)
        .enzyme_parameters()
        .enzyme
        .unwrap();
    assert_eq!(
        crate::digestion::classify_peptide(&enzyme, expanded, false),
        crate::digestion::classify_peptide(&enzyme, find(b"SAMPLEK"), false),
    );
    // The generated decoy of PEPWIDEK, at decoy positions.
    let decoy = find(b"PEDIWPEK");
    assert!(decoy.decoy);
    assert_eq!(decoy.substitutions(), "X5W");
    assert_eq!(find(b"PENTIDQK").substitutions(), "B3N;Z7Q");
    assert_eq!(find(b"SAMPLEK").substitutions(), "");

    // Expansion off: skipped with a warning, as before.
    let peptides = ambiguous_parameters(false).peptides_from_tsv(tsv);
    assert!(peptides
        .iter()
        .all(|peptide| peptide.proteins("rev_", false) == "T3"));
    assert!(!peptides.is_empty());

    // Over the variant cap: the X row is dropped, the B/Z row (4) kept.
    let mut capped = ambiguous_parameters(true);
    capped.max_ambiguous_variants = 4;
    let peptides = capped.peptides_from_tsv(tsv);
    assert!(peptides
        .iter()
        .all(|peptide| peptide.proteins("rev_", false) != "T1"));
    assert_eq!(
        peptides
            .iter()
            .filter(|peptide| !peptide.decoy && peptide.proteins("rev_", false) == "T2")
            .count(),
        4
    );
}
