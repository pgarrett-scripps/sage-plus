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
        }],
        reference,
    }
}

#[test]
fn sequence_coherent_partition_never_splits_terminal_variants() {
    let chunks = Parameters::partition_digests_by_sequence(
        vec![
            digest_group("PEPTIDER", Position::Internal),
            digest_group("SEQUENCEK", Position::Full),
            digest_group("PEPTIDER", Position::Nterm),
        ],
        1,
    );

    assert_eq!(chunks.len(), 2);
    assert!(chunks.iter().any(|chunk| {
        chunk.len() == 2
            && chunk
                .iter()
                .all(|group| group.reference.sequence == "PEPTIDER")
    }));
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
        fasta: "none".into(),
        peptides: None,
        custom_cleavage_sites: None,
        prefilter: false,
        prefilter_chunk_size: 0,
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
        assert!(sequences.iter().any(|(_, decoy)| *decoy));
        for peptide in &peptides {
            assert!(
                peptide
                    .sequence
                    .iter()
                    .all(|residue| crate::mass::VALID_AA.contains(residue)),
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
        .all(|residue| crate::mass::VALID_AA.contains(residue))));
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
