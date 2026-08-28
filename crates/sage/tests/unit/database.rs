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
fn fragment_index_preserves_exact_ids_masses_and_ranges() {
    assert_eq!(std::mem::size_of::<PackedFragment>(), 6);
    let parameters = Builder {
        generate_decoys: Some(false),
        ..Builder::default()
    }
    .make_parameters();
    let mut peptides = parameters
        .peptides_from_tsv("sequence\tprotein\nPEPTIDER\tprotein-a\nSEQUENCEK\tprotein-b\n");
    Parameters::reorder_peptides(&mut peptides);

    let mut expected = BTreeMap::<u32, Vec<Theoretical>>::new();
    let suffix_bits = fragment_suffix_bits(parameters.bucket_size);
    for (peptide_index, peptide) in peptides.iter().enumerate() {
        for mass in preliminary_fragment_masses(&parameters, peptide) {
            expected
                .entry(mass.to_bits() >> suffix_bits)
                .or_default()
                .push(Theoretical {
                    peptide_index: PeptideIx(peptide_index as u32),
                    fragment_mz: mass,
                });
        }
    }

    let database = parameters.build_from_peptides(peptides);
    let expected = expected.into_values().flatten().collect::<Vec<_>>();
    let actual = (0..database.buckets().len())
        .flat_map(|bucket| database.fragments.bucket(bucket))
        .collect::<Vec<_>>();
    assert_eq!(actual, expected);
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
            "C": {
                "mass": 57.0215,
                "name": "Carbamidomethyl"
            }
        },
        "variable_mods": {
            "M": [15.9949],
            "K": [
                {
                    "mass": 42.0106,
                    "max_count": 1,
                    "name": "Acetyl",
                    "neutral_losses": [17.0265],
                    "neutral_loss_mode": "required"
                },
                {"mass": 14.0157}
            ]
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
    let k_entries = &serialized["variable_mods"]["K"];
    assert!(k_entries[0].is_object());
    assert_eq!(k_entries[0]["max_count"], 1);
    assert_eq!(k_entries[0]["name"], "Acetyl");
    assert_eq!(k_entries[0]["neutral_loss_mode"], "required");
    assert!(k_entries[1].is_object());
    assert!(k_entries[1].get("max_count").is_none());
    assert!(serialized["variable_mods"]["M"][0].is_number());
    assert_eq!(serialized["static_mods"]["C"]["name"], "Carbamidomethyl");
}

#[test]
fn channel_offsets_generate_complete_static_channels() {
    let builder: Builder = serde_json::from_value(serde_json::json!({
        "fasta": "none",
        "generate_decoys": false,
        "static_mods": {
            "K": {
                "mass": 0.0,
                "name": "SILAC-K",
                "channel_offsets": {"light": 0.0, "heavy": 8.014199}
            },
            "R": {
                "mass": 0.0,
                "name": "SILAC-R",
                "neutral_losses": [17.026549],
                "neutral_loss_mode": "required",
                "channel_offsets": {"light": 0.0, "heavy": 10.008269}
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
            "K": {
                "mass": 0.0,
                "channel_offsets": {"light": 0.0, "heavy": 8.014199}
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
        "variable_mods": {"K": [{
            "mass": 0.0,
            "name": "SILAC-K",
            "channel_offsets": {"light": 0.0, "heavy": 8.014199}
        }]}
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
            "K": {
                "mass": 229.162932,
                "name": "TMT-SILAC-K",
                "channel_offsets": {"light": 0.0, "heavy": 8.014199}
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
            "K": {
                "mass": 229.16293,
                "name": "Labeled-K",
                "neutral_losses": [17.02655],
                "neutral_loss_mode": "required",
                "channel_offsets": {"light": 0.0, "heavy": 8.014199}
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
            "S": [{
                "mass": 79.96633,
                "name": "Phospho",
                "max_count": 2,
                "site_mode": "both",
                "neutral_losses": [97.9769]
            }]
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
    assert_eq!(serialized["variable_mods"]["S"][0]["site_mode"], "both");
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
                mass: 79.96633,
                max_count: Some(2),
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
            protein: Arc::from("P1"),
            position: 1,
            residue: b'S',
            modification: Arc::from("Phospho"),
        },
        PtmLibrarySite {
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
                mass: 79.96633,
                max_count: Some(2),
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
            protein: Arc::from("P1"),
            position: 1,
            residue: b'S',
            modification: Arc::from("Phospho"),
        },
        PtmLibrarySite {
            protein: Arc::from("P2"),
            position: 2,
            residue: b'S',
            modification: Arc::from("Phospho"),
        },
    ])));
    let fasta = Fasta::parse(">P1\nMSSK\n>P2\nMSSK\n".into(), "rev_", false).unwrap();

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
