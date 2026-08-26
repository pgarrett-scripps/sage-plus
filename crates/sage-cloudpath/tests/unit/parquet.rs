use super::*;
use parquet::file::reader::{FileReader, SerializedFileReader};
use sage_core::database::PeptideIx;
use sage_core::enzyme::ProteinOccurrence;
use sage_core::ion_series::Kind;
use sage_core::peptide::Peptide;
use sage_core::spectral_library::{LibraryFragment, SpectralLibraryEntry, SpectralLibrarySettings};
use std::sync::Arc;

#[test]
fn lfq_preserves_missingness_and_ms2_evidence() -> parquet::errors::Result<()> {
    let mut database = IndexedDatabase::default();
    database.peptides.push(Peptide {
        sequence: Arc::from(&b"PEPTIDE"[..]),
        proteins: vec![Arc::from("P12345")],
        ..Peptide::default()
    });
    let mut areas = HashMap::new();
    areas.insert(
        (PrecursorId::Charged((PeptideIx(0), 2)), false),
        QuantifiedPeak {
            peak: sage_core::lfq::Peak {
                score: 12.5,
                spectral_angle: 0.91,
                q_value: 0.005,
                ..Default::default()
            },
            intensities: vec![Some(42.0), None],
            ms2_confirmed: vec![true, false],
        },
    );

    let bytes = serialize_lfq(&areas, &["run-a".into(), "run-b".into()], &database)?;
    let reader = SerializedFileReader::new(bytes::Bytes::from(bytes))?;
    let metadata = reader
        .metadata()
        .file_metadata()
        .key_value_metadata()
        .expect("LFQ schema metadata");
    assert!(metadata
        .iter()
        .any(|entry| { entry.key == "sage.schema.name" && entry.value.as_deref() == Some("lfq") }));
    assert!(metadata.iter().any(|entry| {
        entry.key == "sage.schema.version" && entry.value.as_deref() == Some("1")
    }));
    let rows = reader
        .get_row_iter(None)?
        .collect::<parquet::errors::Result<Vec<_>>>()?;
    assert_eq!(rows.len(), 2);
    fn values(row: &Row) -> HashMap<&str, &Field> {
        row.get_column_iter()
            .map(|(name, field)| (name.as_str(), field))
            .collect::<HashMap<_, _>>()
    }
    assert_eq!(values(&rows[0])["intensity"], &Field::Double(42.0));
    assert_eq!(values(&rows[0])["ms2_confirmed"], &Field::Bool(true));
    assert_eq!(values(&rows[1])["intensity"], &Field::Null);
    assert_eq!(values(&rows[1])["ms2_confirmed"], &Field::Bool(false));
    Ok(())
}

#[test]
fn lfq_serialization_is_independent_of_hashmap_insertion_order() -> parquet::errors::Result<()> {
    let mut database = IndexedDatabase::default();
    for sequence in [b"PEPTIDE".as_slice(), b"SEQUENCE".as_slice()] {
        database.peptides.push(Peptide {
            sequence: Arc::from(sequence),
            proteins: vec![Arc::from("P12345")],
            ..Peptide::default()
        });
    }

    let first_key = (PrecursorId::Charged((PeptideIx(0), 2)), false);
    let second_key = (PrecursorId::Charged((PeptideIx(1), 3)), false);
    let first_peak = || QuantifiedPeak {
        peak: sage_core::lfq::Peak {
            score: 12.5,
            spectral_angle: 0.91,
            q_value: 0.005,
            ..Default::default()
        },
        intensities: vec![Some(42.0)],
        ms2_confirmed: vec![true],
    };
    let second_peak = || QuantifiedPeak {
        peak: sage_core::lfq::Peak {
            score: 9.0,
            spectral_angle: 0.75,
            q_value: 0.01,
            ..Default::default()
        },
        intensities: vec![Some(21.0)],
        ms2_confirmed: vec![false],
    };

    let mut forward = HashMap::new();
    forward.insert(first_key, first_peak());
    forward.insert(second_key, second_peak());
    let mut reverse = HashMap::new();
    reverse.insert(second_key, second_peak());
    reverse.insert(first_key, first_peak());

    let filenames = ["run-a".into()];
    assert_eq!(
        serialize_lfq(&forward, &filenames, &database)?,
        serialize_lfq(&reverse, &filenames, &database)?
    );
    Ok(())
}

#[test]
fn labeled_lfq_writes_channels_groups_and_reference_ratios() -> parquet::errors::Result<()> {
    let builder: sage_core::database::Builder = serde_json::from_value(serde_json::json!({
        "generate_decoys": false,
        "static_mods": {
            "R": {
                "mass": 0.0,
                "channel_offsets": {"light": 0.0, "heavy": 10.008269}
            }
        }
    }))
    .unwrap();
    let parameters = builder.make_parameters();
    parameters.validate_channels().unwrap();
    let peptides = parameters.peptides_from_tsv("sequence\nPEPTIDER\n");
    let database = parameters.build_from_peptides(peptides);
    let mut areas = HashMap::new();
    for (index, peptide) in database.peptides.iter().enumerate() {
        let intensity = match peptide.label_channel.as_deref() {
            Some("light") => 10.0,
            Some("heavy") => 30.0,
            _ => continue,
        };
        areas.insert(
            (PrecursorId::Charged((PeptideIx(index as u32), 2)), false),
            QuantifiedPeak {
                peak: sage_core::lfq::Peak::default(),
                intensities: vec![Some(intensity)],
                ms2_confirmed: vec![true],
            },
        );
    }

    let bytes = serialize_lfq(&areas, &["run-a".into()], &database)?;
    let reader = SerializedFileReader::new(bytes::Bytes::from(bytes))?;
    let metadata = reader
        .metadata()
        .file_metadata()
        .key_value_metadata()
        .unwrap();
    assert!(metadata.iter().any(|entry| {
        entry.key == "sage.schema.version" && entry.value.as_deref() == Some("2")
    }));
    let rows = reader
        .get_row_iter(None)?
        .collect::<parquet::errors::Result<Vec<_>>>()?;
    assert_eq!(rows.len(), 2);
    for row in rows {
        let values = row
            .get_column_iter()
            .map(|(name, field)| (name.as_str(), field))
            .collect::<HashMap<_, _>>();
        assert_eq!(values["label_group"], &Field::Str("PEPTIDER".into()));
        match values["label_channel"] {
            Field::Str(channel) if channel == "light" => {
                assert_eq!(values["ratio_to_reference"], &Field::Double(1.0));
            }
            Field::Str(channel) if channel == "heavy" => {
                assert_eq!(values["ratio_to_reference"], &Field::Double(3.0));
            }
            channel => panic!("unexpected label channel {channel}"),
        }
    }
    Ok(())
}

#[test]
fn results_preserve_typed_protein_occurrences() -> parquet::errors::Result<()> {
    let mut database = IndexedDatabase::default();
    database.peptides.push(Peptide {
        sequence: Arc::from(&b"PEPTIDE"[..]),
        proteins: vec![Arc::from("P12345"), Arc::from("P67890")],
        protein_sites: Arc::from([
            ProteinOccurrence {
                protein: Arc::from("P12345"),
                start: Some(4),
                prev_aa: Some(b'K'),
                next_aa: Some(b'R'),
            },
            ProteinOccurrence {
                protein: Arc::from("P67890"),
                start: Some(9),
                prev_aa: None,
                next_aa: None,
            },
        ]),
        ..Peptide::default()
    });
    let feature = Feature {
        peptide_idx: PeptideIx(0),
        ..Feature::default()
    };

    let bytes = serialize_features(&[&feature], &[], &["run-a".into()], &database, 1.0)?;
    let reader = SerializedFileReader::new(bytes::Bytes::from(bytes))?;
    let metadata = reader
        .metadata()
        .file_metadata()
        .key_value_metadata()
        .unwrap();
    assert!(metadata.iter().any(|entry| {
        entry.key == "sage.schema.version" && entry.value.as_deref() == Some("3")
    }));
    let rows = reader
        .get_row_iter(None)?
        .collect::<parquet::errors::Result<Vec<_>>>()?;
    let values = rows[0]
        .get_column_iter()
        .map(|(name, field)| (name.as_str(), field))
        .collect::<HashMap<_, _>>();
    assert_eq!(
        values["protein_sites"].to_string(),
        "[{protein: \"P12345\", start: 5, end: 11, prev_aa: \"K\", next_aa: \"R\"}, {protein: \"P67890\", start: 10, end: 16, prev_aa: null, next_aa: null}]"
    );
    Ok(())
}

#[test]
fn labeled_results_write_channel_and_group_columns() -> parquet::errors::Result<()> {
    let builder: sage_core::database::Builder = serde_json::from_value(serde_json::json!({
        "generate_decoys": false,
        "static_mods": {
            "R": {
                "mass": 0.0,
                "channel_offsets": {"light": 0.0, "heavy": 10.008269}
            }
        }
    }))
    .unwrap();
    let parameters = builder.make_parameters();
    let peptides = parameters.peptides_from_tsv("sequence\nPEPTIDER\n");
    let database = parameters.build_from_peptides(peptides);
    let light = database
        .peptides
        .iter()
        .position(|peptide| peptide.label_channel.as_deref() == Some("light"))
        .unwrap();
    let feature = Feature {
        peptide_idx: PeptideIx(light as u32),
        label: 1,
        ..Feature::default()
    };
    let bytes = serialize_features(&[&feature], &[], &["run-a".into()], &database, 1.0)?;
    let reader = SerializedFileReader::new(bytes::Bytes::from(bytes))?;
    let metadata = reader
        .metadata()
        .file_metadata()
        .key_value_metadata()
        .unwrap();
    assert!(metadata.iter().any(|entry| {
        entry.key == "sage.schema.version" && entry.value.as_deref() == Some("4")
    }));
    let rows = reader
        .get_row_iter(None)?
        .collect::<parquet::errors::Result<Vec<_>>>()?;
    let values = rows[0]
        .get_column_iter()
        .map(|(name, field)| (name.as_str(), field))
        .collect::<HashMap<_, _>>();
    assert_eq!(values["label_channel"], &Field::Str("light".into()));
    assert_eq!(values["label_group"], &Field::Str("PEPTIDER".into()));
    Ok(())
}

#[test]
fn spectral_library_has_versioned_long_form_rows() -> parquet::errors::Result<()> {
    let entry = SpectralLibraryEntry {
        library_entry_id: "PEPTIDE/2".into(),
        source_psm_id: 42,
        source_file: "sample.mzML".into(),
        source_spectrum: "scan=42".into(),
        modified_peptide: "PEPTIDE".into(),
        proforma: "PEPTIDE".into(),
        stripped_peptide: "PEPTIDE".into(),
        proteins: "P12345".into(),
        label_channel: None,
        label_group: None,
        label_reference: None,
        precursor_charge: 2,
        precursor_neutral_mass: 798.3854,
        precursor_mz: 400.2,
        retention_time_minutes: 12.5,
        aligned_retention_time_minutes: 11.8,
        ion_mobility: 1.1,
        spectrum_q: 0.001,
        peptide_q: 0.002,
        supporting_psms: 3,
        fragments: vec![
            LibraryFragment {
                kind: Kind::B,
                ordinal: 2,
                charge: 1,
                neutral_loss: 0.0,
                mz: 200.1,
                relative_intensity: 0.5,
            },
            LibraryFragment {
                kind: Kind::Y,
                ordinal: 4,
                charge: 2,
                neutral_loss: 18.010_565,
                mz: 350.2,
                relative_intensity: 1.0,
            },
        ],
    };
    let settings = SpectralLibrarySettings {
        strategy: SpectralLibraryStrategy::Consensus,
        ..SpectralLibrarySettings::default()
    };
    let bytes = serialize_spectral_library(&[entry], &settings)?;
    let search_entries = deserialize_spectral_library(bytes.clone())?;
    assert_eq!(search_entries.len(), 1);
    assert_eq!(search_entries[0].library_entry_id, "PEPTIDE/2");
    assert_eq!(search_entries[0].source_file, "sample.mzML");
    assert_eq!(search_entries[0].source_spectrum, "scan=42");
    assert_eq!(search_entries[0].fragments.len(), 2);
    assert!(!search_entries[0].is_decoy);

    let reader = SerializedFileReader::new(bytes::Bytes::from(bytes))?;
    assert_eq!(reader.metadata().file_metadata().num_rows(), 2);
    assert_eq!(
        reader
            .metadata()
            .file_metadata()
            .schema_descr()
            .num_columns(),
        23
    );
    let metadata = reader
        .metadata()
        .file_metadata()
        .key_value_metadata()
        .expect("spectral-library schema metadata");
    assert!(metadata.iter().any(|entry| {
        entry.key == "sage.schema.name" && entry.value.as_deref() == Some("spectral_library")
    }));
    assert!(metadata.iter().any(|entry| {
        entry.key == "sage.spectral_library.strategy" && entry.value.as_deref() == Some("consensus")
    }));
    assert!(metadata.iter().any(|entry| {
        entry.key == "sage.schema.version" && entry.value.as_deref() == Some("1")
    }));
    let rows = reader
        .get_row_iter(None)?
        .collect::<parquet::errors::Result<Vec<_>>>()?;
    let first = rows[0]
        .get_column_iter()
        .map(|(name, field)| (name.as_str(), field))
        .collect::<HashMap<_, _>>();
    assert_eq!(first["library_entry_id"], &Field::Str("PEPTIDE/2".into()));
    assert_eq!(first["fragment_type"], &Field::Str("b".into()));
    assert_eq!(first["relative_intensity"], &Field::Float(0.5));
    Ok(())
}

#[test]
fn labeled_spectral_library_round_trips_channel_metadata() -> parquet::errors::Result<()> {
    let entry = SpectralLibraryEntry {
        library_entry_id: "PEPTIDER[+10.008269]/2".into(),
        source_psm_id: 1,
        source_file: "sample.mzML".into(),
        source_spectrum: "scan=1".into(),
        modified_peptide: "PEPTIDER[Arg10]".into(),
        proforma: "PEPTIDER[+10.008269]".into(),
        stripped_peptide: "PEPTIDER".into(),
        proteins: "P1".into(),
        label_channel: Some("heavy".into()),
        label_group: Some("PEPTIDER".into()),
        label_reference: Some("light".into()),
        precursor_charge: 2,
        precursor_neutral_mass: 966.0,
        precursor_mz: 484.0,
        retention_time_minutes: 10.0,
        aligned_retention_time_minutes: 10.0,
        ion_mobility: 1.0,
        spectrum_q: 0.001,
        peptide_q: 0.001,
        supporting_psms: 1,
        fragments: vec![LibraryFragment {
            kind: Kind::Y,
            ordinal: 3,
            charge: 1,
            neutral_loss: 0.0,
            mz: 400.0,
            relative_intensity: 1.0,
        }],
    };
    let bytes = serialize_spectral_library(&[entry], &SpectralLibrarySettings::default())?;
    let entries = deserialize_spectral_library(bytes)?;
    assert_eq!(entries[0].label_channel.as_deref(), Some("heavy"));
    assert_eq!(entries[0].label_group.as_deref(), Some("PEPTIDER"));
    assert_eq!(entries[0].label_reference.as_deref(), Some("light"));
    Ok(())
}

#[test]
fn ptm_library_round_trip() {
    let sites = vec![
        PtmLibrarySite {
            protein: Arc::from("P12345"),
            position: 41,
            residue: b'S',
            modification: Arc::from("Phospho"),
        },
        PtmLibrarySite {
            protein: Arc::from("P12345"),
            position: 41,
            residue: b'S',
            modification: Arc::from("Phospho"),
        },
    ];
    let encoded = serialize_ptm_library(&sites).unwrap();
    let decoded = deserialize_ptm_library(encoded).unwrap();
    assert_eq!(decoded.len(), 1);
    assert_eq!(decoded.sites_for("P12345")[0], sites[0]);
}

#[test]
fn deserialize_custom_cleavage_library() -> parquet::errors::Result<()> {
    let schema = parquet::schema::parser::parse_message_type(
        r#"
            message schema {
                required byte_array protein (utf8);
                required int64 position;
                required byte_array context (utf8);
            }
            "#,
    )
    .unwrap();
    let mut writer = SerializedFileWriter::new(
        Vec::new(),
        schema.into(),
        WriterProperties::default().into(),
    )
    .unwrap();
    let mut row_group = writer.next_row_group().unwrap();
    write_required_column!(
        row_group,
        vec![ByteArray::from("P1"), ByteArray::from("P2")],
        ByteArrayType
    );
    write_required_column!(row_group, vec![4_i64, 0_i64], Int64Type);
    write_required_column!(
        row_group,
        vec![ByteArray::from("PEPK|TIDE"), ByteArray::from("A|CDE")],
        ByteArrayType
    );
    row_group.close().unwrap();
    let bytes = writer.into_inner().unwrap();

    let library = deserialize_custom_cleavage_sites(bytes).unwrap();
    let fasta = sage_core::fasta::Fasta::parse(">P1\nMPEPKTIDER\n>P2\nACDE\n".into(), "rev_", true)
        .unwrap();
    let validated = library.validate(&fasta).unwrap();
    assert_eq!(validated.boundaries_for("P1"), &[5]);
    assert_eq!(validated.boundaries_for("P2"), &[1]);
    assert_eq!(validated.sites_without_context, 0);
    Ok(())
}

#[test]
fn serialize_ptm_site_reports() {
    let ptm = serialize_ptm_sites(&[PtmSiteRecord {
        psm_id: 42,
        filename: "sample.mzML".into(),
        scannr: "scan=42".into(),
        peptide: "AAS[+79.966]AATAA".into(),
        proteins: "P12345".into(),
        charge: 2,
        spectrum_q: 0.005,
        peptide_q: 0.006,
        modification: "Phospho".into(),
        modification_mass: 79.96633,
        position: 3,
        residue: "S".into(),
        localization_probability: 0.982,
        delta_localization_score: 18.7,
        target_decoy_score: 21.0,
        localization_q_value: 0.01,
        candidate_sites: 2,
        site_determining_ions_matched: 6,
        site_determining_ions_total: 8,
        site_probabilities: "S3:0.982;T6:0.018".into(),
    }])
    .unwrap();
    let reader = SerializedFileReader::new(bytes::Bytes::from(ptm)).unwrap();
    assert_eq!(reader.metadata().file_metadata().num_rows(), 1);
    assert_eq!(
        reader
            .metadata()
            .file_metadata()
            .schema_descr()
            .num_columns(),
        20
    );

    let protein = serialize_protein_sites(&[ProteinSiteRecord {
        protein: "P12345".into(),
        peptide: "AAS[+79.966]AATAA".into(),
        residue: "S".into(),
        position_in_peptide: 3,
        modification: "Phospho".into(),
        modification_mass: 79.96633,
        num_psms: 2,
        best_localization_probability: 0.982,
        best_delta_localization_score: 18.7,
        best_localization_q_value: 0.01,
        best_spectrum_q: 0.005,
    }])
    .unwrap();
    let reader = SerializedFileReader::new(bytes::Bytes::from(protein)).unwrap();
    assert_eq!(reader.metadata().file_metadata().num_rows(), 1);
    assert_eq!(
        reader
            .metadata()
            .file_metadata()
            .schema_descr()
            .num_columns(),
        11
    );
}
