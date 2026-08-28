use super::*;

pub fn build_spectral_library_schema() -> parquet::errors::Result<Type> {
    parquet::schema::parser::parse_message_type(include_str!(
        "../../../../schemas/spectral_library.sage.v1.parquet.schema"
    ))
}

fn build_spectral_library_schema_version(has_labels: bool) -> parquet::errors::Result<Type> {
    parquet::schema::parser::parse_message_type(if has_labels {
        include_str!("../../../../schemas/spectral_library.sage.v2.parquet.schema")
    } else {
        include_str!("../../../../schemas/spectral_library.sage.v1.parquet.schema")
    })
}

/// Serialize one long-form row per selected transition in the empirical library.
pub fn serialize_spectral_library(
    entries: &[SpectralLibraryEntry],
    settings: &SpectralLibrarySettings,
) -> parquet::errors::Result<Vec<u8>> {
    let has_labels = entries.iter().any(|entry| entry.label_channel.is_some());
    let schema = build_spectral_library_schema_version(has_labels)?;
    let options = WriterProperties::builder()
        .set_compression(parquet::basic::Compression::ZSTD(ZstdLevel::try_new(3)?))
        .set_key_value_metadata(Some(vec![
            KeyValue::new("sage.schema.name".into(), Some("spectral_library".into())),
            KeyValue::new(
                "sage.schema.version".into(),
                Some(if has_labels { "2" } else { "1" }.into()),
            ),
            KeyValue::new(
                "sage.spectral_library.strategy".into(),
                Some(
                    match settings.strategy {
                        SpectralLibraryStrategy::BestPsm => "best_psm",
                        SpectralLibraryStrategy::Consensus => "consensus",
                    }
                    .into(),
                ),
            ),
            KeyValue::new(
                "sage.spectral_library.psm_q_max".into(),
                Some(settings.psm_q_value.to_string()),
            ),
            KeyValue::new(
                "sage.spectral_library.peptide_q_max".into(),
                Some(settings.peptide_q_value.to_string()),
            ),
            KeyValue::new(
                "sage.spectral_library.min_consensus_psms".into(),
                Some(settings.min_consensus_psms.to_string()),
            ),
            KeyValue::new(
                "sage.spectral_library.min_fragment_frequency".into(),
                Some(settings.min_fragment_frequency.to_string()),
            ),
        ]))
        .build();
    let mut writer = SerializedFileWriter::new(Vec::new(), schema.into(), options.into())?;
    let rows = entries
        .iter()
        .flat_map(|entry| {
            entry
                .fragments
                .iter()
                .map(move |fragment| (entry, fragment))
        })
        .collect::<Vec<_>>();

    for rows in rows.chunks(65_536) {
        let mut rg = writer.next_row_group()?;
        write_required_column!(
            rg,
            rows.iter()
                .map(|(entry, _)| entry.library_entry_id.as_str().into())
                .collect::<Vec<ByteArray>>(),
            ByteArrayType
        );
        write_required_column!(
            rg,
            rows.iter()
                .map(|(entry, _)| entry.source_psm_id as i64)
                .collect::<Vec<_>>(),
            Int64Type
        );
        write_required_column!(
            rg,
            rows.iter()
                .map(|(entry, _)| entry.source_file.as_str().into())
                .collect::<Vec<ByteArray>>(),
            ByteArrayType
        );
        write_required_column!(
            rg,
            rows.iter()
                .map(|(entry, _)| entry.source_spectrum.as_str().into())
                .collect::<Vec<ByteArray>>(),
            ByteArrayType
        );
        write_required_column!(
            rg,
            rows.iter()
                .map(|(entry, _)| entry.modified_peptide.as_str().into())
                .collect::<Vec<ByteArray>>(),
            ByteArrayType
        );
        write_required_column!(
            rg,
            rows.iter()
                .map(|(entry, _)| entry.proforma.as_str().into())
                .collect::<Vec<ByteArray>>(),
            ByteArrayType
        );
        write_required_column!(
            rg,
            rows.iter()
                .map(|(entry, _)| entry.stripped_peptide.as_str().into())
                .collect::<Vec<ByteArray>>(),
            ByteArrayType
        );
        write_required_column!(
            rg,
            rows.iter()
                .map(|(entry, _)| entry.proteins.as_str().into())
                .collect::<Vec<ByteArray>>(),
            ByteArrayType
        );
        if has_labels {
            write_required_column!(
                rg,
                rows.iter()
                    .map(|(entry, _)| entry.label_channel.as_deref().unwrap_or("").into())
                    .collect::<Vec<ByteArray>>(),
                ByteArrayType
            );
            write_required_column!(
                rg,
                rows.iter()
                    .map(|(entry, _)| entry.label_group.as_deref().unwrap_or("").into())
                    .collect::<Vec<ByteArray>>(),
                ByteArrayType
            );
            write_required_column!(
                rg,
                rows.iter()
                    .map(|(entry, _)| entry.label_reference.as_deref().unwrap_or("").into())
                    .collect::<Vec<ByteArray>>(),
                ByteArrayType
            );
        }
        write_required_column!(
            rg,
            rows.iter()
                .map(|(entry, _)| i32::from(entry.precursor_charge))
                .collect::<Vec<_>>(),
            Int32Type
        );
        write_required_column!(
            rg,
            rows.iter()
                .map(|(entry, _)| entry.precursor_neutral_mass)
                .collect::<Vec<_>>(),
            FloatType
        );
        write_required_column!(
            rg,
            rows.iter()
                .map(|(entry, _)| entry.precursor_mz)
                .collect::<Vec<_>>(),
            FloatType
        );
        write_required_column!(
            rg,
            rows.iter()
                .map(|(entry, _)| entry.retention_time_minutes)
                .collect::<Vec<_>>(),
            FloatType
        );
        write_required_column!(
            rg,
            rows.iter()
                .map(|(entry, _)| entry.aligned_retention_time_minutes)
                .collect::<Vec<_>>(),
            FloatType
        );
        write_required_column!(
            rg,
            rows.iter()
                .map(|(entry, _)| entry.ion_mobility)
                .collect::<Vec<_>>(),
            FloatType
        );
        write_required_column!(
            rg,
            rows.iter()
                .map(|(entry, _)| entry.spectrum_q)
                .collect::<Vec<_>>(),
            FloatType
        );
        write_required_column!(
            rg,
            rows.iter()
                .map(|(entry, _)| entry.peptide_q)
                .collect::<Vec<_>>(),
            FloatType
        );
        write_required_column!(
            rg,
            rows.iter()
                .map(|(entry, _)| entry.supporting_psms as i32)
                .collect::<Vec<_>>(),
            Int32Type
        );
        write_required_column!(
            rg,
            rows.iter()
                .map(|(_, fragment)| match fragment.kind {
                    Kind::A => "a".into(),
                    Kind::B => "b".into(),
                    Kind::C => "c".into(),
                    Kind::X => "x".into(),
                    Kind::Y => "y".into(),
                    Kind::Z => "z".into(),
                })
                .collect::<Vec<ByteArray>>(),
            ByteArrayType
        );
        write_required_column!(
            rg,
            rows.iter()
                .map(|(_, fragment)| fragment.ordinal)
                .collect::<Vec<_>>(),
            Int32Type
        );
        write_required_column!(
            rg,
            rows.iter()
                .map(|(_, fragment)| fragment.charge)
                .collect::<Vec<_>>(),
            Int32Type
        );
        write_required_column!(
            rg,
            rows.iter()
                .map(|(_, fragment)| fragment.neutral_loss)
                .collect::<Vec<_>>(),
            FloatType
        );
        write_required_column!(
            rg,
            rows.iter()
                .map(|(_, fragment)| fragment.mz)
                .collect::<Vec<_>>(),
            FloatType
        );
        write_required_column!(
            rg,
            rows.iter()
                .map(|(_, fragment)| fragment.relative_intensity)
                .collect::<Vec<_>>(),
            FloatType
        );
        rg.close()?;
    }
    writer.into_inner().map(|bytes| bytes.to_vec())
}
