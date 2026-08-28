use super::*;

pub fn build_lfq_schema() -> parquet::errors::Result<Type> {
    parquet::schema::parser::parse_message_type(include_str!(
        "../../../../schemas/lfq.v1.parquet.schema"
    ))
}

fn build_lfq_schema_version(has_labels: bool) -> parquet::errors::Result<Type> {
    parquet::schema::parser::parse_message_type(if has_labels {
        include_str!("../../../../schemas/lfq.v2.parquet.schema")
    } else {
        include_str!("../../../../schemas/lfq.v1.parquet.schema")
    })
}

pub fn serialize_lfq<H: BuildHasher>(
    areas: &HashMap<(PrecursorId, bool), QuantifiedPeak, H>,
    filenames: &[String],
    database: &IndexedDatabase,
) -> parquet::errors::Result<Vec<u8>> {
    if let Some((_, quantified)) = areas.iter().find(|(_, quantified)| {
        quantified.intensities.len() != filenames.len()
            || quantified.ms2_confirmed.len() != filenames.len()
    }) {
        return Err(ParquetError::General(format!(
            "LFQ row has {} intensities and {} MS2 evidence values for {} files",
            quantified.intensities.len(),
            quantified.ms2_confirmed.len(),
            filenames.len()
        )));
    }
    let mut rows = areas.iter().collect::<Vec<_>>();
    rows.sort_unstable_by(|left, right| left.0.cmp(right.0));
    let has_labels = !database.label_channels.is_empty();
    let schema = build_lfq_schema_version(has_labels)?;
    let reference_intensities = database
        .label_reference
        .as_deref()
        .map(|reference| {
            areas
                .iter()
                .filter_map(|((id, decoy), quantified)| {
                    let (peptide_idx, charge) = match id {
                        PrecursorId::Combined(peptide_idx) => (*peptide_idx, None),
                        PrecursorId::Charged((peptide_idx, charge)) => {
                            (*peptide_idx, Some(*charge))
                        }
                    };
                    let peptide = &database[peptide_idx];
                    (peptide.label_channel.as_deref() == Some(reference)).then(|| {
                        (
                            (peptide.label_group(), charge, *decoy),
                            quantified.intensities.as_slice(),
                        )
                    })
                })
                .collect::<HashMap<_, _>>()
        })
        .unwrap_or_default();

    let options = WriterProperties::builder()
        .set_compression(parquet::basic::Compression::ZSTD(ZstdLevel::try_new(3)?))
        .set_key_value_metadata(Some(vec![
            KeyValue::new("sage.schema.name".into(), Some("lfq".into())),
            KeyValue::new(
                "sage.schema.version".into(),
                Some(if has_labels { "2" } else { "1" }.into()),
            ),
        ]))
        .build();

    let buf = Vec::new();
    let mut writer = SerializedFileWriter::new(buf, schema.into(), options.into())?;
    let mut rg = writer.next_row_group()?;

    if let Some(mut col) = rg.next_column()? {
        let values = rows
            .iter()
            .copied()
            .flat_map(|((id, _), _)| {
                let peptide_idx = match id {
                    PrecursorId::Combined(x) | PrecursorId::Charged((x, _)) => x,
                };
                let val = database[*peptide_idx].to_string().as_bytes().into();
                std::iter::repeat_n(val, filenames.len())
            })
            .collect::<Vec<_>>();

        col.typed::<ByteArrayType>()
            .write_batch(&values, None, None)?;
        col.close()?;
    }

    if let Some(mut col) = rg.next_column()? {
        let values = rows
            .iter()
            .copied()
            .flat_map(|((id, _), _)| {
                let peptide_idx = match id {
                    PrecursorId::Combined(x) | PrecursorId::Charged((x, _)) => x,
                };
                let val = database[*peptide_idx].sequence.as_ref().into();
                std::iter::repeat_n(val, filenames.len())
            })
            .collect::<Vec<_>>();

        col.typed::<ByteArrayType>()
            .write_batch(&values, None, None)?;
        col.close()?;
    }

    if has_labels {
        if let Some(mut col) = rg.next_column()? {
            let values = rows
                .iter()
                .copied()
                .flat_map(|((id, _), _)| {
                    let peptide_idx = match id {
                        PrecursorId::Combined(x) | PrecursorId::Charged((x, _)) => x,
                    };
                    let value: ByteArray = database[*peptide_idx]
                        .label_channel
                        .as_deref()
                        .unwrap_or("")
                        .into();
                    std::iter::repeat_n(value, filenames.len())
                })
                .collect::<Vec<_>>();
            col.typed::<ByteArrayType>()
                .write_batch(&values, None, None)?;
            col.close()?;
        }

        if let Some(mut col) = rg.next_column()? {
            let values = rows
                .iter()
                .copied()
                .flat_map(|((id, _), _)| {
                    let peptide_idx = match id {
                        PrecursorId::Combined(x) | PrecursorId::Charged((x, _)) => x,
                    };
                    let value: ByteArray = database[*peptide_idx].label_group().as_bytes().into();
                    std::iter::repeat_n(value, filenames.len())
                })
                .collect::<Vec<_>>();
            col.typed::<ByteArrayType>()
                .write_batch(&values, None, None)?;
            col.close()?;
        }
    }

    if let Some(mut col) = rg.next_column()? {
        let mut values = Vec::with_capacity(areas.len() * filenames.len());
        let mut def_levels = Vec::with_capacity(areas.len() * filenames.len());

        for ((id, _), _) in rows.iter().copied() {
            match id {
                PrecursorId::Combined(_) => {
                    def_levels.extend(std::iter::repeat_n(0, filenames.len()));
                }
                PrecursorId::Charged((_, charge)) => {
                    values.extend(std::iter::repeat_n(*charge as i32, filenames.len()));
                    def_levels.extend(std::iter::repeat_n(1, filenames.len()));
                }
            }
        }

        col.typed::<Int32Type>()
            .write_batch(&values, Some(&def_levels), None)?;
        col.close()?;
    }

    if let Some(mut col) = rg.next_column()? {
        let values = rows
            .iter()
            .copied()
            .flat_map(|((id, _), _)| {
                let peptide_idx = match id {
                    PrecursorId::Combined(x) | PrecursorId::Charged((x, _)) => x,
                };
                let val = database[*peptide_idx]
                    .proteins(&database.decoy_tag, database.generate_decoys)
                    .as_str()
                    .into();
                std::iter::repeat_n(val, filenames.len())
            })
            .collect::<Vec<_>>();

        col.typed::<ByteArrayType>()
            .write_batch(&values, None, None)?;
        col.close()?;
    }

    if let Some(mut col) = rg.next_column()? {
        let values = rows
            .iter()
            .copied()
            .flat_map(|((_, decoy), _)| std::iter::repeat_n(*decoy, filenames.len()))
            .collect::<Vec<_>>();

        col.typed::<BoolType>().write_batch(&values, None, None)?;
        col.close()?;
    }

    if let Some(mut col) = rg.next_column()? {
        let values = rows
            .iter()
            .copied()
            .flat_map(|(_, quantified)| {
                std::iter::repeat_n(quantified.peak.q_value, filenames.len())
            })
            .collect::<Vec<_>>();

        col.typed::<FloatType>().write_batch(&values, None, None)?;
        col.close()?;
    }

    if let Some(mut col) = rg.next_column()? {
        let values = rows
            .iter()
            .copied()
            .flat_map(|(_, quantified)| std::iter::repeat_n(quantified.peak.score, filenames.len()))
            .collect::<Vec<_>>();

        col.typed::<DoubleType>().write_batch(&values, None, None)?;
        col.close()?;
    }

    if let Some(mut col) = rg.next_column()? {
        let values = rows
            .iter()
            .copied()
            .flat_map(|(_, quantified)| {
                std::iter::repeat_n(quantified.peak.spectral_angle, filenames.len())
            })
            .collect::<Vec<_>>();

        col.typed::<DoubleType>().write_batch(&values, None, None)?;
        col.close()?;
    }

    if let Some(mut col) = rg.next_column()? {
        let values = rows
            .iter()
            .flat_map(|_| filenames.iter().map(|filename| filename.as_bytes().into()))
            .collect::<Vec<_>>();

        col.typed::<ByteArrayType>()
            .write_batch(&values, None, None)?;

        col.close()?;
    }

    if let Some(mut col) = rg.next_column()? {
        let mut values = Vec::with_capacity(areas.len() * filenames.len());
        let mut def_levels = Vec::with_capacity(areas.len() * filenames.len());
        for (_, quantified) in rows.iter().copied() {
            for intensity in &quantified.intensities {
                if let Some(intensity) = intensity {
                    values.push(*intensity);
                    def_levels.push(1);
                } else {
                    def_levels.push(0);
                }
            }
        }

        col.typed::<DoubleType>()
            .write_batch(&values, Some(&def_levels), None)?;
        col.close()?;
    }

    if has_labels {
        if let Some(mut col) = rg.next_column()? {
            let mut values = Vec::new();
            let mut def_levels = Vec::with_capacity(areas.len() * filenames.len());
            for ((id, decoy), quantified) in rows.iter().copied() {
                let (peptide_idx, charge) = match id {
                    PrecursorId::Combined(peptide_idx) => (*peptide_idx, None),
                    PrecursorId::Charged((peptide_idx, charge)) => (*peptide_idx, Some(*charge)),
                };
                let reference = reference_intensities.get(&(
                    database[peptide_idx].label_group(),
                    charge,
                    *decoy,
                ));
                for file_id in 0..filenames.len() {
                    let ratio = quantified.intensities[file_id].zip(
                        reference
                            .and_then(|intensities| intensities.get(file_id))
                            .copied()
                            .flatten(),
                    );
                    if let Some((intensity, reference)) =
                        ratio.filter(|(_, reference)| *reference > 0.0)
                    {
                        values.push(intensity / reference);
                        def_levels.push(1);
                    } else {
                        def_levels.push(0);
                    }
                }
            }
            col.typed::<DoubleType>()
                .write_batch(&values, Some(&def_levels), None)?;
            col.close()?;
        }
    }

    if let Some(mut col) = rg.next_column()? {
        let values = rows
            .iter()
            .map(|(_, quantified)| *quantified)
            .flat_map(|quantified| quantified.ms2_confirmed.iter().copied())
            .collect::<Vec<_>>();

        col.typed::<BoolType>().write_batch(&values, None, None)?;
        col.close()?;
    }

    rg.close()?;
    writer.into_inner()
}
