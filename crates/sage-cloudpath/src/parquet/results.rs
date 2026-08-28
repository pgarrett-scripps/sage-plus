use super::*;

pub fn build_schema() -> Result<Type, parquet::errors::ParquetError> {
    parquet::schema::parser::parse_message_type(include_str!(
        "../../../../schemas/results.sage.v1.parquet.schema"
    ))
}

fn build_results_schema(has_labels: bool) -> Result<Type, parquet::errors::ParquetError> {
    parquet::schema::parser::parse_message_type(if has_labels {
        include_str!("../../../../schemas/results.sage.v2.parquet.schema")
    } else {
        include_str!("../../../../schemas/results.sage.v1.parquet.schema")
    })
}

struct OutputProteinSite {
    protein: ByteArray,
    start: i32,
    end: i32,
    prev_aa: Option<ByteArray>,
    next_aa: Option<ByteArray>,
}

fn output_protein_sites(feature: &Feature, database: &IndexedDatabase) -> Vec<OutputProteinSite> {
    let peptide = &database[feature.peptide_idx];
    peptide
        .protein_sites
        .iter()
        .filter_map(|site| {
            let start = site.start?;
            let protein = if peptide.decoy && database.generate_decoys {
                format!("{}{}", database.decoy_tag, site.protein)
            } else {
                site.protein.to_string()
            };
            Some(OutputProteinSite {
                protein: protein.into_bytes().into(),
                start: start.saturating_add(1) as i32,
                end: start.saturating_add(peptide.sequence.len() as u32) as i32,
                prev_aa: site.prev_aa.map(|aa| vec![aa].into()),
                next_aa: site.next_aa.map(|aa| vec![aa].into()),
            })
        })
        .collect()
}

/// Caller must guarantee that `reporter_ions` is not an empty slice
fn write_reporter_ions(
    mut column: SerializedColumnWriter,
    features: &[&Feature],
    reporter_ions: &[TmtQuant],
) -> parquet::errors::Result<()> {
    let mut scan_map = HashMap::new();

    for r in reporter_ions {
        scan_map.entry((r.file_id, &r.spec_id)).or_insert(r);
    }

    // Caller guarantees `reporter_ions` is not empty
    let channels = reporter_ions[0].peaks.len();

    // https://docs.rs/parquet/44.0.0/parquet/column/index.html
    // Using the low level API here is not very pleasant...
    let def_levels = vec![3; channels];
    let mut rep_levels = vec![1; channels];
    rep_levels[0] = 0;

    let col = column.typed::<FloatType>();
    for feature in features {
        if let Some(rs) = scan_map.get(&(feature.file_id, &feature.spec_id)) {
            col.write_batch(&rs.peaks, Some(&def_levels), Some(&rep_levels))?;
        } else {
            col.write_batch(&[], Some(&[0]), Some(&[0]))?;
        }
    }

    column.close()?;
    Ok(())
}

fn write_null_column(
    mut column: SerializedColumnWriter,
    length: usize,
) -> Result<usize, parquet::errors::ParquetError> {
    let levels = vec![0i16; length];
    let wrote = column
        .typed::<FloatType>()
        .write_batch(&[], Some(&levels), Some(&levels))?;
    column.close().map(|_| wrote)
}

pub fn serialize_features(
    features: &[&Feature],
    reporter_ions: &[TmtQuant],
    filenames: &[String],
    database: &IndexedDatabase,
    output_psm_q_value: f32,
) -> Result<Vec<u8>, parquet::errors::ParquetError> {
    let has_labels = !database.label_channels.is_empty();
    let schema = build_results_schema(has_labels)?;

    let options = WriterProperties::builder()
        .set_compression(parquet::basic::Compression::ZSTD(ZstdLevel::try_new(3)?))
        .set_key_value_metadata(Some(vec![
            KeyValue::new("sage.schema.name".into(), Some("results.sage".into())),
            KeyValue::new(
                "sage.schema.version".into(),
                Some(if has_labels { "2" } else { "1" }.into()),
            ),
            KeyValue::new(
                "sage.output_filter.spectrum_q_max".into(),
                Some(output_psm_q_value.to_string()),
            ),
        ]))
        .build();

    let buf = Vec::new();
    let mut writer = SerializedFileWriter::new(buf, schema.into(), options.into())?;

    for features in features.chunks(65536) {
        let mut rg = writer.next_row_group()?;
        macro_rules! write_col {
            ($field:ident, $ty:ident) => {
                if let Some(mut col) = rg.next_column()? {
                    col.typed::<$ty>().write_batch(
                        &features
                            .iter()
                            .map(|f| f.$field as <$ty as DataType>::T)
                            .collect::<Vec<_>>(),
                        None,
                        None,
                    )?;
                    col.close()?;
                }
            };
            ($lambda:expr, $ty:ident) => {
                if let Some(mut col) = rg.next_column()? {
                    col.typed::<$ty>().write_batch(
                        &features.iter().map($lambda).collect::<Vec<_>>(),
                        None,
                        None,
                    )?;
                    col.close()?;
                }
            };
        }

        write_col!(|f: &&Feature| f.psm_id as i64, Int64Type);
        write_col!(
            |f: &&Feature| filenames[f.file_id].as_str().into(),
            ByteArrayType
        );
        write_col!(|f: &&Feature| f.spec_id.as_str().into(), ByteArrayType);
        write_col!(
            |f: &&Feature| database[f.peptide_idx].to_string().as_bytes().into(),
            ByteArrayType
        );
        write_col!(
            |f: &&Feature| f.ambiguity_sequence.as_str().into(),
            ByteArrayType
        );
        write_col!(mass_shift, FloatType);
        write_col!(
            |f: &&Feature| database[f.peptide_idx].sequence.as_ref().into(),
            ByteArrayType
        );
        if has_labels {
            write_col!(
                |f: &&Feature| database[f.peptide_idx]
                    .label_channel
                    .as_deref()
                    .unwrap_or("")
                    .into(),
                ByteArrayType
            );
            write_col!(
                |f: &&Feature| database[f.peptide_idx].label_group().as_bytes().into(),
                ByteArrayType
            );
        }
        write_col!(
            |f: &&Feature| database[f.peptide_idx]
                .proteins(&database.decoy_tag, database.generate_decoys)
                .as_str()
                .into(),
            ByteArrayType
        );

        let protein_sites = features
            .iter()
            .map(|feature| output_protein_sites(feature, database))
            .collect::<Vec<_>>();

        macro_rules! write_required_protein_site_column {
            ($field:ident, $ty:ident) => {
                if let Some(mut column) = rg.next_column()? {
                    let mut values = Vec::new();
                    let mut definition_levels = Vec::new();
                    let mut repetition_levels = Vec::new();
                    for sites in &protein_sites {
                        if sites.is_empty() {
                            definition_levels.push(0);
                            repetition_levels.push(0);
                        } else {
                            for (index, site) in sites.iter().enumerate() {
                                values.push(site.$field.clone());
                                definition_levels.push(1);
                                repetition_levels.push(i16::from(index > 0));
                            }
                        }
                    }
                    column.typed::<$ty>().write_batch(
                        &values,
                        Some(&definition_levels),
                        Some(&repetition_levels),
                    )?;
                    column.close()?;
                }
            };
        }

        macro_rules! write_optional_protein_site_column {
            ($field:ident) => {
                if let Some(mut column) = rg.next_column()? {
                    let mut values = Vec::new();
                    let mut definition_levels = Vec::new();
                    let mut repetition_levels = Vec::new();
                    for sites in &protein_sites {
                        if sites.is_empty() {
                            definition_levels.push(0);
                            repetition_levels.push(0);
                        } else {
                            for (index, site) in sites.iter().enumerate() {
                                repetition_levels.push(i16::from(index > 0));
                                if let Some(value) = &site.$field {
                                    values.push(value.clone());
                                    definition_levels.push(2);
                                } else {
                                    definition_levels.push(1);
                                }
                            }
                        }
                    }
                    column.typed::<ByteArrayType>().write_batch(
                        &values,
                        Some(&definition_levels),
                        Some(&repetition_levels),
                    )?;
                    column.close()?;
                }
            };
        }

        write_required_protein_site_column!(protein, ByteArrayType);
        write_required_protein_site_column!(start, Int32Type);
        write_required_protein_site_column!(end, Int32Type);
        write_optional_protein_site_column!(prev_aa);
        write_optional_protein_site_column!(next_aa);

        write_col!(
            |f: &&Feature| f.protein_groups.as_deref().unwrap_or("").into(),
            ByteArrayType
        );
        write_col!(
            |f: &&Feature| database[f.peptide_idx].proteins.len() as i32,
            Int32Type
        );
        write_col!(num_protein_groups, Int32Type);
        write_col!(rank, Int32Type);
        write_col!(|f: &&Feature| f.label == -1, BoolType);
        write_col!(expmass, FloatType);
        write_col!(calcmass, FloatType);
        write_col!(charge, Int32Type);
        write_col!(peptide_len, Int32Type);
        write_col!(missed_cleavages, Int32Type);
        write_col!(
            |f: &&Feature| database[f.peptide_idx].semi_enzymatic,
            BoolType
        );
        write_col!(ms2_intensity, FloatType);
        write_col!(isotope_error, FloatType);
        write_col!(delta_mass, FloatType);
        write_col!(average_ppm, FloatType);
        write_col!(aligned_delta_mass, FloatType);
        write_col!(aligned_average_ppm, FloatType);
        write_col!(hyperscore, FloatType);
        write_col!(delta_next, FloatType);
        write_col!(delta_best, FloatType);
        write_col!(rt, FloatType);
        write_col!(aligned_rt, FloatType);
        write_col!(predicted_rt, FloatType);
        write_col!(delta_rt_model, FloatType);
        write_col!(ims, FloatType);
        write_col!(predicted_ims, FloatType);
        write_col!(delta_ims_model, FloatType);
        write_col!(matched_peaks, Int32Type);
        write_col!(longest_b, Int32Type);
        write_col!(longest_y, Int32Type);
        write_col!(longest_y_pct, FloatType);
        write_col!(matched_intensity_pct, FloatType);
        // Retain the former library-search columns as zeros for Parquet schema compatibility.
        write_col!(|_: &&Feature| 0.0_f32, FloatType);
        write_col!(|_: &&Feature| 0.0_f32, FloatType);
        write_col!(|_: &&Feature| 0.0_f32, FloatType);
        write_col!(scored_candidates, Int32Type);
        write_col!(poisson, FloatType);
        write_col!(discriminant_score, FloatType);
        write_col!(posterior_error, FloatType);
        write_col!(spectrum_q, FloatType);
        write_col!(peptide_q, FloatType);
        write_col!(protein_q, FloatType);
        write_col!(protein_group_q, FloatType);

        if let Some(col) = rg.next_column()? {
            if reporter_ions.is_empty() {
                write_null_column(col, features.len())?;
            } else {
                write_reporter_ions(col, features, reporter_ions)?;
            }
        }

        rg.close()?;
    }
    writer.into_inner()
}

pub fn build_matched_fragment_schema() -> parquet::errors::Result<Type> {
    let msg = r#"
        message schema {
            required int64 psm_id;
            required byte_array fragment_type (utf8);
            required int32 fragment_ordinals;
            required int32 fragment_charge;
            required float fragment_mz_experimental;
            required float fragment_mz_calculated;
            required float neutral_loss;
            required float fragment_intensity;
        }
    "#;

    parquet::schema::parser::parse_message_type(msg)
}

pub fn serialize_matched_fragments(
    features: &[&Feature],
    output_psm_q_value: f32,
) -> Result<Vec<u8>, parquet::errors::ParquetError> {
    let schema = build_matched_fragment_schema()?;

    let options = WriterProperties::builder()
        .set_compression(parquet::basic::Compression::ZSTD(ZstdLevel::try_new(3)?))
        .set_key_value_metadata(Some(vec![KeyValue::new(
            "sage.output_filter.spectrum_q_max".into(),
            Some(output_psm_q_value.to_string()),
        )]))
        .build();

    let buf = Vec::new();

    let mut writer = SerializedFileWriter::new(buf, schema.into(), options.into())?;

    for features in features.chunks(65536) {
        let mut rg = writer.next_row_group()?;

        if let Some(mut col) = rg.next_column()? {
            let psm_ids = features
                .iter()
                .flat_map(|f| {
                    std::iter::repeat_n(
                        f.psm_id as i64,
                        f.fragments
                            .as_ref()
                            .map(|fragments| fragments.fragment_ordinals.len())
                            .unwrap_or_default(),
                    )
                })
                .collect::<Vec<_>>();

            col.typed::<Int64Type>().write_batch(&psm_ids, None, None)?;
            col.close()?;
        }

        if let Some(mut col) = rg.next_column()? {
            let fragment_types = features
                .iter()
                .flat_map(|f| {
                    f.fragments
                        .as_ref()
                        .map(|fragments| fragments.kinds.iter().copied())
                })
                .flatten()
                .map(|kind| match kind {
                    Kind::A => "a".as_bytes().into(),
                    Kind::B => "b".as_bytes().into(),
                    Kind::C => "c".as_bytes().into(),
                    Kind::X => "x".as_bytes().into(),
                    Kind::Y => "y".as_bytes().into(),
                    Kind::Z => "z".as_bytes().into(),
                })
                .collect::<Vec<ByteArray>>();

            col.typed::<ByteArrayType>()
                .write_batch(&fragment_types, None, None)?;
            col.close()?;
        }

        if let Some(mut col) = rg.next_column()? {
            let fragment_ordinals = features
                .iter()
                .flat_map(|f| {
                    f.fragments
                        .as_ref()
                        .map(|fragments| fragments.fragment_ordinals.iter().copied())
                })
                .flatten()
                .collect::<Vec<_>>();

            col.typed::<Int32Type>()
                .write_batch(&fragment_ordinals, None, None)?;
            col.close()?;
        }

        if let Some(mut col) = rg.next_column()? {
            let fragment_charge = features
                .iter()
                .flat_map(|f| {
                    f.fragments
                        .as_ref()
                        .map(|fragments| fragments.charges.iter().copied())
                })
                .flatten()
                .collect::<Vec<i32>>();

            col.typed::<Int32Type>()
                .write_batch(&fragment_charge, None, None)?;
            col.close()?;
        }

        if let Some(mut col) = rg.next_column()? {
            let fragment_mz_experimental = features
                .iter()
                .flat_map(|f| {
                    f.fragments
                        .as_ref()
                        .map(|fragments| fragments.mz_experimental.iter().copied())
                })
                .flatten()
                .collect::<Vec<_>>();

            col.typed::<FloatType>()
                .write_batch(&fragment_mz_experimental, None, None)?;
            col.close()?;
        }

        if let Some(mut col) = rg.next_column()? {
            let fragment_mz_calculated = features
                .iter()
                .flat_map(|f| {
                    f.fragments
                        .as_ref()
                        .map(|fragments| fragments.mz_calculated.iter().copied())
                })
                .flatten()
                .collect::<Vec<_>>();

            col.typed::<FloatType>()
                .write_batch(&fragment_mz_calculated, None, None)?;
            col.close()?;
        }

        if let Some(mut col) = rg.next_column()? {
            let neutral_losses = features
                .iter()
                .flat_map(|f| {
                    f.fragments
                        .as_ref()
                        .map(|fragments| fragments.neutral_losses.iter().copied())
                })
                .flatten()
                .collect::<Vec<_>>();

            col.typed::<FloatType>()
                .write_batch(&neutral_losses, None, None)?;
            col.close()?;
        }

        if let Some(mut col) = rg.next_column()? {
            let fragment_intensity = features
                .iter()
                .flat_map(|f| {
                    f.fragments
                        .as_ref()
                        .map(|fragments| fragments.intensities.iter().copied())
                })
                .flatten()
                .collect::<Vec<_>>();

            col.typed::<FloatType>()
                .write_batch(&fragment_intensity, None, None)?;
            col.close()?;
        }

        rg.close()?;
    }

    writer.into_inner()
}
