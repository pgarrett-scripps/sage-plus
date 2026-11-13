use super::*;
use sage_core::database::PeptideIx;
use sage_core::in_source::InSourceAnnotations;

pub fn build_lfq_schema() -> parquet::errors::Result<Type> {
    parquet::schema::parser::parse_message_type(include_str!(
        "../../../../schemas/lfq.v9.parquet.schema"
    ))
}

/// Schema version of `lfq.parquet`: 9 without precursor labels, 10 with them.
/// Versions 9 and 10 are 7 and 8 plus `in_source_parent_peptide` and
/// `in_source_type`, which are always written and null when the annotation is
/// off, as every optional column before them.
fn lfq_schema_version(has_labels: bool) -> &'static str {
    if has_labels {
        "10"
    } else {
        "9"
    }
}

fn build_lfq_schema_version(has_labels: bool) -> parquet::errors::Result<Type> {
    parquet::schema::parser::parse_message_type(if has_labels {
        include_str!("../../../../schemas/lfq.v10.parquet.schema")
    } else {
        include_str!("../../../../schemas/lfq.v9.parquet.schema")
    })
}

pub fn serialize_lfq<H: BuildHasher>(
    areas: &HashMap<(PrecursorId, bool), QuantifiedPeak, H>,
    filenames: &[String],
    database: &IndexedDatabase,
    provenance: &[(String, String)],
    peak_center: &str,
    in_source: Option<&InSourceAnnotations>,
) -> parquet::errors::Result<Vec<u8>> {
    if let Some((_, quantified)) = areas.iter().find(|(_, quantified)| {
        quantified.intensities.len() != filenames.len()
            || quantified.ms2_confirmed.len() != filenames.len()
            || quantified.ms2_confirmed_strict.len() != filenames.len()
            || quantified.file_evidence.len() != filenames.len()
    }) {
        return Err(ParquetError::General(format!(
            "LFQ row has {} intensities, {} MS2 values, {} strict MS2 values and {} file evidence values for {} files",
            quantified.intensities.len(),
            quantified.ms2_confirmed.len(),
            quantified.ms2_confirmed_strict.len(),
            quantified.file_evidence.len(),
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
                    (peptide.label_channel_str() == Some(reference)).then(|| {
                        (
                            (peptide.label_group(), charge, *decoy),
                            quantified.intensities.as_slice(),
                        )
                    })
                })
                .collect::<HashMap<_, _>>()
        })
        .unwrap_or_default();

    let options = writer_properties(
        vec![
            KeyValue::new("sage.schema.name".into(), Some("lfq".into())),
            KeyValue::new(
                "sage.lfq.q_value_scope".into(),
                Some("precursor_across_files".into()),
            ),
            KeyValue::new(
                "sage.lfq.file_score_status".into(),
                Some("experimental_uncalibrated".into()),
            ),
            KeyValue::new(
                "sage.lfq.extraction_q_value_scope".into(),
                Some("precursor_file_all_rows".into()),
            ),
            KeyValue::new("sage.lfq.peak_center".into(), Some(peak_center.into())),
            KeyValue::new(
                "sage.schema.version".into(),
                Some(lfq_schema_version(has_labels).into()),
            ),
        ],
        provenance,
    )?;

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
                        .label_channel_str()
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

    if let Some(mut col) = rg.next_column()? {
        let values = rows
            .iter()
            .flat_map(|(_, peak)| peak.ms2_confirmed_strict.iter().copied())
            .collect::<Vec<_>>();
        col.typed::<BoolType>().write_batch(&values, None, None)?;
        col.close()?;
    }
    for field in 0..3 {
        if let Some(mut col) = rg.next_column()? {
            let mut values = Vec::new();
            let mut levels = Vec::new();
            for evidence in rows.iter().flat_map(|(_, peak)| &peak.file_evidence) {
                if let Some(evidence) = evidence {
                    values.push(match field {
                        0 => evidence.score,
                        1 => evidence.spectral_angle,
                        _ => evidence.trace_cosine,
                    });
                    levels.push(1);
                } else {
                    levels.push(0);
                }
            }
            col.typed::<DoubleType>()
                .write_batch(&values, Some(&levels), None)?;
            col.close()?;
        }
    }
    if let Some(mut col) = rg.next_column()? {
        let mut values = Vec::new();
        let mut levels = Vec::new();
        for evidence in rows.iter().flat_map(|(_, peak)| &peak.file_evidence) {
            if let Some(evidence) = evidence {
                values.push(evidence.rt_shift_bins);
                levels.push(1);
            } else {
                levels.push(0);
            }
        }
        col.typed::<Int32Type>()
            .write_batch(&values, Some(&levels), None)?;
        col.close()?;
    }
    if let Some(mut col) = rg.next_column()? {
        let mut values = Vec::new();
        let mut levels = Vec::new();
        for evidence in rows.iter().flat_map(|(_, peak)| &peak.file_evidence) {
            if let Some(evidence) = evidence {
                values.push(evidence.transfer_candidate);
                levels.push(1);
            } else {
                levels.push(0);
            }
        }
        col.typed::<BoolType>()
            .write_batch(&values, Some(&levels), None)?;
        col.close()?;
    }
    if let Some(mut col) = rg.next_column()? {
        let mut values = Vec::new();
        let mut levels = Vec::new();
        for evidence in rows.iter().flat_map(|(_, peak)| &peak.file_evidence) {
            match evidence
                .as_ref()
                .and_then(|evidence| evidence.extraction_q_value)
            {
                Some(q) => {
                    values.push(q);
                    levels.push(1);
                }
                None => levels.push(0),
            }
        }
        col.typed::<FloatType>()
            .write_batch(&values, Some(&levels), None)?;
        col.close()?;
    }
    for field in 0..5 {
        if let Some(mut col) = rg.next_column()? {
            let mut values = Vec::new();
            let mut levels = Vec::new();
            for evidence in rows.iter().flat_map(|(_, peak)| &peak.file_evidence) {
                let value = evidence.as_ref().and_then(|evidence| match field {
                    0 => Some(evidence.apex_rt),
                    1 => Some(evidence.peak_start_rt),
                    2 => Some(evidence.peak_end_rt),
                    3 => evidence.fwhm,
                    _ => evidence.id_apex_offset,
                });
                match value.filter(|value| value.is_finite()) {
                    Some(value) => {
                        values.push(value);
                        levels.push(1);
                    }
                    None => levels.push(0),
                }
            }
            col.typed::<FloatType>()
                .write_batch(&values, Some(&levels), None)?;
            col.close()?;
        }
    }
    // A quantified row is judged on its own peak: its apex retention time in
    // this file against the parent's apex in the same file, not on a PSM's
    // retention time, so a row transferred between runs is judged too. Decoy
    // rows are never annotated. Both columns are always written and are null
    // when the annotation is off.
    let parent_apexes = in_source.map(|_| {
        let mut apexes: HashMap<(PeptideIx, usize), Vec<f32>> = HashMap::new();
        for ((id, decoy), peak) in areas.iter() {
            if *decoy {
                continue;
            }
            let peptide_idx = match id {
                PrecursorId::Combined(idx) | PrecursorId::Charged((idx, _)) => *idx,
            };
            for (file_id, evidence) in peak.file_evidence.iter().enumerate() {
                if let Some(apex) = evidence
                    .as_ref()
                    .map(|evidence| evidence.apex_rt)
                    .filter(|apex| apex.is_finite())
                {
                    apexes.entry((peptide_idx, file_id)).or_default().push(apex);
                }
            }
        }
        apexes
    });
    let parent_apexes = parent_apexes.as_ref();
    let verdicts = rows
        .iter()
        .flat_map(|((id, decoy), peak)| {
            let peptide_idx = match id {
                PrecursorId::Combined(idx) | PrecursorId::Charged((idx, _)) => *idx,
            };
            let decoy = *decoy;
            peak.file_evidence
                .iter()
                .enumerate()
                .map(move |(file_id, evidence)| {
                    let annotations = in_source?;
                    let apexes = parent_apexes?;
                    if decoy {
                        return None;
                    }
                    let apex = evidence
                        .as_ref()
                        .map(|evidence| evidence.apex_rt)
                        .filter(|apex| apex.is_finite())?;
                    let candidates = annotations.candidates(peptide_idx);
                    annotations.coeluting(file_id, apex, candidates, |parent| {
                        apexes.get(&(parent, file_id)).and_then(|parent_apexes| {
                            parent_apexes.iter().copied().min_by(|left, right| {
                                (left - apex).abs().total_cmp(&(right - apex).abs())
                            })
                        })
                    })
                })
        })
        .collect::<Vec<_>>();
    if let Some(mut col) = rg.next_column()? {
        let mut values = Vec::new();
        let mut levels = Vec::new();
        for verdict in &verdicts {
            match verdict {
                Some(annotation) => {
                    values.push(database[annotation.parent].to_string().into_bytes().into());
                    levels.push(1);
                }
                None => levels.push(0),
            }
        }
        col.typed::<ByteArrayType>()
            .write_batch(&values, Some(&levels), None)?;
        col.close()?;
    }
    if let Some(mut col) = rg.next_column()? {
        let mut values = Vec::new();
        let mut levels = Vec::new();
        for verdict in &verdicts {
            match verdict {
                Some(annotation) => {
                    values.push(annotation.kind.label().into_bytes().into());
                    levels.push(1);
                }
                None => levels.push(0),
            }
        }
        col.typed::<ByteArrayType>()
            .write_batch(&values, Some(&levels), None)?;
        col.close()?;
    }
    rg.close()?;
    writer.into_inner()
}
