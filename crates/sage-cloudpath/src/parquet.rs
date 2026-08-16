//! Use the low-level `parquet` file writer API to serialize Sage results
//!
//! Modifying the file formats here requires some digging into documentation
//! about Dremel definition and repetition levels and the Parquet file format
//! https://akshays-blog.medium.com/wrapping-head-around-repetition-and-definition-levels-in-dremel-powering-bigquery-c1a33c9695da
//! https://blog.twitter.com/engineering/en_us/a/2013/dremel-made-simple-with-parquet
//! https://github.com/apache/parquet-format/blob/master/LogicalTypes.md

#![cfg(feature = "parquet")]

use std::collections::HashMap;
use std::hash::BuildHasher;

use parquet::data_type::{BoolType, ByteArray, FloatType, Int64Type};
use parquet::file::writer::SerializedColumnWriter;
use parquet::{
    basic::ZstdLevel,
    data_type::{ByteArrayType, DataType, Int32Type},
    file::{properties::WriterProperties, writer::SerializedFileWriter},
    schema::types::Type,
};
use sage_core::database::IndexedDatabase;
use sage_core::ion_series::Kind;
use sage_core::lfq::{Peak, PrecursorId};
use sage_core::scoring::Feature;
use sage_core::tmt::TmtQuant;

pub struct PtmSiteRecord {
    pub psm_id: i64,
    pub filename: String,
    pub scannr: String,
    pub peptide: String,
    pub proteins: String,
    pub charge: i32,
    pub spectrum_q: f32,
    pub peptide_q: f32,
    pub modification: String,
    pub modification_mass: f32,
    pub position: i32,
    pub residue: String,
    pub localization_probability: f32,
    pub delta_localization_score: f32,
    pub candidate_sites: i32,
    pub site_determining_ions_matched: i32,
    pub site_determining_ions_total: i32,
    pub site_probabilities: String,
}

pub struct ProteinSiteRecord {
    pub protein: String,
    pub peptide: String,
    pub residue: String,
    pub position_in_peptide: i32,
    pub modification: String,
    pub modification_mass: f32,
    pub num_psms: i32,
    pub best_localization_probability: f32,
    pub best_delta_localization_score: f32,
    pub best_spectrum_q: f32,
}

fn ptm_site_schema() -> parquet::errors::Result<Type> {
    parquet::schema::parser::parse_message_type(
        r#"
        message schema {
            required int64 psm_id;
            required byte_array filename (utf8);
            required byte_array scannr (utf8);
            required byte_array peptide (utf8);
            required byte_array proteins (utf8);
            required int32 charge;
            required float spectrum_q;
            required float peptide_q;
            required byte_array modification (utf8);
            required float modification_mass;
            required int32 position;
            required byte_array residue (utf8);
            required float localization_probability;
            required float delta_localization_score;
            required int32 candidate_sites;
            required int32 site_determining_ions_matched;
            required int32 site_determining_ions_total;
            required byte_array site_probabilities (utf8);
        }
        "#,
    )
}

fn protein_site_schema() -> parquet::errors::Result<Type> {
    parquet::schema::parser::parse_message_type(
        r#"
        message schema {
            required byte_array protein (utf8);
            required byte_array peptide (utf8);
            required byte_array residue (utf8);
            required int32 position_in_peptide;
            required byte_array modification (utf8);
            required float modification_mass;
            required int32 num_psms;
            required float best_localization_probability;
            required float best_delta_localization_score;
            required float best_spectrum_q;
        }
        "#,
    )
}

macro_rules! write_required_column {
    ($row_group:expr, $values:expr, $ty:ident) => {
        if let Some(mut column) = $row_group.next_column()? {
            column.typed::<$ty>().write_batch(&$values, None, None)?;
            column.close()?;
        }
    };
}

pub fn serialize_ptm_sites(records: &[PtmSiteRecord]) -> parquet::errors::Result<Vec<u8>> {
    let schema = ptm_site_schema()?;
    let options = WriterProperties::builder()
        .set_compression(parquet::basic::Compression::ZSTD(ZstdLevel::try_new(3)?))
        .build();
    let mut writer = SerializedFileWriter::new(Vec::new(), schema.into(), options.into())?;

    for records in records.chunks(65536) {
        let mut rg = writer.next_row_group()?;
        write_required_column!(
            rg,
            records.iter().map(|r| r.psm_id).collect::<Vec<_>>(),
            Int64Type
        );
        write_required_column!(
            rg,
            records
                .iter()
                .map(|r| r.filename.as_str().into())
                .collect::<Vec<ByteArray>>(),
            ByteArrayType
        );
        write_required_column!(
            rg,
            records
                .iter()
                .map(|r| r.scannr.as_str().into())
                .collect::<Vec<ByteArray>>(),
            ByteArrayType
        );
        write_required_column!(
            rg,
            records
                .iter()
                .map(|r| r.peptide.as_str().into())
                .collect::<Vec<ByteArray>>(),
            ByteArrayType
        );
        write_required_column!(
            rg,
            records
                .iter()
                .map(|r| r.proteins.as_str().into())
                .collect::<Vec<ByteArray>>(),
            ByteArrayType
        );
        write_required_column!(
            rg,
            records.iter().map(|r| r.charge).collect::<Vec<_>>(),
            Int32Type
        );
        write_required_column!(
            rg,
            records.iter().map(|r| r.spectrum_q).collect::<Vec<_>>(),
            FloatType
        );
        write_required_column!(
            rg,
            records.iter().map(|r| r.peptide_q).collect::<Vec<_>>(),
            FloatType
        );
        write_required_column!(
            rg,
            records
                .iter()
                .map(|r| r.modification.as_str().into())
                .collect::<Vec<ByteArray>>(),
            ByteArrayType
        );
        write_required_column!(
            rg,
            records
                .iter()
                .map(|r| r.modification_mass)
                .collect::<Vec<_>>(),
            FloatType
        );
        write_required_column!(
            rg,
            records.iter().map(|r| r.position).collect::<Vec<_>>(),
            Int32Type
        );
        write_required_column!(
            rg,
            records
                .iter()
                .map(|r| r.residue.as_str().into())
                .collect::<Vec<ByteArray>>(),
            ByteArrayType
        );
        write_required_column!(
            rg,
            records
                .iter()
                .map(|r| r.localization_probability)
                .collect::<Vec<_>>(),
            FloatType
        );
        write_required_column!(
            rg,
            records
                .iter()
                .map(|r| r.delta_localization_score)
                .collect::<Vec<_>>(),
            FloatType
        );
        write_required_column!(
            rg,
            records
                .iter()
                .map(|r| r.candidate_sites)
                .collect::<Vec<_>>(),
            Int32Type
        );
        write_required_column!(
            rg,
            records
                .iter()
                .map(|r| r.site_determining_ions_matched)
                .collect::<Vec<_>>(),
            Int32Type
        );
        write_required_column!(
            rg,
            records
                .iter()
                .map(|r| r.site_determining_ions_total)
                .collect::<Vec<_>>(),
            Int32Type
        );
        write_required_column!(
            rg,
            records
                .iter()
                .map(|r| r.site_probabilities.as_str().into())
                .collect::<Vec<ByteArray>>(),
            ByteArrayType
        );
        rg.close()?;
    }

    writer.into_inner().map(|bytes| bytes.to_vec())
}

pub fn serialize_protein_sites(records: &[ProteinSiteRecord]) -> parquet::errors::Result<Vec<u8>> {
    let schema = protein_site_schema()?;
    let options = WriterProperties::builder()
        .set_compression(parquet::basic::Compression::ZSTD(ZstdLevel::try_new(3)?))
        .build();
    let mut writer = SerializedFileWriter::new(Vec::new(), schema.into(), options.into())?;

    for records in records.chunks(65536) {
        let mut rg = writer.next_row_group()?;
        write_required_column!(
            rg,
            records
                .iter()
                .map(|r| r.protein.as_str().into())
                .collect::<Vec<ByteArray>>(),
            ByteArrayType
        );
        write_required_column!(
            rg,
            records
                .iter()
                .map(|r| r.peptide.as_str().into())
                .collect::<Vec<ByteArray>>(),
            ByteArrayType
        );
        write_required_column!(
            rg,
            records
                .iter()
                .map(|r| r.residue.as_str().into())
                .collect::<Vec<ByteArray>>(),
            ByteArrayType
        );
        write_required_column!(
            rg,
            records
                .iter()
                .map(|r| r.position_in_peptide)
                .collect::<Vec<_>>(),
            Int32Type
        );
        write_required_column!(
            rg,
            records
                .iter()
                .map(|r| r.modification.as_str().into())
                .collect::<Vec<ByteArray>>(),
            ByteArrayType
        );
        write_required_column!(
            rg,
            records
                .iter()
                .map(|r| r.modification_mass)
                .collect::<Vec<_>>(),
            FloatType
        );
        write_required_column!(
            rg,
            records.iter().map(|r| r.num_psms).collect::<Vec<_>>(),
            Int32Type
        );
        write_required_column!(
            rg,
            records
                .iter()
                .map(|r| r.best_localization_probability)
                .collect::<Vec<_>>(),
            FloatType
        );
        write_required_column!(
            rg,
            records
                .iter()
                .map(|r| r.best_delta_localization_score)
                .collect::<Vec<_>>(),
            FloatType
        );
        write_required_column!(
            rg,
            records
                .iter()
                .map(|r| r.best_spectrum_q)
                .collect::<Vec<_>>(),
            FloatType
        );
        rg.close()?;
    }

    writer.into_inner().map(|bytes| bytes.to_vec())
}

pub fn build_schema() -> Result<Type, parquet::errors::ParquetError> {
    let msg = r#"
        message schema {
            required int64 psm_id;
            required byte_array filename (utf8);
            required byte_array scannr (utf8);
            required byte_array peptide (utf8);
            required byte_array ambiguity_sequence (utf8);
            required float mass_shift;
            required byte_array stripped_peptide (utf8);
            required byte_array proteins (utf8);
            required byte_array protein_groups (utf8);
            required int32 num_proteins;
            required int32 num_protein_groups;
            required int32 rank;
            required boolean is_decoy;
            required float expmass;
            required float calcmass;
            required int32 charge;
            required int32 peptide_len;
            required int32 missed_cleavages;
            required boolean semi_enzymatic;
            required float ms2_intensity;
            required float isotope_error;
            required float precursor_ppm;
            required float fragment_ppm;
            required float hyperscore;
            required float delta_next;
            required float delta_best;
            required float rt;
            required float aligned_rt;
            required float predicted_rt;
            required float delta_rt_model;
            required float ion_mobility;
            required float predicted_mobility;
            required float delta_mobility;
            required int32 matched_peaks;
            required int32 longest_b;
            required int32 longest_y;
            required float longest_y_pct;
            required float matched_intensity_pct;
            required int32 scored_candidates;
            required float poisson;
            required float sage_discriminant_score;
            required float posterior_error;
            required float spectrum_q;
            required float peptide_q;
            required float protein_q;
            required float protein_group_q;
            optional group reporter_ion_intensity (LIST) {
                repeated group list {
                    optional float element;
                }
            }
        }
    "#;
    parquet::schema::parser::parse_message_type(msg)
}

/// Caller must guarantee that `reporter_ions` is not an empty slice
fn write_reporter_ions(
    mut column: SerializedColumnWriter,
    features: &[Feature],
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
    features: &[Feature],
    reporter_ions: &[TmtQuant],
    filenames: &[String],
    database: &IndexedDatabase,
) -> Result<Vec<u8>, parquet::errors::ParquetError> {
    let schema = build_schema()?;

    let options = WriterProperties::builder()
        .set_compression(parquet::basic::Compression::ZSTD(ZstdLevel::try_new(3)?))
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

        write_col!(|f: &Feature| f.psm_id as i64, Int64Type);
        write_col!(
            |f: &Feature| filenames[f.file_id].as_str().into(),
            ByteArrayType
        );
        write_col!(|f: &Feature| f.spec_id.as_str().into(), ByteArrayType);
        write_col!(
            |f: &Feature| database[f.peptide_idx].to_string().as_bytes().into(),
            ByteArrayType
        );
        write_col!(
            |f: &Feature| f.ambiguity_sequence.as_str().into(),
            ByteArrayType
        );
        write_col!(mass_shift, FloatType);
        write_col!(
            |f: &Feature| database[f.peptide_idx].sequence.as_ref().into(),
            ByteArrayType
        );
        write_col!(
            |f: &Feature| database[f.peptide_idx]
                .proteins(&database.decoy_tag, database.generate_decoys)
                .as_str()
                .into(),
            ByteArrayType
        );
        write_col!(
            |f: &Feature| f.protein_groups.as_deref().unwrap_or("").into(),
            ByteArrayType
        );
        write_col!(
            |f: &Feature| database[f.peptide_idx].proteins.len() as i32,
            Int32Type
        );
        write_col!(num_protein_groups, Int32Type);
        write_col!(rank, Int32Type);
        write_col!(|f: &Feature| f.label == -1, BoolType);
        write_col!(expmass, FloatType);
        write_col!(calcmass, FloatType);
        write_col!(charge, Int32Type);
        write_col!(peptide_len, Int32Type);
        write_col!(missed_cleavages, Int32Type);
        write_col!(
            |f: &Feature| database[f.peptide_idx].semi_enzymatic,
            BoolType
        );
        write_col!(ms2_intensity, FloatType);
        write_col!(isotope_error, FloatType);
        write_col!(delta_mass, FloatType);
        write_col!(average_ppm, FloatType);
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
            required float fragment_intensity;
        }
    "#;

    parquet::schema::parser::parse_message_type(msg)
}

pub fn serialize_matched_fragments(
    features: &[Feature],
) -> Result<Vec<u8>, parquet::errors::ParquetError> {
    let schema = build_matched_fragment_schema()?;

    let options = WriterProperties::builder()
        .set_compression(parquet::basic::Compression::ZSTD(ZstdLevel::try_new(3)?))
        .build();

    let buf = Vec::new();

    let mut writer = SerializedFileWriter::new(buf, schema.into(), options.into())?;

    for features in features.chunks(65536) {
        let mut rg = writer.next_row_group()?;

        if let Some(mut col) = rg.next_column()? {
            let psm_ids = features
                .iter()
                .flat_map(|f| {
                    std::iter::repeat(f.psm_id as i64).take(
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

pub fn build_lfq_schema() -> parquet::errors::Result<Type> {
    let msg = r#"
        message schema {
            required byte_array peptide (utf8);
            required byte_array stripped_peptide (utf8);
            optional int32 charge;
            required byte_array proteins (utf8);
            required boolean is_decoy;
            required float q_value;
            required byte_array filename (utf8);
            required float intensity;
        }
    "#;
    parquet::schema::parser::parse_message_type(msg)
}

pub fn serialize_lfq<H: BuildHasher>(
    areas: &HashMap<(PrecursorId, bool), (Peak, Vec<f64>), H>,
    filenames: &[String],
    database: &IndexedDatabase,
) -> parquet::errors::Result<Vec<u8>> {
    let schema = build_lfq_schema()?;

    let options = WriterProperties::builder()
        .set_compression(parquet::basic::Compression::ZSTD(ZstdLevel::try_new(3)?))
        .build();

    let buf = Vec::new();
    let mut writer = SerializedFileWriter::new(buf, schema.into(), options.into())?;
    let mut rg = writer.next_row_group()?;

    if let Some(mut col) = rg.next_column()? {
        let values = areas
            .iter()
            .flat_map(|((id, _), _)| {
                let peptide_idx = match id {
                    PrecursorId::Combined(x) | PrecursorId::Charged((x, _)) => x,
                };
                let val = database[*peptide_idx].to_string().as_bytes().into();
                std::iter::repeat(val).take(filenames.len())
            })
            .collect::<Vec<_>>();

        col.typed::<ByteArrayType>()
            .write_batch(&values, None, None)?;
        col.close()?;
    }

    if let Some(mut col) = rg.next_column()? {
        let values = areas
            .iter()
            .flat_map(|((id, _), _)| {
                let peptide_idx = match id {
                    PrecursorId::Combined(x) | PrecursorId::Charged((x, _)) => x,
                };
                let val = database[*peptide_idx].sequence.as_ref().into();
                std::iter::repeat(val).take(filenames.len())
            })
            .collect::<Vec<_>>();

        col.typed::<ByteArrayType>()
            .write_batch(&values, None, None)?;
        col.close()?;
    }

    if let Some(mut col) = rg.next_column()? {
        let mut values = Vec::with_capacity(areas.len() * filenames.len());
        let mut def_levels = Vec::with_capacity(areas.len() * filenames.len());

        for ((id, _), _) in areas.iter() {
            match id {
                PrecursorId::Combined(_) => {
                    def_levels.extend(std::iter::repeat(0).take(filenames.len()));
                }
                PrecursorId::Charged((_, charge)) => {
                    values.extend(std::iter::repeat(*charge as i32).take(filenames.len()));
                    def_levels.extend(std::iter::repeat(1).take(filenames.len()));
                }
            }
        }

        col.typed::<Int32Type>()
            .write_batch(&values, Some(&def_levels), None)?;
        col.close()?;
    }

    if let Some(mut col) = rg.next_column()? {
        let values = areas
            .iter()
            .flat_map(|((id, _), _)| {
                let peptide_idx = match id {
                    PrecursorId::Combined(x) | PrecursorId::Charged((x, _)) => x,
                };
                let val = database[*peptide_idx]
                    .proteins(&database.decoy_tag, database.generate_decoys)
                    .as_str()
                    .into();
                std::iter::repeat(val).take(filenames.len())
            })
            .collect::<Vec<_>>();

        col.typed::<ByteArrayType>()
            .write_batch(&values, None, None)?;
        col.close()?;
    }

    if let Some(mut col) = rg.next_column()? {
        let values = areas
            .iter()
            .flat_map(|((_, decoy), _)| std::iter::repeat(*decoy).take(filenames.len()))
            .collect::<Vec<_>>();

        col.typed::<BoolType>().write_batch(&values, None, None)?;
        col.close()?;
    }

    if let Some(mut col) = rg.next_column()? {
        let values = areas
            .iter()
            .flat_map(|(_, (peak, _))| std::iter::repeat(peak.q_value).take(filenames.len()))
            .collect::<Vec<_>>();

        col.typed::<FloatType>().write_batch(&values, None, None)?;
        col.close()?;
    }

    if let Some(mut col) = rg.next_column()? {
        let values = areas
            .iter()
            .flat_map(|(_, (_, values))| {
                (0..values.len()).map(|idx| filenames[idx].as_bytes().into())
            })
            .collect::<Vec<_>>();

        col.typed::<ByteArrayType>()
            .write_batch(&values, None, None)?;

        col.close()?;
    }

    if let Some(mut col) = rg.next_column()? {
        let values = areas
            .iter()
            .flat_map(|(_, (_, values))| values.iter().copied().map(|v| v as f32))
            .collect::<Vec<_>>();

        col.typed::<FloatType>().write_batch(&values, None, None)?;
        col.close()?;
    }

    rg.close()?;
    writer.into_inner()
}

#[cfg(test)]
mod ptm_tests {
    use super::*;
    use parquet::file::reader::{FileReader, SerializedFileReader};

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
            18
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
            10
        );
    }
}
