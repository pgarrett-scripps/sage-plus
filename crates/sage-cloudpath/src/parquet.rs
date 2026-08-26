//! Use the low-level `parquet` file writer API to serialize Sage results
//!
//! Modifying the file formats here requires some digging into documentation
//! about Dremel definition and repetition levels and the Parquet file format
//! https://akshays-blog.medium.com/wrapping-head-around-repetition-and-definition-levels-in-dremel-powering-bigquery-c1a33c9695da
//! https://blog.twitter.com/engineering/en_us/a/2013/dremel-made-simple-with-parquet
//! https://github.com/apache/parquet-format/blob/master/LogicalTypes.md

#![cfg(feature = "parquet")]

use std::collections::HashMap;
use std::fs::File;
use std::hash::BuildHasher;
use std::path::Path;

use parquet::data_type::{BoolType, ByteArray, DoubleType, FloatType, Int64Type};
use parquet::errors::ParquetError;
use parquet::file::metadata::KeyValue;
use parquet::file::reader::{FileReader, SerializedFileReader};
use parquet::file::writer::SerializedColumnWriter;
use parquet::record::{Field, Row};
use parquet::{
    basic::ZstdLevel,
    data_type::{ByteArrayType, DataType, Int32Type},
    file::{properties::WriterProperties, writer::SerializedFileWriter},
    schema::types::Type,
};
use sage_core::cleavage::CustomCleavageLibrary;
use sage_core::database::IndexedDatabase;
use sage_core::ion_series::Kind;
use sage_core::lfq::{PrecursorId, QuantifiedPeak};
use sage_core::ptm_library::{PtmLibrary, PtmLibrarySite};
use sage_core::scoring::Feature;
use sage_core::spectral_library::{
    LibraryFragment, SpectralLibraryEntry, SpectralLibrarySettings, SpectralLibraryStrategy,
};
use sage_core::spectral_library_search::{DdaLibraryEntry, DdaLibraryIndex};
use sage_core::tmt::TmtQuant;

fn field_to_json(field: &Field) -> serde_json::Value {
    use serde_json::{Number, Value};

    match field {
        Field::Null => Value::Null,
        Field::Bool(value) => Value::Bool(*value),
        Field::Byte(value) => Value::Number(Number::from(*value)),
        Field::Short(value) => Value::Number(Number::from(*value)),
        Field::Int(value) => Value::Number(Number::from(*value)),
        Field::Long(value) => Value::Number(Number::from(*value)),
        Field::UByte(value) => Value::Number(Number::from(*value)),
        Field::UShort(value) => Value::Number(Number::from(*value)),
        Field::UInt(value) => Value::Number(Number::from(*value)),
        Field::ULong(value) => Value::Number(Number::from(*value)),
        Field::Float(value) => Number::from_f64(f64::from(*value))
            .map(Value::Number)
            .unwrap_or(Value::Null),
        Field::Double(value) => Number::from_f64(*value)
            .map(Value::Number)
            .unwrap_or(Value::Null),
        Field::Str(value) => Value::String(value.clone()),
        Field::Bytes(value) => Value::String(base64::encode(value.data())),
        Field::Date(value) => Value::Number(Number::from(*value)),
        Field::TimestampMillis(value) | Field::TimestampMicros(value) => {
            Value::Number(Number::from(*value))
        }
        // Sage's analytical schemas currently use only the scalar types above.
        // Preserve any future logical/nested values as their Parquet display form.
        other => Value::String(other.to_string()),
    }
}

/// Scan a Parquet file as typed JSON objects without materializing the full file.
/// The predicate is evaluated before `limit` is applied.
pub fn scan_json_rows<F>(
    path: &Path,
    scan_limit: usize,
    limit: usize,
    mut predicate: F,
) -> parquet::errors::Result<(Vec<serde_json::Value>, usize, bool)>
where
    F: FnMut(&serde_json::Map<String, serde_json::Value>) -> bool,
{
    let file = File::open(path).map_err(|error| ParquetError::External(Box::new(error)))?;
    let reader = SerializedFileReader::new(file)?;
    let mut rows = Vec::new();
    let mut scanned_rows = 0usize;
    let mut truncated = false;

    for row in reader.get_row_iter(None)?.take(scan_limit) {
        scanned_rows += 1;
        let row = row?;
        let object = row
            .get_column_iter()
            .map(|(name, field)| (name.clone(), field_to_json(field)))
            .collect::<serde_json::Map<_, _>>();
        if !predicate(&object) {
            continue;
        }
        if rows.len() == limit {
            truncated = true;
            break;
        }
        rows.push(serde_json::Value::Object(object));
    }
    if scanned_rows == scan_limit {
        truncated = true;
    }
    Ok((rows, scanned_rows, truncated))
}

/// Read a compact PTM site library. Required columns are `protein`,
/// `position` (one-based), `residue`, and `modification`. Additional evidence
/// columns are intentionally ignored by database construction.
pub fn deserialize_ptm_library(bytes: Vec<u8>) -> parquet::errors::Result<PtmLibrary> {
    use parquet::errors::ParquetError;
    use parquet::file::reader::{FileReader, SerializedFileReader};
    use parquet::record::Field;
    use std::sync::Arc;

    fn text(field: &Field, column: &str) -> parquet::errors::Result<String> {
        match field {
            Field::Str(value) => Ok(value.clone()),
            _ => Err(ParquetError::General(format!(
                "PTM library column `{column}` must be a UTF-8 string"
            ))),
        }
    }

    fn position(field: &Field) -> parquet::errors::Result<u32> {
        let value = match field {
            Field::Int(value) => i64::from(*value),
            Field::Long(value) => *value,
            Field::UInt(value) => i64::from(*value),
            _ => {
                return Err(ParquetError::General(
                    "PTM library column `position` must be an integer".into(),
                ))
            }
        };
        u32::try_from(value)
            .ok()
            .filter(|value| *value > 0)
            .map(|value| value - 1)
            .ok_or_else(|| {
                ParquetError::General(
                    "PTM library positions must be positive, one-based integers".into(),
                )
            })
    }

    let reader = SerializedFileReader::new(bytes::Bytes::from(bytes))?;
    let mut sites = Vec::new();
    for (row_index, row) in reader.get_row_iter(None)?.enumerate() {
        let row = row?;
        let columns = row
            .get_column_iter()
            .map(|(name, field)| (name.as_str(), field))
            .collect::<HashMap<_, _>>();
        let required = |name: &str| {
            columns.get(name).copied().ok_or_else(|| {
                ParquetError::General(format!("PTM library is missing required column `{name}`"))
            })
        };
        let protein = text(required("protein")?, "protein")?;
        let modification = text(required("modification")?, "modification")?;
        let residue = text(required("residue")?, "residue")?;
        let residue = residue.as_bytes();
        if protein.is_empty() || modification.is_empty() || residue.len() != 1 {
            return Err(ParquetError::General(format!(
                "invalid PTM library values in row {}",
                row_index + 1
            )));
        }
        sites.push(PtmLibrarySite {
            protein: Arc::from(protein),
            position: position(required("position")?)?,
            residue: residue[0],
            modification: Arc::from(modification),
        });
    }
    Ok(PtmLibrary::new(sites))
}

fn ptm_library_schema() -> parquet::errors::Result<Type> {
    parquet::schema::parser::parse_message_type(
        r#"
        message schema {
            required byte_array protein (utf8);
            required int32 position;
            required byte_array residue (utf8);
            required byte_array modification (utf8);
        }
        "#,
    )
}

pub fn serialize_ptm_library(sites: &[PtmLibrarySite]) -> parquet::errors::Result<Vec<u8>> {
    if sites.iter().any(|site| site.position >= i32::MAX as u32) {
        return Err(parquet::errors::ParquetError::General(
            "PTM library position exceeds the int32 schema limit".into(),
        ));
    }
    let schema = ptm_library_schema()?;
    let options = WriterProperties::builder()
        .set_compression(parquet::basic::Compression::ZSTD(ZstdLevel::try_new(3)?))
        .build();
    let mut writer = SerializedFileWriter::new(Vec::new(), schema.into(), options.into())?;
    for sites in sites.chunks(65536) {
        let mut rg = writer.next_row_group()?;
        if let Some(mut column) = rg.next_column()? {
            let values = sites
                .iter()
                .map(|site| site.protein.as_ref().into())
                .collect::<Vec<ByteArray>>();
            column
                .typed::<ByteArrayType>()
                .write_batch(&values, None, None)?;
            column.close()?;
        }
        if let Some(mut column) = rg.next_column()? {
            let values = sites
                .iter()
                .map(|site| (site.position + 1) as i32)
                .collect::<Vec<_>>();
            column
                .typed::<Int32Type>()
                .write_batch(&values, None, None)?;
            column.close()?;
        }
        if let Some(mut column) = rg.next_column()? {
            let values = sites
                .iter()
                .map(|site| vec![site.residue].into())
                .collect::<Vec<ByteArray>>();
            column
                .typed::<ByteArrayType>()
                .write_batch(&values, None, None)?;
            column.close()?;
        }
        if let Some(mut column) = rg.next_column()? {
            let values = sites
                .iter()
                .map(|site| site.modification.as_ref().into())
                .collect::<Vec<ByteArray>>();
            column
                .typed::<ByteArrayType>()
                .write_batch(&values, None, None)?;
            column.close()?;
        }
        rg.close()?;
    }
    writer.into_inner().map(|bytes| bytes.to_vec())
}

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
    pub target_decoy_score: f32,
    pub localization_q_value: f32,
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
    pub best_localization_q_value: f32,
    pub best_spectrum_q: f32,
}

/// Read a protein-specific custom cleavage library from Parquet bytes.
/// Required columns are UTF-8 `protein` and integer `position`; optional
/// `context` is UTF-8 and may be null.
pub fn deserialize_custom_cleavage_sites(
    bytes: Vec<u8>,
) -> parquet::errors::Result<CustomCleavageLibrary> {
    let reader = SerializedFileReader::new(bytes::Bytes::from(bytes))?;
    let columns = reader
        .metadata()
        .file_metadata()
        .schema_descr()
        .columns()
        .iter()
        .map(|column| column.name())
        .collect::<Vec<_>>();
    for required in ["protein", "position"] {
        if !columns.contains(&required) {
            return Err(ParquetError::General(format!(
                "custom cleavage-site Parquet is missing required `{required}` column"
            )));
        }
    }

    let mut records = Vec::new();
    for (index, row) in reader.get_row_iter(None)?.enumerate() {
        let row_number = index + 1;
        let row = row?;
        let protein = match row_field(&row, "protein") {
            Some(Field::Str(protein)) => protein.clone(),
            Some(field) => {
                return Err(ParquetError::General(format!(
                    "custom cleavage-site Parquet row {row_number} has non-string `protein` value `{field}`"
                )));
            }
            None => {
                return Err(ParquetError::General(format!(
                    "custom cleavage-site Parquet row {row_number} is missing `protein`"
                )));
            }
        };
        let position = parquet_position(row_field(&row, "position"), row_number)?;
        let context = match row_field(&row, "context") {
            None | Some(Field::Null) => None,
            Some(Field::Str(context)) => Some(context.clone()),
            Some(field) => {
                return Err(ParquetError::General(format!(
                    "custom cleavage-site Parquet row {row_number} has non-string `context` value `{field}`"
                )));
            }
        };
        records.push((protein, position, context));
    }

    CustomCleavageLibrary::from_records(records)
        .map_err(|error| ParquetError::General(error.to_string()))
}

fn row_field<'a>(row: &'a Row, name: &str) -> Option<&'a Field> {
    row.get_column_iter()
        .find_map(|(column, field)| (column == name).then_some(field))
}

fn parquet_position(field: Option<&Field>, row: usize) -> parquet::errors::Result<usize> {
    let value = match field {
        Some(Field::Byte(value)) => i64::from(*value),
        Some(Field::Short(value)) => i64::from(*value),
        Some(Field::Int(value)) => i64::from(*value),
        Some(Field::Long(value)) => *value,
        Some(Field::UByte(value)) => return Ok(usize::from(*value)),
        Some(Field::UShort(value)) => return Ok(usize::from(*value)),
        Some(Field::UInt(value)) => {
            return usize::try_from(*value).map_err(|_| {
                ParquetError::General(format!(
                "custom cleavage-site Parquet row {row} has `position` outside the supported range"
            ))
            })
        }
        Some(Field::ULong(value)) => {
            return usize::try_from(*value).map_err(|_| {
                ParquetError::General(format!(
                    "custom cleavage-site Parquet row {row} has `position` outside the supported range"
                ))
            });
        }
        Some(field) => {
            return Err(ParquetError::General(format!(
                "custom cleavage-site Parquet row {row} has non-integer `position` value `{field}`"
            )));
        }
        None => {
            return Err(ParquetError::General(format!(
                "custom cleavage-site Parquet row {row} is missing `position`"
            )));
        }
    };
    usize::try_from(value).map_err(|_| {
        ParquetError::General(format!(
            "custom cleavage-site Parquet row {row} has negative `position` {value}"
        ))
    })
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
            required float target_decoy_score;
            required float localization_q_value;
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
            required float best_localization_q_value;
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
                .map(|r| r.target_decoy_score)
                .collect::<Vec<_>>(),
            FloatType
        );
        write_required_column!(
            rg,
            records
                .iter()
                .map(|r| r.localization_q_value)
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
                .map(|r| r.best_localization_q_value)
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
    parquet::schema::parser::parse_message_type(include_str!(
        "../../../schemas/results.sage.v1.parquet.schema"
    ))
}

fn build_results_schema(has_labels: bool) -> Result<Type, parquet::errors::ParquetError> {
    parquet::schema::parser::parse_message_type(if has_labels {
        include_str!("../../../schemas/results.sage.v2.parquet.schema")
    } else {
        include_str!("../../../schemas/results.sage.v1.parquet.schema")
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
        write_col!(spectral_angle, FloatType);
        write_col!(explained_library_intensity, FloatType);
        write_col!(explained_query_intensity, FloatType);
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

pub fn build_spectral_library_schema() -> parquet::errors::Result<Type> {
    parquet::schema::parser::parse_message_type(include_str!(
        "../../../schemas/spectral_library.sage.v1.parquet.schema"
    ))
}

fn build_spectral_library_schema_version(has_labels: bool) -> parquet::errors::Result<Type> {
    parquet::schema::parser::parse_message_type(if has_labels {
        include_str!("../../../schemas/spectral_library.sage.v2.parquet.schema")
    } else {
        include_str!("../../../schemas/spectral_library.sage.v1.parquet.schema")
    })
}

/// Read Sage's canonical long-form empirical library into DDA search entries.
///
/// Sage exports contain target entries only. `is_decoy` therefore remains
/// false until a future schema defines portable decoy provenance.
pub fn deserialize_spectral_library(
    bytes: Vec<u8>,
) -> parquet::errors::Result<Vec<DdaLibraryEntry>> {
    fn required<'a>(
        row: &'a Row,
        name: &str,
        row_number: usize,
    ) -> parquet::errors::Result<&'a Field> {
        row_field(row, name).ok_or_else(|| {
            ParquetError::General(format!(
                "spectral-library row {row_number} is missing `{name}`"
            ))
        })
    }

    fn text(field: &Field, name: &str, row_number: usize) -> parquet::errors::Result<String> {
        match field {
            Field::Str(value) => Ok(value.clone()),
            _ => Err(ParquetError::General(format!(
                "spectral-library row {row_number} has non-string `{name}`"
            ))),
        }
    }

    fn float(field: &Field, name: &str, row_number: usize) -> parquet::errors::Result<f32> {
        match field {
            Field::Float(value) => Ok(*value),
            Field::Double(value) => Ok(*value as f32),
            _ => Err(ParquetError::General(format!(
                "spectral-library row {row_number} has non-floating-point `{name}`"
            ))),
        }
    }

    fn integer(field: &Field, name: &str, row_number: usize) -> parquet::errors::Result<i64> {
        match field {
            Field::Byte(value) => Ok(i64::from(*value)),
            Field::Short(value) => Ok(i64::from(*value)),
            Field::Int(value) => Ok(i64::from(*value)),
            Field::Long(value) => Ok(*value),
            Field::UByte(value) => Ok(i64::from(*value)),
            Field::UShort(value) => Ok(i64::from(*value)),
            Field::UInt(value) => Ok(i64::from(*value)),
            Field::ULong(value) => i64::try_from(*value).map_err(|_| {
                ParquetError::General(format!(
                    "spectral-library row {row_number} has out-of-range `{name}`"
                ))
            }),
            _ => Err(ParquetError::General(format!(
                "spectral-library row {row_number} has non-integer `{name}`"
            ))),
        }
    }

    fn charge(field: &Field, name: &str, row_number: usize) -> parquet::errors::Result<u8> {
        u8::try_from(integer(field, name, row_number)?)
            .ok()
            .filter(|charge| *charge > 0)
            .ok_or_else(|| {
                ParquetError::General(format!(
                    "spectral-library row {row_number} has invalid `{name}`"
                ))
            })
    }

    fn fragment_kind(field: &Field, row_number: usize) -> parquet::errors::Result<Kind> {
        match text(field, "fragment_type", row_number)?.as_str() {
            "a" => Ok(Kind::A),
            "b" => Ok(Kind::B),
            "c" => Ok(Kind::C),
            "x" => Ok(Kind::X),
            "y" => Ok(Kind::Y),
            "z" => Ok(Kind::Z),
            kind => Err(ParquetError::General(format!(
                "spectral-library row {row_number} has unsupported fragment type `{kind}`"
            ))),
        }
    }

    let reader = SerializedFileReader::new(bytes::Bytes::from(bytes))?;
    let metadata = reader.metadata().file_metadata().key_value_metadata();
    let has_metadata = |key: &str, value: &str| {
        metadata.is_some_and(|entries| {
            entries
                .iter()
                .any(|entry| entry.key == key && entry.value.as_deref() == Some(value))
        })
    };
    let schema_version = if has_metadata("sage.schema.version", "2") {
        2
    } else if has_metadata("sage.schema.version", "1") {
        1
    } else {
        0
    };
    if !has_metadata("sage.schema.name", "spectral_library") || schema_version == 0 {
        return Err(ParquetError::General(
            "input is not a supported Sage spectral_library Parquet file".into(),
        ));
    }

    let mut entries = Vec::<DdaLibraryEntry>::new();
    let mut entry_indices = HashMap::<String, usize>::new();
    for (row_index, row) in reader.get_row_iter(None)?.enumerate() {
        let row_number = row_index + 1;
        let row = row?;
        let id = text(
            required(&row, "library_entry_id", row_number)?,
            "library_entry_id",
            row_number,
        )?;
        let proforma = text(
            required(&row, "proforma", row_number)?,
            "proforma",
            row_number,
        )?;
        let stripped_peptide = text(
            required(&row, "stripped_peptide", row_number)?,
            "stripped_peptide",
            row_number,
        )?;
        let proteins = text(
            required(&row, "proteins", row_number)?,
            "proteins",
            row_number,
        )?;
        let label_channel = (schema_version >= 2)
            .then(|| {
                text(
                    required(&row, "label_channel", row_number)?,
                    "label_channel",
                    row_number,
                )
            })
            .transpose()?
            .filter(|channel| !channel.is_empty());
        let label_group = (schema_version >= 2)
            .then(|| {
                text(
                    required(&row, "label_group", row_number)?,
                    "label_group",
                    row_number,
                )
            })
            .transpose()?
            .filter(|group| !group.is_empty());
        let label_reference = (schema_version >= 2)
            .then(|| {
                text(
                    required(&row, "label_reference", row_number)?,
                    "label_reference",
                    row_number,
                )
            })
            .transpose()?
            .filter(|reference| !reference.is_empty());
        let source_file = text(
            required(&row, "source_file", row_number)?,
            "source_file",
            row_number,
        )?;
        let source_spectrum = text(
            required(&row, "source_spectrum", row_number)?,
            "source_spectrum",
            row_number,
        )?;
        let precursor_charge = charge(
            required(&row, "precursor_charge", row_number)?,
            "precursor_charge",
            row_number,
        )?;
        let precursor_neutral_mass = float(
            required(&row, "precursor_neutral_mass", row_number)?,
            "precursor_neutral_mass",
            row_number,
        )?;
        let precursor_mz = float(
            required(&row, "precursor_mz", row_number)?,
            "precursor_mz",
            row_number,
        )?;
        let retention_time_minutes = float(
            required(&row, "aligned_retention_time_minutes", row_number)?,
            "aligned_retention_time_minutes",
            row_number,
        )?;
        let ion_mobility = float(
            required(&row, "ion_mobility", row_number)?,
            "ion_mobility",
            row_number,
        )?;
        let source_spectrum_q = float(
            required(&row, "spectrum_q", row_number)?,
            "spectrum_q",
            row_number,
        )?;

        let entry_index = match entry_indices.get(&id).copied() {
            Some(entry_index) => {
                let existing = &entries[entry_index];
                if existing.proforma != proforma
                    || existing.precursor_charge != precursor_charge
                    || existing.precursor_neutral_mass != precursor_neutral_mass
                    || existing.label_channel != label_channel
                    || existing.label_group != label_group
                    || existing.label_reference != label_reference
                {
                    return Err(ParquetError::General(format!(
                        "spectral-library entry `{id}` has inconsistent precursor metadata"
                    )));
                }
                entry_index
            }
            None => {
                let entry_index = entries.len();
                entry_indices.insert(id.clone(), entry_index);
                entries.push(DdaLibraryEntry {
                    library_entry_id: id,
                    source_file,
                    source_spectrum,
                    proforma,
                    stripped_peptide,
                    proteins,
                    label_channel,
                    label_group,
                    label_reference,
                    precursor_charge,
                    precursor_neutral_mass,
                    precursor_mz,
                    retention_time_minutes,
                    ion_mobility,
                    source_spectrum_q,
                    is_decoy: false,
                    fragments: Vec::new(),
                });
                entry_index
            }
        };

        let fragment_charge = charge(
            required(&row, "fragment_charge", row_number)?,
            "fragment_charge",
            row_number,
        )?;
        entries[entry_index].fragments.push(LibraryFragment {
            kind: fragment_kind(required(&row, "fragment_type", row_number)?, row_number)?,
            ordinal: i32::try_from(integer(
                required(&row, "fragment_ordinal", row_number)?,
                "fragment_ordinal",
                row_number,
            )?)
            .map_err(|_| {
                ParquetError::General(format!(
                    "spectral-library row {row_number} has out-of-range `fragment_ordinal`"
                ))
            })?,
            charge: i32::from(fragment_charge),
            neutral_loss: float(
                required(&row, "neutral_loss", row_number)?,
                "neutral_loss",
                row_number,
            )?,
            mz: float(
                required(&row, "fragment_mz", row_number)?,
                "fragment_mz",
                row_number,
            )?,
            relative_intensity: float(
                required(&row, "relative_intensity", row_number)?,
                "relative_intensity",
                row_number,
            )?,
        });
    }

    DdaLibraryIndex::new(entries)
        .map(|index| index.entries().to_vec())
        .map_err(ParquetError::General)
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

pub fn build_lfq_schema() -> parquet::errors::Result<Type> {
    parquet::schema::parser::parse_message_type(include_str!(
        "../../../schemas/lfq.v1.parquet.schema"
    ))
}

fn build_lfq_schema_version(has_labels: bool) -> parquet::errors::Result<Type> {
    parquet::schema::parser::parse_message_type(if has_labels {
        include_str!("../../../schemas/lfq.v2.parquet.schema")
    } else {
        include_str!("../../../schemas/lfq.v1.parquet.schema")
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
                std::iter::repeat(val).take(filenames.len())
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
                std::iter::repeat(val).take(filenames.len())
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
                    std::iter::repeat(value).take(filenames.len())
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
                    std::iter::repeat(value).take(filenames.len())
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
                std::iter::repeat(val).take(filenames.len())
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
            .flat_map(|((_, decoy), _)| std::iter::repeat(*decoy).take(filenames.len()))
            .collect::<Vec<_>>();

        col.typed::<BoolType>().write_batch(&values, None, None)?;
        col.close()?;
    }

    if let Some(mut col) = rg.next_column()? {
        let values = rows
            .iter()
            .copied()
            .flat_map(|(_, quantified)| {
                std::iter::repeat(quantified.peak.q_value).take(filenames.len())
            })
            .collect::<Vec<_>>();

        col.typed::<FloatType>().write_batch(&values, None, None)?;
        col.close()?;
    }

    if let Some(mut col) = rg.next_column()? {
        let values = rows
            .iter()
            .copied()
            .flat_map(|(_, quantified)| {
                std::iter::repeat(quantified.peak.score).take(filenames.len())
            })
            .collect::<Vec<_>>();

        col.typed::<DoubleType>().write_batch(&values, None, None)?;
        col.close()?;
    }

    if let Some(mut col) = rg.next_column()? {
        let values = rows
            .iter()
            .copied()
            .flat_map(|(_, quantified)| {
                std::iter::repeat(quantified.peak.spectral_angle).take(filenames.len())
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

#[cfg(test)]
#[path = "../tests/unit/parquet.rs"]
mod ptm_tests;
