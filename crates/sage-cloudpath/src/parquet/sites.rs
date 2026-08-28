use super::*;

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
