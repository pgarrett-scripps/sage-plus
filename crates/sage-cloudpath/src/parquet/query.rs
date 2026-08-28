use super::*;

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
