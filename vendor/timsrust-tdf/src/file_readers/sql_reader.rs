pub(crate) mod frame_groups;
pub(crate) mod frames;
pub(crate) mod metadata;
pub(crate) mod pasef_frame_msms;
pub(crate) mod precursors;
pub(crate) mod quad_settings;

use std::collections::HashMap;

use serde::Deserialize;

use crate::{TDFPathError, TDFPathLike};

pub(crate) struct SqlReader {
    inner: timsrust_reader::MetadataDatabase,
}

impl SqlReader {
    pub(crate) fn open(path: impl TDFPathLike) -> Result<Self, SqlReaderError> {
        let path = path.to_timstof_path()?;
        let inner = timsrust_reader::MetadataDatabase::open(path.tdf().as_ref())?;
        Ok(Self { inner })
    }

    pub(crate) fn read_table<T: for<'de> Deserialize<'de>>(
        &self,
        table: &str,
    ) -> Result<Vec<T>, SqlReaderError> {
        self.inner
            .table(table)?
            .into_iter()
            .map(|row| {
                serde_json::from_value(serde_json::Value::Object(row)).map_err(SqlReaderError::from)
            })
            .collect()
    }
}

pub(crate) trait ReadableSqlTable: for<'de> Deserialize<'de> {
    fn table_name() -> &'static str;

    fn from_sql_reader(reader: &SqlReader) -> Result<Vec<Self>, SqlReaderError>
    where
        Self: Sized,
    {
        Ok(reader.read_table::<Self>(Self::table_name())?)
    }
}

pub(crate) trait ReadableSqlHashMap {
    fn table_name() -> &'static str;

    fn from_sql_reader(reader: &SqlReader) -> Result<HashMap<String, String>, SqlReaderError>
    where
        Self: Sized,
    {
        #[derive(Deserialize)]
        struct KvRow {
            #[serde(rename = "Key")]
            key: String,
            #[serde(rename = "Value")]
            value: String,
        }
        let rows = reader.read_table::<KvRow>(Self::table_name())?;
        Ok(rows.into_iter().map(|r| (r.key, r.value)).collect())
    }
}

#[derive(Debug, thiserror::Error)]
pub(crate) enum SqlReaderError {
    #[error("{0}")]
    TDFPathError(#[from] TDFPathError),
    #[error("{0}")]
    Native(#[from] timsrust_reader::Error),
    #[error("{0}")]
    Json(#[from] serde_json::Error),
}

impl std::fmt::Debug for SqlReader {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("SqlReader").finish_non_exhaustive()
    }
}
