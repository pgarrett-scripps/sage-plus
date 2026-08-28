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
    SpectralLibraryEntry, SpectralLibrarySettings, SpectralLibraryStrategy,
};
use sage_core::tmt::TmtQuant;

macro_rules! write_required_column {
    ($row_group:expr, $values:expr, $ty:ident) => {
        if let Some(mut column) = $row_group.next_column()? {
            column.typed::<$ty>().write_batch(&$values, None, None)?;
            column.close()?;
        }
    };
}

mod lfq;
mod query;
mod results;
mod sites;
mod spectral_library;

pub use lfq::{build_lfq_schema, serialize_lfq};
pub use query::scan_json_rows;
pub use results::{
    build_matched_fragment_schema, build_schema, serialize_features, serialize_matched_fragments,
};
pub use sites::{
    deserialize_custom_cleavage_sites, deserialize_ptm_library, serialize_protein_sites,
    serialize_ptm_library, serialize_ptm_sites, ProteinSiteRecord, PtmSiteRecord,
};
pub use spectral_library::{build_spectral_library_schema, serialize_spectral_library};

#[cfg(test)]
#[path = "../tests/unit/parquet.rs"]
mod ptm_tests;
