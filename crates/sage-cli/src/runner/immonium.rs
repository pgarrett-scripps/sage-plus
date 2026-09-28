//! Immonium-ion PIN columns. The Parquet columns are written with the rest of
//! `results.sage.parquet`.

use super::*;

/// PIN columns added before `Peptide` when `immonium` is enabled.
pub(super) const PIN_COLUMNS: [&str; 5] = [
    "immonium_explained",
    "immonium_missing",
    "immonium_unexplained",
    "immonium_modified_explained",
    "immonium_modified_unexplained",
];

pub(super) fn push_pin_fields(record: &mut csv::ByteRecord, feature: &Feature) {
    let evidence = feature.immonium.unwrap_or_default();
    for value in evidence.lda_row() {
        record.push_field(itoa::Buffer::new().format(value as u8).as_bytes());
    }
}
