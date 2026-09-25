//! `crosslinks.sage.parquet`: one row per spectrum's best crosslink match.

use crate::search::{peptide_position, protein_position, Csm};
use parquet::basic::{Compression, ZstdLevel};
use parquet::data_type::{BoolType, ByteArray, ByteArrayType, DoubleType, FloatType, Int32Type};
use parquet::file::metadata::KeyValue;
use parquet::file::properties::WriterProperties;
use parquet::file::writer::SerializedFileWriter;
use sage_core::database::IndexedDatabase;
use std::sync::Arc;

pub const SCHEMA: &str = "
message crosslinks {
  required binary filename (STRING);
  required binary spectrum_id (STRING);
  required int32 scannr;
  required float rt;
  required int32 charge;
  required float expmass;
  required float calcmass;
  required int32 isotope_error;
  required float precursor_ppm;
  required binary alpha_peptide (STRING);
  required binary alpha_proteins (STRING);
  required int32 alpha_link_site;
  required int32 alpha_protein_position;
  required boolean alpha_doublet;
  required double alpha_hyperscore;
  required int32 alpha_matched_peaks;
  required binary beta_peptide (STRING);
  required binary beta_proteins (STRING);
  required int32 beta_link_site;
  required int32 beta_protein_position;
  required boolean beta_doublet;
  required double beta_hyperscore;
  required int32 beta_matched_peaks;
  required binary crosslink_class (STRING);
  required boolean is_decoy;
  required binary link_type (STRING);
  required double hyperscore;
  required double delta_next;
  required int32 matched_peaks;
  required float matched_intensity_pct;
  required int32 doublets;
  required double discriminant_score;
  required float csm_q;
  required float residue_pair_q;
}
";

macro_rules! column {
    ($row_group:expr, $ty:ident, $values:expr) => {
        if let Some(mut column) = $row_group.next_column()? {
            column.typed::<$ty>().write_batch(&$values, None, None)?;
            column.close()?;
        }
    };
}

/// Scan number from a Thermo/mzML native ID, or the ID itself when numeric.
fn scan_number(id: &str) -> i32 {
    id.rsplit_once("scan=")
        .map_or(id, |(_, scan)| scan)
        .split(|c: char| !c.is_ascii_digit())
        .next()
        .and_then(|scan| scan.parse().ok())
        .unwrap_or(0)
}

/// Serialize `csms` (already filtered) to parquet bytes.
pub fn serialize(
    csms: &[&Csm],
    db: &IndexedDatabase,
    filenames: &[String],
    linker: &str,
) -> parquet::errors::Result<Vec<u8>> {
    let schema = Arc::new(parquet::schema::parser::parse_message_type(SCHEMA)?);
    let options = WriterProperties::builder()
        .set_compression(Compression::ZSTD(ZstdLevel::try_new(3)?))
        .set_key_value_metadata(Some(vec![
            KeyValue::new("sage.schema.name".into(), Some("crosslinks".into())),
            KeyValue::new("sage.schema.version".into(), Some("1".into())),
            KeyValue::new("sage.crosslink.linker".into(), Some(linker.into())),
            KeyValue::new(
                "sage.crosslink.fdr".into(),
                Some("(TD-DD)/TT, separately for intra and inter links".into()),
            ),
        ]))
        .build();
    let mut writer = SerializedFileWriter::new(Vec::new(), schema, Arc::new(options))?;
    let mut rg = writer.next_row_group()?;

    let text = |f: &dyn Fn(&Csm) -> String| -> Vec<ByteArray> {
        csms.iter()
            .map(|c| ByteArray::from(f(c).as_str()))
            .collect()
    };
    let int = |f: &dyn Fn(&Csm) -> i32| -> Vec<i32> { csms.iter().map(|c| f(c)).collect() };
    let float = |f: &dyn Fn(&Csm) -> f32| -> Vec<f32> { csms.iter().map(|c| f(c)).collect() };
    let double = |f: &dyn Fn(&Csm) -> f64| -> Vec<f64> { csms.iter().map(|c| f(c)).collect() };
    let boolean = |f: &dyn Fn(&Csm) -> bool| -> Vec<bool> { csms.iter().map(|c| f(c)).collect() };

    column!(
        rg,
        ByteArrayType,
        text(&|c| filenames.get(c.file_id).cloned().unwrap_or_default())
    );
    column!(rg, ByteArrayType, text(&|c| c.spectrum_id.clone()));
    column!(rg, Int32Type, int(&|c| scan_number(&c.spectrum_id)));
    column!(rg, FloatType, float(&|c| c.rt));
    column!(rg, Int32Type, int(&|c| c.charge as i32));
    column!(rg, FloatType, float(&|c| c.expmass));
    column!(rg, FloatType, float(&|c| c.calcmass));
    column!(rg, Int32Type, int(&|c| c.isotope_error as i32));
    column!(rg, FloatType, float(&|c| c.precursor_ppm));
    for beta in [false, true] {
        let chain = move |c: &Csm| {
            if beta {
                c.beta.clone()
            } else {
                c.alpha.clone()
            }
        };
        column!(
            rg,
            ByteArrayType,
            text(&|c| db[chain(c).peptide].to_string())
        );
        column!(
            rg,
            ByteArrayType,
            text(&|c| db[chain(c).peptide].proteins.join(";"))
        );
        column!(
            rg,
            Int32Type,
            int(&|c| {
                let m = chain(c);
                peptide_position(&db[m.peptide], m.site) as i32
            })
        );
        column!(
            rg,
            Int32Type,
            int(&|c| {
                let m = chain(c);
                protein_position(&db[m.peptide], m.site).map_or(0, |(_, p)| p as i32)
            })
        );
        column!(rg, BoolType, boolean(&|c| chain(c).doublet));
        column!(rg, DoubleType, double(&|c| chain(c).hyperscore));
        column!(rg, Int32Type, int(&|c| chain(c).matched_peaks as i32));
    }
    column!(rg, ByteArrayType, text(&|c| c.class.as_str().into()));
    column!(
        rg,
        BoolType,
        boolean(&|c| c.class != crate::search::Class::TT)
    );
    column!(
        rg,
        ByteArrayType,
        text(&|c| if c.intra { "intra" } else { "inter" }.into())
    );
    column!(rg, DoubleType, double(&|c| c.hyperscore));
    column!(rg, DoubleType, double(&|c| c.delta_next));
    column!(rg, Int32Type, int(&|c| c.matched_peaks as i32));
    column!(rg, FloatType, float(&|c| c.matched_intensity_pct));
    column!(rg, Int32Type, int(&|c| c.doublets as i32));
    column!(rg, DoubleType, double(&|c| c.discriminant_score));
    column!(rg, FloatType, float(&|c| c.csm_q));
    column!(rg, FloatType, float(&|c| c.residue_pair_q));

    rg.close()?;
    writer.into_inner()
}

#[cfg(test)]
mod tests {
    use super::scan_number;

    #[test]
    fn scan_numbers_from_native_ids() {
        assert_eq!(
            scan_number("controllerType=0 controllerNumber=1 scan=35020"),
            35020
        );
        assert_eq!(scan_number("1234"), 1234);
        assert_eq!(scan_number("index=7"), 0);
    }
}
