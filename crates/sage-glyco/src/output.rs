//! `glyco.sage.parquet`: one row per explained glyco candidate.

use parquet::basic::ZstdLevel;
use parquet::data_type::{BoolType, ByteArray, ByteArrayType, DoubleType, FloatType, Int32Type};
use parquet::file::metadata::KeyValue;
use parquet::file::properties::WriterProperties;
use parquet::file::writer::SerializedFileWriter;
use sage_core::database::IndexedDatabase;
use sage_core::peptide::Site;

use crate::composition::GlycanLibrary;
use crate::fdr::GlycoPsm;

/// Output file name.
pub const FILE_NAME: &str = "glyco.sage.parquet";

pub fn schema() -> parquet::errors::Result<parquet::schema::types::Type> {
    parquet::schema::parser::parse_message_type(include_str!(
        "../../../schemas/glyco.sage.v1.parquet.schema"
    ))
}

macro_rules! column {
    ($row_group:expr, $values:expr, $ty:ident) => {
        if let Some(mut column) = $row_group.next_column()? {
            column.typed::<$ty>().write_batch(&$values, None, None)?;
            column.close()?;
        }
    };
}

fn winner(p: &GlycoPsm) -> &crate::search::Explanation {
    &p.candidate.explanations[p.assignment.explanation]
}

/// Serialize `psms` (in any order) to Parquet bytes. `filenames` is indexed
/// by `Feature::file_id`.
pub fn serialize(
    psms: &[GlycoPsm],
    db: &IndexedDatabase,
    library: &GlycanLibrary,
    filenames: &[String],
) -> parquet::errors::Result<Vec<u8>> {
    let options = WriterProperties::builder()
        .set_compression(parquet::basic::Compression::ZSTD(ZstdLevel::try_new(3)?))
        .set_key_value_metadata(Some(vec![
            KeyValue::new("sage.schema.name".into(), Some("glyco".into())),
            KeyValue::new("sage.schema.version".into(), Some("1".into())),
            KeyValue::new("sage.schema.status".into(), Some("experimental".into())),
        ]))
        .build();
    let mut writer = SerializedFileWriter::new(Vec::new(), schema()?.into(), options.into())?;
    for rows in psms.chunks(65_536) {
        let mut rg = writer.next_row_group()?;
        let peptide = |p: &GlycoPsm| &db.peptides[p.candidate.feature.peptide_idx.0 as usize];
        let counts = |p: &GlycoPsm| {
            let e = winner(p);
            if p.assignment.decoy {
                e.decoy
            } else {
                e.target
            }
        };
        let text = |f: &dyn Fn(&GlycoPsm) -> String| {
            rows.iter()
                .map(|p| f(p).as_str().into())
                .collect::<Vec<ByteArray>>()
        };
        column!(
            rg,
            text(&|p| filenames
                .get(p.candidate.feature.file_id)
                .cloned()
                .unwrap_or_default()),
            ByteArrayType
        );
        column!(
            rg,
            text(&|p| p.candidate.feature.spec_id.clone()),
            ByteArrayType
        );
        column!(
            rg,
            rows.iter()
                .map(|p| p.candidate.feature.label)
                .collect::<Vec<_>>(),
            Int32Type
        );
        column!(rg, text(&|p| peptide(p).to_string()), ByteArrayType);
        column!(
            rg,
            text(&|p| peptide(p).sequence.as_str().to_string()),
            ByteArrayType
        );
        column!(
            rg,
            text(&|p| peptide(p).proteins(&db.decoy_tag, db.generate_decoys)),
            ByteArrayType
        );
        column!(
            rg,
            rows.iter()
                .map(|p| match p.candidate.feature.mass_offset.map(|m| m.site) {
                    Some(Site::Sequence(index)) => index as i32 + 1,
                    _ => 0,
                })
                .collect::<Vec<_>>(),
            Int32Type
        );
        column!(
            rg,
            rows.iter()
                .map(|p| p.candidate.feature.charge as i32)
                .collect::<Vec<_>>(),
            Int32Type
        );
        column!(
            rg,
            rows.iter()
                .map(|p| p.candidate.feature.expmass)
                .collect::<Vec<_>>(),
            FloatType
        );
        column!(
            rg,
            rows.iter()
                .map(|p| p.candidate.peptide_mass)
                .collect::<Vec<_>>(),
            FloatType
        );
        let composition = |p: &GlycoPsm| {
            library
                .get(winner(p).composition as usize)
                .expect("library index")
        };
        column!(rg, text(&|p| composition(p).to_string()), ByteArrayType);
        column!(
            rg,
            rows.iter()
                .map(|p| composition(p).mass())
                .collect::<Vec<_>>(),
            DoubleType
        );
        column!(
            rg,
            rows.iter()
                .map(|p| winner(p).isotope as i32)
                .collect::<Vec<_>>(),
            Int32Type
        );
        column!(
            rg,
            rows.iter()
                .map(|p| winner(p).adducts as i32)
                .collect::<Vec<_>>(),
            Int32Type
        );
        column!(
            rg,
            rows.iter().map(|p| winner(p).error_ppm).collect::<Vec<_>>(),
            FloatType
        );
        column!(
            rg,
            rows.iter()
                .map(|p| p.candidate.explanations.len() as i32)
                .collect::<Vec<_>>(),
            Int32Type
        );
        column!(
            rg,
            rows.iter()
                .map(|p| p.candidate.feature.hyperscore)
                .collect::<Vec<_>>(),
            DoubleType
        );
        column!(
            rg,
            rows.iter().map(|p| p.discriminant).collect::<Vec<_>>(),
            DoubleType
        );
        column!(
            rg,
            rows.iter().map(|p| p.peptide_q).collect::<Vec<_>>(),
            FloatType
        );
        column!(
            rg,
            rows.iter().map(|p| p.assignment.score).collect::<Vec<_>>(),
            DoubleType
        );
        column!(
            rg,
            rows.iter()
                .map(|p| {
                    let margin = p.assignment.score - p.assignment.runner_up;
                    if margin.is_finite() {
                        margin
                    } else {
                        f64::MAX
                    }
                })
                .collect::<Vec<_>>(),
            DoubleType
        );
        column!(
            rg,
            rows.iter().map(|p| p.glycan_q).collect::<Vec<_>>(),
            FloatType
        );
        column!(
            rg,
            rows.iter().map(|p| p.assignment.decoy).collect::<Vec<_>>(),
            BoolType
        );
        column!(
            rg,
            rows.iter()
                .map(|p| counts(p).y_matched() as i32)
                .collect::<Vec<_>>(),
            Int32Type
        );
        column!(
            rg,
            rows.iter()
                .map(|p| counts(p).y_generated() as i32)
                .collect::<Vec<_>>(),
            Int32Type
        );
        column!(
            rg,
            rows.iter()
                .map(|p| p.candidate.oxonium_ions as i32)
                .collect::<Vec<_>>(),
            Int32Type
        );
        column!(
            rg,
            rows.iter()
                .map(|p| p.candidate.feature.rt)
                .collect::<Vec<_>>(),
            FloatType
        );
        column!(
            rg,
            rows.iter().map(|p| p.passes_fdr).collect::<Vec<_>>(),
            BoolType
        );
        rg.close()?;
    }
    writer.into_inner()
}
