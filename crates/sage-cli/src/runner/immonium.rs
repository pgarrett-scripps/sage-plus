//! Immonium-ion outputs: `immonium.tsv` and the extra PIN columns.

use super::*;
use sage_core::immonium::ImmoniumEvidence;

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

impl Runner {
    /// Write `immonium.tsv`: one row per reported PSM. `None` when off.
    pub(super) fn write_immonium(
        &self,
        features: &[&Feature],
        filenames: &[String],
    ) -> anyhow::Result<Option<Url>> {
        let Some(settings) = self.parameters.immonium.as_ref() else {
            return Ok(None);
        };
        let mut writer = csv::WriterBuilder::new()
            .delimiter(b'\t')
            .from_writer(Vec::new());
        writer.write_record([
            "psm_id",
            "file",
            "scannr",
            "peptide",
            "is_decoy",
            "rank",
            "spectrum_q",
            "explained",
            "missing",
            "unexplained",
            "residue_ions",
            "modified_explained",
            "modified_unexplained",
            "modified_ions",
        ])?;
        for feature in features {
            let evidence: ImmoniumEvidence = feature.immonium.unwrap_or_default();
            writer.write_record([
                feature.psm_id.to_string().as_str(),
                filenames[feature.file_id].as_str(),
                feature.spec_id.as_str(),
                self.database[feature.peptide_idx].to_string().as_str(),
                if feature.label == -1 { "true" } else { "false" },
                feature.rank.to_string().as_str(),
                feature.spectrum_q.to_string().as_str(),
                evidence.explained.to_string().as_str(),
                evidence.missing.to_string().as_str(),
                evidence.unexplained.to_string().as_str(),
                evidence.residue_names().as_str(),
                evidence.modified_explained.to_string().as_str(),
                evidence.modified_unexplained.to_string().as_str(),
                settings.modified_names(evidence.modified_observed).as_str(),
            ])?;
        }
        writer.flush()?;
        let path = self.make_path("immonium.tsv");
        sage_cloudpath::write_bytes_sync(&path, writer.into_inner()?)?;
        Ok(Some(path))
    }
}
