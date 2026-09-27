//! Run-level quality-control summaries written next to the search results.

use super::*;
use sage_core::digestion::DigestionSummary;
use sage_core::enzyme::EnzymeParameters;

/// Spectrum- and peptide-level q-value a PSM must pass to count towards the
/// digestion summary, matching the run summary's 1% FDR counts.
const DIGESTION_Q_VALUE: f32 = 0.01;

#[derive(Debug, Clone, Default, serde::Serialize, serde::Deserialize)]
pub struct QcRunStats {
    #[serde(default)]
    pub digestion: DigestionRunStats,
}

/// Digestion summary of rank-1 PSMs at 1% spectrum and peptide q-value,
/// counted over distinct peptide sequences with decoys subtracted.
#[derive(Debug, Clone, Default, serde::Serialize, serde::Deserialize)]
pub struct DigestionRunStats {
    pub q_value: f32,
    pub files: Vec<DigestionFileStats>,
    pub total: DigestionSummary,
}

#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct DigestionFileStats {
    pub file: String,
    #[serde(flatten)]
    pub summary: DigestionSummary,
}

impl Runner {
    /// Summarize missed cleavages and non-enzymatic termini per file and for
    /// the whole run.
    pub(super) fn digestion_stats(
        &self,
        features: &[Feature],
        filenames: &[String],
    ) -> DigestionRunStats {
        let enzyme = EnzymeParameters::from(self.database_parameters.enzyme.clone()).enzyme;
        let passing = features
            .iter()
            .filter(|feature| {
                feature.rank == 1
                    && feature.spectrum_q <= DIGESTION_Q_VALUE
                    && feature.peptide_q <= DIGESTION_Q_VALUE
            })
            .collect::<Vec<_>>();
        let files = filenames
            .iter()
            .enumerate()
            .map(|(file_id, file)| DigestionFileStats {
                file: file.clone(),
                summary: sage_core::digestion::summarize(
                    enzyme.as_ref(),
                    passing
                        .iter()
                        .filter(|feature| feature.file_id == file_id)
                        .map(|feature| &self.database[feature.peptide_idx]),
                ),
            })
            .collect();
        let total = sage_core::digestion::summarize(
            enzyme.as_ref(),
            passing
                .iter()
                .map(|feature| &self.database[feature.peptide_idx]),
        );
        info!(
            "digestion: {} peptides (decoy-corrected), missed cleavage {:.1}%, semi N-terminal {:.1}%, semi C-terminal {:.1}%, non-enzymatic {:.1}%",
            total.peptides,
            total.missed_cleavage_pct,
            total.semi_n_pct,
            total.semi_c_pct,
            total.non_enzymatic_pct
        );
        DigestionRunStats {
            q_value: DIGESTION_Q_VALUE,
            files,
            total,
        }
    }

    pub(super) fn write_digestion(&self, stats: &DigestionRunStats) -> anyhow::Result<Url> {
        let path = self.make_path("digestion.tsv");
        sage_cloudpath::write_bytes_sync(&path, serialize_digestion(stats)?)?;
        Ok(path)
    }
}

pub(super) fn serialize_digestion(stats: &DigestionRunStats) -> anyhow::Result<Vec<u8>> {
    let mut writer = csv::WriterBuilder::new()
        .delimiter(b'\t')
        .from_writer(Vec::new());
    writer.write_record([
        "file",
        "target_peptides",
        "decoy_peptides",
        "peptides",
        "missed_cleavages_0",
        "missed_cleavages_1",
        "missed_cleavages_2_plus",
        "semi_n",
        "semi_c",
        "non_enzymatic",
        "missed_cleavage_pct",
        "semi_n_pct",
        "semi_c_pct",
        "non_enzymatic_pct",
    ])?;
    let rows = stats
        .files
        .iter()
        .map(|file| (file.file.as_str(), &file.summary))
        .chain(std::iter::once(("total", &stats.total)));
    for (file, summary) in rows {
        writer.write_record([
            file.to_string(),
            summary.target_peptides.to_string(),
            summary.decoy_peptides.to_string(),
            summary.peptides.to_string(),
            summary.missed_cleavages_0.to_string(),
            summary.missed_cleavages_1.to_string(),
            summary.missed_cleavages_2_plus.to_string(),
            summary.semi_n.to_string(),
            summary.semi_c.to_string(),
            summary.non_enzymatic.to_string(),
            format!("{:.3}", summary.missed_cleavage_pct),
            format!("{:.3}", summary.semi_n_pct),
            format!("{:.3}", summary.semi_c_pct),
            format!("{:.3}", summary.non_enzymatic_pct),
        ])?;
    }
    writer.flush()?;
    Ok(writer.into_inner()?)
}
