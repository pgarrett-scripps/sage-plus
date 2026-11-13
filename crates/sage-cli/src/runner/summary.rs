//! The `run-summary.json` types: `RunSummary` and the per-stage stats it holds.
//!
//! Moved out of `runner.rs` unchanged; every `runner::*` path is unchanged.

use super::*;

#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct RunSummary {
    #[serde(default = "run_summary_schema_version")]
    pub schema_version: u32,
    pub runtime_secs: u64,
    pub files: usize,
    pub peptides_in_database: usize,
    pub fragments_in_database: usize,
    pub psms_at_one_percent_fdr: usize,
    pub peptides_at_one_percent_fdr: usize,
    pub proteins_at_one_percent_fdr: usize,
    pub protein_groups_at_one_percent_fdr: usize,
    #[serde(default)]
    pub ptm_localization: PtmLocalizationRunStats,
    #[serde(default)]
    pub models: ModelRunStats,
    #[serde(default)]
    pub quantification: QuantificationRunStats,
    #[serde(default)]
    pub execution: ExecutionRunStats,
    #[serde(default)]
    pub inputs: InputRunStats,
    #[serde(default)]
    pub modifications: ModificationRunStats,
    #[serde(default)]
    pub spectral_library: SpectralLibraryRunStats,
    /// Quality-control summaries (digestion, and the optional MS1 and MS2
    /// checks when enabled).
    #[serde(default)]
    pub qc: qc::QcRunStats,
    pub output_paths: Vec<String>,
    #[serde(default)]
    pub warnings: Vec<crate::events::RunWarning>,
    #[serde(default)]
    pub provenance: RunProvenance,
    /// Precursor and fragment tolerances suggested by the confident PSMs.
    #[serde(default)]
    pub recommended_tolerances: ToleranceRecommendation,
    /// In-source fragmentation and neutral-loss annotation; absent when the
    /// annotation is off, so a search without it writes the summary it wrote
    /// before.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub in_source: Option<InSourceRunStats>,
    /// Peak resident set size of the Sage process, in bytes, at the end of
    /// the run; `None` where the platform does not report it.
    #[serde(default)]
    pub peak_rss_bytes: Option<u64>,
}

/// In-source fragmentation and in-source neutral loss annotation of one run.
/// `flagged_psms` is how many PSMs were annotated; `modes` is the built-in
/// check summed over the files, one record per way co-elution was decided,
/// each with its own unit, tolerance and placements; `files` carries each
/// file's own record, including its difference histogram per mode. There is
/// no record pooling the modes: a difference in peak widths is not a
/// difference in minutes.
#[derive(Debug, Clone, Default, PartialEq, serde::Serialize, serde::Deserialize)]
pub struct InSourceRunStats {
    /// Configured co-elution tolerance, when one was set; otherwise each
    /// file's tolerance is in `files`.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub rt_tolerance: Option<f32>,
    pub flagged_psms: usize,
    /// One record per way co-elution was decided: `apex` where quantification
    /// gave both the product and the parent a peak in the file, `rt` where
    /// PSM retention times were compared.
    #[serde(default)]
    pub modes: Vec<sage_core::in_source::InSourceModeStats>,
    pub files: Vec<sage_core::in_source::InSourceFileStats>,
}

/// Search tolerances suggested by the signed mass errors of rank-1 target
/// PSMs at 1% spectrum q-value, pooled over all files. Errors are the raw,
/// uncorrected ones, so the suggestion applies to a search without
/// `mass_recalibration`.
#[derive(Debug, Clone, Default, PartialEq, serde::Serialize, serde::Deserialize)]
pub struct ToleranceRecommendation {
    /// Confident PSMs the estimates were computed from.
    pub psms: usize,
    /// Precursor estimate, from `precursor_ppm`.
    pub precursor: Option<ToleranceEstimate>,
    /// Fragment estimate, from each PSM's signed fragment error and its
    /// ions' spread around it.
    pub fragment: Option<ToleranceEstimate>,
    /// Why no recommendation was made, such as too few PSMs.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub skipped: Option<String>,
}

impl ToleranceRecommendation {
    /// Estimate tolerances from rank-1 target PSMs at 1% spectrum q-value.
    /// Spectrum q-values must already be assigned.
    pub fn from_features(features: &[Feature]) -> Self {
        let confident = features
            .iter()
            .filter(|feature| feature.rank == 1 && feature.label == 1 && feature.spectrum_q <= 0.01)
            .collect::<Vec<_>>();
        let psms = confident.len();
        if psms < MIN_TOLERANCE_PSMS {
            return Self {
                psms,
                skipped: Some(format!(
                    "{psms} confident PSMs, fewer than the {MIN_TOLERANCE_PSMS} needed"
                )),
                ..Self::default()
            };
        }
        let precursor = confident
            .iter()
            .map(|feature| feature.delta_mass)
            .collect::<Vec<_>>();
        let fragment = confident
            .iter()
            .map(|feature| feature.signed_fragment_ppm)
            .collect::<Vec<_>>();
        let fragment_sd = confident
            .iter()
            .map(|feature| feature.fragment_ppm_sd)
            .collect::<Vec<_>>();
        Self {
            psms,
            precursor: ToleranceEstimate::from_errors(&precursor),
            fragment: ToleranceEstimate::from_fragment_errors(&fragment, &fragment_sd),
            skipped: None,
        }
    }

    /// One log line, for example
    /// `recommended tolerances: precursor ±10 ppm, fragment ±20 ppm (from 5321 PSMs)`.
    pub fn log_line(&self) -> String {
        if let Some(reason) = &self.skipped {
            return format!("recommended tolerances: skipped ({reason})");
        }
        let describe = |estimate: &Option<ToleranceEstimate>| match estimate {
            Some(ToleranceEstimate {
                recommended_ppm: Some(ppm),
                ..
            }) => format!("±{ppm} ppm"),
            Some(estimate) => format!("wider than ±100 ppm (needs ±{:.1})", estimate.required_ppm),
            None => "not estimated".into(),
        };
        format!(
            "recommended tolerances: precursor {}, fragment {} (from {} PSMs)",
            describe(&self.precursor),
            describe(&self.fragment),
            self.psms
        )
    }
}

#[derive(Debug, Clone, Default, serde::Serialize, serde::Deserialize)]
pub struct RunProvenance {
    pub software_version: String,
    pub mzmlb_enabled: bool,
    pub input_identity_mode: String,
    pub inputs: Vec<InputIdentity>,
    /// The provenance record written to every Parquet footer, keyed as there
    /// without the `sage.` prefix. Absent in summaries from before Beta 13.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub metadata: Option<crate::provenance::Provenance>,
}

#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct InputIdentity {
    pub path: String,
    pub size_bytes: Option<u64>,
    pub modified_unix_secs: Option<u64>,
}

impl InputIdentity {
    pub(super) fn from_url(url: &Url) -> Self {
        let metadata = url
            .to_file_path()
            .ok()
            .and_then(|path| std::fs::metadata(path).ok());
        Self {
            path: url.to_string(),
            size_bytes: metadata
                .as_ref()
                .filter(|metadata| metadata.is_file())
                .map(std::fs::Metadata::len),
            modified_unix_secs: metadata
                .and_then(|metadata| metadata.modified().ok())
                .and_then(|time| time.duration_since(std::time::UNIX_EPOCH).ok())
                .map(|time| time.as_secs()),
        }
    }
}

#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct MassAlignmentFileStats {
    pub file_id: usize,
    pub calibration_psms: usize,
    pub precursor: Option<sage_core::mass_calibration::CalibrationModel>,
    pub fragment: Option<sage_core::mass_calibration::CalibrationModel>,
    pub precursor_skip_reason: Option<String>,
    pub fragment_skip_reason: Option<String>,
}

/// Search-time mass recalibration of one file: the discovery pass and the
/// precursor and fragment models it selected.
#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct MassRecalibrationFileStats {
    pub file_id: usize,
    /// MS2 spectra searched in the discovery pass.
    pub discovery_spectra: usize,
    /// Rank-1 target PSMs at 1% Poisson spectrum q-value in the discovery
    /// pass, with q-values computed within each acquisition group.
    pub discovery_psms: usize,
    pub discovery_ms: u64,
    pub precursor: sage_core::mass_recalibration::ModelSelection,
    /// Fragment models, one per acquisition group (analyzer and activation).
    pub fragment: Vec<sage_core::mass_recalibration::GroupModelSelection>,
    #[serde(skip)]
    pub correction: sage_core::mass_recalibration::FileMassCorrection,
}

#[derive(Debug, Clone, Default, serde::Serialize, serde::Deserialize)]
pub struct MassRecalibrationRunStats {
    /// Configured mode: `static`, `linear`, or `auto`.
    pub mode: String,
    pub files: Vec<MassRecalibrationFileStats>,
}

const fn run_summary_schema_version() -> u32 {
    1
}

#[derive(Debug, Clone, Default, serde::Serialize, serde::Deserialize)]
pub struct PtmLocalizationRunStats {
    pub enabled: bool,
    pub localized_psms: usize,
    pub psm_q_value: f32,
    pub localization_q_value: f32,
    /// Separating-ion margin a localization needed to compete.
    #[serde(default)]
    pub min_separating_margin: u32,
    /// Target protein sites in the site-level FDR competition.
    #[serde(default)]
    pub target_protein_sites: usize,
    /// Decoy protein sites in the site-level FDR competition.
    #[serde(default)]
    pub decoy_protein_sites: usize,
    /// Target protein sites at site-level q-value <= 0.01.
    #[serde(default)]
    pub protein_sites_at_one_percent_fdr: usize,
}

#[derive(Debug, Clone, Copy, Default)]
pub(super) struct SiteFdrRunStats {
    pub(super) target_protein_sites: usize,
    pub(super) decoy_protein_sites: usize,
    pub(super) protein_sites_at_one_percent_fdr: usize,
}

impl SiteFdrRunStats {
    pub(super) fn from_sites(sites: &[artifacts::ProteinSite]) -> Self {
        let mut stats = Self::default();
        for site in sites {
            if site.decoy {
                stats.decoy_protein_sites += 1;
            } else {
                stats.target_protein_sites += 1;
                stats.protein_sites_at_one_percent_fdr += usize::from(site.q_value <= 0.01);
            }
        }
        stats
    }
}

#[derive(Debug, Clone, Default, serde::Serialize, serde::Deserialize)]
pub struct ModelRunStats {
    pub mass_alignment_applied: bool,
    #[serde(default)]
    pub mass_alignment_files: Vec<MassAlignmentFileStats>,
    /// Search-time mass recalibration; absent when `mass_recalibration` is off.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub mass_recalibration: Option<MassRecalibrationRunStats>,
    pub retention_time_prediction_enabled: bool,
    pub retention_time_model_fitted: bool,
    pub retention_time_features: String,
    pub retention_time_alignment: Option<String>,
    pub ion_mobility_observed: bool,
    pub ion_mobility_model_enabled: bool,
    pub ion_mobility_model_fitted: bool,
    pub ion_mobility_features: String,
    /// Scan-to-1/K0 scale of timsTOF inputs (`calibrated` or `linear`).
    /// Absent when no Bruker input was searched.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub ion_mobility_scale: Option<String>,
}

#[derive(Debug, Clone, Default, serde::Serialize, serde::Deserialize)]
pub struct QuantificationRunStats {
    pub lfq_enabled: bool,
    pub lfq_features: usize,
    pub tmt: Option<String>,
    pub tmt_features: usize,
    /// Per-channel reporter-ion coverage over all quantified spectra. Missing
    /// channels are counted, not averaged in as 0. Absent without TMT.
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub tmt_channels: Vec<sage_core::tmt::ReporterChannelSummary>,
    #[serde(default)]
    pub ms1_label_channels: usize,
    #[serde(default)]
    pub ms1_label_reference: Option<String>,
}

#[derive(Debug, Clone, Default, serde::Serialize, serde::Deserialize)]
pub struct ExecutionRunStats {
    pub batch_size: usize,
    pub parallelism: usize,
    #[serde(default)]
    pub rayon_threads: usize,
    pub max_memory_gb: Option<f64>,
}

#[derive(Debug, Clone, Default, serde::Serialize, serde::Deserialize)]
pub struct InputRunStats {
    pub mzml_files: usize,
    pub thermo_raw_files: usize,
    pub other_files: usize,
}

#[derive(Debug, Clone, Default, serde::Serialize, serde::Deserialize)]
pub struct ModificationRunStats {
    pub static_definitions: usize,
    pub variable_definitions: usize,
    pub max_variable_mods: usize,
    pub max_total_variable_mods: usize,
    pub max_combinations: Option<usize>,
    pub ptm_library_sites: usize,
    #[serde(default)]
    pub label_channels: usize,
    #[serde(default)]
    pub labeled_peptides: usize,
    /// Distinct search-time mass offset definitions.
    #[serde(default)]
    pub mass_offset_definitions: usize,
    /// Reported candidate PSMs whose peptidoform includes a mass offset.
    #[serde(default)]
    pub mass_offset_psms: usize,
    /// Distinct offset peptidoforms not already present in the index.
    #[serde(default)]
    pub mass_offset_peptidoforms: usize,
}

#[derive(Debug, Clone, Default, serde::Serialize, serde::Deserialize)]
pub struct SpectralLibraryRunStats {
    pub enabled: bool,
    pub entries: usize,
    pub transitions: usize,
    pub strategy: String,
    pub psm_q_value: f32,
    pub peptide_q_value: f32,
    pub formats: Vec<String>,
}

#[cfg(test)]
#[path = "../../tests/unit/runner/summary.rs"]
mod test;
