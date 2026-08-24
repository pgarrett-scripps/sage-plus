use crate::memory::MemoryLimits;
use anyhow::{ensure, Context};
use clap::ArgMatches;
use sage_cloudpath::tdf::BrukerProcessingConfig;
use sage_cloudpath::Url;
use sage_core::scoring::ScoreType;
use sage_core::{
    database::{Builder, Parameters},
    lfq::LfqSettings,
    mass::Tolerance,
    ml::mobility_model::IonMobilitySettings,
    ml::retention_alignment::AlignmentMethod,
    ml::retention_model::RetentionTimeSettings,
    spectral_library::SpectralLibrarySettings,
    spectral_library_search::LibrarySearchSettings,
    tmt::Isobaric,
};
use serde::{Deserialize, Serialize};

#[derive(Serialize, Deserialize, Clone, Copy, Debug)]
#[serde(default)]
pub struct PtmLocalizationSettings {
    /// Compute PTM site localization and write site-level reports.
    pub enabled: bool,
    /// Identification q-value cutoff for PSMs included in PTM site reports.
    #[serde(alias = "q_value")]
    pub psm_q_value: f32,
    /// Arrangement-level false-localization-rate q-value cutoff.
    pub localization_q_value: f32,
}

impl Default for PtmLocalizationSettings {
    fn default() -> Self {
        Self {
            enabled: false,
            psm_q_value: 0.01,
            localization_q_value: 0.01,
        }
    }
}

#[derive(Serialize, Clone)]
/// Actual search parameters - may include overrides or default values not set by user
pub struct Search {
    pub version: String,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub database: Option<Parameters>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub library_search: Option<LibrarySearchSettings>,
    pub quant: QuantSettings,
    pub precursor_tol: Tolerance,
    pub fragment_tol: Tolerance,
    pub precursor_charge: (u8, u8),
    pub override_precursor_charge: bool,
    pub isotope_errors: (i8, i8),
    pub deisotope: bool,
    pub chimera: bool,
    pub wide_window: bool,
    pub min_peaks: usize,
    pub max_peaks: usize,
    pub max_fragment_charge: Option<u8>,
    pub min_matched_peaks: u16,
    pub report_psms: usize,
    pub output_filter: OutputFilter,
    pub predict_rt: bool,
    pub retention_time_model: RetentionTimeSettings,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub retention_time_alignment: Option<AlignmentMethod>,
    pub ion_mobility_model: IonMobilitySettings,
    pub mzml_paths: Vec<Url>,
    pub output_paths: Vec<Url>,
    pub bruker_config: BrukerProcessingConfig,
    pub protein_grouping: bool,
    pub protein_grouping_peptide_fdr: f32,
    /// Maximum resident memory Sage may use, in GiB. `None` or zero disables the limit.
    pub max_memory_gb: Option<f64>,
    /// Minimum system memory Sage must leave available, in GiB. `None` or zero disables the limit.
    pub min_free_memory_gb: Option<f64>,
    /// Number of input files to load and search at once.
    pub batch_size: usize,

    pub ptm_localization: PtmLocalizationSettings,
    pub spectral_library: SpectralLibrarySettings,

    /// ppm threshold below which a precursor delta mass is treated as no shift
    /// for sequence-ambiguity annotation (`ambiguity_sequence` / `mass_shift`)
    pub mass_shift_ppm: f32,

    #[serde(skip_serializing)]
    pub output_directory: Url,

    #[serde(skip_serializing)]
    pub write_pin: bool,

    #[serde(skip_serializing)]
    pub write_report: bool,

    #[serde(skip_serializing)]
    pub annotate_matches: bool,

    pub score_type: ScoreType,
    /// Use the bitmap-based preliminary search (requires `deisotope: true`).
    pub use_bitmap: bool,
}

#[derive(Deserialize)]
/// Input search parameters deserialized from JSON file
pub struct Input {
    pub database: Option<Builder>,
    pub library_search: Option<LibrarySearchSettings>,
    pub precursor_tol: Tolerance,
    pub fragment_tol: Tolerance,
    pub report_psms: Option<usize>,
    pub output_filter: Option<OutputFilter>,
    pub chimera: Option<bool>,
    pub wide_window: Option<bool>,
    pub min_peaks: Option<usize>,
    pub max_peaks: Option<usize>,
    pub max_fragment_charge: Option<u8>,
    pub min_matched_peaks: Option<u16>,
    pub precursor_charge: Option<(u8, u8)>,
    pub override_precursor_charge: Option<bool>,
    pub isotope_errors: Option<(i8, i8)>,
    pub deisotope: Option<bool>,
    pub quant: Option<QuantOptions>,
    pub predict_rt: Option<bool>,
    pub retention_time_model: Option<RetentionTimeSettings>,
    pub retention_time_alignment: Option<AlignmentMethod>,
    pub ion_mobility_model: Option<IonMobilitySettings>,
    pub output_directory: Option<String>,
    pub mzml_paths: Option<Vec<String>>,
    pub bruker_config: Option<BrukerProcessingConfig>,
    pub protein_grouping: Option<bool>,
    pub protein_grouping_peptide_fdr: Option<f32>,
    pub max_memory_gb: Option<f64>,
    pub min_free_memory_gb: Option<f64>,
    pub batch_size: Option<usize>,

    pub ptm_localization: Option<PtmLocalizationSettings>,
    pub spectral_library: Option<SpectralLibrarySettings>,
    pub mass_shift_ppm: Option<f32>,

    pub annotate_matches: Option<bool>,
    pub write_pin: Option<bool>,
    pub write_report: Option<bool>,
    pub score_type: Option<ScoreType>,
    pub use_bitmap: Option<bool>,
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
#[serde(default)]
pub struct OutputFilter {
    /// Maximum spectrum-level q-value written to the PSM and matched-fragment tables.
    pub psm_q_value: f32,
}

impl Default for OutputFilter {
    fn default() -> Self {
        Self { psm_q_value: 0.1 }
    }
}

#[derive(Serialize, Deserialize, Debug)]
pub struct LfqOptions {
    pub peak_scoring: Option<sage_core::lfq::PeakScoringStrategy>,
    pub integration: Option<sage_core::lfq::IntegrationStrategy>,
    pub spectral_angle: Option<f64>,
    pub ppm_tolerance: Option<f32>,
    pub rt_pct_tolerance: Option<f32>,
    pub mobility_pct_tolerance: Option<f32>,
    pub combine_charge_states: Option<bool>,
    pub peptide_q_value: Option<f32>,
}

impl From<LfqOptions> for LfqSettings {
    fn from(value: LfqOptions) -> LfqSettings {
        let default = LfqSettings::default();
        let settings = LfqSettings {
            peak_scoring: value.peak_scoring.unwrap_or(default.peak_scoring),
            integration: value.integration.unwrap_or(default.integration),
            spectral_angle: value.spectral_angle.unwrap_or(default.spectral_angle).abs(),
            ppm_tolerance: value.ppm_tolerance.unwrap_or(default.ppm_tolerance).abs(),
            rt_pct_tolerance: value
                .rt_pct_tolerance
                .unwrap_or(default.rt_pct_tolerance)
                .abs(),
            peptide_q_value: value.peptide_q_value.unwrap_or(default.peptide_q_value),
            mobility_pct_tolerance: value
                .mobility_pct_tolerance
                .unwrap_or(default.mobility_pct_tolerance),
            combine_charge_states: value
                .combine_charge_states
                .unwrap_or(default.combine_charge_states),
        };
        if settings.ppm_tolerance > 20.0 {
            log::warn!("lfq_settings.ppm_tolerance is higher than expected");
        }
        if settings.rt_pct_tolerance > 2.0 {
            log::warn!("lfq_settings.rt_pct_tolerance is higher than expected");
        }
        if settings.rt_pct_tolerance < 0.05 {
            log::warn!("lfq_settings.rt_pct_tolerance is smaller than expected");
        }
        if settings.mobility_pct_tolerance > 4.0 {
            log::warn!("lfq_settings.mobility_pct_tolerance is higher than expected");
        }
        if settings.mobility_pct_tolerance < 0.05 {
            log::warn!("lfq_settings.mobility_pct_tolerance is smaller than expected");
        }
        if settings.spectral_angle < 0.50 {
            log::warn!("lfq_settings.spectral_angle is lower than expected");
        }
        if settings.peptide_q_value > 0.01 {
            log::info!("lfq_settings.peptide_q_value is higher than expected, expect increased runtime and memory usage");
        }
        if settings.peptide_q_value < 0.01 {
            log::warn!("lfq_settings.peptide_q_value is lower than expected, not all identified peptides will have MS1 intensities extracted");
        }

        settings
    }
}

#[cfg(test)]
mod tests {
    use super::LfqOptions;
    use sage_core::lfq::LfqSettings;

    #[test]
    fn parses_lfq_rt_percent_tolerance() {
        let options: LfqOptions = serde_json::from_str(r#"{"rt_pct_tolerance": 1.25}"#).unwrap();
        let settings: LfqSettings = options.into();

        assert_eq!(settings.rt_pct_tolerance, 1.25);
    }

    #[test]
    fn defaults_lfq_rt_percent_tolerance() {
        let options: LfqOptions = serde_json::from_str("{}").unwrap();
        let settings: LfqSettings = options.into();

        assert_eq!(settings.rt_pct_tolerance, 0.5);
    }
}

#[derive(Serialize, Deserialize, Debug)]
pub struct TmtOptions {
    pub level: Option<u8>,
    pub sn: Option<bool>,
}

#[derive(Copy, Clone, Serialize, Debug)]
pub struct TmtSettings {
    pub level: u8,
    pub sn: bool,
}

impl From<TmtOptions> for TmtSettings {
    fn from(value: TmtOptions) -> Self {
        let default = Self::default();
        Self {
            level: value.level.unwrap_or(default.level),
            sn: value.sn.unwrap_or(default.sn),
        }
    }
}

impl Default for TmtSettings {
    fn default() -> Self {
        Self {
            level: 3,
            sn: false,
        }
    }
}

#[derive(Serialize, Deserialize, Default, Debug)]
pub struct QuantOptions {
    pub tmt: Option<Isobaric>,
    #[serde(rename = "tmt_settings")]
    pub tmt_options: Option<TmtOptions>,

    pub lfq: Option<bool>,
    #[serde(rename = "lfq_settings")]
    pub lfq_options: Option<LfqOptions>,
}

#[derive(Serialize, Default, Clone)]
pub struct QuantSettings {
    pub tmt: Option<Isobaric>,
    pub tmt_settings: TmtSettings,
    pub lfq: bool,
    pub lfq_settings: LfqSettings,
}

impl From<QuantOptions> for QuantSettings {
    fn from(value: QuantOptions) -> Self {
        Self {
            tmt: value.tmt,
            tmt_settings: value.tmt_options.map(Into::into).unwrap_or_default(),

            lfq: value.lfq.unwrap_or(false),
            lfq_settings: value.lfq_options.map(Into::into).unwrap_or_default(),
        }
    }
}

impl Input {
    pub fn from_arguments(matches: ArgMatches) -> anyhow::Result<Self> {
        let path = matches
            .get_one::<String>("parameters")
            .expect("required parameters");
        let mut input = Input::load(path)
            .with_context(|| format!("Failed to read parameters from `{path}`"))?;

        // Handle JSON configuration overrides
        if let Some(output_directory) = matches.get_one::<String>("output_directory") {
            input.output_directory = Some(output_directory.into());
        }
        if let Some(fasta) = matches.get_one::<String>("fasta") {
            input.database.get_or_insert_default().fasta = Some(fasta.into());
        }
        if let Some(mzml_paths) = matches.get_many::<String>("mzml_paths") {
            input.mzml_paths = Some(mzml_paths.into_iter().map(|p| p.into()).collect());
        }

        if matches.get_flag("write-pin") {
            input.write_pin = Some(true);
        }

        if matches.get_flag("write-report") {
            input.write_report = Some(true);
        }

        if matches.get_flag("annotate-matches") {
            input.annotate_matches = Some(true);
        }
        if let Some(batch_size) = matches.get_one::<u16>("batch-size").copied() {
            input.batch_size = Some(batch_size as usize);
        }
        if let Some(max_memory_gb) = matches.get_one::<f64>("max-memory").copied() {
            input.max_memory_gb = Some(max_memory_gb);
        }

        // Only override the config-file value when the flag is explicitly set.
        if matches.get_flag("localize") {
            input.ptm_localization.get_or_insert_default().enabled = true;
        }
        if matches.get_flag("spectral-library") {
            input.spectral_library.get_or_insert_default().enabled = true;
        }

        input.validate()?;
        Ok(input)
    }

    pub fn load<S: AsRef<str>>(path: S) -> anyhow::Result<Self> {
        sage_cloudpath::util::read_json(path).map_err(anyhow::Error::from)
    }

    fn check_mass_tolerances(tolerance: &Tolerance) {
        let (lo, hi) = match tolerance {
            Tolerance::Ppm(lo, hi) => (*lo, *hi),
            Tolerance::Pct(lo, hi) => {
                log::warn!(
                    "Pct tolerances are very rarely used for mass tolerances, did you mean ppm?"
                );
                (*lo, *hi)
            }
            Tolerance::Da(lo, hi) => (*lo, *hi),
        };
        if hi.abs() > lo.abs() {
            log::warn!(
                "Tolerances are applied to experimental masses, not theoretical: [{}, {}]",
                lo,
                hi
            );
        }
        if lo > 0.0 {
            log::warn!(
                "The `left` tolerance should probably be negative, for example: [{}, {}]",
                -lo,
                hi.abs()
            )
        }
        if hi < 0.0 {
            log::warn!(
                "The `right` tolerance should probably be positive, for example: [{}, {}]",
                -lo.abs(),
                hi
            )
        }
    }

    /// Validate logical configuration constraints without reading inputs or writing outputs.
    pub fn validate(&self) -> anyhow::Result<()> {
        ensure!(
            self.database.is_some() ^ self.library_search.is_some(),
            "exactly one of `database` or `library_search` must be configured"
        );
        if let Some(database) = &self.database {
            ensure!(
                database.fasta.is_some() || database.peptides.is_some(),
                "Either `database.fasta` or `database.peptides` must be set"
            );
            ensure!(
                database.custom_cleavage_sites.is_none() || database.fasta.is_some(),
                "`database.custom_cleavage_sites` requires `database.fasta`"
            );
        }
        if let Some(library) = &self.library_search {
            library.validate().map_err(anyhow::Error::msg)?;
            ensure!(
                self.report_psms.unwrap_or(1) == 1,
                "`report_psms` must be 1 with `library_search` for target-decoy competition"
            );
            ensure!(
                !self.chimera.unwrap_or(false),
                "`chimera` is not supported with `library_search`"
            );
            ensure!(
                !self.wide_window.unwrap_or(false),
                "`wide_window` is not supported with `library_search`"
            );
            ensure!(
                !self.use_bitmap.unwrap_or(false),
                "`use_bitmap` is not supported with `library_search`"
            );
            ensure!(
                !self.ptm_localization.unwrap_or_default().enabled,
                "`ptm_localization` is not supported with `library_search`"
            );
            ensure!(
                !self
                    .spectral_library
                    .as_ref()
                    .is_some_and(|settings| settings.enabled),
                "spectral-library export is not supported during `library_search`"
            );
        }
        ensure!(
            self.mzml_paths.as_ref().map(Vec::len).unwrap_or_default() > 0,
            "`mzml_paths` must contain at least one spectra file"
        );
        if let Some((low, high)) = self.isotope_errors {
            ensure!(
                low <= high,
                "isotope errors must be specified [low, high], received [{low}, {high}]"
            );
        }
        if let Some((low, high)) = self.precursor_charge {
            ensure!(
                low <= high,
                "precursor charges must be specified [low, high], received [{low}, {high}]"
            );
        }
        if let (Some(min), Some(max)) = (self.min_peaks, self.max_peaks) {
            ensure!(min <= max, "`min_peaks` cannot exceed `max_peaks`");
        }
        ensure!(
            self.report_psms.unwrap_or(1) > 0,
            "`report_psms` must be greater than zero"
        );
        let output_psm_q_value = self.output_filter.unwrap_or_default().psm_q_value;
        ensure!(
            output_psm_q_value.is_finite() && (0.0..=1.0).contains(&output_psm_q_value),
            "output_filter.psm_q_value must be between 0 and 1"
        );
        self.memory_limits()?;
        resolve_batch_size(self.batch_size)?;
        if let Some(rt_pct_tolerance) = self
            .quant
            .as_ref()
            .and_then(|quant| quant.lfq_options.as_ref())
            .and_then(|lfq| lfq.rt_pct_tolerance)
        {
            ensure!(
                rt_pct_tolerance != 0.0,
                "`lfq_settings.rt_pct_tolerance` must not be zero"
            );
        }
        if let Some(settings) = self.ptm_localization {
            ensure!(
                settings.psm_q_value.is_finite() && (0.0..=1.0).contains(&settings.psm_q_value),
                "ptm_localization.psm_q_value must be between 0 and 1"
            );
            ensure!(
                settings.localization_q_value.is_finite()
                    && (0.0..=1.0).contains(&settings.localization_q_value),
                "ptm_localization.localization_q_value must be between 0 and 1"
            );
        }
        if let Some(settings) = &self.spectral_library {
            settings.validate().map_err(anyhow::Error::msg)?;
        }
        Ok(())
    }

    pub fn build(self) -> anyhow::Result<Search> {
        self.validate()?;
        let memory_limits = self.memory_limits()?;
        let batch_size = resolve_batch_size(self.batch_size)?;
        let database = self.database.map(Builder::make_parameters);
        if let Some(database) = &database {
            database
                .validate_ptm_library(&sage_core::ptm_library::PtmLibrary::default())
                .map_err(anyhow::Error::msg)?;
        }
        let library_mode = self.library_search.is_some();

        Self::check_mass_tolerances(&self.fragment_tol);
        Self::check_mass_tolerances(&self.precursor_tol);

        let mzml_paths = self
            .mzml_paths
            .expect("'mzml_paths' must be provided!")
            .iter()
            .map(|s| sage_cloudpath::to_url(s))
            .collect::<Result<Vec<_>, _>>()?;

        let output_directory = match self.output_directory {
            Some(path) => {
                match sage_cloudpath::try_parse_url(&path) {
                    Some(mut url) => {
                        // Valid URL, might still be a local directory that doesn't exist
                        if url.scheme() == "file" {
                            let path = url.to_file_path().expect("url scheme is file");
                            std::fs::create_dir_all(path)?;
                        }

                        if !url.path().ends_with("/") {
                            url.set_path(&format!("{}/", url.path()));
                        }
                        url
                    }
                    None => {
                        // Treat as a local path (covers Windows `C:\...` which
                        // otherwise parses as a URL with scheme `c`).
                        let path = std::path::Path::new(&path);
                        std::fs::create_dir_all(path)?;
                        Url::from_directory_path(path.canonicalize()?).expect("valid path")
                    }
                }
            }
            None => {
                let dir = std::env::current_dir()?;
                Url::from_directory_path(dir).expect("valid path")
            }
        };

        let score_type = self.score_type.unwrap_or(ScoreType::SageHyperScore);

        let ptm_localization = self.ptm_localization.unwrap_or_default();
        ensure!(
            ptm_localization.psm_q_value.is_finite()
                && (0.0..=1.0).contains(&ptm_localization.psm_q_value),
            "ptm_localization.psm_q_value must be between 0 and 1"
        );
        ensure!(
            ptm_localization.localization_q_value.is_finite()
                && (0.0..=1.0).contains(&ptm_localization.localization_q_value),
            "ptm_localization.localization_q_value must be between 0 and 1"
        );
        let spectral_library = self.spectral_library.unwrap_or_default();
        spectral_library.validate().map_err(anyhow::Error::msg)?;

        Ok(Search {
            version: clap::crate_version!().into(),
            database,
            library_search: self.library_search,
            quant: self.quant.map(Into::into).unwrap_or_default(),
            mzml_paths,
            output_directory,
            precursor_tol: self.precursor_tol,
            fragment_tol: self.fragment_tol,
            report_psms: self.report_psms.unwrap_or(1),
            output_filter: self.output_filter.unwrap_or_default(),
            max_peaks: self.max_peaks.unwrap_or(150),
            min_peaks: self.min_peaks.unwrap_or(15),
            min_matched_peaks: self.min_matched_peaks.unwrap_or(4),
            max_fragment_charge: self.max_fragment_charge,
            annotate_matches: self.annotate_matches.unwrap_or(false),
            precursor_charge: self.precursor_charge.unwrap_or((2, 4)),
            override_precursor_charge: self.override_precursor_charge.unwrap_or(false),
            isotope_errors: self.isotope_errors.unwrap_or((0, 0)),
            deisotope: self.deisotope.unwrap_or(true),
            chimera: self.chimera.unwrap_or(false),
            wide_window: self.wide_window.unwrap_or(false),
            predict_rt: self.predict_rt.unwrap_or(!library_mode),
            retention_time_model: self.retention_time_model.unwrap_or_default(),
            retention_time_alignment: self.retention_time_alignment,
            ion_mobility_model: self.ion_mobility_model.unwrap_or_default(),
            output_paths: Vec::new(),
            write_pin: self.write_pin.unwrap_or(false),
            bruker_config: self.bruker_config.unwrap_or_default(),
            write_report: self.write_report.unwrap_or(false),
            protein_grouping: self.protein_grouping.unwrap_or(true),
            protein_grouping_peptide_fdr: self.protein_grouping_peptide_fdr.unwrap_or(0.01),
            max_memory_gb: memory_limits.max_gib(),
            min_free_memory_gb: memory_limits.min_free_gib(),
            batch_size,
            ptm_localization,
            spectral_library,
            mass_shift_ppm: self
                .mass_shift_ppm
                .unwrap_or(sage_core::ambiguity::DEFAULT_MASS_SHIFT_PPM),
            score_type,
            use_bitmap: self.use_bitmap.unwrap_or(false),
        })
    }

    /// Validate and convert the configured memory limits.
    pub fn memory_limits(&self) -> anyhow::Result<MemoryLimits> {
        MemoryLimits::from_gib(self.max_memory_gb, self.min_free_memory_gb)
    }
}

fn resolve_batch_size(batch_size: Option<usize>) -> anyhow::Result<usize> {
    match batch_size {
        Some(batch_size) => {
            ensure!(batch_size > 0, "`batch_size` must be greater than zero");
            Ok(batch_size)
        }
        None => Ok((num_cpus::get() / 2).max(1)),
    }
}

#[cfg(test)]
mod test {

    use super::{resolve_batch_size, Input, OutputFilter, PtmLocalizationSettings};
    use sage_core::{
        database::EnzymeBuilder,
        enzyme::EnzymeParameters,
        ml::retention_alignment::AlignmentMethod,
        ml::retention_model::{RetentionTimeFeatureSet, RetentionTimeSettings},
        spectral_library::{
            SpectralLibraryFormat, SpectralLibrarySettings, SpectralLibraryStrategy,
        },
    };

    #[test]
    fn deserialize_enriched_retention_time_settings() -> Result<(), serde_json::Error> {
        let settings: RetentionTimeSettings = serde_json::from_value(serde_json::json!({
            "features": "additive_ptm",
            "folds": 5,
            "seed": 7,
            "ptm_regularization": 12.5
        }))?;

        assert_eq!(settings.features, RetentionTimeFeatureSet::AdditivePtm);
        assert_eq!(settings.folds, 5);
        assert_eq!(settings.seed, 7);
        assert_eq!(settings.ptm_regularization, 12.5);
        Ok(())
    }

    #[test]
    fn deserialize_ptm_localization_settings() -> Result<(), serde_json::Error> {
        let configured: PtmLocalizationSettings = serde_json::from_value(serde_json::json!({
            "enabled": true,
            "psm_q_value": 0.025,
            "localization_q_value": 0.05
        }))?;
        assert!(configured.enabled);
        assert_eq!(configured.psm_q_value, 0.025);
        assert_eq!(configured.localization_q_value, 0.05);

        let partial: PtmLocalizationSettings =
            serde_json::from_value(serde_json::json!({ "enabled": true }))?;
        assert!(partial.enabled);
        assert_eq!(partial.psm_q_value, 0.01);
        assert_eq!(partial.localization_q_value, 0.01);

        let legacy_name: PtmLocalizationSettings =
            serde_json::from_value(serde_json::json!({ "q_value": 0.02 }))?;
        assert_eq!(legacy_name.psm_q_value, 0.02);
        Ok(())
    }

    #[test]
    fn deserialize_spectral_library_settings() -> Result<(), serde_json::Error> {
        let configured: SpectralLibrarySettings = serde_json::from_value(serde_json::json!({
            "enabled": true,
            "strategy": "consensus",
            "max_fragments": 12,
            "min_consensus_psms": 2,
            "min_fragment_frequency": 0.6,
            "formats": ["sage_parquet", "mzspeclib"]
        }))?;
        assert!(configured.enabled);
        assert_eq!(configured.psm_q_value, 0.01);
        assert_eq!(configured.peptide_q_value, 0.01);
        assert_eq!(configured.strategy, SpectralLibraryStrategy::Consensus);
        assert_eq!(configured.max_fragments, 12);
        assert_eq!(configured.min_consensus_psms, 2);
        assert_eq!(configured.min_fragment_frequency, 0.6);
        assert_eq!(
            configured.formats,
            vec![
                SpectralLibraryFormat::SageParquet,
                SpectralLibraryFormat::MzSpecLib
            ]
        );
        Ok(())
    }

    #[test]
    fn spectral_library_settings_are_validated() {
        let input: Input = serde_json::from_value(serde_json::json!({
            "database": { "fasta": "test.fasta" },
            "precursor_tol": { "ppm": [-10, 10] },
            "fragment_tol": { "ppm": [-10, 10] },
            "mzml_paths": ["test.mzML"],
            "spectral_library": { "enabled": true, "max_fragments": 0 }
        }))
        .unwrap();
        let error = input.validate().unwrap_err().to_string();
        assert!(error.contains("spectral_library.max_fragments"));
    }

    fn base_search_space(value: serde_json::Value) -> Input {
        serde_json::from_value(serde_json::json!({
            "precursor_tol": { "ppm": [-10, 10] },
            "fragment_tol": { "ppm": [-10, 10] },
            "mzml_paths": ["test.mzML"],
            "database": value
        }))
        .unwrap()
    }

    #[test]
    fn database_and_library_search_are_mutually_exclusive() {
        let database = base_search_space(serde_json::json!({ "fasta": "test.fasta" }));
        assert!(database.validate().is_ok());

        let library: Input = serde_json::from_value(serde_json::json!({
            "library_search": { "path": "library.mzspeclib.txt" },
            "precursor_tol": { "ppm": [-10, 10] },
            "fragment_tol": { "ppm": [-10, 10] },
            "isotope_errors": [-1, 2],
            "mzml_paths": ["tests/LQSRPAAPPAPGPGQLTLR.mzML"]
        }))
        .unwrap();
        assert!(library.validate().is_ok());

        let both: Input = serde_json::from_value(serde_json::json!({
            "database": { "fasta": "test.fasta" },
            "library_search": { "path": "library.sage.parquet" },
            "precursor_tol": { "ppm": [-10, 10] },
            "fragment_tol": { "ppm": [-10, 10] },
            "mzml_paths": ["test.mzML"]
        }))
        .unwrap();
        assert!(both
            .validate()
            .unwrap_err()
            .to_string()
            .contains("exactly one"));

        let neither: Input = serde_json::from_value(serde_json::json!({
            "precursor_tol": { "ppm": [-10, 10] },
            "fragment_tol": { "ppm": [-10, 10] },
            "mzml_paths": ["test.mzML"]
        }))
        .unwrap();
        assert!(neither
            .validate()
            .unwrap_err()
            .to_string()
            .contains("exactly one"));
    }

    #[test]
    fn library_search_build_does_not_create_database_parameters() {
        let spectra = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("../../tests/LQSRPAAPPAPGPGQLTLR.mzML");
        let input: Input = serde_json::from_value(serde_json::json!({
            "library_search": { "path": "library.sage.parquet" },
            "precursor_tol": { "ppm": [-10, 10] },
            "fragment_tol": { "ppm": [-10, 10] },
            "isotope_errors": [-1, 2],
            "mzml_paths": [spectra]
        }))
        .unwrap();
        let search = input.build().unwrap();
        assert!(search.database.is_none());
        assert!(search.library_search.is_some());
        assert_eq!(search.isotope_errors, (-1, 2));
        assert!(!search.predict_rt);
    }

    #[test]
    fn output_filter_defaults_and_deserializes() -> Result<(), serde_json::Error> {
        let default: OutputFilter = serde_json::from_value(serde_json::json!({}))?;
        assert_eq!(default.psm_q_value, 0.1);

        let configured: OutputFilter =
            serde_json::from_value(serde_json::json!({ "psm_q_value": 0.025 }))?;
        assert_eq!(configured.psm_q_value, 0.025);
        Ok(())
    }

    #[test]
    fn output_filter_q_value_must_be_a_probability() {
        let input: Input = serde_json::from_value(serde_json::json!({
            "database": { "fasta": "test.fasta" },
            "precursor_tol": { "ppm": [-10, 10] },
            "fragment_tol": { "ppm": [-10, 10] },
            "mzml_paths": ["test.mzML"],
            "output_filter": { "psm_q_value": 1.1 }
        }))
        .unwrap();

        let error = input.validate().unwrap_err().to_string();
        assert!(error.contains("output_filter.psm_q_value"));
    }

    #[test]
    fn deserialize_enzyme_builder() -> Result<(), serde_json::Error> {
        let a: EnzymeBuilder = serde_json::from_value(serde_json::json!({
            "cleave_at": "KR",
        }))?;
        let b: EnzymeBuilder = serde_json::from_value(serde_json::json!({
            "cleave_at": "KR",
            "restrict": "P",
        }))?;
        let c: EnzymeBuilder = serde_json::from_value(serde_json::json!({
            "cleave_at": "KR",
            "restrict": "",
        }))?;

        let a: EnzymeParameters = a.into();
        let b: EnzymeParameters = b.into();
        let c: EnzymeParameters = c.into();

        assert_eq!(a.enzyme.map(|e| e.skip_suffix), Some([false; 26]));
        {
            let mut expected = [false; 26];
            expected[(b'P' - b'A') as usize] = true;
            assert_eq!(b.enzyme.map(|e| e.skip_suffix), Some(expected));
        }
        assert_eq!(c.enzyme.map(|e| e.skip_suffix), Some([false; 26]));

        Ok(())
    }

    #[test]
    fn deserialize_custom_cleavage_site_path() -> Result<(), serde_json::Error> {
        let input: Input = serde_json::from_value(serde_json::json!({
            "database": {
                "fasta": "proteome.fasta",
                "custom_cleavage_sites": "cleavage-sites.tsv"
            },
            "precursor_tol": { "ppm": [-10.0, 10.0] },
            "fragment_tol": { "ppm": [-20.0, 20.0] },
            "mzml_paths": ["input.mzML"]
        }))?;

        assert_eq!(
            input
                .database
                .as_ref()
                .and_then(|database| database.custom_cleavage_sites.as_deref()),
            Some("cleavage-sites.tsv")
        );
        assert!(input.validate().is_ok());
        Ok(())
    }

    #[test]
    fn deserialize_runtime_memory_settings() -> Result<(), serde_json::Error> {
        let input: Input = serde_json::from_value(serde_json::json!({
            "database": {},
            "precursor_tol": { "ppm": [-10.0, 10.0] },
            "fragment_tol": { "ppm": [-20.0, 20.0] },
            "max_memory_gb": 12.5,
            "min_free_memory_gb": 2.0,
            "batch_size": 1
        }))?;

        assert_eq!(input.max_memory_gb, Some(12.5));
        assert_eq!(input.min_free_memory_gb, Some(2.0));
        assert_eq!(input.batch_size, Some(1));
        assert!(input.memory_limits().unwrap().is_enabled());
        Ok(())
    }

    #[test]
    fn batch_size_must_be_positive() {
        assert!(resolve_batch_size(Some(0)).is_err());
        assert_eq!(resolve_batch_size(Some(3)).unwrap(), 3);
        assert!(resolve_batch_size(None).unwrap() >= 1);
    }

    #[test]
    fn deserialize_nonlinear_retention_time_alignment() -> Result<(), serde_json::Error> {
        let method: AlignmentMethod = serde_json::from_value(serde_json::json!("nonlinear"))?;
        assert_eq!(method, AlignmentMethod::Nonlinear);
        Ok(())
    }

    #[test]
    fn validation_returns_range_errors_instead_of_exiting() {
        let input: Input = serde_json::from_value(serde_json::json!({
            "database": { "fasta": "test.fasta" },
            "precursor_tol": { "ppm": [-10, 10] },
            "fragment_tol": { "ppm": [-10, 10] },
            "isotope_errors": [3, -1],
            "mzml_paths": ["test.mzML"]
        }))
        .unwrap();

        let error = input.validate().unwrap_err().to_string();
        assert!(error.contains("isotope errors"));
    }
}
