use super::input::Search;
use super::memory::MemoryLimits;
use super::output::SageResults;
use super::telemetry;
use crate::events::{CancellationToken, EventEmitter, EventKind};
use anyhow::Context;
use log::{info, warn};
use rayon::prelude::*;
use sage_cloudpath::{FileFormat, Url};
use sage_core::cleavage::{CustomCleavageLibrary, ValidatedCustomCleavageLibrary};
use sage_core::database::{Builder, IndexedDatabase, Parameters};
use sage_core::fasta::Fasta;
use sage_core::lfq::{PrecursorId, QuantifiedPeak};
use sage_core::mass::Tolerance;
use sage_core::mass_calibration::{
    align_fragment_error, fit as fit_mass_calibration, CalibrationPoint, FitOptions,
};
use sage_core::peptide::Peptide;
use sage_core::scoring::{Feature, Scorer};
use sage_core::spectral_library::{
    self, LibrarySelection, SpectralLibraryFormat, SpectralLibraryStrategy,
};
use sage_core::spectrum::{ProcessedSpectrum, SpectrumProcessor};
use std::cmp::Ordering;
use std::collections::{HashMap, HashSet};
use std::io::{BufWriter, Write};
use std::sync::Arc;
use std::time::Instant;
// HTML report specific imports
use maud::{html, PreEscaped};
use report_builder::{
    plots::{plot_boxplot, plot_pp, plot_scatter, plot_score_histogram},
    Report, ReportSection,
};

mod library_runner;

enum OutputTarget {
    Local(BufWriter<std::fs::File>),
    Remote(Box<BufWriter<sage_cloudpath::CloudWriter>>),
}

impl OutputTarget {
    fn new(path: &Url) -> anyhow::Result<Self> {
        if let Ok(local_path) = path.to_file_path() {
            if let Some(parent) = local_path.parent() {
                std::fs::create_dir_all(parent)?;
            }
            Ok(Self::Local(BufWriter::new(std::fs::File::create(
                local_path,
            )?)))
        } else {
            Ok(Self::Remote(Box::new(BufWriter::with_capacity(
                1024 * 1024,
                sage_cloudpath::CloudWriter::new(path)?,
            ))))
        }
    }

    fn finish(mut self, _path: &Url) -> anyhow::Result<()> {
        self.flush()?;
        if let Self::Remote(writer) = self {
            (*writer)
                .into_inner()
                .map_err(|error| error.into_error())?
                .finish()?;
        }
        Ok(())
    }
}

impl Write for OutputTarget {
    fn write(&mut self, buf: &[u8]) -> std::io::Result<usize> {
        match self {
            Self::Local(writer) => writer.write(buf),
            Self::Remote(writer) => writer.write(buf),
        }
    }

    fn flush(&mut self) -> std::io::Result<()> {
        match self {
            Self::Local(writer) => writer.flush(),
            Self::Remote(writer) => writer.flush(),
        }
    }
}

fn finish_csv_writer(mut writer: csv::Writer<OutputTarget>, path: &Url) -> anyhow::Result<()> {
    writer.flush()?;
    let output = writer
        .into_inner()
        .map_err(|error| anyhow::anyhow!("failed to flush CSV output: {}", error.error()))?;
    output.finish(path)
}

pub struct Runner {
    pub database: IndexedDatabase,
    pub parameters: Search,
    database_parameters: Parameters,
    library_search: Option<library_runner::LibrarySearchRuntime>,
    start: Instant,
    events: EventEmitter,
    cancellation: CancellationToken,
}

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
    #[serde(default)]
    pub library_search: LibrarySearchRunStats,
    pub output_paths: Vec<String>,
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
}

#[derive(Debug, Clone, Default, serde::Serialize, serde::Deserialize)]
pub struct ModelRunStats {
    pub mass_alignment_applied: bool,
    pub retention_time_prediction_enabled: bool,
    pub retention_time_model_fitted: bool,
    pub retention_time_features: String,
    pub retention_time_alignment: Option<String>,
    pub ion_mobility_observed: bool,
    pub ion_mobility_model_enabled: bool,
    pub ion_mobility_model_fitted: bool,
    pub ion_mobility_features: String,
    /// Library-reference RT alignment, separate from database RT prediction.
    pub library_retention_time_alignment: Option<String>,
    pub library_retention_time_files_aligned: usize,
    /// Library-reference ion-mobility alignment, separate from sequence models.
    pub library_ion_mobility_alignment: Option<String>,
    pub library_ion_mobility_files_aligned: usize,
    /// Final scoring model used for library-search PSMs.
    pub library_rescoring: Option<String>,
}

#[derive(Debug, Clone, Default, serde::Serialize, serde::Deserialize)]
pub struct QuantificationRunStats {
    pub lfq_enabled: bool,
    pub lfq_features: usize,
    pub tmt: Option<String>,
    pub tmt_features: usize,
    #[serde(default)]
    pub ms1_label_channels: usize,
    #[serde(default)]
    pub ms1_label_reference: Option<String>,
}

#[derive(Debug, Clone, Default, serde::Serialize, serde::Deserialize)]
pub struct ExecutionRunStats {
    pub batch_size: usize,
    pub parallelism: usize,
    pub bitmap_search: bool,
    pub max_memory_gb: Option<f64>,
    pub min_free_memory_gb: Option<f64>,
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

#[derive(Debug, Clone, Default, serde::Serialize, serde::Deserialize)]
pub struct LibrarySearchRunStats {
    pub enabled: bool,
    pub path: Option<String>,
    pub target_entries: usize,
    pub decoy_entries: usize,
    pub transitions: usize,
}

/// A single localized modification site for one PSM, used to build the
/// PTM-site and protein-site reports.
struct SiteRow {
    psm_id: usize,
    filename: String,
    scannr: String,
    peptide: String,
    peptide_sequence: String,
    proteins: String,
    charge: u8,
    spectrum_q: f32,
    peptide_q: f32,
    modification: String,
    modification_mass: f32,
    /// 1-based position within the peptide.
    position: usize,
    residue: u8,
    localization_probability: f32,
    delta_score: f32,
    target_decoy_score: f32,
    localization_q_value: f32,
    candidate_sites: usize,
    site_determining_matched: u32,
    site_determining_total: u32,
    site_probabilities: String,
}

fn passes_localization_filter(feature: &Feature, psm_q_value: f32) -> bool {
    feature.label == 1 && feature.spectrum_q <= psm_q_value
}

fn passes_output_filter(feature: &Feature, psm_q_value: f32) -> bool {
    feature.spectrum_q <= psm_q_value
}

fn feature_identity_cmp(left: &Feature, right: &Feature) -> Ordering {
    left.file_id
        .cmp(&right.file_id)
        .then_with(|| left.spec_id.cmp(&right.spec_id))
        .then_with(|| left.rank.cmp(&right.rank))
        .then_with(|| left.peptide_idx.cmp(&right.peptide_idx))
        .then_with(|| left.charge.cmp(&right.charge))
}

fn sort_features_by_discriminant(features: &mut [Feature]) {
    features.par_sort_unstable_by(|left, right| {
        right
            .discriminant_score
            .total_cmp(&left.discriminant_score)
            .then_with(|| feature_identity_cmp(left, right))
    });
}

fn assign_psm_ids(features: &mut [Feature]) {
    for (index, feature) in features.iter_mut().enumerate() {
        feature.psm_id = index + 1;
    }
}

#[derive(Default)]
struct SpectrumAccumulator {
    pub ms1: Vec<ProcessedSpectrum>,
    pub msn: Vec<ProcessedSpectrum>,
}

#[derive(Default)]
struct PostprocessStats {
    annotated_psms: usize,
    annotated_fragments: usize,
    localized_psms: usize,
    library_selections: Vec<LibrarySelection>,
}

impl SpectrumAccumulator {
    pub fn fold_op(mut self, rhs: ProcessedSpectrum) -> Self {
        if rhs.level == 1 {
            self.ms1.push(rhs);
        } else {
            self.msn.push(rhs);
        }
        self
    }

    pub fn reduce(mut self, other: Self) -> Self {
        self.ms1.extend(other.ms1);
        self.msn.extend(other.msn);
        self
    }
}

impl FromParallelIterator<ProcessedSpectrum> for SpectrumAccumulator {
    fn from_par_iter<I>(par_iter: I) -> Self
    where
        I: IntoParallelIterator<Item = ProcessedSpectrum>,
    {
        par_iter
            .into_par_iter()
            .fold(SpectrumAccumulator::default, SpectrumAccumulator::fold_op)
            .reduce(SpectrumAccumulator::default, SpectrumAccumulator::reduce)
    }
}

impl FromIterator<ProcessedSpectrum> for SpectrumAccumulator {
    fn from_iter<I>(iter: I) -> Self
    where
        I: IntoIterator<Item = ProcessedSpectrum>,
    {
        iter.into_iter()
            .fold(SpectrumAccumulator::default(), SpectrumAccumulator::fold_op)
    }
}

impl Runner {
    pub fn new(parameters: Search, parallel: usize) -> anyhow::Result<Self> {
        Self::new_with_control(
            parameters,
            parallel,
            EventEmitter::disabled(),
            CancellationToken::default(),
        )
    }

    pub fn new_with_control(
        parameters: Search,
        parallel: usize,
        events: EventEmitter,
        cancellation: CancellationToken,
    ) -> anyhow::Result<Self> {
        let parameters = parameters.clone();
        let mut database_parameters = parameters
            .database
            .clone()
            .unwrap_or_else(|| Builder::default().make_parameters());
        let start = Instant::now();
        cancellation.check()?;
        if let Some(settings) = parameters.library_search.as_ref() {
            events.emit(EventKind::LibrarySearchStarted {
                path: settings.path.clone(),
            });
            let (database, library_search) = library_runner::load_library_search(settings)?;
            let overlapping_sources =
                library_search.overlapping_source_files(&parameters.mzml_paths);
            if !overlapping_sources.is_empty() {
                let message = format!(
                    "library search input overlaps recorded library source file(s): {}; same-source searches are unsuitable for FDR validation",
                    overlapping_sources.join(", ")
                );
                warn!("{message}");
                events.emit(EventKind::Warning {
                    code: "library_source_overlap".into(),
                    message,
                });
            }
            events.emit(EventKind::LibrarySearchBuilt {
                target_entries: library_search.target_entries,
                decoy_entries: library_search.target_entries,
                transitions: library_search.target_transitions,
            });
            events.check()?;
            return Ok(Self {
                database,
                parameters,
                database_parameters,
                library_search: Some(library_search),
                start,
                events,
                cancellation,
            });
        }
        database_parameters.use_bitmap = parameters.use_bitmap;
        if let Some(settings) = database_parameters.ptm_library.clone() {
            let library = if sage_core::ptm_library::is_tsv_path(&settings.path) {
                let contents = sage_cloudpath::util::read_text(&settings.path)
                    .with_context(|| format!("Failed to read PTM library `{}`", settings.path))?;
                sage_core::ptm_library::PtmLibrary::from_tsv(&contents).map_err(anyhow::Error::msg)
            } else {
                let bytes = sage_cloudpath::util::read_bytes(&settings.path)
                    .with_context(|| format!("Failed to read PTM library `{}`", settings.path))?;
                sage_cloudpath::parquet::deserialize_ptm_library(bytes).map_err(anyhow::Error::from)
            }
            .with_context(|| format!("Failed to parse PTM library `{}`", settings.path))?;
            database_parameters
                .validate_ptm_library(&library)
                .map_err(anyhow::Error::msg)?;
            info!("loaded {} unique PTM library sites", library.len());
            database_parameters.loaded_ptm_library = Some(Arc::new(library));
        }
        events.emit(EventKind::DatabaseStarted);
        let limits =
            MemoryLimits::from_gib(parameters.max_memory_gb, parameters.min_free_memory_gb)?;
        // Collect peptides from FASTA (if configured).
        let mut all_peptides: Vec<Peptide> = if !database_parameters.fasta.is_empty() {
            let fasta_url = sage_cloudpath::to_url(&database_parameters.fasta)?;
            let fasta = sage_cloudpath::util::read_fasta(
                &fasta_url,
                &database_parameters.decoy_tag,
                database_parameters.generate_decoys,
            )
            .with_context(|| {
                format!(
                    "Failed to build database from `{}`",
                    database_parameters.fasta
                )
            })?;
            let custom_cleavages = if let Some(path) =
                database_parameters.custom_cleavage_sites.as_deref()
            {
                let library = if path.to_ascii_lowercase().ends_with(".parquet") {
                    let content = sage_cloudpath::util::read_bytes(path).with_context(|| {
                        format!("Failed to read custom cleavage-site file `{path}`")
                    })?;
                    sage_cloudpath::parquet::deserialize_custom_cleavage_sites(content)
                        .with_context(|| {
                            format!("Failed to parse custom cleavage-site file `{path}`")
                        })?
                } else {
                    let content = sage_cloudpath::util::read_text(path).with_context(|| {
                        format!("Failed to read custom cleavage-site file `{path}`")
                    })?;
                    CustomCleavageLibrary::from_tsv(&content).with_context(|| {
                        format!("Failed to parse custom cleavage-site file `{path}`")
                    })?
                };
                let validated = library.validate(&fasta).with_context(|| {
                    format!("Failed to validate custom cleavage-site file `{path}`")
                })?;
                info!(
                    "custom cleavage sites: {} matched, {} unmatched",
                    validated.matched_sites, validated.unmatched_sites
                );
                if validated.unmatched_sites > 0 {
                    warn!(
                        "{} custom cleavage sites refer to proteins absent from the FASTA",
                        validated.unmatched_sites
                    );
                }
                if validated.sites_without_context > 0 {
                    warn!(
                        "{} custom cleavage sites have no sequence context and were validated by coordinate only",
                        validated.sites_without_context
                    );
                }
                Some(validated)
            } else {
                None
            };

            if let (Some(settings), Some(library)) = (
                database_parameters.ptm_library.as_ref(),
                database_parameters.loaded_ptm_library.as_deref(),
            ) {
                let proteins = fasta
                    .targets
                    .iter()
                    .map(|(accession, sequence)| (accession.as_ref(), sequence.as_bytes()))
                    .collect::<HashMap<_, _>>();
                let mut invalid = Vec::new();
                for site in library.iter() {
                    match proteins.get(site.protein.as_ref()) {
                        None => invalid.push(format!(
                            "{}:{} references a protein absent from the FASTA",
                            site.protein,
                            site.position + 1
                        )),
                        Some(sequence)
                            if sequence.get(site.position as usize) != Some(&site.residue) =>
                        {
                            invalid.push(format!(
                                "{}:{} expects residue {}",
                                site.protein,
                                site.position + 1,
                                site.residue as char
                            ));
                        }
                        _ => {}
                    }
                }
                if !invalid.is_empty() {
                    let message = format!(
                        "{} PTM library sites did not match the FASTA; first issue: {}",
                        invalid.len(),
                        invalid[0]
                    );
                    if settings.strict {
                        anyhow::bail!(message);
                    }
                    log::warn!("{message}");
                }
            }

            let needs_estimate = limits.is_enabled()
                || (database_parameters.prefilter && database_parameters.prefilter_chunk_size == 0);
            if needs_estimate {
                let full_estimate = database_parameters
                    .estimate_memory_with_custom_cleavages(&fasta, custom_cleavages.as_ref());
                events.emit(EventKind::DatabaseEstimated {
                    unmodified_peptides: full_estimate.unmodified_peptides,
                    modified_peptides: full_estimate.modified_peptides,
                    fragments: full_estimate.fragments,
                    peak_bytes: full_estimate
                        .unmodified_peak_bytes
                        .max(full_estimate.modified_peak_bytes)
                        .max(full_estimate.fragment_peak_bytes),
                });
                if database_parameters.prefilter && database_parameters.prefilter_chunk_size == 0 {
                    database_parameters.auto_calculate_prefilter_chunk_size(
                        &fasta,
                        full_estimate.modified_peptides,
                    );
                }

                if limits.is_enabled() {
                    info!(
                        "database preflight: {} unmodified peptides ({:.2} GiB peak), up to {} modified peptides ({:.2} GiB peak), up to {} fragments ({:.2} GiB index peak)",
                        full_estimate.unmodified_peptides,
                        full_estimate.unmodified_peak_bytes as f64 / (1024.0 * 1024.0 * 1024.0),
                        full_estimate.modified_peptides,
                        full_estimate.modified_peak_bytes as f64 / (1024.0 * 1024.0 * 1024.0),
                        full_estimate.fragments,
                        full_estimate.fragment_peak_bytes as f64 / (1024.0 * 1024.0 * 1024.0),
                    );
                    limits.check_estimate(
                        "unmodified-peptide",
                        full_estimate.unmodified_peak_bytes,
                    )?;

                    if database_parameters.prefilter {
                        let mut modified_peak = 0u64;
                        let mut fragment_peak = 0u64;
                        for chunk in fasta.iter_chunks(database_parameters.prefilter_chunk_size) {
                            let estimate = database_parameters
                                .estimate_memory_with_custom_cleavages(
                                    &chunk,
                                    custom_cleavages.as_ref(),
                                );
                            modified_peak = modified_peak.max(estimate.modified_peak_bytes);
                            fragment_peak = fragment_peak.max(estimate.fragment_peak_bytes);
                        }
                        limits.check_estimate("modified-peptide", modified_peak)?;
                        limits.check_estimate("fragment-index", fragment_peak)?;
                    }
                }
            }

            match database_parameters.prefilter {
                false => {
                    let digests = database_parameters
                        .digest_unmodified_with_custom_cleavages(&fasta, custom_cleavages.as_ref());
                    if limits.is_enabled() {
                        let estimate = database_parameters.estimate_modified_memory(&digests);
                        info!(
                            "modification preflight: {} unmodified peptides may expand to {} modified peptides ({:.2} GiB additional peak)",
                            estimate.unmodified_peptides,
                            estimate.modified_peptides,
                            estimate.modified_peak_bytes as f64 / (1024.0 * 1024.0 * 1024.0),
                        );
                        limits.check_estimate("modified-peptide", estimate.modified_peak_bytes)?;
                    }
                    database_parameters.modify_digests(digests)
                }
                true => {
                    if database_parameters.prefilter_chunk_size >= fasta.targets.len() {
                        database_parameters
                            .digest_with_custom_cleavages(&fasta, custom_cleavages.as_ref())
                    } else {
                        info!(
                            "using {} db chunks of size {}",
                            (fasta.targets.len() + database_parameters.prefilter_chunk_size - 1)
                                / database_parameters.prefilter_chunk_size,
                            database_parameters.prefilter_chunk_size,
                        );
                        let mini_runner = Self {
                            database: IndexedDatabase::default(),
                            parameters: parameters.clone(),
                            database_parameters: database_parameters.clone(),
                            library_search: None,
                            start,
                            events: events.clone(),
                            cancellation: cancellation.clone(),
                        };
                        mini_runner.prefilter_peptides(parallel, fasta, custom_cleavages)
                    }
                }
            }
        } else {
            if database_parameters.loaded_ptm_library.is_some() {
                anyhow::bail!("database.ptm_library requires database.fasta");
            }
            vec![]
        };

        // Append peptides from TSV file (if configured), additive with FASTA.
        if let Some(peptides_path) = database_parameters.peptides.clone() {
            let content = sage_cloudpath::util::read_text(&peptides_path)
                .with_context(|| format!("Failed to read peptide file `{peptides_path}`"))?;
            all_peptides.extend(database_parameters.peptides_from_tsv(&content));
        }

        // Merge, deduplicate, and build the index.
        Parameters::reorder_peptides(&mut all_peptides);
        if limits.is_enabled() {
            let estimate = database_parameters.estimate_index_memory(&all_peptides);
            info!(
                "final database preflight: {} peptides, {} fragments, estimated {:.2} GiB peak",
                estimate.modified_peptides,
                estimate.fragments,
                estimate.fragment_peak_bytes as f64 / (1024.0 * 1024.0 * 1024.0),
            );
            limits.check_estimate(
                "final fragment-index",
                estimate
                    .fragment_peak_bytes
                    .saturating_sub(estimate.modified_peak_bytes),
            )?;
        }
        let database = database_parameters
            .clone()
            .build_from_peptides(all_peptides);

        cancellation.check()?;
        events.emit(EventKind::DatabaseBuilt {
            peptides: database.peptides.len(),
            fragments: database.fragments.len(),
        });
        events.check()?;

        info!(
            "generated {} fragments, {} peptides in {:#?}",
            database.fragments.len(),
            database.peptides.len(),
            (start.elapsed())
        );
        Ok(Self {
            database,
            parameters,
            database_parameters,
            library_search: None,
            start,
            events,
            cancellation,
        })
    }

    pub fn prefilter_peptides(
        self,
        parallel: usize,
        fasta: Fasta,
        custom_cleavages: Option<ValidatedCustomCleavageLibrary>,
    ) -> Vec<Peptide> {
        let spectra: Option<Vec<ProcessedSpectrum>> =
            match parallel >= self.parameters.mzml_paths.len() {
                true => Some(
                    self.read_processed_spectra(&self.parameters.mzml_paths, 0, 0)
                        .1,
                ),
                false => None,
            };

        let db_params = self.database_parameters.clone();
        // TODO: Don't generate decoys for fast searching
        // * if `generate_decoys` is used, we should re-generate at the end
        //  to ensure that picked-peptide conditions are used, otherwise,
        //  if the user supplied decoys in the fasta file, then we should retain them
        //
        // db_params.generate_decoys = false;

        let mut all_peptides: Vec<Peptide> = fasta
            .iter_chunks(self.database_parameters.prefilter_chunk_size)
            .enumerate()
            .flat_map(|(chunk_id, fasta_chunk)| {
                let start = Instant::now();
                info!("pre-filtering fasta chunk {}", chunk_id,);
                let mut db = db_params
                    .clone()
                    .build_with_custom_cleavages(fasta_chunk, custom_cleavages.as_ref());

                info!(
                    "generated {} fragments, {} peptides in {}ms",
                    db.fragments.len(),
                    db.peptides.len(),
                    (Instant::now() - start).as_millis()
                );

                let scorer = Scorer {
                    db: &db,
                    precursor_tol: self.parameters.precursor_tol,
                    fragment_tol: self.parameters.fragment_tol,
                    min_matched_peaks: self.parameters.min_matched_peaks,
                    min_isotope_err: self.parameters.isotope_errors.0,
                    max_isotope_err: self.parameters.isotope_errors.1,
                    min_precursor_charge: self.parameters.precursor_charge.0,
                    max_precursor_charge: self.parameters.precursor_charge.1,
                    override_precursor_charge: self.parameters.override_precursor_charge,
                    max_fragment_charge: self.parameters.max_fragment_charge,
                    chimera: self.parameters.chimera,
                    report_psms: self.parameters.report_psms + 1, // Q: Why is 1 being added here? (JSPP: Feb 2024)
                    wide_window: self.parameters.wide_window,
                    annotate_matches: false,
                    mass_shift_ppm: self.parameters.mass_shift_ppm,
                    score_type: self.parameters.score_type,
                    use_bitmap: self.parameters.use_bitmap,
                };

                // Allocate an array of booleans indicating whether a peptide was identified in a
                // preliminary pass of the data
                let keep = (0..db.peptides.len())
                    .map(|_| std::sync::atomic::AtomicBool::new(false))
                    .collect::<Vec<_>>();

                match &spectra {
                    Some(spectra) => self.peptide_filter_processed_spectra(&scorer, spectra, &keep),
                    None => self
                        .parameters
                        .mzml_paths
                        .chunks(parallel)
                        .enumerate()
                        .for_each(|(chunk_idx, chunk)| {
                            let spectra_chunk =
                                self.read_processed_spectra(chunk, chunk_idx, parallel).1;
                            self.peptide_filter_processed_spectra(&scorer, &spectra_chunk, &keep)
                        }),
                };

                // Retain only peptides where `keep[ix] = true`
                let peptides = db
                    .peptides
                    .drain(..)
                    .enumerate()
                    .filter_map(|(ix, peptide)| {
                        let val = keep[ix].load(std::sync::atomic::Ordering::Relaxed);
                        if val {
                            Some(peptide)
                        } else {
                            None
                        }
                    })
                    .collect::<Vec<_>>();

                info!(
                    "found {} pre-filtered peptides for fasta chunk {}",
                    peptides.len(),
                    chunk_id,
                );
                peptides
            })
            .collect();

        Parameters::reorder_peptides(&mut all_peptides);
        all_peptides
    }

    fn peptide_filter_processed_spectra(
        &self,
        scorer: &Scorer,
        spectra: &[ProcessedSpectrum],
        keep: &[std::sync::atomic::AtomicBool],
    ) {
        use std::sync::atomic::{AtomicUsize, Ordering};
        let counter = AtomicUsize::new(0);
        let start = Instant::now();

        spectra
            .par_iter()
            .filter(|spec| spec.masses.len() >= self.parameters.min_peaks && spec.level == 2)
            .for_each(|spectrum| {
                let prev = counter.fetch_add(1, Ordering::Relaxed);
                if prev > 0 && prev % 10_000 == 0 {
                    let duration = Instant::now().duration_since(start).as_millis() as usize;

                    let rate = prev * 1000 / (duration + 1);
                    log::trace!("- searched {} spectra ({} spectra/s)", prev, rate);
                }
                scorer.quick_score(
                    spectrum,
                    self.database_parameters.prefilter_low_memory,
                    keep,
                )
            });

        let duration = Instant::now().duration_since(start).as_millis() as usize;
        let prev = counter.load(Ordering::Relaxed);
        let rate = prev * 1000 / (duration + 1);
        log::info!(
            "- prefilter search:  {:8} ms ({} spectra/s)",
            duration,
            rate
        );
    }

    fn spectrum_fdr(&self, features: &mut Vec<Feature>) -> usize {
        if sage_core::ml::linear_discriminant::score_psms(features, self.parameters.precursor_tol)
            .is_none()
        {
            log::warn!("linear model fitting failed, falling back to heuristic discriminant score");
            self.events.emit(EventKind::Warning {
                code: "discriminant_model_fallback".into(),
                message: "linear model fitting failed; using heuristic discriminant score".into(),
            });
            features.par_iter_mut().for_each(|feat| {
                feat.discriminant_score = (-feat.poisson as f32).ln_1p() + feat.longest_y_pct / 3.0
            });
        }
        sort_features_by_discriminant(features);
        sage_core::ml::qvalue::spectrum_q_value(features)
    }

    /// Align systematic precursor and fragment mass errors per raw file before
    /// fitting the final FDR model. Models are trained only on provisional 1%
    /// spectrum-q rank-1 targets and are then applied equally to targets and
    /// decoys. Raw output errors remain unchanged.
    fn align_mass_errors(&self, features: &mut [Feature]) {
        features.par_iter_mut().for_each(|feature| {
            feature.aligned_delta_mass = feature.delta_mass;
            feature.aligned_average_ppm = feature.average_ppm;
        });

        let fit_options = FitOptions {
            // A line is useful even for modest drift; reject it only when its
            // held-in robust residual is worse than the static center.
            min_linear_improvement: 0.0,
            ..FitOptions::default()
        };

        for file_id in 0..self.parameters.mzml_paths.len() {
            let calibration_psms = features
                .iter()
                .filter(|feature| {
                    feature.file_id == file_id
                        && feature.rank == 1
                        && feature.label == 1
                        && feature.spectrum_q <= 0.01
                })
                .collect::<Vec<_>>();

            let precursor_points = calibration_psms
                .iter()
                .map(|feature| CalibrationPoint {
                    rt_minutes: feature.rt,
                    error_ppm: feature.delta_mass,
                })
                .collect::<Vec<_>>();
            let fragment_points = calibration_psms
                .iter()
                .map(|feature| CalibrationPoint {
                    rt_minutes: feature.rt,
                    error_ppm: feature.signed_fragment_ppm,
                })
                .collect::<Vec<_>>();

            let precursor_fit = matches!(self.parameters.precursor_tol, Tolerance::Ppm(_, _))
                .then(|| fit_mass_calibration(&precursor_points, fit_options))
                .flatten();
            let fragment_fit = fit_mass_calibration(&fragment_points, fit_options);

            if let Some(fit) = precursor_fit {
                log::info!(
                    "- file {} precursor mass alignment: {:?}, offset={:.3} ppm, slope={:.4} ppm/min, n={}",
                    file_id,
                    fit.model.kind,
                    fit.model.intercept_ppm,
                    fit.model.slope_ppm_per_min,
                    fit.inliers,
                );
            }
            if let Some(fit) = fragment_fit {
                log::info!(
                    "- file {} fragment mass alignment: {:?}, offset={:.3} ppm, slope={:.4} ppm/min, n={}",
                    file_id,
                    fit.model.kind,
                    fit.model.intercept_ppm,
                    fit.model.slope_ppm_per_min,
                    fit.inliers,
                );
            }

            features
                .iter_mut()
                .filter(|feature| feature.file_id == file_id)
                .for_each(|feature| {
                    if let Some(fit) = precursor_fit {
                        feature.aligned_delta_mass =
                            feature.delta_mass - fit.model.predict_ppm(feature.rt);
                    }
                    if let Some(fit) = fragment_fit {
                        let predicted = fit.model.predict_ppm(feature.rt);
                        // Preserve the within-PSM absolute-error spread while
                        // translating its signed center to the fitted baseline.
                        feature.aligned_average_ppm = align_fragment_error(
                            feature.average_ppm,
                            feature.signed_fragment_ppm,
                            predicted,
                        );
                    }
                });
        }
    }

    // Create a path for `file_name` in the specified output directory, if it exists,
    // otherwise, write to current directory
    fn make_path<S: AsRef<str>>(&self, file_name: S) -> Url {
        self.parameters
            .output_directory
            .join(file_name.as_ref())
            .expect("valid path segment")
    }

    fn search_processed_spectra(
        &self,
        scorer: &Scorer,
        msn_spectra: &[ProcessedSpectrum],
    ) -> Vec<Feature> {
        use std::sync::atomic::{AtomicUsize, Ordering};
        let counter = AtomicUsize::new(0);
        let start = Instant::now();

        let features: Vec<_> = msn_spectra
            .par_iter()
            .filter(|spec| {
                !self.cancellation.is_cancelled()
                    && spec.masses.len() >= self.parameters.min_peaks
                    && spec.level == 2
            })
            .map(|x| {
                let prev = counter.fetch_add(1, Ordering::Relaxed);
                if prev > 0 && prev % 10_000 == 0 {
                    let duration = Instant::now().duration_since(start).as_millis() as usize;

                    let rate = prev * 1000 / (duration + 1);
                    log::trace!("- searched {} spectra ({} spectra/s)", prev, rate);
                }
                x
            })
            .flat_map(|spec| scorer.score(spec))
            .collect();

        let duration = Instant::now().duration_since(start).as_millis() as usize;
        let prev = counter.load(Ordering::Relaxed);
        let rate = prev * 1000 / (duration + 1);
        log::info!("- search:  {:8} ms ({} spectra/s)", duration, rate);
        features
    }

    fn complete_features(
        &self,
        msn_spectra: Vec<ProcessedSpectrum>,
        ms1_spectra: Vec<ProcessedSpectrum>,
        features: Vec<Feature>,
    ) -> SageResults {
        let quant = self
            .parameters
            .quant
            .tmt
            .as_ref()
            .map(|isobaric| {
                let level = self.parameters.quant.tmt_settings.level;
                if level != 2 && level != 3 {
                    log::warn!("TMT quant level set at {}, is this correct?", level);
                }
                sage_core::tmt::quantify(&msn_spectra, isobaric, Tolerance::Ppm(-20.0, 20.0), level)
            })
            .unwrap_or_default();

        SageResults {
            features,
            quant,
            ms1: ms1_spectra,
        }
    }

    fn requires_ms1(&self) -> bool {
        self.parameters.quant.lfq
    }

    fn process_chunk(
        &self,
        scorer: &Scorer,
        chunk: &[Url],
        chunk_idx: usize,
        batch_size: usize,
    ) -> SageResults {
        let spectra = self.read_processed_spectra(chunk, chunk_idx, batch_size);
        let features = self.search_processed_spectra(scorer, &spectra.1);
        self.complete_features(spectra.1, spectra.0, features)
    }

    fn read_processed_spectra(
        &self,
        chunk: &[Url],
        chunk_idx: usize,
        batch_size: usize,
    ) -> (Vec<ProcessedSpectrum>, Vec<ProcessedSpectrum>) {
        self.read_processed_spectra_with_ms1(chunk, chunk_idx, batch_size, self.requires_ms1())
    }

    fn read_processed_spectra_with_ms1(
        &self,
        chunk: &[Url],
        chunk_idx: usize,
        batch_size: usize,
        requires_ms1: bool,
    ) -> (Vec<ProcessedSpectrum>, Vec<ProcessedSpectrum>) {
        // Read all of the spectra at once - this can help prevent memory over-consumption issues
        info!(
            "processing files {} .. {} ",
            batch_size * chunk_idx,
            batch_size * chunk_idx + chunk.len()
        );
        let start = Instant::now();

        let sn = self
            .parameters
            .quant
            .tmt_settings
            .sn
            .then_some(self.parameters.quant.tmt_settings.level);

        let min_deisotope_mz = match &self.parameters.quant.tmt {
            Some(i) => match self.parameters.quant.tmt_settings.level {
                2 => i.reporter_masses().last().map(|x| x * (1.0 + 20E-6)),
                _ => None,
            },
            None => None,
        };

        let sp = SpectrumProcessor::new(
            self.parameters.max_peaks,
            self.parameters.deisotope,
            min_deisotope_mz.unwrap_or(0.0),
        );

        // If the file format supports parallel reading, then we can read
        // then it is faster to read each file in series. (since each spectra
        // will be processed internally in parallel).
        let file_serial_read = chunk
            .iter()
            .all(|path| FileFormat::from(path.as_ref()).within_file_parallel());
        log::trace!("file serial read: {}", file_serial_read);
        let inner_closure = |(idx, path): (usize, &Url)| {
            let file_id = chunk_idx * batch_size + idx;
            self.events.emit(EventKind::FileStarted {
                file_id,
                path: path.to_string(),
            });
            let res = sage_cloudpath::util::read_spectra(
                path,
                file_id,
                sn,
                self.parameters.bruker_config,
                requires_ms1,
            );

            match res {
                Ok(s) => {
                    log::trace!("- {}: read {} spectra", path, s.len());
                    let spectra = s
                        .into_par_iter()
                        .map(|spectrum| sp.process(spectrum))
                        .collect::<SpectrumAccumulator>();
                    self.events.emit(EventKind::FileCompleted {
                        file_id,
                        path: path.to_string(),
                        spectra: spectra.ms1.len() + spectra.msn.len(),
                    });
                    Ok(spectra)
                }
                Err(e) => {
                    log::error!("- {}: {}", path, e);
                    Err(e)
                }
            }
        };

        let spectra: SpectrumAccumulator = if file_serial_read {
            chunk
                .iter()
                .enumerate()
                .flat_map(inner_closure)
                .fold(SpectrumAccumulator::default(), SpectrumAccumulator::reduce)
        } else {
            chunk
                .par_iter()
                .enumerate()
                .flat_map(inner_closure)
                .reduce(SpectrumAccumulator::default, SpectrumAccumulator::reduce)
        };

        let has_ims = spectra.ms1.iter().any(|x| !x.mobilities.is_empty());
        if spectra.ms1.is_empty() {
            log::trace!("no MS1 spectra found");
        } else {
            if has_ims {
                log::trace!("Processing MS1 spectra with IMS columns");
            } else {
                log::trace!("Processing MS1 spectra without IMS");
            }
        }

        self.events.emit(EventKind::SpectraProcessed {
            ms1_spectra: spectra.ms1.len(),
            msn_spectra: spectra.msn.len(),
        });

        let io_time = Instant::now() - start;
        info!("- file IO: {:8} ms", io_time.as_millis());

        (spectra.ms1, spectra.msn)
    }

    /// Re-read only the MS2 files needed for post-FDR work. Detailed matched
    /// fragments and PTM localization share this pass so enabling both does
    /// not double the input I/O.
    fn postprocess_features(
        &self,
        scorer: &Scorer,
        features: &mut [Feature],
        batch_size: usize,
    ) -> anyhow::Result<PostprocessStats> {
        #[derive(Default)]
        struct SpectrumWork {
            /// Every reported PSM for this spectrum, required to replay
            /// rank-ordered peak removal in chimera mode.
            feature_indices: Vec<usize>,
            annotate: bool,
            localization_indices: Vec<usize>,
        }

        let annotate_matches = self.parameters.annotate_matches;
        let localize = self.parameters.ptm_localization.enabled;
        let library_selections = spectral_library::select_psms(
            features,
            &self.database,
            &self.parameters.spectral_library,
        );
        let library_annotation_indices = library_selections
            .iter()
            .flat_map(|selection| selection.feature_indices.iter().copied())
            .collect::<HashSet<_>>();
        let annotate_fragments = annotate_matches || !library_annotation_indices.is_empty();
        if !annotate_fragments && !localize {
            return Ok(PostprocessStats {
                library_selections,
                ..PostprocessStats::default()
            });
        }

        let output_psm_q_value = self.parameters.output_filter.psm_q_value;
        let annotation_indices = features
            .iter()
            .enumerate()
            .filter_map(|(idx, feature)| {
                ((annotate_matches && passes_output_filter(feature, output_psm_q_value))
                    || library_annotation_indices.contains(&idx))
                .then_some(idx)
            })
            .collect::<HashSet<_>>();
        let mut work: HashMap<usize, HashMap<String, SpectrumWork>> = HashMap::new();

        if annotate_fragments {
            for (idx, feature) in features.iter().enumerate() {
                work.entry(feature.file_id)
                    .or_default()
                    .entry(feature.spec_id.clone())
                    .or_default()
                    .feature_indices
                    .push(idx);
            }
            for idx in &annotation_indices {
                let feature = &features[*idx];
                if let Some(spectrum) = work
                    .get_mut(&feature.file_id)
                    .and_then(|file| file.get_mut(feature.spec_id.as_str()))
                {
                    spectrum.annotate = true;
                }
            }
        }

        if localize {
            for (idx, feature) in features.iter().enumerate() {
                if passes_localization_filter(feature, self.parameters.ptm_localization.psm_q_value)
                    && sage_core::ptm::has_localizable_modification(
                        &self.database[feature.peptide_idx],
                        &self.database.potential_mods,
                    )
                {
                    work.entry(feature.file_id)
                        .or_default()
                        .entry(feature.spec_id.clone())
                        .or_default()
                        .localization_indices
                        .push(idx);
                }
            }
        }

        for file in work.values_mut() {
            file.retain(|_, spectrum| {
                spectrum.annotate || !spectrum.localization_indices.is_empty()
            });
            for spectrum in file.values_mut() {
                spectrum
                    .feature_indices
                    .sort_unstable_by_key(|idx| features[*idx].rank);
            }
        }
        work.retain(|_, file| !file.is_empty());

        let expected_annotations = annotation_indices.len();
        if work.is_empty() {
            anyhow::ensure!(
                expected_annotations == 0,
                "internal error: selected PSMs were not scheduled for fragment annotation"
            );
            return Ok(PostprocessStats {
                library_selections,
                ..PostprocessStats::default()
            });
        }

        let start = Instant::now();
        let mut annotations = Vec::with_capacity(expected_annotations);
        let mut localizations = Vec::new();

        for (chunk_idx, chunk) in self.parameters.mzml_paths.chunks(batch_size).enumerate() {
            self.cancellation.check()?;
            let first_file_id = chunk_idx * batch_size;
            if !(first_file_id..first_file_id + chunk.len())
                .any(|file_id| work.contains_key(&file_id))
            {
                continue;
            }

            let spectra = self
                .read_processed_spectra_with_ms1(chunk, chunk_idx, batch_size, false)
                .1;
            let results = spectra
                .par_iter()
                .filter_map(|spectrum| {
                    let spectrum_work = work.get(&spectrum.file_id)?.get(spectrum.id.as_str())?;
                    let mut annotated = Vec::new();

                    if spectrum_work.annotate {
                        let ranked_features = spectrum_work
                            .feature_indices
                            .iter()
                            .map(|&idx| &features[idx])
                            .collect::<Vec<_>>();
                        let selected = spectrum_work
                            .feature_indices
                            .iter()
                            .map(|idx| annotation_indices.contains(idx))
                            .collect::<Vec<_>>();
                        annotated.extend(
                            spectrum_work
                                .feature_indices
                                .iter()
                                .copied()
                                .zip(scorer.annotate_ranked_candidates(
                                    spectrum,
                                    &ranked_features,
                                    &selected,
                                ))
                                .filter_map(|(idx, fragments)| {
                                    fragments.map(|fragments| (idx, fragments))
                                }),
                        );
                    }

                    let localized = spectrum_work
                        .localization_indices
                        .iter()
                        .copied()
                        .filter_map(|idx| {
                            let feature = &features[idx];
                            let peptide = &self.database[feature.peptide_idx];
                            let localization = sage_core::ptm::localize(
                                peptide,
                                spectrum,
                                &self.database.ion_kinds,
                                &self.database.potential_mods,
                                self.parameters.fragment_tol,
                                self.parameters.max_fragment_charge,
                                feature.charge,
                            );
                            (!localization.mods.is_empty()).then_some((idx, localization))
                        })
                        .collect::<Vec<_>>();

                    Some((annotated, localized))
                })
                .collect::<Vec<_>>();

            for (mut batch_annotations, mut batch_localizations) in results {
                annotations.append(&mut batch_annotations);
                localizations.append(&mut batch_localizations);
            }
        }

        let mut annotated_indices = HashSet::with_capacity(annotations.len());
        let mut annotated_fragments = 0usize;
        for (idx, fragments) in annotations {
            anyhow::ensure!(
                annotated_indices.insert(idx),
                "spectrum {} in file {} was encountered more than once during deferred annotation",
                features[idx].spec_id,
                features[idx].file_id
            );
            annotated_fragments += fragments.fragment_ordinals.len();
            features[idx].fragments = Some(fragments);
        }

        if annotated_indices.len() != expected_annotations {
            let missing = features
                .iter()
                .enumerate()
                .find(|(idx, _feature)| {
                    annotation_indices.contains(idx) && !annotated_indices.contains(idx)
                })
                .map(|(_, feature)| {
                    format!(
                        "file {} spectrum {} (PSM {})",
                        feature.file_id, feature.spec_id, feature.psm_id
                    )
                })
                .unwrap_or_else(|| "unknown PSM".into());
            anyhow::bail!(
                "deferred fragment annotation completed {}/{} selected PSMs; missing {}",
                annotated_indices.len(),
                expected_annotations,
                missing
            );
        }

        let localized_psms = localizations.len();
        for (idx, localization) in localizations {
            features[idx].localization = Some(localization);
        }

        let localization_indices = features
            .iter()
            .enumerate()
            .flat_map(|(feature_idx, feature)| {
                feature.localization.iter().flat_map(move |localization| {
                    localization
                        .mods
                        .iter()
                        .enumerate()
                        .filter(|(_, modification)| modification.competition_eligible)
                        .map(move |(mod_idx, _)| (feature_idx, mod_idx))
                })
            })
            .collect::<Vec<_>>();
        let evidence = localization_indices
            .iter()
            .map(|&(feature_idx, mod_idx)| {
                let modification =
                    &features[feature_idx].localization.as_ref().unwrap().mods[mod_idx];
                (modification.target_decoy_score, modification.decoy_winner)
            })
            .collect::<Vec<_>>();
        let q_values = sage_core::ptm::target_decoy_q_values(&evidence);
        for ((feature_idx, mod_idx), q_value) in localization_indices.into_iter().zip(q_values) {
            features[feature_idx].localization.as_mut().unwrap().mods[mod_idx]
                .localization_q_value = q_value;
        }

        log::info!(
            "- post-FDR MS2 pass: {} annotated PSMs ({} fragments), {} localized PSMs in {} ms",
            annotated_indices.len(),
            annotated_fragments,
            localized_psms,
            start.elapsed().as_millis()
        );

        Ok(PostprocessStats {
            annotated_psms: annotated_indices.len(),
            annotated_fragments,
            localized_psms,
            library_selections,
        })
    }

    pub fn batch_files(&self, scorer: &Scorer, batch_size: usize) -> SageResults {
        self.parameters
            .mzml_paths
            .chunks(batch_size)
            .enumerate()
            .map(|(chunk_idx, chunk)| {
                let results = self.process_chunk(scorer, chunk, chunk_idx, batch_size);
                self.events.emit(EventKind::SearchProgress {
                    files_completed: (chunk_idx * batch_size + chunk.len())
                        .min(self.parameters.mzml_paths.len()),
                    files_total: self.parameters.mzml_paths.len(),
                });
                results
            })
            .collect::<SageResults>()
    }

    pub fn run(self, parallel: usize) -> anyhow::Result<telemetry::Telemetry> {
        self.run_with_summary(parallel)
            .map(|(telemetry, _summary)| telemetry)
    }

    pub fn run_with_summary(
        mut self,
        parallel: usize,
    ) -> anyhow::Result<(telemetry::Telemetry, RunSummary)> {
        anyhow::ensure!(parallel > 0, "batch size must be greater than zero");
        if self.library_search.is_some() {
            return self.run_library_with_summary(parallel);
        }
        self.cancellation.check()?;
        self.events.check()?;
        let scorer = Scorer {
            db: &self.database,
            precursor_tol: self.parameters.precursor_tol,
            fragment_tol: self.parameters.fragment_tol,
            min_matched_peaks: self.parameters.min_matched_peaks,
            min_isotope_err: self.parameters.isotope_errors.0,
            max_isotope_err: self.parameters.isotope_errors.1,
            min_precursor_charge: self.parameters.precursor_charge.0,
            max_precursor_charge: self.parameters.precursor_charge.1,
            override_precursor_charge: self.parameters.override_precursor_charge,
            max_fragment_charge: self.parameters.max_fragment_charge,
            chimera: self.parameters.chimera,
            report_psms: self.parameters.report_psms,
            wide_window: self.parameters.wide_window,
            annotate_matches: false,
            mass_shift_ppm: self.parameters.mass_shift_ppm,
            score_type: self.parameters.score_type,
            use_bitmap: self.parameters.use_bitmap,
        };

        //Collect all results into a single container
        let mut outputs = self.batch_files(&scorer, parallel);
        self.cancellation.check()?;
        self.events.check()?;

        // Establish provisional q-values from the search-only Poisson feature,
        // then use confident PSMs for mass-error alignment and property-model
        // training before final FDR fitting.
        outputs.features.par_sort_unstable_by(|left, right| {
            left.poisson
                .total_cmp(&right.poisson)
                .then_with(|| feature_identity_cmp(left, right))
        });
        assign_psm_ids(&mut outputs.features);
        sage_core::ml::qvalue::spectrum_q_value_by(&mut outputs.features, |feature| {
            feature.poisson
        });
        self.align_mass_errors(&mut outputs.features);
        self.events.emit(EventKind::MassAlignmentCompleted {
            files: self.parameters.mzml_paths.len(),
        });

        let has_ion_mobility = outputs
            .features
            .iter()
            .any(|feature| feature.ims.is_finite() && feature.ims > 0.0);

        let needs_alignment = self.parameters.predict_rt
            || self.parameters.quant.lfq
            || self.parameters.retention_time_alignment.is_some();
        let alignments = if needs_alignment {
            // Alignment landmarks are observed target PSMs passing 1% spectrum FDR.
            let alignment_method = self.parameters.retention_time_alignment.unwrap_or_default();
            let alignments = sage_core::ml::retention_alignment::global_alignment_with_method(
                &mut outputs.features,
                self.parameters.mzml_paths.len(),
                alignment_method,
            );
            self.events.emit(EventKind::RetentionAlignmentCompleted {
                method: format!("{alignment_method:?}").to_lowercase(),
                files: self.parameters.mzml_paths.len(),
            });
            Some(alignments)
        } else {
            None
        };

        let retention_time_model_fitted = if self.parameters.predict_rt {
            if sage_core::ml::retention_model::predict(
                &self.database,
                &mut outputs.features,
                &self.parameters.retention_time_model,
            )
            .is_some()
            {
                self.events.emit(EventKind::RtModelFitted);
                true
            } else {
                self.events.emit(EventKind::RtModelSkipped {
                    reason: "insufficient high-confidence observations or poor model fit".into(),
                });
                false
            }
        } else {
            self.events.emit(EventKind::RtModelSkipped {
                reason: "retention-time prediction disabled".into(),
            });
            false
        };
        let ion_mobility_model_fitted = if has_ion_mobility {
            if sage_core::ml::mobility_model::predict(
                &self.database,
                &mut outputs.features,
                &self.parameters.ion_mobility_model,
            )
            .is_some()
            {
                self.events.emit(EventKind::MobilityModelFitted);
                true
            } else {
                self.events.emit(EventKind::MobilityModelSkipped {
                    reason: "ion mobility unavailable or model fitting failed".into(),
                });
                false
            }
        } else {
            self.events.emit(EventKind::MobilityModelSkipped {
                reason: "no ion-mobility observations present".into(),
            });
            false
        };

        let q_spectrum = self.spectrum_fdr(&mut outputs.features);
        self.events.emit(EventKind::LdaScoringCompleted);
        let q_peptide = sage_core::fdr::picked_peptide(&self.database, &mut outputs.features);
        // Protein FDR is based exclusively on proteotypic (unique, non-shared) peptides. Shared peptides
        // are reported with protein FDR = 1.0
        let q_protein = sage_core::fdr::picked_protein(&self.database, &mut outputs.features);
        // Conducts "IDPicker-based protein grouping at 1% peptide FDR"
        sage_core::protein_grouping::generate_protein_groups(
            &self.database,
            &mut outputs.features,
            self.parameters.protein_grouping,
            Some(self.parameters.protein_grouping_peptide_fdr),
        );
        // Uses the "Picked Group FDR" approach to compute protein group FDR for the IDPicker groups,
        // including rescued subset grouping (rsG). Shared peptides (between different groups)
        // are reported with protein group FDR = 1.0
        let q_protein_group =
            sage_core::fdr::picked_protein_group(&self.database, &mut outputs.features);
        self.events.emit(EventKind::FdrCompleted {
            psms: q_spectrum,
            peptides: q_peptide,
            proteins: q_protein,
            protein_groups: q_protein_group,
        });
        self.cancellation.check()?;

        let postprocess = self.postprocess_features(&scorer, &mut outputs.features, parallel)?;
        if self.parameters.annotate_matches {
            self.events.emit(EventKind::FragmentAnnotationCompleted {
                psms: postprocess.annotated_psms,
                fragments: postprocess.annotated_fragments,
            });
        }
        let localized_psms = if self.parameters.ptm_localization.enabled {
            self.events.emit(EventKind::PtmLocalizationCompleted {
                psms: postprocess.localized_psms,
            });
            postprocess.localized_psms
        } else {
            0
        };

        let filenames = self
            .parameters
            .mzml_paths
            .iter()
            .map(|url| {
                sage_cloudpath::filename(url)
                    .unwrap_or_else(|| url.as_str())
                    .to_string()
            })
            .collect::<Vec<_>>();

        let library_entries = spectral_library::build_entries(
            &outputs.features,
            &self.database,
            &filenames,
            &postprocess.library_selections,
            &self.parameters.spectral_library,
        )
        .map_err(anyhow::Error::msg)?;
        let library_transitions = library_entries
            .iter()
            .map(|entry| entry.fragments.len())
            .sum::<usize>();

        let areas = alignments.and_then(|alignments| {
            if self.parameters.quant.lfq {
                log::trace!("performing LFQ");
                let mut areas = sage_core::lfq::build_feature_map(
                    self.parameters.quant.lfq_settings,
                    self.parameters.precursor_charge,
                    &outputs.features,
                    &self.database,
                )
                .quantify(&self.database, &outputs.ms1, &alignments);

                let q_precursor = sage_core::fdr::picked_precursor(&mut areas);
                self.events.emit(EventKind::QuantificationCompleted {
                    kind: "lfq".into(),
                    features: areas.len(),
                });

                log::info!("discovered {} target MS1 peaks at 5% FDR", q_precursor);
                Some(areas)
            } else {
                None
            }
        });

        if !outputs.quant.is_empty() {
            self.events.emit(EventKind::QuantificationCompleted {
                kind: "tmt".into(),
                features: outputs.quant.len(),
            });
        }
        let lfq_features = areas.as_ref().map(|areas| areas.len()).unwrap_or_default();
        let tmt_features = outputs.quant.len();

        log::info!(
            "discovered {} target peptide-spectrum matches at 1% FDR",
            q_spectrum
        );
        log::info!("discovered {} target peptides at 1% FDR", q_peptide);
        log::info!(
            "discovered {} target proteins (supported by proteotypic peptides only) at 1% FDR",
            q_protein
        );
        log::info!(
            "discovered {} target protein groups (supported by proteotypic peptides only) at 1% FDR",
            q_protein_group
        );
        log::trace!("writing outputs");

        let output_psm_q_value = self.parameters.output_filter.psm_q_value;
        let output_features = outputs
            .features
            .iter()
            .filter(|feature| passes_output_filter(feature, output_psm_q_value))
            .collect::<Vec<_>>();
        log::info!(
            "writing {} of {} PSMs at spectrum q-value <= {}",
            output_features.len(),
            outputs.features.len(),
            output_psm_q_value
        );

        let bytes = sage_cloudpath::parquet::serialize_features(
            &output_features,
            &outputs.quant,
            &filenames,
            &self.database,
            output_psm_q_value,
        )?;

        let path = self.make_path("results.sage.parquet");
        sage_cloudpath::write_bytes_sync(&path, bytes)?;
        self.parameters.output_paths.push(path);

        if self.parameters.annotate_matches {
            let bytes = sage_cloudpath::parquet::serialize_matched_fragments(
                &output_features,
                output_psm_q_value,
            )?;
            let path = self.make_path("matched_fragments.sage.parquet");
            sage_cloudpath::write_bytes_sync(&path, bytes)?;
            self.parameters.output_paths.push(path);
        }

        if self
            .parameters
            .spectral_library
            .writes(SpectralLibraryFormat::SageParquet)
        {
            let bytes = sage_cloudpath::parquet::serialize_spectral_library(
                &library_entries,
                &self.parameters.spectral_library,
            )?;
            let path = self.make_path("spectral_library.sage.parquet");
            sage_cloudpath::write_bytes_sync(&path, bytes)?;
            self.parameters.output_paths.push(path);
        }
        if self
            .parameters
            .spectral_library
            .writes(SpectralLibraryFormat::MzSpecLib)
        {
            let bytes = spectral_library::serialize_mzspeclib(
                &library_entries,
                self.parameters.version.as_str(),
                self.parameters.spectral_library.strategy,
            );
            let path = self.make_path("spectral_library.mzspeclib.txt");
            sage_cloudpath::write_bytes_sync(&path, bytes)?;
            self.parameters.output_paths.push(path);
        }
        if self.parameters.spectral_library.enabled {
            self.events.emit(EventKind::SpectralLibraryCompleted {
                entries: library_entries.len(),
                transitions: library_transitions,
                formats: self.parameters.spectral_library.formats.len(),
            });
        }

        if let Some(areas) = &areas {
            let bytes = sage_cloudpath::parquet::serialize_lfq(areas, &filenames, &self.database)?;

            let path = self.make_path("lfq.parquet");
            sage_cloudpath::write_bytes_sync(&path, bytes)?;
            self.parameters.output_paths.push(path);
        }

        // PTM site reports follow the selected main output format.
        if self.parameters.ptm_localization.enabled {
            self.parameters
                .output_paths
                .push(self.write_ptm_sites(&outputs.features, &filenames)?);
            self.parameters
                .output_paths
                .push(self.write_protein_sites(&outputs.features, &filenames)?);
            self.parameters
                .output_paths
                .extend(self.write_ptm_library(&outputs.features, &filenames)?);
        }

        // Write percolator input file if requested
        if self.parameters.write_pin {
            self.parameters
                .output_paths
                .push(self.write_pin(&outputs.features, &filenames)?);
        }

        // Write an html report if requested
        if self.parameters.write_report {
            self.parameters.output_paths.push(self.write_report(
                &outputs.features,
                areas,
                &filenames,
            )?);
        }

        let path = self.make_path("results.json");
        if !self.events.is_enabled() {
            println!("{}", serde_json::to_string_pretty(&self.parameters)?);
        }

        let bytes = serde_json::to_vec_pretty(&self.parameters)?;
        sage_cloudpath::write_bytes_sync(&path, bytes)?;
        self.parameters.output_paths.push(path);

        let run_time = (Instant::now() - self.start).as_secs();
        info!("finished in {}s", run_time);
        info!("cite: \"Sage: An Open-Source Tool for Fast Proteomics Searching and Quantification at Scale\" https://doi.org/10.1021/acs.jproteome.3c00486");

        let summary_path = self.make_path("run-summary.json");
        let mut output_paths = self
            .parameters
            .output_paths
            .iter()
            .map(ToString::to_string)
            .collect::<Vec<_>>();
        output_paths.push(summary_path.to_string());
        let mut input_stats = InputRunStats::default();
        for path in &self.parameters.mzml_paths {
            let path = path.path().trim_end_matches('/').to_ascii_lowercase();
            if path.ends_with(".raw") {
                input_stats.thermo_raw_files += 1;
            } else if path.ends_with(".mzml") || path.ends_with(".mzml.gz") {
                input_stats.mzml_files += 1;
            } else {
                input_stats.other_files += 1;
            }
        }
        let summary = RunSummary {
            schema_version: 7,
            runtime_secs: run_time,
            files: self.parameters.mzml_paths.len(),
            peptides_in_database: self.database.peptides.len(),
            fragments_in_database: self.database.fragments.len(),
            psms_at_one_percent_fdr: q_spectrum,
            peptides_at_one_percent_fdr: q_peptide,
            proteins_at_one_percent_fdr: q_protein,
            protein_groups_at_one_percent_fdr: q_protein_group,
            ptm_localization: PtmLocalizationRunStats {
                enabled: self.parameters.ptm_localization.enabled,
                localized_psms,
                psm_q_value: self.parameters.ptm_localization.psm_q_value,
                localization_q_value: self.parameters.ptm_localization.localization_q_value,
            },
            models: ModelRunStats {
                mass_alignment_applied: true,
                retention_time_prediction_enabled: self.parameters.predict_rt,
                retention_time_model_fitted,
                retention_time_features: format!(
                    "{:?}",
                    self.parameters.retention_time_model.features
                )
                .to_lowercase(),
                retention_time_alignment: needs_alignment.then(|| {
                    format!(
                        "{:?}",
                        self.parameters.retention_time_alignment.unwrap_or_default()
                    )
                    .to_lowercase()
                }),
                ion_mobility_observed: has_ion_mobility,
                ion_mobility_model_enabled: self.parameters.ion_mobility_model.enabled,
                ion_mobility_model_fitted,
                ion_mobility_features: format!("{:?}", self.parameters.ion_mobility_model.features)
                    .to_lowercase(),
                ..ModelRunStats::default()
            },
            quantification: QuantificationRunStats {
                lfq_enabled: self.parameters.quant.lfq,
                lfq_features,
                tmt: self
                    .parameters
                    .quant
                    .tmt
                    .as_ref()
                    .map(|tmt| format!("{tmt:?}").to_lowercase()),
                tmt_features,
                ms1_label_channels: self.database.label_channels.len(),
                ms1_label_reference: self.database.label_reference.as_deref().map(str::to_owned),
            },
            execution: ExecutionRunStats {
                batch_size: self.parameters.batch_size,
                parallelism: parallel,
                bitmap_search: self.parameters.use_bitmap,
                max_memory_gb: self.parameters.max_memory_gb,
                min_free_memory_gb: self.parameters.min_free_memory_gb,
            },
            inputs: input_stats,
            modifications: ModificationRunStats {
                static_definitions: self.database_parameters.static_mods.len(),
                variable_definitions: self
                    .database_parameters
                    .variable_mods
                    .values()
                    .map(Vec::len)
                    .sum(),
                max_variable_mods: self.database_parameters.max_variable_mods,
                max_total_variable_mods: self.database_parameters.max_total_variable_mods,
                max_combinations: self.database_parameters.max_combinations,
                ptm_library_sites: self
                    .database_parameters
                    .loaded_ptm_library
                    .as_deref()
                    .map_or(0, |library| library.len()),
                label_channels: self.database.label_channels.len(),
                labeled_peptides: self
                    .database
                    .peptides
                    .iter()
                    .filter(|peptide| !peptide.decoy && peptide.label_channel.is_some())
                    .count(),
            },
            spectral_library: SpectralLibraryRunStats {
                enabled: self.parameters.spectral_library.enabled,
                entries: library_entries.len(),
                transitions: library_transitions,
                strategy: match self.parameters.spectral_library.strategy {
                    SpectralLibraryStrategy::BestPsm => "best_psm".into(),
                    SpectralLibraryStrategy::Consensus => "consensus".into(),
                },
                psm_q_value: self.parameters.spectral_library.psm_q_value,
                peptide_q_value: self.parameters.spectral_library.peptide_q_value,
                formats: self
                    .parameters
                    .spectral_library
                    .formats
                    .iter()
                    .map(|format| match format {
                        SpectralLibraryFormat::SageParquet => "sage_parquet".into(),
                        SpectralLibraryFormat::MzSpecLib => "mzspeclib".into(),
                    })
                    .collect(),
            },
            library_search: LibrarySearchRunStats::default(),
            output_paths,
        };
        sage_cloudpath::write_bytes_sync(&summary_path, serde_json::to_vec_pretty(&summary)?)?;
        self.parameters.output_paths.push(summary_path);

        for path in &self.parameters.output_paths {
            self.events.emit(EventKind::OutputWritten {
                path: path.to_string(),
            });
        }

        let telemetry = telemetry::Telemetry::new(
            self.parameters,
            self.database.peptides.len(),
            self.database.fragments.len(),
            run_time,
        );

        self.events.emit(EventKind::JobCompleted {
            runtime_secs: run_time,
            outputs: summary.output_paths.len(),
        });
        self.events.check()?;

        Ok((telemetry, summary))
    }
    /// Flatten FDR-passing target PSMs into one [`SiteRow`] per localized
    /// modification site. Shared by the PSM-site and protein-site reports.
    fn collect_site_rows(&self, features: &[Feature], filenames: &[String]) -> Vec<SiteRow> {
        let mut rows = Vec::new();
        for feature in features {
            // Only confidently-identified target PSMs.
            if !passes_localization_filter(feature, self.parameters.ptm_localization.psm_q_value) {
                continue;
            }
            let localization = match &feature.localization {
                Some(loc) => loc,
                None => continue,
            };
            let peptide = &self.database[feature.peptide_idx];
            let peptide_str = peptide.to_string();
            let proteins =
                peptide.proteins(&self.database.decoy_tag, self.database.generate_decoys);
            let filename = filenames.get(feature.file_id).cloned().unwrap_or_default();

            for m in &localization.mods {
                if m.decoy_winner
                    || m.localization_q_value
                        > self.parameters.ptm_localization.localization_q_value
                {
                    continue;
                }
                let modification = m.label.clone().unwrap_or_else(|| format!("{:+}", m.mass));
                let site_probabilities = m
                    .all_sites
                    .iter()
                    .map(|s| {
                        format!(
                            "{}{}:{:.4}",
                            s.residue as char,
                            s.position + 1,
                            s.probability
                        )
                    })
                    .collect::<Vec<_>>()
                    .join(";");

                for site in &m.best_sites {
                    rows.push(SiteRow {
                        psm_id: feature.psm_id,
                        filename: filename.clone(),
                        scannr: feature.spec_id.clone(),
                        peptide: peptide_str.clone(),
                        peptide_sequence: String::from_utf8_lossy(&peptide.sequence).into_owned(),
                        proteins: proteins.clone(),
                        charge: feature.charge,
                        spectrum_q: feature.spectrum_q,
                        peptide_q: feature.peptide_q,
                        modification: modification.clone(),
                        modification_mass: m.mass,
                        position: site.position + 1,
                        residue: site.residue,
                        localization_probability: site.probability,
                        delta_score: m.delta_score,
                        target_decoy_score: m.target_decoy_score,
                        localization_q_value: m.localization_q_value,
                        candidate_sites: m.candidate_sites,
                        site_determining_matched: m.site_determining_matched,
                        site_determining_total: m.site_determining_ions,
                        site_probabilities: site_probabilities.clone(),
                    });
                }
            }
        }
        rows
    }

    /// Write a per-PSM-site PTM localization report (one row per localized
    /// modification site of each FDR-passing PSM).
    pub fn write_ptm_sites(
        &self,
        features: &[Feature],
        filenames: &[String],
    ) -> anyhow::Result<Url> {
        let rows = self.collect_site_rows(features, filenames);

        use sage_cloudpath::parquet::PtmSiteRecord;
        let records = rows
            .iter()
            .map(|row| PtmSiteRecord {
                psm_id: row.psm_id as i64,
                filename: row.filename.clone(),
                scannr: row.scannr.clone(),
                peptide: row.peptide.clone(),
                proteins: row.proteins.clone(),
                charge: row.charge as i32,
                spectrum_q: row.spectrum_q,
                peptide_q: row.peptide_q,
                modification: row.modification.clone(),
                modification_mass: row.modification_mass,
                position: row.position as i32,
                residue: (row.residue as char).to_string(),
                localization_probability: row.localization_probability,
                delta_localization_score: row.delta_score,
                target_decoy_score: row.target_decoy_score,
                localization_q_value: row.localization_q_value,
                candidate_sites: row.candidate_sites as i32,
                site_determining_ions_matched: row.site_determining_matched as i32,
                site_determining_ions_total: row.site_determining_total as i32,
                site_probabilities: row.site_probabilities.clone(),
            })
            .collect::<Vec<_>>();
        let path = self.make_path("results.sage.ptm-sites.parquet");
        let bytes = sage_cloudpath::parquet::serialize_ptm_sites(&records)?;
        sage_cloudpath::write_bytes_sync(&path, bytes)?;
        Ok(path)
    }

    /// Write a collapsed protein-site report: the best localization for each
    /// (protein, modified peptide site) aggregated across all supporting PSMs.
    pub fn write_protein_sites(
        &self,
        features: &[Feature],
        filenames: &[String],
    ) -> anyhow::Result<Url> {
        let rows = self.collect_site_rows(features, filenames);

        // Key on (protein, peptide, position, mod mass). Protein coordinates
        // are not resolved (the FASTA is consumed during indexing), so a row
        // represents a localized site on a peptide, attributed to each protein
        // the peptide maps to.
        #[derive(Clone)]
        struct Agg {
            protein: String,
            peptide: String,
            residue: u8,
            position: usize,
            modification: String,
            modification_mass: f32,
            n_psms: u32,
            best_probability: f32,
            best_delta_score: f32,
            best_localization_q_value: f32,
            best_spectrum_q: f32,
        }

        let mut map: HashMap<(String, String, usize, i64), Agg> = HashMap::new();
        for row in &rows {
            for protein in row.proteins.split(';').filter(|p| !p.is_empty()) {
                let mass_key = (row.modification_mass * 1e3).round() as i64;
                let key = (
                    protein.to_string(),
                    row.peptide.clone(),
                    row.position,
                    mass_key,
                );
                let entry = map.entry(key).or_insert_with(|| Agg {
                    protein: protein.to_string(),
                    peptide: row.peptide.clone(),
                    residue: row.residue,
                    position: row.position,
                    modification: row.modification.clone(),
                    modification_mass: row.modification_mass,
                    n_psms: 0,
                    best_probability: 0.0,
                    best_delta_score: f32::MIN,
                    best_localization_q_value: 1.0,
                    best_spectrum_q: f32::MAX,
                });
                entry.n_psms += 1;
                entry.best_probability = entry.best_probability.max(row.localization_probability);
                entry.best_delta_score = entry.best_delta_score.max(row.delta_score);
                entry.best_localization_q_value = entry
                    .best_localization_q_value
                    .min(row.localization_q_value);
                entry.best_spectrum_q = entry.best_spectrum_q.min(row.spectrum_q);
            }
        }

        let mut aggregated: Vec<Agg> = map.into_values().collect();
        aggregated.sort_by(|a, b| {
            a.protein
                .cmp(&b.protein)
                .then_with(|| a.peptide.cmp(&b.peptide))
                .then_with(|| a.position.cmp(&b.position))
        });

        use sage_cloudpath::parquet::ProteinSiteRecord;
        let records = aggregated
            .iter()
            .map(|agg| ProteinSiteRecord {
                protein: agg.protein.clone(),
                peptide: agg.peptide.clone(),
                residue: (agg.residue as char).to_string(),
                position_in_peptide: agg.position as i32,
                modification: agg.modification.clone(),
                modification_mass: agg.modification_mass,
                num_psms: agg.n_psms as i32,
                best_localization_probability: agg.best_probability,
                best_delta_localization_score: agg.best_delta_score,
                best_localization_q_value: agg.best_localization_q_value,
                best_spectrum_q: agg.best_spectrum_q,
            })
            .collect::<Vec<_>>();
        let path = self.make_path("results.sage.protein-sites.parquet");
        let bytes = sage_cloudpath::parquet::serialize_protein_sites(&records)?;
        sage_cloudpath::write_bytes_sync(&path, bytes)?;
        Ok(path)
    }

    /// Emit a compact, reusable protein-coordinate site library from passing
    /// localized PSMs. Only names defined by this search's `variable_mods` are
    /// included, so every emitted row can be resolved by the same config.
    fn write_ptm_library(
        &self,
        features: &[Feature],
        filenames: &[String],
    ) -> anyhow::Result<Vec<Url>> {
        if self.database_parameters.fasta.is_empty() {
            return Ok(Vec::new());
        }

        let known_names = self
            .database_parameters
            .variable_mods
            .values()
            .flatten()
            .filter_map(|entry| entry.definition().name.map(|name| name.to_string()))
            .collect::<HashSet<_>>();
        let fasta_url = sage_cloudpath::to_url(&self.database_parameters.fasta)?;
        let fasta = sage_cloudpath::util::read_fasta(
            &fasta_url,
            &self.database_parameters.decoy_tag,
            self.database_parameters.generate_decoys,
        )?;
        let proteins = fasta
            .targets
            .iter()
            .map(|(accession, sequence)| (accession.as_ref(), sequence.as_str()))
            .collect::<HashMap<_, _>>();

        let mut sites = HashSet::new();
        let mut skipped_unnamed = 0usize;
        for row in self.collect_site_rows(features, filenames) {
            if !known_names.contains(&row.modification) {
                skipped_unnamed += 1;
                continue;
            }
            let peptide_position = row.position - 1;
            for protein in row
                .proteins
                .split(';')
                .filter(|protein| !protein.is_empty())
            {
                let Some(sequence) = proteins.get(protein) else {
                    continue;
                };
                for (start, _) in sequence.match_indices(&row.peptide_sequence) {
                    let protein_position = start + peptide_position;
                    if sequence.as_bytes().get(protein_position) == Some(&row.residue) {
                        sites.insert(sage_core::ptm_library::PtmLibrarySite {
                            protein: Arc::from(protein),
                            position: protein_position as u32,
                            residue: row.residue,
                            modification: Arc::from(row.modification.as_str()),
                        });
                    }
                }
            }
        }
        if skipped_unnamed > 0 {
            log::warn!(
                "skipped {} localized sites without a matching configured modification name",
                skipped_unnamed
            );
        }

        let mut sites = sites.into_iter().collect::<Vec<_>>();
        sites.sort_unstable_by(|a, b| {
            a.protein
                .cmp(&b.protein)
                .then_with(|| a.position.cmp(&b.position))
                .then_with(|| a.modification.cmp(&b.modification))
        });
        let parquet_path = self.make_path("results.sage.ptm-library.parquet");
        let bytes = sage_cloudpath::parquet::serialize_ptm_library(&sites)?;
        sage_cloudpath::write_bytes_sync(&parquet_path, bytes)?;

        let tsv_path = self.make_path("results.sage.ptm-library.tsv");
        let mut writer = csv::WriterBuilder::new()
            .delimiter(b'\t')
            .from_writer(Vec::new());
        writer.write_record(["protein", "position", "residue", "modification"])?;
        for site in &sites {
            writer.write_record([
                site.protein.as_ref(),
                &(site.position + 1).to_string(),
                std::str::from_utf8(&[site.residue])?,
                site.modification.as_ref(),
            ])?;
        }
        writer.flush()?;
        sage_cloudpath::write_bytes_sync(&tsv_path, writer.into_inner()?)?;
        Ok(vec![parquet_path, tsv_path])
    }

    fn serialize_pin(
        &self,
        re: &regex::Regex,
        feature: &Feature,
        filenames: &[String],
    ) -> csv::ByteRecord {
        let scannr = re
            .captures_iter(&feature.spec_id)
            .last()
            .and_then(|cap| cap.get(1).map(|cap| cap.as_str()))
            .unwrap_or(&feature.spec_id);

        let mut record = csv::ByteRecord::new();
        let peptide = &self.database[feature.peptide_idx];
        record.push_field(itoa::Buffer::new().format(feature.psm_id).as_bytes());
        record.push_field(itoa::Buffer::new().format(feature.label).as_bytes());
        record.push_field(scannr.as_bytes());
        record.push_field(ryu::Buffer::new().format(feature.expmass).as_bytes());
        record.push_field(ryu::Buffer::new().format(feature.calcmass).as_bytes());
        record.push_field(filenames[feature.file_id].as_bytes());
        record.push_field(ryu::Buffer::new().format(feature.rt).as_bytes());
        record.push_field(ryu::Buffer::new().format(feature.ims).as_bytes());
        record.push_field(itoa::Buffer::new().format(feature.rank).as_bytes());
        record.push_field(
            itoa::Buffer::new()
                .format((feature.charge == 2) as i32)
                .as_bytes(),
        );
        record.push_field(
            itoa::Buffer::new()
                .format((feature.charge == 3) as i32)
                .as_bytes(),
        );
        record.push_field(
            itoa::Buffer::new()
                .format((feature.charge == 4) as i32)
                .as_bytes(),
        );
        record.push_field(
            itoa::Buffer::new()
                .format((feature.charge == 5) as i32)
                .as_bytes(),
        );
        record.push_field(
            itoa::Buffer::new()
                .format((feature.charge == 6) as i32)
                .as_bytes(),
        );
        record.push_field(
            itoa::Buffer::new()
                .format(if feature.charge < 2 || feature.charge > 6 {
                    feature.charge
                } else {
                    0
                })
                .as_bytes(),
        );
        record.push_field(itoa::Buffer::new().format(feature.peptide_len).as_bytes());
        record.push_field(
            itoa::Buffer::new()
                .format(feature.missed_cleavages)
                .as_bytes(),
        );
        record.push_field(
            itoa::Buffer::new()
                .format(peptide.semi_enzymatic as u8)
                .as_bytes(),
        );
        record.push_field(ryu::Buffer::new().format(feature.isotope_error).as_bytes());
        record.push_field(
            ryu::Buffer::new()
                .format(feature.delta_mass.abs().ln_1p())
                .as_bytes(),
        );
        record.push_field(ryu::Buffer::new().format(feature.average_ppm).as_bytes());
        record.push_field(
            ryu::Buffer::new()
                .format(feature.hyperscore.ln_1p())
                .as_bytes(),
        );
        record.push_field(
            ryu::Buffer::new()
                .format(feature.delta_next.ln_1p())
                .as_bytes(),
        );
        record.push_field(
            ryu::Buffer::new()
                .format(feature.delta_best.ln_1p())
                .as_bytes(),
        );
        record.push_field(ryu::Buffer::new().format(feature.aligned_rt).as_bytes());
        record.push_field(ryu::Buffer::new().format(feature.predicted_rt).as_bytes());
        record.push_field(
            ryu::Buffer::new()
                .format(feature.delta_rt_model.clamp(0.001, 1.0).sqrt())
                .as_bytes(),
        );
        record.push_field(ryu::Buffer::new().format(feature.predicted_ims).as_bytes());
        record.push_field(
            ryu::Buffer::new()
                .format(feature.delta_ims_model)
                .as_bytes(),
        );
        record.push_field(itoa::Buffer::new().format(feature.matched_peaks).as_bytes());
        record.push_field(itoa::Buffer::new().format(feature.longest_b).as_bytes());
        record.push_field(itoa::Buffer::new().format(feature.longest_y).as_bytes());
        record.push_field(ryu::Buffer::new().format(feature.longest_y_pct).as_bytes());
        record.push_field(
            ryu::Buffer::new()
                .format(feature.matched_intensity_pct.ln_1p())
                .as_bytes(),
        );
        record.push_field(
            itoa::Buffer::new()
                .format(feature.scored_candidates)
                .as_bytes(),
        );
        record.push_field(
            ryu::Buffer::new()
                .format((-feature.poisson).ln_1p())
                .as_bytes(),
        );
        record.push_field(
            ryu::Buffer::new()
                .format(feature.posterior_error)
                .as_bytes(),
        );
        record.push_field(peptide.to_string().as_bytes());
        record.push_field(
            peptide
                .proteins(&self.database.decoy_tag, self.database.generate_decoys)
                .as_bytes(),
        );
        record
    }

    pub fn write_pin(&self, features: &[Feature], filenames: &[String]) -> anyhow::Result<Url> {
        let path = self.make_path("results.sage.pin");

        let mut wtr = csv::WriterBuilder::new()
            .delimiter(b'\t')
            .from_writer(OutputTarget::new(&path)?);

        let headers = csv::ByteRecord::from(vec![
            "SpecId",
            "Label",
            "ScanNr",
            "ExpMass",
            "CalcMass",
            "FileName",
            "retentiontime",
            "ion_mobility",
            "rank",
            "z=2",
            "z=3",
            "z=4",
            "z=5",
            "z=6",
            "z=other",
            "peptide_len",
            "missed_cleavages",
            "semi_enzymatic",
            "isotope_error",
            "ln(precursor_ppm)",
            "fragment_ppm",
            "ln(hyperscore)",
            "ln(delta_next)",
            "ln(delta_best)",
            "aligned_rt",
            "predicted_rt",
            "sqrt(delta_rt_model)",
            "predicted_mobility",
            "sqrt(delta_mobility)",
            "matched_peaks",
            "longest_b",
            "longest_y",
            "longest_y_pct",
            "ln(matched_intensity_pct)",
            "scored_candidates",
            "ln(-poisson)",
            "posterior_error",
            "Peptide",
            "Proteins",
        ]);

        let re = regex::Regex::new(r"scan=(\d+)").expect("This is valid regex");

        wtr.write_byte_record(&headers)?;
        for chunk in features.chunks(1024) {
            for record in chunk
                .par_iter()
                .map(|feat| self.serialize_pin(&re, feat, filenames))
                .collect::<Vec<_>>()
            {
                wtr.write_byte_record(&record)?;
            }
        }

        finish_csv_writer(wtr, &path)?;
        Ok(path)
    }

    fn write_report(
        &self,
        features: &[Feature],
        areas: Option<HashMap<(PrecursorId, bool), QuantifiedPeak, fnv::FnvBuildHasher>>,
        filenames: &[String],
    ) -> anyhow::Result<Url> {
        let path = self.make_path("results.sage.report.html");

        let global_q_value_filter = 0.01;
        let predict_section_q_value_filter = 0.01;

        // Create a new report
        let mut report = Report::new(
            "Sage",
            &self.parameters.version,
            Some(
                "https://github.com/pgarrett-scripps/sage-plus/blob/main/figures/logo.png?raw=true",
            ),
            "Sage Report",
        );

        /* Section 1: Introduction */
        {
            let mut intro_section = ReportSection::new("Results Overview");
            intro_section.add_content(html! {
                "The following files were processed:"
                ul {
                    @for filename in filenames {
                        li { (filename) }
                    }
                }
            });

            // Number of targets identified at global q-value filter at spectrum level per file
            let num_psm_targets_per_file: Vec<usize> = filenames
                .iter()
                .map(|filename| {
                    features
                        .iter()
                        .filter(|f| {
                            f.label == 1
                                && f.spectrum_q <= global_q_value_filter
                                && filenames[f.file_id] == *filename
                        })
                        .count()
                })
                .collect();

            // Number of peptides identified at global q-value filter at peptide level per file
            let mut num_peptide_targets_per_file: Vec<usize> = Vec::new();
            for filename in filenames {
                let mut peptides = HashSet::new();
                for feature in features.iter().filter(|f| {
                    f.label == 1
                        && f.peptide_q <= global_q_value_filter
                        && filenames[f.file_id] == *filename
                }) {
                    peptides.insert(self.database[feature.peptide_idx].to_string());
                }
                num_peptide_targets_per_file.push(peptides.len());
            }

            // Number of proteins identified at global q-value filter at protein level per file
            let mut num_protein_targets_per_file: Vec<usize> = Vec::new();
            for filename in filenames {
                let mut proteins = HashSet::new();
                for feature in features.iter().filter(|f| {
                    f.label == 1
                        && f.protein_q <= global_q_value_filter
                        && filenames[f.file_id] == *filename
                }) {
                    proteins.insert(
                        self.database[feature.peptide_idx]
                            .proteins(&self.database.decoy_tag, self.database.generate_decoys),
                    );
                }
                num_protein_targets_per_file.push(proteins.len());
            }

            // Total MS2 intensity at global q-value filter at each level per file
            let total_ms2_intensity_per_file: Vec<f32> = filenames
                .iter()
                .map(|filename| {
                    features
                        .iter()
                        .filter(|f| {
                            f.label == 1
                                && f.spectrum_q <= global_q_value_filter
                                && f.peptide_q <= global_q_value_filter
                                && f.protein_q <= global_q_value_filter
                                && filenames[f.file_id] == *filename
                        })
                        .map(|f| f.ms2_intensity)
                        .sum()
                })
                .collect();

            // Total LFQ (MS1) intensity at global q-value filter per file (if LFQ is enabled)
            let total_lfq_intensity_per_file: Vec<f32> = if let Some(areas) = &areas {
                let mut total_lfq_intensities = Vec::new();
                for i in 0..filenames.len() {
                    let mut intensities = Vec::new();
                    for ((_id, decoy), quantified) in areas {
                        if !decoy && quantified.peak.q_value <= global_q_value_filter {
                            if let Some(intensity) = quantified.intensities[i] {
                                intensities.push(intensity as f32);
                            }
                        }
                    }
                    total_lfq_intensities.push(intensities.iter().sum());
                }
                total_lfq_intensities
            } else {
                vec![0.0; filenames.len()]
            };

            // Mmedian MS1 mass accuracy for each file, using feature.delta_mass
            let median_ms1_mass_accuracy_per_file: Vec<f32> = filenames
                .iter()
                .map(|filename| {
                    let mut accuracies = Vec::new();
                    for feature in features.iter().filter(|f| {
                        filenames[f.file_id] == *filename
                            && f.label == 1
                            && f.spectrum_q <= global_q_value_filter
                    }) {
                        accuracies.push(feature.delta_mass);
                    }
                    accuracies.sort_by(|a, b| a.partial_cmp(b).unwrap());
                    let mid = accuracies.len() / 2;

                    if accuracies.is_empty() {
                        return f32::NAN;
                    }

                    if accuracies.len() % 2 == 0 {
                        if mid > 0 {
                            (accuracies[mid - 1] + accuracies[mid]) / 2.0
                        } else {
                            accuracies[mid]
                        }
                    } else {
                        accuracies[mid]
                    }
                })
                .collect();

            // Median MS2 mass accuracy for each file, using feature.average_ppm
            let median_ms2_mass_accuracy_per_file: Vec<f32> = filenames
                .iter()
                .map(|filename| {
                    let mut accuracies = Vec::new();
                    for feature in features.iter().filter(|f| {
                        filenames[f.file_id] == *filename
                            && f.label == 1
                            && f.spectrum_q <= global_q_value_filter
                    }) {
                        accuracies.push(feature.average_ppm);
                    }
                    accuracies.sort_by(|a, b| a.partial_cmp(b).unwrap());
                    let mid = accuracies.len() / 2;

                    if accuracies.is_empty() {
                        return f32::NAN;
                    }

                    if accuracies.len() % 2 == 0 {
                        if mid > 0 {
                            (accuracies[mid - 1] + accuracies[mid]) / 2.0
                        } else {
                            accuracies[mid]
                        }
                    } else {
                        accuracies[mid]
                    }
                })
                .collect();

            // Median RT deviation for each file, using feature.delta_rt_model
            let median_rt_deviation_per_file: Vec<f32> = filenames
                .iter()
                .map(|filename| {
                    let mut deviations = Vec::new();
                    for feature in features.iter().filter(|f| {
                        filenames[f.file_id] == *filename
                            && f.label == 1
                            && f.spectrum_q <= global_q_value_filter
                    }) {
                        deviations.push(feature.delta_rt_model);
                    }
                    deviations.sort_by(|a, b| a.partial_cmp(b).unwrap());
                    let mid = deviations.len() / 2;

                    if deviations.is_empty() {
                        return f32::NAN;
                    }

                    if deviations.len() % 2 == 0 {
                        if mid > 0 {
                            (deviations[mid - 1] + deviations[mid]) / 2.0
                        } else {
                            deviations[mid]
                        }
                    } else {
                        deviations[mid]
                    }
                })
                .collect();

            // Median IM deviation for each file, using feature.delta_ims_model
            let median_im_deviation_per_file: Vec<f32> = filenames
                .iter()
                .map(|filename| {
                    let mut deviations = Vec::new();
                    for feature in features.iter().filter(|f| {
                        filenames[f.file_id] == *filename
                            && f.label == 1
                            && f.spectrum_q <= global_q_value_filter
                    }) {
                        deviations.push(feature.delta_ims_model);
                    }
                    deviations.sort_by(|a, b| a.partial_cmp(b).unwrap());
                    let mid = deviations.len() / 2;

                    if deviations.is_empty() {
                        return f32::NAN;
                    }

                    if deviations.len() % 2 == 0 {
                        if mid > 0 {
                            (deviations[mid - 1] + deviations[mid]) / 2.0
                        } else {
                            deviations[mid]
                        }
                    } else {
                        deviations[mid]
                    }
                })
                .collect();

            // Average peptide length for each file
            let avg_peptide_length_per_file: Vec<f32> = filenames
                .iter()
                .map(|filename| {
                    let mut lengths = Vec::new();
                    for feature in features.iter().filter(|f| {
                        filenames[f.file_id] == *filename
                            && f.label == 1
                            && f.spectrum_q <= global_q_value_filter
                    }) {
                        lengths.push(feature.peptide_len as f32);
                    }
                    lengths.iter().sum::<f32>() / lengths.len() as f32
                })
                .collect();

            // Average peptide charge for each file
            let avg_peptide_charge_per_file: Vec<f32> = filenames
                .iter()
                .map(|filename| {
                    let mut charges = Vec::new();
                    for feature in features.iter().filter(|f| {
                        filenames[f.file_id] == *filename
                            && f.label == 1
                            && f.spectrum_q <= global_q_value_filter
                    }) {
                        charges.push(feature.charge as f32);
                    }
                    charges.iter().sum::<f32>() / charges.len() as f32
                })
                .collect();

            // Average number of matched peaks for each file
            let avg_matched_peaks_per_file: Vec<f32> = filenames
                .iter()
                .map(|filename| {
                    let mut peaks = Vec::new();
                    for feature in features.iter().filter(|f| {
                        filenames[f.file_id] == *filename
                            && f.label == 1
                            && f.spectrum_q <= global_q_value_filter
                    }) {
                        peaks.push(feature.matched_peaks as f32);
                    }
                    peaks.iter().sum::<f32>() / peaks.len() as f32
                })
                .collect();

            // Prepare html table to add to the report
            let table = html! {
                div class="table-container" {
                    table id="dataTable"  class="display" {
                        thead {
                            tr {
                                th { "File" }
                                th { "PSMs" }
                                th { "Peptides" }
                                th { "Proteins" }
                                th { "Total MS1 Intensity" }
                                th { "Total MS2 Intensity" }
                                th { "Median MS1 Delta Mass" }
                                th { "Median MS2 Delta Mass" }
                                th { "Median RT Deviation" }
                                th { "Median IM Deviation" }
                                th { "Average Peptide Length" }
                                th { "Average Peptide Charge" }
                                th { "Average Matched Peaks" }
                            }
                        }
                        tbody {
                            @for (i, filename) in filenames.iter().enumerate() {
                                tr {
                                    td { (filename) }
                                    td { (num_psm_targets_per_file[i]) }
                                    td { (num_peptide_targets_per_file[i]) }
                                    td { (num_protein_targets_per_file[i]) }
                                    td { (total_lfq_intensity_per_file[i]) }
                                    td { (total_ms2_intensity_per_file[i]) }
                                    td { (median_ms1_mass_accuracy_per_file[i]) }
                                    td { (median_ms2_mass_accuracy_per_file[i]) }
                                    td { (median_rt_deviation_per_file[i]) }
                                    td { (median_im_deviation_per_file[i]) }
                                    td { (avg_peptide_length_per_file[i]) }
                                    td { (avg_peptide_charge_per_file[i]) }
                                    td { (avg_matched_peaks_per_file[i]) }
                                }
                            }
                        }
                    }
                    button id="downloadCsv" { "Download as CSV" }
                }
            };

            intro_section.add_content(table);

            // Add boxplot of the LFQ intensities from areas if available
            if let Some(areas) = areas {
                let mut lfq_intensities: Vec<Vec<f64>> = Vec::new();
                for i in 0..filenames.len() {
                    let mut intensities = Vec::new();
                    for ((_id, decoy), quantified) in &areas {
                        if !decoy && quantified.peak.q_value <= global_q_value_filter {
                            if let Some(intensity) = quantified.intensities[i] {
                                intensities.push(intensity.log2());
                            }
                        }
                    }
                    lfq_intensities.push(intensities);
                }

                let lfq_boxplot = plot_boxplot(
                    &lfq_intensities,
                    filenames.to_vec(),
                    &format!("LFQ Intensities ({:?}% Q-value)", global_q_value_filter),
                    "Run",
                    "Log2(Intensity)",
                )
                .unwrap();
                intro_section.add_plot(lfq_boxplot);
            }

            report.add_section(intro_section);
        }

        /* Section 2: Scoring QC */
        {
            let mut scoring_section = ReportSection::new("Scoring Quality Control");

            scoring_section.add_content(html! {
                "It is important to assess the quality of the scoring model to ensure that the model is performing as expected, and that we're not overfitting or violating any assumptions of the Target-Decoy approach. The plot below shows the distribution of discriminant scores for each PSM, colored by whether the PSM is a target or decoy. We would expect the target distributions to be bimodal, where the first mode represents false targets that should align with the decoy distribution, and the second mode represents true targets."
            });

            // Extract sage_discriminant_score and label from features
            let (scores, labels): (Vec<f64>, Vec<i32>) = features
                .iter()
                .map(|f| (f.discriminant_score as f64, f.label))
                .unzip();

            if !scores.is_empty() && scores.len() > 100 {
                let score_histogram =
                    plot_score_histogram(&scores, &labels, "LDA Score", "Score").unwrap();

                scoring_section.add_plot(score_histogram);

                let pp_plot = plot_pp(&scores, &labels, "PP Plot").unwrap();

                scoring_section.add_content(html! {
                    "The Probability-Probability (PP) plot is a diagnostic tool that can be used to assess the quality of the scoring model. It plots the empirical cumulative distribution function (ECDF) of the target distribution against the ECDF of the decoy distribution. See: Debrie, E. et. al. (2023) Journal of Proteome Research. for more information."
                });
                scoring_section.add_plot(pp_plot);

                let spectrum_q_histogram = plot_score_histogram(
                    &features
                        .iter()
                        .map(|f| f.spectrum_q as f64)
                        .collect::<Vec<f64>>(),
                    &labels,
                    "Spectrum Q-value",
                    "Q-value",
                )
                .unwrap();
                scoring_section.add_plot(spectrum_q_histogram);

                let peptide_q_histogram = plot_score_histogram(
                    &features
                        .iter()
                        .map(|f| f.peptide_q as f64)
                        .collect::<Vec<f64>>(),
                    &labels,
                    "Peptide Q-value",
                    "Q-value",
                )
                .unwrap();
                scoring_section.add_plot(peptide_q_histogram);

                let protein_q_histogram = plot_score_histogram(
                    &features
                        .iter()
                        .map(|f| f.protein_q as f64)
                        .collect::<Vec<f64>>(),
                    &labels,
                    "Protein Q-value",
                    "Q-value",
                )
                .unwrap();
                scoring_section.add_plot(protein_q_histogram);
            } else {
                scoring_section.add_content(html! {
                    div style="margin-top: 10px; margin-bottom: 10px; padding: 15px; background-color: #ffe6e6; border: 1px solid #ff9999; color: #cc0000; border-radius: 5px; white-space: pre-line;" {
                        p {
                            "There are not enough scores to plot the scoring quality control plots."
                        }
                    }
                });
            }

            report.add_section(scoring_section);
        }

        /* Section 3: Predicted Properties */
        {
            let mut predicted_properties_section = ReportSection::new("Predicted Properties");

            predicted_properties_section.add_content(html! {
                "The following plots show the predicted properties of target peptides. The plots show the predicted retention time and ion mobility if present. The predicted properties are used to assess the quality of the model and to identify potential outliers."
            });

            // Normalized experimental RT per file
            let mut rt_per_file: Vec<Vec<f64>> = Vec::new();
            for i in 0..filenames.len() {
                let mut rts = Vec::new();
                for feature in features.iter().filter(|f| {
                    f.label == 1
                        && f.spectrum_q <= predict_section_q_value_filter
                        && filenames[f.file_id] == filenames[i]
                }) {
                    rts.push(feature.rt as f64);
                }

                let min_rt = rts.iter().cloned().fold(f64::INFINITY, f64::min);
                let max_rt = rts.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
                rts = rts
                    .iter()
                    .map(|rt| (rt - min_rt) / (max_rt - min_rt))
                    .collect();

                rt_per_file.push(rts);
            }

            // Predicted RT per file
            let mut predicted_rt_per_file: Vec<Vec<f64>> = Vec::new();
            for i in 0..filenames.len() {
                let mut predicted_rts = Vec::new();
                for feature in features.iter().filter(|f| {
                    f.label == 1
                        && f.spectrum_q <= predict_section_q_value_filter
                        && filenames[f.file_id] == filenames[i]
                }) {
                    predicted_rts.push(feature.predicted_rt as f64);
                }
                predicted_rt_per_file.push(predicted_rts);
            }

            let rt_scatter = plot_scatter(
                &rt_per_file,
                &predicted_rt_per_file,
                filenames.to_vec(),
                "Retention Time LR Model",
                "Retention Time",
                "Predicted Retention Time",
            )
            .unwrap();
            predicted_properties_section.add_plot(rt_scatter);

            // Experimental IMS per file
            let mut ims_per_file: Vec<Vec<f64>> = Vec::new();
            for i in 0..filenames.len() {
                let mut imss = Vec::new();
                for feature in features.iter().filter(|f| {
                    f.label == 1
                        && f.spectrum_q <= predict_section_q_value_filter
                        && filenames[f.file_id] == filenames[i]
                }) {
                    imss.push(feature.ims as f64);
                }

                ims_per_file.push(imss);
            }

            // Predicted IMS per file
            let mut predicted_ims_per_file: Vec<Vec<f64>> = Vec::new();
            for i in 0..filenames.len() {
                let mut predicted_imss = Vec::new();
                for feature in features.iter().filter(|f| {
                    f.label == 1
                        && f.spectrum_q <= predict_section_q_value_filter
                        && filenames[f.file_id] == filenames[i]
                }) {
                    predicted_imss.push(feature.predicted_ims as f64);
                }
                predicted_ims_per_file.push(predicted_imss);
            }

            if !ims_per_file.is_empty() && !predicted_ims_per_file.is_empty() {
                let ims_scatter = plot_scatter(
                    &ims_per_file,
                    &predicted_ims_per_file,
                    filenames.to_vec(),
                    "Ion Mobility LR Model",
                    "Ion Mobility",
                    "Predicted Ion Mobility",
                )
                .unwrap();
                predicted_properties_section.add_plot(ims_scatter);
            }

            report.add_section(predicted_properties_section);
        }

        /* Section 4: Configuration */
        {
            let mut config_section = ReportSection::new("Configuration");
            config_section.add_content(html! {
                style {
                    ".code-container {
                        background-color: #f5f5f5;
                        padding: 10px;
                        border-radius: 5px;
                        overflow-x: auto;
                        font-family: monospace;
                        white-space: pre-wrap;
                    }"
                }
                div class="code-container" {
                    pre {
                        code { (PreEscaped(serde_json::to_string_pretty(&self.parameters)?)) }
                    }
                }
            });
            report.add_section(config_section);
        }

        let bytes = report.to_string().into_bytes();
        sage_cloudpath::write_bytes_sync(&path, bytes)?;

        Ok(path)
    }
}

#[cfg(test)]
mod tests {
    use super::{
        assign_psm_ids, passes_localization_filter, passes_output_filter,
        sort_features_by_discriminant, RunSummary,
    };
    use sage_core::database::PeptideIx;
    use sage_core::scoring::Feature;

    #[test]
    fn tied_features_receive_repeatable_psm_ids() {
        let feature = |file_id, spec_id: &str, peptide_idx| Feature {
            file_id,
            spec_id: spec_id.into(),
            rank: 1,
            peptide_idx: PeptideIx(peptide_idx),
            discriminant_score: 5.0,
            ..Feature::default()
        };
        let mut forward = vec![feature(1, "scan=2", 2), feature(0, "scan=1", 1)];
        let mut reversed = forward.iter().cloned().rev().collect::<Vec<_>>();

        sort_features_by_discriminant(&mut forward);
        assign_psm_ids(&mut forward);
        sort_features_by_discriminant(&mut reversed);
        assign_psm_ids(&mut reversed);

        let identities = |features: &[Feature]| {
            features
                .iter()
                .map(|feature| (feature.file_id, feature.spec_id.clone(), feature.psm_id))
                .collect::<Vec<_>>()
        };
        assert_eq!(identities(&forward), identities(&reversed));
    }

    #[test]
    fn localization_filter_requires_passing_target_psm() {
        let passing = Feature {
            label: 1,
            spectrum_q: 0.01,
            ..Default::default()
        };
        assert!(passes_localization_filter(&passing, 0.01));

        let failing = Feature {
            spectrum_q: 0.011,
            ..passing.clone()
        };
        assert!(!passes_localization_filter(&failing, 0.01));

        let decoy = Feature {
            label: -1,
            ..passing
        };
        assert!(!passes_localization_filter(&decoy, 0.01));
    }

    #[test]
    fn output_filter_is_inclusive_and_applies_to_targets_and_decoys() {
        let target = Feature {
            label: 1,
            spectrum_q: 0.1,
            ..Default::default()
        };
        assert!(passes_output_filter(&target, 0.1));

        let decoy = Feature {
            label: -1,
            ..target.clone()
        };
        assert!(passes_output_filter(&decoy, 0.1));

        let failing = Feature {
            spectrum_q: 0.100_001,
            ..target
        };
        assert!(!passes_output_filter(&failing, 0.1));
    }

    #[test]
    fn older_run_summaries_receive_compatible_defaults() {
        let summary: RunSummary = serde_json::from_value(serde_json::json!({
            "runtime_secs": 1,
            "files": 1,
            "peptides_in_database": 10,
            "fragments_in_database": 20,
            "psms_at_one_percent_fdr": 2,
            "peptides_at_one_percent_fdr": 1,
            "proteins_at_one_percent_fdr": 1,
            "protein_groups_at_one_percent_fdr": 1,
            "output_paths": []
        }))
        .unwrap();

        assert_eq!(summary.schema_version, 1);
        assert!(!summary.ptm_localization.enabled);
        assert_eq!(summary.models.library_retention_time_alignment, None);
        assert_eq!(summary.models.library_retention_time_files_aligned, 0);
        assert_eq!(summary.models.library_ion_mobility_alignment, None);
        assert_eq!(summary.models.library_ion_mobility_files_aligned, 0);
        assert_eq!(summary.models.library_rescoring, None);
        assert_eq!(summary.quantification.lfq_features, 0);
    }
}
