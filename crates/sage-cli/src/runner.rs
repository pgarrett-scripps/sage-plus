use super::input::Search;
use super::memory::{trim_allocator, AllocatorTrimResult, MemoryLimits};
use super::output::SageResults;
use super::telemetry;
use crate::events::{CancellationToken, EventEmitter, EventKind};
use anyhow::Context;
use log::{info, warn};
use rayon::prelude::*;
use sage_cloudpath::{FileFormat, Url};
use sage_core::cleavage::{CustomCleavageLibrary, ValidatedCustomCleavageLibrary};
use sage_core::database::{IndexedDatabase, Parameters};
use sage_core::fasta::Fasta;
use sage_core::lfq::{PrecursorId, QuantifiedPeak};
use sage_core::mass::Tolerance;
use sage_core::mass_calibration::{
    align_fragment_error, fit as fit_mass_calibration, CalibrationPoint, FitOptions,
};
use sage_core::peptide::Peptide;
use sage_core::scoring::{AtomicBitSet, Feature, Scorer};
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

enum OutputTarget {
    Local(BufWriter<std::fs::File>),
    Remote(Box<BufWriter<sage_cloudpath::CloudWriter>>),
}

/// Numeric channel-group membership used only during prefilter closure. Group
/// strings are discarded after construction, leaving compact CSR arrays.
struct LabelGroupIndex {
    peptide_groups: Vec<u32>,
    offsets: Vec<u32>,
    members: Vec<sage_core::database::PeptideIx>,
}

impl LabelGroupIndex {
    fn new(peptides: &[Peptide]) -> Self {
        let mut groups = HashMap::<String, u32>::new();
        let mut peptide_groups = vec![u32::MAX; peptides.len()];
        let mut counts = Vec::<u32>::new();
        for (index, peptide) in peptides.iter().enumerate() {
            if peptide.label_channel.is_none() && peptide.label_group_override.is_none() {
                continue;
            }
            let next = groups.len() as u32;
            let group = *groups.entry(peptide.label_group()).or_insert(next);
            peptide_groups[index] = group;
            if group as usize == counts.len() {
                counts.push(0);
            }
            counts[group as usize] += 1;
        }

        let mut offsets = Vec::with_capacity(counts.len() + 1);
        offsets.push(0);
        for count in counts {
            offsets.push(offsets.last().copied().unwrap() + count);
        }
        let mut positions = offsets[..offsets.len().saturating_sub(1)].to_vec();
        let mut members = vec![
            sage_core::database::PeptideIx::default();
            offsets.last().copied().unwrap_or(0) as usize
        ];
        for (index, &group) in peptide_groups.iter().enumerate() {
            if group == u32::MAX {
                continue;
            }
            let position = &mut positions[group as usize];
            members[*position as usize] = sage_core::database::PeptideIx(index as u32);
            *position += 1;
        }
        Self {
            peptide_groups,
            offsets,
            members,
        }
    }

    fn close(&self, keep: &AtomicBitSet) {
        let selected_groups = AtomicBitSet::new(self.offsets.len().saturating_sub(1));
        for (index, &group) in self.peptide_groups.iter().enumerate() {
            if group != u32::MAX && keep.contains(index) {
                selected_groups.insert(group as usize);
            }
        }
        for group in 0..selected_groups.len() {
            if !selected_groups.contains(group) {
                continue;
            }
            let start = self.offsets[group] as usize;
            let end = self.offsets[group + 1] as usize;
            for peptide in &self.members[start..end] {
                keep.insert(peptide.0 as usize);
            }
        }
    }
}

fn close_prefilter_pairs(database: &IndexedDatabase, keep: &AtomicBitSet) {
    for index in 0..keep.len() {
        if !keep.contains(index) {
            continue;
        }
        if let Some(pair) =
            database.paired_peptide_index(sage_core::database::PeptideIx(index as u32))
        {
            keep.insert(pair.0 as usize);
        }
    }
}

impl OutputTarget {
    fn new(path: &Url) -> anyhow::Result<Self> {
        if path.scheme() == "file" {
            let local_path = path
                .to_file_path()
                .map_err(|_| anyhow::anyhow!("invalid local output URL `{path}`"))?;
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

fn median_finite(values: impl IntoIterator<Item = f32>) -> Option<f32> {
    let mut values = values
        .into_iter()
        .filter(|value| value.is_finite())
        .collect::<Vec<_>>();
    values.sort_by(f32::total_cmp);
    let middle = values.len() / 2;
    match values.len() {
        0 => None,
        length if length % 2 == 0 => Some((values[middle - 1] + values[middle]) / 2.0),
        _ => Some(values[middle]),
    }
}

fn average_finite(values: impl IntoIterator<Item = f32>) -> Option<f32> {
    let values = values
        .into_iter()
        .filter(|value| value.is_finite())
        .collect::<Vec<_>>();
    (!values.is_empty()).then(|| values.iter().sum::<f32>() / values.len() as f32)
}

fn normalize_finite(values: Vec<f64>) -> Vec<f64> {
    if values.is_empty() {
        return values;
    }
    let min = values.iter().copied().fold(f64::INFINITY, f64::min);
    let max = values.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    let range = max - min;
    if !range.is_finite() || range <= f64::EPSILON {
        return vec![0.0; values.len()];
    }
    values
        .into_iter()
        .map(|value| (value - min) / range)
        .collect()
}

fn labeled_finite_values(
    features: &[Feature],
    value: impl Fn(&Feature) -> f64,
) -> (Vec<f64>, Vec<i32>) {
    features
        .iter()
        .filter_map(|feature| {
            let value = value(feature);
            (value.is_finite() && matches!(feature.label, -1 | 1)).then_some((value, feature.label))
        })
        .unzip()
}

pub struct Runner {
    pub database: IndexedDatabase,
    pub parameters: Search,
    database_parameters: Parameters,
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

fn missing_decoy_warning(
    generate_decoys: bool,
    decoy_labels: impl IntoIterator<Item = bool>,
) -> Option<String> {
    if decoy_labels.into_iter().any(|decoy| decoy) {
        return None;
    }
    let remedy = if generate_decoys {
        "Check that target peptides can produce non-colliding reversed decoys or provide explicit decoys"
    } else {
        "Add decoys to the input database or set database.generate_decoys to true"
    };
    Some(format!(
        "the peptide database contains no decoys. FDR, q-values, rescoring, protein grouping, and LFQ filtering cannot be estimated reliably. {remedy}"
    ))
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
        let mut database_parameters = parameters.database.clone();
        let start = Instant::now();
        cancellation.check()?;
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
                            fasta
                                .targets
                                .len()
                                .div_ceil(database_parameters.prefilter_chunk_size),
                            database_parameters.prefilter_chunk_size,
                        );
                        let mini_runner = Self {
                            database: IndexedDatabase::default(),
                            parameters: parameters.clone(),
                            database_parameters: database_parameters.clone(),
                            start,
                            events: events.clone(),
                            cancellation: cancellation.clone(),
                        };
                        mini_runner.prefilter_peptides(parallel, fasta, custom_cleavages)?
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

        if let Some(message) = missing_decoy_warning(
            database_parameters.generate_decoys,
            database.peptides.iter().map(|peptide| peptide.decoy),
        ) {
            warn!("{message}");
            events.emit(EventKind::Warning {
                code: "database_without_decoys".into(),
                message,
            });
        }

        let trim_started = Instant::now();
        let trim_result = trim_allocator();
        let trim_elapsed = trim_started.elapsed();
        match trim_result {
            AllocatorTrimResult::Released => {
                log::debug!(
                    "released unused allocator pages after database construction in {trim_elapsed:#?}"
                );
            }
            AllocatorTrimResult::NoRelease => {
                log::trace!(
                    "allocator had no unused pages to release after database construction in {trim_elapsed:#?}"
                );
            }
            AllocatorTrimResult::Unsupported => {}
        }

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
        info!(
            "fragment index allocated {} bytes ({:.3} bytes per fragment)",
            database.fragments.allocated_bytes(),
            database.fragments.allocated_bytes() as f64 / database.fragments.len().max(1) as f64
        );
        Ok(Self {
            database,
            parameters,
            database_parameters,
            start,
            events,
            cancellation,
        })
    }
}
mod artifacts;
mod execution;
mod postprocess;
mod prefilter;
mod search;

#[cfg(test)]
#[path = "../tests/unit/runner.rs"]
mod tests;
