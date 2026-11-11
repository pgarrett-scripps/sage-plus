use super::input::Search;
use super::memory::{trim_allocator, AllocatorTrimResult, MemoryLimits};
use super::output::SageResults;
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
use sage_core::mass::PROTON;
use sage_core::mass_calibration::{
    align_fragment_error, fit as fit_mass_calibration, CalibrationPoint, FitOptions,
    ToleranceEstimate, MIN_TOLERANCE_PSMS,
};
use sage_core::mass_recalibration::{
    confident_per_group, select_group_models, select_model, stable_hash, stratified_sample,
    FileMassCorrection, GroupMassCorrection, MassErrorPoint, MassModelKind, MassRecalibration,
    MassRecalibrationMode, ModelSelection, RecalibrationOptions,
};
use sage_core::peptide::Peptide;
use sage_core::scoring::{AtomicBitSet, Feature, Scorer};
use sage_core::spectral_library::{
    self, LibrarySelection, SpectralLibraryFormat, SpectralLibraryStrategy,
};
use sage_core::spectrum::{AcquisitionGroup, ProcessedSpectrum, SpectrumProcessor};
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

#[cfg(test)]
fn close_prefilter_pairs(database: &IndexedDatabase, keep: &AtomicBitSet) {
    close_pairs(keep, |index| {
        database
            .paired_peptide_index(sage_core::database::PeptideIx(index as u32))
            .map(|pair| pair.0 as usize)
    });
}

/// Add the target or decoy partner, given by `pair_of`, of every kept peptide.
#[cfg(test)]
fn close_pairs(keep: &AtomicBitSet, pair_of: impl Fn(usize) -> Option<usize> + Sync) {
    // Pairs are resolved in parallel and repeated until no partner is added,
    // so the closure does not depend on visit order.
    let mut frontier = (0..keep.len())
        .into_par_iter()
        .filter(|&index| keep.contains(index))
        .collect::<Vec<_>>();
    while !frontier.is_empty() {
        frontier = frontier
            .into_par_iter()
            .filter_map(|index| {
                let pair = pair_of(index)?;
                (!keep.contains(pair)).then_some(pair)
            })
            .collect::<Vec<_>>();
        frontier.sort_unstable();
        frontier.dedup();
        for &pair in &frontier {
            keep.insert(pair);
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

pub struct Runner {
    pub database: IndexedDatabase,
    pub parameters: Search,
    database_parameters: Parameters,
    start: Instant,
    events: EventEmitter,
    cancellation: CancellationToken,
    /// Spectra read by the prefilter, kept for the search so each file is
    /// read and processed once.
    retained_spectra: std::sync::Mutex<RetainedSpectra>,
    /// Search-time mass corrections selected for each searched file.
    mass_recalibration: std::sync::Mutex<Vec<MassRecalibrationFileStats>>,
    /// Quality-control scans of each file's raw spectra, by file id.
    file_qc: Arc<std::sync::Mutex<std::collections::BTreeMap<usize, qc::FileQc>>>,
    /// Custom cleavage sites of the searched FASTA, for the digestion summary.
    custom_cleavages: Option<ValidatedCustomCleavageLibrary>,
    /// Provenance key-value pairs for Parquet footers, set before outputs are
    /// written.
    parquet_provenance: Vec<(String, String)>,
}

/// Processed MS1 and MSn spectra of one file batch.
pub(crate) type SpectrumBatch = (Vec<ProcessedSpectrum>, Vec<ProcessedSpectrum>);

/// File batches retained from the prefilter, indexed by batch. They are only
/// valid for a search that splits the files into batches of the same size.
#[derive(Default)]
pub(crate) struct RetainedSpectra {
    batch_size: usize,
    batches: Vec<Option<SpectrumBatch>>,
}

/// A single localized modification site for one PSM, used to build the
/// PTM-site and protein-site reports.
struct SiteRow {
    /// Site from a decoy PSM: counted in site-level FDR, never written.
    decoy: bool,
    /// Discriminant score of the supporting PSM.
    discriminant_score: f32,
    ambiguous: bool,
    protein_sites: Arc<[sage_core::enzyme::ProteinOccurrence]>,
    attachment: sage_core::ptm_library::Attachment,
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
    /// `modification` is a name, not a signed mass.
    named: bool,
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

/// PSMs entering PTM localization and site-level FDR. Targets and decoys pass
/// the same identification cutoff so decoy sites can be counted.
fn passes_site_psm_filter(feature: &Feature, psm_q_value: f32) -> bool {
    feature.spectrum_q <= psm_q_value
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

/// Zero-based occurrence of each spectrum's `(file_id, id)` among the
/// preceding spectra. Unique IDs are always 0. Readers return spectra in file
/// order, so the numbering is identical between the search pass and the
/// post-FDR reread of the same files.
fn spectrum_id_occurrences(spectra: &[ProcessedSpectrum]) -> Vec<usize> {
    let mut seen = HashMap::with_capacity(spectra.len());
    spectra
        .iter()
        .map(|spectrum| {
            let count = seen
                .entry((spectrum.file_id, spectrum.id.as_str()))
                .or_insert(0usize);
            *count += 1;
            *count - 1
        })
        .collect()
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

/// Load `database.ptm_library`, if configured, into the database parameters.
/// Does a PTM library record's residue match the FASTA residue at its
/// position? With ambiguous-residue expansion on, a B, Z or X in the FASTA
/// also matches any residue it expands to, since the record is applied to
/// that expanded peptide.
fn library_residue_matches(fasta: Option<u8>, record: u8, expand_ambiguous: bool) -> bool {
    fasta.is_some_and(|written| {
        written == record
            || (expand_ambiguous
                && sage_core::ambiguous_residues::is_expansion_of(&[written], &[record]))
    })
}

fn load_ptm_library(database_parameters: &mut Parameters) -> anyhow::Result<()> {
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
    Ok(())
}

fn load_fasta(database_parameters: &Parameters) -> anyhow::Result<Fasta> {
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
    database_parameters.log_ambiguous_proteins(&fasta);
    Ok(fasta)
}

fn load_custom_cleavages(
    database_parameters: &Parameters,
    fasta: &Fasta,
) -> anyhow::Result<Option<ValidatedCustomCleavageLibrary>> {
    let Some(path) = database_parameters.custom_cleavage_sites.as_deref() else {
        return Ok(None);
    };
    let library = if path.to_ascii_lowercase().ends_with(".parquet") {
        let content = sage_cloudpath::util::read_bytes(path)
            .with_context(|| format!("Failed to read custom cleavage-site file `{path}`"))?;
        sage_cloudpath::parquet::deserialize_custom_cleavage_sites(content)
            .with_context(|| format!("Failed to parse custom cleavage-site file `{path}`"))?
    } else {
        let content = sage_cloudpath::util::read_text(path)
            .with_context(|| format!("Failed to read custom cleavage-site file `{path}`"))?;
        CustomCleavageLibrary::from_tsv(&content)
            .with_context(|| format!("Failed to parse custom cleavage-site file `{path}`"))?
    };
    let validated = library
        .validate(fasta)
        .with_context(|| format!("Failed to validate custom cleavage-site file `{path}`"))?;
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
    Ok(Some(validated))
}

fn empty_database_message(parameters: &sage_core::database::Parameters) -> String {
    let enzyme = &parameters.enzyme;
    let mut message = format!(
        "the database contains no target peptides. Check that `database.fasta` or \
         `database.peptides` contains target sequences, and that peptides of length {}-{} \
         and mass {}-{} Da can be digested from them",
        enzyme.min_len.unwrap_or(5),
        enzyme.effective_max_len(),
        parameters.peptide_min_mass,
        parameters.peptide_max_mass,
    );
    if parameters.prefilter {
        message.push_str(
            ", or that some peptides match the spectra closely enough to pass the prefilter",
        );
    }
    message
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
        load_ptm_library(&mut database_parameters)?;
        events.emit(EventKind::DatabaseStarted);
        let limits = MemoryLimits::from_gib(parameters.max_memory_gb)?;
        let mut retained_spectra = RetainedSpectra::default();
        let file_qc = Arc::new(std::sync::Mutex::new(std::collections::BTreeMap::new()));
        // Prefilter survivors are already mass-ordered and deduplicated.
        let mut reordered = false;
        let mut custom_cleavage_sites = None;
        // Collect peptides from FASTA (if configured).
        let mut all_peptides: Vec<Peptide> = if !database_parameters.fasta.is_empty() {
            let fasta = load_fasta(&database_parameters)?;
            let custom_cleavages = load_custom_cleavages(&database_parameters, &fasta)?;
            custom_cleavage_sites = custom_cleavages.clone();

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
                            if !library_residue_matches(
                                sequence.get(site.position as usize).copied(),
                                site.residue,
                                database_parameters.expand_ambiguous_residues,
                            ) =>
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

            database_parameters.log_ambiguous_expansion(&fasta, custom_cleavages.as_ref());
            match database_parameters.prefilter {
                false => {
                    let digests = database_parameters
                        .digest_unmodified_with_custom_cleavages(&fasta, custom_cleavages.as_ref());
                    database_parameters.modify_digests(digests)
                }
                true => {
                    let mini_runner = Self {
                        database: IndexedDatabase::default(),
                        parameters: parameters.clone(),
                        database_parameters: database_parameters.clone(),
                        start,
                        events: events.clone(),
                        cancellation: cancellation.clone(),
                        retained_spectra: Default::default(),
                        mass_recalibration: Default::default(),
                        file_qc: file_qc.clone(),
                        custom_cleavages: None,
                        parquet_provenance: Vec::new(),
                    };
                    let (peptides, retained) =
                        mini_runner.prefilter_peptides(parallel, fasta, custom_cleavages)?;
                    retained_spectra = retained;
                    reordered = true;
                    peptides
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
            reordered = false;
        }

        // Merge, deduplicate, and build the index.
        if !reordered {
            database_parameters.reorder_merged_peptides(&mut all_peptides);
        }
        // Prefilter spectra are kept to skip rereading them, unless the
        // fragment index looks unlikely to fit beside them. The search
        // rereads released spectra.
        if limits.is_enabled() && !retained_spectra.batches.is_empty() {
            let estimate = database_parameters.estimate_index_memory(&all_peptides);
            let additional = estimate
                .fragment_peak_bytes
                .saturating_sub(estimate.modified_peak_bytes);
            if !limits.estimate_fits(additional) {
                info!("releasing prefilter spectra to make room for the fragment index");
                retained_spectra = RetainedSpectra::default();
                let _ = trim_allocator();
            }
        }
        let database = database_parameters
            .clone()
            .build_from_peptides(all_peptides);
        anyhow::ensure!(
            database.peptides.iter().any(|peptide| !peptide.decoy),
            "{}",
            empty_database_message(&database_parameters)
        );

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
            retained_spectra: std::sync::Mutex::new(retained_spectra),
            mass_recalibration: Default::default(),
            file_qc,
            custom_cleavages: custom_cleavage_sites,
            parquet_provenance: Vec::new(),
        })
    }
}
mod artifacts;
pub mod estimate;
mod execution;
mod immonium;
mod postprocess;
pub(crate) mod prefilter;
pub mod qc;
mod search;
mod summary;
mod util;

pub use summary::*;
use util::*;

#[cfg(test)]
#[path = "../tests/unit/runner.rs"]
mod tests;
