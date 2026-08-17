use super::input::Search;
use super::memory::MemoryLimits;
use super::output::SageResults;
use super::telemetry;
use anyhow::Context;
use csv::ByteRecord;
use log::info;
use rayon::prelude::*;
use sage_cloudpath::{FileFormat, Url};
use sage_core::database::{IndexedDatabase, Parameters};
use sage_core::fasta::Fasta;
use sage_core::ion_series::Kind;
use sage_core::lfq::{Peak, PrecursorId};
use sage_core::mass::Tolerance;
use sage_core::peptide::Peptide;
use sage_core::scoring::Fragments;
use sage_core::scoring::{Feature, Scorer};
use sage_core::spectrum::{ProcessedSpectrum, SpectrumProcessor};
use sage_core::tmt::TmtQuant;
use std::collections::{HashMap, HashSet};
use std::io::{BufWriter, Write};
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
    start: Instant,
}

/// A single localized modification site for one PSM, used to build the
/// PTM-site and protein-site reports.
struct SiteRow {
    psm_id: usize,
    filename: String,
    scannr: String,
    peptide: String,
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
    candidate_sites: usize,
    site_determining_matched: u32,
    site_determining_total: u32,
    site_probabilities: String,
}

fn passes_localization_filter(feature: &Feature, psm_q_value: f32) -> bool {
    feature.label == 1 && feature.spectrum_q <= psm_q_value
}

#[derive(Default)]
struct SpectrumAccumulator {
    pub ms1: Vec<ProcessedSpectrum>,
    pub msn: Vec<ProcessedSpectrum>,
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
        let mut parameters = parameters.clone();
        parameters.database.use_bitmap = parameters.use_bitmap;
        let start = Instant::now();
        let limits =
            MemoryLimits::from_gib(parameters.max_memory_gb, parameters.min_free_memory_gb)?;
        // Collect peptides from FASTA (if configured).
        let mut all_peptides: Vec<Peptide> = if !parameters.database.fasta.is_empty() {
            let fasta_url = sage_cloudpath::to_url(&parameters.database.fasta)?;
            let fasta = sage_cloudpath::util::read_fasta(
                &fasta_url,
                &parameters.database.decoy_tag,
                parameters.database.generate_decoys,
            )
            .with_context(|| {
                format!(
                    "Failed to build database from `{}`",
                    parameters.database.fasta
                )
            })?;

            let needs_estimate = limits.is_enabled()
                || (parameters.database.prefilter && parameters.database.prefilter_chunk_size == 0);
            if needs_estimate {
                let full_estimate = parameters.database.estimate_memory(&fasta);
                if parameters.database.prefilter && parameters.database.prefilter_chunk_size == 0 {
                    parameters.database.auto_calculate_prefilter_chunk_size(
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

                    if parameters.database.prefilter {
                        let mut modified_peak = 0u64;
                        let mut fragment_peak = 0u64;
                        for chunk in fasta.iter_chunks(parameters.database.prefilter_chunk_size) {
                            let estimate = parameters.database.estimate_memory(&chunk);
                            modified_peak = modified_peak.max(estimate.modified_peak_bytes);
                            fragment_peak = fragment_peak.max(estimate.fragment_peak_bytes);
                        }
                        limits.check_estimate("modified-peptide", modified_peak)?;
                        limits.check_estimate("fragment-index", fragment_peak)?;
                    }
                }
            }

            match parameters.database.prefilter {
                false => {
                    let digests = parameters.database.digest_unmodified(&fasta);
                    if limits.is_enabled() {
                        let estimate = parameters.database.estimate_modified_memory(&digests);
                        info!(
                            "modification preflight: {} unmodified peptides may expand to {} modified peptides ({:.2} GiB additional peak)",
                            estimate.unmodified_peptides,
                            estimate.modified_peptides,
                            estimate.modified_peak_bytes as f64 / (1024.0 * 1024.0 * 1024.0),
                        );
                        limits.check_estimate("modified-peptide", estimate.modified_peak_bytes)?;
                    }
                    parameters.database.modify_digests(digests)
                }
                true => {
                    if parameters.database.prefilter_chunk_size >= fasta.targets.len() {
                        parameters.database.digest(&fasta)
                    } else {
                        info!(
                            "using {} db chunks of size {}",
                            (fasta.targets.len() + parameters.database.prefilter_chunk_size - 1)
                                / parameters.database.prefilter_chunk_size,
                            parameters.database.prefilter_chunk_size,
                        );
                        let mini_runner = Self {
                            database: IndexedDatabase::default(),
                            parameters: parameters.clone(),
                            start,
                        };
                        mini_runner.prefilter_peptides(parallel, fasta)
                    }
                }
            }
        } else {
            vec![]
        };

        // Append peptides from TSV file (if configured), additive with FASTA.
        if let Some(peptides_path) = parameters.database.peptides.clone() {
            let content = sage_cloudpath::util::read_text(&peptides_path)
                .with_context(|| format!("Failed to read peptide file `{peptides_path}`"))?;
            all_peptides.extend(parameters.database.peptides_from_tsv(&content));
        }

        // Merge, deduplicate, and build the index.
        Parameters::reorder_peptides(&mut all_peptides);
        if limits.is_enabled() {
            let estimate = parameters.database.estimate_index_memory(&all_peptides);
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
        let database = parameters
            .database
            .clone()
            .build_from_peptides(all_peptides);

        info!(
            "generated {} fragments, {} peptides in {:#?}",
            database.fragments.len(),
            database.peptides.len(),
            (start.elapsed())
        );
        Ok(Self {
            database,
            parameters,
            start,
        })
    }

    pub fn prefilter_peptides(self, parallel: usize, fasta: Fasta) -> Vec<Peptide> {
        let spectra: Option<Vec<ProcessedSpectrum>> =
            match parallel >= self.parameters.mzml_paths.len() {
                true => Some(
                    self.read_processed_spectra(&self.parameters.mzml_paths, 0, 0)
                        .1,
                ),
                false => None,
            };

        let db_params = self.parameters.database.clone();
        // TODO: Don't generate decoys for fast searching
        // * if `generate_decoys` is used, we should re-generate at the end
        //  to ensure that picked-peptide conditions are used, otherwise,
        //  if the user supplied decoys in the fasta file, then we should retain them
        //
        // db_params.generate_decoys = false;

        let mut all_peptides: Vec<Peptide> = fasta
            .iter_chunks(self.parameters.database.prefilter_chunk_size)
            .enumerate()
            .flat_map(|(chunk_id, fasta_chunk)| {
                let start = Instant::now();
                info!("pre-filtering fasta chunk {}", chunk_id,);
                let mut db = db_params.clone().build(fasta_chunk);

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
                    annotate_matches: self.parameters.annotate_matches,
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
                    self.parameters.database.prefilter_low_memory,
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

    fn spectrum_fdr(&self, features: &mut Vec<Feature>) -> (usize, bool) {
        let percolator_succeeded = self.parameters.percolator.enabled
            && match sage_core::ml::percolator::score_psms(
                features,
                &self.parameters.percolator,
                self.parameters.precursor_tol,
            ) {
                Ok(()) => {
                    log::info!(
                        "- Percolator-style {:?} PSM rescoring completed",
                        self.parameters.percolator.model
                    );
                    true
                }
                Err(error) => {
                    log::warn!("Percolator-style rescoring failed ({error}); falling back to LDA");
                    false
                }
            };

        if !percolator_succeeded
            && sage_core::ml::linear_discriminant::score_psms(
                features,
                self.parameters.precursor_tol,
            )
            .is_none()
        {
            log::warn!("linear model fitting failed, falling back to heuristic discriminant score");
            features.par_iter_mut().for_each(|feat| {
                feat.discriminant_score = (-feat.poisson as f32).ln_1p() + feat.longest_y_pct / 3.0
            });
        }
        features.par_sort_unstable_by(|a, b| b.discriminant_score.total_cmp(&a.discriminant_score));
        if percolator_succeeded {
            // Rescoring can rerank multiple candidates from one spectrum. Keep
            // exactly one target-or-decoy winner before confidence estimation
            // and downstream picked-peptide/protein calculations.
            let mut spectra = HashSet::new();
            features.retain(|feature| spectra.insert((feature.file_id, feature.spec_id.clone())));
        }
        (
            sage_core::ml::qvalue::spectrum_q_value(features),
            percolator_succeeded,
        )
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
            .filter(|spec| spec.masses.len() >= self.parameters.min_peaks && spec.level == 2)
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
        let inner_closure = |(idx, path)| {
            let file_id = chunk_idx * batch_size + idx;
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
                    Ok(s.into_par_iter()
                        .map(|spectrum| sp.process(spectrum))
                        .collect::<SpectrumAccumulator>())
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

        let io_time = Instant::now() - start;
        info!("- file IO: {:8} ms", io_time.as_millis());

        (spectra.ms1, spectra.msn)
    }

    /// Re-read MS2 spectra and localize only target PSMs that passed the
    /// configured identification q-value. This keeps localization out of the
    /// search hot path without retaining every processed spectrum in memory.
    fn localize_features(&self, features: &mut [Feature], batch_size: usize) {
        let mut feature_indices: HashMap<usize, HashMap<String, Vec<usize>>> = HashMap::new();
        for (idx, feature) in features.iter().enumerate() {
            if passes_localization_filter(feature, self.parameters.ptm_localization.psm_q_value)
                && sage_core::ptm::has_localizable_modification(
                    &self.database[feature.peptide_idx],
                    &self.database.potential_mods,
                )
            {
                feature_indices
                    .entry(feature.file_id)
                    .or_default()
                    .entry(feature.spec_id.clone())
                    .or_default()
                    .push(idx);
            }
        }

        if feature_indices.is_empty() {
            log::info!("- PTM localization: no passing target PSMs");
            return;
        }

        let start = Instant::now();
        let mut localized = 0usize;
        for (chunk_idx, chunk) in self.parameters.mzml_paths.chunks(batch_size).enumerate() {
            let first_file_id = chunk_idx * batch_size;
            if !(first_file_id..first_file_id + chunk.len())
                .any(|file_id| feature_indices.contains_key(&file_id))
            {
                continue;
            }
            let spectra = self
                .read_processed_spectra_with_ms1(chunk, chunk_idx, batch_size, false)
                .1;
            let results = spectra
                .par_iter()
                .map(|spectrum| {
                    feature_indices
                        .get(&spectrum.file_id)
                        .and_then(|file| file.get(spectrum.id.as_str()))
                        .into_iter()
                        .flatten()
                        .filter_map(|&idx| {
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
                        .collect::<Vec<_>>()
                })
                .flatten()
                .collect::<Vec<_>>();

            localized += results.len();
            for (idx, localization) in results {
                features[idx].localization = Some(localization);
            }
        }

        log::info!(
            "- PTM localization: {} PSMs in {} ms",
            localized,
            start.elapsed().as_millis()
        );
    }

    pub fn batch_files(&self, scorer: &Scorer, batch_size: usize) -> SageResults {
        self.parameters
            .mzml_paths
            .chunks(batch_size)
            .enumerate()
            .map(|(chunk_idx, chunk)| self.process_chunk(scorer, chunk, chunk_idx, batch_size))
            .collect::<SageResults>()
    }

    pub fn run(mut self, parallel: usize, parquet: bool) -> anyhow::Result<telemetry::Telemetry> {
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
            annotate_matches: self.parameters.annotate_matches,
            mass_shift_ppm: self.parameters.mass_shift_ppm,
            score_type: self.parameters.score_type,
            use_bitmap: self.parameters.use_bitmap,
        };

        //Collect all results into a single container
        let mut outputs = self.batch_files(&scorer, parallel);

        let alignments = if self.parameters.predict_rt {
            // Poisson probability is usually the best single feature for refining FDR.
            // Take our set of 1% FDR filtered PSMs, and use them to train a linear
            // regression model for predicting retention time
            outputs
                .features
                .par_sort_unstable_by(|a, b| a.poisson.total_cmp(&b.poisson));
            sage_core::ml::qvalue::spectrum_q_value(&mut outputs.features);

            let alignments = sage_core::ml::retention_alignment::global_alignment(
                &mut outputs.features,
                self.parameters.mzml_paths.len(),
            );
            let _ = sage_core::ml::retention_model::predict(&self.database, &mut outputs.features);
            let _ = sage_core::ml::mobility_model::predict(&self.database, &mut outputs.features);
            Some(alignments)
        } else {
            None
        };

        let (q_spectrum, percolator_succeeded) = self.spectrum_fdr(&mut outputs.features);
        let q_peptide = if percolator_succeeded {
            sage_core::fdr::picked_peptide_tdc(&self.database, &mut outputs.features)
        } else {
            sage_core::fdr::picked_peptide(&self.database, &mut outputs.features)
        };
        // Protein FDR is based exclusively on proteotypic (unique, non-shared) peptides. Shared peptides
        // are reported with protein FDR = 1.0
        let q_protein = if percolator_succeeded {
            sage_core::fdr::picked_protein_tdc(&self.database, &mut outputs.features)
        } else {
            sage_core::fdr::picked_protein(&self.database, &mut outputs.features)
        };
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
        let q_protein_group = if percolator_succeeded {
            sage_core::fdr::picked_protein_group_tdc(&self.database, &mut outputs.features)
        } else {
            sage_core::fdr::picked_protein_group(&self.database, &mut outputs.features)
        };

        if self.parameters.ptm_localization.enabled {
            self.localize_features(&mut outputs.features, parallel);
        }

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

        let areas = alignments.and_then(|alignments| {
            if self.parameters.quant.lfq {
                log::trace!("performing LFQ");
                let mut areas = sage_core::lfq::build_feature_map(
                    self.parameters.quant.lfq_settings,
                    self.parameters.precursor_charge,
                    &outputs.features,
                )
                .quantify(&self.database, &outputs.ms1, &alignments);

                let q_precursor = sage_core::fdr::picked_precursor(&mut areas);

                log::info!("discovered {} target MS1 peaks at 5% FDR", q_precursor);
                Some(areas)
            } else {
                None
            }
        });

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

        // Write either a single parquet file, or multiple tsv files
        if parquet {
            log::warn!("parquet output format is currently unstable! There may be failures or schema changes!");

            let bytes = sage_cloudpath::parquet::serialize_features(
                &outputs.features,
                &outputs.quant,
                &filenames,
                &self.database,
            )?;

            let path = self.make_path("results.sage.parquet");
            sage_cloudpath::write_bytes_sync(&path, bytes)?;
            self.parameters.output_paths.push(path);

            if self.parameters.annotate_matches {
                let bytes =
                    sage_cloudpath::parquet::serialize_matched_fragments(&outputs.features)?;
                let path = self.make_path("matched_fragments.sage.parquet");
                sage_cloudpath::write_bytes_sync(&path, bytes)?;
                self.parameters.output_paths.push(path);
            }

            if let Some(areas) = &areas {
                let bytes =
                    sage_cloudpath::parquet::serialize_lfq(areas, &filenames, &self.database)?;

                let path = self.make_path("lfq.parquet");
                sage_cloudpath::write_bytes_sync(&path, bytes)?;
                self.parameters.output_paths.push(path);
            }
        } else {
            self.parameters
                .output_paths
                .push(self.write_features(&outputs.features, &filenames)?);

            if self.parameters.annotate_matches {
                self.parameters
                    .output_paths
                    .push(self.write_fragments(&outputs.features)?);
            }

            if !outputs.quant.is_empty() {
                self.parameters
                    .output_paths
                    .push(self.write_tmt(&outputs.quant, &filenames)?);
            }
            if let Some(areas) = &areas {
                self.parameters
                    .output_paths
                    .push(self.write_lfq(areas, &filenames)?);
            }
        }

        // PTM site reports follow the selected main output format.
        if self.parameters.ptm_localization.enabled {
            self.parameters.output_paths.push(self.write_ptm_sites(
                &outputs.features,
                &filenames,
                parquet,
            )?);
            self.parameters.output_paths.push(self.write_protein_sites(
                &outputs.features,
                &filenames,
                parquet,
            )?);
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
        println!("{}", serde_json::to_string_pretty(&self.parameters)?);

        let bytes = serde_json::to_vec_pretty(&self.parameters)?;
        sage_cloudpath::write_bytes_sync(&path, bytes)?;
        self.parameters.output_paths.push(path);

        let run_time = (Instant::now() - self.start).as_secs();
        info!("finished in {}s", run_time);
        info!("cite: \"Sage: An Open-Source Tool for Fast Proteomics Searching and Quantification at Scale\" https://doi.org/10.1021/acs.jproteome.3c00486");

        let telemetry = telemetry::Telemetry::new(
            self.parameters,
            self.database.peptides.len(),
            self.database.fragments.len(),
            parquet,
            run_time,
        );

        Ok(telemetry)
    }
    pub fn serialize_feature(&self, feature: &Feature, filenames: &[String]) -> csv::ByteRecord {
        let mut record = csv::ByteRecord::new();

        record.push_field(itoa::Buffer::new().format(feature.psm_id).as_bytes());

        let peptide = &self.database[feature.peptide_idx];
        record.push_field(peptide.to_string().as_bytes());
        record.push_field(feature.ambiguity_sequence.as_bytes());
        record.push_field(ryu::Buffer::new().format(feature.mass_shift).as_bytes());
        record.push_field(
            peptide
                .proteins(&self.database.decoy_tag, self.database.generate_decoys)
                .as_bytes(),
        );
        record.push_field(feature.protein_groups.as_deref().unwrap_or("").as_bytes());
        record.push_field(
            itoa::Buffer::new()
                .format(peptide.proteins.len())
                .as_bytes(),
        );
        record.push_field(
            itoa::Buffer::new()
                .format(feature.num_protein_groups)
                .as_bytes(),
        );
        record.push_field(filenames[feature.file_id].as_bytes());
        record.push_field(feature.spec_id.as_bytes());
        record.push_field(itoa::Buffer::new().format(feature.rank).as_bytes());
        record.push_field(itoa::Buffer::new().format(feature.label).as_bytes());
        record.push_field(ryu::Buffer::new().format(feature.expmass).as_bytes());
        record.push_field(ryu::Buffer::new().format(feature.calcmass).as_bytes());
        record.push_field(itoa::Buffer::new().format(feature.charge).as_bytes());
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
        record.push_field(ryu::Buffer::new().format(feature.delta_mass).as_bytes());
        record.push_field(ryu::Buffer::new().format(feature.average_ppm).as_bytes());
        record.push_field(ryu::Buffer::new().format(feature.hyperscore).as_bytes());
        record.push_field(ryu::Buffer::new().format(feature.delta_next).as_bytes());
        record.push_field(ryu::Buffer::new().format(feature.delta_best).as_bytes());
        record.push_field(ryu::Buffer::new().format(feature.rt).as_bytes());
        record.push_field(ryu::Buffer::new().format(feature.aligned_rt).as_bytes());
        record.push_field(ryu::Buffer::new().format(feature.predicted_rt).as_bytes());
        record.push_field(ryu::Buffer::new().format(feature.delta_rt_model).as_bytes());
        record.push_field(ryu::Buffer::new().format(feature.ims).as_bytes());
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
                .format(feature.matched_intensity_pct)
                .as_bytes(),
        );
        record.push_field(
            itoa::Buffer::new()
                .format(feature.scored_candidates)
                .as_bytes(),
        );
        record.push_field(ryu::Buffer::new().format(feature.poisson).as_bytes());
        record.push_field(
            ryu::Buffer::new()
                .format(feature.discriminant_score)
                .as_bytes(),
        );
        record.push_field(
            ryu::Buffer::new()
                .format(feature.posterior_error)
                .as_bytes(),
        );
        record.push_field(ryu::Buffer::new().format(feature.spectrum_q).as_bytes());
        record.push_field(ryu::Buffer::new().format(feature.peptide_q).as_bytes());
        record.push_field(ryu::Buffer::new().format(feature.protein_q).as_bytes());
        record.push_field(
            ryu::Buffer::new()
                .format(feature.protein_group_q)
                .as_bytes(),
        );
        record.push_field(ryu::Buffer::new().format(feature.ms2_intensity).as_bytes());
        record
    }

    pub fn serialize_fragments(
        &self,
        psm_id: usize,
        fragments_: &Option<Fragments>,
    ) -> Vec<ByteRecord> {
        let mut frag_records = vec![];

        if let Some(fragments) = fragments_ {
            for id in 0..fragments.fragment_ordinals.len() {
                let mut record = ByteRecord::new();
                record.push_field(itoa::Buffer::new().format(psm_id).as_bytes());
                let ion_type = match fragments.kinds[id] {
                    Kind::A => "a",
                    Kind::B => "b",
                    Kind::C => "c",
                    Kind::X => "x",
                    Kind::Y => "y",
                    Kind::Z => "z",
                };
                record.push_field(ion_type.as_bytes());
                record.push_field(
                    itoa::Buffer::new()
                        .format(fragments.fragment_ordinals[id])
                        .as_bytes(),
                );
                record.push_field(itoa::Buffer::new().format(fragments.charges[id]).as_bytes());
                record.push_field(
                    ryu::Buffer::new()
                        .format(fragments.mz_calculated[id])
                        .as_bytes(),
                );
                record.push_field(
                    ryu::Buffer::new()
                        .format(fragments.mz_experimental[id])
                        .as_bytes(),
                );
                record.push_field(
                    ryu::Buffer::new()
                        .format(fragments.neutral_losses[id])
                        .as_bytes(),
                );
                record.push_field(
                    ryu::Buffer::new()
                        .format(fragments.intensities[id])
                        .as_bytes(),
                );
                frag_records.push(record);
            }
        }

        frag_records
    }

    pub fn write_features(
        &self,
        features: &[Feature],
        filenames: &[String],
    ) -> anyhow::Result<Url> {
        let path = self.make_path("results.sage.tsv");

        let mut wtr = csv::WriterBuilder::new()
            .delimiter(b'\t')
            .from_writer(OutputTarget::new(&path)?);

        let csv_headers = vec![
            "psm_id",
            "peptide",
            "ambiguity_sequence",
            "mass_shift",
            "proteins",
            "protein_groups",
            "num_proteins",
            "num_protein_groups",
            "filename",
            "scannr",
            "rank",
            "label",
            "expmass",
            "calcmass",
            "charge",
            "peptide_len",
            "missed_cleavages",
            "semi_enzymatic",
            "isotope_error",
            "precursor_ppm",
            "fragment_ppm",
            "hyperscore",
            "delta_next",
            "delta_best",
            "rt",
            "aligned_rt",
            "predicted_rt",
            "delta_rt_model",
            "ion_mobility",
            "predicted_mobility",
            "delta_mobility",
            "matched_peaks",
            "longest_b",
            "longest_y",
            "longest_y_pct",
            "matched_intensity_pct",
            "scored_candidates",
            "poisson",
            "sage_discriminant_score",
            "posterior_error",
            "spectrum_q",
            "peptide_q",
            "protein_q",
            "protein_group_q",
            "ms2_intensity",
        ];

        let headers = csv::ByteRecord::from(csv_headers);

        wtr.write_byte_record(&headers)?;
        for chunk in features.chunks(1024) {
            for record in chunk
                .par_iter()
                .map(|feat| self.serialize_feature(feat, filenames))
                .collect::<Vec<_>>()
            {
                wtr.write_byte_record(&record)?;
            }
        }

        finish_csv_writer(wtr, &path)?;
        Ok(path)
    }

    pub fn write_fragments(&self, features: &[Feature]) -> anyhow::Result<Url> {
        let path = self.make_path("matched_fragments.sage.tsv");

        let mut wtr = csv::WriterBuilder::new()
            .delimiter(b'\t')
            .from_writer(OutputTarget::new(&path)?);

        let headers = csv::ByteRecord::from(vec![
            "psm_id",
            "fragment_type",
            "fragment_ordinals",
            "fragment_charge",
            "fragment_mz_calculated",
            "fragment_mz_experimental",
            "neutral_loss",
            "fragment_intensity",
        ]);

        wtr.write_byte_record(&headers)?;

        for chunk in features.chunks(1024) {
            for record in chunk
                .par_iter()
                .map(|feat| self.serialize_fragments(feat.psm_id, &feat.fragments))
                .flatten()
                .collect::<Vec<_>>()
            {
                wtr.write_byte_record(&record)?;
            }
        }

        finish_csv_writer(wtr, &path)?;
        Ok(path)
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
        parquet: bool,
    ) -> anyhow::Result<Url> {
        let rows = self.collect_site_rows(features, filenames);

        if parquet {
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
                    candidate_sites: row.candidate_sites as i32,
                    site_determining_ions_matched: row.site_determining_matched as i32,
                    site_determining_ions_total: row.site_determining_total as i32,
                    site_probabilities: row.site_probabilities.clone(),
                })
                .collect::<Vec<_>>();
            let path = self.make_path("results.sage.ptm-sites.parquet");
            let bytes = sage_cloudpath::parquet::serialize_ptm_sites(&records)?;
            sage_cloudpath::write_bytes_sync(&path, bytes)?;
            return Ok(path);
        }

        let path = self.make_path("results.sage.ptm-sites.tsv");

        let mut wtr = csv::WriterBuilder::new()
            .delimiter(b'\t')
            .from_writer(vec![]);

        wtr.write_byte_record(&csv::ByteRecord::from(vec![
            "psm_id",
            "filename",
            "scannr",
            "peptide",
            "proteins",
            "charge",
            "spectrum_q",
            "peptide_q",
            "modification",
            "modification_mass",
            "position",
            "residue",
            "localization_probability",
            "delta_localization_score",
            "candidate_sites",
            "site_determining_ions_matched",
            "site_determining_ions_total",
            "site_probabilities",
        ]))?;

        for row in &rows {
            let mut record = ByteRecord::new();
            record.push_field(itoa::Buffer::new().format(row.psm_id).as_bytes());
            record.push_field(row.filename.as_bytes());
            record.push_field(row.scannr.as_bytes());
            record.push_field(row.peptide.as_bytes());
            record.push_field(row.proteins.as_bytes());
            record.push_field(itoa::Buffer::new().format(row.charge).as_bytes());
            record.push_field(ryu::Buffer::new().format(row.spectrum_q).as_bytes());
            record.push_field(ryu::Buffer::new().format(row.peptide_q).as_bytes());
            record.push_field(row.modification.as_bytes());
            record.push_field(ryu::Buffer::new().format(row.modification_mass).as_bytes());
            record.push_field(itoa::Buffer::new().format(row.position).as_bytes());
            record.push_field([row.residue].as_slice());
            record.push_field(
                ryu::Buffer::new()
                    .format(row.localization_probability)
                    .as_bytes(),
            );
            record.push_field(ryu::Buffer::new().format(row.delta_score).as_bytes());
            record.push_field(itoa::Buffer::new().format(row.candidate_sites).as_bytes());
            record.push_field(
                itoa::Buffer::new()
                    .format(row.site_determining_matched)
                    .as_bytes(),
            );
            record.push_field(
                itoa::Buffer::new()
                    .format(row.site_determining_total)
                    .as_bytes(),
            );
            record.push_field(row.site_probabilities.as_bytes());
            wtr.write_byte_record(&record)?;
        }

        wtr.flush()?;
        let bytes = wtr.into_inner()?;
        sage_cloudpath::write_bytes_sync(&path, bytes)?;
        Ok(path)
    }

    /// Write a collapsed protein-site report: the best localization for each
    /// (protein, modified peptide site) aggregated across all supporting PSMs.
    pub fn write_protein_sites(
        &self,
        features: &[Feature],
        filenames: &[String],
        parquet: bool,
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
                    best_spectrum_q: f32::MAX,
                });
                entry.n_psms += 1;
                entry.best_probability = entry.best_probability.max(row.localization_probability);
                entry.best_delta_score = entry.best_delta_score.max(row.delta_score);
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

        if parquet {
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
                    best_spectrum_q: agg.best_spectrum_q,
                })
                .collect::<Vec<_>>();
            let path = self.make_path("results.sage.protein-sites.parquet");
            let bytes = sage_cloudpath::parquet::serialize_protein_sites(&records)?;
            sage_cloudpath::write_bytes_sync(&path, bytes)?;
            return Ok(path);
        }

        let path = self.make_path("results.sage.protein-sites.tsv");

        let mut wtr = csv::WriterBuilder::new()
            .delimiter(b'\t')
            .from_writer(vec![]);

        wtr.write_byte_record(&csv::ByteRecord::from(vec![
            "protein",
            "peptide",
            "residue",
            "position_in_peptide",
            "modification",
            "modification_mass",
            "num_psms",
            "best_localization_probability",
            "best_delta_localization_score",
            "best_spectrum_q",
        ]))?;

        for agg in &aggregated {
            let mut record = ByteRecord::new();
            record.push_field(agg.protein.as_bytes());
            record.push_field(agg.peptide.as_bytes());
            record.push_field([agg.residue].as_slice());
            record.push_field(itoa::Buffer::new().format(agg.position).as_bytes());
            record.push_field(agg.modification.as_bytes());
            record.push_field(ryu::Buffer::new().format(agg.modification_mass).as_bytes());
            record.push_field(itoa::Buffer::new().format(agg.n_psms).as_bytes());
            record.push_field(ryu::Buffer::new().format(agg.best_probability).as_bytes());
            record.push_field(ryu::Buffer::new().format(agg.best_delta_score).as_bytes());
            record.push_field(ryu::Buffer::new().format(agg.best_spectrum_q).as_bytes());
            wtr.write_byte_record(&record)?;
        }

        wtr.flush()?;
        let bytes = wtr.into_inner()?;
        sage_cloudpath::write_bytes_sync(&path, bytes)?;
        Ok(path)
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

    pub fn write_tmt(&self, quant: &[TmtQuant], filenames: &[String]) -> anyhow::Result<Url> {
        let path = self.make_path("tmt.tsv");

        let mut wtr = csv::WriterBuilder::new()
            .delimiter(b'\t')
            .from_writer(OutputTarget::new(&path)?);
        let mut headers = csv::ByteRecord::from(vec!["filename", "scannr", "ion_injection_time"]);
        headers.extend(
            self.parameters
                .quant
                .tmt
                .as_ref()
                .map(|tmt| tmt.headers())
                .expect("TMT quant cannot be performed without setting this parameter"),
        );

        wtr.write_byte_record(&headers)?;

        for chunk in quant.chunks(1024) {
            for record in chunk
                .par_iter()
                .map(|q| {
                    let mut record = csv::ByteRecord::new();
                    record.push_field(filenames[q.file_id].as_bytes());
                    record.push_field(q.spec_id.as_bytes());
                    record.push_field(ryu::Buffer::new().format(q.ion_injection_time).as_bytes());
                    for peak in &q.peaks {
                        record.push_field(ryu::Buffer::new().format(*peak).as_bytes());
                    }
                    record
                })
                .collect::<Vec<csv::ByteRecord>>()
            {
                wtr.write_record(&record)?;
            }
        }
        finish_csv_writer(wtr, &path)?;
        Ok(path)
    }

    pub fn write_lfq(
        &self,
        areas: &HashMap<(PrecursorId, bool), (Peak, Vec<f64>), fnv::FnvBuildHasher>,
        filenames: &[String],
    ) -> anyhow::Result<Url> {
        let path = self.make_path("lfq.tsv");

        let mut wtr = csv::WriterBuilder::new()
            .delimiter(b'\t')
            .from_writer(OutputTarget::new(&path)?);
        let mut headers = csv::ByteRecord::from(vec![
            "peptide",
            "charge",
            "proteins",
            "q_value",
            "score",
            "spectral_angle",
        ]);
        headers.extend(filenames);

        wtr.write_byte_record(&headers)?;

        for ((id, decoy), (peak, data)) in areas {
            if *decoy {
                continue;
            }
            let mut record = csv::ByteRecord::new();
            let (peptide_ix, charge) = match id {
                PrecursorId::Combined(x) => (*x, None),
                PrecursorId::Charged((x, charge)) => (*x, Some(*charge as i32)),
            };
            record.push_field(self.database[peptide_ix].to_string().as_bytes());
            record.push_field(itoa::Buffer::new().format(charge.unwrap_or(-1)).as_bytes());
            record.push_field(
                self.database[peptide_ix]
                    .proteins(&self.database.decoy_tag, self.database.generate_decoys)
                    .as_bytes(),
            );
            record.push_field(ryu::Buffer::new().format(peak.q_value).as_bytes());
            record.push_field(ryu::Buffer::new().format(peak.score).as_bytes());
            record.push_field(ryu::Buffer::new().format(peak.spectral_angle).as_bytes());
            for x in data {
                record.push_field(ryu::Buffer::new().format(*x).as_bytes());
            }
            wtr.write_record(&record)?;
        }
        finish_csv_writer(wtr, &path)?;
        Ok(path)
    }

    fn write_report(
        &self,
        features: &[Feature],
        areas: Option<HashMap<(PrecursorId, bool), (Peak, Vec<f64>), fnv::FnvBuildHasher>>,
        filenames: &[String],
    ) -> anyhow::Result<Url> {
        let path = self.make_path("results.sage.report.html");

        let global_q_value_filter = 0.01;
        let predict_section_q_value_filter = 0.01;

        // Create a new report
        let mut report = Report::new(
            "Sage",
            &self.parameters.version,
            Some("https://github.com/lazear/sage/blob/master/figures/logo.png?raw=true"),
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
                    for ((_id, decoy), (peak, data)) in areas {
                        if !decoy && peak.q_value <= global_q_value_filter {
                            intensities.push(data[i] as f32);
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
                    for ((_id, decoy), (peak, data)) in &areas {
                        if !decoy && peak.q_value <= global_q_value_filter {
                            intensities.push(data[i].log2());
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
    use super::passes_localization_filter;
    use sage_core::scoring::Feature;

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
}
