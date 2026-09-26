use super::*;
use sage_core::enzyme::DigestGroup;
use sage_core::sequence::PeptideSequence;
use sage_core::spectrum_index::{SpectrumIndex, SpectrumIndexBuilder, SpectrumIndexSettings};

/// Spectrum index size, as a fraction of `max_memory_gb`, at which a batch of
/// spectra is searched before more files are read.
const INDEX_MEMORY_FRACTION: f64 = 0.25;
const DEFAULT_INDEX_BUDGET_GIB: f64 = 8.0;

impl Runner {
    /// Retain every peptide that can contribute a preliminary fragment match
    /// to any spectrum. Spectra are indexed, and each sequence-coherent digest
    /// chunk is expanded and streamed through the spectrum index, so no
    /// fragment index is built until the final survivor database. When the
    /// spectra exceed the index budget, they are indexed in batches and the
    /// database is streamed once per batch.
    ///
    /// Processed spectra are also returned for the search, from the first
    /// file batch until they would exceed the index budget.
    pub(crate) fn prefilter_peptides(
        self,
        parallel: usize,
        fasta: Fasta,
        custom_cleavages: Option<ValidatedCustomCleavageLibrary>,
        unmodified_bytes: Option<u64>,
    ) -> anyhow::Result<(Vec<Peptide>, RetainedSpectra)> {
        let db_params = self.database_parameters.clone();
        let unmodified_bytes = unmodified_bytes.unwrap_or_else(|| {
            db_params.estimate_unmodified_bytes(&fasta, custom_cleavages.as_ref())
        });
        let budgets = PrefilterBudgets::for_search(&self.parameters);
        let passes = digest_passes(budgets.digest_bytes, unmodified_bytes);
        let mut digests = DigestSource {
            db_params: &db_params,
            fasta: &fasta,
            custom_cleavages: custom_cleavages.as_ref(),
            passes,
            chunks_per_pass: chunks_per_pass(
                fasta.targets.len(),
                db_params.prefilter_chunk_size,
                passes,
            ),
            cached: None,
        };
        if passes > 1 {
            info!(
                "streaming the digest in {} sequence buckets ({:.2} GiB unmodified digest estimated)",
                passes,
                unmodified_bytes as f64 / (1024.0 * 1024.0 * 1024.0),
            );
        }

        let budget = budgets.index_bytes;
        let mut pass = SurvivorPass {
            db_params: &db_params,
            keeps: Vec::new(),
            peptides: Vec::new(),
        };
        let mut builder = self.spectrum_index_builder(&db_params);
        let batches = self
            .parameters
            .mzml_paths
            .chunks(parallel.max(1))
            .collect::<Vec<_>>();
        let mut retained = RetainedSpectra {
            batch_size: parallel.max(1),
            batches: Vec::with_capacity(batches.len()),
        };
        let mut retained_bytes = 0u64;
        let mut retaining = true;
        for (batch_idx, batch) in batches.iter().enumerate() {
            // The search reports per-file progress when it takes these spectra
            // or reads them again.
            let spectra = self.read_processed_spectra_with_ms1(
                batch,
                batch_idx,
                parallel.max(1),
                self.requires_ms1(),
                false,
            )?;
            builder.add(&spectra.1);
            let bytes = spectra
                .0
                .iter()
                .chain(&spectra.1)
                .map(spectrum_bytes)
                .fold(0u64, u64::saturating_add);
            retaining &= retained_bytes.saturating_add(bytes) <= budget;
            if retaining {
                retained_bytes += bytes;
                retained.batches.push(Some(spectra));
            } else {
                retained.batches.push(None);
                drop(spectra);
            }
            let last = batch_idx + 1 == batches.len();
            if !last && (builder.allocated_bytes() as u64) < budget {
                continue;
            }
            let index = Self::finish_spectrum_index(builder);
            pass.run(&index, &mut digests, last);
            builder = self.spectrum_index_builder(&db_params);
        }

        let mut all_peptides = pass.peptides;
        Parameters::reorder_peptides(&mut all_peptides);
        let kept = retained.batches.iter().flatten().count();
        if kept > 0 {
            info!(
                "retained {} of {} file batches ({:.1} MiB of processed spectra) for the search",
                kept,
                batches.len(),
                retained_bytes as f64 / (1024.0 * 1024.0)
            );
        }
        Ok((all_peptides, retained))
    }

    fn spectrum_index_builder(&self, db_params: &Parameters) -> SpectrumIndexBuilder {
        let settings = SpectrumIndexSettings {
            // Precursors are corrected only under ppm tolerances.
            precursor_tol: if self.mass_recalibration_enabled()
                && matches!(self.parameters.precursor_tol, Tolerance::Ppm(_, _))
            {
                super::search::widen_for_recalibration(self.parameters.precursor_tol)
            } else {
                self.parameters.precursor_tol
            },
            fragment_tol: if self.mass_recalibration_enabled() {
                super::search::widen_for_recalibration(self.parameters.fragment_tol)
            } else {
                self.parameters.fragment_tol
            },
            min_isotope_err: self.parameters.isotope_errors.0,
            max_isotope_err: self.parameters.isotope_errors.1,
            min_precursor_charge: self.parameters.precursor_charge.0,
            max_precursor_charge: self.parameters.precursor_charge.1,
            override_precursor_charge: self.parameters.override_precursor_charge,
            max_fragment_charge: self.parameters.max_fragment_charge,
            wide_window: self.parameters.wide_window,
            min_peaks: self.parameters.min_peaks,
            // Never stricter than the search itself.
            min_matched_peaks: db_params
                .prefilter_min_matched_peaks
                .min(self.parameters.min_matched_peaks.max(1)),
            max_peaks: db_params.prefilter_max_peaks,
        };
        SpectrumIndexBuilder::new(settings, &db_params.mass_offset_modifications())
    }

    fn finish_spectrum_index(builder: SpectrumIndexBuilder) -> SpectrumIndex {
        let start = Instant::now();
        let index = builder.finish();
        info!(
            "indexed {} spectra: {} peak sets, {} peaks, {:.1} MiB, window depth {:.1}{} in {}ms",
            index.spectra(),
            index.probes(),
            index.peaks(),
            index.allocated_bytes() as f64 / (1024.0 * 1024.0),
            index.depth(),
            if index.uses_global_index() {
                ", global peak index"
            } else {
                ""
            },
            start.elapsed().as_millis(),
        );
        index
    }
}

/// Approximate heap and inline size of a processed spectrum.
fn spectrum_bytes(spectrum: &ProcessedSpectrum) -> u64 {
    (std::mem::size_of::<ProcessedSpectrum>()
        + spectrum.id.capacity()
        + spectrum.precursors.capacity() * std::mem::size_of::<sage_core::spectrum::Precursor>()
        + spectrum.masses.capacity() * std::mem::size_of::<f32>()
        + spectrum.intensities.capacity() * std::mem::size_of::<f32>()
        + spectrum.charges.capacity()
        + spectrum.charge_is_known.capacity()
        + spectrum.mobilities.capacity() * std::mem::size_of::<f32>()) as u64
}

/// Unmodified digest size, as a fraction of `max_memory_gb`, held by one
/// digest pass of the prefilter.
const DIGEST_MEMORY_FRACTION: f64 = 0.125;
const DEFAULT_DIGEST_BUDGET_GIB: f64 = 2.0;
/// Each pass digests the whole FASTA, so a tiny budget must not multiply
/// that work without bound.
const MAX_DIGEST_PASSES: u64 = 256;

/// Byte budgets the prefilter works within, derived from `max_memory_gb`.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) struct PrefilterBudgets {
    /// Unmodified digest held by one digest pass.
    pub digest_bytes: u64,
    /// Spectrum index built from one batch of files.
    pub index_bytes: u64,
}

impl PrefilterBudgets {
    pub(crate) fn from_max_memory(max_memory_gb: Option<f64>) -> Self {
        let limit = max_memory_gb.filter(|gib| *gib > 0.0);
        let bytes = |gib: f64| ((gib * 1024.0 * 1024.0 * 1024.0) as u64).max(1);
        Self {
            digest_bytes: bytes(limit.map_or(DEFAULT_DIGEST_BUDGET_GIB, |gib| {
                gib * DIGEST_MEMORY_FRACTION
            })),
            index_bytes: bytes(
                limit.map_or(DEFAULT_INDEX_BUDGET_GIB, |gib| gib * INDEX_MEMORY_FRACTION),
            ),
        }
    }

    /// Budgets for `parameters`, honoring a test override.
    pub(crate) fn for_search(parameters: &Search) -> Self {
        parameters
            .prefilter_budgets
            .unwrap_or_else(|| Self::from_max_memory(parameters.max_memory_gb))
    }
}

/// Number of sequence buckets the prefilter digests one at a time so the
/// unmodified digest of each stays within `budget_bytes`.
pub(crate) fn digest_passes(budget_bytes: u64, unmodified_bytes: u64) -> u64 {
    unmodified_bytes
        .div_ceil(budget_bytes.max(1))
        .clamp(1, MAX_DIGEST_PASSES)
}

/// Chunks each digest pass is split into, so all passes together give about
/// one chunk per `chunk_size` FASTA proteins.
pub(crate) fn chunks_per_pass(proteins: usize, chunk_size: usize, passes: u64) -> usize {
    proteins
        .div_ceil(chunk_size.max(1))
        .max(1)
        .div_ceil(passes.max(1) as usize)
        .max(1)
}

/// One sequence bucket of the digest, split into chunks, with the target
/// sequences its decoys are checked against. The set is shared so a cached
/// pass is not copied for every spectrum batch.
#[derive(Clone, Default)]
struct DigestPass {
    chunks: Vec<Vec<DigestGroup>>,
    target_sequences: Arc<HashSet<PeptideSequence>>,
}

/// Digest passes, regenerated for every spectrum batch. When the whole digest
/// is one pass it is kept between batches instead.
struct DigestSource<'a> {
    db_params: &'a Parameters,
    fasta: &'a Fasta,
    custom_cleavages: Option<&'a ValidatedCustomCleavageLibrary>,
    passes: u64,
    chunks_per_pass: usize,
    cached: Option<DigestPass>,
}

impl DigestSource<'_> {
    /// Digest bucket `pass`. Digestion and chunking are deterministic, so
    /// every spectrum batch sees the same chunks in the same order.
    fn pass(&mut self, pass: u64, last: bool) -> DigestPass {
        if let Some(cached) = &mut self.cached {
            return match last {
                true => std::mem::take(cached),
                false => cached.clone(),
            };
        }
        let start = Instant::now();
        let digests = self.db_params.digest_unmodified_bucket(
            self.fasta,
            self.custom_cleavages,
            pass,
            self.passes,
        );
        let target_sequences = Arc::new(
            digests
                .iter()
                .filter(|digest| !digest.reference.decoy)
                .map(|digest| digest.reference.sequence.clone())
                .collect::<HashSet<_>>(),
        );
        let chunk_size = digests.len().div_ceil(self.chunks_per_pass).max(1);
        let groups = digests.len();
        let chunks = Parameters::partition_digests_by_sequence(digests, chunk_size);
        info!(
            "digest pass {}: {} peptide groups in {} sequence-coherent chunks ({}ms)",
            pass,
            groups,
            chunks.len(),
            start.elapsed().as_millis(),
        );
        let digested = DigestPass {
            chunks,
            target_sequences,
        };
        if self.passes == 1 && !last {
            self.cached = Some(digested.clone());
        }
        digested
    }
}

/// Survivor state shared by every spectrum batch.
struct SurvivorPass<'a> {
    db_params: &'a Parameters,
    /// One survivor set per digest chunk. Chunk expansion is deterministic,
    /// so peptide positions agree between batches.
    keeps: Vec<AtomicBitSet>,
    peptides: Vec<Peptide>,
}

impl SurvivorPass<'_> {
    /// Stream every digest chunk through `index`. The final batch also closes
    /// label and decoy partners and collects the survivors.
    fn run(&mut self, index: &SpectrumIndex, digests: &mut DigestSource, last: bool) {
        let search_start = Instant::now();
        let mut streamed = 0usize;
        let mut retained = 0usize;
        let mut chunk_id = 0usize;
        for pass in 0..digests.passes {
            let DigestPass {
                chunks,
                target_sequences,
            } = digests.pass(pass, last);
            for digest_chunk in chunks {
                let (chunk_streamed, chunk_retained) =
                    self.run_chunk(index, digest_chunk, &target_sequences, chunk_id, last);
                streamed += chunk_streamed;
                retained += chunk_retained;
                chunk_id += 1;
            }
        }
        match last {
            true => info!(
                "- prefilter search:  {:8} ms ({} peptides streamed, {} retained)",
                search_start.elapsed().as_millis(),
                streamed,
                retained,
            ),
            false => info!(
                "- prefilter batch:   {:8} ms ({} peptides streamed)",
                search_start.elapsed().as_millis(),
                streamed,
            ),
        }
    }

    /// Stream one chunk; returns the peptides streamed and retained.
    fn run_chunk(
        &mut self,
        index: &SpectrumIndex,
        digest_chunk: Vec<DigestGroup>,
        target_sequences: &HashSet<PeptideSequence>,
        chunk_id: usize,
        last: bool,
    ) -> (usize, usize) {
        let start = Instant::now();
        let peptides = self
            .db_params
            .clone()
            .modify_digests_with_target_sequences(digest_chunk, target_sequences);
        let generated = Instant::now();
        if self.keeps.len() == chunk_id {
            self.keeps.push(AtomicBitSet::new(peptides.len()));
        }
        let keep = &self.keeps[chunk_id];
        assert_eq!(
            keep.len(),
            peptides.len(),
            "prefilter chunk expansion changed between spectrum batches"
        );
        index.filter(self.db_params, &peptides, keep);
        let filtered = Instant::now();
        if !last {
            info!(
                "prefilter chunk {}: streamed {} peptides (generate {}ms, stream {}ms)",
                chunk_id,
                peptides.len(),
                (generated - start).as_millis(),
                (filtered - generated).as_millis(),
            );
            return (peptides.len(), 0);
        }

        // Closure needs mass-ordered peptides and decoy pairing, but no
        // fragment index.
        let db = self.db_params.clone().build_peptide_table(peptides);
        LabelGroupIndex::new(&db.peptides).close(keep);
        close_prefilter_pairs(&db, keep);

        // Discarded peptides are released in parallel.
        let total = db.peptides.len();
        let peptides = db
            .peptides
            .into_par_iter()
            .enumerate()
            .filter_map(|(ix, peptide)| keep.contains(ix).then_some(peptide))
            .collect::<Vec<_>>();

        info!(
            "prefilter chunk {}: kept {} of {} peptides (generate {}ms, stream {}ms, closure {}ms)",
            chunk_id,
            peptides.len(),
            total,
            (generated - start).as_millis(),
            (filtered - generated).as_millis(),
            filtered.elapsed().as_millis(),
        );
        let kept = peptides.len();
        self.peptides.extend(peptides);
        (total, kept)
    }
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn digest_passes_fit_the_budget_and_are_capped() {
        const GIB: u64 = 1024 * 1024 * 1024;
        let passes = |max_memory_gb, bytes| {
            digest_passes(
                PrefilterBudgets::from_max_memory(max_memory_gb).digest_bytes,
                bytes,
            )
        };
        // An eighth of 16 GiB is 2 GiB per pass.
        assert_eq!(passes(Some(16.0), 0), 1);
        assert_eq!(passes(Some(16.0), 2 * GIB), 1);
        assert_eq!(passes(Some(16.0), 2 * GIB + 1), 2);
        assert_eq!(passes(Some(16.0), 38 * GIB), 19);
        // No or invalid limit falls back to the 2 GiB default.
        assert_eq!(passes(None, 5 * GIB), 3);
        assert_eq!(passes(Some(0.0), 5 * GIB), 3);
        assert_eq!(passes(Some(1e-9), u64::MAX), MAX_DIGEST_PASSES);
        // A quarter of 16 GiB indexes spectra; 8 GiB without a limit.
        assert_eq!(
            PrefilterBudgets::from_max_memory(Some(16.0)).index_bytes,
            4 * GIB
        );
        assert_eq!(PrefilterBudgets::from_max_memory(None).index_bytes, 8 * GIB);
    }

    #[test]
    fn passes_share_the_requested_chunks() {
        // 12 chunks of 100 proteins, over 1, 5, or more passes than chunks.
        assert_eq!(chunks_per_pass(1200, 100, 1), 12);
        assert_eq!(chunks_per_pass(1200, 100, 5), 3);
        assert_eq!(chunks_per_pass(1200, 100, 20), 1);
        assert_eq!(chunks_per_pass(0, 0, 1), 1);
    }
}
