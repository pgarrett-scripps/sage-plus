use super::*;
use sage_core::enzyme::DigestGroup;
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
    ) -> anyhow::Result<(Vec<Peptide>, RetainedSpectra)> {
        let db_params = self.database_parameters.clone();
        let digests =
            db_params.digest_unmodified_with_custom_cleavages(&fasta, custom_cleavages.as_ref());
        let target_sequences = digests
            .iter()
            .filter(|digest| !digest.reference.decoy)
            .map(|digest| digest.reference.sequence.clone())
            .collect::<HashSet<_>>();
        let requested_chunks = fasta
            .targets
            .len()
            .div_ceil(db_params.prefilter_chunk_size.max(1))
            .max(1);
        let digest_chunk_size = digests.len().div_ceil(requested_chunks).max(1);
        let mut digest_chunks =
            Parameters::partition_digests_by_sequence(digests, digest_chunk_size);
        info!(
            "using {} sequence-coherent prefilter chunks",
            digest_chunks.len()
        );

        let budget = self.spectrum_index_budget();
        let mut pass = SurvivorPass {
            db_params: &db_params,
            target_sequences: &target_sequences,
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
            match last {
                true => pass.run(&index, std::mem::take(&mut digest_chunks), true),
                false => pass.run(&index, digest_chunks.clone(), false),
            }
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

    fn spectrum_index_budget(&self) -> u64 {
        let gib = std::env::var("SAGE_PREFILTER_INDEX_GB")
            .ok()
            .and_then(|value| value.parse::<f64>().ok())
            .or_else(|| {
                self.parameters
                    .max_memory_gb
                    .filter(|gib| *gib > 0.0)
                    .map(|gib| gib * INDEX_MEMORY_FRACTION)
            })
            .unwrap_or(DEFAULT_INDEX_BUDGET_GIB);
        (gib * 1024.0 * 1024.0 * 1024.0) as u64
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

/// Survivor state shared by every spectrum batch.
struct SurvivorPass<'a> {
    db_params: &'a Parameters,
    target_sequences: &'a HashSet<sage_core::sequence::PeptideSequence>,
    /// One survivor set per digest chunk. Chunk expansion is deterministic,
    /// so peptide positions agree between batches.
    keeps: Vec<AtomicBitSet>,
    peptides: Vec<Peptide>,
}

impl SurvivorPass<'_> {
    /// Stream every digest chunk through `index`. The final batch also closes
    /// label and decoy partners and collects the survivors.
    fn run(&mut self, index: &SpectrumIndex, digest_chunks: Vec<Vec<DigestGroup>>, last: bool) {
        let search_start = Instant::now();
        let mut streamed = 0usize;
        let mut retained = 0usize;
        for (chunk_id, digest_chunk) in digest_chunks.into_iter().enumerate() {
            let start = Instant::now();
            let peptides = self
                .db_params
                .clone()
                .modify_digests_with_target_sequences(digest_chunk, self.target_sequences);
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
            streamed += peptides.len();
            let filtered = Instant::now();
            if !last {
                info!(
                    "prefilter chunk {}: streamed {} peptides (generate {}ms, stream {}ms)",
                    chunk_id,
                    peptides.len(),
                    (generated - start).as_millis(),
                    (filtered - generated).as_millis(),
                );
                continue;
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
            retained += peptides.len();
            self.peptides.extend(peptides);
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
}
