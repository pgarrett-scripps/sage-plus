use super::*;

impl Runner {
    pub fn prefilter_peptides(
        self,
        parallel: usize,
        fasta: Fasta,
        custom_cleavages: Option<ValidatedCustomCleavageLibrary>,
    ) -> anyhow::Result<Vec<Peptide>> {
        let spectra: Option<Vec<ProcessedSpectrum>> =
            match parallel >= self.parameters.mzml_paths.len() {
                true => Some(
                    self.read_processed_spectra(&self.parameters.mzml_paths, 0, 0)?
                        .1,
                ),
                false => None,
            };

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
        let digest_chunks = Parameters::partition_digests_by_sequence(digests, digest_chunk_size);
        info!(
            "using {} sequence-coherent prefilter chunks",
            digest_chunks.len()
        );
        let mut all_peptides = Vec::new();
        for (chunk_id, digest_chunk) in digest_chunks.into_iter().enumerate() {
            let start = Instant::now();
            info!("pre-filtering fasta chunk {}", chunk_id,);
            let peptides = db_params
                .clone()
                .modify_digests_with_target_sequences(digest_chunk, &target_sequences);
            let mut db = db_params.clone().build_from_peptides(peptides);

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
            };

            let keep = AtomicBitSet::new(db.peptides.len());

            match &spectra {
                Some(spectra) => self.peptide_filter_processed_spectra(&scorer, spectra, &keep),
                None => {
                    for (chunk_idx, chunk) in
                        self.parameters.mzml_paths.chunks(parallel).enumerate()
                    {
                        let spectra_chunk =
                            self.read_processed_spectra(chunk, chunk_idx, parallel)?.1;
                        self.peptide_filter_processed_spectra(&scorer, &spectra_chunk, &keep);
                    }
                }
            };

            LabelGroupIndex::new(&db.peptides).close(&keep);
            close_prefilter_pairs(&db, &keep);

            // Retain only peptides where `keep[ix] = true`
            let peptides = db
                .peptides
                .drain(..)
                .enumerate()
                .filter_map(|(ix, peptide)| {
                    if keep.contains(ix) {
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
            all_peptides.extend(peptides);
        }

        Parameters::reorder_peptides(&mut all_peptides);
        Ok(all_peptides)
    }

    pub(super) fn peptide_filter_processed_spectra(
        &self,
        scorer: &Scorer,
        spectra: &[ProcessedSpectrum],
        keep: &AtomicBitSet,
    ) {
        use std::sync::atomic::{AtomicUsize, Ordering};
        let counter = AtomicUsize::new(0);
        let start = Instant::now();

        spectra
            .par_iter()
            .filter(|spec| spec.masses.len() >= self.parameters.min_peaks && spec.level == 2)
            .for_each(|spectrum| {
                let prev = counter.fetch_add(1, Ordering::Relaxed);
                if prev > 0 && prev.is_multiple_of(10_000) {
                    let duration = Instant::now().duration_since(start).as_millis() as usize;

                    let rate = prev * 1000 / (duration + 1);
                    log::trace!("- searched {} spectra ({} spectra/s)", prev, rate);
                }
                scorer.exact_prefilter(spectrum, keep)
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
}
