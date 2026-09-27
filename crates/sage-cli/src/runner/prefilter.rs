use super::*;
use sage_core::database::{sequence_hashes, DigestExpander};
use sage_core::enzyme::{
    group_protein_digests, Digest, DigestGroup, EnzymeParameters, ProteinOccurrence,
};
use sage_core::spectrum_index::{SpectrumIndex, SpectrumIndexBuilder, SpectrumIndexSettings};

/// Spectrum index size, as a fraction of `max_memory_gb`, at which a batch of
/// spectra is searched before more files are read.
const INDEX_MEMORY_FRACTION: f64 = 0.25;
const DEFAULT_INDEX_BUDGET_GIB: f64 = 8.0;

impl Runner {
    /// Retain every peptide that can contribute a preliminary fragment match
    /// to any spectrum. Spectra are indexed, then proteins are streamed in
    /// parallel: each protein is digested, modified, and checked against the
    /// spectrum index on its own, and only its survivors are kept. When the
    /// spectra exceed the index budget, they are indexed in batches and the
    /// proteins are streamed once per batch.
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
        let started = Instant::now();
        let enzyme: EnzymeParameters = db_params.enzyme.clone().into();
        let shared = SharedSequences::scan(
            &fasta,
            &enzyme,
            custom_cleavages.as_ref(),
            db_params.generate_decoys,
        );
        let scanned = started.elapsed();

        let budget = PrefilterBudgets::for_search(&self.parameters).index_bytes;
        let mut stream = ProteinStream {
            db_params: &db_params,
            fasta: &fasta,
            enzyme: &enzyme,
            custom_cleavages: custom_cleavages.as_ref(),
            shared,
            keeps: (0..fasta.targets.len())
                .map(|_| AtomicBitSet::new(0))
                .collect(),
            units: None,
            peptides: Vec::new(),
            checked: 0,
            deferred_digests: 0,
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
        let mut spectrum_batches = 0usize;
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
            stream.run(&index, last);
            spectrum_batches += 1;
            builder = self.spectrum_index_builder(&db_params);
        }

        let (checked, deferred_digests) = (stream.checked, stream.deferred_digests);
        let mut all_peptides = stream.peptides;
        db_params.reorder_peptides_with_labels(&mut all_peptides);
        info!(
            "- prefilter search:  {:8} ms ({} proteins streamed{}, {} peptides checked, {} kept; {} shared sequences handled together; sharing scan {}ms)",
            started.elapsed().as_millis(),
            fasta.targets.len(),
            match spectrum_batches {
                1 => String::new(),
                n => format!(" {n} times"),
            },
            checked,
            all_peptides.len(),
            deferred_digests,
            scanned.as_millis(),
        );
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

    /// Whole-digest prefilter the streamed one must reproduce: every
    /// spectrum in one index, the whole database digested and expanded at
    /// once, then filtered and closed.
    #[cfg(test)]
    pub(crate) fn prefilter_whole_digest(
        &self,
        parallel: usize,
        fasta: &Fasta,
        custom_cleavages: Option<&ValidatedCustomCleavageLibrary>,
    ) -> anyhow::Result<Vec<Peptide>> {
        let db_params = &self.database_parameters;
        let mut builder = self.spectrum_index_builder(db_params);
        for (batch_idx, batch) in self.parameters.mzml_paths.chunks(parallel).enumerate() {
            let spectra = self.read_processed_spectra_with_ms1(
                batch,
                batch_idx,
                parallel,
                self.requires_ms1(),
                false,
            )?;
            builder.add(&spectra.1);
        }
        let index = Self::finish_spectrum_index(builder);
        let groups = db_params.digest_unmodified_with_custom_cleavages(fasta, custom_cleavages);
        let db = db_params
            .clone()
            .build_peptide_table(db_params.modify_digests(groups));
        let keep = AtomicBitSet::new(db.peptides.len());
        index.filter(db_params, &db.peptides, &keep);
        LabelGroupIndex::new(&db.peptides).close(&keep);
        close_prefilter_pairs(&db, &keep);
        Ok(db
            .peptides
            .into_iter()
            .enumerate()
            .filter_map(|(ix, peptide)| keep.contains(ix).then_some(peptide))
            .collect())
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

/// Byte budget the prefilter works within, derived from `max_memory_gb`.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) struct PrefilterBudgets {
    /// Spectrum index built from one batch of files.
    pub index_bytes: u64,
}

impl PrefilterBudgets {
    pub(crate) fn from_max_memory(max_memory_gb: Option<f64>) -> Self {
        let limit = max_memory_gb.filter(|gib| *gib > 0.0);
        let gib = limit.map_or(DEFAULT_INDEX_BUDGET_GIB, |gib| gib * INDEX_MEMORY_FRACTION);
        Self {
            index_bytes: ((gib * 1024.0 * 1024.0 * 1024.0) as u64).max(1),
        }
    }

    /// Budget for `parameters`, honoring a test override.
    pub(crate) fn for_search(parameters: &Search) -> Self {
        parameters
            .prefilter_budgets
            .unwrap_or_else(|| Self::from_max_memory(parameters.max_memory_gb))
    }
}

/// Canonical keys of digest sequences that occur more than once in the
/// database, counting each generated decoy's sequence as an occurrence too.
///
/// A digest's key is the smaller hash of its sequence and, with generated
/// decoys, of its decoy sequence (`reversed_internal` is an involution). A
/// sequence, its decoy, and any target equal to that decoy share a key.
///
/// A digest whose key is unique yields peptides that no other digest shares:
/// no other protein merges into them, no decoy collides with them, and their
/// decoy and label partners come from the same digest. Such digests are
/// expanded and filtered one protein at a time. The rest are grouped by key;
/// each key is an independent small database, filtered exactly as the whole
/// database digest would filter it. A hash collision only moves a digest to
/// the shared set or merges two keys.
struct SharedSequences {
    keys: Vec<u64>,
    generate_decoys: bool,
}

impl SharedSequences {
    fn scan(
        fasta: &Fasta,
        enzyme: &EnzymeParameters,
        custom_cleavages: Option<&ValidatedCustomCleavageLibrary>,
        generate_decoys: bool,
    ) -> Self {
        let mut keys = (0..fasta.targets.len())
            .into_par_iter()
            .flat_map_iter(|index| {
                fasta
                    .digest_protein(index, enzyme, custom_cleavages)
                    .into_iter()
                    .map(|digest| {
                        let (forward, reversed) = sequence_hashes(digest.sequence.as_bytes());
                        match generate_decoys {
                            // A palindrome's decoy collides with itself.
                            true if forward == reversed => PALINDROME,
                            true => forward.min(reversed),
                            false => forward,
                        }
                    })
            })
            .collect::<Vec<_>>();
        keys.par_sort_unstable();
        let mut shared = keys
            .windows(2)
            .filter(|pair| pair[0] == pair[1])
            .map(|pair| pair[0])
            .collect::<Vec<_>>();
        drop(keys);
        shared.dedup();
        Self {
            keys: shared,
            generate_decoys,
        }
    }

    /// Canonical key of `digest`, when it is shared.
    fn shared_key(&self, digest: &Digest) -> Option<u64> {
        let (forward, reversed) = sequence_hashes(digest.sequence.as_bytes());
        let key = match self.generate_decoys {
            true if forward == reversed => return Some(PALINDROME),
            true => forward.min(reversed),
            false => forward,
        };
        self.keys.binary_search(&key).is_ok().then_some(key)
    }
}

/// Key under which palindromic digests are always shared.
const PALINDROME: u64 = u64::MAX;

/// Digests sharing one canonical key, with the peptides they kept so far.
struct SharedUnit {
    groups: Vec<DigestGroup>,
    keep: AtomicBitSet,
}

/// Survivor state shared by every spectrum batch.
struct ProteinStream<'a> {
    db_params: &'a Parameters,
    fasta: &'a Fasta,
    enzyme: &'a EnzymeParameters,
    custom_cleavages: Option<&'a ValidatedCustomCleavageLibrary>,
    shared: SharedSequences,
    /// Survivors of each protein's own peptides. Expansion is deterministic,
    /// so peptide positions agree between spectrum batches.
    keeps: Vec<AtomicBitSet>,
    /// Shared digests by key, collected while streaming the first batch.
    units: Option<Vec<SharedUnit>>,
    peptides: Vec<Peptide>,
    checked: usize,
    deferred_digests: usize,
}

impl ProteinStream<'_> {
    /// Stream every protein through `index`. The final batch also closes
    /// label and decoy partners and collects the survivors.
    fn run(&mut self, index: &SpectrumIndex, last: bool) {
        let collect_shared = self.units.is_none();
        let (fasta, enzyme, custom_cleavages, shared) =
            (self.fasta, self.enzyme, self.custom_cleavages, &self.shared);
        let generate_decoys = self.db_params.generate_decoys;
        let results = self.db_params.with_digest_expander(|expander| {
            self.keeps
                .par_iter_mut()
                .enumerate()
                .map_init(
                    || index.scratch(),
                    |scratch, (protein, keep)| {
                        let mut deferred = Vec::new();
                        let mut peptides = Vec::new();
                        let first = keep.is_empty();
                        for digest in fasta.digest_protein(protein, enzyme, custom_cleavages) {
                            if let Some(key) = shared.shared_key(&digest) {
                                if collect_shared {
                                    deferred.push((key, digest));
                                }
                                continue;
                            }
                            let group = DigestGroup {
                                origins: vec![ProteinOccurrence::of_protein_digest(&digest)],
                                reference: digest,
                            };
                            peptides.extend(expander.expand(group, None));
                        }
                        let checked = filter_unit(
                            expander,
                            index,
                            scratch,
                            &mut peptides,
                            keep,
                            last,
                            generate_decoys,
                        );
                        (peptides, deferred, if first { checked } else { 0 })
                    },
                )
                .collect::<Vec<_>>()
        });

        let mut deferred = Vec::new();
        for (peptides, protein_deferred, checked) in results {
            self.peptides.extend(peptides);
            deferred.extend(protein_deferred);
            self.checked += checked;
        }
        if collect_shared {
            self.deferred_digests = deferred.len();
            self.units = Some(shared_units(deferred));
        }
        self.run_shared(index, last);
    }

    /// Filter each shared key's digests as their own small database, as the
    /// whole-database digest would.
    fn run_shared(&mut self, index: &SpectrumIndex, last: bool) {
        let generate_decoys = self.db_params.generate_decoys;
        let units = self.units.as_mut().expect("shared digests are collected");
        let results = self.db_params.with_digest_expander(|expander| {
            units
                .par_iter_mut()
                .map_init(
                    || index.scratch(),
                    |scratch, unit| {
                        let first = unit.keep.is_empty();
                        let targets = unit
                            .groups
                            .iter()
                            .filter(|group| !group.reference.decoy)
                            .map(|group| group.reference.sequence.clone())
                            .collect::<HashSet<_>>();
                        let groups = match last {
                            true => std::mem::take(&mut unit.groups),
                            false => unit.groups.clone(),
                        };
                        let mut peptides = groups
                            .into_iter()
                            .flat_map(|group| expander.expand(group, Some(&targets)))
                            .collect::<Vec<_>>();
                        let checked = filter_unit(
                            expander,
                            index,
                            scratch,
                            &mut peptides,
                            &mut unit.keep,
                            last,
                            generate_decoys,
                        );
                        (peptides, if first { checked } else { 0 })
                    },
                )
                .collect::<Vec<_>>()
        });
        for (peptides, checked) in results {
            self.peptides.extend(peptides);
            self.checked += checked;
        }
        if last {
            self.units = Some(Vec::new());
        }
    }
}

/// Group shared digests into one unit per canonical key.
fn shared_units(mut deferred: Vec<(u64, Digest)>) -> Vec<SharedUnit> {
    deferred.par_sort_unstable_by_key(|(key, _)| *key);
    let mut units = Vec::new();
    let mut digests = deferred.into_iter().peekable();
    while let Some((key, digest)) = digests.next() {
        let mut unit = vec![digest];
        while let Some((_, digest)) = digests.next_if(|(next, _)| *next == key) {
            unit.push(digest);
        }
        units.push(unit);
    }
    units
        .into_par_iter()
        .map(|digests| SharedUnit {
            groups: group_protein_digests(digests),
            keep: AtomicBitSet::new(0),
        })
        .collect()
}

/// Mark the peptides of one protein or shared unit that match `index`,
/// adding to `keep` from earlier batches. On the last batch, close label
/// groups and decoy pairs and leave only the survivors in `peptides`;
/// otherwise leave it empty. Returns the number of peptides checked.
fn filter_unit(
    expander: &DigestExpander,
    index: &SpectrumIndex,
    scratch: &mut sage_core::spectrum_index::Scratch,
    peptides: &mut Vec<Peptide>,
    keep: &mut AtomicBitSet,
    last: bool,
    generate_decoys: bool,
) -> usize {
    expander.reorder(peptides);
    let checked = peptides.len();
    if keep.is_empty() {
        *keep = AtomicBitSet::new(peptides.len());
    }
    assert_eq!(
        keep.len(),
        peptides.len(),
        "prefilter expansion changed between spectrum batches"
    );
    for (ix, peptide) in peptides.iter().enumerate() {
        if !keep.contains(ix) && index.matches(expander.parameters(), peptide, scratch) {
            keep.insert(ix);
        }
    }
    if !last {
        // Drop the buffer too; only the keep bits outlive this batch.
        *peptides = Vec::new();
        return checked;
    }
    LabelGroupIndex::new(peptides).close(keep);
    if generate_decoys {
        // A partner's partner is the peptide itself, so one pass closes pairs.
        let partners = (0..peptides.len())
            .filter(|&ix| keep.contains(ix))
            .filter_map(|ix| sorted_pair(peptides, ix))
            .collect::<Vec<_>>();
        for ix in partners {
            keep.insert(ix);
        }
    }
    let mut ix = 0;
    peptides.retain(|_| {
        ix += 1;
        keep.contains(ix - 1)
    });
    // Survivors are held until the whole database is streamed.
    peptides.shrink_to_fit();
    // Release the survivor set; it is not needed after the last batch.
    *keep = AtomicBitSet::new(0);
    checked
}

/// Index of the generated target or decoy partner of `peptides[index]` in
/// mass-sorted `peptides`, as [`IndexedDatabase::paired_peptide_index`] finds
/// it.
fn sorted_pair(peptides: &[Peptide], index: usize) -> Option<usize> {
    let peptide = &peptides[index];
    let paired = peptide.reverse();
    let found = peptides
        .binary_search_by(|candidate| {
            candidate
                .monoisotopic
                .total_cmp(&paired.monoisotopic)
                .then_with(|| candidate.initial_sort(&paired))
        })
        .ok()?;
    (peptides[found].decoy != peptide.decoy).then_some(found)
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn index_budget_follows_max_memory() {
        const GIB: u64 = 1024 * 1024 * 1024;
        // A quarter of 16 GiB indexes spectra; 8 GiB without a limit.
        assert_eq!(
            PrefilterBudgets::from_max_memory(Some(16.0)).index_bytes,
            4 * GIB
        );
        assert_eq!(PrefilterBudgets::from_max_memory(None).index_bytes, 8 * GIB);
        assert_eq!(
            PrefilterBudgets::from_max_memory(Some(0.0)).index_bytes,
            8 * GIB
        );
    }
}
