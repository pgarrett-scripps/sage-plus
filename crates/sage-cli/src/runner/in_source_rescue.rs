//! The in-source rescue pass: per file, search the MS2 spectra the main search
//! left unidentified for in-source products of the file's confident peptides,
//! each only where its parent elutes. See [`sage_core::in_source_rescue`].

use super::*;
use fnv::{FnvHashMap, FnvHashSet};
use sage_core::database::PeptideIx;
use sage_core::enzyme::EnzymeParameters;
use sage_core::in_source::{
    file_retention, InSourceAnnotation, InSourceAnnotations, InSourceKind, InSourceSettings,
    Placement, PEPTIDE_Q_VALUE, PSM_Q_VALUE,
};
use sage_core::in_source_rescue::{
    candidates_of, competition_q_values, parent_windows, precursor_reach_problem,
    InSourceRescueFailure, InSourceRescueFileStats, Window, RESCUE_Q_VALUE, RESCUE_REPORT_PSMS,
};
use sage_core::ml::retention_alignment::Alignment;

/// Where a candidate came from: every confident parent that gives it.
#[derive(Clone, Copy, Debug, PartialEq)]
pub(super) struct Origin {
    pub(super) parent: PeptideIx,
    pub(super) kind: InSourceKind,
    pub(super) parent_len: usize,
}

/// One rescue hit before it is numbered and registered: the PSM, the parent
/// it is a product of, and its peptide in the file's rescue database.
pub(super) type RescueHit = (Feature, Origin, sage_core::peptide::Peptide);

/// What the rescue pass produced for the whole run. It is all or nothing: when
/// any file failed, `features` is empty for every file, nothing was added to
/// the database or the annotations, and `error` names the file.
#[derive(Debug, Default)]
pub(super) struct RescueOutcome {
    pub(super) features: Vec<Feature>,
    /// Each rescued PSM's occurrence among its file's spectra of that ID,
    /// where it is not the first, keyed by `psm_id` as the search keys its own.
    pub(super) occurrences: HashMap<usize, usize>,
    pub(super) stats: Vec<InSourceRescueFileStats>,
    pub(super) error: Option<InSourceRescueFailure>,
}

#[cfg(test)]
thread_local! {
    /// Test hook: PSMs added to the search's own after its FDR, as if it had
    /// found them. Compiled into test builds only.
    pub(super) static INJECT_PSMS: std::cell::RefCell<Vec<Feature>> =
        const { std::cell::RefCell::new(Vec::new()) };
    /// Test hook: make the pass fail on this file id, as any per-file error
    /// would, so a test can drive a whole run through a rescue failure.
    /// Compiled into test builds only.
    pub(super) static FAIL_REREAD: std::cell::Cell<Option<usize>> =
        const { std::cell::Cell::new(None) };
}

/// The null-control shift of every co-elution window, in widths. Only a
/// binary built with the `in-source-rescue-null-control` feature reads
/// [`sage_core::in_source_rescue::NULL_SHIFT_ENV`], so no released build can
/// have its results changed by the environment; a value that is set but is
/// not a finite number stops the run instead of being ignored.
#[cfg(feature = "in-source-rescue-null-control")]
pub(super) fn null_shift() -> anyhow::Result<Option<f32>> {
    use sage_core::in_source_rescue::NULL_SHIFT_ENV;
    match std::env::var(NULL_SHIFT_ENV) {
        Err(std::env::VarError::NotPresent) => Ok(None),
        Err(error) => anyhow::bail!("{NULL_SHIFT_ENV}: {error}"),
        Ok(value) => match value.trim().parse::<f32>() {
            Ok(widths) if widths.is_finite() => Ok(Some(widths)),
            _ => anyhow::bail!("{NULL_SHIFT_ENV} must be a finite number of widths, not `{value}`"),
        },
    }
}

#[cfg(not(feature = "in-source-rescue-null-control"))]
pub(super) fn null_shift() -> anyhow::Result<Option<f32>> {
    Ok(None)
}

impl Runner {
    /// Run the rescue pass over every file and return its PSMs, numbered after
    /// every PSM in `features`, with each one's parent recorded in
    /// `annotations`. The rescued peptides are added to the database beside the
    /// materialized mass-offset peptides, so the writers resolve them like any
    /// other. Nothing in `features` changes: rescued PSMs never enter the main
    /// FDR, protein inference or quantification, which have all run by now.
    /// All or nothing: when any file fails, no rescue rows are returned for
    /// any file and `error` names the file.
    pub(super) fn in_source_rescue(
        &mut self,
        features: &[Feature],
        settings: &InSourceSettings,
        annotations: &mut InSourceAnnotations,
        alignments: Option<&[Alignment]>,
        occurrences: &HashMap<usize, usize>,
    ) -> anyhow::Result<RescueOutcome> {
        let files = self.parameters.mzml_paths.len();
        // The configuration is refused up front for these (`Input::build`);
        // a runner built some other way fails the same way here.
        let refusal = if self.parameters.wide_window || !self.parameters.dia.is_off() {
            Some("it searches DDA spectra only, not wide-window or DIA searches".to_string())
        } else if self.parameters.override_precursor_charge {
            Some(
                "it needs each scan's own precursor charge, not `override_precursor_charge`"
                    .to_string(),
            )
        } else {
            precursor_reach_problem(
                self.parameters.precursor_tol,
                self.parameters.isotope_errors,
                self.database_parameters.peptide_max_mass,
            )
        };
        if let Some(reason) = refusal {
            return Ok(RescueOutcome {
                error: Some(InSourceRescueFailure {
                    message: reason,
                    ..Default::default()
                }),
                ..Default::default()
            });
        }
        let null_shift = null_shift()?;
        if let Some(widths) = null_shift {
            log::warn!(
                "in_source.rescue null control: every co-elution window is shifted by {widths} widths"
            );
        }
        // Every file is searched before anything is kept, so a failure on one
        // leaves no rescue rows from the others.
        let mut per_file = Vec::with_capacity(files);
        let mut error = None;
        for file_id in 0..files {
            self.cancellation.check()?;
            let started = Instant::now();
            let (hits, file_stats) = self.rescue_file(
                file_id,
                features,
                settings,
                annotations,
                occurrences,
                null_shift,
            );
            match &file_stats.error {
                Some(reason) => {
                    log::error!(
                        "- in-source rescue failed on file {file_id} (`{}`): {reason}",
                        self.parameters.mzml_paths[file_id]
                    );
                    error.get_or_insert_with(|| InSourceRescueFailure {
                        file_id: Some(file_id),
                        path: Some(self.parameters.mzml_paths[file_id].to_string()),
                        message: reason.clone(),
                    });
                }
                None => log::info!(
                    "- in-source rescue, file {file_id}: {} of {} unidentified MS2 spectra rescued at {} q-value ({} targets, {} decoys competed; {} candidates from {} parents) in {} ms",
                    file_stats.rescued_psms,
                    file_stats.searched_spectra,
                    RESCUE_Q_VALUE,
                    file_stats.target_hits,
                    file_stats.decoy_hits,
                    file_stats.candidates,
                    file_stats.parents,
                    started.elapsed().as_millis()
                ),
            }
            per_file.push((hits, file_stats));
        }
        if error.is_some() {
            log::error!("in_source.rescue wrote no rescue rows for any file, because it did not complete on every file");
            // Nothing of the files that finished is kept, not even counts.
            return Ok(RescueOutcome {
                error,
                ..Default::default()
            });
        }
        let stats = per_file.iter().map(|(_, stats)| stats.clone()).collect();
        let mut next_psm_id = features.iter().map(|f| f.psm_id + 1).max().unwrap_or(0);
        let mut registered: HashMap<(String, bool), PeptideIx> = HashMap::new();
        // A product rescued in several files is one peptide, in the proteins
        // of every parent that gave it in any file.
        let mut merged: HashMap<(String, bool), sage_core::peptide::Peptide> = HashMap::new();
        for ((_, _, peptide), _) in per_file.iter().flat_map(|(hits, _)| hits) {
            match merged.entry((peptide.to_string(), peptide.decoy)) {
                std::collections::hash_map::Entry::Occupied(mut entry) => {
                    merge_proteins(entry.get_mut(), peptide)
                }
                std::collections::hash_map::Entry::Vacant(entry) => {
                    entry.insert(peptide.clone());
                }
            }
        }
        let mut rescued = Vec::new();
        let mut rescued_occurrences = HashMap::new();
        for (hits, _) in per_file {
            for ((mut feature, origin, peptide), occurrence) in hits {
                let key = (peptide.to_string(), peptide.decoy);
                feature.peptide_idx = *registered.entry(key.clone()).or_insert_with(|| {
                    self.database
                        .push_unindexed_peptide(merged.remove(&key).unwrap_or(peptide))
                });
                feature.psm_id = next_psm_id;
                next_psm_id += 1;
                if occurrence > 0 {
                    rescued_occurrences.insert(feature.psm_id, occurrence);
                }
                // Retention times are written on the run's aligned scale, as
                // every other row's are, when the run was aligned.
                feature.aligned_rt = alignments
                    .and_then(|alignments| alignments.get(feature.file_id))
                    .map_or(feature.rt, |alignment| alignment.transform(feature.rt));
                annotations.insert_rescued(
                    feature.psm_id,
                    InSourceAnnotation {
                        parent: origin.parent,
                        kind: origin.kind,
                    },
                );
                rescued.push(feature);
            }
        }
        Ok(RescueOutcome {
            features: rescued,
            occurrences: rescued_occurrences,
            stats,
            error: None,
        })
    }

    /// The rescue hits of one file and what the pass did there. A file whose
    /// spectra cannot be read again is returned with `error` set and no hits.
    pub(super) fn rescue_file(
        &self,
        file_id: usize,
        features: &[Feature],
        settings: &InSourceSettings,
        annotations: &InSourceAnnotations,
        occurrences: &HashMap<usize, usize>,
        null_shift: Option<f32>,
    ) -> (Vec<(RescueHit, usize)>, InSourceRescueFileStats) {
        let mut stats = InSourceRescueFileStats {
            file_id,
            null_shift_widths: null_shift,
            ..Default::default()
        };
        #[cfg(test)]
        if FAIL_REREAD.with(|fail| fail.get()) == Some(file_id) {
            stats.error = Some("injected by the test harness".into());
            return (Vec::new(), stats);
        }
        // Spectra with a confident PSM are identified; every other MS2
        // spectrum is searched. A spectrum is its ID and its occurrence among
        // the file's spectra of that ID, as the search numbers them, so a
        // repeated ID is told apart. Parents are the confident target PSMs, by
        // the thresholds the in-source annotation uses.
        let occurrence_of =
            |feature: &Feature| occurrences.get(&feature.psm_id).copied().unwrap_or(0);
        let mut identified = FnvHashSet::default();
        let mut rts: FnvHashMap<PeptideIx, Vec<f32>> = FnvHashMap::default();
        for feature in features.iter().filter(|f| f.file_id == file_id) {
            if feature.spectrum_q <= PSM_Q_VALUE && feature.label == 1 {
                identified.insert((feature.spec_id.as_str(), occurrence_of(feature)));
            }
            if feature.label == 1
                && feature.spectrum_q <= PSM_Q_VALUE
                && feature.peptide_q <= PEPTIDE_Q_VALUE
                && !self.database[feature.peptide_idx].decoy
            {
                rts.entry(feature.peptide_idx).or_default().push(feature.rt);
            }
        }
        if rts.is_empty() {
            return (Vec::new(), stats);
        }
        // The medians and tolerance the annotation derives for this file.
        let (medians, rt_tolerance, _) = file_retention(&mut rts, settings);
        let windows = parent_windows(
            file_id,
            &medians,
            rt_tolerance,
            annotations.apexes(),
            null_shift,
        );
        stats.parents = medians.len();

        let (origins, peptides) = self.rescue_candidates(&medians, &mut stats);
        if peptides.is_empty() {
            return (Vec::new(), stats);
        }
        let database = self.rescue_database(peptides);
        // Each rescue peptide's parents, by its index in the rescue database.
        let origins = database
            .peptides
            .iter()
            .map(|peptide| origins.get(&(peptide.to_string(), peptide.decoy)))
            .collect::<Vec<_>>();

        // The unidentified MS2 spectra, read and processed as the main search
        // reads and processes them, searchable by the same rule.
        let url = &self.parameters.mzml_paths[file_id];
        let (sn, processor) = self.spectrum_reading();
        let reread = sage_cloudpath::util::read_spectra(
            url,
            file_id,
            sn,
            self.parameters.bruker_config,
            false,
        );
        let raw = match reread {
            Ok(raw) => raw,
            Err(error) => {
                stats.error = Some(format!("could not read `{url}` again: {error}"));
                return (Vec::new(), stats);
            }
        };
        // Occurrences are numbered over the MS2-and-higher spectra in file
        // order, as the search numbers them.
        let mut seen: HashMap<String, usize> = HashMap::new();
        let numbered = raw
            .into_iter()
            .filter(|spectrum| spectrum.ms_level != 1)
            .map(|spectrum| {
                let count = seen.entry(spectrum.id.clone()).or_insert(0);
                *count += 1;
                (spectrum, *count - 1)
            })
            .collect::<Vec<_>>();
        let spectra = numbered
            .into_par_iter()
            .filter(|(spectrum, occurrence)| {
                spectrum.ms_level == 2 && !identified.contains(&(spectrum.id.as_str(), *occurrence))
            })
            .map(|(spectrum, occurrence)| (processor.process(spectrum), occurrence))
            .filter(|(spectrum, _)| spectrum.is_searchable(self.parameters.min_peaks))
            .collect::<Vec<_>>();
        // A scan without a precursor charge would be searched at every
        // charge, which reads a parent's own spectrum at a fraction of its
        // mass, where its truncations lie; decoys cannot model that.
        let searched = spectra.len();
        let spectra = spectra
            .into_iter()
            .filter(|(spectrum, _)| {
                spectrum
                    .precursors
                    .first()
                    .is_some_and(|precursor| precursor.charge.is_some())
            })
            .collect::<Vec<_>>();
        stats.skipped_without_charge = searched - spectra.len();
        stats.searched_spectra = spectra.len();

        // Each spectrum is searched against only the candidates whose parent
        // co-elutes with it, dropped before any candidate is cut, so a
        // candidate of a parent eluting elsewhere can never take its place.
        let mut scorer = self.scorer();
        scorer.db = &database;
        scorer.report_psms = scorer.report_psms.max(RESCUE_REPORT_PSMS);
        scorer.chimera = false;
        let hits = spectra
            .par_iter()
            .filter_map(|(spectrum, occurrence)| {
                let rt = spectrum.scan_start_time;
                let keep = |peptide: PeptideIx| {
                    origins[peptide.0 as usize]
                        .is_some_and(|origins| coeluting_origin(origins, &windows, rt).is_some())
                };
                let mut scored = scorer.score_restricted(spectrum, &keep).into_iter();
                let mut feature = scored.next()?;
                let origin = coeluting_origin(
                    origins[feature.peptide_idx.0 as usize]?,
                    &windows,
                    feature.rt,
                )?;
                // Every reported candidate co-elutes, so delta_next and
                // delta_best are among the candidates this spectrum was
                // searched against, as for any search.
                debug_assert!(feature.delta_best == 0.0);
                let peptide = database[feature.peptide_idx].clone();
                if let Some(immonium) = &self.parameters.immonium {
                    feature.immonium = Some(immonium.evaluate(spectrum, &peptide));
                }
                feature.rank = 1;
                Some(((feature, origin, peptide), *occurrence))
            })
            .collect::<Vec<_>>();
        let (hits, hit_occurrences): (Vec<_>, Vec<_>) = hits.into_iter().unzip();
        let hits = finish_hits(hits, &mut stats);
        (hits.into_iter().zip(hit_occurrences).collect(), stats)
    }

    /// The database a file's rescue candidates are searched in: the main
    /// search's settings without its mass offsets, because a rescue candidate
    /// is a product as it is and an offset hypothesis would only take scorer
    /// slots from it.
    pub(super) fn rescue_database(
        &self,
        mut peptides: Vec<sage_core::peptide::Peptide>,
    ) -> sage_core::database::IndexedDatabase {
        sage_core::database::Parameters::reorder_peptides(&mut peptides);
        let mut database = self
            .database_parameters
            .clone()
            .build_from_peptides(peptides);
        database.mass_offsets.clear();
        database.offset_library = None;
        database
    }

    /// The candidates of a file's parents, deduplicated across parents, without
    /// those the database already holds, each with its decoy. A pair whose
    /// decoy is its own target (a palindromic interior), another candidate's
    /// target, or a database peptide is dropped, as the database digest drops
    /// such decoys; otherwise the pair would merge into one target and lose its
    /// decoy. Returns every peptide's parents, keyed by peptide and decoy flag.
    #[allow(clippy::type_complexity)]
    pub(super) fn rescue_candidates(
        &self,
        parents: &FnvHashMap<PeptideIx, f32>,
        stats: &mut InSourceRescueFileStats,
    ) -> (
        HashMap<(String, bool), Vec<Origin>>,
        Vec<sage_core::peptide::Peptide>,
    ) {
        let min_len = self.database_parameters.enzyme.min_len.unwrap_or(5);
        let enzyme = EnzymeParameters::from(self.database_parameters.enzyme.clone()).enzyme;
        let mut parents = parents.keys().copied().collect::<Vec<_>>();
        parents.sort_unstable();
        let mut targets: HashMap<String, (sage_core::peptide::Peptide, Vec<Origin>)> =
            HashMap::new();
        let mut order = Vec::new();
        for parent in parents {
            for candidate in candidates_of(&self.database[parent], parent, min_len, enzyme.as_ref())
            {
                if self
                    .database
                    .indexed_peptidoform(&candidate.peptide)
                    .is_some()
                {
                    continue;
                }
                let origin = Origin {
                    parent: candidate.parent,
                    kind: candidate.kind,
                    parent_len: candidate.parent_len,
                };
                let key = candidate.peptide.to_string();
                match targets.entry(key.clone()) {
                    std::collections::hash_map::Entry::Occupied(mut entry) => {
                        // A product of several parents lies in all of their
                        // proteins.
                        merge_proteins(&mut entry.get_mut().0, &candidate.peptide);
                        entry.get_mut().1.push(origin);
                    }
                    std::collections::hash_map::Entry::Vacant(entry) => {
                        order.push(key);
                        entry.insert((candidate.peptide, vec![origin]));
                    }
                }
            }
        }
        let target_sequences = targets
            .values()
            .map(|(peptide, _)| {
                self.database_parameters
                    .decoy_collision_key(&peptide.sequence)
            })
            .collect::<HashSet<_>>();
        let mut origins = HashMap::new();
        let mut peptides = Vec::new();
        for key in order {
            let (target, parents) = targets.remove(&key).expect("every key was inserted");
            let decoy = target.reverse();
            let mut as_target = decoy.clone();
            as_target.decoy = false;
            let collides = decoy.to_string() == key
                || target_sequences.contains(
                    &self
                        .database_parameters
                        .decoy_collision_key(&decoy.sequence),
                )
                || self.database.indexed_peptidoform(&as_target).is_some();
            if collides {
                stats.dropped_decoy_collisions += 1;
                continue;
            }
            stats.candidates += 1;
            origins.insert((decoy.to_string(), true), parents.clone());
            origins.insert((key, false), parents);
            peptides.push(target);
            peptides.push(decoy);
        }
        (origins, peptides)
    }
}

/// Add `other`'s proteins and protein sites to `peptide`'s, as the database
/// digest merges a peptide found in several proteins.
pub(super) fn merge_proteins(
    peptide: &mut sage_core::peptide::Peptide,
    other: &sage_core::peptide::Peptide,
) {
    peptide.proteins.extend(other.proteins.iter().cloned());
    peptide.proteins.sort_unstable();
    peptide.proteins.dedup();
    let mut sites = peptide.protein_sites.to_vec();
    sites.extend(other.protein_sites.iter().cloned());
    sites.sort_unstable();
    sites.dedup();
    peptide.protein_sites = sites.into();
}

/// The longest of a peptide's parents that elutes at `rt`, the lowest index
/// breaking a tie, or none.
pub(super) fn coeluting_origin(
    origins: &[Origin],
    windows: &FnvHashMap<PeptideIx, Vec<Window>>,
    rt: f32,
) -> Option<Origin> {
    origins
        .iter()
        .filter(|origin| {
            windows
                .get(&origin.parent)
                .is_some_and(|windows| windows.iter().any(|window| window.contains(rt)))
        })
        .max_by(|a, b| {
            a.parent_len
                .cmp(&b.parent_len)
                .then_with(|| b.parent.cmp(&a.parent))
        })
        .copied()
}

/// Give the kept hits of one file their q-values from target-decoy
/// competition among themselves, count them, and set the fields a rescue row
/// carries.
pub(super) fn finish_hits(
    hits: Vec<RescueHit>,
    stats: &mut InSourceRescueFileStats,
) -> Vec<RescueHit> {
    let evidence = hits
        .iter()
        .map(|(feature, _, peptide)| (feature.hyperscore, peptide.decoy))
        .collect::<Vec<_>>();
    let q_values = competition_q_values(&evidence);
    let mut out = Vec::with_capacity(hits.len());
    for ((mut feature, origin, peptide), q) in hits.into_iter().zip(q_values) {
        if peptide.decoy {
            stats.decoy_hits += 1;
        } else {
            stats.target_hits += 1;
            if q <= RESCUE_Q_VALUE {
                stats.rescued_psms += 1;
                match (origin.kind.placement, origin.kind.loss) {
                    (_, Some(_)) => stats.rescued_loss += 1,
                    (Placement::Prefix, None) => stats.rescued_prefix += 1,
                    _ => stats.rescued_suffix += 1,
                }
            }
        }
        // A rescued PSM carries the q-value of the rescue competition and no
        // peptide- or protein-level confidence: it is a product, not a sample
        // peptide, and it was never in those calculations.
        feature.rank = 1;
        feature.label = if peptide.decoy { -1 } else { 1 };
        feature.spectrum_q = q;
        feature.peptide_q = 1.0;
        feature.protein_q = 1.0;
        feature.protein_group_q = 1.0;
        feature.discriminant_score = 0.0;
        feature.posterior_error = 1.0;
        feature.missed_cleavages = peptide.missed_cleavages;
        out.push((feature, origin, peptide));
    }
    out
}

#[cfg(test)]
#[path = "../../tests/unit/runner/in_source_rescue.rs"]
mod test;
