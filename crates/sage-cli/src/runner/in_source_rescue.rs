//! The in-source rescue pass: per file, search the MS2 spectra the main search
//! left unidentified for in-source products of the file's confident peptides,
//! each only where its parent elutes. See [`sage_core::in_source_rescue`].

use super::*;
use fnv::{FnvHashMap, FnvHashSet};
use sage_core::database::PeptideIx;
use sage_core::enzyme::EnzymeParameters;
use sage_core::in_source::{
    derive_rt_tolerance, pairwise_half_distances, InSourceAnnotation, InSourceAnnotations,
    InSourceKind, InSourceSettings, Placement, PEPTIDE_Q_VALUE, PSM_Q_VALUE,
};
use sage_core::in_source_rescue::{
    candidates_of, competition_q_values, parent_windows, precursor_reach_problem,
    InSourceRescueFileStats, Window, RESCUE_Q_VALUE, RESCUE_REPORT_PSMS,
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

/// The null-control shift of every co-elution window, in widths. Only a
/// binary built with the `in-source-rescue-null-control` feature reads
/// [`sage_core::in_source_rescue::NULL_SHIFT_ENV`], so no released build can
/// have its results changed by the environment; a value that is set but is
/// not a finite number stops the run instead of being ignored.
#[cfg(feature = "in-source-rescue-null-control")]
fn null_shift() -> anyhow::Result<Option<f32>> {
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
fn null_shift() -> anyhow::Result<Option<f32>> {
    Ok(None)
}

impl Runner {
    /// Run the rescue pass over every file and return its PSMs, numbered after
    /// every PSM in `features`, with each one's parent recorded in
    /// `annotations`. The rescued peptides are added to the database beside the
    /// materialized mass-offset peptides, so the writers resolve them like any
    /// other. Nothing in `features` changes: rescued PSMs never enter the main
    /// FDR, protein inference or quantification, which have all run by now.
    /// A file that cannot be searched is recorded and skipped; the run goes on.
    pub(super) fn in_source_rescue(
        &mut self,
        features: &[Feature],
        settings: &InSourceSettings,
        annotations: &mut InSourceAnnotations,
        alignments: Option<&[Alignment]>,
    ) -> anyhow::Result<(Vec<Feature>, Vec<InSourceRescueFileStats>)> {
        let files = self.parameters.mzml_paths.len();
        let skip_all = |reason: String| {
            log::warn!("in_source.rescue skipped: {reason}");
            (0..files)
                .map(|file_id| InSourceRescueFileStats {
                    file_id,
                    skipped: Some(reason.clone()),
                    ..Default::default()
                })
                .collect::<Vec<_>>()
        };
        if self.parameters.wide_window || !self.parameters.dia.is_off() {
            let reason =
                "it searches DDA spectra only, not wide-window or DIA searches".to_string();
            return Ok((Vec::new(), skip_all(reason)));
        }
        if let Some(problem) = precursor_reach_problem(
            self.parameters.precursor_tol,
            self.parameters.isotope_errors,
            self.database_parameters.peptide_max_mass,
        ) {
            return Ok((Vec::new(), skip_all(problem)));
        }
        let null_shift = null_shift()?;
        if let Some(widths) = null_shift {
            log::warn!(
                "in_source.rescue null control: every co-elution window is shifted by {widths} widths"
            );
        }
        let mut next_psm_id = features.iter().map(|f| f.psm_id + 1).max().unwrap_or(0);
        let mut registered: HashMap<(String, bool), PeptideIx> = HashMap::new();
        let mut rescued = Vec::new();
        let mut stats = Vec::new();
        for file_id in 0..files {
            self.cancellation.check()?;
            let started = Instant::now();
            let (hits, file_stats) =
                self.rescue_file(file_id, features, settings, annotations, null_shift);
            match &file_stats.skipped {
                Some(reason) => log::warn!("- in-source rescue, file {file_id} skipped: {reason}"),
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
            for (mut feature, origin, peptide) in hits {
                let key = (peptide.to_string(), peptide.decoy);
                feature.peptide_idx = *registered
                    .entry(key)
                    .or_insert_with(|| self.database.push_unindexed_peptide(peptide));
                feature.psm_id = next_psm_id;
                next_psm_id += 1;
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
            stats.push(file_stats);
        }
        Ok((rescued, stats))
    }

    /// The rescue hits of one file and what the pass did there. A file whose
    /// spectra cannot be read again is returned with `skipped` set and no hits.
    pub(super) fn rescue_file(
        &self,
        file_id: usize,
        features: &[Feature],
        settings: &InSourceSettings,
        annotations: &InSourceAnnotations,
        null_shift: Option<f32>,
    ) -> (Vec<RescueHit>, InSourceRescueFileStats) {
        let mut stats = InSourceRescueFileStats {
            file_id,
            null_shift_widths: null_shift,
            ..Default::default()
        };
        // Spectra with a confident PSM are identified; every other MS2
        // spectrum is searched. Parents are the confident target PSMs, by the
        // thresholds the in-source annotation uses.
        let mut identified = FnvHashSet::default();
        let mut rts: FnvHashMap<PeptideIx, Vec<f32>> = FnvHashMap::default();
        for feature in features.iter().filter(|f| f.file_id == file_id) {
            if feature.spectrum_q <= PSM_Q_VALUE && feature.label == 1 {
                identified.insert(feature.spec_id.clone());
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
        // The retention-time tolerance the annotation derives for this file.
        let mut spreads = Vec::new();
        let (mut rt_min, mut rt_max) = (f32::MAX, f32::MIN);
        let mut medians = FnvHashMap::default();
        for (&peptide, values) in rts.iter_mut() {
            values.sort_by(f32::total_cmp);
            rt_min = rt_min.min(values[0]);
            rt_max = rt_max.max(values[values.len() - 1]);
            pairwise_half_distances(values, &mut spreads);
            let mid = values.len() / 2;
            let median = if values.len() % 2 == 0 {
                (values[mid - 1] + values[mid]) / 2.0
            } else {
                values[mid]
            };
            medians.insert(peptide, median);
        }
        let rt_tolerance = match settings.rt_tolerance {
            Some(configured) => configured,
            None => derive_rt_tolerance(&mut spreads, (rt_max - rt_min).max(0.0)).0,
        };
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
        let mut peptides = peptides;
        sage_core::database::Parameters::reorder_peptides(&mut peptides);
        let database = self
            .database_parameters
            .clone()
            .build_from_peptides(peptides);

        // The unidentified MS2 spectra, read and processed as the main search
        // reads and processes them, searchable by the same rule.
        let url = &self.parameters.mzml_paths[file_id];
        let sn = self
            .parameters
            .quant
            .tmt_settings
            .sn
            .then_some(self.parameters.quant.tmt_settings.level);
        let min_deisotope_mz = match &self.parameters.quant.tmt {
            Some(i) if self.parameters.quant.tmt_settings.level == 2 => {
                i.reporter_masses().last().map(|x| x * (1.0 + 20E-6))
            }
            _ => None,
        };
        let processor = SpectrumProcessor::with_deisotope_settings(
            self.parameters.max_peaks,
            self.parameters.deisotope,
            min_deisotope_mz.unwrap_or(0.0),
        );
        let raw = match sage_cloudpath::util::read_spectra(
            url,
            file_id,
            sn,
            self.parameters.bruker_config,
            false,
        ) {
            Ok(raw) => raw,
            Err(error) => {
                stats.skipped = Some(format!("could not read `{url}` again: {error}"));
                return (Vec::new(), stats);
            }
        };
        let spectra = raw
            .into_par_iter()
            .filter(|spectrum| spectrum.ms_level == 2 && !identified.contains(&spectrum.id))
            .map(|spectrum| processor.process(spectrum))
            .filter(|spectrum| spectrum.is_searchable(self.parameters.min_peaks))
            .collect::<Vec<_>>();
        stats.searched_spectra = spectra.len();

        let mut scorer = self.scorer();
        scorer.db = &database;
        scorer.report_psms = scorer.report_psms.max(RESCUE_REPORT_PSMS);
        scorer.chimera = false;
        // Each spectrum keeps its best hit among the candidates whose parent
        // co-elutes with it, and reports the longest such parent.
        let hits = spectra
            .par_iter()
            .filter_map(|spectrum| {
                let hit = scorer
                    .score(spectrum)
                    .into_iter()
                    .filter(|feature| feature.mass_offset.is_none())
                    .find_map(|feature| {
                        let peptide = &database[feature.peptide_idx];
                        let key = (peptide.to_string(), peptide.decoy);
                        coeluting_origin(origins.get(&key)?, &windows, feature.rt)
                            .map(|origin| (feature, origin, peptide.clone()))
                    });
                hit.map(|(mut feature, origin, peptide)| {
                    if let Some(immonium) = &self.parameters.immonium {
                        feature.immonium = Some(immonium.evaluate(spectrum, &peptide));
                    }
                    (feature, origin, peptide)
                })
            })
            .collect::<Vec<_>>();
        (finish_hits(hits, &mut stats), stats)
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
                let entry = targets.entry(key.clone()).or_insert_with(|| {
                    order.push(key);
                    (candidate.peptide, Vec::new())
                });
                entry.1.push(origin);
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
