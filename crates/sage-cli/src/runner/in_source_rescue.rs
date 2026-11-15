//! The in-source rescue pass: per file, search the MS2 spectra the main search
//! left unidentified for in-source products of the file's confident peptides,
//! each only where its parent elutes. See [`sage_core::in_source_rescue`].

use super::*;
use fnv::{FnvHashMap, FnvHashSet};
use sage_core::database::PeptideIx;
use sage_core::in_source::{
    derive_rt_tolerance, pairwise_half_distances, InSourceAnnotation, InSourceAnnotations,
    InSourceKind, InSourceSettings, PEPTIDE_Q_VALUE, PSM_Q_VALUE,
};
use sage_core::in_source_rescue::{
    candidates_of, competition_q_values, parent_windows, InSourceRescueFileStats, Window,
    NULL_SHIFT_ENV, RESCUE_Q_VALUE, RESCUE_REPORT_PSMS,
};

/// Where a candidate came from: every confident parent that gives it.
#[derive(Clone, Copy)]
struct Origin {
    parent: PeptideIx,
    kind: InSourceKind,
    parent_len: usize,
}

impl Runner {
    /// Run the rescue pass over every file and return its PSMs, numbered after
    /// `next_psm_id`, with each one's parent recorded in `annotations`. The
    /// rescued peptides are added to the database beside the materialized
    /// mass-offset peptides, so the writers resolve them like any other.
    /// Nothing in `features` changes: rescued PSMs never enter the main FDR,
    /// protein inference or quantification, which have all run by now.
    pub(super) fn in_source_rescue(
        &mut self,
        features: &[Feature],
        settings: &InSourceSettings,
        annotations: &mut InSourceAnnotations,
    ) -> anyhow::Result<(Vec<Feature>, Vec<InSourceRescueFileStats>)> {
        if self.parameters.wide_window || !self.parameters.dia.is_off() {
            log::warn!(
                "in_source.rescue searches DDA spectra only; it is skipped for wide-window and DIA searches"
            );
            return Ok((Vec::new(), Vec::new()));
        }
        let null_shift = std::env::var(NULL_SHIFT_ENV)
            .ok()
            .and_then(|value| value.parse::<f32>().ok());
        if let Some(widths) = null_shift {
            log::warn!("in_source.rescue null control: every co-elution window is shifted by {widths} widths");
        }
        let min_len = self.database_parameters.enzyme.min_len.unwrap_or(5);
        let mut next_psm_id = features.iter().map(|f| f.psm_id + 1).max().unwrap_or(0);
        let mut registered: HashMap<(String, bool), PeptideIx> = HashMap::new();
        let mut rescued = Vec::new();
        let mut stats = Vec::new();
        for file_id in 0..self.parameters.mzml_paths.len() {
            self.cancellation.check()?;
            let started = Instant::now();
            let (mut file_features, file_stats) = self.rescue_file(
                file_id,
                features,
                settings,
                annotations,
                min_len,
                null_shift,
            )?;
            log::info!(
                "- in-source rescue, file {file_id}: {} of {} unidentified MS2 spectra rescued at {} q-value ({} targets, {} decoys competed; {} candidates from {} parents) in {} ms",
                file_stats.rescued_psms,
                file_stats.searched_spectra,
                RESCUE_Q_VALUE,
                file_stats.target_hits,
                file_stats.decoy_hits,
                file_stats.candidates,
                file_stats.parents,
                started.elapsed().as_millis()
            );
            for (feature, origin, peptide) in file_features.drain(..) {
                let key = (peptide.to_string(), peptide.decoy);
                let index = *registered
                    .entry(key)
                    .or_insert_with(|| self.database.push_unindexed_peptide(peptide));
                let mut feature = feature;
                feature.peptide_idx = index;
                feature.psm_id = next_psm_id;
                next_psm_id += 1;
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

    #[allow(clippy::type_complexity)]
    fn rescue_file(
        &self,
        file_id: usize,
        features: &[Feature],
        settings: &InSourceSettings,
        annotations: &InSourceAnnotations,
        min_len: usize,
        null_shift: Option<f32>,
    ) -> anyhow::Result<(
        Vec<(Feature, Origin, sage_core::peptide::Peptide)>,
        InSourceRescueFileStats,
    )> {
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
            return Ok((Vec::new(), stats));
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

        // Candidates, deduplicated across parents, without those the main
        // search already scored, each with its decoy.
        let mut origins: HashMap<(String, bool), Vec<Origin>> = HashMap::new();
        let mut peptides = Vec::new();
        let mut parents = medians.keys().copied().collect::<Vec<_>>();
        parents.sort_unstable();
        for parent in parents {
            for candidate in candidates_of(&self.database[parent], parent, min_len) {
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
                let decoy = candidate.peptide.reverse();
                for peptide in [candidate.peptide, decoy] {
                    let key = (peptide.to_string(), peptide.decoy);
                    let entry = origins.entry(key).or_default();
                    if entry.is_empty() {
                        peptides.push(peptide);
                    }
                    entry.push(origin);
                }
            }
        }
        stats.candidates = peptides.iter().filter(|p| !p.decoy).count();
        if peptides.is_empty() {
            return Ok((Vec::new(), stats));
        }
        sage_core::database::Parameters::reorder_peptides(&mut peptides);
        let database = self
            .database_parameters
            .clone()
            .build_from_peptides(peptides);

        // The unidentified MS2 spectra, processed as the main search processes
        // them.
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
        let spectra = sage_cloudpath::util::read_spectra(
            url,
            file_id,
            sn,
            self.parameters.bruker_config,
            false,
        )
        .with_context(|| format!("in_source.rescue: failed to reread `{url}`"))?
        .into_par_iter()
        .filter(|spectrum| spectrum.ms_level == 2 && !identified.contains(&spectrum.id))
        .map(|spectrum| processor.process(spectrum))
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
                scorer
                    .score(spectrum)
                    .into_iter()
                    .filter(|feature| feature.mass_offset.is_none())
                    .find_map(|feature| {
                        let peptide = &database[feature.peptide_idx];
                        let key = (peptide.to_string(), peptide.decoy);
                        origins
                            .get(&key)?
                            .iter()
                            .filter(|origin| {
                                windows.get(&origin.parent).is_some_and(|w| {
                                    w.iter().any(|w: &Window| w.contains(feature.rt))
                                })
                            })
                            .max_by(|a, b| {
                                a.parent_len
                                    .cmp(&b.parent_len)
                                    .then_with(|| b.parent.cmp(&a.parent))
                            })
                            .copied()
                            .map(|origin| (feature, origin, peptide.clone()))
                    })
            })
            .collect::<Vec<_>>();

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
                        (sage_core::in_source::Placement::Prefix, None) => {
                            stats.rescued_prefix += 1
                        }
                        _ => stats.rescued_suffix += 1,
                    }
                }
            }
            // A rescued PSM carries the q-value of the rescue competition and
            // no peptide- or protein-level confidence: it is a product, not a
            // sample peptide, and it was never in those calculations.
            feature.rank = 1;
            feature.label = if peptide.decoy { -1 } else { 1 };
            feature.spectrum_q = q;
            feature.peptide_q = 1.0;
            feature.protein_q = 1.0;
            feature.protein_group_q = 1.0;
            feature.discriminant_score = 0.0;
            feature.posterior_error = 1.0;
            feature.aligned_rt = feature.rt;
            out.push((feature, origin, peptide));
        }
        Ok((out, stats))
    }
}
