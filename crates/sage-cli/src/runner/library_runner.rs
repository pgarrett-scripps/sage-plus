use super::*;
use sage_core::database::PeptideIx;
use sage_core::enzyme::Position;
use sage_core::mass_calibration::CalibrationModel;
use sage_core::ml::retention_alignment::{fit_reference_alignment, Alignment};
use sage_core::spectral_library_search::{
    generate_decoys, parse_proforma, DdaLibraryIndex, DdaLibrarySearchParameters,
    LibrarySearchSettings,
};

fn source_basename(source: &str) -> Option<String> {
    source
        .rsplit(['/', '\\'])
        .next()
        .filter(|source| !source.is_empty())
        .map(str::to_ascii_lowercase)
}

pub(super) struct LibrarySearchRuntime {
    index: DdaLibraryIndex,
    peptide_indices: Vec<PeptideIx>,
    pub(super) target_entries: usize,
    pub(super) target_transitions: usize,
    max_retention_time_minutes: f32,
    source_files: HashSet<String>,
}

impl LibrarySearchRuntime {
    pub(super) fn overlapping_source_files(&self, query_paths: &[Url]) -> Vec<String> {
        let mut overlaps = query_paths
            .iter()
            .filter_map(sage_cloudpath::filename)
            .filter(|filename| self.source_files.contains(&filename.to_ascii_lowercase()))
            .map(str::to_owned)
            .collect::<Vec<_>>();
        overlaps.sort_unstable();
        overlaps.dedup();
        overlaps
    }
}

#[derive(Clone, Copy, Debug, Default)]
struct LibraryMassCalibration {
    precursor: Option<CalibrationModel>,
    fragment: Option<CalibrationModel>,
}

fn assign_entity_q_values(
    features: &mut [Feature],
    key: impl Fn(&Feature) -> Option<String>,
    assign: impl Fn(&mut Feature, f32),
) -> usize {
    let mut best = HashMap::<String, [f32; 2]>::new();
    for feature in features.iter() {
        if let Some(key) = key(feature) {
            let scores = best.entry(key).or_insert([f32::MIN; 2]);
            let decoy = usize::from(feature.label == -1);
            scores[decoy] = scores[decoy].max(feature.discriminant_score);
        }
    }
    let mut rows = best
        .into_iter()
        .map(|(key, scores)| {
            let decoy = scores[1] >= scores[0];
            (key, decoy, scores[usize::from(decoy)], 1.0f32)
        })
        .collect::<Vec<_>>();
    rows.sort_unstable_by(|left, right| right.2.total_cmp(&left.2));
    let mut targets = 0usize;
    let mut decoys = 1usize;
    for (_, decoy, _, q) in &mut rows {
        if *decoy {
            decoys += 1;
        } else {
            targets += 1;
        }
        *q = (decoys as f32 / targets as f32).min(1.0);
    }
    let mut minimum = 1.0f32;
    for row in rows.iter_mut().rev() {
        minimum = minimum.min(row.3);
        row.3 = minimum;
    }
    let values = rows
        .into_iter()
        .map(|(key, decoy, _, q)| ((key, decoy), q))
        .collect::<HashMap<_, _>>();
    let mut passing = HashSet::new();
    for feature in features {
        let Some(key) = key(feature) else {
            continue;
        };
        let q = values
            .get(&(key.clone(), feature.label == -1))
            .copied()
            .unwrap_or(1.0);
        assign(feature, q);
        if feature.label == 1 && q <= 0.01 {
            passing.insert(key);
        }
    }
    passing.len()
}

fn library_peptide(
    entry: &sage_core::spectral_library_search::DdaLibraryEntry,
) -> anyhow::Result<Peptide> {
    let parsed = parse_proforma(&entry.proforma).map_err(anyhow::Error::msg)?;
    let proteins = entry
        .proteins
        .split(';')
        .map(str::trim)
        .filter(|protein| !protein.is_empty())
        .map(Arc::<str>::from)
        .collect::<Vec<_>>();
    anyhow::ensure!(
        !proteins.is_empty(),
        "library entry `{}` has no protein accessions",
        entry.library_entry_id
    );
    Ok(Peptide {
        decoy: entry.is_decoy,
        sequence: Arc::from(parsed.sequence.into_boxed_slice()),
        modifications: parsed.modifications,
        applied_modifications: Arc::default(),
        label_channel: entry.label_channel.as_deref().map(Arc::from),
        label_group_override: entry.label_group.as_deref().map(Arc::from),
        nterm: parsed.nterm,
        cterm: parsed.cterm,
        monoisotopic: entry.precursor_neutral_mass,
        missed_cleavages: 0,
        semi_enzymatic: false,
        position: Position::Full,
        proteins,
        protein_sites: Arc::default(),
    })
}

pub(super) fn load_library_search(
    settings: &LibrarySearchSettings,
) -> anyhow::Result<(IndexedDatabase, LibrarySearchRuntime)> {
    let bytes = sage_cloudpath::util::read_bytes(&settings.path)
        .with_context(|| format!("Failed to read spectral library `{}`", settings.path))?;
    let targets = if settings.path.to_ascii_lowercase().ends_with(".parquet") {
        sage_cloudpath::parquet::deserialize_spectral_library(bytes).map_err(anyhow::Error::from)?
    } else {
        let text = String::from_utf8(bytes)
            .with_context(|| format!("spectral library `{}` is not UTF-8", settings.path))?;
        sage_core::spectral_library_search::deserialize_mzspeclib(&text)
            .map_err(anyhow::Error::msg)?
    };
    anyhow::ensure!(!targets.is_empty(), "spectral library contains no entries");
    let label_reference = targets
        .iter()
        .find_map(|entry| entry.label_reference.as_deref())
        .map(Arc::from);
    let mut label_channels = Vec::<Arc<str>>::new();
    for channel in targets
        .iter()
        .filter_map(|entry| entry.label_channel.as_deref())
    {
        if !label_channels
            .iter()
            .any(|configured| configured.as_ref() == channel)
        {
            label_channels.push(Arc::from(channel));
        }
    }
    if !label_channels.is_empty() {
        anyhow::ensure!(
            targets.iter().all(|entry| {
                if entry.label_channel.is_some() {
                    entry.label_group.is_some()
                        && entry.label_reference.as_deref() == label_reference.as_deref()
                } else {
                    entry.label_group.is_none() && entry.label_reference.is_none()
                }
            }),
            "spectral library contains inconsistent label metadata"
        );
    }
    let target_entries = targets.len();
    let target_transitions = targets.iter().map(|entry| entry.fragments.len()).sum();
    let source_files = targets
        .iter()
        .filter_map(|entry| source_basename(&entry.source_file))
        .collect::<HashSet<_>>();
    let max_retention_time_minutes = targets
        .iter()
        .map(|entry| entry.retention_time_minutes)
        .filter(|rt| rt.is_finite() && *rt > 0.0)
        .fold(0.0, f32::max);
    let entries = generate_decoys(targets, settings).map_err(anyhow::Error::msg)?;
    let index = DdaLibraryIndex::new(entries).map_err(anyhow::Error::msg)?;

    let mut peptides = Vec::<Peptide>::new();
    let mut peptide_indices = Vec::with_capacity(index.entries().len());
    let mut peptidoforms = HashMap::<(String, bool), PeptideIx>::new();
    for entry in index.entries() {
        // Targets with the same peptidoform (for example, multiple charge
        // states) share identity. Keep decoys pair-specific even if two
        // shuffles happen to produce the same sequence.
        let identity = if entry.is_decoy {
            entry.library_entry_id.clone()
        } else {
            entry.proforma.clone()
        };
        let key = (identity, entry.is_decoy);
        let peptide_idx = match peptidoforms.get(&key).copied() {
            Some(peptide_idx) => peptide_idx,
            None => {
                let peptide_idx = PeptideIx(peptides.len() as u32);
                peptides.push(library_peptide(entry)?);
                peptidoforms.insert(key, peptide_idx);
                peptide_idx
            }
        };
        peptide_indices.push(peptide_idx);
    }

    let target_by_id = index
        .entries()
        .iter()
        .zip(peptide_indices.iter().copied())
        .filter(|(entry, _)| !entry.is_decoy)
        .map(|(entry, peptide_idx)| (entry.library_entry_id.as_str(), peptide_idx))
        .collect::<HashMap<_, _>>();
    let mut decoy_pairing = (0..peptides.len())
        .map(|index| PeptideIx(index as u32))
        .collect::<Vec<_>>();
    for (entry, peptide_idx) in index.entries().iter().zip(peptide_indices.iter().copied()) {
        if entry.is_decoy {
            let target_id = entry
                .library_entry_id
                .strip_prefix(&settings.decoy_tag)
                .ok_or_else(|| anyhow::anyhow!("invalid generated decoy identifier"))?;
            let target_idx = target_by_id.get(target_id).copied().ok_or_else(|| {
                anyhow::anyhow!("generated decoy `{}` has no target", entry.library_entry_id)
            })?;
            decoy_pairing[peptide_idx.0 as usize] = target_idx;
        }
    }
    drop(target_by_id);
    let database = IndexedDatabase {
        peptides,
        generate_decoys: true,
        decoy_tag: settings.decoy_tag.clone(),
        decoy_pairing,
        label_reference,
        label_channels,
        ..IndexedDatabase::default()
    };
    Ok((
        database,
        LibrarySearchRuntime {
            index,
            peptide_indices,
            target_entries,
            target_transitions,
            max_retention_time_minutes,
            source_files,
        },
    ))
}

impl Runner {
    fn search_library_spectrum(
        &self,
        spectrum: &ProcessedSpectrum,
        calibration: LibraryMassCalibration,
    ) -> Vec<Feature> {
        let Some(runtime) = self.library_search.as_ref() else {
            return Vec::new();
        };
        let Some(precursor) = spectrum.precursors.first() else {
            return Vec::new();
        };
        let charges = match precursor.charge {
            Some(charge) if !self.parameters.override_precursor_charge => charge..=charge,
            _ => self.parameters.precursor_charge.0..=self.parameters.precursor_charge.1,
        };
        let search_parameters = DdaLibrarySearchParameters {
            precursor_tolerance: self.parameters.precursor_tol,
            fragment_tolerance: self.parameters.fragment_tol,
            min_matched_peaks: usize::from(self.parameters.min_matched_peaks),
            max_hits: self.parameters.report_psms,
            min_isotope_error: self.parameters.isotope_errors.0,
            max_isotope_error: self.parameters.isotope_errors.1,
            annotate_matches: self.parameters.annotate_matches,
            precursor_offset_ppm: calibration
                .precursor
                .map(|model| model.predict_ppm(spectrum.scan_start_time))
                .unwrap_or_default(),
            fragment_offset_ppm: calibration
                .fragment
                .map(|model| model.predict_ppm(spectrum.scan_start_time))
                .unwrap_or_default(),
        };
        let mut matches = charges
            .flat_map(|charge| {
                let neutral_mass = (precursor.mz - sage_core::mass::PROTON) * f32::from(charge);
                runtime
                    .index
                    .search(spectrum, neutral_mass, charge, search_parameters)
            })
            .collect::<Vec<_>>();
        matches.sort_unstable_by(|left, right| {
            right
                .spectral_angle
                .total_cmp(&left.spectral_angle)
                .then_with(|| right.matched_peaks.cmp(&left.matched_peaks))
                .then_with(|| {
                    left.precursor_ppm
                        .abs()
                        .total_cmp(&right.precursor_ppm.abs())
                })
                .then_with(|| {
                    left.isotope_error
                        .unsigned_abs()
                        .cmp(&right.isotope_error.unsigned_abs())
                })
                .then_with(|| left.entry_index.cmp(&right.entry_index))
        });
        let mut seen = HashSet::new();
        matches.retain(|matched| seen.insert(matched.entry_index));
        matches.truncate(self.parameters.report_psms);
        let best = matches
            .first()
            .map(|matched| matched.spectral_angle)
            .unwrap_or_default();
        let scored_candidates = matches.len() as u32;
        matches
            .iter()
            .enumerate()
            .map(|(rank, matched)| {
                let entry = &runtime.index.entries()[matched.entry_index];
                let peptide_idx = runtime.peptide_indices[matched.entry_index];
                let peptide = &self.database[peptide_idx];
                let next = matches
                    .get(rank + 1)
                    .map(|next| next.spectral_angle)
                    .unwrap_or_default();
                let ion_mobility = precursor.inverse_ion_mobility.unwrap_or_default();
                Feature {
                    peptide_idx,
                    psm_id: 0,
                    peptide_len: peptide.sequence.len(),
                    spec_id: spectrum.id.clone(),
                    file_id: spectrum.file_id,
                    rank: rank as u32 + 1,
                    label: peptide.label(),
                    expmass: (precursor.mz - sage_core::mass::PROTON)
                        * f32::from(entry.precursor_charge),
                    calcmass: entry.precursor_neutral_mass,
                    charge: entry.precursor_charge,
                    rt: spectrum.scan_start_time,
                    aligned_rt: spectrum.scan_start_time,
                    predicted_rt: if runtime.max_retention_time_minutes > 0.0 {
                        entry.retention_time_minutes / runtime.max_retention_time_minutes
                    } else {
                        0.0
                    },
                    ims: ion_mobility,
                    predicted_ims: entry.ion_mobility,
                    delta_mass: matched.raw_precursor_ppm,
                    aligned_delta_mass: matched.precursor_ppm,
                    hyperscore: f64::from(matched.spectral_angle),
                    delta_next: f64::from(matched.spectral_angle - next),
                    delta_best: f64::from(best - matched.spectral_angle),
                    matched_peaks: matched.matched_peaks as u32,
                    matched_intensity_pct: matched.explained_query_intensity * 100.0,
                    spectral_angle: matched.spectral_angle,
                    explained_library_intensity: matched.explained_library_intensity,
                    explained_query_intensity: matched.explained_query_intensity,
                    average_ppm: matched.average_fragment_ppm,
                    signed_fragment_ppm: matched.signed_fragment_ppm,
                    aligned_average_ppm: matched.aligned_average_fragment_ppm,
                    isotope_error: f32::from(matched.isotope_error) * sage_core::mass::NEUTRON,
                    scored_candidates,
                    poisson: -f64::from(matched.spectral_angle),
                    discriminant_score: matched.spectral_angle,
                    posterior_error: 1.0,
                    spectrum_q: 1.0,
                    peptide_q: 1.0,
                    protein_q: 1.0,
                    protein_group_q: 1.0,
                    ms2_intensity: matched.explained_query_intensity * spectrum.total_ion_current,
                    ambiguity_sequence: entry.proforma.clone(),
                    delta_rt_model: matched.retention_time_delta_minutes,
                    delta_ims_model: matched.ion_mobility_delta.unwrap_or_default(),
                    fragments: matched.fragments.clone(),
                    ..Feature::default()
                }
            })
            .collect()
    }

    fn batch_library_files(
        &self,
        batch_size: usize,
        calibrations: &[LibraryMassCalibration],
        collect_quantification: bool,
    ) -> anyhow::Result<SageResults> {
        let results = self
            .parameters
            .mzml_paths
            .chunks(batch_size)
            .enumerate()
            .map(|(chunk_idx, chunk)| {
                let spectra = self.read_processed_spectra_with_ms1(
                    chunk,
                    chunk_idx,
                    batch_size,
                    collect_quantification && self.requires_ms1(),
                )?;
                let features = spectra
                    .1
                    .par_iter()
                    .filter(|spectrum| {
                        spectrum.level == 2
                            && spectrum.masses.len() >= self.parameters.min_peaks
                            && !self.cancellation.is_cancelled()
                    })
                    .flat_map(|spectrum| {
                        self.search_library_spectrum(
                            spectrum,
                            calibrations
                                .get(spectrum.file_id)
                                .copied()
                                .unwrap_or_default(),
                        )
                    })
                    .collect();
                self.events.emit(EventKind::SearchProgress {
                    files_completed: (chunk_idx * batch_size + chunk.len())
                        .min(self.parameters.mzml_paths.len()),
                    files_total: self.parameters.mzml_paths.len(),
                });
                let result = if collect_quantification {
                    self.complete_features(spectra.1, spectra.0, features)
                } else {
                    SageResults {
                        features,
                        ..SageResults::default()
                    }
                };
                Ok(result)
            })
            .collect::<anyhow::Result<Vec<_>>>()?;
        Ok(results.into_iter().collect())
    }

    fn fit_library_mass_calibrations(&self, features: &[Feature]) -> Vec<LibraryMassCalibration> {
        let fit_options = FitOptions {
            min_linear_improvement: 0.0,
            ..FitOptions::default()
        };
        (0..self.parameters.mzml_paths.len())
            .map(|file_id| {
                let calibration_psms = features
                    .iter()
                    .filter(|feature| {
                        feature.file_id == file_id
                            && feature.rank == 1
                            && feature.label == 1
                            && feature.spectrum_q <= 0.01
                    })
                    .collect::<Vec<_>>();
                let precursor_points = calibration_psms
                    .iter()
                    .map(|feature| CalibrationPoint {
                        rt_minutes: feature.rt,
                        error_ppm: feature.delta_mass,
                    })
                    .collect::<Vec<_>>();
                let fragment_points = calibration_psms
                    .iter()
                    .map(|feature| CalibrationPoint {
                        rt_minutes: feature.rt,
                        error_ppm: feature.signed_fragment_ppm,
                    })
                    .collect::<Vec<_>>();
                let precursor = matches!(self.parameters.precursor_tol, Tolerance::Ppm(_, _))
                    .then(|| fit_mass_calibration(&precursor_points, fit_options))
                    .flatten()
                    .map(|fit| fit.model);
                let fragment = fit_mass_calibration(&fragment_points, fit_options)
                    .map(|fit| fit.model);
                if let Some(model) = precursor {
                    log::info!(
                        "- file {} library precursor calibration: {:?}, offset={:.3} ppm, slope={:.4} ppm/min, n={}",
                        file_id,
                        model.kind,
                        model.intercept_ppm,
                        model.slope_ppm_per_min,
                        calibration_psms.len(),
                    );
                }
                if let Some(model) = fragment {
                    log::info!(
                        "- file {} library fragment calibration: {:?}, offset={:.3} ppm, slope={:.4} ppm/min, n={}",
                        file_id,
                        model.kind,
                        model.intercept_ppm,
                        model.slope_ppm_per_min,
                        calibration_psms.len(),
                    );
                }
                LibraryMassCalibration {
                    precursor,
                    fragment,
                }
            })
            .collect()
    }

    fn align_library_properties(&self, features: &mut [Feature]) -> (Vec<Alignment>, usize, usize) {
        let mut rt_alignments = Vec::with_capacity(self.parameters.mzml_paths.len());
        let mut mobility_alignments = Vec::with_capacity(self.parameters.mzml_paths.len());
        let mut rt_files_aligned = 0usize;
        let mut mobility_files_aligned = 0usize;

        for file_id in 0..self.parameters.mzml_paths.len() {
            let landmarks = features.iter().filter(|feature| {
                feature.file_id == file_id
                    && feature.rank == 1
                    && feature.label == 1
                    && feature.spectrum_q <= 0.01
            });
            let rt_points = landmarks
                .clone()
                .map(|feature| (feature.rt, feature.predicted_rt))
                .collect::<Vec<_>>();
            let mobility_points = landmarks
                .map(|feature| (feature.ims, feature.predicted_ims))
                .collect::<Vec<_>>();
            let rt = fit_reference_alignment(&rt_points, 16);
            let mobility = fit_reference_alignment(&mobility_points, 16);
            if let Some(alignment) = rt {
                rt_files_aligned += 1;
                log::info!(
                    "- file {} library RT alignment: slope={:.5}, intercept={:.3}, inliers={}/{}",
                    file_id,
                    alignment.slope,
                    alignment.intercept,
                    alignment.inliers,
                    alignment.points,
                );
            }
            if let Some(alignment) = mobility {
                mobility_files_aligned += 1;
                log::info!(
                "- file {} library mobility alignment: slope={:.5}, intercept={:.5}, inliers={}/{}",
                file_id,
                alignment.slope,
                alignment.intercept,
                alignment.inliers,
                alignment.points,
            );
            }
            let max_observed_rt = features
                .iter()
                .filter(|feature| feature.file_id == file_id)
                .map(|feature| feature.rt)
                .filter(|rt| rt.is_finite() && *rt > 0.0)
                .fold(1.0, f32::max);
            rt_alignments.push(rt.unwrap_or(
                sage_core::ml::retention_alignment::ReferenceAlignment {
                    slope: 1.0 / max_observed_rt,
                    intercept: 0.0,
                    points: rt_points.len(),
                    inliers: 0,
                },
            ));
            mobility_alignments.push(mobility);
        }

        features.par_iter_mut().for_each(|feature| {
            feature.aligned_rt = rt_alignments
                .get(feature.file_id)
                .copied()
                .map(|alignment| alignment.transform(feature.rt))
                .unwrap_or(feature.rt);
            feature.delta_rt_model =
                if feature.predicted_rt.is_finite() && feature.predicted_rt > 0.0 {
                    (feature.aligned_rt - feature.predicted_rt).abs()
                } else {
                    0.0
                };
            feature.delta_ims_model = mobility_alignments
                .get(feature.file_id)
                .and_then(|alignment| *alignment)
                .filter(|_| feature.ims.is_finite() && feature.ims > 0.0)
                .map(|alignment| alignment.transform(feature.ims))
                .filter(|_| feature.predicted_ims.is_finite() && feature.predicted_ims > 0.0)
                .map(|aligned| (aligned - feature.predicted_ims).abs())
                .unwrap_or_default();
        });

        let lfq_alignments = rt_alignments
            .iter()
            .enumerate()
            .map(|(file_id, alignment)| alignment.for_lfq(file_id))
            .collect();
        (lfq_alignments, rt_files_aligned, mobility_files_aligned)
    }

    pub(super) fn run_library_with_summary(
        mut self,
        parallel: usize,
    ) -> anyhow::Result<(telemetry::Telemetry, RunSummary)> {
        let (target_entries, target_transitions) = self
            .library_search
            .as_ref()
            .map(|runtime| (runtime.target_entries, runtime.target_transitions))
            .expect("library runtime checked before dispatch");
        let uncalibrated =
            vec![LibraryMassCalibration::default(); self.parameters.mzml_paths.len()];
        let mut outputs = self.batch_library_files(parallel, &uncalibrated, false)?;
        self.cancellation.check()?;
        self.events.check()?;
        sort_features_by_discriminant(&mut outputs.features);
        sage_core::ml::qvalue::spectrum_q_value(&mut outputs.features);
        let calibrations = self.fit_library_mass_calibrations(&outputs.features);
        let mass_alignment_applied = calibrations
            .iter()
            .any(|calibration| calibration.precursor.is_some() || calibration.fragment.is_some());
        let quantification_enabled =
            self.parameters.quant.lfq || self.parameters.quant.tmt.is_some();
        if mass_alignment_applied || quantification_enabled {
            outputs = self.batch_library_files(parallel, &calibrations, true)?;
            self.cancellation.check()?;
            self.events.check()?;
        }
        self.events.emit(EventKind::MassAlignmentCompleted {
            files: self.parameters.mzml_paths.len(),
        });
        sort_features_by_discriminant(&mut outputs.features);
        assign_psm_ids(&mut outputs.features);
        sage_core::ml::qvalue::spectrum_q_value(&mut outputs.features);
        let (library_rt_alignments, rt_files_aligned, mobility_files_aligned) =
            self.align_library_properties(&mut outputs.features);
        let library_rescoring_fitted =
            sage_core::ml::linear_discriminant::score_library_psms(&mut outputs.features).is_some();
        if library_rescoring_fitted {
            self.events.emit(EventKind::LdaScoringCompleted);
        } else {
            let message =
                "insufficient target-decoy evidence for library rescoring, using spectral angle"
                    .to_string();
            log::warn!("{message}");
            self.events.emit(EventKind::Warning {
                code: "library_rescoring_fallback".into(),
                message,
            });
        }
        sort_features_by_discriminant(&mut outputs.features);
        let q_spectrum = sage_core::ml::qvalue::spectrum_q_value(&mut outputs.features);

        let q_peptide = assign_entity_q_values(
            &mut outputs.features,
            |feature| {
                let peptide_idx = self
                    .database
                    .decoy_pairing
                    .get(feature.peptide_idx.0 as usize)
                    .copied()
                    .unwrap_or(feature.peptide_idx);
                Some(self.database[peptide_idx].label_group())
            },
            |feature, q| feature.peptide_q = q,
        );
        let q_protein = assign_entity_q_values(
            &mut outputs.features,
            |feature| {
                let peptide = &self.database[feature.peptide_idx];
                (peptide.proteins.len() == 1).then(|| peptide.proteins[0].to_string())
            },
            |feature, q| feature.protein_q = q,
        );
        sage_core::protein_grouping::generate_protein_groups(
            &self.database,
            &mut outputs.features,
            self.parameters.protein_grouping,
            Some(self.parameters.protein_grouping_peptide_fdr),
        );
        let decoy_tag = self.database.decoy_tag.clone();
        let q_protein_group = assign_entity_q_values(
            &mut outputs.features,
            |feature| {
                (feature.num_protein_groups == 1).then(|| {
                    feature
                        .protein_groups
                        .as_deref()
                        .unwrap_or_default()
                        .split(';')
                        .map(|protein| protein.strip_prefix(&decoy_tag).unwrap_or(protein))
                        .collect::<Vec<_>>()
                        .join(";")
                })
            },
            |feature, q| feature.protein_group_q = q,
        );
        self.events.emit(EventKind::FdrCompleted {
            psms: q_spectrum,
            peptides: q_peptide,
            proteins: q_protein,
            protein_groups: q_protein_group,
        });

        let areas = if self.parameters.quant.lfq {
            let mut areas = sage_core::lfq::build_feature_map(
                self.parameters.quant.lfq_settings,
                self.parameters.precursor_charge,
                &outputs.features,
                &self.database,
            )
            .quantify(&self.database, &outputs.ms1, &library_rt_alignments);
            let q_precursor = sage_core::fdr::picked_precursor(&mut areas);
            log::info!(
                "discovered {} target library-search MS1 peaks at 5% FDR",
                q_precursor
            );
            self.events.emit(EventKind::QuantificationCompleted {
                kind: "lfq".into(),
                features: areas.len(),
            });
            Some(areas)
        } else {
            None
        };
        if !outputs.quant.is_empty() {
            self.events.emit(EventKind::QuantificationCompleted {
                kind: "tmt".into(),
                features: outputs.quant.len(),
            });
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
        let output_psm_q_value = self.parameters.output_filter.psm_q_value;
        let output_features = outputs
            .features
            .iter()
            .filter(|feature| passes_output_filter(feature, output_psm_q_value))
            .collect::<Vec<_>>();
        let bytes = sage_cloudpath::parquet::serialize_features(
            &output_features,
            &outputs.quant,
            &filenames,
            &self.database,
            output_psm_q_value,
        )?;
        let path = self.make_path("results.sage.parquet");
        sage_cloudpath::write_bytes_sync(&path, bytes)?;
        self.parameters.output_paths.push(path);

        if self.parameters.annotate_matches {
            let bytes = sage_cloudpath::parquet::serialize_matched_fragments(
                &output_features,
                output_psm_q_value,
            )?;
            let path = self.make_path("matched_fragments.sage.parquet");
            sage_cloudpath::write_bytes_sync(&path, bytes)?;
            self.parameters.output_paths.push(path);
            self.events.emit(EventKind::FragmentAnnotationCompleted {
                psms: output_features.len(),
                fragments: output_features
                    .iter()
                    .filter_map(|feature| feature.fragments.as_ref())
                    .map(|fragments| fragments.fragment_ordinals.len())
                    .sum(),
            });
        }

        if let Some(areas) = &areas {
            let bytes = sage_cloudpath::parquet::serialize_lfq(areas, &filenames, &self.database)?;
            let path = self.make_path("lfq.parquet");
            sage_cloudpath::write_bytes_sync(&path, bytes)?;
            self.parameters.output_paths.push(path);
        }

        let path = self.make_path("results.json");
        let bytes = serde_json::to_vec_pretty(&self.parameters)?;
        sage_cloudpath::write_bytes_sync(&path, bytes)?;
        self.parameters.output_paths.push(path);

        let runtime_secs = self.start.elapsed().as_secs();
        let summary_path = self.make_path("run-summary.json");
        let mut output_paths = self
            .parameters
            .output_paths
            .iter()
            .map(ToString::to_string)
            .collect::<Vec<_>>();
        output_paths.push(summary_path.to_string());
        let mut input_stats = InputRunStats::default();
        for path in &self.parameters.mzml_paths {
            let path = path.path().trim_end_matches('/').to_ascii_lowercase();
            if path.ends_with(".raw") {
                input_stats.thermo_raw_files += 1;
            } else if path.ends_with(".mzml") || path.ends_with(".mzml.gz") {
                input_stats.mzml_files += 1;
            } else {
                input_stats.other_files += 1;
            }
        }
        let library_path = self
            .parameters
            .library_search
            .as_ref()
            .map(|settings| settings.path.clone());
        let summary = RunSummary {
            schema_version: 7,
            runtime_secs,
            files: self.parameters.mzml_paths.len(),
            peptides_in_database: 0,
            fragments_in_database: 0,
            psms_at_one_percent_fdr: q_spectrum,
            peptides_at_one_percent_fdr: q_peptide,
            proteins_at_one_percent_fdr: q_protein,
            protein_groups_at_one_percent_fdr: q_protein_group,
            ptm_localization: PtmLocalizationRunStats::default(),
            models: ModelRunStats {
                mass_alignment_applied,
                ion_mobility_observed: outputs
                    .features
                    .iter()
                    .any(|feature| feature.ims.is_finite() && feature.ims > 0.0),
                library_retention_time_alignment: (rt_files_aligned > 0).then(|| "linear".into()),
                library_retention_time_files_aligned: rt_files_aligned,
                library_ion_mobility_alignment: (mobility_files_aligned > 0)
                    .then(|| "linear".into()),
                library_ion_mobility_files_aligned: mobility_files_aligned,
                library_rescoring: Some(if library_rescoring_fitted {
                    "linear_discriminant".into()
                } else {
                    "spectral_angle_fallback".into()
                }),
                ..ModelRunStats::default()
            },
            quantification: QuantificationRunStats {
                lfq_enabled: self.parameters.quant.lfq,
                lfq_features: areas.as_ref().map(HashMap::len).unwrap_or_default(),
                tmt: self
                    .parameters
                    .quant
                    .tmt
                    .as_ref()
                    .map(|tmt| format!("{tmt:?}").to_lowercase()),
                tmt_features: outputs.quant.len(),
                ms1_label_channels: self.database.label_channels.len(),
                ms1_label_reference: self.database.label_reference.as_deref().map(str::to_owned),
            },
            execution: ExecutionRunStats {
                batch_size: self.parameters.batch_size,
                parallelism: parallel,
                max_memory_gb: self.parameters.max_memory_gb,
                min_free_memory_gb: self.parameters.min_free_memory_gb,
            },
            inputs: input_stats,
            modifications: ModificationRunStats {
                label_channels: self.database.label_channels.len(),
                labeled_peptides: self
                    .database
                    .peptides
                    .iter()
                    .filter(|peptide| !peptide.decoy && peptide.label_channel.is_some())
                    .count(),
                ..ModificationRunStats::default()
            },
            spectral_library: SpectralLibraryRunStats::default(),
            library_search: LibrarySearchRunStats {
                enabled: true,
                path: library_path,
                target_entries,
                decoy_entries: target_entries,
                transitions: target_transitions,
            },
            output_paths,
        };
        sage_cloudpath::write_bytes_sync(&summary_path, serde_json::to_vec_pretty(&summary)?)?;
        self.parameters.output_paths.push(summary_path);
        for path in &self.parameters.output_paths {
            self.events.emit(EventKind::OutputWritten {
                path: path.to_string(),
            });
        }
        let telemetry = telemetry::Telemetry::new(
            self.parameters,
            target_entries,
            target_transitions,
            runtime_secs,
        );
        self.events.emit(EventKind::JobCompleted {
            runtime_secs,
            outputs: summary.output_paths.len(),
        });
        self.events.check()?;
        Ok((telemetry, summary))
    }
}

#[cfg(test)]
#[path = "../../tests/unit/runner/library_runner.rs"]
mod tests;
