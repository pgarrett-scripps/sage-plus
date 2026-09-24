use super::*;

/// Most MS2 spectra per file searched in the discovery pass. Each
/// acquisition group keeps its share of the sample, taken at even positions
/// across its own spectra, so the sample spans the whole gradient and does
/// not alias with interleaved scan cycles.
const MAX_DISCOVERY_SPECTRA: usize = 25_000;

/// Largest correction, in ppm, a model may apply for `tolerance`: the
/// window half-width for ppm tolerances, or the Da half-width at m/z 2000.
pub(super) fn recalibration_cap_ppm(tolerance: Tolerance) -> Option<f32> {
    match tolerance {
        Tolerance::Ppm(lo, hi) => Some(lo.abs().max(hi.abs())),
        Tolerance::Da(lo, hi) => Some(lo.abs().max(hi.abs()) * 1e6 / 2000.0),
        Tolerance::Pct(_, _) => None,
    }
}

/// Widen a tolerance by its own half-width on both sides, so a window that
/// is later moved by at most the recalibration cap stays inside it.
pub(super) fn widen_for_recalibration(tolerance: Tolerance) -> Tolerance {
    match tolerance {
        Tolerance::Ppm(lo, hi) => {
            let m = lo.abs().max(hi.abs());
            Tolerance::Ppm(lo - m, hi + m)
        }
        Tolerance::Da(lo, hi) => {
            let m = lo.abs().max(hi.abs());
            Tolerance::Da(lo - m, hi + m)
        }
        other => other,
    }
}

impl Runner {
    /// Whether the search runs a discovery pass and recalibrates masses.
    pub(super) fn mass_recalibration_enabled(&self) -> bool {
        self.parameters.mass_recalibration != MassRecalibrationMode::Off
            && !self.parameters.wide_window
    }

    /// Corrections selected so far, indexed by file.
    pub(super) fn mass_recalibration_models(&self) -> Option<Arc<MassRecalibration>> {
        let stats = self
            .mass_recalibration
            .lock()
            .expect("mass recalibration lock");
        if stats.iter().all(|file| file.correction.is_identity()) {
            return None;
        }
        let mut files = vec![FileMassCorrection::default(); self.parameters.mzml_paths.len()];
        for file in stats.iter() {
            if let Some(slot) = files.get_mut(file.file_id) {
                *slot = file.correction.clone();
            }
        }
        Some(Arc::new(MassRecalibration { files }))
    }

    pub(super) fn mass_recalibration_stats(&self) -> Vec<MassRecalibrationFileStats> {
        let mut stats = self
            .mass_recalibration
            .lock()
            .expect("mass recalibration lock")
            .clone();
        stats.sort_by_key(|file| file.file_id);
        stats
    }

    /// Discovery pass for every file in a batch: search a sample of each
    /// file's MS2 spectra, stratified by acquisition group, with observed
    /// masses, keep rank-1 targets at 1% Poisson spectrum q-value computed
    /// within each acquisition group, and select a per-file precursor model
    /// and per-group fragment models. Models are fitted on targets only and are applied to
    /// every spectrum of the file, so targets and decoys are treated alike.
    fn discover_mass_corrections(&self, scorer: &Scorer, spectra: &[ProcessedSpectrum]) {
        let mut file_ids = spectra
            .iter()
            .filter(|spectrum| spectrum.level == 2)
            .map(|spectrum| spectrum.file_id)
            .collect::<Vec<_>>();
        file_ids.sort_unstable();
        file_ids.dedup();
        let max_kind = self.parameters.mass_recalibration.max_kind();

        for file_id in file_ids {
            let start = Instant::now();
            let candidates = spectra
                .iter()
                .filter(|spectrum| {
                    spectrum.file_id == file_id && spectrum.is_searchable(self.parameters.min_peaks)
                })
                .collect::<Vec<_>>();
            let searchable = candidates.len();
            let keys = candidates
                .iter()
                .map(|spectrum| spectrum.acquisition)
                .collect::<Vec<_>>();
            let sample = stratified_sample(&keys, MAX_DISCOVERY_SPECTRA)
                .into_iter()
                .map(|index| candidates[index])
                .collect::<Vec<_>>();
            let features = sample
                .par_iter()
                .enumerate()
                .flat_map_iter(|(index, spectrum)| {
                    scorer
                        .score(spectrum)
                        .into_iter()
                        .filter(|feature| feature.rank == 1)
                        .map(move |feature| (index, feature))
                })
                .collect::<Vec<_>>();
            // Q-values are computed within each acquisition group: pooled,
            // a low-accuracy group searched with a wide fragment window
            // floods the decoy counts and leaves almost no confident PSMs in
            // an accurate group. The precursor model uses the union; each
            // fragment model sees only its own group's PSMs.
            let confident = confident_per_group(
                features
                    .into_iter()
                    .map(|(index, feature)| (sample[index].acquisition, index, feature))
                    .collect(),
                0.01,
            );

            let precursor_points = confident
                .iter()
                .filter(|(_, feature)| {
                    feature.isotope_error == 0.0 && feature.mass_offset.is_none()
                })
                .map(|(_, feature)| {
                    let z = feature.charge.max(1) as f32;
                    let observed = feature.expmass / z + PROTON;
                    let theoretical = feature.calcmass / z + PROTON;
                    MassErrorPoint {
                        rt_minutes: feature.rt,
                        mz: observed,
                        error_ppm: (observed - theoretical) * 1e6 / theoretical,
                        group: stable_hash(&feature.spec_id),
                    }
                })
                .collect::<Vec<_>>();
            let fragment_points = confident
                .par_iter()
                .filter(|(_, feature)| feature.mass_offset.is_none())
                .flat_map_iter(|(index, feature)| {
                    let spectrum = sample[*index];
                    let fragments = scorer.annotate_candidate(spectrum, feature);
                    let acquisition = spectrum.acquisition;
                    let group = stable_hash(&feature.spec_id);
                    let rt = feature.rt;
                    fragments
                        .mz_experimental
                        .into_iter()
                        .zip(fragments.mz_calculated)
                        .map(move |(observed, theoretical)| {
                            let point = MassErrorPoint {
                                rt_minutes: rt,
                                mz: observed,
                                error_ppm: (observed - theoretical) * 1e6 / theoretical,
                                group,
                            };
                            (acquisition, point)
                        })
                })
                .collect::<Vec<_>>();
            let mut groups: Vec<(AcquisitionGroup, usize)> = Vec::new();
            for spectrum in &sample {
                match groups.iter_mut().find(|(g, _)| *g == spectrum.acquisition) {
                    Some((_, count)) => *count += 1,
                    None => groups.push((spectrum.acquisition, 1)),
                }
            }

            let options = |cap: Option<f32>, min_mz_span: f32| RecalibrationOptions {
                max_kind: if cap.is_some() {
                    max_kind
                } else {
                    MassModelKind::None
                },
                max_abs_ppm: cap.unwrap_or(0.0),
                min_mz_span,
                ..RecalibrationOptions::default()
            };
            // Precursor corrections require a ppm window, as in post-hoc alignment.
            let precursor_cap = match self.parameters.precursor_tol {
                Tolerance::Ppm(_, _) => recalibration_cap_ppm(self.parameters.precursor_tol),
                _ => None,
            };
            let precursor = select_model(&precursor_points, options(precursor_cap, 100.0));
            let fragment = select_group_models(
                &fragment_points,
                &groups,
                options(recalibration_cap_ppm(self.parameters.fragment_tol), 200.0),
            );
            let correction = FileMassCorrection {
                precursor: precursor.model.clone(),
                fragment: fragment
                    .iter()
                    .map(|group| GroupMassCorrection {
                        group: group.group,
                        model: group.selection.model.clone(),
                    })
                    .collect(),
            };
            let describe = |selection: &ModelSelection| match &selection.model {
                Some(model) => format!(
                    "{:?}/{:?} ({} parameters, offset {:.2} ppm)",
                    model.kind,
                    model.axes,
                    model.parameters(),
                    model.intercept_ppm
                )
                .to_lowercase(),
                None => format!(
                    "none ({})",
                    selection
                        .skipped
                        .as_deref()
                        .unwrap_or("no validated improvement")
                ),
            };
            info!(
                "- file {} mass recalibration: {} of {} spectra searched, {} confident PSMs; precursor {}; fragment {} [{} ms]",
                file_id,
                sample.len(),
                searchable,
                confident.len(),
                describe(&precursor),
                fragment
                    .iter()
                    .map(|group| format!(
                        "{} ({} spectra, {} PSMs) {}",
                        group.group.label(),
                        group.spectra,
                        group.selection.psms,
                        describe(&group.selection)
                    ))
                    .collect::<Vec<_>>()
                    .join(", "),
                start.elapsed().as_millis(),
            );
            self.mass_recalibration
                .lock()
                .expect("mass recalibration lock")
                .push(MassRecalibrationFileStats {
                    file_id,
                    discovery_spectra: sample.len(),
                    discovery_psms: confident.len(),
                    discovery_ms: start.elapsed().as_millis() as u64,
                    precursor,
                    fragment,
                    correction,
                });
        }
    }
}

impl Runner {
    pub(super) fn spectrum_fdr(&self, features: &mut Vec<Feature>) -> usize {
        if sage_core::ml::linear_discriminant::score_psms(features, self.parameters.precursor_tol)
            .is_none()
        {
            log::warn!("linear model fitting failed, falling back to heuristic discriminant score");
            self.events.emit(EventKind::Warning {
                code: "discriminant_model_fallback".into(),
                message: "linear model fitting failed; using heuristic discriminant score".into(),
            });
            features.par_iter_mut().for_each(|feat| {
                feat.discriminant_score = (-feat.poisson as f32).ln_1p() + feat.longest_y_pct / 3.0
            });
        }
        sort_features_by_discriminant(features);
        sage_core::ml::qvalue::spectrum_q_value(features)
    }

    /// Align systematic precursor and fragment mass errors per raw file before
    /// fitting the final FDR model. Models are trained only on provisional 1%
    /// spectrum-q rank-1 targets and are then applied equally to targets and
    /// decoys. Raw output errors remain unchanged.
    pub(super) fn align_mass_errors(
        &self,
        features: &mut [Feature],
    ) -> Vec<MassAlignmentFileStats> {
        let mut diagnostics = Vec::with_capacity(self.parameters.mzml_paths.len());
        // Files corrected at search time already carry residual errors in the
        // aligned columns; post-hoc alignment only handles the rest.
        let recalibration = self.mass_recalibration_models();
        let corrected = |file_id: usize| {
            let correction = recalibration.as_deref().and_then(|r| r.file(file_id));
            (
                correction.is_some_and(|c| c.precursor.as_ref().is_some_and(|m| !m.is_identity())),
                correction.is_some_and(|c| c.corrects_fragments()),
            )
        };
        features.par_iter_mut().for_each(|feature| {
            let (precursor, fragment) = corrected(feature.file_id);
            if !precursor {
                feature.aligned_delta_mass = feature.delta_mass;
            }
            if !fragment {
                feature.aligned_average_ppm = feature.average_ppm;
            }
        });

        let fit_options = FitOptions {
            // A line is useful even for modest drift; reject it only when its
            // held-in robust residual is worse than the static center.
            min_linear_improvement: 0.0,
            ..FitOptions::default()
        };

        for file_id in 0..self.parameters.mzml_paths.len() {
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

            let (precursor_corrected, fragment_corrected) = corrected(file_id);
            let precursor_fit = (!precursor_corrected
                && matches!(self.parameters.precursor_tol, Tolerance::Ppm(_, _)))
            .then(|| fit_mass_calibration(&precursor_points, fit_options))
            .flatten();
            let fragment_fit = (!fragment_corrected)
                .then(|| fit_mass_calibration(&fragment_points, fit_options))
                .flatten();
            diagnostics.push(MassAlignmentFileStats {
                file_id,
                calibration_psms: calibration_psms.len(),
                precursor: precursor_fit.map(|fit| fit.model),
                fragment: fragment_fit.map(|fit| fit.model),
                precursor_skip_reason: precursor_fit.is_none().then(|| {
                    if precursor_corrected {
                        "corrected by search-time mass recalibration".into()
                    } else if matches!(self.parameters.precursor_tol, Tolerance::Ppm(_, _)) {
                        "insufficient finite high-confidence observations".into()
                    } else {
                        "precursor alignment requires ppm tolerance".into()
                    }
                }),
                fragment_skip_reason: fragment_fit.is_none().then(|| {
                    if fragment_corrected {
                        "corrected by search-time mass recalibration".into()
                    } else {
                        "insufficient finite high-confidence observations".into()
                    }
                }),
            });

            if let Some(fit) = precursor_fit {
                log::info!(
                    "- file {} precursor mass alignment: {:?}, offset={:.3} ppm, slope={:.4} ppm/min, n={}",
                    file_id,
                    fit.model.kind,
                    fit.model.intercept_ppm,
                    fit.model.slope_ppm_per_min,
                    fit.inliers,
                );
            }
            if let Some(fit) = fragment_fit {
                log::info!(
                    "- file {} fragment mass alignment: {:?}, offset={:.3} ppm, slope={:.4} ppm/min, n={}",
                    file_id,
                    fit.model.kind,
                    fit.model.intercept_ppm,
                    fit.model.slope_ppm_per_min,
                    fit.inliers,
                );
            }

            features
                .iter_mut()
                .filter(|feature| feature.file_id == file_id)
                .for_each(|feature| {
                    if let Some(fit) = precursor_fit {
                        feature.aligned_delta_mass =
                            feature.delta_mass - fit.model.predict_ppm(feature.rt);
                    }
                    if let Some(fit) = fragment_fit {
                        let predicted = fit.model.predict_ppm(feature.rt);
                        // Preserve the within-PSM absolute-error spread while
                        // translating its signed center to the fitted baseline.
                        feature.aligned_average_ppm = align_fragment_error(
                            feature.average_ppm,
                            feature.signed_fragment_ppm,
                            predicted,
                        );
                    }
                });
        }
        diagnostics
    }

    // Create a path for `file_name` in the specified output directory, if it exists,
    // otherwise, write to current directory
    pub(super) fn make_path<S: AsRef<str>>(&self, file_name: S) -> Url {
        self.parameters
            .output_directory
            .join(file_name.as_ref())
            .expect("valid path segment")
    }

    /// Score MS2 spectra. Also returns search-time `psm_id` -> occurrence for
    /// PSMs from spectra whose ID repeats within a file, so the post-FDR pass
    /// can match each PSM to the spectrum it was scored against.
    pub(super) fn search_processed_spectra(
        &self,
        scorer: &Scorer,
        msn_spectra: &[ProcessedSpectrum],
    ) -> (Vec<Feature>, HashMap<usize, usize>) {
        use std::sync::atomic::{AtomicUsize, Ordering};
        let counter = AtomicUsize::new(0);
        let start = Instant::now();
        let occurrences = spectrum_id_occurrences(msn_spectra);
        let repeated = std::sync::Mutex::new(HashMap::new());

        let features: Vec<_> = msn_spectra
            .par_iter()
            .zip(occurrences.par_iter())
            .filter(|(spec, _)| {
                !self.cancellation.is_cancelled() && spec.is_searchable(self.parameters.min_peaks)
            })
            .map(|x| {
                let prev = counter.fetch_add(1, Ordering::Relaxed);
                if prev > 0 && prev.is_multiple_of(10_000) {
                    let duration = Instant::now().duration_since(start).as_millis() as usize;

                    let rate = prev * 1000 / (duration + 1);
                    log::trace!("- searched {} spectra ({} spectra/s)", prev, rate);
                }
                x
            })
            .flat_map(|(spec, &occurrence)| {
                let features = scorer.score(spec);
                if occurrence > 0 {
                    repeated
                        .lock()
                        .expect("repeated spectrum map")
                        .extend(features.iter().map(|feature| (feature.psm_id, occurrence)));
                }
                features
            })
            .collect();

        let duration = Instant::now().duration_since(start).as_millis() as usize;
        let prev = counter.load(Ordering::Relaxed);
        let rate = prev * 1000 / (duration + 1);
        log::info!("- search:  {:8} ms ({} spectra/s)", duration, rate);
        (
            features,
            repeated.into_inner().expect("repeated spectrum map"),
        )
    }

    pub(super) fn complete_features(
        &self,
        msn_spectra: Vec<ProcessedSpectrum>,
        ms1_spectra: Vec<ProcessedSpectrum>,
        features: Vec<Feature>,
        repeated_spectrum_psms: HashMap<usize, usize>,
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
            repeated_spectrum_psms,
        }
    }

    pub(super) fn requires_ms1(&self) -> bool {
        self.parameters.quant.lfq
    }

    pub(super) fn process_chunk(
        &self,
        scorer: &Scorer,
        chunk: &[Url],
        chunk_idx: usize,
        batch_size: usize,
    ) -> anyhow::Result<SageResults> {
        let spectra = match self.take_retained_spectra(chunk_idx, batch_size) {
            Some(spectra) => {
                self.emit_retained_batch_events(chunk, chunk_idx, batch_size, &spectra);
                spectra
            }
            None => self.read_processed_spectra(chunk, chunk_idx, batch_size)?,
        };
        let (features, repeated_spectrum_psms) = if self.mass_recalibration_enabled() {
            self.discover_mass_corrections(scorer, &spectra.1);
            let recalibrated = Scorer {
                mass_recalibration: self.mass_recalibration_models(),
                ..self.scorer()
            };
            self.search_processed_spectra(&recalibrated, &spectra.1)
        } else {
            self.search_processed_spectra(scorer, &spectra.1)
        };
        Ok(self.complete_features(spectra.1, spectra.0, features, repeated_spectrum_psms))
    }

    /// Take the prefilter's spectra for this file batch. Spectra retained
    /// for a different batch size are released unused.
    fn take_retained_spectra(&self, chunk_idx: usize, batch_size: usize) -> Option<SpectrumBatch> {
        let mut retained = self.retained_spectra.lock().expect("retained spectra lock");
        if retained.batch_size != batch_size {
            retained.batches = Vec::new();
            return None;
        }
        retained.batches.get_mut(chunk_idx).and_then(Option::take)
    }

    /// Emit the progress events that reading this file batch would emit.
    fn emit_retained_batch_events(
        &self,
        chunk: &[Url],
        chunk_idx: usize,
        batch_size: usize,
        spectra: &SpectrumBatch,
    ) {
        info!(
            "processing files {} .. {} (read during prefiltering)",
            batch_size * chunk_idx,
            batch_size * chunk_idx + chunk.len()
        );
        for (idx, path) in chunk.iter().enumerate() {
            let file_id = chunk_idx * batch_size + idx;
            self.events.emit(EventKind::FileStarted {
                file_id,
                path: path.to_string(),
            });
            self.events.emit(EventKind::FileCompleted {
                file_id,
                path: path.to_string(),
                spectra: spectra
                    .0
                    .iter()
                    .chain(&spectra.1)
                    .filter(|spectrum| spectrum.file_id == file_id)
                    .count(),
            });
        }
        self.events.emit(EventKind::SpectraProcessed {
            ms1_spectra: spectra.0.len(),
            msn_spectra: spectra.1.len(),
        });
    }

    pub(super) fn read_processed_spectra(
        &self,
        chunk: &[Url],
        chunk_idx: usize,
        batch_size: usize,
    ) -> anyhow::Result<(Vec<ProcessedSpectrum>, Vec<ProcessedSpectrum>)> {
        self.read_processed_spectra_with_ms1(
            chunk,
            chunk_idx,
            batch_size,
            self.requires_ms1(),
            true,
        )
    }

    /// `search_events` controls the per-file progress events
    /// (`file_started`, `file_completed`, `spectra_processed`). Rereads after
    /// the search has reported completion pass `false`; read failures are
    /// always reported.
    pub(super) fn read_processed_spectra_with_ms1(
        &self,
        chunk: &[Url],
        chunk_idx: usize,
        batch_size: usize,
        requires_ms1: bool,
        search_events: bool,
    ) -> anyhow::Result<(Vec<ProcessedSpectrum>, Vec<ProcessedSpectrum>)> {
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

        let sp = SpectrumProcessor::with_deisotope_settings(
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
        let inner_closure = |(idx, path): (usize, &Url)| {
            let file_id = chunk_idx * batch_size + idx;
            if search_events {
                self.events.emit(EventKind::FileStarted {
                    file_id,
                    path: path.to_string(),
                });
            }
            let res = sage_cloudpath::util::read_spectra(
                path,
                file_id,
                sn,
                self.parameters.bruker_config,
                requires_ms1,
            );

            match res {
                Ok(s) => {
                    if s.is_empty() {
                        let message = "input contains no spectra".to_string();
                        self.events.emit(EventKind::FileFailed {
                            file_id,
                            path: path.to_string(),
                            message: message.clone(),
                        });
                        anyhow::bail!("failed to read spectra file `{path}`: {message}");
                    }
                    log::trace!("- {}: read {} spectra", path, s.len());
                    let spectra = s
                        .into_par_iter()
                        .map(|spectrum| sp.process(spectrum))
                        .collect::<SpectrumAccumulator>();
                    if search_events {
                        self.events.emit(EventKind::FileCompleted {
                            file_id,
                            path: path.to_string(),
                            spectra: spectra.ms1.len() + spectra.msn.len(),
                        });
                    }
                    Ok(spectra)
                }
                Err(e) => {
                    log::error!("- {}: {}", path, e);
                    self.events.emit(EventKind::FileFailed {
                        file_id,
                        path: path.to_string(),
                        message: e.to_string(),
                    });
                    Err(anyhow::Error::new(e)
                        .context(format!("failed to read spectra file `{path}`")))
                }
            }
        };

        let spectra: SpectrumAccumulator = if file_serial_read {
            chunk.iter().enumerate().map(inner_closure).try_fold(
                SpectrumAccumulator::default(),
                |accumulator, spectra| {
                    Ok::<_, anyhow::Error>(SpectrumAccumulator::reduce(accumulator, spectra?))
                },
            )?
        } else {
            chunk
                .par_iter()
                .enumerate()
                .map(inner_closure)
                .try_reduce(SpectrumAccumulator::default, |left, right| {
                    Ok(SpectrumAccumulator::reduce(left, right))
                })?
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

        if search_events {
            self.events.emit(EventKind::SpectraProcessed {
                ms1_spectra: spectra.ms1.len(),
                msn_spectra: spectra.msn.len(),
            });
        }

        let io_time = Instant::now() - start;
        info!("- file IO: {:8} ms", io_time.as_millis());

        Ok((spectra.ms1, spectra.msn))
    }
}
