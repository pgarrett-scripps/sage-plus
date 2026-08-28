use super::*;

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
    pub(super) fn align_mass_errors(&self, features: &mut [Feature]) {
        features.par_iter_mut().for_each(|feature| {
            feature.aligned_delta_mass = feature.delta_mass;
            feature.aligned_average_ppm = feature.average_ppm;
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

            let precursor_fit = matches!(self.parameters.precursor_tol, Tolerance::Ppm(_, _))
                .then(|| fit_mass_calibration(&precursor_points, fit_options))
                .flatten();
            let fragment_fit = fit_mass_calibration(&fragment_points, fit_options);

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
    }

    // Create a path for `file_name` in the specified output directory, if it exists,
    // otherwise, write to current directory
    pub(super) fn make_path<S: AsRef<str>>(&self, file_name: S) -> Url {
        self.parameters
            .output_directory
            .join(file_name.as_ref())
            .expect("valid path segment")
    }

    pub(super) fn search_processed_spectra(
        &self,
        scorer: &Scorer,
        msn_spectra: &[ProcessedSpectrum],
    ) -> Vec<Feature> {
        use std::sync::atomic::{AtomicUsize, Ordering};
        let counter = AtomicUsize::new(0);
        let start = Instant::now();

        let features: Vec<_> = msn_spectra
            .par_iter()
            .filter(|spec| {
                !self.cancellation.is_cancelled()
                    && spec.masses.len() >= self.parameters.min_peaks
                    && spec.level == 2
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
            .flat_map(|spec| scorer.score(spec))
            .collect();

        let duration = Instant::now().duration_since(start).as_millis() as usize;
        let prev = counter.load(Ordering::Relaxed);
        let rate = prev * 1000 / (duration + 1);
        log::info!("- search:  {:8} ms ({} spectra/s)", duration, rate);
        features
    }

    pub(super) fn complete_features(
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
        let spectra = self.read_processed_spectra(chunk, chunk_idx, batch_size)?;
        let features = self.search_processed_spectra(scorer, &spectra.1);
        Ok(self.complete_features(spectra.1, spectra.0, features))
    }

    pub(super) fn read_processed_spectra(
        &self,
        chunk: &[Url],
        chunk_idx: usize,
        batch_size: usize,
    ) -> anyhow::Result<(Vec<ProcessedSpectrum>, Vec<ProcessedSpectrum>)> {
        self.read_processed_spectra_with_ms1(chunk, chunk_idx, batch_size, self.requires_ms1())
    }

    pub(super) fn read_processed_spectra_with_ms1(
        &self,
        chunk: &[Url],
        chunk_idx: usize,
        batch_size: usize,
        requires_ms1: bool,
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
            self.events.emit(EventKind::FileStarted {
                file_id,
                path: path.to_string(),
            });
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
                    self.events.emit(EventKind::FileCompleted {
                        file_id,
                        path: path.to_string(),
                        spectra: spectra.ms1.len() + spectra.msn.len(),
                    });
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

        self.events.emit(EventKind::SpectraProcessed {
            ms1_spectra: spectra.ms1.len(),
            msn_spectra: spectra.msn.len(),
        });

        let io_time = Instant::now() - start;
        info!("- file IO: {:8} ms", io_time.as_millis());

        Ok((spectra.ms1, spectra.msn))
    }
}
