use super::*;

impl Runner {
    pub fn batch_files(&self, scorer: &Scorer, batch_size: usize) -> anyhow::Result<SageResults> {
        let results = self
            .parameters
            .mzml_paths
            .chunks(batch_size)
            .enumerate()
            .map(|(chunk_idx, chunk)| {
                let results = self.process_chunk(scorer, chunk, chunk_idx, batch_size);
                self.events.emit(EventKind::SearchProgress {
                    files_completed: (chunk_idx * batch_size + chunk.len())
                        .min(self.parameters.mzml_paths.len()),
                    files_total: self.parameters.mzml_paths.len(),
                });
                results
            })
            .collect::<anyhow::Result<Vec<_>>>()?;
        Ok(results.into_iter().collect())
    }

    pub fn run(self, parallel: usize) -> anyhow::Result<telemetry::Telemetry> {
        self.run_with_summary(parallel)
            .map(|(telemetry, _summary)| telemetry)
    }

    pub fn run_with_summary(
        mut self,
        parallel: usize,
    ) -> anyhow::Result<(telemetry::Telemetry, RunSummary)> {
        anyhow::ensure!(parallel > 0, "batch size must be greater than zero");
        self.cancellation.check()?;
        self.events.check()?;
        let scorer = Scorer {
            db: &self.database,
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
            report_psms: self.parameters.report_psms,
            wide_window: self.parameters.wide_window,
            annotate_matches: false,
            mass_shift_ppm: self.parameters.mass_shift_ppm,
            score_type: self.parameters.score_type,
        };

        //Collect all results into a single container
        let mut outputs = self.batch_files(&scorer, parallel)?;
        self.cancellation.check()?;
        self.events.check()?;

        // Establish provisional q-values from the search-only Poisson feature,
        // then use confident PSMs for mass-error alignment and property-model
        // training before final FDR fitting.
        outputs.features.par_sort_unstable_by(|left, right| {
            left.poisson
                .total_cmp(&right.poisson)
                .then_with(|| feature_identity_cmp(left, right))
        });
        assign_psm_ids(&mut outputs.features);
        sage_core::ml::qvalue::spectrum_q_value_by(&mut outputs.features, |feature| {
            feature.poisson
        });
        self.align_mass_errors(&mut outputs.features);
        self.events.emit(EventKind::MassAlignmentCompleted {
            files: self.parameters.mzml_paths.len(),
        });

        let has_ion_mobility = outputs
            .features
            .iter()
            .any(|feature| feature.ims.is_finite() && feature.ims > 0.0);

        let needs_alignment = self.parameters.predict_rt
            || self.parameters.quant.lfq
            || self.parameters.retention_time_alignment.is_some();
        let alignments = if needs_alignment {
            // Alignment landmarks are observed target PSMs passing 1% spectrum FDR.
            let alignment_method = self.parameters.retention_time_alignment.unwrap_or_default();
            let alignments = sage_core::ml::retention_alignment::global_alignment_with_method(
                &mut outputs.features,
                self.parameters.mzml_paths.len(),
                alignment_method,
            );
            self.events.emit(EventKind::RetentionAlignmentCompleted {
                method: format!("{alignment_method:?}").to_lowercase(),
                files: self.parameters.mzml_paths.len(),
            });
            Some(alignments)
        } else {
            None
        };

        let retention_time_model_fitted = if self.parameters.predict_rt {
            if sage_core::ml::retention_model::predict(
                &self.database,
                &mut outputs.features,
                &self.parameters.retention_time_model,
            )
            .is_some()
            {
                self.events.emit(EventKind::RtModelFitted);
                true
            } else {
                self.events.emit(EventKind::RtModelSkipped {
                    reason: "insufficient high-confidence observations or poor model fit".into(),
                });
                false
            }
        } else {
            self.events.emit(EventKind::RtModelSkipped {
                reason: "retention-time prediction disabled".into(),
            });
            false
        };
        let ion_mobility_model_fitted =
            if has_ion_mobility && self.parameters.ion_mobility_model.enabled {
                if sage_core::ml::mobility_model::predict(
                    &self.database,
                    &mut outputs.features,
                    &self.parameters.ion_mobility_model,
                )
                .is_some()
                {
                    self.events.emit(EventKind::MobilityModelFitted);
                    true
                } else {
                    self.events.emit(EventKind::MobilityModelSkipped {
                        reason: "ion mobility unavailable or model fitting failed".into(),
                    });
                    false
                }
            } else if !self.parameters.ion_mobility_model.enabled {
                self.events.emit(EventKind::MobilityModelSkipped {
                    reason: "ion-mobility prediction disabled".into(),
                });
                false
            } else {
                self.events.emit(EventKind::MobilityModelSkipped {
                    reason: "no ion-mobility observations present".into(),
                });
                false
            };

        let q_spectrum = self.spectrum_fdr(&mut outputs.features);
        self.events.emit(EventKind::LdaScoringCompleted);
        let q_peptide = sage_core::fdr::picked_peptide(&self.database, &mut outputs.features);
        // Protein FDR is based exclusively on proteotypic (unique, non-shared) peptides. Shared peptides
        // are reported with protein FDR = 1.0
        let q_protein = sage_core::fdr::picked_protein(&self.database, &mut outputs.features);
        // Conducts "IDPicker-based protein grouping at 1% peptide FDR"
        sage_core::protein_grouping::generate_protein_groups(
            &self.database,
            &mut outputs.features,
            self.parameters.protein_grouping,
            Some(self.parameters.protein_grouping_peptide_fdr),
        );
        // Uses the "Picked Group FDR" approach to compute protein group FDR for the IDPicker groups,
        // including rescued subset grouping (rsG). Shared peptides (between different groups)
        // are reported with protein group FDR = 1.0
        let q_protein_group =
            sage_core::fdr::picked_protein_group(&self.database, &mut outputs.features);
        self.events.emit(EventKind::FdrCompleted {
            psms: q_spectrum,
            peptides: q_peptide,
            proteins: q_protein,
            protein_groups: q_protein_group,
        });
        self.cancellation.check()?;

        let postprocess = self.postprocess_features(&scorer, &mut outputs.features, parallel)?;
        if self.parameters.annotate_matches {
            self.events.emit(EventKind::FragmentAnnotationCompleted {
                psms: postprocess.annotated_psms,
                fragments: postprocess.annotated_fragments,
            });
        }
        let localized_psms = if self.parameters.ptm_localization.enabled {
            self.events.emit(EventKind::PtmLocalizationCompleted {
                psms: postprocess.localized_psms,
            });
            postprocess.localized_psms
        } else {
            0
        };

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

        let library_entries = spectral_library::build_entries(
            &outputs.features,
            &self.database,
            &filenames,
            &postprocess.library_selections,
            &self.parameters.spectral_library,
        )
        .map_err(anyhow::Error::msg)?;
        let library_transitions = library_entries
            .iter()
            .map(|entry| entry.fragments.len())
            .sum::<usize>();

        let areas = alignments.and_then(|alignments| {
            if self.parameters.quant.lfq {
                log::trace!("performing LFQ");
                let mut areas = sage_core::lfq::build_feature_map(
                    self.parameters.quant.lfq_settings,
                    self.parameters.precursor_charge,
                    &outputs.features,
                    &self.database,
                )
                .quantify(&self.database, &outputs.ms1, &alignments);

                let q_precursor = sage_core::fdr::picked_precursor(&mut areas);
                self.events.emit(EventKind::QuantificationCompleted {
                    kind: "lfq".into(),
                    features: areas.len(),
                });

                log::info!("discovered {} target MS1 peaks at 5% FDR", q_precursor);
                Some(areas)
            } else {
                None
            }
        });

        if !outputs.quant.is_empty() {
            self.events.emit(EventKind::QuantificationCompleted {
                kind: "tmt".into(),
                features: outputs.quant.len(),
            });
        }
        let lfq_features = areas.as_ref().map(|areas| areas.len()).unwrap_or_default();
        let tmt_features = outputs.quant.len();

        log::info!(
            "discovered {} target peptide-spectrum matches at 1% FDR",
            q_spectrum
        );
        log::info!("discovered {} target peptides at 1% FDR", q_peptide);
        log::info!(
            "discovered {} target proteins (supported by proteotypic peptides only) at 1% FDR",
            q_protein
        );
        log::info!(
            "discovered {} target protein groups (supported by proteotypic peptides only) at 1% FDR",
            q_protein_group
        );
        log::trace!("writing outputs");

        let output_psm_q_value = self.parameters.output_filter.psm_q_value;
        let output_features = outputs
            .features
            .iter()
            .filter(|feature| passes_output_filter(feature, output_psm_q_value))
            .collect::<Vec<_>>();
        log::info!(
            "writing {} of {} PSMs at spectrum q-value <= {}",
            output_features.len(),
            outputs.features.len(),
            output_psm_q_value
        );

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
        }

        if self
            .parameters
            .spectral_library
            .writes(SpectralLibraryFormat::SageParquet)
        {
            let bytes = sage_cloudpath::parquet::serialize_spectral_library(
                &library_entries,
                &self.parameters.spectral_library,
            )?;
            let path = self.make_path("spectral_library.sage.parquet");
            sage_cloudpath::write_bytes_sync(&path, bytes)?;
            self.parameters.output_paths.push(path);
        }
        if self
            .parameters
            .spectral_library
            .writes(SpectralLibraryFormat::MzSpecLib)
        {
            let bytes = spectral_library::serialize_mzspeclib(
                &library_entries,
                self.parameters.version.as_str(),
                self.parameters.spectral_library.strategy,
            );
            let path = self.make_path("spectral_library.mzspeclib.txt");
            sage_cloudpath::write_bytes_sync(&path, bytes)?;
            self.parameters.output_paths.push(path);
        }
        if self.parameters.spectral_library.enabled {
            self.events.emit(EventKind::SpectralLibraryCompleted {
                entries: library_entries.len(),
                transitions: library_transitions,
                formats: self.parameters.spectral_library.formats.len(),
            });
        }

        if let Some(areas) = &areas {
            let bytes = sage_cloudpath::parquet::serialize_lfq(areas, &filenames, &self.database)?;

            let path = self.make_path("lfq.parquet");
            sage_cloudpath::write_bytes_sync(&path, bytes)?;
            self.parameters.output_paths.push(path);
        }

        // PTM site reports follow the selected main output format.
        if self.parameters.ptm_localization.enabled {
            self.parameters
                .output_paths
                .push(self.write_ptm_sites(&outputs.features, &filenames)?);
            self.parameters
                .output_paths
                .push(self.write_protein_sites(&outputs.features, &filenames)?);
            self.parameters
                .output_paths
                .extend(self.write_ptm_library(&outputs.features, &filenames)?);
        }

        // Write percolator input file if requested
        if self.parameters.write_pin {
            self.parameters
                .output_paths
                .push(self.write_pin(&outputs.features, &filenames)?);
        }

        // Write an html report if requested
        if self.parameters.write_report {
            self.parameters.output_paths.push(self.write_report(
                &outputs.features,
                areas,
                &filenames,
            )?);
        }

        let path = self.make_path("results.json");
        if !self.events.is_enabled() {
            println!("{}", serde_json::to_string_pretty(&self.parameters)?);
        }

        let bytes = serde_json::to_vec_pretty(&self.parameters)?;
        sage_cloudpath::write_bytes_sync(&path, bytes)?;
        self.parameters.output_paths.push(path);

        let run_time = (Instant::now() - self.start).as_secs();
        info!("finished in {}s", run_time);
        info!("cite: \"Sage: An Open-Source Tool for Fast Proteomics Searching and Quantification at Scale\" https://doi.org/10.1021/acs.jproteome.3c00486");

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
        let summary = RunSummary {
            schema_version: 8,
            runtime_secs: run_time,
            files: self.parameters.mzml_paths.len(),
            peptides_in_database: self.database.peptides.len(),
            fragments_in_database: self.database.fragments.len(),
            psms_at_one_percent_fdr: q_spectrum,
            peptides_at_one_percent_fdr: q_peptide,
            proteins_at_one_percent_fdr: q_protein,
            protein_groups_at_one_percent_fdr: q_protein_group,
            ptm_localization: PtmLocalizationRunStats {
                enabled: self.parameters.ptm_localization.enabled,
                localized_psms,
                psm_q_value: self.parameters.ptm_localization.psm_q_value,
                localization_q_value: self.parameters.ptm_localization.localization_q_value,
            },
            models: ModelRunStats {
                mass_alignment_applied: true,
                retention_time_prediction_enabled: self.parameters.predict_rt,
                retention_time_model_fitted,
                retention_time_features: format!(
                    "{:?}",
                    self.parameters.retention_time_model.features
                )
                .to_lowercase(),
                retention_time_alignment: needs_alignment.then(|| {
                    format!(
                        "{:?}",
                        self.parameters.retention_time_alignment.unwrap_or_default()
                    )
                    .to_lowercase()
                }),
                ion_mobility_observed: has_ion_mobility,
                ion_mobility_model_enabled: self.parameters.ion_mobility_model.enabled,
                ion_mobility_model_fitted,
                ion_mobility_features: format!("{:?}", self.parameters.ion_mobility_model.features)
                    .to_lowercase(),
            },
            quantification: QuantificationRunStats {
                lfq_enabled: self.parameters.quant.lfq,
                lfq_features,
                tmt: self
                    .parameters
                    .quant
                    .tmt
                    .as_ref()
                    .map(|tmt| format!("{tmt:?}").to_lowercase()),
                tmt_features,
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
                static_definitions: self.database_parameters.static_mods.len(),
                variable_definitions: self
                    .database_parameters
                    .variable_mods
                    .values()
                    .map(Vec::len)
                    .sum(),
                max_variable_mods: self.database_parameters.max_variable_mods,
                max_total_variable_mods: self.database_parameters.max_total_variable_mods,
                max_combinations: self.database_parameters.max_combinations,
                ptm_library_sites: self
                    .database_parameters
                    .loaded_ptm_library
                    .as_deref()
                    .map_or(0, |library| library.len()),
                label_channels: self.database.label_channels.len(),
                labeled_peptides: self
                    .database
                    .peptides
                    .iter()
                    .filter(|peptide| !peptide.decoy && peptide.label_channel.is_some())
                    .count(),
            },
            spectral_library: SpectralLibraryRunStats {
                enabled: self.parameters.spectral_library.enabled,
                entries: library_entries.len(),
                transitions: library_transitions,
                strategy: match self.parameters.spectral_library.strategy {
                    SpectralLibraryStrategy::BestPsm => "best_psm".into(),
                    SpectralLibraryStrategy::Consensus => "consensus".into(),
                },
                psm_q_value: self.parameters.spectral_library.psm_q_value,
                peptide_q_value: self.parameters.spectral_library.peptide_q_value,
                formats: self
                    .parameters
                    .spectral_library
                    .formats
                    .iter()
                    .map(|format| match format {
                        SpectralLibraryFormat::SageParquet => "sage_parquet".into(),
                        SpectralLibraryFormat::MzSpecLib => "mzspeclib".into(),
                    })
                    .collect(),
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
            self.database.peptides.len(),
            self.database.fragments.len(),
            run_time,
        );

        self.events.emit(EventKind::JobCompleted {
            runtime_secs: run_time,
            outputs: summary.output_paths.len(),
        });
        self.events.check()?;

        Ok((telemetry, summary))
    }
}
