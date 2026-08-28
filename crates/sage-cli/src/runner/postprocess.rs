use super::*;

impl Runner {
    /// Re-read only the MS2 files needed for post-FDR work. Detailed matched
    /// fragments and PTM localization share this pass so enabling both does
    /// not double the input I/O.
    pub(super) fn postprocess_features(
        &self,
        scorer: &Scorer,
        features: &mut [Feature],
        batch_size: usize,
    ) -> anyhow::Result<PostprocessStats> {
        #[derive(Default)]
        struct SpectrumWork {
            /// Every reported PSM for this spectrum, required to replay
            /// rank-ordered peak removal in chimera mode.
            feature_indices: Vec<usize>,
            annotate: bool,
            localization_indices: Vec<usize>,
        }

        let annotate_matches = self.parameters.annotate_matches;
        let localize = self.parameters.ptm_localization.enabled;
        let library_selections = spectral_library::select_psms(
            features,
            &self.database,
            &self.parameters.spectral_library,
        );
        let library_annotation_indices = library_selections
            .iter()
            .flat_map(|selection| selection.feature_indices.iter().copied())
            .collect::<HashSet<_>>();
        let annotate_fragments = annotate_matches || !library_annotation_indices.is_empty();
        if !annotate_fragments && !localize {
            return Ok(PostprocessStats {
                library_selections,
                ..PostprocessStats::default()
            });
        }

        let output_psm_q_value = self.parameters.output_filter.psm_q_value;
        let annotation_indices = features
            .iter()
            .enumerate()
            .filter_map(|(idx, feature)| {
                ((annotate_matches && passes_output_filter(feature, output_psm_q_value))
                    || library_annotation_indices.contains(&idx))
                .then_some(idx)
            })
            .collect::<HashSet<_>>();
        let mut work: HashMap<usize, HashMap<String, SpectrumWork>> = HashMap::new();

        if annotate_fragments {
            for (idx, feature) in features.iter().enumerate() {
                work.entry(feature.file_id)
                    .or_default()
                    .entry(feature.spec_id.clone())
                    .or_default()
                    .feature_indices
                    .push(idx);
            }
            for idx in &annotation_indices {
                let feature = &features[*idx];
                if let Some(spectrum) = work
                    .get_mut(&feature.file_id)
                    .and_then(|file| file.get_mut(feature.spec_id.as_str()))
                {
                    spectrum.annotate = true;
                }
            }
        }

        if localize {
            for (idx, feature) in features.iter().enumerate() {
                if passes_localization_filter(feature, self.parameters.ptm_localization.psm_q_value)
                    && sage_core::ptm::has_localizable_modification(
                        &self.database[feature.peptide_idx],
                        &self.database.potential_mods,
                    )
                {
                    work.entry(feature.file_id)
                        .or_default()
                        .entry(feature.spec_id.clone())
                        .or_default()
                        .localization_indices
                        .push(idx);
                }
            }
        }

        for file in work.values_mut() {
            file.retain(|_, spectrum| {
                spectrum.annotate || !spectrum.localization_indices.is_empty()
            });
            for spectrum in file.values_mut() {
                spectrum
                    .feature_indices
                    .sort_unstable_by_key(|idx| features[*idx].rank);
            }
        }
        work.retain(|_, file| !file.is_empty());

        let expected_annotations = annotation_indices.len();
        if work.is_empty() {
            anyhow::ensure!(
                expected_annotations == 0,
                "internal error: selected PSMs were not scheduled for fragment annotation"
            );
            return Ok(PostprocessStats {
                library_selections,
                ..PostprocessStats::default()
            });
        }

        let start = Instant::now();
        let mut annotations = Vec::with_capacity(expected_annotations);
        let mut localizations = Vec::new();

        for (chunk_idx, chunk) in self.parameters.mzml_paths.chunks(batch_size).enumerate() {
            self.cancellation.check()?;
            let first_file_id = chunk_idx * batch_size;
            if !(first_file_id..first_file_id + chunk.len())
                .any(|file_id| work.contains_key(&file_id))
            {
                continue;
            }

            let spectra = self
                .read_processed_spectra_with_ms1(chunk, chunk_idx, batch_size, false)?
                .1;
            let results = spectra
                .par_iter()
                .filter_map(|spectrum| {
                    let spectrum_work = work.get(&spectrum.file_id)?.get(spectrum.id.as_str())?;
                    let mut annotated = Vec::new();

                    if spectrum_work.annotate {
                        let ranked_features = spectrum_work
                            .feature_indices
                            .iter()
                            .map(|&idx| &features[idx])
                            .collect::<Vec<_>>();
                        let selected = spectrum_work
                            .feature_indices
                            .iter()
                            .map(|idx| annotation_indices.contains(idx))
                            .collect::<Vec<_>>();
                        annotated.extend(
                            spectrum_work
                                .feature_indices
                                .iter()
                                .copied()
                                .zip(scorer.annotate_ranked_candidates(
                                    spectrum,
                                    &ranked_features,
                                    &selected,
                                ))
                                .filter_map(|(idx, fragments)| {
                                    fragments.map(|fragments| (idx, fragments))
                                }),
                        );
                    }

                    let localized = spectrum_work
                        .localization_indices
                        .iter()
                        .copied()
                        .filter_map(|idx| {
                            let feature = &features[idx];
                            let peptide = &self.database[feature.peptide_idx];
                            let localization = sage_core::ptm::localize(
                                peptide,
                                spectrum,
                                &self.database.ion_kinds,
                                &self.database.potential_mods,
                                self.parameters.fragment_tol,
                                self.parameters.max_fragment_charge,
                                feature.charge,
                            );
                            (!localization.mods.is_empty()).then_some((idx, localization))
                        })
                        .collect::<Vec<_>>();

                    Some((annotated, localized))
                })
                .collect::<Vec<_>>();

            for (mut batch_annotations, mut batch_localizations) in results {
                annotations.append(&mut batch_annotations);
                localizations.append(&mut batch_localizations);
            }
        }

        let mut annotated_indices = HashSet::with_capacity(annotations.len());
        let mut annotated_fragments = 0usize;
        for (idx, fragments) in annotations {
            anyhow::ensure!(
                annotated_indices.insert(idx),
                "spectrum {} in file {} was encountered more than once during deferred annotation",
                features[idx].spec_id,
                features[idx].file_id
            );
            annotated_fragments += fragments.fragment_ordinals.len();
            features[idx].fragments = Some(fragments);
        }

        if annotated_indices.len() != expected_annotations {
            let missing = features
                .iter()
                .enumerate()
                .find(|(idx, _feature)| {
                    annotation_indices.contains(idx) && !annotated_indices.contains(idx)
                })
                .map(|(_, feature)| {
                    format!(
                        "file {} spectrum {} (PSM {})",
                        feature.file_id, feature.spec_id, feature.psm_id
                    )
                })
                .unwrap_or_else(|| "unknown PSM".into());
            anyhow::bail!(
                "deferred fragment annotation completed {}/{} selected PSMs; missing {}",
                annotated_indices.len(),
                expected_annotations,
                missing
            );
        }

        let localized_psms = localizations.len();
        for (idx, localization) in localizations {
            features[idx].localization = Some(localization);
        }

        let localization_indices = features
            .iter()
            .enumerate()
            .flat_map(|(feature_idx, feature)| {
                feature.localization.iter().flat_map(move |localization| {
                    localization
                        .mods
                        .iter()
                        .enumerate()
                        .filter(|(_, modification)| modification.competition_eligible)
                        .map(move |(mod_idx, _)| (feature_idx, mod_idx))
                })
            })
            .collect::<Vec<_>>();
        let evidence = localization_indices
            .iter()
            .map(|&(feature_idx, mod_idx)| {
                let modification =
                    &features[feature_idx].localization.as_ref().unwrap().mods[mod_idx];
                (modification.target_decoy_score, modification.decoy_winner)
            })
            .collect::<Vec<_>>();
        let q_values = sage_core::ptm::target_decoy_q_values(&evidence);
        for ((feature_idx, mod_idx), q_value) in localization_indices.into_iter().zip(q_values) {
            features[feature_idx].localization.as_mut().unwrap().mods[mod_idx]
                .localization_q_value = q_value;
        }

        log::info!(
            "- post-FDR MS2 pass: {} annotated PSMs ({} fragments), {} localized PSMs in {} ms",
            annotated_indices.len(),
            annotated_fragments,
            localized_psms,
            start.elapsed().as_millis()
        );

        Ok(PostprocessStats {
            annotated_psms: annotated_indices.len(),
            annotated_fragments,
            localized_psms,
            library_selections,
        })
    }
}
