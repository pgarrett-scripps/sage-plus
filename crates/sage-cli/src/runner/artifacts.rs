use super::*;

impl Runner {
    /// Flatten FDR-passing target PSMs into one [`SiteRow`] per localized
    /// modification site. Shared by the PSM-site and protein-site reports.
    pub(super) fn collect_site_rows(
        &self,
        features: &[Feature],
        filenames: &[String],
    ) -> Vec<SiteRow> {
        let mut rows = Vec::new();
        for feature in features {
            // Only confidently-identified target PSMs.
            if !passes_localization_filter(feature, self.parameters.ptm_localization.psm_q_value) {
                continue;
            }
            let localization = match &feature.localization {
                Some(loc) => loc,
                None => continue,
            };
            let peptide = &self.database[feature.peptide_idx];
            let peptide_str = peptide.to_string();
            let proteins =
                peptide.proteins(&self.database.decoy_tag, self.database.generate_decoys);
            let filename = filenames.get(feature.file_id).cloned().unwrap_or_default();

            for m in &localization.mods {
                if m.decoy_winner
                    || m.localization_q_value
                        > self.parameters.ptm_localization.localization_q_value
                {
                    continue;
                }
                let modification = m.label.clone().unwrap_or_else(|| format!("{:+}", m.mass));
                let site_probabilities = m
                    .all_sites
                    .iter()
                    .map(|s| {
                        format!(
                            "{}{}:{:.4}",
                            s.residue as char,
                            s.position + 1,
                            s.probability
                        )
                    })
                    .collect::<Vec<_>>()
                    .join(";");

                for site in &m.best_sites {
                    rows.push(SiteRow {
                        psm_id: feature.psm_id,
                        filename: filename.clone(),
                        scannr: feature.spec_id.clone(),
                        peptide: peptide_str.clone(),
                        peptide_sequence: String::from_utf8_lossy(&peptide.sequence).into_owned(),
                        proteins: proteins.clone(),
                        charge: feature.charge,
                        spectrum_q: feature.spectrum_q,
                        peptide_q: feature.peptide_q,
                        modification: modification.clone(),
                        modification_mass: m.mass,
                        position: site.position + 1,
                        residue: site.residue,
                        localization_probability: site.probability,
                        delta_score: m.delta_score,
                        target_decoy_score: m.target_decoy_score,
                        localization_q_value: m.localization_q_value,
                        candidate_sites: m.candidate_sites,
                        site_determining_matched: m.site_determining_matched,
                        site_determining_total: m.site_determining_ions,
                        site_probabilities: site_probabilities.clone(),
                    });
                }
            }
        }
        rows
    }

    /// Write a per-PSM-site PTM localization report (one row per localized
    /// modification site of each FDR-passing PSM).
    pub fn write_ptm_sites(
        &self,
        features: &[Feature],
        filenames: &[String],
    ) -> anyhow::Result<Url> {
        let rows = self.collect_site_rows(features, filenames);

        use sage_cloudpath::parquet::PtmSiteRecord;
        let records = rows
            .iter()
            .map(|row| PtmSiteRecord {
                psm_id: row.psm_id as i64,
                filename: row.filename.clone(),
                scannr: row.scannr.clone(),
                peptide: row.peptide.clone(),
                proteins: row.proteins.clone(),
                charge: row.charge as i32,
                spectrum_q: row.spectrum_q,
                peptide_q: row.peptide_q,
                modification: row.modification.clone(),
                modification_mass: row.modification_mass,
                position: row.position as i32,
                residue: (row.residue as char).to_string(),
                localization_probability: row.localization_probability,
                delta_localization_score: row.delta_score,
                target_decoy_score: row.target_decoy_score,
                localization_q_value: row.localization_q_value,
                candidate_sites: row.candidate_sites as i32,
                site_determining_ions_matched: row.site_determining_matched as i32,
                site_determining_ions_total: row.site_determining_total as i32,
                site_probabilities: row.site_probabilities.clone(),
            })
            .collect::<Vec<_>>();
        let path = self.make_path("results.sage.ptm-sites.parquet");
        let bytes = sage_cloudpath::parquet::serialize_ptm_sites(&records)?;
        sage_cloudpath::write_bytes_sync(&path, bytes)?;
        Ok(path)
    }

    /// Write a collapsed protein-site report: the best localization for each
    /// (protein, modified peptide site) aggregated across all supporting PSMs.
    pub fn write_protein_sites(
        &self,
        features: &[Feature],
        filenames: &[String],
    ) -> anyhow::Result<Url> {
        let rows = self.collect_site_rows(features, filenames);

        // Key on (protein, peptide, position, mod mass). Protein coordinates
        // are not resolved (the FASTA is consumed during indexing), so a row
        // represents a localized site on a peptide, attributed to each protein
        // the peptide maps to.
        #[derive(Clone)]
        struct Agg {
            protein: String,
            peptide: String,
            residue: u8,
            position: usize,
            modification: String,
            modification_mass: f32,
            n_psms: u32,
            best_probability: f32,
            best_delta_score: f32,
            best_localization_q_value: f32,
            best_spectrum_q: f32,
        }

        let mut map: HashMap<(String, String, usize, i64), Agg> = HashMap::new();
        for row in &rows {
            for protein in row.proteins.split(';').filter(|p| !p.is_empty()) {
                let mass_key = (row.modification_mass * 1e3).round() as i64;
                let key = (
                    protein.to_string(),
                    row.peptide.clone(),
                    row.position,
                    mass_key,
                );
                let entry = map.entry(key).or_insert_with(|| Agg {
                    protein: protein.to_string(),
                    peptide: row.peptide.clone(),
                    residue: row.residue,
                    position: row.position,
                    modification: row.modification.clone(),
                    modification_mass: row.modification_mass,
                    n_psms: 0,
                    best_probability: 0.0,
                    best_delta_score: f32::MIN,
                    best_localization_q_value: 1.0,
                    best_spectrum_q: f32::MAX,
                });
                entry.n_psms += 1;
                entry.best_probability = entry.best_probability.max(row.localization_probability);
                entry.best_delta_score = entry.best_delta_score.max(row.delta_score);
                entry.best_localization_q_value = entry
                    .best_localization_q_value
                    .min(row.localization_q_value);
                entry.best_spectrum_q = entry.best_spectrum_q.min(row.spectrum_q);
            }
        }

        let mut aggregated: Vec<Agg> = map.into_values().collect();
        aggregated.sort_by(|a, b| {
            a.protein
                .cmp(&b.protein)
                .then_with(|| a.peptide.cmp(&b.peptide))
                .then_with(|| a.position.cmp(&b.position))
        });

        use sage_cloudpath::parquet::ProteinSiteRecord;
        let records = aggregated
            .iter()
            .map(|agg| ProteinSiteRecord {
                protein: agg.protein.clone(),
                peptide: agg.peptide.clone(),
                residue: (agg.residue as char).to_string(),
                position_in_peptide: agg.position as i32,
                modification: agg.modification.clone(),
                modification_mass: agg.modification_mass,
                num_psms: agg.n_psms as i32,
                best_localization_probability: agg.best_probability,
                best_delta_localization_score: agg.best_delta_score,
                best_localization_q_value: agg.best_localization_q_value,
                best_spectrum_q: agg.best_spectrum_q,
            })
            .collect::<Vec<_>>();
        let path = self.make_path("results.sage.protein-sites.parquet");
        let bytes = sage_cloudpath::parquet::serialize_protein_sites(&records)?;
        sage_cloudpath::write_bytes_sync(&path, bytes)?;
        Ok(path)
    }

    /// Emit a compact, reusable protein-coordinate site library from passing
    /// localized PSMs. Only names defined by this search's `variable_mods` are
    /// included, so every emitted row can be resolved by the same config.
    pub(super) fn write_ptm_library(
        &self,
        features: &[Feature],
        filenames: &[String],
    ) -> anyhow::Result<Vec<Url>> {
        if self.database_parameters.fasta.is_empty() {
            return Ok(Vec::new());
        }

        let known_names = self
            .database_parameters
            .variable_mods
            .values()
            .flatten()
            .filter_map(|entry| entry.definition().name.map(|name| name.to_string()))
            .collect::<HashSet<_>>();
        let fasta_url = sage_cloudpath::to_url(&self.database_parameters.fasta)?;
        let fasta = sage_cloudpath::util::read_fasta(
            &fasta_url,
            &self.database_parameters.decoy_tag,
            self.database_parameters.generate_decoys,
        )?;
        let proteins = fasta
            .targets
            .iter()
            .map(|(accession, sequence)| (accession.as_ref(), sequence.as_str()))
            .collect::<HashMap<_, _>>();

        let mut sites = HashSet::new();
        let mut skipped_unnamed = 0usize;
        for row in self.collect_site_rows(features, filenames) {
            if !known_names.contains(&row.modification) {
                skipped_unnamed += 1;
                continue;
            }
            let peptide_position = row.position - 1;
            for protein in row
                .proteins
                .split(';')
                .filter(|protein| !protein.is_empty())
            {
                let Some(sequence) = proteins.get(protein) else {
                    continue;
                };
                for (start, _) in sequence.match_indices(&row.peptide_sequence) {
                    let protein_position = start + peptide_position;
                    if sequence.as_bytes().get(protein_position) == Some(&row.residue) {
                        sites.insert(sage_core::ptm_library::PtmLibrarySite {
                            protein: Arc::from(protein),
                            position: protein_position as u32,
                            residue: row.residue,
                            modification: Arc::from(row.modification.as_str()),
                        });
                    }
                }
            }
        }
        if skipped_unnamed > 0 {
            log::warn!(
                "skipped {} localized sites without a matching configured modification name",
                skipped_unnamed
            );
        }

        let mut sites = sites.into_iter().collect::<Vec<_>>();
        sites.sort_unstable_by(|a, b| {
            a.protein
                .cmp(&b.protein)
                .then_with(|| a.position.cmp(&b.position))
                .then_with(|| a.modification.cmp(&b.modification))
        });
        let parquet_path = self.make_path("results.sage.ptm-library.parquet");
        let bytes = sage_cloudpath::parquet::serialize_ptm_library(&sites)?;
        sage_cloudpath::write_bytes_sync(&parquet_path, bytes)?;

        let tsv_path = self.make_path("results.sage.ptm-library.tsv");
        let mut writer = csv::WriterBuilder::new()
            .delimiter(b'\t')
            .from_writer(Vec::new());
        writer.write_record(["protein", "position", "residue", "modification"])?;
        for site in &sites {
            writer.write_record([
                site.protein.as_ref(),
                &(site.position + 1).to_string(),
                std::str::from_utf8(&[site.residue])?,
                site.modification.as_ref(),
            ])?;
        }
        writer.flush()?;
        sage_cloudpath::write_bytes_sync(&tsv_path, writer.into_inner()?)?;
        Ok(vec![parquet_path, tsv_path])
    }

    pub(super) fn serialize_pin(
        &self,
        re: &regex::Regex,
        feature: &Feature,
        filenames: &[String],
    ) -> csv::ByteRecord {
        let scannr = re
            .captures_iter(&feature.spec_id)
            .last()
            .and_then(|cap| cap.get(1).map(|cap| cap.as_str()))
            .unwrap_or(&feature.spec_id);

        let mut record = csv::ByteRecord::new();
        let peptide = &self.database[feature.peptide_idx];
        record.push_field(itoa::Buffer::new().format(feature.psm_id).as_bytes());
        record.push_field(itoa::Buffer::new().format(feature.label).as_bytes());
        record.push_field(scannr.as_bytes());
        record.push_field(ryu::Buffer::new().format(feature.expmass).as_bytes());
        record.push_field(ryu::Buffer::new().format(feature.calcmass).as_bytes());
        record.push_field(filenames[feature.file_id].as_bytes());
        record.push_field(ryu::Buffer::new().format(feature.rt).as_bytes());
        record.push_field(ryu::Buffer::new().format(feature.ims).as_bytes());
        record.push_field(itoa::Buffer::new().format(feature.rank).as_bytes());
        record.push_field(
            itoa::Buffer::new()
                .format((feature.charge == 2) as i32)
                .as_bytes(),
        );
        record.push_field(
            itoa::Buffer::new()
                .format((feature.charge == 3) as i32)
                .as_bytes(),
        );
        record.push_field(
            itoa::Buffer::new()
                .format((feature.charge == 4) as i32)
                .as_bytes(),
        );
        record.push_field(
            itoa::Buffer::new()
                .format((feature.charge == 5) as i32)
                .as_bytes(),
        );
        record.push_field(
            itoa::Buffer::new()
                .format((feature.charge == 6) as i32)
                .as_bytes(),
        );
        record.push_field(
            itoa::Buffer::new()
                .format(if feature.charge < 2 || feature.charge > 6 {
                    feature.charge
                } else {
                    0
                })
                .as_bytes(),
        );
        record.push_field(itoa::Buffer::new().format(feature.peptide_len).as_bytes());
        record.push_field(
            itoa::Buffer::new()
                .format(feature.missed_cleavages)
                .as_bytes(),
        );
        record.push_field(
            itoa::Buffer::new()
                .format(peptide.semi_enzymatic as u8)
                .as_bytes(),
        );
        record.push_field(ryu::Buffer::new().format(feature.isotope_error).as_bytes());
        record.push_field(
            ryu::Buffer::new()
                .format(feature.delta_mass.abs().ln_1p())
                .as_bytes(),
        );
        record.push_field(ryu::Buffer::new().format(feature.average_ppm).as_bytes());
        record.push_field(
            ryu::Buffer::new()
                .format(feature.hyperscore.ln_1p())
                .as_bytes(),
        );
        record.push_field(
            ryu::Buffer::new()
                .format(feature.delta_next.ln_1p())
                .as_bytes(),
        );
        record.push_field(
            ryu::Buffer::new()
                .format(feature.delta_best.ln_1p())
                .as_bytes(),
        );
        record.push_field(ryu::Buffer::new().format(feature.aligned_rt).as_bytes());
        record.push_field(ryu::Buffer::new().format(feature.predicted_rt).as_bytes());
        record.push_field(
            ryu::Buffer::new()
                .format(feature.delta_rt_model.clamp(0.001, 1.0).sqrt())
                .as_bytes(),
        );
        record.push_field(ryu::Buffer::new().format(feature.predicted_ims).as_bytes());
        record.push_field(
            ryu::Buffer::new()
                .format(feature.delta_ims_model)
                .as_bytes(),
        );
        record.push_field(itoa::Buffer::new().format(feature.matched_peaks).as_bytes());
        record.push_field(itoa::Buffer::new().format(feature.longest_b).as_bytes());
        record.push_field(itoa::Buffer::new().format(feature.longest_y).as_bytes());
        record.push_field(ryu::Buffer::new().format(feature.longest_y_pct).as_bytes());
        record.push_field(
            ryu::Buffer::new()
                .format(feature.matched_intensity_pct.ln_1p())
                .as_bytes(),
        );
        record.push_field(
            itoa::Buffer::new()
                .format(feature.scored_candidates)
                .as_bytes(),
        );
        record.push_field(
            ryu::Buffer::new()
                .format((-feature.poisson).ln_1p())
                .as_bytes(),
        );
        record.push_field(
            ryu::Buffer::new()
                .format(feature.posterior_error)
                .as_bytes(),
        );
        record.push_field(peptide.to_string().as_bytes());
        record.push_field(
            peptide
                .proteins(&self.database.decoy_tag, self.database.generate_decoys)
                .as_bytes(),
        );
        record
    }

    pub fn write_pin(&self, features: &[Feature], filenames: &[String]) -> anyhow::Result<Url> {
        let path = self.make_path("results.sage.pin");

        let mut wtr = csv::WriterBuilder::new()
            .delimiter(b'\t')
            .from_writer(OutputTarget::new(&path)?);

        let headers = csv::ByteRecord::from(vec![
            "SpecId",
            "Label",
            "ScanNr",
            "ExpMass",
            "CalcMass",
            "FileName",
            "retentiontime",
            "ion_mobility",
            "rank",
            "z=2",
            "z=3",
            "z=4",
            "z=5",
            "z=6",
            "z=other",
            "peptide_len",
            "missed_cleavages",
            "semi_enzymatic",
            "isotope_error",
            "ln(precursor_ppm)",
            "fragment_ppm",
            "ln(hyperscore)",
            "ln(delta_next)",
            "ln(delta_best)",
            "aligned_rt",
            "predicted_rt",
            "sqrt(delta_rt_model)",
            "predicted_mobility",
            "sqrt(delta_mobility)",
            "matched_peaks",
            "longest_b",
            "longest_y",
            "longest_y_pct",
            "ln(matched_intensity_pct)",
            "scored_candidates",
            "ln(-poisson)",
            "posterior_error",
            "Peptide",
            "Proteins",
        ]);

        let re = regex::Regex::new(r"scan=(\d+)").expect("This is valid regex");

        wtr.write_byte_record(&headers)?;
        for chunk in features.chunks(1024) {
            for record in chunk
                .par_iter()
                .map(|feat| self.serialize_pin(&re, feat, filenames))
                .collect::<Vec<_>>()
            {
                wtr.write_byte_record(&record)?;
            }
        }

        finish_csv_writer(wtr, &path)?;
        Ok(path)
    }

    pub(super) fn write_report(
        &self,
        features: &[Feature],
        areas: Option<HashMap<(PrecursorId, bool), QuantifiedPeak, fnv::FnvBuildHasher>>,
        filenames: &[String],
    ) -> anyhow::Result<Url> {
        let path = self.make_path("results.sage.report.html");

        let global_q_value_filter = 0.01;
        let predict_section_q_value_filter = 0.01;

        // Create a new report
        let mut report = Report::new(
            "Sage",
            &self.parameters.version,
            Some(
                "https://github.com/pgarrett-scripps/sage-plus/blob/main/figures/logo.png?raw=true",
            ),
            "Sage Report",
        );

        /* Section 1: Introduction */
        {
            let mut intro_section = ReportSection::new("Results Overview");
            intro_section.add_content(html! {
                "The following files were processed:"
                ul {
                    @for filename in filenames {
                        li { (filename) }
                    }
                }
            });

            // Number of targets identified at global q-value filter at spectrum level per file
            let num_psm_targets_per_file: Vec<usize> = filenames
                .iter()
                .map(|filename| {
                    features
                        .iter()
                        .filter(|f| {
                            f.label == 1
                                && f.spectrum_q <= global_q_value_filter
                                && filenames[f.file_id] == *filename
                        })
                        .count()
                })
                .collect();

            // Number of peptides identified at global q-value filter at peptide level per file
            let mut num_peptide_targets_per_file: Vec<usize> = Vec::new();
            for filename in filenames {
                let mut peptides = HashSet::new();
                for feature in features.iter().filter(|f| {
                    f.label == 1
                        && f.peptide_q <= global_q_value_filter
                        && filenames[f.file_id] == *filename
                }) {
                    peptides.insert(self.database[feature.peptide_idx].to_string());
                }
                num_peptide_targets_per_file.push(peptides.len());
            }

            // Number of proteins identified at global q-value filter at protein level per file
            let mut num_protein_targets_per_file: Vec<usize> = Vec::new();
            for filename in filenames {
                let mut proteins = HashSet::new();
                for feature in features.iter().filter(|f| {
                    f.label == 1
                        && f.protein_q <= global_q_value_filter
                        && filenames[f.file_id] == *filename
                }) {
                    proteins.insert(
                        self.database[feature.peptide_idx]
                            .proteins(&self.database.decoy_tag, self.database.generate_decoys),
                    );
                }
                num_protein_targets_per_file.push(proteins.len());
            }

            // Total MS2 intensity at global q-value filter at each level per file
            let total_ms2_intensity_per_file: Vec<f32> = filenames
                .iter()
                .map(|filename| {
                    features
                        .iter()
                        .filter(|f| {
                            f.label == 1
                                && f.spectrum_q <= global_q_value_filter
                                && f.peptide_q <= global_q_value_filter
                                && f.protein_q <= global_q_value_filter
                                && filenames[f.file_id] == *filename
                        })
                        .map(|f| f.ms2_intensity)
                        .sum()
                })
                .collect();

            // Total LFQ (MS1) intensity at global q-value filter per file (if LFQ is enabled)
            let total_lfq_intensity_per_file: Vec<f32> = if let Some(areas) = &areas {
                let mut total_lfq_intensities = Vec::new();
                for i in 0..filenames.len() {
                    let mut intensities = Vec::new();
                    for ((_id, decoy), quantified) in areas {
                        if !decoy && quantified.peak.q_value <= global_q_value_filter {
                            if let Some(intensity) = quantified.intensities[i] {
                                intensities.push(intensity as f32);
                            }
                        }
                    }
                    total_lfq_intensities.push(intensities.iter().sum());
                }
                total_lfq_intensities
            } else {
                vec![0.0; filenames.len()]
            };

            // Mmedian MS1 mass accuracy for each file, using feature.delta_mass
            let median_ms1_mass_accuracy_per_file: Vec<f32> = filenames
                .iter()
                .map(|filename| {
                    median_finite(features.iter().filter_map(|feature| {
                        (filenames[feature.file_id] == *filename
                            && feature.label == 1
                            && feature.spectrum_q <= global_q_value_filter)
                            .then_some(feature.delta_mass)
                    }))
                    .unwrap_or(f32::NAN)
                })
                .collect();

            // Median MS2 mass accuracy for each file, using feature.average_ppm
            let median_ms2_mass_accuracy_per_file: Vec<f32> = filenames
                .iter()
                .map(|filename| {
                    median_finite(features.iter().filter_map(|feature| {
                        (filenames[feature.file_id] == *filename
                            && feature.label == 1
                            && feature.spectrum_q <= global_q_value_filter)
                            .then_some(feature.average_ppm)
                    }))
                    .unwrap_or(f32::NAN)
                })
                .collect();

            // Median RT deviation for each file, using feature.delta_rt_model
            let median_rt_deviation_per_file: Vec<f32> = filenames
                .iter()
                .map(|filename| {
                    median_finite(features.iter().filter_map(|feature| {
                        (filenames[feature.file_id] == *filename
                            && feature.label == 1
                            && feature.spectrum_q <= global_q_value_filter)
                            .then_some(feature.delta_rt_model)
                    }))
                    .unwrap_or(f32::NAN)
                })
                .collect();

            // Median IM deviation for each file, using feature.delta_ims_model
            let median_im_deviation_per_file: Vec<f32> = filenames
                .iter()
                .map(|filename| {
                    median_finite(features.iter().filter_map(|feature| {
                        (filenames[feature.file_id] == *filename
                            && feature.label == 1
                            && feature.spectrum_q <= global_q_value_filter)
                            .then_some(feature.delta_ims_model)
                    }))
                    .unwrap_or(f32::NAN)
                })
                .collect();

            // Average peptide length for each file
            let avg_peptide_length_per_file: Vec<f32> = filenames
                .iter()
                .map(|filename| {
                    average_finite(features.iter().filter_map(|feature| {
                        (filenames[feature.file_id] == *filename
                            && feature.label == 1
                            && feature.spectrum_q <= global_q_value_filter)
                            .then_some(feature.peptide_len as f32)
                    }))
                    .unwrap_or(f32::NAN)
                })
                .collect();

            // Average peptide charge for each file
            let avg_peptide_charge_per_file: Vec<f32> = filenames
                .iter()
                .map(|filename| {
                    average_finite(features.iter().filter_map(|feature| {
                        (filenames[feature.file_id] == *filename
                            && feature.label == 1
                            && feature.spectrum_q <= global_q_value_filter)
                            .then_some(feature.charge as f32)
                    }))
                    .unwrap_or(f32::NAN)
                })
                .collect();

            // Average number of matched peaks for each file
            let avg_matched_peaks_per_file: Vec<f32> = filenames
                .iter()
                .map(|filename| {
                    average_finite(features.iter().filter_map(|feature| {
                        (filenames[feature.file_id] == *filename
                            && feature.label == 1
                            && feature.spectrum_q <= global_q_value_filter)
                            .then_some(feature.matched_peaks as f32)
                    }))
                    .unwrap_or(f32::NAN)
                })
                .collect();

            // Prepare html table to add to the report
            let table = html! {
                div class="table-container" {
                    table id="dataTable"  class="display" {
                        thead {
                            tr {
                                th { "File" }
                                th { "PSMs" }
                                th { "Peptides" }
                                th { "Proteins" }
                                th { "Total MS1 Intensity" }
                                th { "Total MS2 Intensity" }
                                th { "Median MS1 Delta Mass" }
                                th { "Median MS2 Delta Mass" }
                                th { "Median RT Deviation" }
                                th { "Median IM Deviation" }
                                th { "Average Peptide Length" }
                                th { "Average Peptide Charge" }
                                th { "Average Matched Peaks" }
                            }
                        }
                        tbody {
                            @for (i, filename) in filenames.iter().enumerate() {
                                tr {
                                    td { (filename) }
                                    td { (num_psm_targets_per_file[i]) }
                                    td { (num_peptide_targets_per_file[i]) }
                                    td { (num_protein_targets_per_file[i]) }
                                    td { (total_lfq_intensity_per_file[i]) }
                                    td { (total_ms2_intensity_per_file[i]) }
                                    td { (median_ms1_mass_accuracy_per_file[i]) }
                                    td { (median_ms2_mass_accuracy_per_file[i]) }
                                    td { (median_rt_deviation_per_file[i]) }
                                    td { (median_im_deviation_per_file[i]) }
                                    td { (avg_peptide_length_per_file[i]) }
                                    td { (avg_peptide_charge_per_file[i]) }
                                    td { (avg_matched_peaks_per_file[i]) }
                                }
                            }
                        }
                    }
                    button id="downloadCsv" { "Download as CSV" }
                }
            };

            intro_section.add_content(table);

            // Add boxplot of the LFQ intensities from areas if available
            if let Some(areas) = areas {
                let mut lfq_intensities: Vec<Vec<f64>> = Vec::new();
                for i in 0..filenames.len() {
                    let mut intensities = Vec::new();
                    for ((_id, decoy), quantified) in &areas {
                        if !decoy && quantified.peak.q_value <= global_q_value_filter {
                            if let Some(intensity) = quantified.intensities[i] {
                                if intensity.is_finite() && intensity > 0.0 {
                                    intensities.push(intensity.log2());
                                }
                            }
                        }
                    }
                    lfq_intensities.push(intensities);
                }

                match plot_boxplot(
                    &lfq_intensities,
                    filenames.to_vec(),
                    &format!("LFQ Intensities ({:?}% Q-value)", global_q_value_filter),
                    "Run",
                    "Log2(Intensity)",
                ) {
                    Ok(lfq_boxplot) => intro_section.add_plot(lfq_boxplot),
                    Err(error) => log::warn!("skipping LFQ report plot: {error}"),
                }
            }

            report.add_section(intro_section);
        }

        /* Section 2: Scoring QC */
        {
            let mut scoring_section = ReportSection::new("Scoring Quality Control");

            scoring_section.add_content(html! {
                "It is important to assess the quality of the scoring model to ensure that the model is performing as expected, and that we're not overfitting or violating any assumptions of the Target-Decoy approach. The plot below shows the distribution of discriminant scores for each PSM, colored by whether the PSM is a target or decoy. We would expect the target distributions to be bimodal, where the first mode represents false targets that should align with the decoy distribution, and the second mode represents true targets."
            });

            // Extract sage_discriminant_score and label from features
            let (scores, labels) =
                labeled_finite_values(features, |feature| feature.discriminant_score as f64);
            let has_targets = labels.contains(&1);
            let has_decoys = labels.contains(&-1);

            if scores.len() > 100 && has_targets && has_decoys {
                match plot_score_histogram(&scores, &labels, "LDA Score", "Score") {
                    Ok(score_histogram) => scoring_section.add_plot(score_histogram),
                    Err(error) => log::warn!("skipping LDA score report plot: {error}"),
                }

                match plot_pp(&scores, &labels, "PP Plot") {
                    Ok(pp_plot) => {
                        scoring_section.add_content(html! {
                            "The Probability-Probability (PP) plot is a diagnostic tool that can be used to assess the quality of the scoring model. It plots the empirical cumulative distribution function (ECDF) of the target distribution against the ECDF of the decoy distribution. See: Debrie, E. et. al. (2023) Journal of Proteome Research. for more information."
                        });
                        scoring_section.add_plot(pp_plot);
                    }
                    Err(error) => log::warn!("skipping PP report plot: {error}"),
                }

                for (title, values) in [
                    (
                        "Spectrum Q-value",
                        labeled_finite_values(features, |feature| feature.spectrum_q as f64),
                    ),
                    (
                        "Peptide Q-value",
                        labeled_finite_values(features, |feature| feature.peptide_q as f64),
                    ),
                    (
                        "Protein Q-value",
                        labeled_finite_values(features, |feature| feature.protein_q as f64),
                    ),
                ] {
                    match plot_score_histogram(&values.0, &values.1, title, "Q-value") {
                        Ok(histogram) => scoring_section.add_plot(histogram),
                        Err(error) => log::warn!("skipping {title} report plot: {error}"),
                    }
                }
            } else {
                scoring_section.add_content(html! {
                    div style="margin-top: 10px; margin-bottom: 10px; padding: 15px; background-color: #ffe6e6; border: 1px solid #ff9999; color: #cc0000; border-radius: 5px; white-space: pre-line;" {
                        p {
                            "Scoring quality control plots require more than 100 finite scores with both target and decoy observations."
                        }
                    }
                });
            }

            report.add_section(scoring_section);
        }

        /* Section 3: Predicted Properties */
        {
            let mut predicted_properties_section = ReportSection::new("Predicted Properties");

            predicted_properties_section.add_content(html! {
                "The following plots show the predicted properties of target peptides. The plots show the predicted retention time and ion mobility if present. The predicted properties are used to assess the quality of the model and to identify potential outliers."
            });

            // Normalized experimental RT per file
            let mut rt_per_file: Vec<Vec<f64>> = Vec::new();
            let mut predicted_rt_per_file: Vec<Vec<f64>> = Vec::new();
            for i in 0..filenames.len() {
                let (rts, predicted_rts): (Vec<_>, Vec<_>) = features
                    .iter()
                    .filter(|feature| {
                        feature.label == 1
                            && feature.spectrum_q <= predict_section_q_value_filter
                            && filenames[feature.file_id] == filenames[i]
                            && feature.rt.is_finite()
                            && feature.predicted_rt.is_finite()
                    })
                    .map(|feature| (feature.rt as f64, feature.predicted_rt as f64))
                    .unzip();
                rt_per_file.push(normalize_finite(rts));
                predicted_rt_per_file.push(predicted_rts);
            }

            if !filenames.is_empty() && rt_per_file.iter().any(|values| !values.is_empty()) {
                match plot_scatter(
                    &rt_per_file,
                    &predicted_rt_per_file,
                    filenames.to_vec(),
                    "Retention Time LR Model",
                    "Retention Time",
                    "Predicted Retention Time",
                ) {
                    Ok(rt_scatter) => predicted_properties_section.add_plot(rt_scatter),
                    Err(error) => log::warn!("skipping retention-time report plot: {error}"),
                }
            }

            // Experimental IMS per file
            let mut ims_per_file: Vec<Vec<f64>> = Vec::new();
            let mut predicted_ims_per_file: Vec<Vec<f64>> = Vec::new();
            for i in 0..filenames.len() {
                let (imss, predicted_imss): (Vec<_>, Vec<_>) = features
                    .iter()
                    .filter(|feature| {
                        feature.label == 1
                            && feature.spectrum_q <= predict_section_q_value_filter
                            && filenames[feature.file_id] == filenames[i]
                            && feature.ims.is_finite()
                            && feature.predicted_ims.is_finite()
                    })
                    .map(|feature| (feature.ims as f64, feature.predicted_ims as f64))
                    .unzip();
                ims_per_file.push(imss);
                predicted_ims_per_file.push(predicted_imss);
            }

            if !filenames.is_empty() && ims_per_file.iter().any(|values| !values.is_empty()) {
                match plot_scatter(
                    &ims_per_file,
                    &predicted_ims_per_file,
                    filenames.to_vec(),
                    "Ion Mobility LR Model",
                    "Ion Mobility",
                    "Predicted Ion Mobility",
                ) {
                    Ok(ims_scatter) => predicted_properties_section.add_plot(ims_scatter),
                    Err(error) => log::warn!("skipping ion-mobility report plot: {error}"),
                }
            }

            report.add_section(predicted_properties_section);
        }

        /* Section 4: Configuration */
        {
            let mut config_section = ReportSection::new("Configuration");
            config_section.add_content(html! {
                style {
                    ".code-container {
                        background-color: #f5f5f5;
                        padding: 10px;
                        border-radius: 5px;
                        overflow-x: auto;
                        font-family: monospace;
                        white-space: pre-wrap;
                    }"
                }
                div class="code-container" {
                    pre {
                        code { (PreEscaped(serde_json::to_string_pretty(&self.parameters)?)) }
                    }
                }
            });
            report.add_section(config_section);
        }

        let bytes = report.to_string().into_bytes();
        sage_cloudpath::write_bytes_sync(&path, bytes)?;

        Ok(path)
    }
}
