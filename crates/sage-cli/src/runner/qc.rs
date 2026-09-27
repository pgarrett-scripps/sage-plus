//! Run-level quality-control summaries written next to the search results.

use super::*;
use sage_core::diagnostic::{find_ions, DiagnosticIon};
use sage_core::digestion::DigestionSummary;
use sage_core::enzyme::EnzymeParameters;
use sage_core::polymer::{PolymerFileStats, PolymerScanner};
use sage_core::spectrum::{RawSpectrum, Representation};

/// Spectrum- and peptide-level q-value a PSM must pass to count towards the
/// digestion summary, matching the run summary's 1% FDR counts.
const DIGESTION_Q_VALUE: f32 = 0.01;

/// Percent of MS1 TIC on one polymer above which a run warning is raised.
const POLYMER_WARNING_TIC_PCT: f64 = 5.0;

/// Quality-control scans of one file's raw spectra.
#[derive(Debug, Clone, Default)]
pub struct FileQc {
    /// `None` when the file's MS1 spectra were not read.
    pub polymers: Option<PolymerFileStats>,
    /// `None` when diagnostic ions are off.
    pub diagnostic_ions: Option<DiagnosticFileScan>,
    /// Time spent scanning the raw spectra.
    pub scan_ms: u128,
}

#[derive(Debug, Clone, Default, serde::Serialize, serde::Deserialize)]
pub struct QcRunStats {
    #[serde(default)]
    pub digestion: DigestionRunStats,
    /// Polymer contamination of each file whose MS1 spectra were read.
    #[serde(default)]
    pub polymers: Vec<PolymerRunFileStats>,
    /// Present when `diagnostic_ions` is enabled.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub diagnostic_ions: Option<DiagnosticIonRunStats>,
}

/// Diagnostic ions found in one file's MS2 spectra.
#[derive(Debug, Clone, Default)]
pub struct DiagnosticFileScan {
    pub ms2_spectra: usize,
    /// Spectra containing each ion, parallel to the configured ions.
    pub spectra_with_ion: Vec<usize>,
    /// `(spectrum id, ion index, observed m/z, relative intensity)`, in
    /// file order.
    pub hits: Vec<(String, usize, f32, f32)>,
}

/// Share of MS2 spectra containing each diagnostic ion.
#[derive(Debug, Clone, Default, serde::Serialize, serde::Deserialize)]
pub struct DiagnosticIonRunStats {
    pub ms2_spectra: usize,
    pub ions: Vec<DiagnosticIonShare>,
    pub files: Vec<DiagnosticIonFileStats>,
}

#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct DiagnosticIonFileStats {
    pub file: String,
    pub ms2_spectra: usize,
    pub ions: Vec<DiagnosticIonShare>,
}

#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct DiagnosticIonShare {
    pub name: String,
    pub mz: f32,
    /// MS2 spectra containing the ion.
    pub spectra: usize,
    /// Percent of MS2 spectra containing the ion.
    pub spectra_pct: f64,
}

fn ion_shares(
    ions: &[DiagnosticIon],
    ms2_spectra: usize,
    counts: &[usize],
) -> Vec<DiagnosticIonShare> {
    ions.iter()
        .zip(counts)
        .map(|(ion, &spectra)| DiagnosticIonShare {
            name: ion.name.clone(),
            mz: ion.mz,
            spectra,
            spectra_pct: if ms2_spectra == 0 {
                0.0
            } else {
                100.0 * spectra as f64 / ms2_spectra as f64
            },
        })
        .collect()
}

/// Diagnostic ions in a file's raw MS2 spectra.
fn scan_diagnostic_ions(ions: &[DiagnosticIon], spectra: &[RawSpectrum]) -> DiagnosticFileScan {
    let per_spectrum = spectra
        .par_iter()
        .filter(|spectrum| spectrum.ms_level >= 2)
        .map(|spectrum| {
            let (mz, intensity) = sorted_peaks(spectrum);
            (spectrum.id.as_str(), find_ions(ions, &mz, &intensity))
        })
        .collect::<Vec<_>>();
    let mut scan = DiagnosticFileScan {
        ms2_spectra: per_spectrum.len(),
        spectra_with_ion: vec![0; ions.len()],
        hits: Vec::new(),
    };
    for (id, hits) in per_spectrum {
        for hit in hits {
            scan.spectra_with_ion[hit.ion] += 1;
            scan.hits
                .push((id.to_string(), hit.ion, hit.mz, hit.relative_intensity));
        }
    }
    scan
}

#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct PolymerRunFileStats {
    pub file: String,
    #[serde(flatten)]
    pub stats: PolymerFileStats,
}

/// `(mz, intensity)` of a spectrum sorted by m/z, borrowed when already
/// sorted.
fn sorted_peaks(
    spectrum: &RawSpectrum,
) -> (std::borrow::Cow<'_, [f32]>, std::borrow::Cow<'_, [f32]>) {
    use std::borrow::Cow;
    if spectrum.mz.windows(2).all(|pair| pair[0] <= pair[1]) {
        return (
            Cow::Borrowed(&spectrum.mz),
            Cow::Borrowed(&spectrum.intensity),
        );
    }
    let mut order = (0..spectrum.mz.len()).collect::<Vec<_>>();
    order.sort_by(|&left, &right| spectrum.mz[left].total_cmp(&spectrum.mz[right]));
    (
        Cow::Owned(order.iter().map(|&index| spectrum.mz[index]).collect()),
        Cow::Owned(
            order
                .iter()
                .map(|&index| spectrum.intensity[index])
                .collect(),
        ),
    )
}

/// Polymer ladders in a file's centroided MS1 spectra. `None` when the file
/// has no MS1 spectra.
fn scan_polymers(spectra: &[RawSpectrum]) -> Option<PolymerFileStats> {
    let scanner = PolymerScanner::default();
    let empty = PolymerFileStats::new(&scanner.polymers);
    let mut stats = spectra
        .par_iter()
        .filter(|spectrum| spectrum.ms_level == 1)
        .fold(
            || None::<PolymerFileStats>,
            |stats, spectrum| {
                let mut stats = stats.unwrap_or_else(|| empty.clone());
                if spectrum.representation == Representation::Profile {
                    stats.skipped_profile_spectra += 1;
                } else {
                    let (mz, intensity) = sorted_peaks(spectrum);
                    let tic = intensity.iter().map(|&value| value as f64).sum();
                    stats.add(tic, &scanner.scan(&mz, &intensity));
                }
                Some(stats)
            },
        )
        .reduce(
            || None,
            |left, right| match (left, right) {
                (Some(mut left), Some(right)) => {
                    left.merge(&right);
                    Some(left)
                }
                (left, right) => left.or(right),
            },
        )?;
    stats.finish();
    Some(stats)
}

/// Digestion summary of rank-1 PSMs at 1% spectrum and peptide q-value,
/// counted over distinct peptide sequences with decoys subtracted.
#[derive(Debug, Clone, Default, serde::Serialize, serde::Deserialize)]
pub struct DigestionRunStats {
    pub q_value: f32,
    pub files: Vec<DigestionFileStats>,
    pub total: DigestionSummary,
}

#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct DigestionFileStats {
    pub file: String,
    #[serde(flatten)]
    pub summary: DigestionSummary,
}

impl Runner {
    /// Run the quality-control scans on one file's raw spectra, once per
    /// file: rereads of a file already scanned are skipped.
    pub(super) fn collect_file_qc(&self, file_id: usize, spectra: &[RawSpectrum]) {
        if self
            .file_qc
            .lock()
            .expect("file QC lock")
            .contains_key(&file_id)
        {
            return;
        }
        let start = Instant::now();
        let polymers = scan_polymers(spectra);
        let diagnostic_ions = self
            .parameters
            .diagnostic_ions
            .as_ref()
            .map(|ions| scan_diagnostic_ions(ions, spectra));
        let qc = FileQc {
            polymers,
            diagnostic_ions,
            scan_ms: start.elapsed().as_millis(),
        };
        log::debug!("- file {file_id}: QC scan {} ms", qc.scan_ms);
        self.file_qc
            .lock()
            .expect("file QC lock")
            .insert(file_id, qc);
    }

    /// Report polymer contamination per file, warning when one polymer
    /// carries more than 5% of a file's MS1 TIC.
    pub(super) fn polymer_stats(&self, filenames: &[String]) -> Vec<PolymerRunFileStats> {
        let file_qc = self.file_qc.lock().expect("file QC lock");
        let mut files = Vec::new();
        for (file_id, file) in filenames.iter().enumerate() {
            let Some(stats) = file_qc.get(&file_id).and_then(|qc| qc.polymers.clone()) else {
                continue;
            };
            let shares = stats
                .polymers
                .iter()
                .map(|share| format!("{} {:.2}%", share.name, share.tic_pct))
                .collect::<Vec<_>>()
                .join(", ");
            info!(
                "polymers: {file}: {shares} of MS1 TIC ({} centroided MS1 spectra{})",
                stats.ms1_spectra,
                match stats.skipped_profile_spectra {
                    0 => String::new(),
                    n => format!(", {n} profile spectra skipped"),
                }
            );
            for share in &stats.polymers {
                if share.tic_pct > POLYMER_WARNING_TIC_PCT {
                    let message = format!(
                        "{file}: {} ladders carry {:.1}% of the MS1 TIC (above {POLYMER_WARNING_TIC_PCT}%); the sample may be contaminated",
                        share.name, share.tic_pct
                    );
                    warn!("{message}");
                    self.events.emit(EventKind::Warning {
                        code: "polymer_contamination".into(),
                        message,
                    });
                }
            }
            files.push(PolymerRunFileStats {
                file: file.clone(),
                stats,
            });
        }
        let scan_ms = file_qc.values().map(|qc| qc.scan_ms).sum::<u128>();
        info!("- QC spectrum scans: {scan_ms:8} ms (summed over files)");
        files
    }

    /// Summarize missed cleavages and non-enzymatic termini per file and for
    /// the whole run.
    pub(super) fn digestion_stats(
        &self,
        features: &[Feature],
        filenames: &[String],
    ) -> DigestionRunStats {
        let enzyme = EnzymeParameters::from(self.database_parameters.enzyme.clone()).enzyme;
        let clip = self.database_parameters.clip_n_term_met;
        let passing = features
            .iter()
            .filter(|feature| {
                feature.rank == 1
                    && feature.spectrum_q <= DIGESTION_Q_VALUE
                    && feature.peptide_q <= DIGESTION_Q_VALUE
            })
            .collect::<Vec<_>>();
        let files = filenames
            .iter()
            .enumerate()
            .map(|(file_id, file)| DigestionFileStats {
                file: file.clone(),
                summary: sage_core::digestion::summarize(
                    enzyme.as_ref(),
                    clip,
                    passing
                        .iter()
                        .filter(|feature| feature.file_id == file_id)
                        .map(|feature| &self.database[feature.peptide_idx]),
                ),
            })
            .collect();
        let total = sage_core::digestion::summarize(
            enzyme.as_ref(),
            clip,
            passing
                .iter()
                .map(|feature| &self.database[feature.peptide_idx]),
        );
        info!(
            "digestion: {} peptides (decoy-corrected), missed cleavage {:.1}%, semi N-terminal {:.1}%, semi C-terminal {:.1}%, non-enzymatic {:.1}%",
            total.peptides,
            total.missed_cleavage_pct,
            total.semi_n_pct,
            total.semi_c_pct,
            total.non_enzymatic_pct
        );
        DigestionRunStats {
            q_value: DIGESTION_Q_VALUE,
            files,
            total,
        }
    }

    /// Log and summarize the share of MS2 spectra containing each
    /// diagnostic ion; `None` when diagnostic ions are off.
    pub(super) fn diagnostic_ion_stats(
        &self,
        filenames: &[String],
    ) -> Option<DiagnosticIonRunStats> {
        let ions = self.parameters.diagnostic_ions.as_ref()?;
        let file_qc = self.file_qc.lock().expect("file QC lock");
        let mut ms2_spectra = 0;
        let mut counts = vec![0; ions.len()];
        let mut files = Vec::new();
        for (file_id, file) in filenames.iter().enumerate() {
            let Some(scan) = file_qc
                .get(&file_id)
                .and_then(|qc| qc.diagnostic_ions.as_ref())
            else {
                continue;
            };
            ms2_spectra += scan.ms2_spectra;
            for (count, file_count) in counts.iter_mut().zip(&scan.spectra_with_ion) {
                *count += file_count;
            }
            files.push(DiagnosticIonFileStats {
                file: file.clone(),
                ms2_spectra: scan.ms2_spectra,
                ions: ion_shares(ions, scan.ms2_spectra, &scan.spectra_with_ion),
            });
        }
        let ions = ion_shares(ions, ms2_spectra, &counts);
        info!(
            "diagnostic ions in {ms2_spectra} MS2 spectra: {}",
            ions.iter()
                .map(|ion| format!("{} {:.2}%", ion.name, ion.spectra_pct))
                .collect::<Vec<_>>()
                .join(", ")
        );
        Some(DiagnosticIonRunStats {
            ms2_spectra,
            ions,
            files,
        })
    }

    /// Write `diagnostic_ions.tsv`: one row per ion found in a spectrum.
    pub(super) fn write_diagnostic_ions(
        &self,
        filenames: &[String],
    ) -> anyhow::Result<Option<Url>> {
        let Some(ions) = self.parameters.diagnostic_ions.as_ref() else {
            return Ok(None);
        };
        let mut writer = csv::WriterBuilder::new()
            .delimiter(b'\t')
            .from_writer(Vec::new());
        writer.write_record(["file", "scannr", "ion", "mz", "relative_intensity"])?;
        let file_qc = self.file_qc.lock().expect("file QC lock");
        for (file_id, file) in filenames.iter().enumerate() {
            let Some(scan) = file_qc
                .get(&file_id)
                .and_then(|qc| qc.diagnostic_ions.as_ref())
            else {
                continue;
            };
            for (id, ion, mz, relative_intensity) in &scan.hits {
                writer.write_record([
                    file.as_str(),
                    id.as_str(),
                    ions[*ion].name.as_str(),
                    &format!("{mz:.4}"),
                    &format!("{relative_intensity:.6}"),
                ])?;
            }
        }
        writer.flush()?;
        let path = self.make_path("diagnostic_ions.tsv");
        sage_cloudpath::write_bytes_sync(&path, writer.into_inner()?)?;
        Ok(Some(path))
    }

    pub(super) fn write_digestion(&self, stats: &DigestionRunStats) -> anyhow::Result<Url> {
        let path = self.make_path("digestion.tsv");
        sage_cloudpath::write_bytes_sync(&path, serialize_digestion(stats)?)?;
        Ok(path)
    }
}

pub(super) fn serialize_digestion(stats: &DigestionRunStats) -> anyhow::Result<Vec<u8>> {
    let mut writer = csv::WriterBuilder::new()
        .delimiter(b'\t')
        .from_writer(Vec::new());
    writer.write_record([
        "file",
        "target_peptides",
        "decoy_peptides",
        "peptides",
        "missed_cleavages_0",
        "missed_cleavages_1",
        "missed_cleavages_2_plus",
        "semi_n",
        "semi_c",
        "non_enzymatic",
        "missed_cleavage_pct",
        "semi_n_pct",
        "semi_c_pct",
        "non_enzymatic_pct",
    ])?;
    let rows = stats
        .files
        .iter()
        .map(|file| (file.file.as_str(), &file.summary))
        .chain(std::iter::once(("total", &stats.total)));
    for (file, summary) in rows {
        writer.write_record([
            file.to_string(),
            summary.target_peptides.to_string(),
            summary.decoy_peptides.to_string(),
            summary.peptides.to_string(),
            summary.missed_cleavages_0.to_string(),
            summary.missed_cleavages_1.to_string(),
            summary.missed_cleavages_2_plus.to_string(),
            summary.semi_n.to_string(),
            summary.semi_c.to_string(),
            summary.non_enzymatic.to_string(),
            format!("{:.3}", summary.missed_cleavage_pct),
            format!("{:.3}", summary.semi_n_pct),
            format!("{:.3}", summary.semi_c_pct),
            format!("{:.3}", summary.non_enzymatic_pct),
        ])?;
    }
    writer.flush()?;
    Ok(writer.into_inner()?)
}

#[cfg(test)]
mod test {
    use super::*;

    fn ms1(mz: Vec<f32>, intensity: Vec<f32>, representation: Representation) -> RawSpectrum {
        RawSpectrum {
            ms_level: 1,
            representation,
            mz,
            intensity,
            ..Default::default()
        }
    }

    #[test]
    fn polymer_scan_counts_centroided_ms1_only() {
        // Singly protonated PEG n = 10..=15, listed out of order, plus one
        // unrelated peak of the same total intensity.
        let mut mz = (10..=15)
            .rev()
            .map(|n| (18.010_565 + n as f64 * 44.026_215 + 1.007_276) as f32)
            .collect::<Vec<_>>();
        mz.push(300.0);
        let mut intensity = vec![1.0; 6];
        intensity.push(6.0);
        let spectra = vec![
            ms1(mz.clone(), intensity.clone(), Representation::Centroid),
            ms1(mz, intensity, Representation::Profile),
            RawSpectrum {
                ms_level: 2,
                ..Default::default()
            },
        ];
        let stats = scan_polymers(&spectra).expect("MS1 present");
        assert_eq!(stats.ms1_spectra, 1);
        assert_eq!(stats.skipped_profile_spectra, 1);
        let peg = &stats.polymers[0];
        assert_eq!(peg.name, "peg");
        assert!((peg.tic_pct - 50.0).abs() < 1e-9, "{peg:?}");
        assert!(stats.polymers[1..].iter().all(|share| share.tic_pct == 0.0));
    }

    #[test]
    fn polymer_scan_needs_ms1() {
        let spectra = vec![RawSpectrum {
            ms_level: 2,
            ..Default::default()
        }];
        assert!(scan_polymers(&spectra).is_none());
    }
}
