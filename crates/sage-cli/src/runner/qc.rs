//! Run-level quality-control summaries written next to the search results.

use super::*;
use sage_core::diagnostic::{find_ions, DiagnosticIon};
use sage_core::digestion::DigestionSummary;
use sage_core::enzyme::EnzymeParameters;
use sage_core::polymer::{PolymerFileStats, PolymerScanner};
use sage_core::spectrum::{RawSpectrum, Representation};
use std::collections::HashSet;

/// Spectrum-, peptide- and protein-level q-value a PSM must pass to count
/// towards any QC table, matching the run summary's 1% FDR counts.
const QC_Q_VALUE: f32 = 0.01;

/// Spectrum- and peptide-level q-value a PSM must pass to count towards the
/// digestion summary, matching the run summary's 1% FDR counts.
const DIGESTION_Q_VALUE: f32 = QC_Q_VALUE;

/// Equal-width retention-time bins each file's gradient is split into for
/// `qc.gradient.tsv`.
const GRADIENT_BINS: usize = 10;

/// Percent of MS1 TIC on one polymer above which a run warning is raised.
const POLYMER_WARNING_TIC_PCT: f64 = 5.0;

/// Quality-control scans of one file's raw spectra.
#[derive(Debug, Clone, Default)]
pub struct FileQc {
    /// `None` when the file's MS1 spectra were not read.
    pub polymers: Option<PolymerFileStats>,
    /// `None` when diagnostic ions are off.
    pub diagnostic_ions: Option<DiagnosticFileScan>,
    /// What the file acquired: the identification rate's denominator and the
    /// gradient bins.
    pub census: SpectrumCensus,
    /// Time spent scanning the raw spectra.
    pub scan_ms: u128,
}

/// What one file acquired, counted once when the file is read. The
/// identification rate's denominator and the retention-time bins both come
/// from here, so `qc.tsv` and `qc.gradient.tsv` agree with each other.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct SpectrumCensus {
    pub ms1_spectra: usize,
    /// Spectra of MS level 2 or above.
    pub ms2_spectra: usize,
    /// First and last scan start time read, in minutes.
    pub rt_start: f32,
    pub rt_end: f32,
    /// False when the file's MS1 spectra were not read, which makes every
    /// MS1 count and intensity zero rather than measured.
    pub ms1_read: bool,
    /// [`GRADIENT_BINS`] equal-width bins over `rt_start..=rt_end`; empty
    /// when the file held no spectra.
    pub bins: Vec<GradientBinCensus>,
}

/// One retention-time bin of a file's acquisition.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct GradientBinCensus {
    pub ms2_acquired: usize,
    pub ms1_spectra: usize,
    /// Summed MS1 total ion current of the bin.
    pub ms1_tic: f64,
}

/// Bin of a retention time among [`GRADIENT_BINS`] equal-width bins over
/// `start..=end`, the last bin inclusive of `end`. Everything lands in bin 0
/// when the file spans no time at all.
fn gradient_bin(rt: f32, start: f32, end: f32) -> usize {
    if end <= start {
        return 0;
    }
    let width = (end - start) / GRADIENT_BINS as f32;
    let index = ((rt - start) / width).floor();
    if index < 1.0 {
        // Covers a negative index and a NaN retention time alike.
        return 0;
    }
    (index as usize).min(GRADIENT_BINS - 1)
}

/// A spectrum's total ion current: the value the acquisition states, or the
/// summed peak intensity where it states none (the TDF reader sums the peaks
/// itself, an mzML without the `total ion current` cvParam states nothing).
fn total_ion_current(spectrum: &RawSpectrum) -> f64 {
    if spectrum.total_ion_current > 0.0 {
        return f64::from(spectrum.total_ion_current);
    }
    spectrum
        .intensity
        .iter()
        .map(|&value| f64::from(value))
        .sum()
}

/// Count what a file acquired, in one pass over its raw spectra.
fn census(spectra: &[RawSpectrum]) -> SpectrumCensus {
    let mut census = SpectrumCensus::default();
    if spectra.is_empty() {
        return census;
    }
    let (mut start, mut end) = (f32::INFINITY, f32::NEG_INFINITY);
    for spectrum in spectra {
        if spectrum.scan_start_time.is_finite() {
            start = start.min(spectrum.scan_start_time);
            end = end.max(spectrum.scan_start_time);
        }
    }
    if !start.is_finite() {
        // No file states no retention time, but a synthetic or converted one
        // can; one bin holding everything is then the honest answer.
        start = 0.0;
        end = 0.0;
    }
    census.rt_start = start;
    census.rt_end = end;
    census.bins = vec![GradientBinCensus::default(); GRADIENT_BINS];
    for spectrum in spectra {
        let bin = &mut census.bins[gradient_bin(spectrum.scan_start_time, start, end)];
        match spectrum.ms_level {
            1 => {
                census.ms1_spectra += 1;
                bin.ms1_spectra += 1;
                bin.ms1_tic += total_ion_current(spectrum);
            }
            level if level >= 2 => {
                census.ms2_spectra += 1;
                bin.ms2_acquired += 1;
            }
            // MS level 0 is not a spectrum either side of this counts.
            _ => {}
        }
    }
    census.ms1_read = census.ms1_spectra > 0;
    census
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
            census: census(spectra),
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
            let Some(qc) = file_qc.get(&file_id) else {
                // Every file read records a QC entry; a gap means a read path
                // skipped the scan, so polymer and diagnostic-ion QC miss it.
                warn!("{file}: no quality-control scan was recorded; polymer and diagnostic-ion QC omit this file");
                continue;
            };
            let Some(stats) = qc.polymers.clone() else {
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
                summary: sage_core::digestion::summarize_with_custom_cleavages(
                    enzyme.as_ref(),
                    clip,
                    self.custom_cleavages.as_ref(),
                    passing
                        .iter()
                        .filter(|feature| feature.file_id == file_id)
                        .map(|feature| &self.database[feature.peptide_idx]),
                ),
            })
            .collect();
        let total = sage_core::digestion::summarize_with_custom_cleavages(
            enzyme.as_ref(),
            clip,
            self.custom_cleavages.as_ref(),
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

/// Identification rate of one file and of each of its gradient bins. Counts
/// are of rank-1 target PSMs at [`QC_Q_VALUE`]; a file is identified where
/// its MS2 spectra carry such a PSM, so the share is of spectra, not of PSMs.
#[derive(Debug, Clone, PartialEq)]
pub struct IdentificationFileStats {
    pub file: String,
    /// What the file acquired.
    pub census: SpectrumCensus,
    /// Rank-1 target PSMs at [`QC_Q_VALUE`] spectrum q-value.
    pub psms: usize,
    /// MS2 spectra carrying at least one of those PSMs.
    pub ms2_identified: usize,
    /// Distinct target peptide sequences (modifications ignored) among PSMs
    /// that also pass the peptide q-value, counted as `digestion.tsv` counts
    /// them but without the decoy correction.
    pub peptides: usize,
    /// Distinct protein accessions of PSMs that also pass the protein
    /// q-value, which only proteotypic peptides can.
    pub proteins: usize,
    /// MS1 total ion current summed over the file; zero when the file's MS1
    /// spectra were not read.
    pub ms1_tic: f64,
    pub bins: Vec<GradientBinStats>,
}

/// One gradient bin's acquisition and identification.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct GradientBinStats {
    /// Zero-based bin index, in retention-time order.
    pub bin: usize,
    pub rt_start: f32,
    pub rt_end: f32,
    pub ms2_acquired: usize,
    pub ms2_identified: usize,
    pub ms1_spectra: usize,
    pub ms1_tic: f64,
}

/// Percent of `acquired` spectra that were identified, zero when the
/// denominator is zero.
fn identified_pct(identified: usize, acquired: usize) -> f64 {
    if acquired == 0 {
        return 0.0;
    }
    100.0 * identified as f64 / acquired as f64
}

impl GradientBinStats {
    pub fn identified_pct(&self) -> f64 {
        identified_pct(self.ms2_identified, self.ms2_acquired)
    }
}

impl IdentificationFileStats {
    /// Percent of acquired MS2 spectra that carry a confident PSM. In DIA
    /// pseudo mode the denominator is the acquired isolation windows, while
    /// the PSMs come from pseudo-spectra built from them, so the share can
    /// exceed 100%.
    pub fn identified_pct(&self) -> f64 {
        identified_pct(self.ms2_identified, self.census.ms2_spectra)
    }

    /// Bins that acquired MS2 spectra and identified none: a dead region of
    /// the gradient, or a spray dropout.
    pub fn empty_bins(&self) -> usize {
        self.bins
            .iter()
            .filter(|bin| bin.ms2_acquired > 0 && bin.ms2_identified == 0)
            .count()
    }

    /// Lowest and highest identification rate among the bins that acquired
    /// MS2 spectra; `None` when no bin did.
    pub fn bin_pct_range(&self) -> Option<(f64, f64)> {
        let mut shares = self
            .bins
            .iter()
            .filter(|bin| bin.ms2_acquired > 0)
            .map(|bin| bin.identified_pct());
        let first = shares.next()?;
        Some(shares.fold((first, first), |(low, high), share| {
            (low.min(share), high.max(share))
        }))
    }
}

impl Runner {
    /// Identification rate per file and per gradient bin, from the PSMs that
    /// pass [`QC_Q_VALUE`] and the spectrum census taken when each file was
    /// read.
    pub(super) fn identification_stats(
        &self,
        features: &[Feature],
        filenames: &[String],
    ) -> Vec<IdentificationFileStats> {
        let file_qc = self.file_qc.lock().expect("file QC lock");
        let mut stats = filenames
            .iter()
            .enumerate()
            .map(|(file_id, file)| {
                let census = file_qc
                    .get(&file_id)
                    .map(|qc| qc.census.clone())
                    .unwrap_or_default();
                let bins = census
                    .bins
                    .iter()
                    .enumerate()
                    .map(|(index, bin)| {
                        let width = (census.rt_end - census.rt_start) / GRADIENT_BINS as f32;
                        GradientBinStats {
                            bin: index,
                            rt_start: census.rt_start + width * index as f32,
                            rt_end: census.rt_start + width * (index + 1) as f32,
                            ms2_acquired: bin.ms2_acquired,
                            ms2_identified: 0,
                            ms1_spectra: bin.ms1_spectra,
                            ms1_tic: bin.ms1_tic,
                        }
                    })
                    .collect::<Vec<_>>();
                IdentificationFileStats {
                    file: file.clone(),
                    ms1_tic: census.bins.iter().map(|bin| bin.ms1_tic).sum(),
                    census,
                    psms: 0,
                    ms2_identified: 0,
                    peptides: 0,
                    proteins: 0,
                    bins,
                }
            })
            .collect::<Vec<_>>();
        drop(file_qc);

        // One pass over the PSMs. Spectra are counted through their IDs
        // because a chimeric search can report several rank-1 PSMs for one
        // spectrum, and the share the design asks for is of MS2 spectra.
        let mut identified_spectra = vec![HashSet::new(); stats.len()];
        let mut bin_spectra = vec![vec![HashSet::new(); GRADIENT_BINS]; stats.len()];
        let mut sequences = vec![HashSet::new(); stats.len()];
        let mut proteins = vec![HashSet::new(); stats.len()];
        for feature in features {
            if feature.rank != 1 || feature.label != 1 || feature.spectrum_q > QC_Q_VALUE {
                continue;
            }
            let Some(file) = stats.get_mut(feature.file_id) else {
                continue;
            };
            file.psms += 1;
            identified_spectra[feature.file_id].insert(feature.spec_id.as_str());
            let bin = gradient_bin(feature.rt, file.census.rt_start, file.census.rt_end);
            if let Some(spectra) = bin_spectra[feature.file_id].get_mut(bin) {
                spectra.insert(feature.spec_id.as_str());
            }
            let peptide = &self.database[feature.peptide_idx];
            if feature.peptide_q <= QC_Q_VALUE {
                sequences[feature.file_id].insert(peptide.sequence.as_bytes());
            }
            if feature.protein_q <= QC_Q_VALUE {
                for protein in &peptide.proteins {
                    proteins[feature.file_id].insert(protein.as_ref());
                }
            }
        }
        for (file_id, file) in stats.iter_mut().enumerate() {
            file.ms2_identified = identified_spectra[file_id].len();
            file.peptides = sequences[file_id].len();
            file.proteins = proteins[file_id].len();
            for (bin, spectra) in file.bins.iter_mut().zip(&bin_spectra[file_id]) {
                bin.ms2_identified = spectra.len();
            }
            info!(
                "identification: {}: {} of {} MS2 spectra identified ({:.1}%), {} PSMs, {} peptides, {} proteins at {}% FDR",
                file.file,
                file.ms2_identified,
                file.census.ms2_spectra,
                file.identified_pct(),
                file.psms,
                file.peptides,
                file.proteins,
                QC_Q_VALUE * 100.0
            );
            if file.census.ms2_spectra > 0 && file.ms2_identified == 0 {
                warn!(
                    "{}: none of its {} MS2 spectra were identified; the file may have failed",
                    file.file, file.census.ms2_spectra
                );
            }
        }
        stats
    }

    /// Write `qc.tsv` and `qc.gradient.tsv`.
    pub(super) fn write_qc_tables(
        &self,
        identification: &[IdentificationFileStats],
        digestion: &DigestionRunStats,
        polymers: &[PolymerRunFileStats],
    ) -> anyhow::Result<Vec<Url>> {
        let table = self.make_path("qc.tsv");
        sage_cloudpath::write_bytes_sync(
            &table,
            serialize_qc(identification, digestion, polymers)?,
        )?;
        let gradient = self.make_path("qc.gradient.tsv");
        sage_cloudpath::write_bytes_sync(&gradient, serialize_qc_gradient(identification)?)?;
        Ok(vec![table, gradient])
    }
}

/// `qc.tsv`: one row per input file, repeating the headline number of every
/// QC output so one table shows a whole run. Schema:
/// `schemas/qc.v1.tsv.schema.json`.
pub(super) fn serialize_qc(
    identification: &[IdentificationFileStats],
    digestion: &DigestionRunStats,
    polymers: &[PolymerRunFileStats],
) -> anyhow::Result<Vec<u8>> {
    let mut writer = csv::WriterBuilder::new()
        .delimiter(b'\t')
        .from_writer(Vec::new());
    writer.write_record([
        "file",
        "ms1_spectra",
        "ms2_spectra",
        "rt_start_min",
        "rt_end_min",
        "ms1_read",
        "ms1_tic",
        "psms",
        "ms2_identified",
        "ms2_identified_pct",
        "peptides",
        "proteins",
        "q_value",
        "digestion_peptides",
        "missed_cleavage_pct",
        "semi_n_pct",
        "semi_c_pct",
        "non_enzymatic_pct",
        "top_polymer",
        "top_polymer_tic_pct",
        "gradient_bins",
        "empty_gradient_bins",
        "min_bin_identified_pct",
        "max_bin_identified_pct",
    ])?;
    for file in identification {
        let digested = digestion
            .files
            .iter()
            .find(|row| row.file == file.file)
            .map(|row| &row.summary);
        let polymer = polymers
            .iter()
            .find(|row| row.file == file.file)
            .and_then(|row| {
                row.stats
                    .polymers
                    .iter()
                    .max_by(|left, right| left.tic_pct.total_cmp(&right.tic_pct))
            });
        let (min_pct, max_pct) = file.bin_pct_range().unwrap_or((0.0, 0.0));
        writer.write_record([
            file.file.clone(),
            file.census.ms1_spectra.to_string(),
            file.census.ms2_spectra.to_string(),
            format!("{:.4}", file.census.rt_start),
            format!("{:.4}", file.census.rt_end),
            u8::from(file.census.ms1_read).to_string(),
            format!("{:.1}", file.ms1_tic),
            file.psms.to_string(),
            file.ms2_identified.to_string(),
            format!("{:.3}", file.identified_pct()),
            file.peptides.to_string(),
            file.proteins.to_string(),
            format!("{QC_Q_VALUE}"),
            digested.map_or(0, |summary| summary.peptides).to_string(),
            format!("{:.3}", digested.map_or(0.0, |s| s.missed_cleavage_pct)),
            format!("{:.3}", digested.map_or(0.0, |s| s.semi_n_pct)),
            format!("{:.3}", digested.map_or(0.0, |s| s.semi_c_pct)),
            format!("{:.3}", digested.map_or(0.0, |s| s.non_enzymatic_pct)),
            polymer.map_or("-", |share| share.name.as_str()).to_string(),
            format!("{:.3}", polymer.map_or(0.0, |share| share.tic_pct)),
            file.bins.len().to_string(),
            file.empty_bins().to_string(),
            format!("{min_pct:.3}"),
            format!("{max_pct:.3}"),
        ])?;
    }
    writer.flush()?;
    Ok(writer.into_inner()?)
}

/// `qc.gradient.tsv`: one row per file and retention-time bin. Schema:
/// `schemas/qc_gradient.v1.tsv.schema.json`.
pub(super) fn serialize_qc_gradient(
    identification: &[IdentificationFileStats],
) -> anyhow::Result<Vec<u8>> {
    let mut writer = csv::WriterBuilder::new()
        .delimiter(b'\t')
        .from_writer(Vec::new());
    writer.write_record([
        "file",
        "bin",
        "rt_start_min",
        "rt_end_min",
        "ms2_acquired",
        "ms2_identified",
        "ms2_identified_pct",
        "ms1_spectra",
        "ms1_tic",
    ])?;
    for file in identification {
        for bin in &file.bins {
            writer.write_record([
                file.file.clone(),
                bin.bin.to_string(),
                format!("{:.4}", bin.rt_start),
                format!("{:.4}", bin.rt_end),
                bin.ms2_acquired.to_string(),
                bin.ms2_identified.to_string(),
                format!("{:.3}", bin.identified_pct()),
                bin.ms1_spectra.to_string(),
                format!("{:.1}", bin.ms1_tic),
            ])?;
        }
    }
    writer.flush()?;
    Ok(writer.into_inner()?)
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
#[path = "../../tests/unit/runner/qc.rs"]
mod test;
