//! Run-level quality-control summaries written next to the search results.

use super::*;
use fnv::FnvHashMap;
use sage_core::chromatography::{
    quantification_by_file, repeat_psm_spread, ChromatographyStats, PeakWidth, RepeatRts,
};
use sage_core::composition::{CompositionProfile, CompositionSummary};
use sage_core::diagnostic::{find_ions, DiagnosticIon};
use sage_core::digestion::DigestionSummary;
use sage_core::enzyme::EnzymeParameters;
use sage_core::lfq::{PrecursorId, QuantifiedPeak};
use sage_core::polymer::{PolymerFileStats, PolymerScanner};
use sage_core::spectrum::{RawSpectrum, Representation};

/// Spectrum- and peptide-level q-value a PSM must pass to count towards the
/// digestion summary, matching the run summary's 1% FDR counts.
const DIGESTION_Q_VALUE: f32 = 0.01;

/// Percent of MS1 TIC on one polymer above which a run warning is raised.
const POLYMER_WARNING_TIC_PCT: f64 = 5.0;

/// Spectrum- and peptide-level q-value a PSM must pass to count towards the
/// amino acid composition profile, matching the digestion summary.
const COMPOSITION_Q_VALUE: f32 = 0.01;

/// The same for the charge and chromatography summary.
const CHROMATOGRAPHY_Q_VALUE: f32 = 0.01;

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
    /// Headline numbers of the amino acid composition profile. The per-residue
    /// table is `qc.composition.tsv`.
    #[serde(default)]
    pub composition: CompositionRunStats,
    /// Charge states and peak width per file; the same rows are
    /// `qc.chromatography.tsv`.
    #[serde(default)]
    pub chromatography: ChromatographyRunStats,
}

/// Charge-state shares, peak width and quantification counts per file.
#[derive(Debug, Clone, Default, serde::Serialize, serde::Deserialize)]
pub struct ChromatographyRunStats {
    pub q_value: f32,
    pub files: Vec<ChromatographyFileStats>,
    /// Every file together. The charge shares and the quantification counts
    /// pool; the peak width does not, because a width in one file's
    /// retention-time unit is not comparable with another's.
    pub total: ChromatographyStats,
}

#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct ChromatographyFileStats {
    pub file: String,
    #[serde(flatten)]
    pub stats: ChromatographyStats,
}

/// Headline numbers of the amino acid composition comparison, per file and for
/// the whole run. The per-residue rows go to `qc.composition.tsv`, which would
/// make this block twenty times longer for no reader's benefit.
#[derive(Debug, Clone, Default, serde::Serialize, serde::Deserialize)]
pub struct CompositionRunStats {
    pub q_value: f32,
    pub files: Vec<CompositionHeadline>,
    pub total: CompositionHeadline,
}

#[derive(Debug, Clone, Default, serde::Serialize, serde::Deserialize)]
pub struct CompositionHeadline {
    pub file: String,
    /// Distinct confident target peptide sequences observed.
    pub observed_peptides: u64,
    pub observed_residues: u64,
    /// Observed residues whose missed-cleavage class the database could match,
    /// which are the ones weighting the expected profile.
    pub matched_observed_residues: u64,
    /// Distinct database peptide sequences of the proteins identified here.
    pub expected_peptides: u64,
    pub expected_residues: u64,
    /// Residue whose observed share differs most from the expected one, and
    /// that difference as a log2 ratio. Absent when no residue has a finite
    /// ratio, for example when nothing was identified.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub largest_shift_residue: Option<char>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub largest_shift_log2_ratio: Option<f64>,
    /// Jensen-Shannon distance over the whole profile; absent when either
    /// profile counted no residues.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub jensen_shannon_distance: Option<f64>,
}

impl CompositionHeadline {
    fn new(file: &str, summary: &CompositionSummary) -> Self {
        Self {
            file: file.to_string(),
            observed_peptides: summary.observed_peptides,
            observed_residues: summary.observed_residues,
            matched_observed_residues: summary.matched_observed_residues,
            expected_peptides: summary.expected_peptides,
            expected_residues: summary.expected_residues,
            largest_shift_residue: summary.largest_shift.map(|(residue, _)| residue),
            largest_shift_log2_ratio: summary.largest_shift.map(|(_, ratio)| ratio),
            jensen_shannon_distance: summary.jensen_shannon_distance,
        }
    }
}

/// The full per-file composition profiles, including the per-residue rows.
#[derive(Debug, Clone, Default)]
pub struct CompositionProfiles {
    pub q_value: f32,
    /// One entry per searched file, in input order.
    pub files: Vec<CompositionFileProfile>,
    /// Every file together: distinct sequences over the run, against the
    /// database peptides of every protein identified anywhere in it.
    pub total: CompositionSummary,
}

#[derive(Debug, Clone)]
pub struct CompositionFileProfile {
    pub file: String,
    pub summary: CompositionSummary,
}

impl CompositionProfiles {
    /// The small block written to `run-summary.json`.
    pub fn run_stats(&self) -> CompositionRunStats {
        CompositionRunStats {
            q_value: self.q_value,
            files: self
                .files
                .iter()
                .map(|file| CompositionHeadline::new(&file.file, &file.summary))
                .collect(),
            total: CompositionHeadline::new("total", &self.total),
        }
    }
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

    /// Interior amino acid composition of each file's confident target
    /// peptides, against the composition of the database peptides of the
    /// proteins identified in that file.
    ///
    /// The expected profile is taken from the database as built, so a
    /// prefiltered or limited database makes a prefiltered expectation.
    pub(super) fn composition_stats(
        &self,
        features: &[Feature],
        filenames: &[String],
    ) -> CompositionProfiles {
        use std::collections::{HashMap, HashSet};

        let start = Instant::now();
        let enzyme = EnzymeParameters::from(self.database_parameters.enzyme.clone()).enzyme;
        let clip = self.database_parameters.clip_n_term_met;
        // One slot per file, and one more for the run as a whole.
        let run = filenames.len();
        let mut observed = vec![CompositionProfile::default(); run + 1];
        let mut observed_seen = vec![HashSet::new(); run + 1];
        let mut proteins = vec![HashSet::new(); run + 1];
        for feature in features {
            if feature.rank != 1
                || feature.spectrum_q > COMPOSITION_Q_VALUE
                || feature.peptide_q > COMPOSITION_Q_VALUE
            {
                continue;
            }
            let peptide = &self.database[feature.peptide_idx];
            if peptide.decoy {
                continue;
            }
            let sequence = peptide.sequence.as_str();
            for slot in [feature.file_id, run] {
                // Distinct sequences, so an abundant peptide counts once.
                if observed_seen[slot].insert(sequence) {
                    observed[slot].add_peptide(enzyme.as_ref(), clip, peptide);
                }
                for protein in &peptide.proteins {
                    proteins[slot].insert(&**protein);
                }
            }
        }

        // A protein identified in several files credits its peptides to each.
        let mut protein_files: HashMap<&str, Vec<usize>> = HashMap::new();
        for (slot, identified) in proteins.iter().enumerate() {
            for protein in identified {
                protein_files.entry(protein).or_default().push(slot);
            }
        }
        let mut expected = vec![CompositionProfile::default(); run + 1];
        let mut expected_seen: HashSet<&str> = HashSet::new();
        let mut credited: Vec<usize> = Vec::new();
        for peptide in &self.database.peptides {
            if peptide.decoy {
                continue;
            }
            credited.clear();
            for protein in &peptide.proteins {
                let slots = protein_files
                    .get(&**protein)
                    .map_or(&[][..], |slots| slots.as_slice());
                for &slot in slots {
                    if !credited.contains(&slot) {
                        credited.push(slot);
                    }
                }
            }
            if credited.is_empty() {
                continue;
            }
            // Every modification variant of a sequence shares its proteins, so
            // the first occurrence settles which files it is counted in.
            if !expected_seen.insert(peptide.sequence.as_str()) {
                continue;
            }
            let mut counts = CompositionProfile::default();
            counts.add_peptide(enzyme.as_ref(), clip, peptide);
            for &slot in &credited {
                expected[slot].merge(&counts);
            }
        }

        let files = filenames
            .iter()
            .enumerate()
            .map(|(slot, file)| CompositionFileProfile {
                file: file.clone(),
                summary: CompositionSummary::new(&observed[slot], &expected[slot]),
            })
            .collect();
        let total = CompositionSummary::new(&observed[run], &expected[run]);
        info!(
            "composition: {} peptides, {} interior residues ({} matched on missed cleavages) against {} database peptides of {} proteins; largest shift {}, Jensen-Shannon distance {} ({} ms)",
            total.observed_peptides,
            total.observed_residues,
            total.matched_observed_residues,
            total.expected_peptides,
            proteins[run].len(),
            match total.largest_shift {
                Some((residue, ratio)) => format!("{residue} {ratio:+.2} log2"),
                None => "none".to_string(),
            },
            match total.jensen_shannon_distance {
                Some(distance) => format!("{distance:.4}"),
                None => "none".to_string(),
            },
            start.elapsed().as_millis(),
        );
        CompositionProfiles {
            q_value: COMPOSITION_Q_VALUE,
            files,
            total,
        }
    }

    /// Charge states of each file's confident precursors, a labelled peak
    /// width, and what quantification did in that file.
    ///
    /// `areas` is the LFQ feature map when quantification ran. With it the peak
    /// width is the median full width at half maximum of the file's integrated
    /// peaks; without it, or for a file whose peaks never fell to half height,
    /// the width is the spread of that file's repeat PSMs instead, and the kind
    /// column says which of the two a row carries.
    pub(super) fn chromatography_stats(
        &self,
        features: &[Feature],
        filenames: &[String],
        areas: Option<&FnvHashMap<(PrecursorId, bool), QuantifiedPeak>>,
    ) -> ChromatographyRunStats {
        let start = Instant::now();
        let mut files = vec![ChromatographyStats::default(); filenames.len()];
        let mut total = ChromatographyStats::default();
        let mut rts = vec![RepeatRts::default(); filenames.len()];
        for feature in features {
            if feature.rank != 1
                || feature.label != 1
                || feature.spectrum_q > CHROMATOGRAPHY_Q_VALUE
                || feature.peptide_q > CHROMATOGRAPHY_Q_VALUE
            {
                continue;
            }
            let Some(file) = files.get_mut(feature.file_id) else {
                continue;
            };
            file.add_charge(feature.charge);
            total.add_charge(feature.charge);
            rts[feature.file_id]
                .entry(feature.peptide_idx)
                .or_default()
                .push(feature.rt);
        }

        let quantification = areas.map(|areas| quantification_by_file(areas, filenames.len()));
        for (file_id, stats) in files.iter_mut().enumerate() {
            let fwhm = quantification.as_ref().and_then(|by_file| {
                let (quantified, transferred, fwhm) = by_file[file_id];
                stats.quantified_precursors = Some(quantified);
                stats.transferred_precursors = Some(transferred);
                fwhm
            });
            match fwhm {
                Some(fwhm) => {
                    stats.peak_width = Some(fwhm);
                    stats.peak_width_kind = Some(PeakWidth::LfqMedianFwhm);
                }
                // No peaks to measure: the repeat-PSM spread stands in, and the
                // kind column keeps the two apart.
                None => {
                    if let Some((spread, kind)) = repeat_psm_spread(&rts[file_id]) {
                        stats.peak_width = Some(spread);
                        stats.peak_width_kind = Some(kind);
                    }
                }
            }
        }
        if let Some(by_file) = quantification.as_ref() {
            total.quantified_precursors =
                Some(by_file.iter().map(|&(quantified, _, _)| quantified).sum());
            total.transferred_precursors =
                Some(by_file.iter().map(|&(_, transferred, _)| transferred).sum());
        }
        info!(
            "chromatography: charges 1+ {:.1}%, 2+ {:.1}%, 3+ {:.1}%, 4+ and above {:.1}% of {} confident precursors{} ({} ms)",
            total.charge_pct(1),
            total.charge_pct(2),
            total.charge_pct(3),
            total.charge_pct(4),
            total.confident_precursors,
            match (total.quantified_precursors, total.transferred_pct()) {
                (Some(quantified), Some(transferred)) => format!(
                    "; {quantified} quantified precursor/file pairs, {transferred:.1}% without direct MS2 evidence in their file"
                ),
                (Some(quantified), None) => format!("; {quantified} quantified precursor/file pairs"),
                _ => String::new(),
            },
            start.elapsed().as_millis(),
        );
        for (file, stats) in filenames.iter().zip(&files) {
            if let (Some(width), Some(kind)) = (stats.peak_width, stats.peak_width_kind) {
                info!("- {file}: {} {width:.4}", kind.label());
            }
        }
        ChromatographyRunStats {
            q_value: CHROMATOGRAPHY_Q_VALUE,
            files: filenames
                .iter()
                .cloned()
                .zip(files)
                .map(|(file, stats)| ChromatographyFileStats { file, stats })
                .collect(),
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

    pub(super) fn write_chromatography(
        &self,
        stats: &ChromatographyRunStats,
    ) -> anyhow::Result<Url> {
        let path = self.make_path("qc.chromatography.tsv");
        sage_cloudpath::write_bytes_sync(&path, serialize_chromatography(stats)?)?;
        Ok(path)
    }

    pub(super) fn write_composition(&self, profiles: &CompositionProfiles) -> anyhow::Result<Url> {
        let path = self.make_path("qc.composition.tsv");
        sage_cloudpath::write_bytes_sync(&path, serialize_composition(profiles)?)?;
        Ok(path)
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

pub(super) fn serialize_composition(profiles: &CompositionProfiles) -> anyhow::Result<Vec<u8>> {
    let mut writer = csv::WriterBuilder::new()
        .delimiter(b'\t')
        .from_writer(Vec::new());
    writer.write_record([
        "file",
        "residue",
        "observed",
        "observed_fraction",
        "expected",
        "expected_fraction",
        "log2_ratio",
    ])?;
    let rows = profiles
        .files
        .iter()
        .map(|file| (file.file.as_str(), &file.summary))
        .chain(std::iter::once(("total", &profiles.total)));
    for (file, summary) in rows {
        for residue in &summary.residues {
            writer.write_record([
                file.to_string(),
                residue.residue.to_string(),
                residue.observed.to_string(),
                format!("{:.6}", residue.observed_fraction),
                residue.expected.to_string(),
                format!("{:.6}", residue.expected_fraction),
                // Empty when one side has none of the residue, where the
                // ratio has no finite value.
                residue
                    .log2_ratio
                    .map(|ratio| format!("{ratio:.4}"))
                    .unwrap_or_default(),
            ])?;
        }
    }
    writer.flush()?;
    Ok(writer.into_inner()?)
}

pub(super) fn serialize_chromatography(stats: &ChromatographyRunStats) -> anyhow::Result<Vec<u8>> {
    let mut writer = csv::WriterBuilder::new()
        .delimiter(b'\t')
        .from_writer(Vec::new());
    writer.write_record([
        "file",
        "confident_precursors",
        "charge_1_pct",
        "charge_2_pct",
        "charge_3_pct",
        "charge_4_plus_pct",
        "charge_unknown_pct",
        "peak_width",
        "peak_width_kind",
        "quantified_precursors",
        "transferred_precursors",
        "transferred_pct",
    ])?;
    let rows = stats
        .files
        .iter()
        .map(|file| (file.file.as_str(), &file.stats))
        .chain(std::iter::once(("total", &stats.total)));
    let optional = |value: Option<u64>| value.map(|v| v.to_string()).unwrap_or_default();
    for (file, stats) in rows {
        writer.write_record([
            file.to_string(),
            stats.confident_precursors.to_string(),
            format!("{:.3}", stats.charge_pct(1)),
            format!("{:.3}", stats.charge_pct(2)),
            format!("{:.3}", stats.charge_pct(3)),
            format!("{:.3}", stats.charge_pct(4)),
            format!("{:.3}", stats.charge_unknown_pct()),
            // Empty when the file had neither peaks nor repeat PSMs to measure.
            stats
                .peak_width
                .map(|width| format!("{width:.4}"))
                .unwrap_or_default(),
            stats
                .peak_width_kind
                .map(|kind| kind.label().to_string())
                .unwrap_or_default(),
            optional(stats.quantified_precursors),
            optional(stats.transferred_precursors),
            stats
                .transferred_pct()
                .map(|pct| format!("{pct:.3}"))
                .unwrap_or_default(),
        ])?;
    }
    writer.flush()?;
    Ok(writer.into_inner()?)
}

#[cfg(test)]
#[path = "../../tests/unit/runner/qc.rs"]
mod test;
