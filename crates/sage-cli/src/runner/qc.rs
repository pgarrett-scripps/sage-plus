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
use sage_core::occupancy::{OccupancyCounts, OccupancyRule};
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

/// Share of a file's MS2 spectra a gradient bin must have acquired before it
/// counts towards the headline bin columns of `qc.tsv`
/// (`empty_gradient_bins`, `min_bin_identified_pct`,
/// `max_bin_identified_pct`). The edge bins of a gradient can hold a handful
/// of scans, and an identification rate over two spectra is noise that would
/// otherwise set a whole file's reported range. `qc.gradient.tsv` keeps every
/// bin with its own counts, so nothing is hidden.
const HEADLINE_BIN_MS2_SHARE: f64 = 0.01;

/// Spectrum- and peptide-level q-value a PSM must pass to count towards the
/// modification occupancy table. The same threshold as
/// [`DIGESTION_Q_VALUE`], so both tables describe the same identifications.
const OCCUPANCY_Q_VALUE: f32 = 0.01;

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
            // An MS1 scan a reader kept without decoding its peaks, which
            // is what a Thermo `.raw` searched without quantification gives,
            // still dates the file above, and so sets the retention-time span
            // and the bin edges. Nothing was measured in it, so it counts as
            // no MS1 spectrum and no TIC, and `ms1_read` stays false.
            1 if !spectrum.metadata_only => {
                census.ms1_spectra += 1;
                bin.ms1_spectra += 1;
                bin.ms1_tic += total_ion_current(spectrum);
            }
            1 => {}
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
    /// Present when the search has at least one variable modification.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub modifications: Option<ModificationOccupancyRunStats>,
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
/// has no MS1 spectra whose peaks were read: a polymer share is a share of
/// MS1 intensity, so a run that never decoded its MS1 peaks has nothing to
/// scan.
fn scan_polymers(spectra: &[RawSpectrum]) -> Option<PolymerFileStats> {
    let scanner = PolymerScanner::default();
    let empty = PolymerFileStats::new(&scanner.polymers);
    let mut stats = spectra
        .par_iter()
        .filter(|spectrum| spectrum.ms_level == 1 && !spectrum.metadata_only)
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

/// Occupancy of every searched variable modification, counted over distinct
/// confident target peptidoforms. Same rows as `qc.modifications.tsv`.
#[derive(Debug, Clone, Default, serde::Serialize, serde::Deserialize)]
pub struct ModificationOccupancyRunStats {
    pub q_value: f32,
    /// Every file's peptidoforms together, each counted once.
    pub modifications: Vec<ModificationOccupancy>,
    pub files: Vec<ModificationOccupancyFileStats>,
}

#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct ModificationOccupancyFileStats {
    pub file: String,
    pub modifications: Vec<ModificationOccupancy>,
}

/// One searched definition and site rule in one file, or in the whole run.
#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct ModificationOccupancy {
    /// The site rule as a configuration key writes it, e.g. `M`,
    /// `peptide_n_term` or `internal_residue:K`.
    pub site_rule: String,
    pub mass: f32,
    /// The configured name; absent for a definition given as a bare mass.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub name: Option<String>,
    /// `database` or `mass_offset`.
    pub search_mode: String,
    /// Peptidoforms with at least one eligible site.
    pub peptidoforms: usize,
    pub modified_peptidoforms: usize,
    pub eligible_sites: usize,
    pub modified_sites: usize,
    /// `modified_sites` as a percent of `eligible_sites`.
    pub occupancy_pct: f64,
}

fn occupancy_rows(
    rules: &[OccupancyRule],
    counts: &[OccupancyCounts],
) -> Vec<ModificationOccupancy> {
    rules
        .iter()
        .zip(counts)
        .map(|(rule, counts)| ModificationOccupancy {
            site_rule: rule.site_rule_name(),
            mass: rule.mass(),
            name: rule.name().map(str::to_string),
            search_mode: rule.search_mode_name().to_string(),
            peptidoforms: counts.peptidoforms,
            modified_peptidoforms: counts.modified_peptidoforms,
            eligible_sites: counts.eligible_sites,
            modified_sites: counts.modified_sites,
            occupancy_pct: counts.occupancy_pct(),
        })
        .collect()
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
        // `polymer_check` off means the file's MS1 peaks may never have been
        // decoded, and means the user asked for no polymer scan either way.
        let polymers = self
            .parameters
            .polymer_check
            .then(|| scan_polymers(spectra))
            .flatten();
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
    /// Occupancy of every searched variable modification per file; `None`
    /// when the search has no variable modification.
    ///
    /// Each distinct peptidoform of a file counts once, however many PSMs
    /// support it, and the run total counts each peptidoform once across
    /// files, so it is not the sum of the file rows. Decoys are excluded:
    /// a decoy peptidoform's sites say nothing about the sample.
    pub(super) fn modification_occupancy_stats(
        &self,
        features: &[Feature],
        filenames: &[String],
    ) -> Option<ModificationOccupancyRunStats> {
        let rules = sage_core::occupancy::rules(&self.database_parameters);
        if rules.is_empty() {
            return None;
        }
        let mut per_file = vec![Vec::new(); filenames.len()];
        let mut run = Vec::new();
        let mut seen = HashSet::new();
        let mut seen_in_run = HashSet::new();
        for feature in features {
            if feature.rank != 1
                || feature.spectrum_q > OCCUPANCY_Q_VALUE
                || feature.peptide_q > OCCUPANCY_Q_VALUE
            {
                continue;
            }
            let peptide = &self.database[feature.peptide_idx];
            if peptide.decoy {
                continue;
            }
            if let Some(peptidoforms) = per_file.get_mut(feature.file_id) {
                if seen.insert((feature.file_id, feature.peptide_idx)) {
                    peptidoforms.push(peptide);
                }
            }
            if seen_in_run.insert(feature.peptide_idx) {
                run.push(peptide);
            }
        }
        let files = filenames
            .iter()
            .zip(&per_file)
            .map(|(file, peptidoforms)| ModificationOccupancyFileStats {
                file: file.clone(),
                modifications: occupancy_rows(
                    &rules,
                    &sage_core::occupancy::count(
                        &self.database_parameters,
                        &rules,
                        peptidoforms.iter().copied(),
                    ),
                ),
            })
            .collect();
        let modifications = occupancy_rows(
            &rules,
            &sage_core::occupancy::count(&self.database_parameters, &rules, run),
        );
        info!(
            "modification occupancy: {}",
            modifications
                .iter()
                .map(|row| format!(
                    "{} {:+.4} {:.2}% of {} sites",
                    row.site_rule, row.mass, row.occupancy_pct, row.eligible_sites
                ))
                .collect::<Vec<_>>()
                .join(", ")
        );
        Some(ModificationOccupancyRunStats {
            q_value: OCCUPANCY_Q_VALUE,
            modifications,
            files,
        })
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

    /// Write `qc.modifications.tsv`, or nothing when the search has no
    /// variable modification.
    pub(super) fn write_modification_occupancy(
        &self,
        stats: Option<&ModificationOccupancyRunStats>,
    ) -> anyhow::Result<Option<Url>> {
        let Some(stats) = stats else {
            return Ok(None);
        };
        let path = self.make_path("qc.modifications.tsv");
        sage_cloudpath::write_bytes_sync(&path, serialize_modification_occupancy(stats)?)?;
        Ok(Some(path))
    }

    pub(super) fn write_digestion(&self, stats: &DigestionRunStats) -> anyhow::Result<Url> {
        let path = self.make_path("digestion.tsv");
        sage_cloudpath::write_bytes_sync(&path, serialize_digestion(stats)?)?;
        Ok(path)
    }
}

pub(super) fn serialize_modification_occupancy(
    stats: &ModificationOccupancyRunStats,
) -> anyhow::Result<Vec<u8>> {
    let mut writer = csv::WriterBuilder::new()
        .delimiter(b'\t')
        .from_writer(Vec::new());
    writer.write_record([
        "file",
        "site_rule",
        "mass",
        "name",
        "search_mode",
        "peptidoforms",
        "modified_peptidoforms",
        "eligible_sites",
        "modified_sites",
        "occupancy_pct",
    ])?;
    let rows = stats
        .files
        .iter()
        .map(|file| (file.file.as_str(), &file.modifications))
        .chain(std::iter::once(("total", &stats.modifications)));
    for (file, modifications) in rows {
        for row in modifications {
            writer.write_record([
                file.to_string(),
                row.site_rule.clone(),
                format!("{:.6}", row.mass),
                row.name.clone().unwrap_or_default(),
                row.search_mode.clone(),
                row.peptidoforms.to_string(),
                row.modified_peptidoforms.to_string(),
                row.eligible_sites.to_string(),
                row.modified_sites.to_string(),
                format!("{:.3}", row.occupancy_pct),
            ])?;
        }
    }
    writer.flush()?;
    Ok(writer.into_inner()?)
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

    /// Bins carrying enough of the file's MS2 spectra for their
    /// identification rate to mean something: at least
    /// [`HEADLINE_BIN_MS2_SHARE`] of them, and at least one. A file that
    /// acquired no MS2 spectrum has none.
    fn headline_bins(&self) -> impl Iterator<Item = &GradientBinStats> {
        let floor = HEADLINE_BIN_MS2_SHARE * self.census.ms2_spectra as f64;
        self.bins
            .iter()
            .filter(move |bin| bin.ms2_acquired > 0 && bin.ms2_acquired as f64 >= floor)
    }

    /// Bins that acquired a meaningful share of the file's MS2 spectra and
    /// identified none: a dead region of the gradient, or a spray dropout. A
    /// bin that acquired no MS2 spectrum is not one of these, because nothing
    /// was tried there; those are counted by [`Self::no_ms2_bins`]. A bin
    /// holding a handful of scans is not one either, because a rate over a
    /// handful of spectra is noise; see [`HEADLINE_BIN_MS2_SHARE`].
    pub fn empty_bins(&self) -> usize {
        self.headline_bins()
            .filter(|bin| bin.ms2_identified == 0)
            .count()
    }

    /// Bins that acquired no MS2 spectrum at all, which is ordinary at the
    /// start and end of a gradient and says nothing about identification.
    /// Literal: a bin with one MS2 spectrum is not counted here, however
    /// little that is.
    pub fn no_ms2_bins(&self) -> usize {
        self.bins.iter().filter(|bin| bin.ms2_acquired == 0).count()
    }

    /// Lowest and highest identification rate among the bins that carry the
    /// headline columns ([`Self::headline_bins`]); `None` when no bin does.
    pub fn bin_pct_range(&self) -> Option<(f64, f64)> {
        let mut shares = self.headline_bins().map(|bin| bin.identified_pct());
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
    pub(super) fn write_qc_tables(&self, sources: &QcSources<'_>) -> anyhow::Result<Vec<Url>> {
        let table = self.make_path("qc.tsv");
        sage_cloudpath::write_bytes_sync(&table, serialize_qc(sources)?)?;
        let gradient = self.make_path("qc.gradient.tsv");
        sage_cloudpath::write_bytes_sync(
            &gradient,
            serialize_qc_gradient(sources.identification)?,
        )?;
        Ok(vec![table, gradient])
    }
}

/// `qc.tsv`: one row per input file, repeating the headline number of every
/// QC output so one table shows a whole run. Schema:
/// `schemas/qc.v1.tsv.schema.json`.
/// Every QC pass `qc.tsv` repeats a headline number from, borrowed. Nothing
/// here is recomputed: each field is the stats object that pass already built,
/// so a number in `qc.tsv` is the same number as in its own output.
pub(super) struct QcSources<'a> {
    pub identification: &'a [IdentificationFileStats],
    pub digestion: &'a DigestionRunStats,
    pub polymers: &'a [PolymerRunFileStats],
    pub composition: &'a CompositionProfiles,
    pub chromatography: &'a ChromatographyRunStats,
    /// `None` when no modification was searched.
    pub modifications: Option<&'a ModificationOccupancyRunStats>,
    /// `None` without a `contaminant_fasta`.
    pub contaminants: Option<&'a ContaminantRunStats>,
    /// `None` when the in-source annotation is off.
    pub in_source: Option<&'a sage_core::in_source::InSourceAnnotations>,
}

pub(super) fn serialize_qc(sources: &QcSources<'_>) -> anyhow::Result<Vec<u8>> {
    let QcSources {
        identification,
        digestion,
        polymers,
        composition,
        chromatography,
        modifications,
        contaminants,
        in_source,
    } = *sources;
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
        "largest_shift_residue",
        "largest_shift_log2_ratio",
        "jensen_shannon_distance",
        "charge_1_pct",
        "charge_2_pct",
        "charge_3_pct",
        "charge_4_plus_pct",
        "charge_unknown_pct",
        "peak_width",
        "peak_width_kind",
        "quantified_precursors",
        "transferred_pct",
        "modifications_searched",
        "top_modification",
        "top_modification_occupancy_pct",
        "contaminant_psms",
        "contaminant_psm_pct",
        "in_source_flagged_psms",
        "gradient_bins",
        "empty_gradient_bins",
        "no_ms2_gradient_bins",
        "min_bin_identified_pct",
        "max_bin_identified_pct",
    ])?;
    for (file_id, file) in identification.iter().enumerate() {
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
        let shift = composition
            .files
            .iter()
            .find(|profile| profile.file == file.file)
            .map(|profile| &profile.summary);
        let largest_shift = shift.and_then(|summary| summary.largest_shift);
        let chromatography = chromatography
            .files
            .iter()
            .find(|row| row.file == file.file)
            .map(|row| &row.stats);
        // The definition with the highest occupancy in this file, named as
        // qc.modifications.tsv names it: the configured name, or the site rule
        // and mass for a definition given as a bare mass.
        let top_modification = modifications
            .and_then(|stats| stats.files.iter().find(|row| row.file == file.file))
            .map(|row| row.modifications.as_slice());
        let top = top_modification.and_then(|rows| {
            rows.iter()
                .max_by(|left, right| left.occupancy_pct.total_cmp(&right.occupancy_pct))
        });
        let contaminant =
            contaminants.and_then(|stats| stats.files.iter().find(|row| row.filename == file.file));
        let flagged = in_source.map(|annotations| {
            annotations
                .files()
                .iter()
                .filter(|row| row.file_id == file_id)
                .map(|row| row.flagged_psms)
                .sum::<usize>()
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
            // Empty, as in qc.composition.tsv, where no residue has a finite
            // ratio or neither profile counted residues: an absent shift is
            // not a shift of zero.
            largest_shift.map_or(String::new(), |(residue, _)| residue.to_string()),
            largest_shift.map_or(String::new(), |(_, ratio)| format!("{ratio:.4}")),
            shift
                .and_then(|summary| summary.jensen_shannon_distance)
                .map_or(String::new(), |distance| format!("{distance:.4}")),
            // Chromatography: the charge shares of this file's confident
            // precursors, its peak width with the quantity that width is, and
            // what quantification found. Empty where that pass measured
            // nothing, which is not a zero.
            chromatography.map_or(String::new(), |stats| format!("{:.3}", stats.charge_pct(1))),
            chromatography.map_or(String::new(), |stats| format!("{:.3}", stats.charge_pct(2))),
            chromatography.map_or(String::new(), |stats| format!("{:.3}", stats.charge_pct(3))),
            chromatography.map_or(String::new(), |stats| format!("{:.3}", stats.charge_pct(4))),
            chromatography.map_or(String::new(), |stats| {
                format!("{:.3}", stats.charge_unknown_pct())
            }),
            chromatography
                .and_then(|stats| stats.peak_width)
                .map_or(String::new(), |width| format!("{width:.4}")),
            chromatography
                .and_then(|stats| stats.peak_width_kind)
                .map_or(String::new(), |kind| kind.label().to_string()),
            chromatography
                .and_then(|stats| stats.quantified_precursors)
                .map_or(String::new(), |count| count.to_string()),
            chromatography
                .and_then(|stats| stats.transferred_pct())
                .map_or(String::new(), |pct| format!("{pct:.3}")),
            // Modification occupancy: how many definitions were searched and
            // which one is most occupied in this file.
            top_modification.map_or(String::new(), |rows| rows.len().to_string()),
            top.map_or(String::new(), modification_label),
            top.map_or(String::new(), |row| format!("{:.3}", row.occupancy_pct)),
            // Contaminants and in-source products, empty when those are off.
            contaminant.map_or(String::new(), |row| row.contaminant_psms.to_string()),
            contaminant.map_or(String::new(), |row| {
                format!("{:.3}", row.contaminant_psm_pct)
            }),
            flagged.map_or(String::new(), |count| count.to_string()),
            file.bins.len().to_string(),
            file.empty_bins().to_string(),
            file.no_ms2_bins().to_string(),
            format!("{min_pct:.3}"),
            format!("{max_pct:.3}"),
        ])?;
    }
    writer.flush()?;
    Ok(writer.into_inner()?)
}

/// A searched modification named as `qc.modifications.tsv` names it: the
/// configured name, or the site rule and mass for a definition given as a bare
/// mass, so one modification reads the same in both tables.
fn modification_label(occupancy: &ModificationOccupancy) -> String {
    match &occupancy.name {
        Some(name) => name.clone(),
        None => format!("{} {:+.4}", occupancy.site_rule, occupancy.mass),
    }
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
