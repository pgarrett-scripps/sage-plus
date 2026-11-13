//! In-source fragmentation and in-source neutral loss annotation.
//!
//! A peptide can reach the mass analyser as a fragment of a longer peptide, or
//! having lost a small neutral molecule, because it fragmented in the ion
//! source rather than in the collision cell. Such a product is a real spectrum
//! of a real ion, so the search identifies it, but it is not an experimental
//! peptide: it is a product of one that was in the sample.
//!
//! This module annotates those PSMs. It changes no score, removes no row and
//! reweights nothing: every annotated PSM keeps the identification it had. Two
//! things are reported per flagged PSM, the *parent* peptide the product came
//! from and what kind of product it is, so a downstream analysis can filter or
//! group on them.
//!
//! A confident target PSM is flagged when both hold in the same file:
//!
//! 1. its peptide occupies a span of a longer confident peptide of the same
//!    protein, or the same span with a lower mass; the modifications on the
//!    shared residues are equal except at one site; and its mass is the mass a
//!    product of that span would have (see [`Placement`]); and
//! 2. it co-elutes with the parent in that file. Where quantification gave
//!    both a peak there and the PSM fired inside the product's own peak, the
//!    two peaks' apexes decide, within a fraction of the peak width;
//!    otherwise the PSM's retention time is compared with the parent's median
//!    retention time, within a tolerance derived from the file itself.
//!
//! The flag names the parent: the longest co-eluting peptide that contains the
//! product. A product of a product therefore points at the intact peptide,
//! with no chain to follow.
//!
//! The decision is per PSM, not per peptide. The same peptide can be an
//! experimental peptide at one retention time and an in-source product at
//! another in the same file, which is what happens when a labile modification
//! is lost in the source.
//!
//! Containment is decided in protein coordinates, from
//! [`Peptide::protein_sites`], not by searching one peptide's sequence inside
//! another: a product and its parent then come from the same place in the same
//! protein, which is what in-source fragmentation means.

use crate::database::{IndexedDatabase, PeptideIx};
use crate::lfq::{PrecursorId, QuantifiedPeak};
use crate::mass::{H2O, NH3};
use crate::peptide::{ModificationKind, Peptide, Site};
use crate::scoring::Feature;
use fnv::FnvHashMap;
use serde::{Deserialize, Serialize};
use std::collections::HashMap;
use std::hash::BuildHasher;

/// Maximum spectrum q-value for a PSM to be a product or a parent.
pub const PSM_Q_VALUE: f32 = 0.01;
/// Maximum peptide q-value for a PSM to be a product or a parent.
pub const PEPTIDE_Q_VALUE: f32 = 0.01;
/// Window in daltons for matching a mass difference to a loss. Peptide and
/// modification masses are `f32`, so a difference around 2,000 Da carries
/// about 1e-4 Da of representation error.
pub const LOSS_EPSILON: f32 = 0.005;
/// Percentile of a peptide's retention-time spread taken as the co-elution
/// tolerance of the retention-time mode. The median, not a high percentile:
/// the window decides how much of the flagged set is coincidence, and the
/// chance estimate falls with the window and nothing else (see
/// [`BACKGROUND_BAND`]). On a semi-enzymatic search of one PXD028735 file the
/// median is 0.255 min, which admits 265 of 4,168 candidate pairs with 36.6
/// expected by chance (13.8%) where no quantification runs; a window of 1 min
/// admits 308 of the quantified file's 2,513 with 90.2 expected (29.3%).
pub const RT_PERCENTILE: f32 = 50.0;
/// Fewest spread observations needed to derive a tolerance from a file.
pub const MIN_SPREAD_OBSERVATIONS: usize = 20;
/// Fallback co-elution tolerance, as a fraction of a file's retention-time
/// span, when the file has too few repeated peptides to derive one.
pub const FALLBACK_RT_FRACTION: f32 = 0.002;
/// Fallback co-elution tolerance when a file's retention-time span is unusable.
pub const FALLBACK_RT_TOLERANCE: f32 = 0.1;
/// Fraction of a peak width within which two apexes are the same peak. An
/// apex is a property of the peak, so the window is a fraction of the peak's
/// own full width at half maximum rather than a fixed time: see
/// [`apex_tolerance`]. On a semi-enzymatic search of one PXD028735 file, where
/// the median peak is 0.26 min wide, a quarter of a width admits 153 of the
/// 2,513 candidate pairs the apexes judge, with 6.2 expected by chance (4.0%),
/// while the retention times of the same pairs cannot do better than 12.3% at
/// any window. Narrowing to a tenth of a width takes the chance estimate to
/// 2.9% and gives up 69 of those pairs, of which chance explains 3.7.
pub const APEX_FWHM_FRACTION: f32 = 0.25;
/// Width of a retention-time difference histogram bin, in the file's
/// retention-time units. The grid is fixed rather than a function of the
/// tolerance, so the distribution can be read, and a tolerance chosen, after
/// the fact.
pub const HISTOGRAM_BIN: f32 = 0.05;
/// Half-width of the histogram, in the file's retention-time units. Pairs
/// further out are counted in [`InSourcePlacementStats::outside_pairs`].
pub const HISTOGRAM_HALF_WIDTH: f32 = 5.0;
/// Retention-time differences, in the file's units, whose pair density
/// estimates the rate at which unrelated pairs fall inside the co-elution
/// window. Products co-elute, so a band this far from zero holds pairs that
/// are related in sequence and not in time.
pub const BACKGROUND_BAND: std::ops::Range<f32> = 1.0..HISTOGRAM_HALF_WIDTH;

/// `in_source` configuration: `true` or `false` to override the automatic
/// choice, or an object to set the co-elution tolerance.
#[derive(Clone, Copy, Debug, PartialEq, Deserialize, schemars::JsonSchema)]
#[serde(untagged)]
pub enum InSourceConfig {
    Enabled(bool),
    Settings(InSourceSettings),
}

impl InSourceConfig {
    /// Settings to use; `None` when the annotation is off. `automatic` is what
    /// an unset key means: the annotation is on for a semi-enzymatic or
    /// non-specific digest, whose own shorter peptides are exactly what
    /// in-source fragmentation produces, and off otherwise.
    pub fn resolve(config: Option<Self>, automatic: bool) -> Option<InSourceSettings> {
        match config {
            None => automatic.then(InSourceSettings::default),
            Some(Self::Enabled(false)) => None,
            Some(Self::Enabled(true)) => Some(InSourceSettings::default()),
            Some(Self::Settings(settings)) => Some(settings),
        }
    }
}

/// In-source annotation settings.
#[derive(Clone, Copy, Debug, Default, PartialEq, Serialize, Deserialize, schemars::JsonSchema)]
#[serde(default, deny_unknown_fields)]
pub struct InSourceSettings {
    /// Co-elution tolerance in the file's retention-time units. When unset,
    /// each file's tolerance is derived from the spread of its own repeated
    /// peptides; see [`derive_rt_tolerance`].
    #[serde(skip_serializing_if = "Option::is_none")]
    #[schemars(extend("exclusiveMinimum" = 0.0))]
    pub rt_tolerance: Option<f32>,
    /// Fraction of a peak's full width at half maximum within which the
    /// product's and the parent's quantified apexes are the same peak
    /// (default 0.15, [`APEX_FWHM_FRACTION`]). Used only where quantification
    /// gives both of them a peak in the file.
    #[serde(skip_serializing_if = "Option::is_none")]
    #[schemars(extend("exclusiveMinimum" = 0.0))]
    pub apex_fwhm_fraction: Option<f32>,
}

impl InSourceSettings {
    /// Why these settings cannot be used, if so.
    pub fn validate(&self) -> Result<(), String> {
        for (name, value) in [
            ("rt_tolerance", self.rt_tolerance),
            ("apex_fwhm_fraction", self.apex_fwhm_fraction),
        ] {
            if let Some(value) = value {
                if !(value.is_finite() && value > 0.0) {
                    return Err(format!(
                        "`in_source.{name}` must be a finite value greater than zero"
                    ));
                }
            }
        }
        Ok(())
    }

    /// The fraction of a peak width within which two apexes are one peak.
    pub fn apex_tolerance(&self) -> f32 {
        self.apex_fwhm_fraction.unwrap_or(APEX_FWHM_FRACTION)
    }
}

/// Where a product sits inside its parent, which fixes the mass it must have.
///
/// Backbone cleavage in the source gives the same two ion types as cleavage in
/// the collision cell. The C-terminal piece keeps the parent's C-terminal OH
/// and gains a hydrogen at its new N-terminus, so it has the ordinary peptide
/// mass of its own sequence, like a y ion. The N-terminal piece does not: it
/// is a b-type species, 18.010565 Da (water) lighter than the peptide of the
/// same sequence. An internal piece has both of those new termini, so it too
/// is water-light. A peptide at the ordinary peptide mass of an N-terminal or
/// internal span is therefore not an in-source product of that parent, and is
/// not flagged.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub enum Placement {
    /// The product covers the parent's own span: only a loss can make it
    /// lighter.
    SameSpan,
    /// The product ends where the parent ends (a y-type product).
    Suffix,
    /// The product starts where the parent starts (a b-type product).
    Prefix,
    /// The product shares neither terminus with the parent.
    Internal,
}

impl Placement {
    /// Mass the product lacks relative to the ordinary peptide mass of its
    /// span of the parent, before any neutral loss.
    pub fn structural_deficit(self) -> f32 {
        match self {
            Placement::SameSpan | Placement::Suffix => 0.0,
            Placement::Prefix | Placement::Internal => H2O,
        }
    }

    /// The name used in `in_source_type` and in the run summary.
    pub fn label(self) -> &'static str {
        match self {
            Placement::SameSpan => "same_span",
            Placement::Suffix => "suffix",
            Placement::Prefix => "prefix",
            Placement::Internal => "internal",
        }
    }
}

/// What kind of in-source product a PSM is: where it sits in its parent, and
/// the neutral loss it carries on top of that, if any.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct InSourceKind {
    pub placement: Placement,
    /// Neutral mass lost on top of the placement's own deficit.
    pub loss: Option<f32>,
}

impl InSourceKind {
    /// The value written to the `in_source_type` column. A loss is written as
    /// its mass to four decimal places, because a mass is unambiguous where a
    /// name is not: 18.0106 is water however it was configured, and a loss
    /// declared by a modification cannot be confused with the modification.
    pub fn label(&self) -> String {
        match self.loss {
            None => format!("fragment:{}", self.placement.label()),
            Some(loss) if self.placement == Placement::SameSpan => format!("loss:{loss:.4}"),
            Some(loss) => format!("fragment:{}+loss:{loss:.4}", self.placement.label()),
        }
    }
}

/// The parent peptide a PSM's ion came from, and how.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct InSourceAnnotation {
    pub parent: PeptideIx,
    pub kind: InSourceKind,
}

/// A product/parent pair that passed the sequence and mass test. Whether the
/// two co-elute is a per-file question, asked separately for a PSM (against
/// the parent's PSMs) and for a quantified row (against the parent's own peak
/// in that file), so a candidate itself is file-independent.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct InSourceCandidate {
    pub parent: PeptideIx,
    pub kind: InSourceKind,
    /// Residues of the parent; the first tie-break for the longest parent.
    pub parent_len: usize,
}

/// A quantified peak of one precursor in one file, as far as co-elution is
/// concerned: where it apexes and how wide it is.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PeakShape {
    pub apex_rt: f32,
    /// Integration bounds of the peak, in the file's retention-time units.
    pub start_rt: f32,
    pub end_rt: f32,
    /// Full width at half maximum, where the trace fell to half height.
    pub fwhm: Option<f32>,
}

impl PeakShape {
    /// Whether a retention time lies inside the integrated peak. A PSM
    /// outside it fired in some other elution and is not this peak's PSM.
    pub fn contains(&self, rt: f32) -> bool {
        rt >= self.start_rt && rt <= self.end_rt
    }
}

/// The quantified peaks of a run, by peptide and file. An apex is a property
/// of the peak; a PSM's retention time is where MS2 happened to fire inside
/// that peak, which is why apexes make the better co-elution test when
/// quantification has provided them.
#[derive(Clone, Debug, Default)]
pub struct PeakApexes {
    peaks: FnvHashMap<(PeptideIx, usize), Vec<PeakShape>>,
    rows: FnvHashMap<(PrecursorId, usize), PeakShape>,
    median_fwhm: FnvHashMap<usize, f32>,
}

impl PeakApexes {
    /// Collect the target precursors' peaks from an LFQ result. Decoy rows are
    /// left out: they are never annotated.
    pub fn from_areas<H: BuildHasher>(
        areas: &HashMap<(PrecursorId, bool), QuantifiedPeak, H>,
    ) -> Self {
        let mut peaks: FnvHashMap<(PeptideIx, usize), Vec<PeakShape>> = FnvHashMap::default();
        let mut rows: FnvHashMap<(PrecursorId, usize), PeakShape> = FnvHashMap::default();
        let mut widths: FnvHashMap<usize, Vec<f32>> = FnvHashMap::default();
        for ((id, decoy), peak) in areas.iter() {
            if *decoy {
                continue;
            }
            let peptide_idx = match id {
                PrecursorId::Combined(idx) | PrecursorId::Charged((idx, _)) => *idx,
            };
            for (file_id, evidence) in peak.file_evidence.iter().enumerate() {
                let Some(evidence) = evidence.as_ref() else {
                    continue;
                };
                if !evidence.apex_rt.is_finite() {
                    continue;
                }
                let fwhm = evidence.fwhm.filter(|fwhm| fwhm.is_finite() && *fwhm > 0.0);
                if let Some(fwhm) = fwhm {
                    widths.entry(file_id).or_default().push(fwhm);
                }
                let shape = PeakShape {
                    apex_rt: evidence.apex_rt,
                    start_rt: evidence.peak_start_rt.min(evidence.peak_end_rt),
                    end_rt: evidence.peak_start_rt.max(evidence.peak_end_rt),
                    fwhm,
                };
                peaks.entry((peptide_idx, file_id)).or_default().push(shape);
                rows.insert((*id, file_id), shape);
            }
        }
        let median_fwhm = widths
            .into_iter()
            .map(|(file_id, mut values)| (file_id, median(&mut values)))
            .collect();
        Self {
            peaks,
            rows,
            median_fwhm,
        }
    }

    /// The peak of one quantified row: that precursor, in that file. A row
    /// without a peak there has none.
    pub fn row(&self, id: &PrecursorId, file_id: usize) -> Option<&PeakShape> {
        self.rows.get(&(*id, file_id))
    }

    /// The peaks of a precursor in one file, in no particular order.
    pub fn peaks(&self, peptide: PeptideIx, file_id: usize) -> &[PeakShape] {
        self.peaks
            .get(&(peptide, file_id))
            .map(Vec::as_slice)
            .unwrap_or_default()
    }

    /// Median full width at half maximum of the file's quantified peaks, used
    /// where a peak's own width is unknown.
    pub fn median_fwhm(&self, file_id: usize) -> Option<f32> {
        self.median_fwhm.get(&file_id).copied()
    }
}

/// One retention-time difference bin, [`HISTOGRAM_BIN`] wide.
#[derive(Clone, Copy, Debug, PartialEq, Serialize, Deserialize, schemars::JsonSchema)]
pub struct HistogramBin {
    /// Bin index: the bin is centred on `bin * HISTOGRAM_BIN` and is one
    /// [`HISTOGRAM_BIN`] wide, so bin 0 is centred on zero and bin `-k`
    /// mirrors bin `k`.
    pub bin: i32,
    /// Lower edge of the bin, in the mode's own unit.
    pub start: f32,
    /// Upper edge of the bin.
    pub end: f32,
    pub pairs: usize,
}

/// The built-in check for one placement within one file: how many candidate
/// pairs co-eluted, and how many would be expected to by chance.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize, schemars::JsonSchema)]
pub struct InSourcePlacementStats {
    /// `same_span`, `suffix`, `prefix`, `internal`, or `all`.
    pub placement: String,
    /// Product/parent pairs that passed the sequence and mass test.
    pub candidate_pairs: usize,
    /// Candidate pairs inside the co-elution window.
    pub coeluting_pairs: usize,
    /// Co-eluting pairs expected from the density of candidate pairs
    /// [`BACKGROUND_BAND`] away from zero, scaled to the window's width.
    pub expected_by_chance: f32,
    /// PSMs flagged with this placement.
    pub flagged_psms: usize,
    /// Candidate pairs further from zero than [`HISTOGRAM_HALF_WIDTH`].
    pub outside_pairs: usize,
    /// Retention-time differences of the candidate pairs, on the fixed
    /// [`HISTOGRAM_BIN`] grid. Empty bins are left out.
    pub rt_difference_histogram: Vec<HistogramBin>,
}

impl InSourcePlacementStats {
    /// Sum the same placement's records from several files into one. The
    /// histogram grid is fixed, so the bins add; each file contributes its own
    /// chance estimate, computed with its own tolerance.
    pub fn merge<'a>(label: &str, parts: impl Iterator<Item = &'a Self>) -> Self {
        let mut merged = Self {
            placement: label.into(),
            candidate_pairs: 0,
            coeluting_pairs: 0,
            expected_by_chance: 0.0,
            flagged_psms: 0,
            outside_pairs: 0,
            rt_difference_histogram: Vec::new(),
        };
        let mut bins: FnvHashMap<i32, usize> = FnvHashMap::default();
        for part in parts {
            merged.candidate_pairs += part.candidate_pairs;
            merged.coeluting_pairs += part.coeluting_pairs;
            merged.expected_by_chance += part.expected_by_chance;
            merged.flagged_psms += part.flagged_psms;
            merged.outside_pairs += part.outside_pairs;
            for entry in &part.rt_difference_histogram {
                *bins.entry(entry.bin).or_default() += entry.pairs;
            }
        }
        merged.rt_difference_histogram = histogram_bins(bins.into_iter());
        merged
    }
}

/// The built-in check for one way of deciding co-elution.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize, schemars::JsonSchema)]
pub struct InSourceModeStats {
    /// `apex` where quantification gave both the product and the parent a
    /// peak in the file, `rt` where the product's PSM retention time was
    /// compared with the parent's.
    pub mode: String,
    /// `peak_fwhm` for `apex`, where a difference is measured in widths of
    /// the narrower peak, or `retention_time` for `rt`.
    pub unit: String,
    /// Tolerance in that unit; absent at the run level when the files used
    /// different ones.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub tolerance: Option<f32>,
    /// `configured`, `derived`, `fallback` or `peak_width`; absent at the run
    /// level when the files differ.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub tolerance_source: Option<String>,
    /// The check for this mode as a whole, then per placement.
    pub totals: InSourcePlacementStats,
    pub placements: Vec<InSourcePlacementStats>,
}

/// The built-in check, per file.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize, schemars::JsonSchema)]
pub struct InSourceFileStats {
    pub file_id: usize,
    /// Confident target PSMs that could be a product or a parent.
    pub confident_psms: usize,
    /// PSMs flagged in this file.
    pub flagged_psms: usize,
    /// Median full width at half maximum of the file's quantified peaks, the
    /// fallback width when a peak's own is unknown.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub median_fwhm: Option<f32>,
    /// One record per way of deciding co-elution that was used.
    pub modes: Vec<InSourceModeStats>,
}

/// Per-PSM in-source annotations of a whole run.
#[derive(Clone, Debug, Default)]
pub struct InSourceAnnotations {
    by_psm: FnvHashMap<usize, InSourceAnnotation>,
    candidates: FnvHashMap<PeptideIx, Vec<InSourceCandidate>>,
    apex_fraction: f32,
    apexes: Option<PeakApexes>,
    files: Vec<InSourceFileStats>,
}

impl InSourceAnnotations {
    /// The annotation of a PSM, by `psm_id`.
    pub fn psm(&self, psm_id: usize) -> Option<&InSourceAnnotation> {
        self.by_psm.get(&psm_id)
    }

    /// The product/parent pairs of a peptide that passed the sequence and mass
    /// test, in any file. A caller with retention times of its own — a
    /// quantified row's apex and its parent's apex in the same file — asks
    /// [`InSourceAnnotations::coeluting`] which of them co-elute, so a row
    /// transferred between runs is judged on its peak and not on a PSM.
    pub fn candidates(&self, peptide: PeptideIx) -> &[InSourceCandidate] {
        self.candidates
            .get(&peptide)
            .map(Vec::as_slice)
            .unwrap_or_default()
    }

    /// The per-file record of the built-in check.
    pub fn files(&self) -> &[InSourceFileStats] {
        &self.files
    }

    /// Flagged PSMs across the run.
    pub fn flagged_psms(&self) -> usize {
        self.by_psm.len()
    }

    /// Whether any PSM was flagged.
    pub fn is_empty(&self) -> bool {
        self.by_psm.is_empty()
    }

    /// The built-in check for the whole run, one record per mode, with each
    /// mode's placements. A mode's tolerance is reported only where every
    /// file that used it agreed.
    pub fn run_modes(&self) -> Vec<InSourceModeStats> {
        let mut modes: Vec<&str> = self
            .files
            .iter()
            .flat_map(|file| file.modes.iter().map(|mode| mode.mode.as_str()))
            .collect();
        modes.sort_unstable();
        modes.dedup();
        modes
            .into_iter()
            .map(|mode| {
                let parts = self
                    .files
                    .iter()
                    .flat_map(|file| file.modes.iter())
                    .filter(|stats| stats.mode == mode)
                    .collect::<Vec<_>>();
                let tolerances = parts
                    .iter()
                    .filter_map(|stats| stats.tolerance)
                    .collect::<Vec<_>>();
                let agreed = tolerances
                    .first()
                    .copied()
                    .filter(|first| tolerances.iter().all(|value| value == first));
                let mut placements: Vec<&str> = parts
                    .iter()
                    .flat_map(|stats| {
                        stats
                            .placements
                            .iter()
                            .map(|placement| placement.placement.as_str())
                    })
                    .collect();
                placements.sort_unstable_by_key(|label| placement_order(label));
                placements.dedup();
                InSourceModeStats {
                    mode: mode.into(),
                    unit: parts
                        .first()
                        .map(|stats| stats.unit.clone())
                        .unwrap_or_default(),
                    tolerance: agreed,
                    tolerance_source: agreed.and_then(|_| {
                        parts
                            .first()
                            .and_then(|stats| stats.tolerance_source.clone())
                    }),
                    totals: InSourcePlacementStats::merge(
                        "all",
                        parts.iter().map(|stats| &stats.totals),
                    ),
                    placements: placements
                        .into_iter()
                        .map(|label| {
                            InSourcePlacementStats::merge(
                                label,
                                parts
                                    .iter()
                                    .flat_map(|stats| stats.placements.iter())
                                    .filter(|placement| placement.placement == label),
                            )
                        })
                        .collect(),
                }
            })
            .collect()
    }

    /// The annotation of one quantified row — that precursor, in that file —
    /// judged on the peak quantification gave it against the parent's peaks in
    /// the same file. A row flag asks less of the parent than a PSM flag does:
    /// the parent must be confident somewhere in the run, which every
    /// candidate's parent is, and quantified in this file, but it need not
    /// have a PSM in this file, so a row transferred between runs is judged
    /// too. `None` without quantification, where the row or the parent has no
    /// peak in the file, or where nothing co-elutes. Decoy rows are not in the
    /// peak table, so they are never annotated.
    pub fn precursor_row(&self, id: &PrecursorId, file_id: usize) -> Option<InSourceAnnotation> {
        let apexes = self.apexes.as_ref()?;
        let row = std::slice::from_ref(apexes.row(id, file_id)?);
        let file_fwhm = apexes.median_fwhm(file_id);
        let peptide = match id {
            PrecursorId::Combined(idx) | PrecursorId::Charged((idx, _)) => *idx,
        };
        self.candidates(peptide)
            .iter()
            .filter(|candidate| {
                closest_apex_ratio(row, apexes.peaks(candidate.parent, file_id), file_fwhm)
                    .is_some_and(|(ratio, _)| ratio.abs() < self.apex_fraction)
            })
            .reduce(preferred_parent)
            .map(|candidate| InSourceAnnotation {
                parent: candidate.parent,
                kind: candidate.kind,
            })
    }
}

/// Of every pairing of a product's and a parent's peaks in one file, the one
/// whose apexes agree best: the signed difference in widths of the narrower
/// peak, and the product peak it was measured on. `None` when either side has
/// no peak, or no width is known.
fn closest_apex_ratio<'a>(
    product: &'a [PeakShape],
    parent: &[PeakShape],
    file_fwhm: Option<f32>,
) -> Option<(f32, &'a PeakShape)> {
    product
        .iter()
        .flat_map(|product| {
            parent.iter().filter_map(move |parent| {
                let width = match (product.fwhm, parent.fwhm) {
                    (Some(left), Some(right)) => left.min(right),
                    (Some(width), None) | (None, Some(width)) => width,
                    (None, None) => file_fwhm?,
                };
                (width > 0.0).then(|| ((product.apex_rt - parent.apex_rt) / width, product))
            })
        })
        .reduce(|best, candidate| {
            if candidate.0.abs() < best.0.abs() {
                candidate
            } else {
                best
            }
        })
}

/// The histogram entries of a bin count, sorted by bin. Each bin is one
/// [`HISTOGRAM_BIN`] wide and centred on `bin * HISTOGRAM_BIN`, so bin `-k`
/// mirrors bin `k`.
fn histogram_bins(bins: impl Iterator<Item = (i32, usize)>) -> Vec<HistogramBin> {
    let mut entries = bins
        .map(|(bin, pairs)| HistogramBin {
            bin,
            start: (bin as f32 - 0.5) * HISTOGRAM_BIN,
            end: (bin as f32 + 0.5) * HISTOGRAM_BIN,
            pairs,
        })
        .collect::<Vec<_>>();
    entries.sort_unstable_by_key(|entry| entry.bin);
    entries
}

/// Sort key putting placement records in the order [`Placement`] declares.
fn placement_order(label: &str) -> usize {
    [
        Placement::SameSpan,
        Placement::Suffix,
        Placement::Prefix,
        Placement::Internal,
    ]
    .iter()
    .position(|placement| placement.label() == label)
    .unwrap_or(usize::MAX)
}

/// The parent to report of two candidates: the longer peptide, then the lower
/// peptide index, so the choice never depends on iteration order.
fn preferred_parent<'a>(
    best: &'a InSourceCandidate,
    candidate: &'a InSourceCandidate,
) -> &'a InSourceCandidate {
    if (candidate.parent_len, std::cmp::Reverse(candidate.parent.0))
        > (best.parent_len, std::cmp::Reverse(best.parent.0))
    {
        candidate
    } else {
        best
    }
}

/// A site where a product and its span of the parent carry different
/// modification mass.
#[derive(Clone, Copy, Debug, PartialEq)]
struct SiteDifference {
    /// The site within the parent, where a declared loss is looked up.
    site: Site,
    /// Product mass minus parent mass at that site.
    difference: f32,
}

/// Modification differences between a product and its span of a parent,
/// compared site by site, so a modification moved to another residue is a
/// difference rather than an equal total. `Ok(None)` is a product whose
/// modifications match the parent's everywhere; `Err(())` is more than one
/// differing site, which one neutral loss cannot explain.
#[allow(clippy::result_unit_err)]
fn site_differences(
    product: &Peptide,
    parent: &Peptide,
    offset: usize,
) -> Result<Option<SiteDifference>, ()> {
    let length = product.sequence.len();
    let mut found: Option<SiteDifference> = None;
    let mut record = |site: Site, difference: f32| -> bool {
        if difference.abs() <= LOSS_EPSILON {
            return true;
        }
        if found.is_some() {
            return false;
        }
        found = Some(SiteDifference { site, difference });
        true
    };
    for index in 0..length {
        let difference = product.modification_at(index) - parent.modification_at(index + offset);
        if !record(Site::Sequence((index + offset) as u32), difference) {
            return Err(());
        }
    }
    // A terminal modification of the parent is on the product only where the
    // product keeps that terminus.
    let parent_nterm = if offset == 0 {
        parent.nterm.unwrap_or(0.0)
    } else {
        0.0
    };
    let parent_cterm = if offset + length == parent.sequence.len() {
        parent.cterm.unwrap_or(0.0)
    } else {
        0.0
    };
    if !record(Site::Nterm, product.nterm.unwrap_or(0.0) - parent_nterm) {
        return Err(());
    }
    if !record(Site::Cterm, product.cterm.unwrap_or(0.0) - parent_cterm) {
        return Err(());
    }
    Ok(found)
}

/// Neutral losses a product of `parent` may carry: water and ammonia, every
/// loss configured under `database.fragment_losses`, and the losses the
/// parent's own modifications declare at the site in question. There is no
/// implicit whole-modification loss: a modification's full mass counts only
/// where the user declared it among that modification's neutral losses, and a
/// labelling modification contributes nothing, so the light channel of a
/// labelling pair is never a product of its heavy partner.
fn declared_losses(parent: &Peptide, database: &IndexedDatabase, site: Option<Site>) -> Vec<f32> {
    let mut losses = vec![H2O, NH3];
    if let Some(configured) = database.fragment_losses.as_ref() {
        losses.extend(configured.losses.iter().map(|loss| loss.mass));
    }
    if let Some(site) = site {
        let sequence = parent.sequence.as_bytes();
        for applied in parent.applied_modifications() {
            if applied.kind == ModificationKind::Label {
                continue;
            }
            losses.extend(
                applied
                    .modification
                    .losses_at(sequence, site)
                    .iter()
                    .copied()
                    .filter(|loss| *loss > 0.0),
            );
        }
    }
    losses
}

/// Classify a product against a parent whose residues it shares from `offset`,
/// or `None` when it is not a product of that parent.
fn classify(
    product: &Peptide,
    parent: &Peptide,
    offset: usize,
    placement: Placement,
    database: &IndexedDatabase,
) -> Option<InSourceKind> {
    let difference = site_differences(product, parent, offset).ok()?;
    let deficit = placement.structural_deficit();
    // Mass the product lacks against the ordinary peptide mass of its span of
    // the parent. A product is only ever lighter.
    let missing = -difference.map_or(0.0, |site| site.difference);
    if (missing - deficit).abs() <= LOSS_EPSILON {
        // The placement's own deficit and nothing else. The same span with the
        // same modifications is the same peptidoform, not a product of it.
        return (placement != Placement::SameSpan).then_some(InSourceKind {
            placement,
            loss: None,
        });
    }
    let lost = missing - deficit;
    if lost <= LOSS_EPSILON {
        return None;
    }
    let losses = declared_losses(parent, database, difference.map(|site| site.site));
    let matched = losses
        .iter()
        .copied()
        .filter(|loss| (loss - lost).abs() <= LOSS_EPSILON)
        .min_by(|left, right| (left - lost).abs().total_cmp(&(right - lost).abs()))?;
    Some(InSourceKind {
        placement,
        loss: Some(matched),
    })
}

/// Co-elution tolerance of one file: the [`RT_PERCENTILE`]th percentile (the
/// median) of how far each of two PSMs of one peptide sits from the centre of
/// the pair, over every pair of PSMs of every peptide seen more than once in
/// the file. Different charge
/// states of one peptide count as repeats, so the estimate needs neither MS1
/// data nor quantification. A file with fewer than
/// [`MIN_SPREAD_OBSERVATIONS`] such pairs, or whose repeats all sit at one
/// retention time, falls back to [`FALLBACK_RT_FRACTION`] of its
/// retention-time span.
pub fn derive_rt_tolerance(spreads: &mut [f32], rt_span: f32) -> (f32, &'static str) {
    if spreads.len() >= MIN_SPREAD_OBSERVATIONS {
        spreads.sort_unstable_by(f32::total_cmp);
        let rank = (RT_PERCENTILE / 100.0 * (spreads.len() - 1) as f32).round() as usize;
        let tolerance = spreads[rank.min(spreads.len() - 1)];
        if tolerance > 0.0 {
            return (tolerance, "derived");
        }
    }
    let fallback = rt_span * FALLBACK_RT_FRACTION;
    if fallback > 0.0 {
        (fallback, "fallback")
    } else {
        (FALLBACK_RT_TOLERANCE, "fallback")
    }
}

/// One confident PSM taking part in the annotation.
struct ConfidentPsm {
    psm_id: usize,
    peptide_idx: PeptideIx,
    rt: f32,
}

/// A confident peptide's occurrence in a protein.
struct Occurrence {
    peptide_idx: PeptideIx,
    start: u32,
    end: u32,
    len: usize,
}

fn is_confident(feature: &Feature, database: &IndexedDatabase) -> bool {
    feature.label == 1
        && feature.spectrum_q <= PSM_Q_VALUE
        && feature.peptide_q <= PEPTIDE_Q_VALUE
        && !database[feature.peptide_idx].decoy
}

fn median(values: &mut [f32]) -> f32 {
    values.sort_unstable_by(f32::total_cmp);
    let middle = values.len() / 2;
    if values.len().is_multiple_of(2) {
        (values[middle - 1] + values[middle]) / 2.0
    } else {
        values[middle]
    }
}

/// The bin of a difference on the fixed [`HISTOGRAM_BIN`] grid: the nearest
/// multiple of the bin width, so every bin is one bin wide and bin `-k`
/// mirrors bin `k`. `None` when the difference lies further out than
/// [`HISTOGRAM_HALF_WIDTH`].
fn histogram_bin(difference: f32) -> Option<i32> {
    if !difference.is_finite() || difference.abs() >= HISTOGRAM_HALF_WIDTH {
        return None;
    }
    Some((difference / HISTOGRAM_BIN).round() as i32)
}

/// Candidate pairs per retention-time unit in [`BACKGROUND_BAND`], either side
/// of zero, from which the chance rate inside the window is estimated.
fn background_density(histogram: &FnvHashMap<i32, usize>) -> f32 {
    let low = (BACKGROUND_BAND.start / HISTOGRAM_BIN).round() as i32;
    let high = (BACKGROUND_BAND.end / HISTOGRAM_BIN).round() as i32;
    let pairs: usize = histogram
        .iter()
        .filter(|(bin, _)| {
            let distance = bin.saturating_abs();
            (low..high).contains(&distance)
        })
        .map(|(_, pairs)| *pairs)
        .sum();
    let width = 2.0 * (BACKGROUND_BAND.end - BACKGROUND_BAND.start);
    if width > 0.0 {
        pairs as f32 / width
    } else {
        0.0
    }
}

/// The counts of one placement, built while a file is scanned.
#[derive(Default)]
struct PlacementCounts {
    candidate_pairs: usize,
    coeluting_pairs: usize,
    flagged_psms: usize,
    outside_pairs: usize,
    histogram: FnvHashMap<i32, usize>,
}

impl PlacementCounts {
    fn record(&mut self, difference: f32, tolerance: f32) {
        self.candidate_pairs += 1;
        if difference.abs() < tolerance {
            self.coeluting_pairs += 1;
        }
        match histogram_bin(difference) {
            Some(bin) => *self.histogram.entry(bin).or_default() += 1,
            None => self.outside_pairs += 1,
        }
    }

    fn finish(&self, placement: &str, tolerance: f32) -> InSourcePlacementStats {
        let rt_difference_histogram =
            histogram_bins(self.histogram.iter().map(|(bin, pairs)| (*bin, *pairs)));
        InSourcePlacementStats {
            placement: placement.into(),
            candidate_pairs: self.candidate_pairs,
            coeluting_pairs: self.coeluting_pairs,
            // The window is two tolerances wide, one either side of zero.
            expected_by_chance: background_density(&self.histogram) * 2.0 * tolerance,
            flagged_psms: self.flagged_psms,
            outside_pairs: self.outside_pairs,
            rt_difference_histogram,
        }
    }
}

/// Annotate the in-source products among `features`, per file.
///
/// Only confident target PSMs ([`PSM_Q_VALUE`], [`PEPTIDE_Q_VALUE`]) can be a
/// product or a parent. Nothing in `features` is modified: the annotations are
/// a side table keyed by `psm_id`.
pub fn annotate(
    features: &[Feature],
    database: &IndexedDatabase,
    settings: &InSourceSettings,
    apexes: Option<PeakApexes>,
) -> InSourceAnnotations {
    let mut annotations = InSourceAnnotations {
        apex_fraction: settings.apex_tolerance(),
        apexes,
        ..InSourceAnnotations::default()
    };
    let mut by_file: FnvHashMap<usize, Vec<ConfidentPsm>> = FnvHashMap::default();
    for feature in features {
        if is_confident(feature, database) {
            by_file
                .entry(feature.file_id)
                .or_default()
                .push(ConfidentPsm {
                    psm_id: feature.psm_id,
                    peptide_idx: feature.peptide_idx,
                    rt: feature.rt,
                });
        }
    }

    // Every peptide confident in any file, in protein coordinates. A pair is
    // a candidate on sequence and mass alone, which is file-independent; a
    // parent still has to be confident in the file where the product is
    // flagged, which the per-file pass checks.
    let mut by_protein: FnvHashMap<&str, Vec<Occurrence>> = FnvHashMap::default();
    let mut indexed: Vec<PeptideIx> = by_file
        .values()
        .flat_map(|psms| psms.iter().map(|psm| psm.peptide_idx))
        .collect();
    indexed.sort_unstable();
    indexed.dedup();
    for peptide_idx in indexed {
        let peptide = &database[peptide_idx];
        let len = peptide.sequence.len();
        for site in peptide.protein_sites.iter() {
            let Some(start) = site.start else { continue };
            by_protein
                .entry(site.protein.as_ref())
                .or_default()
                .push(Occurrence {
                    peptide_idx,
                    start,
                    end: start + len as u32,
                    len,
                });
        }
    }

    let mut file_ids = by_file.keys().copied().collect::<Vec<_>>();
    file_ids.sort_unstable();
    for file_id in file_ids {
        annotate_file(
            file_id,
            &by_file[&file_id],
            database,
            settings,
            &by_protein,
            &mut annotations,
        );
    }
    annotations
}

fn annotate_file(
    file_id: usize,
    psms: &[ConfidentPsm],
    database: &IndexedDatabase,
    settings: &InSourceSettings,
    by_protein: &FnvHashMap<&str, Vec<Occurrence>>,
    annotations: &mut InSourceAnnotations,
) {
    // Retention times of every peptide in this file, which give both the
    // parents' retention times and the spread the tolerance comes from.
    let mut rts: FnvHashMap<PeptideIx, Vec<f32>> = FnvHashMap::default();
    let mut rt_min = f32::MAX;
    let mut rt_max = f32::MIN;
    for psm in psms {
        rts.entry(psm.peptide_idx).or_default().push(psm.rt);
        rt_min = rt_min.min(psm.rt);
        rt_max = rt_max.max(psm.rt);
    }
    let mut medians: FnvHashMap<PeptideIx, f32> = FnvHashMap::default();
    let mut spreads = Vec::new();
    for (peptide_idx, values) in rts.iter_mut() {
        medians.insert(*peptide_idx, median(values));
        // Half the distance between two PSMs of one peptide, each pair
        // counted once: how far each of the two sits from the centre of the
        // pair, which is the quantity the window applies to a PSM and its
        // parent's median. Distances from the peptide's own median would be
        // biased low instead: the median PSM is zero from itself by
        // construction, and a peptide seen twice contributes the same half
        // distance twice.
        for (index, left) in values.iter().enumerate() {
            spreads.extend(
                values[index + 1..]
                    .iter()
                    .map(|right| (left - right).abs() / 2.0),
            );
        }
    }
    let (rt_tolerance, rt_tolerance_source) = match settings.rt_tolerance {
        Some(configured) => (configured, "configured"),
        None => derive_rt_tolerance(&mut spreads, (rt_max - rt_min).max(0.0)),
    };

    // One set of counts per mode: a difference in peak widths is not a
    // difference in minutes, so the two are never pooled.
    let mut counts: FnvHashMap<(Mode, Placement), PlacementCounts> = FnvHashMap::default();
    let mut totals: FnvHashMap<Mode, PlacementCounts> = FnvHashMap::default();
    let mut flagged_psms = 0usize;
    let apex_fraction = settings.apex_tolerance();
    let file_fwhm = annotations
        .apexes
        .as_ref()
        .and_then(|apexes| apexes.median_fwhm(file_id));

    // A peptide's candidates depend only on the peptide, so the first file
    // that sees it pays for the scan and every PSM afterwards borrows them.
    for psm in psms {
        annotations
            .candidates
            .entry(psm.peptide_idx)
            .or_insert_with(|| candidates_of(psm.peptide_idx, database, by_protein));
    }
    let InSourceAnnotations {
        by_psm,
        candidates: by_peptide,
        apexes,
        ..
    } = annotations;

    for psm in psms {
        let candidates = by_peptide
            .get(&psm.peptide_idx)
            .map(Vec::as_slice)
            .unwrap_or_default();
        if candidates.is_empty() {
            continue;
        }
        let product_peaks = apexes
            .as_ref()
            .map(|apexes| apexes.peaks(psm.peptide_idx, file_id))
            .unwrap_or_default();
        let mut best: Option<(&InSourceCandidate, Mode)> = None;
        for candidate in candidates {
            // A parent is a peptide the file itself identified confidently.
            // One confident only in another file is no parent here, however
            // its peaks look.
            let Some(parent_rt) = medians.get(&candidate.parent) else {
                continue;
            };
            // Apexes decide where quantification found both peaks in this
            // file and this PSM fired inside the product's peak: an apex is a
            // property of the peak, while a PSM's retention time is only
            // where MS2 happened to fire inside it. A PSM outside that peak
            // belongs to another elution of the same peptide, so it is judged
            // on its own retention time, as it is in a search without
            // quantification.
            let parent_peaks = apexes
                .as_ref()
                .map(|apexes| apexes.peaks(candidate.parent, file_id))
                .unwrap_or_default();
            let apex = closest_apex_ratio(product_peaks, parent_peaks, file_fwhm)
                .filter(|(_, peak)| peak.contains(psm.rt));
            let (mode, difference, tolerance) = match apex {
                Some((ratio, _)) => (Mode::Apex, ratio, apex_fraction),
                None => (Mode::Rt, psm.rt - parent_rt, rt_tolerance),
            };
            counts
                .entry((mode, candidate.kind.placement))
                .or_default()
                .record(difference, tolerance);
            totals
                .entry(mode)
                .or_default()
                .record(difference, tolerance);
            if difference.abs() < tolerance {
                best = match best {
                    Some((incumbent, incumbent_mode)) => {
                        let winner = preferred_parent(incumbent, candidate);
                        if std::ptr::eq(winner, candidate) {
                            Some((candidate, mode))
                        } else {
                            Some((incumbent, incumbent_mode))
                        }
                    }
                    None => Some((candidate, mode)),
                };
            }
        }
        if let Some((candidate, mode)) = best {
            by_psm.insert(
                psm.psm_id,
                InSourceAnnotation {
                    parent: candidate.parent,
                    kind: candidate.kind,
                },
            );
            flagged_psms += 1;
            totals.entry(mode).or_default().flagged_psms += 1;
            counts
                .entry((mode, candidate.kind.placement))
                .or_default()
                .flagged_psms += 1;
        }
    }

    let mut modes = totals.keys().copied().collect::<Vec<_>>();
    modes.sort_unstable();
    let modes = modes
        .into_iter()
        .map(|mode| {
            let mut placements = counts
                .keys()
                .filter(|(counted, _)| *counted == mode)
                .map(|(_, placement)| *placement)
                .collect::<Vec<_>>();
            placements.sort_unstable();
            let (tolerance, source) = match mode {
                Mode::Apex => (apex_fraction, "peak_width"),
                Mode::Rt => (rt_tolerance, rt_tolerance_source),
            };
            InSourceModeStats {
                mode: mode.label().into(),
                unit: mode.unit().into(),
                tolerance: Some(tolerance),
                tolerance_source: Some(source.into()),
                totals: totals[&mode].finish("all", tolerance),
                placements: placements
                    .into_iter()
                    .map(|placement| {
                        counts[&(mode, placement)].finish(placement.label(), tolerance)
                    })
                    .collect(),
            }
        })
        .collect::<Vec<_>>();

    annotations.files.push(InSourceFileStats {
        file_id,
        confident_psms: psms.len(),
        flagged_psms,
        median_fwhm: file_fwhm,
        modes,
    });
}

/// How a product and a parent were compared.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, PartialOrd, Ord)]
enum Mode {
    /// Quantified apexes, in widths of the narrower peak.
    Apex,
    /// PSM retention times, in the file's retention-time units.
    Rt,
}

impl Mode {
    fn label(self) -> &'static str {
        match self {
            Mode::Apex => "apex",
            Mode::Rt => "rt",
        }
    }

    fn unit(self) -> &'static str {
        match self {
            Mode::Apex => "peak_fwhm",
            Mode::Rt => "retention_time",
        }
    }
}

/// Every parent of `peptide` that passes the sequence and mass test.
fn candidates_of(
    peptide: PeptideIx,
    database: &IndexedDatabase,
    by_protein: &FnvHashMap<&str, Vec<Occurrence>>,
) -> Vec<InSourceCandidate> {
    let product = &database[peptide];
    let length = product.sequence.len();
    // One candidate per parent. A parent that contains the product at more
    // than one offset — a repeated region — is classified at every offset and
    // keeps the placement [`Placement`] declares first, then the offset
    // nearest the parent's start, so the result never depends on the order of
    // [`Peptide::protein_sites`].
    let mut best: FnvHashMap<PeptideIx, (Placement, usize, InSourceCandidate)> =
        FnvHashMap::default();
    for site in product.protein_sites.iter() {
        let Some(product_start) = site.start else {
            continue;
        };
        let Some(occurrences) = by_protein.get(site.protein.as_ref()) else {
            continue;
        };
        let product_end = product_start + length as u32;
        for parent in occurrences {
            if parent.peptide_idx == peptide
                || parent.start > product_start
                || parent.end < product_end
            {
                continue;
            }
            let parent_peptide = &database[parent.peptide_idx];
            let offset = (product_start - parent.start) as usize;
            // Expanded ambiguous residues, and I/L twins when
            // `merge_isoleucine_leucine` is on, can give two peptides the same
            // protein span with different residues written, so the residues
            // are compared before the coordinates are trusted.
            if !crate::ambiguous_residues::isoleucine_leucine_eq(
                product.sequence.as_bytes(),
                &parent_peptide.sequence.as_bytes()[offset..offset + length],
            ) {
                continue;
            }
            let placement = match (offset == 0, offset + length == parent.len) {
                (true, true) => Placement::SameSpan,
                (false, true) => Placement::Suffix,
                (true, false) => Placement::Prefix,
                (false, false) => Placement::Internal,
            };
            let Some(kind) = classify(product, parent_peptide, offset, placement, database) else {
                continue;
            };
            let candidate = InSourceCandidate {
                parent: parent.peptide_idx,
                kind,
                parent_len: parent.len,
            };
            match best.entry(parent.peptide_idx) {
                std::collections::hash_map::Entry::Occupied(mut slot) => {
                    if (placement, offset) < (slot.get().0, slot.get().1) {
                        slot.insert((placement, offset, candidate));
                    }
                }
                std::collections::hash_map::Entry::Vacant(slot) => {
                    slot.insert((placement, offset, candidate));
                }
            }
        }
    }
    let mut candidates = best
        .into_values()
        .map(|(_, _, candidate)| candidate)
        .collect::<Vec<_>>();
    candidates.sort_unstable_by_key(|candidate| candidate.parent.0);
    candidates
}

#[cfg(test)]
#[path = "../tests/unit/in_source.rs"]
mod test;
