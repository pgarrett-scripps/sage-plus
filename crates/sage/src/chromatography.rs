//! Chromatography and charge of one file's confident identifications.
//!
//! Four numbers per file, all from data the run already holds: the charge-state
//! shares of its confident precursors, a chromatographic peak width, and, where
//! quantification ran, how many precursors it quantified and how many of those
//! had no direct MS2 evidence in that file.
//!
//! The peak width has two sources and is always labelled with the one used.
//! With label-free quantification the width is the median full width at half
//! maximum of the file's integrated peaks, a real peak width. Without it there
//! are no peaks to measure, so the width is the spread of repeat PSMs of one
//! peptide, the same quantity the in-source check derives its co-elution
//! tolerance from ([`crate::in_source::derive_rt_tolerance`]) and computed by
//! the same code. The two are not interchangeable: a half-distance between two
//! identifications is not a peak width, and nothing here pools them.

use crate::database::PeptideIx;
use crate::in_source::{derive_rt_tolerance, pairwise_half_distances};
use crate::lfq::{PrecursorId, QuantifiedPeak};
use fnv::FnvHashMap;
use std::collections::HashMap;
use std::hash::BuildHasher;

/// Charge states counted separately before the rest are pooled; charge 4 and
/// above share the last bucket, as the QC design asks.
pub const CHARGE_BUCKETS: usize = 4;

/// Which quantity a file's `peak_width` is.
#[derive(Clone, Copy, Debug, PartialEq, Eq, serde::Serialize, serde::Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum PeakWidth {
    /// Median full width at half maximum of the file's integrated LFQ peaks.
    LfqMedianFwhm,
    /// Median spread of repeat PSMs of one peptide, from the file's own
    /// repeats: half the distance between each pair, as the in-source
    /// co-elution tolerance is derived.
    RepeatPsmSpread,
    /// The same, where the file had too few repeats to derive one and fell back
    /// to a fraction of its retention-time span.
    RepeatPsmSpreadFallback,
}

impl PeakWidth {
    /// The value written to `qc.chromatography.tsv` and `run-summary.json`.
    pub fn label(self) -> &'static str {
        match self {
            Self::LfqMedianFwhm => "lfq_median_fwhm",
            Self::RepeatPsmSpread => "repeat_psm_spread",
            Self::RepeatPsmSpreadFallback => "repeat_psm_spread_fallback",
        }
    }
}

/// Charge and chromatography of one file.
#[derive(Clone, Debug, Default, PartialEq, serde::Serialize, serde::Deserialize)]
pub struct ChromatographyStats {
    /// Confident precursors counted: rank-1 target PSMs passing the run's QC
    /// q-value, counted as PSMs because a charge state is a property of a PSM.
    pub confident_precursors: u64,
    /// Precursors at charge 1, 2, 3, and 4 or above, in that order.
    pub charges: [u64; CHARGE_BUCKETS],
    /// Precursors whose charge the file did not state.
    pub charge_unknown: u64,
    /// Peak width in the file's retention-time unit, with the quantity it is.
    /// `None` when neither source had anything to measure.
    pub peak_width: Option<f32>,
    pub peak_width_kind: Option<PeakWidth>,
    /// Precursors with an integrated intensity in this file; `None` without
    /// quantification.
    pub quantified_precursors: Option<u64>,
    /// Of those, the ones with no direct MS2 evidence in this file, which is
    /// how match-between-runs carries a quantity into a file that did not
    /// identify it. `None` without quantification.
    pub transferred_precursors: Option<u64>,
}

impl ChromatographyStats {
    /// Percent of the counted precursors at `charge`, where charge 4 stands for
    /// 4 and above. Zero when nothing was counted.
    pub fn charge_pct(&self, charge: usize) -> f64 {
        let total = self.confident_precursors;
        match (total, self.charges.get(charge.saturating_sub(1))) {
            (0, _) | (_, None) => 0.0,
            (total, Some(&count)) => 100.0 * count as f64 / total as f64,
        }
    }

    /// Percent of the counted precursors whose charge was not stated.
    pub fn charge_unknown_pct(&self) -> f64 {
        if self.confident_precursors == 0 {
            return 0.0;
        }
        100.0 * self.charge_unknown as f64 / self.confident_precursors as f64
    }

    /// Percent of this file's quantified precursors that came from a transfer;
    /// `None` without quantification or with nothing quantified.
    pub fn transferred_pct(&self) -> Option<f64> {
        match (self.quantified_precursors?, self.transferred_precursors?) {
            (0, _) => None,
            (quantified, transferred) => Some(100.0 * transferred as f64 / quantified as f64),
        }
    }

    /// Count one confident PSM's charge.
    pub fn add_charge(&mut self, charge: u8) {
        self.confident_precursors += 1;
        match charge {
            0 => self.charge_unknown += 1,
            charge => {
                let bucket = usize::from(charge).min(CHARGE_BUCKETS) - 1;
                self.charges[bucket] += 1;
            }
        }
    }
}

/// Retention times of the confident PSMs of one peptide in one file.
pub type RepeatRts = FnvHashMap<PeptideIx, Vec<f32>>;

/// The repeat-PSM spread of one file, from the retention times of its confident
/// PSMs grouped by peptide: the same percentile of the same pairwise
/// half-distances the in-source co-elution tolerance uses.
///
/// `None` when the file has no PSM at all. A file with too few repeats falls
/// back to a fraction of its retention-time span, reported as
/// [`PeakWidth::RepeatPsmSpreadFallback`].
pub fn repeat_psm_spread(rts: &RepeatRts) -> Option<(f32, PeakWidth)> {
    let (mut low, mut high) = (f32::INFINITY, f32::NEG_INFINITY);
    let mut spreads = Vec::new();
    for values in rts.values() {
        for &rt in values {
            low = low.min(rt);
            high = high.max(rt);
        }
        pairwise_half_distances(values, &mut spreads);
    }
    if !low.is_finite() || !high.is_finite() {
        return None;
    }
    let (spread, source) = derive_rt_tolerance(&mut spreads, (high - low).max(0.0));
    Some(match source {
        "derived" => (spread, PeakWidth::RepeatPsmSpread),
        _ => (spread, PeakWidth::RepeatPsmSpreadFallback),
    })
}

/// Per-file quantification counts and median peak width from the LFQ areas.
///
/// Returns, for each file index, the precursors with an integrated intensity,
/// the ones among them with no direct MS2 evidence in that file, and the median
/// full width at half maximum of the file's peaks. Decoys are skipped.
pub fn quantification_by_file<H: BuildHasher>(
    areas: &HashMap<(PrecursorId, bool), QuantifiedPeak, H>,
    files: usize,
) -> Vec<(u64, u64, Option<f32>)> {
    let mut quantified = vec![0u64; files];
    let mut transferred = vec![0u64; files];
    let mut widths = vec![Vec::new(); files];
    for ((_, decoy), peak) in areas {
        if *decoy {
            continue;
        }
        for file_id in 0..files {
            if peak.intensities.get(file_id).copied().flatten().is_none() {
                continue;
            }
            quantified[file_id] += 1;
            let Some(Some(evidence)) = peak.file_evidence.get(file_id) else {
                continue;
            };
            if evidence.transfer_candidate {
                transferred[file_id] += 1;
            }
            if let Some(fwhm) = evidence.fwhm.filter(|fwhm| fwhm.is_finite() && *fwhm > 0.0) {
                widths[file_id].push(fwhm);
            }
        }
    }
    quantified
        .into_iter()
        .zip(transferred)
        .zip(widths)
        .map(|((quantified, transferred), mut values)| {
            (quantified, transferred, median(&mut values))
        })
        .collect()
}

/// Median of `values`, or `None` when empty.
fn median(values: &mut [f32]) -> Option<f32> {
    if values.is_empty() {
        return None;
    }
    values.sort_unstable_by(f32::total_cmp);
    let middle = values.len() / 2;
    Some(if values.len().is_multiple_of(2) {
        0.5 * (values[middle - 1] + values[middle])
    } else {
        values[middle]
    })
}

#[cfg(test)]
#[path = "../tests/unit/chromatography.rs"]
mod test;
