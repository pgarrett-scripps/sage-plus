//! In-source fragmentation and in-source neutral loss annotation.
//!
//! A peptide can reach the mass analyser as a piece of a longer peptide, or
//! having lost a small neutral molecule, because it fragmented in the ion
//! source rather than in the collision cell. Such a product is a real
//! spectrum of a real ion, so the search identifies it, but it is not an
//! experimental peptide: it is a piece of one that was in the sample.
//!
//! This module annotates those PSMs. It changes no score, removes no row and
//! reweights nothing: every annotated PSM keeps the identification it had.
//! Two things are reported per flagged PSM, the *root* peptide the product
//! came from and what kind of product it is, so a downstream analysis can
//! filter or group on them.
//!
//! A confident target PSM is flagged when both hold in the same file:
//!
//! 1. its peptide occupies a sub-span of a longer confident peptide of the
//!    same protein (a prefix, a suffix or an internal piece), or the same span
//!    with a lower mass, and the modifications on the shared residues are
//!    equal or differ by exactly one candidate loss; and
//! 2. its retention time co-elutes with the longer or heavier peptide's own
//!    PSMs, within a tolerance derived from the file itself.
//!
//! The flag names the root: the longest co-eluting peptide that contains the
//! product. A piece of a piece therefore points at the intact peptide, with no
//! chain to follow.
//!
//! The decision is per PSM, not per peptide. The same peptide can be an
//! experimental peptide at one retention time and an in-source product at
//! another in the same file, which is what happens when a labile modification
//! is lost in the source.
//!
//! Containment is decided in protein coordinates, from
//! [`Peptide::protein_sites`], not by searching one peptide's sequence inside
//! another: a product and its root then come from the same place in the same
//! protein, which is what in-source fragmentation means.

use crate::database::{IndexedDatabase, PeptideIx};
use crate::mass::{H2O, NH3};
use crate::peptide::Peptide;
use crate::scoring::Feature;
use fnv::FnvHashMap;
use serde::{Deserialize, Serialize};
use std::sync::Arc;

/// Maximum spectrum q-value for a PSM to be a product or a parent.
pub const PSM_Q_VALUE: f32 = 0.01;
/// Maximum peptide q-value for a PSM to be a product or a parent.
pub const PEPTIDE_Q_VALUE: f32 = 0.01;
/// Window in daltons for matching a modification-mass difference to a loss.
/// Peptide and modification masses are `f32`, so a difference around 2,000 Da
/// carries about 1e-4 Da of representation error.
pub const LOSS_EPSILON: f32 = 0.005;
/// Percentile of a peptide's retention-time spread taken as the co-elution
/// tolerance.
pub const RT_PERCENTILE: f32 = 95.0;
/// Fewest spread observations needed to derive a tolerance from a file.
pub const MIN_SPREAD_OBSERVATIONS: usize = 20;
/// Fallback co-elution tolerance, as a fraction of a file's retention-time
/// span, when the file has too few repeated peptides to derive one.
pub const FALLBACK_RT_FRACTION: f32 = 0.002;
/// Fallback co-elution tolerance when a file's retention-time span is unusable.
pub const FALLBACK_RT_TOLERANCE: f32 = 0.1;
/// Half-width of the retention-time difference histogram, in tolerances.
pub const HISTOGRAM_HALF_WIDTH_BINS: i32 = 25;
/// Distance from zero, in tolerances, of the histogram bins that estimate the
/// rate at which unrelated pairs fall inside the co-elution window.
pub const BACKGROUND_BIN_RANGE: std::ops::Range<i32> = 5..HISTOGRAM_HALF_WIDTH_BINS;

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
    /// non-specific digest, whose own pieces of longer peptides are exactly
    /// what in-source fragmentation produces, and off otherwise.
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
    #[schemars(range(min = 0.0))]
    pub rt_tolerance: Option<f32>,
}

/// What kind of in-source product a PSM is.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum InSourceKind {
    /// A piece of the root peptide, carrying the root's modifications on every
    /// shared residue.
    Fragment,
    /// The root's own span, having lost one neutral molecule.
    Loss(Arc<str>),
    /// A piece of the root peptide that also lost one neutral molecule.
    FragmentWithLoss(Arc<str>),
}

impl InSourceKind {
    /// The value written to the `in_source_type` column.
    pub fn label(&self) -> String {
        match self {
            InSourceKind::Fragment => "fragment".into(),
            InSourceKind::Loss(name) => format!("loss:{name}"),
            InSourceKind::FragmentWithLoss(name) => format!("fragment+loss:{name}"),
        }
    }
}

/// The root peptide a PSM's ion came from, and how.
#[derive(Clone, Debug, PartialEq)]
pub struct InSourceAnnotation {
    pub root: PeptideIx,
    pub kind: InSourceKind,
}

/// A structural product/parent pair within one file, before co-elution is
/// checked. Keeping these is what lets a quantified LFQ row be judged on its
/// own apex retention time instead of on a PSM's.
#[derive(Clone, Debug, PartialEq)]
struct ParentWindow {
    root: PeptideIx,
    kind: InSourceKind,
    /// Median retention time of the root's confident PSMs in this file.
    rt: f32,
    /// Residues of the root; the first tie-break for the longest root.
    root_len: usize,
}

/// One retention-time difference bin of a file's built-in check.
#[derive(Clone, Copy, Debug, PartialEq, Serialize, Deserialize, schemars::JsonSchema)]
pub struct HistogramBin {
    /// Bin index in tolerances from zero; bin 0 spans `[0, tolerance)` and
    /// bin -1 spans `[-tolerance, 0)`.
    pub bin: i32,
    /// Lower edge of the bin, in the file's retention-time units.
    pub start: f32,
    /// Upper edge of the bin.
    pub end: f32,
    pub pairs: usize,
}

/// The built-in check, per file: how many candidate pairs co-eluted, and how
/// many would be expected to by chance.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize, schemars::JsonSchema)]
pub struct InSourceFileStats {
    pub file_id: usize,
    /// Co-elution tolerance used for this file.
    pub rt_tolerance: f32,
    /// `configured`, `derived` or `fallback`.
    pub rt_tolerance_source: String,
    /// Confident target PSMs that could be a product or a parent.
    pub confident_psms: usize,
    /// Product/parent pairs that passed the sequence and mass test, before
    /// co-elution.
    pub candidate_pairs: usize,
    /// Candidate pairs inside the co-elution window.
    pub coeluting_pairs: usize,
    /// Co-eluting pairs expected from the rate at which candidate pairs fall
    /// `BACKGROUND_BIN_RANGE` tolerances away from zero. A co-eluting count
    /// far above this is the evidence that the flag finds in-source products
    /// rather than coincidences.
    pub expected_by_chance: f32,
    /// PSMs flagged in this file.
    pub flagged_psms: usize,
    /// Retention-time differences of every candidate pair, in bins one
    /// tolerance wide.
    pub rt_difference_histogram: Vec<HistogramBin>,
}

/// Per-PSM in-source annotations of a whole run.
#[derive(Clone, Debug, Default)]
pub struct InSourceAnnotations {
    by_psm: FnvHashMap<usize, InSourceAnnotation>,
    by_precursor: FnvHashMap<(PeptideIx, usize), Vec<ParentWindow>>,
    rt_tolerance: FnvHashMap<usize, f32>,
    files: Vec<InSourceFileStats>,
}

impl InSourceAnnotations {
    /// The annotation of a PSM, by `psm_id`.
    pub fn psm(&self, psm_id: usize) -> Option<&InSourceAnnotation> {
        self.by_psm.get(&psm_id)
    }

    /// The annotation of a quantified precursor in one file, judged on the
    /// retention time the caller measured (an LFQ peak's apex) rather than on
    /// any PSM's. The longest co-eluting root wins, as it does for a PSM.
    pub fn precursor(
        &self,
        peptide: PeptideIx,
        file_id: usize,
        rt: f32,
    ) -> Option<InSourceAnnotation> {
        let tolerance = *self.rt_tolerance.get(&file_id)?;
        let windows = self.by_precursor.get(&(peptide, file_id))?;
        let window = best_window(windows.iter().filter(|w| (w.rt - rt).abs() <= tolerance))?;
        Some(InSourceAnnotation {
            root: window.root,
            kind: window.kind.clone(),
        })
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
}

/// Pick the root among candidate windows: the longest peptide, then the lower
/// peptide index, so the choice never depends on iteration order.
fn best_window<'a>(windows: impl Iterator<Item = &'a ParentWindow>) -> Option<&'a ParentWindow> {
    windows.reduce(|best, window| {
        if (window.root_len, std::cmp::Reverse(window.root.0))
            > (best.root_len, std::cmp::Reverse(best.root.0))
        {
            window
        } else {
            best
        }
    })
}

/// A candidate loss: the name written to the output and the neutral mass.
#[derive(Clone, Debug, PartialEq)]
struct CandidateLoss {
    name: Arc<str>,
    mass: f32,
}

/// The losses a product of `parent` may carry: water and ammonia, the losses
/// configured under `database.fragment_losses`, the neutral losses the
/// parent's own modifications declare, and the full mass of each modification
/// the parent carries, which is what a product that lost a whole modification
/// differs by.
fn candidate_losses(parent: &Peptide, database: &IndexedDatabase) -> Vec<CandidateLoss> {
    let mut losses = vec![
        CandidateLoss {
            name: "water".into(),
            mass: H2O,
        },
        CandidateLoss {
            name: "ammonia".into(),
            mass: NH3,
        },
    ];
    if let Some(configured) = database.fragment_losses.as_ref() {
        for loss in &configured.losses {
            losses.push(CandidateLoss {
                name: loss.name.clone(),
                mass: loss.mass,
            });
        }
    }
    for applied in parent.applied_modifications() {
        let definition = applied.modification;
        let name: Arc<str> = definition
            .name
            .clone()
            .unwrap_or_else(|| format!("{:.4}", definition.mass).into());
        let declared = definition
            .neutral_losses
            .iter()
            .copied()
            .chain(
                definition
                    .site_losses
                    .iter()
                    .flat_map(|rules| rules.iter())
                    .flat_map(|(_, masses)| masses.iter().copied()),
            )
            .chain(std::iter::once(definition.mass));
        for mass in declared {
            if mass > 0.0 {
                losses.push(CandidateLoss {
                    name: name.clone(),
                    mass,
                });
            }
        }
    }
    losses
}

/// Total modification mass a product carries that its span of the parent does
/// not. Both peptides cover the same residues, so this difference is the whole
/// mass difference between the product and that span.
fn modification_delta(piece: &Peptide, parent: &Peptide, offset: usize) -> f32 {
    let piece_len = piece.sequence.len();
    let mut delta = piece.nterm.unwrap_or(0.0) + piece.cterm.unwrap_or(0.0);
    for index in 0..piece_len {
        delta += piece.modification_at(index) - parent.modification_at(index + offset);
    }
    if offset == 0 {
        delta -= parent.nterm.unwrap_or(0.0);
    }
    if offset + piece_len == parent.sequence.len() {
        delta -= parent.cterm.unwrap_or(0.0);
    }
    delta
}

/// Classify a product against a parent whose residues it shares from `offset`,
/// or `None` when the mass difference is not zero and not one candidate loss.
fn classify(
    piece: &Peptide,
    parent: &Peptide,
    offset: usize,
    losses: &[CandidateLoss],
) -> Option<InSourceKind> {
    let same_span = piece.sequence.len() == parent.sequence.len();
    let delta = modification_delta(piece, parent, offset);
    if delta.abs() <= LOSS_EPSILON {
        // The same span with the same modifications is the same peptidoform,
        // not a product of it.
        return (!same_span).then_some(InSourceKind::Fragment);
    }
    if delta > 0.0 {
        // A product cannot be heavier than the span of the parent it came from.
        return None;
    }
    let lost = -delta;
    let matched = losses
        .iter()
        .filter(|candidate| (candidate.mass - lost).abs() <= LOSS_EPSILON)
        .min_by(|left, right| {
            (left.mass - lost)
                .abs()
                .total_cmp(&(right.mass - lost).abs())
        })?;
    Some(if same_span {
        InSourceKind::Loss(matched.name.clone())
    } else {
        InSourceKind::FragmentWithLoss(matched.name.clone())
    })
}

/// Co-elution tolerance of one file: the [`RT_PERCENTILE`]th percentile of how
/// far a peptide's repeated PSMs sit from that peptide's median retention
/// time, over every peptide seen more than once in the file. Different charge
/// states of one peptide count as repeats, so the estimate needs neither MS1
/// data nor quantification. A file with fewer than
/// [`MIN_SPREAD_OBSERVATIONS`] repeats, or whose repeats all sit at one
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

/// A confident peptide's occurrence in a protein, within one file.
struct Occurrence {
    peptide_idx: PeptideIx,
    start: u32,
    end: u32,
    rt: f32,
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

/// Annotate the in-source products among `features`, per file.
///
/// Only confident target PSMs ([`PSM_Q_VALUE`], [`PEPTIDE_Q_VALUE`]) can be a
/// product or a parent. Nothing in `features` is modified: the annotations are
/// a side table keyed by `psm_id`.
pub fn annotate(
    features: &[Feature],
    database: &IndexedDatabase,
    settings: &InSourceSettings,
) -> InSourceAnnotations {
    let mut annotations = InSourceAnnotations::default();
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

    let mut file_ids = by_file.keys().copied().collect::<Vec<_>>();
    file_ids.sort_unstable();
    for file_id in file_ids {
        let psms = &by_file[&file_id];
        annotate_file(file_id, psms, database, settings, &mut annotations);
    }
    annotations
}

fn annotate_file(
    file_id: usize,
    psms: &[ConfidentPsm],
    database: &IndexedDatabase,
    settings: &InSourceSettings,
    annotations: &mut InSourceAnnotations,
) {
    // Retention times of every peptide seen in this file, which give both the
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
        let median_rt = median(values);
        medians.insert(*peptide_idx, median_rt);
        if values.len() > 1 {
            spreads.extend(values.iter().map(|rt| (rt - median_rt).abs()));
        }
    }
    let (rt_tolerance, rt_tolerance_source) = match settings.rt_tolerance {
        Some(configured) => (configured, "configured"),
        None => derive_rt_tolerance(&mut spreads, (rt_max - rt_min).max(0.0)),
    };
    annotations.rt_tolerance.insert(file_id, rt_tolerance);

    // Confident peptides by protein, in protein coordinates.
    let mut by_protein: FnvHashMap<Arc<str>, Vec<Occurrence>> = FnvHashMap::default();
    for (peptide_idx, median_rt) in &medians {
        let peptide = &database[*peptide_idx];
        let len = peptide.sequence.len();
        for site in peptide.protein_sites.iter() {
            let Some(start) = site.start else { continue };
            by_protein
                .entry(site.protein.clone())
                .or_default()
                .push(Occurrence {
                    peptide_idx: *peptide_idx,
                    start,
                    end: start + len as u32,
                    rt: *median_rt,
                    len,
                });
        }
    }

    // One parent is the parent of many products; its candidate losses depend
    // only on the parent, so they are built once.
    let mut losses_by_parent: FnvHashMap<PeptideIx, Vec<CandidateLoss>> = FnvHashMap::default();
    let mut histogram: FnvHashMap<i32, usize> = FnvHashMap::default();
    let mut candidate_pairs = 0usize;
    let mut coeluting_pairs = 0usize;
    let mut flagged_psms = 0usize;

    for psm in psms {
        let piece = &database[psm.peptide_idx];
        let piece_len = piece.sequence.len();
        let mut windows: Vec<ParentWindow> = Vec::new();
        for site in piece.protein_sites.iter() {
            let Some(piece_start) = site.start else {
                continue;
            };
            let Some(occurrences) = by_protein.get(&site.protein) else {
                continue;
            };
            let piece_end = piece_start + piece_len as u32;
            for parent in occurrences {
                if parent.peptide_idx == psm.peptide_idx
                    || parent.start > piece_start
                    || parent.end < piece_end
                {
                    continue;
                }
                let parent_peptide = &database[parent.peptide_idx];
                let offset = (piece_start - parent.start) as usize;
                // Expanded ambiguous residues can give two peptides the same
                // protein span with different residues; compare the residues
                // before trusting the coordinates.
                if piece.sequence.as_bytes()
                    != &parent_peptide.sequence.as_bytes()[offset..offset + piece_len]
                {
                    continue;
                }
                let losses = losses_by_parent
                    .entry(parent.peptide_idx)
                    .or_insert_with(|| candidate_losses(parent_peptide, database));
                let Some(kind) = classify(piece, parent_peptide, offset, losses) else {
                    continue;
                };
                if windows
                    .iter()
                    .any(|window| window.root == parent.peptide_idx)
                {
                    continue;
                }
                windows.push(ParentWindow {
                    root: parent.peptide_idx,
                    kind,
                    rt: parent.rt,
                    root_len: parent.len,
                });
            }
        }
        if windows.is_empty() {
            continue;
        }
        candidate_pairs += windows.len();
        for window in &windows {
            let difference = psm.rt - window.rt;
            let bin = (difference / rt_tolerance).floor() as i32;
            if bin.abs() <= HISTOGRAM_HALF_WIDTH_BINS {
                *histogram.entry(bin).or_default() += 1;
            }
            if difference.abs() <= rt_tolerance {
                coeluting_pairs += 1;
            }
        }
        if let Some(window) = best_window(
            windows
                .iter()
                .filter(|window| (psm.rt - window.rt).abs() <= rt_tolerance),
        ) {
            annotations.by_psm.insert(
                psm.psm_id,
                InSourceAnnotation {
                    root: window.root,
                    kind: window.kind.clone(),
                },
            );
            flagged_psms += 1;
        }
        // The windows depend only on the peptide, so the first PSM of a
        // peptide in this file records them for the quantified rows.
        annotations
            .by_precursor
            .entry((psm.peptide_idx, file_id))
            .or_insert(windows);
    }

    let background_bins = BACKGROUND_BIN_RANGE.len() * 2;
    let background_pairs: usize = histogram
        .iter()
        .filter(|(bin, _)| BACKGROUND_BIN_RANGE.contains(&bin.abs()))
        .map(|(_, pairs)| *pairs)
        .sum();
    // The co-elution window is two bins wide, one either side of zero.
    let expected_by_chance = if background_bins == 0 {
        0.0
    } else {
        2.0 * background_pairs as f32 / background_bins as f32
    };

    let mut rt_difference_histogram = histogram
        .into_iter()
        .map(|(bin, pairs)| HistogramBin {
            bin,
            start: bin as f32 * rt_tolerance,
            end: (bin + 1) as f32 * rt_tolerance,
            pairs,
        })
        .collect::<Vec<_>>();
    rt_difference_histogram.sort_unstable_by_key(|entry| entry.bin);

    annotations.files.push(InSourceFileStats {
        file_id,
        rt_tolerance,
        rt_tolerance_source: rt_tolerance_source.into(),
        confident_psms: psms.len(),
        candidate_pairs,
        coeluting_pairs,
        expected_by_chance,
        flagged_psms,
        rt_difference_histogram,
    });
}

#[cfg(test)]
#[path = "../tests/unit/in_source.rs"]
mod test;
