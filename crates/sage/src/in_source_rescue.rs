//! In-source rescue: a second pass that searches the MS2 spectra a run left
//! unidentified for in-source products of the peptides it did identify.
//!
//! [`crate::in_source`] labels in-source products that the database already
//! holds, which a fully enzymatic search mostly does not: a truncated tryptic
//! peptide is not tryptic. This pass builds those products from the confident
//! peptides of one file and searches them against that file's unidentified MS2
//! spectra, each candidate only where its parent elutes.
//!
//! The pass is separate from the main search in every way that matters for
//! error control. Every candidate has its own decoy, q-values are computed
//! inside the rescue set only, and a rescued PSM never enters the main
//! spectrum, peptide or protein FDR, protein inference or quantification.
//! Rescued PSMs are products, not sample peptides.
//!
//! This module holds the parts that depend only on peptides and numbers: the
//! candidates of a parent, the co-elution windows, and the q-values. Reading
//! spectra and scoring them is the runner's.

use crate::database::PeptideIx;
use crate::enzyme::{Position, ProteinOccurrence};
use crate::in_source::{InSourceKind, PeakApexes, Placement};
use crate::mass::{monoisotopic, H2O, NH3};
use crate::modification::ModificationDefinition;
use crate::peptide::{AppliedModification, CompactModifications, Peptide, Site};
use fnv::FnvHashMap;
use serde::{Deserialize, Serialize};
use std::sync::Arc;

/// Maximum rescue q-value for a rescued PSM to be counted as rescued.
pub const RESCUE_Q_VALUE: f32 = 0.01;

/// Hits the rescue scorer reports per spectrum. The best hit whose parent
/// co-elutes with the spectrum is kept, so a candidate that outscores it
/// outside its window cannot hide it.
pub const RESCUE_REPORT_PSMS: usize = 10;

/// Environment variable that shifts every co-elution window by this many of
/// its own widths. A test harness only, for the null control: a product cannot
/// elute several peak widths away from its parent, so with the windows shifted
/// the rescued count must fall to about the claimed error rate. It is not a
/// configuration key and is never set in a real search.
pub const NULL_SHIFT_ENV: &str = "SAGE_IN_SOURCE_RESCUE_NULL_SHIFT";

/// One product a parent can give in the source, as a peptide to search.
#[derive(Clone, Debug)]
pub struct RescueCandidate {
    pub peptide: Peptide,
    pub parent: PeptideIx,
    pub kind: InSourceKind,
    /// Residues of the parent; the longest co-eluting parent is reported.
    pub parent_len: usize,
}

/// The products a confident parent can give: every N-terminal truncation
/// (a C-terminal piece, y-type, at the ordinary peptide mass of its span) and
/// every C-terminal truncation (an N-terminal piece, b-type, water-light) down
/// to `min_len` residues, and the parent itself after losing water or ammonia.
///
/// A b-type piece is represented as its span carrying a C-terminal
/// -18.010565 Da, so its precursor mass and its y-series fragments are both
/// those of the b-type ion. A loss is placed on the N-terminus, where the
/// dominant in-source losses of water and ammonia happen (pyroglutamate from
/// an N-terminal E or Q); its precursor mass is exact wherever it happened,
/// and only the b-series of a loss elsewhere is placed wrongly. Internal
/// pieces, and losses declared by modifications, are not generated.
///
/// The modifications on the residues a piece keeps stay where they were; a
/// terminal modification of the parent stays only on a piece that keeps that
/// terminus. A parent whose sequence is not a view of a stored sequence, or
/// that is a decoy, gives nothing.
pub fn candidates_of(
    parent: &Peptide,
    parent_idx: PeptideIx,
    min_len: usize,
) -> Vec<RescueCandidate> {
    if parent.decoy {
        return Vec::new();
    }
    let mut out: Vec<(Peptide, InSourceKind)> = Vec::new();
    let len = parent.sequence.len();
    let min_len = min_len.max(1);
    let kind = |placement, loss| InSourceKind { placement, loss };
    // Suffixes: the piece ends where the parent ends.
    for start in 1..len {
        if len - start < min_len {
            break;
        }
        if let Some(peptide) = truncate(parent, start..len, None, None) {
            out.push((peptide, kind(Placement::Suffix, None)));
        }
    }
    // Prefixes: the piece starts where the parent starts, water-light.
    for end in (1..len).rev() {
        if end < min_len {
            break;
        }
        if let Some(peptide) = truncate(parent, 0..end, None, Some(-H2O)) {
            out.push((peptide, kind(Placement::Prefix, None)));
        }
    }
    for loss in [H2O, NH3] {
        if let Some(peptide) = truncate(parent, 0..len, Some(-loss), None) {
            out.push((peptide, kind(Placement::SameSpan, Some(loss))));
        }
    }
    out.into_iter()
        .map(|(peptide, kind)| RescueCandidate {
            peptide,
            parent: parent_idx,
            kind,
            parent_len: len,
        })
        .collect()
}

/// The span `keep` of `parent`, with its modifications moved to the new
/// positions, an optional extra N-terminal and C-terminal mass, and the mass
/// recomputed from the residues and the modifications it carries.
pub fn truncate(
    parent: &Peptide,
    keep: std::ops::Range<usize>,
    nterm_delta: Option<f32>,
    cterm_delta: Option<f32>,
) -> Option<Peptide> {
    let len = parent.sequence.len();
    if keep.start >= keep.end || keep.end > len {
        return None;
    }
    let (storage, offset) = parent.sequence.source();
    let offset = offset as usize;
    let sequence = storage.peptide(offset + keep.start..offset + keep.end)?;
    let mut applied = Vec::new();
    for modification in parent.applied_modifications() {
        let site = match modification.site {
            Site::Nterm if keep.start == 0 => Site::Nterm,
            Site::Cterm if keep.end == len => Site::Cterm,
            Site::Sequence(i) if keep.contains(&(i as usize)) => {
                Site::Sequence(i - keep.start as u32)
            }
            _ => continue,
        };
        applied.push(AppliedModification {
            site,
            modification: Arc::new(modification.modification.clone()),
            kind: modification.kind,
        });
    }
    for (site, delta) in [(Site::Nterm, nterm_delta), (Site::Cterm, cterm_delta)] {
        if let Some(delta) = delta {
            applied.push(AppliedModification {
                site,
                modification: Arc::new(ModificationDefinition::bare(delta)),
                kind: Default::default(),
            });
        }
    }
    let mut mass = H2O;
    for &residue in sequence.as_bytes() {
        mass += monoisotopic(residue);
    }
    let (mut nterm, mut cterm) = (None, None);
    for modification in &applied {
        mass += modification.modification.mass;
        match modification.site {
            Site::Nterm => nterm = Some(nterm.unwrap_or(0.0) + modification.modification.mass),
            Site::Cterm => cterm = Some(cterm.unwrap_or(0.0) + modification.modification.mass),
            Site::Sequence(_) => {}
        }
    }
    let modifications = CompactModifications::from_applied(applied).ok()?;
    let protein_sites = parent
        .protein_sites
        .iter()
        .map(|site| ProteinOccurrence {
            start: site.start.map(|start| start + keep.start as u32),
            prev_aa: None,
            next_aa: None,
            ..site.clone()
        })
        .collect::<Vec<_>>();
    Some(Peptide {
        decoy: false,
        sequence,
        modifications,
        label: parent.label.clone(),
        nterm,
        cterm,
        monoisotopic: mass,
        missed_cleavages: parent.missed_cleavages,
        semi_enzymatic: keep.start != 0 || keep.end != len || parent.semi_enzymatic,
        position: Position::Internal,
        proteins: parent.proteins.clone(),
        protein_sites: protein_sites.into(),
    })
}

/// Where a parent elutes in one file, in the file's retention-time units.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Window {
    pub start: f32,
    pub end: f32,
}

impl Window {
    pub fn contains(&self, rt: f32) -> bool {
        self.start <= rt && rt <= self.end
    }

    /// This window moved by `widths` of its own width, for the null control.
    pub fn shifted(self, widths: f32) -> Self {
        let width = self.end - self.start;
        Self {
            start: self.start + widths * width,
            end: self.end + widths * width,
        }
    }
}

/// The co-elution windows of every parent in one file: each quantified peak of
/// the parent there when quantification gave it one, otherwise its median
/// retention time plus and minus `rt_tolerance`, the tolerance
/// [`crate::in_source::annotate`] uses.
pub fn parent_windows(
    file_id: usize,
    parents: &FnvHashMap<PeptideIx, f32>,
    rt_tolerance: f32,
    apexes: Option<&PeakApexes>,
    null_shift: Option<f32>,
) -> FnvHashMap<PeptideIx, Vec<Window>> {
    parents
        .iter()
        .map(|(&parent, &median_rt)| {
            let peaks = apexes
                .map(|apexes| apexes.peaks(parent, file_id))
                .unwrap_or_default();
            let mut windows = if peaks.is_empty() {
                vec![Window {
                    start: median_rt - rt_tolerance,
                    end: median_rt + rt_tolerance,
                }]
            } else {
                peaks
                    .iter()
                    .map(|peak| Window {
                        start: peak.start_rt,
                        end: peak.end_rt,
                    })
                    .collect()
            };
            if let Some(widths) = null_shift {
                windows = windows.into_iter().map(|w| w.shifted(widths)).collect();
            }
            (parent, windows)
        })
        .collect()
}

/// Target-decoy competition q-values over `(score, is_decoy)`, higher scores
/// better: the q-value of a row is the lowest (decoys + 1) / targets of any
/// score threshold at or below its own, as the main spectrum FDR computes it.
/// Equal scores share one threshold.
pub fn competition_q_values(evidence: &[(f64, bool)]) -> Vec<f32> {
    let mut order = (0..evidence.len()).collect::<Vec<_>>();
    order.sort_by(|&a, &b| evidence[b].0.total_cmp(&evidence[a].0));
    let mut fdr = vec![1.0f32; evidence.len()];
    let (mut targets, mut decoys) = (0usize, 0usize);
    let mut start = 0;
    while start < order.len() {
        let score = evidence[order[start]].0;
        let mut end = start;
        while end < order.len() && evidence[order[end]].0 == score {
            if evidence[order[end]].1 {
                decoys += 1;
            } else {
                targets += 1;
            }
            end += 1;
        }
        let value = ((decoys + 1) as f32 / targets.max(1) as f32).min(1.0);
        fdr[start..end].fill(value);
        start = end;
    }
    let mut minimum = 1.0f32;
    for rank in (0..fdr.len()).rev() {
        minimum = minimum.min(fdr[rank]);
        fdr[rank] = minimum;
    }
    let mut out = vec![1.0f32; evidence.len()];
    for (rank, &index) in order.iter().enumerate() {
        out[index] = fdr[rank];
    }
    out
}

/// What the rescue pass did in one file, for `run-summary.json`.
#[derive(Clone, Debug, Default, PartialEq, Serialize, Deserialize, schemars::JsonSchema)]
pub struct InSourceRescueFileStats {
    pub file_id: usize,
    /// Confident target peptides of the file that were used as parents.
    pub parents: usize,
    /// Target candidates searched, after dropping those already in the
    /// database; each has one decoy.
    pub candidates: usize,
    /// MS2 spectra of the file with no confident PSM, which were searched.
    pub searched_spectra: usize,
    /// Searched spectra whose best co-eluting hit was a target / a decoy.
    pub target_hits: usize,
    pub decoy_hits: usize,
    /// Target hits at rescue q-value <= [`RESCUE_Q_VALUE`].
    pub rescued_psms: usize,
    /// Of those, how many are N-terminal pieces (suffix), C-terminal pieces
    /// (prefix), and losses.
    pub rescued_suffix: usize,
    pub rescued_prefix: usize,
    pub rescued_loss: usize,
    /// Set when the null-control harness shifted the windows.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub null_shift_widths: Option<f32>,
}

#[cfg(test)]
#[path = "../tests/unit/in_source_rescue.rs"]
mod test;
