//! False discovery proportion estimated from an entrapment database.
//!
//! An entrapment experiment searches the sample's proteins alongside proteins
//! that cannot be in the sample. Every discovery of an entrapment protein is
//! false, so the entrapment discoveries measure how well a search controls
//! its false discovery rate, which is what a reported q-value claims.
//!
//! Wen, Noble and colleagues (*Assessment of false discovery rate control in
//! tandem mass spectrometry analysis using entrapment*, Nature Methods 2025,
//! doi 10.1038/s41592-025-02719-x) compare four ways of turning those counts
//! into an estimate, and only two of them bound the false discovery
//! proportion from above:
//!
//! - **combined** (their equation 1), implemented here:
//!   `(N_E * (1 + 1/r)) / (N_T + N_E)`, an upper bound on the proportion of
//!   false discoveries among all reported discoveries;
//! - **paired** (equation 4) is tighter but needs each entrapment sequence to
//!   be paired with an original target sequence, which a foreign-proteome
//!   entrapment database does not give;
//! - **lower bound** (equation 2), the combined formula without the `1/r`
//!   term, can only show that a search fails to control the FDR, never that
//!   it controls it;
//! - **sample** (equation 3), the proportion among original targets alone, is
//!   called inherently flawed and shows neither.
//!
//! So the estimate reported here is the combined one, named in the output, and
//! `r` is the ratio of the entrapment to the original target database size in
//! the units of the level being estimated: peptides for the peptide level,
//! proteins for the protein level.

use serde::{Deserialize, Serialize};

/// The method this module implements, named in the output so a reader knows
/// which of the paper's estimates a number is.
pub const METHOD: &str = "combined (Wen et al. 2025, Nature Methods, eq. 1)";

/// One level's entrapment estimate: the discoveries of each class, the
/// database size ratio, and the combined estimate.
#[derive(Clone, Copy, Debug, Default, PartialEq, Serialize, Deserialize, schemars::JsonSchema)]
pub struct EntrapmentEstimate {
    /// `N_T`: reported discoveries of original (sample) targets.
    pub target_discoveries: usize,
    /// `N_E`: reported discoveries of entrapment targets.
    pub entrapment_discoveries: usize,
    /// `r`: entrapment over original target database size, in this level's
    /// own units. Absent when the database held no original targets of this
    /// kind, which leaves the estimate undefined.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub database_ratio: Option<f32>,
    /// The combined estimate, `N_E * (1 + 1/r) / (N_T + N_E)`, clamped to 1:
    /// a proportion cannot exceed one. Absent without a ratio, or with no
    /// discoveries at all.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub estimated_fdp: Option<f32>,
}

impl EntrapmentEstimate {
    /// The combined estimate of this level's false discovery proportion.
    pub fn combined(
        target_discoveries: usize,
        entrapment_discoveries: usize,
        database_ratio: Option<f32>,
    ) -> Self {
        Self {
            target_discoveries,
            entrapment_discoveries,
            database_ratio,
            estimated_fdp: combined_fdp(target_discoveries, entrapment_discoveries, database_ratio),
        }
    }
}

/// The combined estimate of the false discovery proportion among all reported
/// discoveries: `N_E * (1 + 1/r) / (N_T + N_E)`. `None` where it is not
/// defined: no ratio, a ratio of zero or worse, or nothing discovered.
pub fn combined_fdp(
    target_discoveries: usize,
    entrapment_discoveries: usize,
    database_ratio: Option<f32>,
) -> Option<f32> {
    let ratio = database_ratio.filter(|ratio| ratio.is_finite() && *ratio > 0.0)?;
    let discoveries = target_discoveries + entrapment_discoveries;
    if discoveries == 0 {
        return None;
    }
    let estimate = entrapment_discoveries as f32 * (1.0 + 1.0 / ratio) / discoveries as f32;
    Some(estimate.min(1.0))
}

/// The ratio `r` of entrapment to original target database size, `None` when
/// there are no original targets to divide by.
pub fn database_ratio(original: usize, entrapment: usize) -> Option<f32> {
    (original > 0).then(|| entrapment as f32 / original as f32)
}

#[cfg(test)]
#[path = "../tests/unit/entrapment.rs"]
mod test;
