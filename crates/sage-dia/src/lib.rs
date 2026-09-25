//! Hill-based DIA rescoring for Sage Plus.
//!
//! Exploration crate; the plan is in `docs/explore/DIA_SEARCH.md`. A
//! spectrum-centric wide-window search proposes candidate peptides, koth
//! detects chromatographic hills for MS1 precursors and for MS2 fragments in
//! each isolation window, and [`coelution`] measures how well a candidate's
//! fragment hills co-elute with each other and with its precursor hill.

pub mod coelution;
pub mod hills;
pub mod pseudo;

pub use coelution::{CoelutionFeatures, CoelutionSettings};
pub use hills::{Channel, CompactHill};
pub use pseudo::{PrecursorTrace, PseudoSettings, PseudoSpectrum};
