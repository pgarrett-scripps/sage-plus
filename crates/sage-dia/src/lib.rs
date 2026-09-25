//! DIA search for Sage Plus.
//!
//! [`pipeline::pseudo_spectra`] backs `dia.mode = "pseudo"` in sage-cli: it
//! turns a DIA run into MS1-anchored pseudo-MS2 spectra. The plan is in `docs/explore/DIA_SEARCH.md`. A
//! spectrum-centric wide-window search proposes candidate peptides, koth
//! detects chromatographic hills for MS1 precursors and for MS2 fragments in
//! each isolation window, and [`coelution`] measures how well a candidate's
//! fragment hills co-elute with each other and with its precursor hill.

pub use koth_core;

pub mod coelution;
pub mod hills;
pub mod pipeline;
pub mod pseudo;
pub mod tims;

pub use coelution::{CoelutionFeatures, CoelutionSettings};
pub use hills::{Channel, CompactHill};
pub use pipeline::{prepare, pseudo_spectra, pseudo_spectra_tdf, DiaMode, DiaSettings};
pub use pseudo::{PrecursorTrace, PseudoSettings, PseudoSpectrum};
