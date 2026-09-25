//! Experimental intact N-glycopeptide search for Sage Plus.
//!
//! A second, glyco-specific pass that runs next to the normal search:
//!
//! 1. spectra are gated on oxonium ions ([`composition::oxonium_evidence`]),
//! 2. gated spectra are searched once, with an open precursor window, against
//!    an index of sequon peptides carrying a labile HexNAc mass offset
//!    ([`search::GlycoSearch`]),
//! 3. the precursor delta of each candidate is explained by glycan
//!    compositions, isotope errors and ammonium adducts, with Y-ion and
//!    oxonium evidence for each explanation and its decoy twin
//!    ([`evidence`]),
//! 4. peptide and glycan FDR are estimated separately ([`fdr`]) and written
//!    to `glyco.sage.parquet` ([`output`]).
//!
//! Design notes and measurements: `docs/explore/GLYCO_SEARCH.md`.

pub mod composition;
pub mod config;
pub mod evidence;
pub mod fdr;
pub mod output;
pub mod search;

pub use config::GlycoConfig;
pub use fdr::{FdrSummary, GlycoPsm};
pub use search::{GlycoCandidate, GlycoSearch, SearchCounts, SearchSettings};
