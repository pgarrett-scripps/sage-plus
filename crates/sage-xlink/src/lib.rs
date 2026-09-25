//! Search for peptides joined by an MS-cleavable crosslinker (DSSO, DSBU).
//!
//! Enabled in the `sage` binary by the `crosslink` Cargo feature and, at run
//! time, by a `"crosslink"` block in the search configuration. See
//! `docs/explore/CROSSLINK_SEARCH.md`.

pub mod doublets;
pub mod fdr;
pub mod linker;
pub mod output;
pub mod search;

pub use fdr::{assign_q_values, FdrSummary};
pub use linker::{CleavableLinker, CrosslinkSettings, LinkerConfig};
pub use search::{ChainMatch, Class, CrosslinkSearch, Csm};
