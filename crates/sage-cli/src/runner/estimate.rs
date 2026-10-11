//! `sage --estimate`: a rough database memory preview that runs no search.
//!
//! Estimates are counted from the FASTA digest and modification rules. They
//! can be far from real usage, so they are only reported here and never stop
//! a run; `max_memory_gb` is enforced on measured memory during the search.

use super::{load_custom_cleavages, load_fasta, load_ptm_library};
use crate::input::Search;
use sage_core::database::DatabaseMemoryEstimate;
use std::fmt;

const GIB: f64 = 1024.0 * 1024.0 * 1024.0;

#[derive(Debug, Clone)]
pub struct MemoryEstimateReport {
    pub proteins: usize,
    pub database: DatabaseMemoryEstimate,
    pub prefilter: bool,
    pub max_memory_gb: Option<f64>,
    /// Configured inputs whose peptides the estimate does not count.
    pub not_counted: Vec<String>,
}

/// Estimate database memory for `parameters` without reading spectra.
pub fn estimate_memory(parameters: &Search) -> anyhow::Result<MemoryEstimateReport> {
    let mut database_parameters = parameters.database.clone();
    anyhow::ensure!(
        !database_parameters.fasta.is_empty(),
        "`--estimate` needs database.fasta"
    );
    load_ptm_library(&mut database_parameters)?;
    let fasta = load_fasta(&database_parameters)?;
    let custom_cleavages = load_custom_cleavages(&database_parameters, &fasta)?;
    let database = database_parameters
        .estimate_memory_with_custom_cleavages(&fasta, custom_cleavages.as_ref());

    let mut not_counted = Vec::new();
    if let Some(peptides) = &database_parameters.peptides {
        not_counted.push(format!("peptide list `{peptides}`"));
    }

    Ok(MemoryEstimateReport {
        proteins: fasta.targets.len(),
        database,
        prefilter: database_parameters.prefilter,
        max_memory_gb: parameters.max_memory_gb,
        not_counted,
    })
}

impl fmt::Display for MemoryEstimateReport {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let db = &self.database;
        writeln!(
            f,
            "Rough memory estimate (real usage can differ a lot; this never limits a run)"
        )?;
        writeln!(f, "  proteins:             {}", self.proteins)?;
        writeln!(
            f,
            "  unmodified peptides:  {:>14}  ~{:.2} GiB",
            db.unmodified_peptides,
            db.unmodified_peak_bytes as f64 / GIB
        )?;
        writeln!(
            f,
            "  modified peptides:    {:>14}  ~{:.2} GiB",
            db.modified_peptides,
            db.modified_peak_bytes as f64 / GIB
        )?;
        writeln!(
            f,
            "  fragments:            {:>14}  ~{:.2} GiB fragment index",
            db.fragments,
            db.fragment_peak_bytes as f64 / GIB
        )?;
        if self.prefilter {
            writeln!(
                f,
                "  prefilter:            on; the final index holds only prefilter survivors, so \
                 it is usually far smaller than the fragment figure above"
            )?;
        } else {
            writeln!(f, "  prefilter:            off")?;
        }
        match self.max_memory_gb {
            Some(gib) => writeln!(
                f,
                "  max_memory_gb:        {gib:.2} GiB (enforced on measured memory)"
            )?,
            None => writeln!(f, "  max_memory_gb:        not set")?,
        }
        writeln!(f, "  not counted:          spectra and search buffers")?;
        for input in &self.not_counted {
            writeln!(f, "                        {input}")?;
        }
        Ok(())
    }
}
