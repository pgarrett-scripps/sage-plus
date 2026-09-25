use rayon::prelude::*;
use sage_core::spectrum::ProcessedSpectrum;
use sage_core::{scoring::Feature, tmt::TmtQuant};
use std::collections::HashMap;

pub(crate) fn prepare_local_directory(
    url: &sage_cloudpath::Url,
    overwrite: bool,
) -> anyhow::Result<()> {
    let Ok(directory) = url.to_file_path() else {
        return Ok(());
    };
    let names = [
        "run-summary.json",
        "results.json",
        "results.sage.parquet",
        "matched_fragments.sage.parquet",
        "spectral_library.sage.parquet",
        "spectral_library.mzspeclib.txt",
        "lfq.parquet",
        "results.sage.ptm-sites.parquet",
        "results.sage.protein-sites.parquet",
        "results.sage.ptm-library.parquet",
        "results.sage.ptm-library.tsv",
        "results.sage.pin",
        "results.sage.report.html",
        "crosslinks.sage.parquet",
    ];
    let existing = names
        .iter()
        .map(|name| directory.join(name))
        .filter(|path| path.symlink_metadata().is_ok())
        .collect::<Vec<_>>();
    anyhow::ensure!(overwrite || existing.is_empty(), "output directory contains Sage artifacts, choose a fresh directory or explicitly use --overwrite");
    for path in existing {
        std::fs::remove_file(path)?;
    }
    Ok(())
}

#[derive(Default)]
pub struct SageResults {
    pub ms1: Vec<ProcessedSpectrum>,
    pub features: Vec<Feature>,
    pub quant: Vec<TmtQuant>,
    /// Search-time `psm_id` -> zero-based occurrence of the PSM's spectrum ID
    /// within its file, recorded only for repeated IDs (e.g. MGF files with
    /// repeated `TITLE=` lines). Empty when every spectrum ID is unique.
    pub repeated_spectrum_psms: HashMap<usize, usize>,
    /// Best crosslink-spectrum match per spectrum, before FDR.
    #[cfg(feature = "crosslink")]
    pub crosslinks: Vec<sage_xlink::Csm>,
}

impl SageResults {
    fn fold(mut self, other: SageResults) -> Self {
        self.ms1.extend(other.ms1);
        self.features.extend(other.features);
        self.quant.extend(other.quant);
        self.repeated_spectrum_psms
            .extend(other.repeated_spectrum_psms);
        #[cfg(feature = "crosslink")]
        self.crosslinks.extend(other.crosslinks);
        self
    }
}

impl FromParallelIterator<SageResults> for SageResults {
    fn from_par_iter<I>(par_iter: I) -> Self
    where
        I: IntoParallelIterator<Item = SageResults>,
    {
        par_iter
            .into_par_iter()
            .reduce(SageResults::default, SageResults::fold)
    }
}

impl FromIterator<SageResults> for SageResults {
    fn from_iter<I>(par_iter: I) -> Self
    where
        I: IntoIterator<Item = SageResults>,
    {
        par_iter
            .into_iter()
            .fold(SageResults::default(), SageResults::fold)
    }
}

#[cfg(test)]
#[path = "../tests/unit/output.rs"]
mod tests;
