//! Optional glycopeptide pass (`"glyco": {...}` in the config).
//!
//! The configuration block is accepted by every build so config files and
//! the committed JSON schema do not depend on build features, but only
//! builds with the `glyco` Cargo feature can run it. Everything here is
//! inert when the block is absent: normal searches are unchanged.

use sage_core::database::Builder;
use serde::Serialize;

/// The resolved glyco configuration carried by [`crate::input::Search`].
#[derive(Serialize, Clone)]
pub struct GlycoSetup {
    /// The configuration block with defaults filled in.
    pub config: serde_json::Value,
    /// The user's database settings, from which the glyco index is derived.
    #[serde(skip)]
    pub builder: Builder,
}

impl GlycoSetup {
    /// Validate the block. Fails on builds without the `glyco` feature.
    pub fn new(value: serde_json::Value, builder: &Builder) -> anyhow::Result<Self> {
        imp::validate(value, builder)
    }
}

#[cfg(feature = "glyco")]
mod imp {
    use super::GlycoSetup;
    use anyhow::Context;
    use sage_core::database::Builder;
    use sage_core::fasta::Fasta;
    use sage_core::spectrum::ProcessedSpectrum;
    use sage_glyco::{GlycoCandidate, GlycoConfig, GlycoSearch, SearchCounts, SearchSettings};
    use std::sync::Mutex;
    use std::time::{Duration, Instant};

    pub(super) fn validate(
        value: serde_json::Value,
        builder: &Builder,
    ) -> anyhow::Result<GlycoSetup> {
        let config = GlycoConfig::from_value(value).map_err(anyhow::Error::msg)?;
        sage_glyco::search::database_parameters(&config, builder).map_err(anyhow::Error::msg)?;
        Ok(GlycoSetup {
            config: serde_json::to_value(&config)?,
            builder: builder.clone(),
        })
    }

    /// The glyco index and the candidates collected over all file batches.
    pub struct GlycoState {
        search: GlycoSearch,
        candidates: Mutex<(Vec<GlycoCandidate>, SearchCounts, Duration)>,
    }

    impl GlycoState {
        pub fn build(
            setup: &GlycoSetup,
            search: &crate::input::Search,
            fasta: &Fasta,
        ) -> anyhow::Result<Self> {
            let start = Instant::now();
            let config: GlycoConfig = serde_json::from_value(setup.config.clone())?;
            let mut files = Vec::new();
            for path in &config.glycan_files {
                files.push(
                    sage_cloudpath::util::read_text(path)
                        .with_context(|| format!("Failed to read glycan list `{path}`"))?,
                );
            }
            let library = config.library(&files).map_err(anyhow::Error::msg)?;
            let parameters = sage_glyco::search::database_parameters(&config, &setup.builder)
                .map_err(anyhow::Error::msg)?;
            let settings = SearchSettings {
                fragment_tol: search.fragment_tol,
                explain_ppm: SearchSettings::explain_ppm(search.precursor_tol),
                precursor_charge: search.precursor_charge,
                override_precursor_charge: search.override_precursor_charge,
                max_fragment_charge: search.max_fragment_charge,
                min_matched_peaks: search.min_matched_peaks,
                score_type: search.score_type,
            };
            let glycans = library.len();
            let search = GlycoSearch::build(config, library, parameters, fasta, settings)
                .map_err(anyhow::Error::msg)?;
            log::info!(
                "glyco: {} glycan compositions, {} sequon peptides, {} fragments in {:#?}",
                glycans,
                search.db.peptides.len(),
                search.db.fragments.len(),
                start.elapsed()
            );
            Ok(GlycoState {
                search,
                candidates: Mutex::new(Default::default()),
            })
        }

        pub fn search(&self, spectra: &[ProcessedSpectrum]) {
            let start = Instant::now();
            let (candidates, counts) = self.search.search(spectra);
            let mut state = self.candidates.lock().expect("glyco candidates lock");
            state.0.extend(candidates);
            state.1 += counts;
            state.2 += start.elapsed();
        }

        /// Run glyco FDR and serialize `glyco.sage.parquet`.
        pub fn finish(&self, filenames: &[String]) -> anyhow::Result<Vec<u8>> {
            let start = Instant::now();
            let (candidates, counts, elapsed) =
                std::mem::take(&mut *self.candidates.lock().expect("glyco candidates lock"));
            let config = &self.search.config;
            let (psms, _, summary) = sage_glyco::fdr::score(
                candidates,
                config.isotope_errors,
                config.ammonium_adducts,
                self.search.settings.explain_ppm,
                config.peptide_fdr,
                config.glycan_fdr,
            );
            log::info!(
                "glyco: {} MS2 spectra, {} oxonium-gated, {} explained; searched in {:#?}",
                counts.spectra,
                counts.gated,
                counts.explained,
                elapsed
            );
            log::info!(
                "glyco: {} explained peptide candidates, {} selected; peptide model `{}`, {} target PSMs at peptide q <= {}, {} glycan decoy winners among them, {} glycoPSMs pass both FDRs",
                summary.explained,
                summary.candidates,
                summary.model,
                summary.peptide_passing,
                config.peptide_fdr,
                summary.glycan_decoys_passing_peptide,
                summary.passing
            );
            let bytes = sage_glyco::output::serialize(
                &psms,
                &self.search.db,
                &self.search.library,
                filenames,
            )?;
            log::debug!("glyco: FDR and output in {:#?}", start.elapsed());
            Ok(bytes)
        }
    }
}

#[cfg(not(feature = "glyco"))]
mod imp {
    use super::GlycoSetup;
    use sage_core::database::Builder;
    use sage_core::fasta::Fasta;
    use sage_core::spectrum::ProcessedSpectrum;

    pub(super) fn validate(_: serde_json::Value, _: &Builder) -> anyhow::Result<GlycoSetup> {
        anyhow::bail!("`glyco` requires a Sage build with the `glyco` feature")
    }

    /// Uninhabited without the `glyco` feature.
    pub enum GlycoState {}

    impl GlycoState {
        pub fn build(_: &GlycoSetup, _: &crate::input::Search, _: &Fasta) -> anyhow::Result<Self> {
            unreachable!("glyco setups are rejected without the glyco feature")
        }

        pub fn search(&self, _: &[ProcessedSpectrum]) {
            match *self {}
        }

        pub fn finish(&self, _: &[String]) -> anyhow::Result<Vec<u8>> {
            match *self {}
        }
    }
}

pub use imp::GlycoState;

/// Output file name of the glyco pass.
pub const FILE_NAME: &str = "glyco.sage.parquet";
