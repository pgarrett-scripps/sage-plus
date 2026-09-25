//! The glyco pass: one open-window search of oxonium-gated spectra against
//! sequon peptides, then glycan explanation of the precursor delta.

use rayon::prelude::*;
use sage_core::database::{Builder, IndexedDatabase, Parameters};
use sage_core::fasta::Fasta;
use sage_core::mass::{Tolerance, NEUTRON};
use sage_core::scoring::{Feature, ScoreType, Scorer};
use sage_core::spectrum::ProcessedSpectrum;

use crate::composition::{oxonium_evidence, GlycanLibrary};
use crate::config::{GlycoConfig, HEXNAC};
use crate::evidence::{ion_set, ClassCounts, IonSet, SpectrumEvidence};

/// NH3, the neutral mass of a noncovalent ammonium adduct.
pub const AMMONIA: f64 = 17.026_549;

/// Search settings shared with the main search.
#[derive(Copy, Clone, Debug)]
pub struct SearchSettings {
    pub fragment_tol: Tolerance,
    /// Tolerance, in ppm of the precursor mass, for explaining the glycan
    /// mass.
    pub explain_ppm: f32,
    pub precursor_charge: (u8, u8),
    pub override_precursor_charge: bool,
    pub max_fragment_charge: Option<u8>,
    pub min_matched_peaks: u16,
    pub score_type: ScoreType,
}

impl SearchSettings {
    /// Explanation tolerance from the main precursor tolerance: its widest
    /// ppm bound, or the Da bound converted at 2000 Da.
    pub fn explain_ppm(precursor_tol: Tolerance) -> f32 {
        match precursor_tol {
            Tolerance::Ppm(lo, hi) => lo.abs().max(hi.abs()),
            Tolerance::Da(lo, hi) => lo.abs().max(hi.abs()) / 2000.0 * 1e6,
            Tolerance::Pct(lo, hi) => lo.abs().max(hi.abs()) * 1e4,
        }
    }
}

/// Database parameters of the glyco index: the user's database settings
/// plus the labile HexNAc mass offset at the sequon.
pub fn database_parameters(config: &GlycoConfig, builder: &Builder) -> Result<Parameters, String> {
    let extra: Builder = serde_json::from_value(serde_json::json!({
        "variable_mods": { "HexNAc": config.hexnac_offset() }
    }))
    .map_err(|error| format!("glyco: invalid HexNAc offset: {error}"))?;
    let mut builder = builder.clone();
    if !config.variable_mods {
        builder.variable_mods = None;
    }
    let mods = builder.variable_mods.get_or_insert_with(Default::default);
    if mods.contains_key("HexNAc") {
        return Err("glyco: `database.variable_mods` must not already define `HexNAc`".into());
    }
    mods.extend(extra.variable_mods.unwrap_or_default());
    builder.bucket_size = Some(config.bucket_size);
    builder.validate_modification_keys()?;
    Ok(builder.make_parameters())
}

/// One way to explain a precursor delta: a composition, an isotope error
/// and a number of ammonium adducts, with the fragment evidence for it and
/// for its decoy twin.
#[derive(Clone, Debug)]
pub struct Explanation {
    /// Index into the [`GlycanLibrary`].
    pub composition: u32,
    pub isotope: i8,
    pub adducts: u8,
    pub error_ppm: f32,
    pub target: ClassCounts,
    pub decoy: ClassCounts,
    pub y_intensity: f32,
}

/// The best explained peptide candidate of one gated spectrum.
#[derive(Clone, Debug)]
pub struct GlycoCandidate {
    pub feature: Feature,
    /// Bare (unglycosylated) peptide mass.
    pub peptide_mass: f32,
    pub oxonium_ions: u8,
    pub oxonium_fraction: f32,
    /// Y0 or Y1 (peptide + HexNAc) was matched.
    pub anchored: bool,
    /// Random match rates of Y-ion and oxonium lookups in this spectrum.
    pub random_y: f32,
    pub random_oxonium: f32,
    pub explanations: Vec<Explanation>,
}

/// Counters for the run summary.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct SearchCounts {
    pub spectra: usize,
    pub gated: usize,
    pub explained: usize,
}

impl std::ops::AddAssign for SearchCounts {
    fn add_assign(&mut self, rhs: Self) {
        self.spectra += rhs.spectra;
        self.gated += rhs.gated;
        self.explained += rhs.explained;
    }
}

pub struct GlycoSearch {
    pub config: GlycoConfig,
    pub library: GlycanLibrary,
    ions: Vec<IonSet>,
    pub db: IndexedDatabase,
    pub settings: SearchSettings,
    precursor_tol: Tolerance,
}

impl GlycoSearch {
    /// Build the sequon peptide index from `fasta`.
    pub fn build(
        config: GlycoConfig,
        library: GlycanLibrary,
        parameters: Parameters,
        fasta: &Fasta,
        settings: SearchSettings,
    ) -> Result<Self, String> {
        let (gmin, gmax) = library
            .mass_range()
            .ok_or("glyco: the glycan library is empty")?;
        let offsets: Vec<_> = parameters
            .mass_offset_modifications()
            .into_iter()
            .filter(|offset| (offset.mass() as f64 - HEXNAC).abs() < 1e-3)
            .collect();
        // Only digests with an asparagine can carry a sequon; filtering
        // before modification expansion keeps the transient digest small.
        let mut digests = parameters.digest_unmodified(fasta);
        digests.retain(|digest| digest.reference.sequence.as_bytes().contains(&b'N'));
        let mut sites = Vec::new();
        let peptides: Vec<_> = parameters
            .modify_digests(digests)
            .into_iter()
            .filter(|peptide| {
                sites.clear();
                for offset in &offsets {
                    for specificity in &offset.specificities {
                        peptide.compatible_sites(*specificity, &mut sites);
                    }
                }
                !sites.is_empty()
            })
            .collect();
        let mut db = parameters.build_from_peptides(peptides);
        // Every glycopeptide carries the innermost HexNAc: the HexNAc
        // hypothesis counts fragments with and without it, so the bare
        // hypothesis would repeat the same work.
        db.offsets_only = db.mass_offsets.len() == 1;

        // The HexNAc hypothesis searches bare peptides at precursor - HexNAc,
        // so the window covers explained deltas of glycan mass, isotope
        // errors and adducts.
        let dmin = gmin + config.isotope_errors.0.min(0) as f64 * NEUTRON as f64;
        let dmax = gmax
            + config.isotope_errors.1.max(0) as f64 * NEUTRON as f64
            + config.ammonium_adducts as f64 * AMMONIA;
        let slack = 0.2;
        let precursor_tol = Tolerance::Da(
            (HEXNAC - dmax - slack) as f32,
            (HEXNAC - dmin + slack) as f32,
        );
        let ions = library.iter().map(|(_, c)| ion_set(c)).collect();
        Ok(GlycoSearch {
            config,
            library,
            ions,
            db,
            settings,
            precursor_tol,
        })
    }

    pub fn scorer(&self) -> Scorer<'_> {
        Scorer {
            db: &self.db,
            precursor_tol: self.precursor_tol,
            fragment_tol: self.settings.fragment_tol,
            min_matched_peaks: self.settings.min_matched_peaks,
            min_isotope_err: 0,
            max_isotope_err: 0,
            min_precursor_charge: self.settings.precursor_charge.0,
            max_precursor_charge: self.settings.precursor_charge.1,
            override_precursor_charge: self.settings.override_precursor_charge,
            max_fragment_charge: self.settings.max_fragment_charge.or(Some(3)),
            chimera: false,
            report_psms: self.config.report_candidates,
            wide_window: false,
            annotate_matches: false,
            mass_shift_ppm: 20.0,
            score_type: self.settings.score_type,
            mass_recalibration: None,
        }
    }

    /// Search `spectra` (MS2 only are considered) in parallel.
    pub fn search(&self, spectra: &[ProcessedSpectrum]) -> (Vec<GlycoCandidate>, SearchCounts) {
        let scorer = self.scorer();
        let results: Vec<(bool, Option<GlycoCandidate>)> = spectra
            .par_iter()
            .filter(|query| query.level == 2 && !query.precursors.is_empty())
            .map(|query| self.search_spectrum(&scorer, query))
            .collect();
        let counts = SearchCounts {
            spectra: results.len(),
            gated: results.iter().filter(|r| r.0).count(),
            explained: results.iter().filter(|r| r.1.is_some()).count(),
        };
        let candidates = results.into_iter().filter_map(|r| r.1).collect();
        (candidates, counts)
    }

    /// Gate one spectrum, search it, and explain its best-ranked candidate
    /// whose precursor delta matches a glycan. Returns (gated, candidate).
    pub fn search_spectrum(
        &self,
        scorer: &Scorer,
        query: &ProcessedSpectrum,
    ) -> (bool, Option<GlycoCandidate>) {
        let tolerance = self.settings.fragment_tol;
        let oxonium = oxonium_evidence(query, tolerance);
        if !oxonium.is_glyco(self.config.min_oxonium_ions) {
            return (false, None);
        }
        for feature in scorer.score(query) {
            let base = &self.db.peptides[feature.peptide_idx.0 as usize];
            let peptide_mass = base.monoisotopic;
            let delta = feature.expmass as f64 - peptide_mass as f64;
            let tol_da = feature.expmass as f64 * self.settings.explain_ppm as f64 * 1e-6;
            let mut found = Vec::new();
            for adducts in 0..=self.config.ammonium_adducts {
                let isotopes = self.config.isotope_errors.0..=self.config.isotope_errors.1;
                for (index, isotope) in
                    self.library
                        .explain(delta - adducts as f64 * AMMONIA, tol_da, isotopes)
                {
                    found.push((index, isotope, adducts));
                }
            }
            if found.is_empty() {
                continue;
            }
            let evidence = SpectrumEvidence::new(query, tolerance, feature.charge);
            let mut anchored = false;
            let explanations = found
                .into_iter()
                .map(|(index, isotope, adducts)| {
                    let mass = self.library.get(index).expect("library index").mass();
                    let (target, decoy, y_intensity, anchor) =
                        evidence.count(peptide_mass, &self.ions[index]);
                    anchored |= anchor;
                    let error =
                        delta - isotope as f64 * NEUTRON as f64 - adducts as f64 * AMMONIA - mass;
                    Explanation {
                        composition: index as u32,
                        isotope,
                        adducts,
                        error_ppm: (error / feature.expmass as f64 * 1e6) as f32,
                        target,
                        decoy,
                        y_intensity,
                    }
                })
                .collect::<Vec<_>>();
            let (random_y, random_oxonium) = evidence.random_rates(peptide_mass, delta as f32);
            return (
                true,
                Some(GlycoCandidate {
                    peptide_mass,
                    oxonium_ions: oxonium.count() as u8,
                    oxonium_fraction: oxonium.intensity_fraction,
                    anchored,
                    random_y,
                    random_oxonium,
                    explanations,
                    feature,
                }),
            );
        }
        (true, None)
    }
}
