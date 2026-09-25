//! The glyco pass: one open-window search of oxonium-gated spectra against
//! sequon peptides, then glycan explanation of the precursor delta.

use std::sync::Arc;

use rayon::prelude::*;
use sage_core::database::{Builder, IndexedDatabase, Parameters, PeakExclusion};
use sage_core::fasta::Fasta;
use sage_core::ion_series::Kind;
use sage_core::mass::{Tolerance, NEUTRON, PROTON};
use sage_core::peptide::{Peptide, Site};
use sage_core::scoring::{Feature, ScoreType, Scorer};
use sage_core::spectrum::ProcessedSpectrum;

use crate::composition::{oxonium_evidence, GlycanLibrary, Monosaccharide};
use crate::config::{GlycoConfig, HEXNAC};
use crate::evidence::{ion_set, ClassCounts, IonSet, SpectrumEvidence};

/// Glycan residue masses added to the bare peptide for the Y ions put in the
/// fragment index: Y0, Y1, Y1+Fuc, and the chitobiose core with one to three
/// Hex. The HexNAc hypothesis also looks every indexed fragment up at
/// +HexNAc, which adds Y2 and Y2+Fuc.
pub const INDEXED_Y_IONS: [f64; 6] = [
    0.0,
    HEXNAC,
    HEXNAC + FUC,
    2.0 * HEXNAC + HEX,
    2.0 * HEXNAC + 2.0 * HEX,
    2.0 * HEXNAC + 3.0 * HEX,
];
const HEX: f64 = Monosaccharide::Hex.mass();
const FUC: f64 = Monosaccharide::Fuc.mass();

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
    /// Y1 (peptide + HexNAc) was matched.
    pub y1: bool,
    /// Oxonium Hex/HexNAc intensity ratio, see
    /// [`crate::composition::OxoniumEvidence::hex_ratio`].
    pub hex_ratio: f32,
    /// Random match rates of Y-ion and oxonium lookups in this spectrum.
    pub random_y: f32,
    pub random_oxonium: f32,
    pub explanations: Vec<Explanation>,
    /// Site-spanning fragment counts, when `glyco.site_features` is on.
    pub site: SiteFragments,
}

/// Matched b/y ions of a glyco candidate, from its fragment annotation.
#[derive(Copy, Clone, Debug, Default, PartialEq)]
pub struct SiteFragments {
    /// Matched b/y ions that span the glycosylation site, in any form.
    pub matched: u16,
    /// The ones matched without the HexNAc.
    pub bare: u16,
    /// Site-spanning ion groups (cleavages) of the peptide.
    pub possible: u16,
    /// Matched b/y ions (any) at fragment charge 3 or more.
    pub high_charge: u16,
}

impl SiteFragments {
    fn new(scorer: &Scorer, query: &ProcessedSpectrum, feature: &Feature) -> Self {
        let Some(Site::Sequence(site)) = feature.mass_offset.map(|assignment| assignment.site)
        else {
            return Self::default();
        };
        let site = site as i32;
        let len = feature.peptide_len as i32;
        let fragments = scorer.annotate_candidate(query, feature);
        let mut out = SiteFragments {
            possible: (len - 1).max(0) as u16,
            ..Default::default()
        };
        for i in 0..fragments.kinds.len() {
            let ordinal = fragments.fragment_ordinals[i];
            let spans = match fragments.kinds[i] {
                Kind::A | Kind::B | Kind::C => site < ordinal,
                Kind::X | Kind::Y | Kind::Z | Kind::ZDot => ordinal >= len - site,
            };
            if spans {
                out.matched += 1;
                if fragments.neutral_losses[i] > 0.0 {
                    out.bare += 1;
                }
            }
            if fragments.charges[i] >= 3 {
                out.high_charge += 1;
            }
        }
        out
    }
}

/// Glycan residue masses of the Y ions excluded from b/y matching: Y0, the
/// 0,2X cross-ring ion, HexNAc(1) with and without fucose, and every
/// HexNAc(2..=6)Hex(0..)Fuc(0..=1)NeuAc(0..=2) core-like fragment (up to 20
/// Hex on the oligomannose trunk, 12 otherwise; sialic acid only with at
/// least three HexNAc and three Hex). Sorted.
fn excluded_y_deltas() -> Vec<f32> {
    use crate::composition::Monosaccharide;
    let (hex, fuc, neuac) = (
        Monosaccharide::Hex.mass(),
        Monosaccharide::Fuc.mass(),
        Monosaccharide::NeuAc.mass(),
    );
    let mut deltas = vec![
        0.0,
        crate::composition::HEXNAC_CROSS_RING,
        HEXNAC,
        HEXNAC + fuc,
    ];
    for hexnac in 2..=6u8 {
        let max_hex = if hexnac == 2 { 20 } else { 12 };
        for h in 0..=max_hex {
            for f in 0..=1u8 {
                let max_neuac = if hexnac >= 3 && h >= 3 { 2 } else { 0 };
                for s in 0..=max_neuac {
                    deltas.push(
                        hexnac as f64 * HEXNAC + h as f64 * hex + f as f64 * fuc + s as f64 * neuac,
                    );
                }
            }
        }
    }
    let mut deltas: Vec<f32> = deltas.into_iter().map(|d| d as f32).collect();
    deltas.sort_unstable_by(f32::total_cmp);
    deltas.dedup_by(|a, b| (*a - *b).abs() < 1e-4);
    deltas
}

/// Oxonium m/z excluded from b/y matching: [`crate::composition::OXONIUM_IONS`]
/// plus HexHexNAc (366.140) and NeuAcHexHexNAc (657.235).
fn excluded_oxonium() -> Vec<f32> {
    let mut mz: Vec<f32> = crate::composition::OXONIUM_IONS
        .iter()
        .map(|ion| ion.mz)
        .collect();
    mz.extend([366.139_5, 657.235]);
    mz
}

/// Peak exclusion for full scoring: peaks at a Y ion of the candidate's
/// bare peptide (any charge up to the precursor's) or at an oxonium ion.
pub fn glycan_peak_exclusion(tolerance: Tolerance) -> Arc<PeakExclusion> {
    let deltas = excluded_y_deltas();
    let oxonium = excluded_oxonium();
    let slack = 2.0 * NEUTRON + AMMONIA as f32;
    Arc::new(
        move |peptide: &Peptide, charge: u8, query: &ProcessedSpectrum, excluded: &mut [bool]| {
            let bare = peptide.monoisotopic - HEXNAC as f32;
            let max_delta = query
                .precursors
                .first()
                .map(|precursor| (precursor.mz - PROTON) * charge as f32 - bare + slack)
                .unwrap_or(f32::MAX);
            for (idx, excluded) in excluded.iter_mut().enumerate() {
                let mz = query.peak_mz(idx);
                let known = query.has_known_charge(idx);
                if !known || query.charges[idx] <= 1 {
                    let (lo, hi) = tolerance.bounds(mz - PROTON);
                    if oxonium.iter().any(|&o| (lo..=hi).contains(&(o - PROTON))) {
                        *excluded = true;
                        continue;
                    }
                }
                let charges = if known {
                    query.charges[idx]..=query.charges[idx]
                } else {
                    1..=charge
                };
                for z in charges {
                    let neutral = (mz - PROTON) * z as f32;
                    let (lo, hi) = tolerance.bounds(neutral);
                    let (lo, hi) = (lo - bare, hi - bare);
                    if hi < 0.0 || lo > max_delta {
                        continue;
                    }
                    let start = deltas.partition_point(|&d| d < lo);
                    if deltas.get(start).is_some_and(|&d| d <= hi) {
                        *excluded = true;
                        break;
                    }
                }
            }
        },
    )
}

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
        let mut db = if config.index_y_ions {
            let y_ions = |peptide: &sage_core::peptide::Peptide, masses: &mut Vec<f32>| {
                let bare = peptide.monoisotopic as f64;
                masses.extend(INDEXED_Y_IONS.iter().map(|glycan| (bare + glycan) as f32));
            };
            parameters.build_from_peptides_with_extra_fragments(peptides, &y_ions)
        } else {
            parameters.build_from_peptides(peptides)
        };
        // Every glycopeptide carries the innermost HexNAc: the HexNAc
        // hypothesis counts fragments with and without it, so the bare
        // hypothesis would repeat the same work.
        db.offsets_only = db.mass_offsets.len() == 1;
        // Site-spanning fragments are matched in each configured form.
        for offset in db
            .mass_offsets
            .iter_mut()
            .filter(|offset| (offset.mass() as f64 - HEXNAC).abs() < 1e-3)
        {
            let mut definition = (*offset.definition).clone();
            definition.neutral_losses = config.site_losses().into();
            offset.definition = Arc::new(definition);
        }
        if config.exclude_glycan_peaks {
            db.peak_exclusion = Some(glycan_peak_exclusion(settings.fragment_tol));
        }

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
            max_fragment_charge: self
                .settings
                .max_fragment_charge
                .or(self.config.max_fragment_charge),
            chimera: false,
            report_psms: self.config.report_candidates,
            wide_window: false,
            annotate_matches: false,
            mass_shift_ppm: 20.0,
            score_type: self.settings.score_type,
            mass_recalibration: None,
        }
    }

    /// Search `spectra` (MS2 only are considered) in parallel. The
    /// candidates of one spectrum are contiguous and in peptide rank order.
    pub fn search(&self, spectra: &[ProcessedSpectrum]) -> (Vec<GlycoCandidate>, SearchCounts) {
        let scorer = self.scorer();
        let results: Vec<(bool, Vec<GlycoCandidate>)> = spectra
            .par_iter()
            .filter(|query| query.level == 2 && !query.precursors.is_empty())
            .map(|query| self.search_spectrum(&scorer, query))
            .collect();
        let counts = SearchCounts {
            spectra: results.len(),
            gated: results.iter().filter(|r| r.0).count(),
            explained: results.iter().filter(|r| !r.1.is_empty()).count(),
        };
        let candidates = results.into_iter().flat_map(|r| r.1).collect();
        (candidates, counts)
    }

    /// Gate one spectrum, search it, and explain up to
    /// `explain_candidates` of its ranked peptide candidates whose precursor
    /// delta matches a glycan. Returns (gated, candidates).
    pub fn search_spectrum(
        &self,
        scorer: &Scorer,
        query: &ProcessedSpectrum,
    ) -> (bool, Vec<GlycoCandidate>) {
        let tolerance = self.settings.fragment_tol;
        let oxonium = oxonium_evidence(query, tolerance);
        if !oxonium.is_glyco(self.config.min_oxonium_ions) {
            return (false, Vec::new());
        }
        let mut candidates = Vec::new();
        for feature in scorer.score(query) {
            if candidates.len() >= self.config.explain_candidates {
                break;
            }
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
            let site = if self.config.site_features {
                SiteFragments::new(scorer, query, &feature)
            } else {
                SiteFragments::default()
            };
            candidates.push(GlycoCandidate {
                site,
                y1: evidence.y_peak(peptide_mass, HEXNAC as f32).is_some(),
                hex_ratio: oxonium.hex_ratio,
                peptide_mass,
                oxonium_ions: oxonium.count() as u8,
                oxonium_fraction: oxonium.intensity_fraction,
                anchored,
                random_y,
                random_oxonium,
                explanations,
                feature,
            });
        }
        (true, candidates)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use sage_core::spectrum::Precursor;

    #[test]
    fn glycan_peaks_are_excluded_and_peptide_peaks_kept() {
        let mut placed = Peptide::try_from(sage_core::enzyme::Digest {
            sequence: "AANGTK".into(),
            ..Default::default()
        })
        .unwrap();
        let bare = placed.monoisotopic;
        placed.monoisotopic += HEXNAC as f32;
        let hex = crate::composition::Monosaccharide::Hex.mass() as f32;
        let glycan = 2.0 * HEXNAC as f32 + 5.0 * hex;
        // Peaks of undetermined charge, stored as singly charged neutral masses.
        let masses = [
            204.086_6 - PROTON,                     // HexNAc oxonium
            300.0,                                  // unrelated
            (bare + HEXNAC as f32) / 2.0,           // Y1 at 2+
            bare + 2.0 * HEXNAC as f32 + 3.0 * hex, // HexNAc(2)Hex(3) Y ion at 1+
        ];
        let query = ProcessedSpectrum {
            level: 2,
            masses: masses.to_vec(),
            charges: vec![1; masses.len()],
            charge_is_known: vec![false; masses.len()],
            intensities: vec![1.0; masses.len()],
            precursors: vec![Precursor {
                mz: (bare + glycan) / 3.0 + PROTON,
                charge: Some(3),
                ..Default::default()
            }],
            ..Default::default()
        };
        let exclusion = glycan_peak_exclusion(Tolerance::Ppm(-10.0, 10.0));
        let mut excluded = vec![false; masses.len()];
        exclusion(&placed, 3, &query, &mut excluded);
        assert_eq!(excluded, vec![true, false, true, true]);
    }
}
