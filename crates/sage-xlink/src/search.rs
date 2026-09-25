//! Per-spectrum crosslink search.
//!
//! 1. Find signature doublets and turn them into chain-mass pair hypotheses.
//! 2. Look each chain up as an ordinary peptide carrying a spectrum-specific
//!    mass offset: its partner plus the linker, placed on a linkable residue.
//!    The offset's neutral losses leave each linker stub, so fragments match
//!    whether the linker stayed intact or cleaved.
//! 3. Combine the best candidates of both chains, counting a peak matched by
//!    both chains once.

use crate::doublets::{
    find_doublets, pair_hypotheses, PairHypothesis, PairingSettings, PrecursorMass,
};
use crate::linker::{CleavableLinker, CrosslinkSettings};
use sage_core::database::{IndexedDatabase, MassOffset, PeptideIx};
use sage_core::ion_series::Kind;
use sage_core::mass::{Tolerance, NEUTRON, PROTON};
use sage_core::modification::{
    ModificationDefinition, ModificationSpecificity, NeutralLossMode, SiteMode,
};
use sage_core::peptide::{Peptide, Site};
use sage_core::scoring::{Fragments, OffsetHypothesis, OffsetMatch, Scorer};
use sage_core::spectrum::ProcessedSpectrum;
use std::sync::Arc;

/// Decoy status of the two chains.
#[derive(Copy, Clone, Debug, PartialEq, Eq, Hash)]
pub enum Class {
    /// Both chains target.
    TT,
    /// One chain decoy.
    TD,
    /// Both chains decoy.
    DD,
}

impl Class {
    pub fn of(alpha_decoy: bool, beta_decoy: bool) -> Self {
        match (alpha_decoy, beta_decoy) {
            (false, false) => Self::TT,
            (true, true) => Self::DD,
            _ => Self::TD,
        }
    }

    pub fn as_str(self) -> &'static str {
        match self {
            Self::TT => "TT",
            Self::TD => "TD",
            Self::DD => "DD",
        }
    }
}

/// One chain of a crosslink-spectrum match.
#[derive(Clone, Debug, PartialEq)]
pub struct ChainMatch {
    pub peptide: PeptideIx,
    /// Linked site on the peptide.
    pub site: Site,
    pub hyperscore: f64,
    pub matched_peaks: u16,
    /// Whether the chain's signature doublet was observed.
    pub doublet: bool,
    pub length: u16,
}

/// The best crosslink-spectrum match of one spectrum.
#[derive(Clone, Debug, PartialEq)]
pub struct Csm {
    pub file_id: usize,
    pub spectrum_id: String,
    pub rt: f32,
    pub charge: u8,
    pub expmass: f32,
    pub calcmass: f32,
    pub isotope_error: i8,
    /// Precursor error after the isotope correction.
    pub precursor_ppm: f32,
    /// The better-scoring chain.
    pub alpha: ChainMatch,
    pub beta: ChainMatch,
    /// Hyperscore over the union of both chains' matched peaks.
    pub hyperscore: f64,
    /// Combined hyperscore minus the next distinct candidate's.
    pub delta_next: f64,
    pub matched_peaks: u16,
    pub matched_intensity_pct: f32,
    /// Signature doublets found in the spectrum.
    pub doublets: u16,
    pub class: Class,
    /// Both chains share a protein (decoy tags ignored).
    pub intra: bool,
    pub discriminant_score: f64,
    pub csm_q: f32,
    pub residue_pair_q: f32,
}

/// A configured crosslink search.
#[derive(Clone, Debug)]
pub struct CrosslinkSearch {
    pub settings: CrosslinkSettings,
    pub linker: CleavableLinker,
    specificities: Vec<ModificationSpecificity>,
    definition_name: Arc<str>,
}

struct Candidate {
    chains: [(usize, usize); 2],
    pair: usize,
    hyperscore: f64,
    matched_peaks: u16,
    matched_intensity: f32,
}

impl CrosslinkSearch {
    pub fn new(settings: CrosslinkSettings) -> Result<Self, String> {
        settings.validate()?;
        let linker = settings.linker.resolve()?;
        let mut specificities = Vec::new();
        for residue in settings.residues.bytes() {
            // A linked residue blocks cleavage after it, so it cannot be the
            // peptide C-terminus unless that is the protein C-terminus.
            specificities.push(ModificationSpecificity::Internal(residue));
            specificities.push(ModificationSpecificity::PeptideNTerm(residue));
            specificities.push(ModificationSpecificity::ProteinCTerm(residue));
        }
        if settings.protein_n_term {
            specificities.push(ModificationSpecificity::ProteinN(None));
        }
        let definition_name = Arc::from(format!("{}_crosslinked", linker.name));
        Ok(Self {
            settings,
            linker,
            specificities,
            definition_name,
        })
    }

    /// The scorer used for chain lookups: the search scorer with the
    /// crosslink precursor tolerance and per-chain peak minimum.
    pub fn chain_scorer<'db>(&self, scorer: &Scorer<'db>) -> Scorer<'db> {
        Scorer {
            precursor_tol: self.settings.precursor_tol.unwrap_or(scorer.precursor_tol),
            min_matched_peaks: self.settings.min_chain_matched_peaks,
            ..scorer.clone()
        }
    }

    /// Offset carried by a chain whose partner has `partner_mass`.
    fn chain_offset(&self, partner_mass: f32) -> MassOffset {
        let mass = partner_mass + self.linker.crosslink_mass;
        let losses: Vec<f32> = self
            .linker
            .fragment_stubs()
            .into_iter()
            .map(|stub| mass - stub)
            .collect();
        MassOffset {
            definition: Arc::new(ModificationDefinition {
                mass,
                name: Some(self.definition_name.clone()),
                neutral_losses: Arc::from(losses),
                neutral_loss_mode: NeutralLossMode::Optional,
                channel_offsets: Arc::default(),
            }),
            specificities: self.specificities.clone(),
            site_mode: SiteMode::Exhaustive,
        }
    }

    fn chain_hypothesis(&self, chain: f32, partner: f32, charge: u8) -> OffsetHypothesis {
        let offset = self.chain_offset(partner);
        // Preliminary matching uses only the cleaved stub shifts. Adding the
        // intact-linker shift cost ~13% CPU and did not change true FDR or
        // correct identifications on the Beveridge library.
        let preliminary_shifts = self.linker.fragment_stubs();
        OffsetHypothesis {
            precursor_mass: chain + offset.mass(),
            precursor_charge: charge,
            offset,
            preliminary_shifts,
        }
    }

    fn precursors(&self, query: &ProcessedSpectrum) -> Vec<PrecursorMass> {
        let Some(precursor) = query.precursors.first() else {
            return Vec::new();
        };
        let half_width = match precursor.isolation_window {
            Some(Tolerance::Da(lo, hi)) if hi > lo => (hi - lo) / 2.0,
            _ => self.settings.isolation_half_width,
        };
        let charges = match precursor.charge {
            Some(charge) => charge..=charge,
            None => self.settings.missing_charges.0..=self.settings.missing_charges.1,
        };
        charges
            .map(|charge| PrecursorMass {
                mass: (precursor.mz - PROTON) * charge as f32,
                charge,
                half_width,
            })
            .collect()
    }

    /// Search one MS2 spectrum. `scorer` is the ordinary search scorer;
    /// search-time recalibration configured on it is applied here.
    pub fn search(&self, scorer: &Scorer, spectrum: &ProcessedSpectrum) -> Option<Csm> {
        let scorer = self.chain_scorer(scorer);
        let query = scorer.recalibrated(spectrum);
        let query = query.as_ref();
        let precursors = self.precursors(query);
        let max_precursor_charge = precursors.iter().map(|p| p.charge).max()?;
        let fragment_charge =
            fragment_charge_limit(scorer.max_fragment_charge, max_precursor_charge);

        let doublets = find_doublets(query, &self.linker, fragment_charge, scorer.fragment_tol);
        if doublets.is_empty() {
            return None;
        }
        let mut pairs = pair_hypotheses(
            &doublets,
            &precursors,
            &self.linker,
            PairingSettings {
                precursor_tol: scorer.precursor_tol,
                isotope_errors: self.settings.isotope_errors,
                min_chain_mass: self.settings.min_chain_mass,
            },
        );
        pairs.truncate(self.settings.max_pairs);
        if pairs.is_empty() {
            return None;
        }

        let hypotheses: Vec<OffsetHypothesis> = pairs
            .iter()
            .flat_map(|pair| {
                [
                    self.chain_hypothesis(pair.alpha_mass, pair.beta_mass, pair.charge),
                    self.chain_hypothesis(pair.beta_mass, pair.alpha_mass, pair.charge),
                ]
            })
            .collect();
        let matches = scorer.score_offset_hypotheses(
            query,
            &hypotheses,
            self.settings.preliminary_candidates,
            self.settings.chain_candidates,
        );

        let mut seen = vec![false; query.masses.len()];
        let mut candidates = Vec::new();
        for (pair_index, _) in pairs.iter().enumerate() {
            let (first, second) = (2 * pair_index, 2 * pair_index + 1);
            for (i, a) in matches[first].iter().enumerate() {
                for (j, b) in matches[second].iter().enumerate() {
                    let (matched_peaks, summed_b, summed_y, intensity) =
                        union_peaks(&mut seen, &a.fragments, &b.fragments);
                    candidates.push(Candidate {
                        chains: [(first, i), (second, j)],
                        pair: pair_index,
                        hyperscore: scorer.score_type.score(
                            matched_peaks.0,
                            matched_peaks.1,
                            summed_b,
                            summed_y,
                        ),
                        matched_peaks: matched_peaks.0 + matched_peaks.1,
                        matched_intensity: intensity,
                    });
                }
            }
        }
        let key = |candidate: &Candidate| {
            let chain = |(h, i): (usize, usize)| {
                let m: &OffsetMatch = &matches[h][i];
                (m.peptide_idx, m.site)
            };
            let (x, y) = (chain(candidate.chains[0]), chain(candidate.chains[1]));
            if x <= y {
                (x, y)
            } else {
                (y, x)
            }
        };
        candidates.sort_by(|a, b| {
            b.hyperscore
                .total_cmp(&a.hyperscore)
                .then_with(|| key(a).cmp(&key(b)))
        });
        let best = candidates.first()?;
        let best_key = key(best);
        let next = candidates
            .iter()
            .find(|c| key(c) != best_key)
            .map_or(0.0, |c| c.hyperscore);

        let pair: &PairHypothesis = &pairs[best.pair];
        let (a, b) = (
            &matches[best.chains[0].0][best.chains[0].1],
            &matches[best.chains[1].0][best.chains[1].1],
        );
        let db = scorer.db;
        let chain = |m: &OffsetMatch, doublet: bool| ChainMatch {
            peptide: m.peptide_idx,
            site: m.site,
            hyperscore: m.hyperscore,
            matched_peaks: m.matched_b + m.matched_y,
            doublet,
            length: db[m.peptide_idx].sequence.len() as u16,
        };
        let mut alpha = chain(a, pair.alpha_observed);
        let mut beta = chain(b, pair.beta_observed);
        if beta.hyperscore > alpha.hyperscore {
            std::mem::swap(&mut alpha, &mut beta);
        }

        let calcmass = db[alpha.peptide].monoisotopic
            + db[beta.peptide].monoisotopic
            + self.linker.crosslink_mass;
        let expmass = pair.observed_mass;
        let isotope = ((expmass - calcmass) / NEUTRON).round();
        let precursor_ppm = (expmass - isotope * NEUTRON - calcmass) / calcmass * 1e6;
        let class = Class::of(db[alpha.peptide].decoy, db[beta.peptide].decoy);
        let intra = shares_protein(db, &db[alpha.peptide], &db[beta.peptide]);
        Some(Csm {
            file_id: spectrum.file_id,
            spectrum_id: spectrum.id.clone(),
            rt: spectrum.scan_start_time,
            charge: pair.charge,
            expmass,
            calcmass,
            isotope_error: isotope.clamp(i8::MIN as f32, i8::MAX as f32) as i8,
            precursor_ppm,
            alpha,
            beta,
            hyperscore: best.hyperscore,
            delta_next: best.hyperscore - next,
            matched_peaks: best.matched_peaks,
            matched_intensity_pct: 100.0 * best.matched_intensity
                / query.total_ion_current.max(f32::MIN_POSITIVE),
            doublets: doublets.len().min(u16::MAX as usize) as u16,
            class,
            intra,
            discriminant_score: best.hyperscore,
            csm_q: 1.0,
            residue_pair_q: 1.0,
        })
    }
}

fn fragment_charge_limit(configured: Option<u8>, precursor_charge: u8) -> u8 {
    precursor_charge
        .min(configured.map(|c| c + 1).unwrap_or(precursor_charge))
        .max(2)
}

/// Count matched peaks of two fragment sets, each peak once. Returns
/// `((n-terminal, c-terminal), summed n-terminal, summed c-terminal, total
/// intensity)`. `seen` must be all false on entry and is left all false.
fn union_peaks(seen: &mut [bool], a: &Fragments, b: &Fragments) -> ((u16, u16), f32, f32, f32) {
    let (mut nb, mut ny, mut sb, mut sy) = (0u16, 0u16, 0.0f32, 0.0f32);
    for fragments in [a, b] {
        for (index, &peak) in fragments.peak_indices.iter().enumerate() {
            let peak = peak as usize;
            if seen[peak] {
                continue;
            }
            seen[peak] = true;
            let intensity = fragments.intensities[index];
            match fragments.kinds[index] {
                Kind::A | Kind::B | Kind::C => {
                    nb += 1;
                    sb += intensity;
                }
                _ => {
                    ny += 1;
                    sy += intensity;
                }
            }
        }
    }
    for fragments in [a, b] {
        for &peak in &fragments.peak_indices {
            seen[peak as usize] = false;
        }
    }
    ((nb, ny), sb, sy, sb + sy)
}

/// Strip the decoy tag from a protein accession.
pub fn target_accession<'a>(db: &IndexedDatabase, accession: &'a str) -> &'a str {
    if db.decoy_tag.is_empty() {
        accession
    } else {
        accession
            .strip_prefix(db.decoy_tag.as_str())
            .unwrap_or(accession)
    }
}

fn shares_protein(db: &IndexedDatabase, a: &Peptide, b: &Peptide) -> bool {
    a.proteins.iter().any(|x| {
        let x = target_accession(db, x);
        b.proteins.iter().any(|y| target_accession(db, y) == x)
    })
}

/// 1-based protein position of the linked residue, from the first protein
/// occurrence, and that protein's accession.
pub fn protein_position(peptide: &Peptide, site: Site) -> Option<(Arc<str>, u32)> {
    let occurrence = peptide.protein_sites.first()?;
    let offset = match site {
        Site::Nterm => 0,
        Site::Cterm => peptide.sequence.len().saturating_sub(1) as u32,
        Site::Sequence(index) => index as u32,
    };
    Some((occurrence.protein.clone(), occurrence.start? + offset + 1))
}

/// 1-based position of the linked residue in the peptide.
pub fn peptide_position(peptide: &Peptide, site: Site) -> u32 {
    match site {
        Site::Nterm => 1,
        Site::Cterm => peptide.sequence.len() as u32,
        Site::Sequence(index) => index as u32 + 1,
    }
}
