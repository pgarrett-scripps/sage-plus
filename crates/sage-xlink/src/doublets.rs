//! Signature doublets and chain-mass pair hypotheses.
//!
//! When an MS-cleavable linker (DSSO, DSBU) breaks in HCD, each released chain
//! appears twice: once carrying the short stub and once carrying the long
//! stub. The two peaks share a charge and differ by a fixed mass (31.9721 Da
//! for DSSO alkene/thiol, 25.9792 Da for DSBU Bu/BuUr). A doublet gives the
//! unmodified mass of one chain directly, and the precursor mass then fixes
//! the mass of its partner:
//!
//! `precursor = alpha + beta + crosslink_mass`
//!
//! Each pair hypothesis turns an n^2 pair search into two closed-window
//! lookups against the ordinary fragment index, one per chain mass.

use crate::linker::CleavableLinker;
use sage_core::mass::{Tolerance, NEUTRON};
use sage_core::spectrum::ProcessedSpectrum;

/// Doublets considered when pairing two observed chains, strongest first.
/// Bounds the quadratic pairing step on noisy spectra.
const MAX_PAIRED_DOUBLETS: usize = 64;

/// One chain released by linker cleavage, seen as a same-charge doublet.
#[derive(Copy, Clone, Debug, PartialEq)]
pub struct ChainDoublet {
    /// Neutral mass of the chain without any linker stub.
    pub chain_mass: f32,
    pub charge: u8,
    pub light_peak: usize,
    pub heavy_peak: usize,
    /// Summed intensity of both peaks.
    pub intensity: f32,
}

/// An observed precursor under one charge assumption.
#[derive(Copy, Clone, Debug, PartialEq)]
pub struct PrecursorMass {
    /// Neutral mass from the reported monoisotopic m/z.
    pub mass: f32,
    pub charge: u8,
    /// Half width of the isolation window, m/z.
    pub half_width: f32,
}

/// A candidate pair of chain masses.
#[derive(Copy, Clone, Debug, PartialEq)]
pub struct PairHypothesis {
    /// Heavier chain mass.
    pub alpha_mass: f32,
    /// Lighter chain mass.
    pub beta_mass: f32,
    /// True when both chains were observed as doublets, false when one chain
    /// mass is inferred from the precursor.
    pub alpha_observed: bool,
    pub beta_observed: bool,
    /// Summed doublet intensity supporting the hypothesis.
    pub intensity: f32,
    pub charge: u8,
    /// Observed precursor mass the hypothesis was derived from.
    pub observed_mass: f32,
}

impl PairHypothesis {
    pub fn both_observed(&self) -> bool {
        self.alpha_observed && self.beta_observed
    }

    /// `alpha + beta + crosslink`.
    pub fn pair_mass(&self, linker: &CleavableLinker) -> f32 {
        self.alpha_mass + self.beta_mass + linker.crosslink_mass
    }
}

#[derive(Copy, Clone)]
struct ChargedPeak {
    neutral: f32,
    charge: u8,
    index: usize,
}

/// Expand peaks to neutral masses at each plausible charge. Deisotoped peaks
/// keep their assigned charge; peaks of unknown charge (stored as `m/z -
/// proton` with charge 1) are tried at every charge below `max_charge`, as the
/// fragment matcher does.
fn charged_peaks(spectrum: &ProcessedSpectrum, max_charge: u8) -> Vec<ChargedPeak> {
    let mut peaks = Vec::new();
    for (index, &mass) in spectrum.masses.iter().enumerate() {
        if spectrum.has_known_charge(index) {
            peaks.push(ChargedPeak {
                neutral: mass,
                charge: spectrum.charges[index],
                index,
            });
        } else {
            for charge in 1..max_charge.max(2) {
                peaks.push(ChargedPeak {
                    neutral: mass * charge as f32,
                    charge,
                    index,
                });
            }
        }
    }
    peaks.sort_unstable_by(|a, b| a.neutral.total_cmp(&b.neutral));
    peaks
}

/// Find every same-charge peak pair separated by the linker's doublet
/// spacing. `max_charge` is exclusive, matching the fragment matcher.
pub fn find_doublets(
    spectrum: &ProcessedSpectrum,
    linker: &CleavableLinker,
    max_charge: u8,
    fragment_tol: Tolerance,
) -> Vec<ChainDoublet> {
    let peaks = charged_peaks(spectrum, max_charge);
    let spacing = linker.doublet_spacing();
    let mut doublets = Vec::new();
    for light in &peaks {
        let expected = light.neutral + spacing;
        // Tolerances are specified per m/z; scale Da windows by charge.
        let tolerance = match fragment_tol {
            Tolerance::Da(lo, hi) => {
                Tolerance::Da(lo * light.charge as f32, hi * light.charge as f32)
            }
            other => other,
        };
        let (lo, hi) = tolerance.bounds(expected);
        let start = peaks.partition_point(|peak| peak.neutral < lo);
        for heavy in peaks[start..].iter().take_while(|peak| peak.neutral <= hi) {
            if heavy.charge != light.charge || heavy.index == light.index {
                continue;
            }
            doublets.push(ChainDoublet {
                chain_mass: light.neutral - linker.light_stub,
                charge: light.charge,
                light_peak: light.index,
                heavy_peak: heavy.index,
                intensity: spectrum.intensities[light.index] + spectrum.intensities[heavy.index],
            });
        }
    }
    doublets
}

/// Settings for [`pair_hypotheses`].
#[derive(Copy, Clone, Debug)]
pub struct PairingSettings {
    pub precursor_tol: Tolerance,
    /// Isotope errors tried when a chain is inferred from the precursor.
    pub isotope_errors: (i8, i8),
    pub min_chain_mass: f32,
}

/// Combine doublets with the precursor into chain-mass pairs.
///
/// * Two observed doublets pair when their summed mass plus the linker lies
///   inside the precursor window: within the isolation half width of the
///   reported mass, extended by the isotope range. This tolerates a
///   misassigned monoisotopic peak; the remaining precursor error is left to
///   scoring. With a zero half width the sum must instead match the
///   precursor within tolerance at one of the isotope errors.
/// * Each doublet also proposes `precursor - isotope * neutron - crosslink -
///   chain` as its partner, for every isotope error.
///
/// Hypotheses with matching masses are merged, keeping the strongest
/// evidence. Pairs with both chains observed come first, then by intensity.
pub fn pair_hypotheses(
    doublets: &[ChainDoublet],
    precursors: &[PrecursorMass],
    linker: &CleavableLinker,
    settings: PairingSettings,
) -> Vec<PairHypothesis> {
    let tol = settings.precursor_tol;
    let mut strongest: Vec<&ChainDoublet> = doublets
        .iter()
        .filter(|d| d.chain_mass >= settings.min_chain_mass)
        .collect();
    strongest.sort_by(|a, b| b.intensity.total_cmp(&a.intensity));
    strongest.truncate(MAX_PAIRED_DOUBLETS);

    let (iso_lo, iso_hi) = settings.isotope_errors;
    let mut pairs: Vec<PairHypothesis> = Vec::new();
    let mut push = |candidate: PairHypothesis| match pairs.iter_mut().find(|existing| {
        existing.charge == candidate.charge
            && tol.contains(existing.alpha_mass, candidate.alpha_mass)
            && tol.contains(existing.beta_mass, candidate.beta_mass)
    }) {
        Some(existing) => {
            existing.alpha_observed |= candidate.alpha_observed;
            existing.beta_observed |= candidate.beta_observed;
            existing.intensity = existing.intensity.max(candidate.intensity);
        }
        None => pairs.push(candidate),
    };

    for precursor in precursors {
        let window = precursor.half_width * precursor.charge as f32;
        let lo = iso_lo.min(0) as f32 * NEUTRON - window;
        let hi = iso_hi.max(0) as f32 * NEUTRON + window;
        for (i, first) in strongest.iter().enumerate() {
            for second in &strongest[i + 1..] {
                let sum = first.chain_mass + second.chain_mass + linker.crosslink_mass;
                let delta = precursor.mass - sum;
                if delta < lo || delta > hi {
                    continue;
                }
                // With no isolation slack the pair must match the precursor
                // at one of the isotope errors.
                if precursor.half_width == 0.0
                    && !(iso_lo..=iso_hi)
                        .any(|k| tol.contains(precursor.mass - k as f32 * NEUTRON, sum))
                {
                    continue;
                }
                let (alpha, beta) = ordered(first.chain_mass, second.chain_mass);
                push(PairHypothesis {
                    alpha_mass: alpha,
                    beta_mass: beta,
                    alpha_observed: true,
                    beta_observed: true,
                    intensity: first.intensity + second.intensity,
                    charge: precursor.charge,
                    observed_mass: precursor.mass,
                });
            }
        }
        for doublet in &strongest {
            for isotope in iso_lo..=iso_hi {
                let target = precursor.mass - isotope as f32 * NEUTRON;
                let partner = target - linker.crosslink_mass - doublet.chain_mass;
                if partner < settings.min_chain_mass {
                    continue;
                }
                let partner_seen = strongest
                    .iter()
                    .any(|other| tol.contains(partner, other.chain_mass));
                let (alpha, beta) = ordered(doublet.chain_mass, partner);
                let doublet_is_alpha = alpha == doublet.chain_mass;
                push(PairHypothesis {
                    alpha_mass: alpha,
                    beta_mass: beta,
                    alpha_observed: doublet_is_alpha || partner_seen,
                    beta_observed: !doublet_is_alpha || partner_seen,
                    intensity: doublet.intensity,
                    charge: precursor.charge,
                    observed_mass: precursor.mass,
                });
            }
        }
    }

    pairs.sort_by(|a, b| {
        b.both_observed()
            .cmp(&a.both_observed())
            .then_with(|| b.intensity.total_cmp(&a.intensity))
            .then_with(|| {
                let da = (a.observed_mass - a.pair_mass(linker)).abs();
                let db = (b.observed_mass - b.pair_mass(linker)).abs();
                da.total_cmp(&db)
            })
    });
    pairs
}

fn ordered(a: f32, b: f32) -> (f32, f32) {
    if a >= b {
        (a, b)
    } else {
        (b, a)
    }
}

#[cfg(test)]
#[path = "../tests/unit/doublets.rs"]
mod test;
