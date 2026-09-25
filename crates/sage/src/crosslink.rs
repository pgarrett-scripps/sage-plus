//! Exploratory prototype: MS-cleavable crosslinker signature doublets.
//!
//! Not wired into the search. See `docs/explore/CROSSLINK_SEARCH.md`.
//!
//! When an MS-cleavable linker (DSSO, DSBU) breaks in HCD, each released chain
//! appears twice: once carrying the short stub and once carrying the long
//! stub. The two peaks share a charge and differ by a fixed mass (31.9721 Da
//! for DSSO alkene/thiol, 25.9792 Da for DSBU Bu/BuUr). Finding a doublet
//! gives the unmodified mass of one chain directly, and the precursor mass then
//! fixes the mass of its partner:
//!
//! `precursor = alpha + beta + crosslink_mass`
//!
//! Each pair hypothesis turns an n^2 pair search into two closed-window
//! lookups against the ordinary fragment index, one per chain mass.

use crate::mass::Tolerance;
use crate::spectrum::ProcessedSpectrum;

/// An MS-cleavable crosslinker. Masses are monoisotopic and neutral.
#[derive(Copy, Clone, Debug, PartialEq)]
pub struct CleavableLinker {
    pub name: &'static str,
    /// Mass added when two residues are linked.
    pub crosslink_mass: f32,
    /// Stub left on a chain by the lighter cleavage product.
    pub light_stub: f32,
    /// Stub left on a chain by the heavier cleavage product.
    pub heavy_stub: f32,
    /// Mass added by a hydrolyzed dead-end (monolink).
    pub hydrolyzed_monolink_mass: f32,
}

/// DSSO. Stubs are alkene (C3H2O) and unsaturated thiol (C3H2OS); the
/// sulfenic acid stub (C3H4O2S, +103.9932) is the thiol plus water.
pub const DSSO: CleavableLinker = CleavableLinker {
    name: "DSSO",
    crosslink_mass: 158.003_77,
    light_stub: 54.010_56,
    heavy_stub: 85.982_63,
    hydrolyzed_monolink_mass: 176.014_33,
};

/// DSBU (BuUrBu). Stubs are Bu (C4H7NO) and BuUr (C5H5NO2).
pub const DSBU: CleavableLinker = CleavableLinker {
    name: "DSBU",
    crosslink_mass: 196.084_8,
    light_stub: 85.052_76,
    heavy_stub: 111.032_03,
    hydrolyzed_monolink_mass: 214.095_36,
};

impl CleavableLinker {
    /// Spacing between the two members of a signature doublet.
    pub fn doublet_spacing(&self) -> f32 {
        self.heavy_stub - self.light_stub
    }
}

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

/// A candidate pair of chain masses consistent with the precursor.
#[derive(Copy, Clone, Debug, PartialEq)]
pub struct PairHypothesis {
    /// Heavier chain mass.
    pub alpha_mass: f32,
    /// Lighter chain mass.
    pub beta_mass: f32,
    /// True when both chains were observed as doublets, false when one chain
    /// mass is only inferred from the precursor.
    pub both_observed: bool,
    /// Summed doublet intensity supporting the hypothesis.
    pub intensity: f32,
}

#[derive(Copy, Clone)]
struct ChargedPeak {
    neutral: f32,
    charge: u8,
    index: usize,
}

/// Expand peaks to neutral masses at each plausible charge. Deisotoped peaks
/// keep their assigned charge; peaks of unknown charge (stored as `m/z -
/// proton` with charge 1) are tried at every charge below `max_charge`, as
/// `FragmentMatchIndex` does.
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
/// spacing. `max_charge` is exclusive, matching `max_fragment_charge`.
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

/// Combine doublets with the precursor mass into chain-mass pairs. Each
/// doublet proposes itself as one chain and `precursor - crosslink - chain` as
/// the other. Hypotheses whose two masses match an existing one within
/// `precursor_tol` are merged, keeping the strongest evidence.
pub fn pair_hypotheses(
    doublets: &[ChainDoublet],
    precursor_mass: f32,
    linker: &CleavableLinker,
    precursor_tol: Tolerance,
    min_chain_mass: f32,
) -> Vec<PairHypothesis> {
    let mut pairs: Vec<PairHypothesis> = Vec::new();
    for doublet in doublets {
        let partner = precursor_mass - linker.crosslink_mass - doublet.chain_mass;
        if partner < min_chain_mass || doublet.chain_mass < min_chain_mass {
            continue;
        }
        let partner_doublet = doublets
            .iter()
            .filter(|other| precursor_tol.contains(partner, other.chain_mass))
            .max_by(|a, b| a.intensity.total_cmp(&b.intensity));
        let (alpha, beta) = if doublet.chain_mass >= partner {
            (doublet.chain_mass, partner)
        } else {
            (partner, doublet.chain_mass)
        };
        let candidate = PairHypothesis {
            alpha_mass: alpha,
            beta_mass: beta,
            both_observed: partner_doublet.is_some(),
            intensity: doublet.intensity + partner_doublet.map_or(0.0, |d| d.intensity),
        };
        match pairs.iter_mut().find(|existing| {
            precursor_tol.contains(existing.alpha_mass, candidate.alpha_mass)
                && precursor_tol.contains(existing.beta_mass, candidate.beta_mass)
        }) {
            Some(existing) => {
                existing.both_observed |= candidate.both_observed;
                existing.intensity = existing.intensity.max(candidate.intensity);
            }
            None => pairs.push(candidate),
        }
    }
    // Fully observed pairs first, then by supporting intensity.
    pairs.sort_by(|a, b| {
        b.both_observed
            .cmp(&a.both_observed)
            .then_with(|| b.intensity.total_cmp(&a.intensity))
    });
    pairs
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::mass::PROTON;

    const ALPHA: f32 = 1500.7342;
    const BETA: f32 = 1010.5127;

    /// Build a spectrum from (neutral mass, charge, known charge, intensity).
    fn spectrum(peaks: &[(f32, u8, bool, f32)]) -> ProcessedSpectrum {
        let mut stored: Vec<(f32, u8, bool, f32)> = peaks
            .iter()
            .map(|&(neutral, charge, known, intensity)| {
                if known {
                    (neutral, charge, true, intensity)
                } else {
                    // Unknown-charge peaks are stored as m/z - proton, charge 1.
                    let mz = neutral / charge as f32 + PROTON;
                    (mz - PROTON, 1, false, intensity)
                }
            })
            .collect();
        stored.sort_by(|a, b| a.0.total_cmp(&b.0));
        ProcessedSpectrum {
            level: 2,
            masses: stored.iter().map(|p| p.0).collect(),
            charges: stored.iter().map(|p| p.1).collect(),
            charge_is_known: stored.iter().map(|p| p.2).collect(),
            intensities: stored.iter().map(|p| p.3).collect(),
            total_ion_current: stored.iter().map(|p| p.3).sum(),
            ..Default::default()
        }
    }

    fn signature(
        linker: &CleavableLinker,
        chain: f32,
        charge: u8,
        known: bool,
    ) -> [(f32, u8, bool, f32); 2] {
        [
            (chain + linker.light_stub, charge, known, 100.0),
            (chain + linker.heavy_stub, charge, known, 80.0),
        ]
    }

    const TOL: Tolerance = Tolerance::Ppm(-10.0, 10.0);

    #[test]
    fn four_peak_dsso_signature_yields_one_observed_pair() {
        let mut peaks = Vec::new();
        peaks.extend(signature(&DSSO, ALPHA, 2, true));
        peaks.extend(signature(&DSSO, BETA, 3, false));
        // Ordinary backbone fragments that must not form doublets.
        peaks.extend([
            (500.25, 1, true, 50.0),
            (620.31, 1, true, 40.0),
            (731.4, 2, true, 30.0),
        ]);
        let spectrum = spectrum(&peaks);
        let precursor = ALPHA + BETA + DSSO.crosslink_mass;

        let doublets = find_doublets(&spectrum, &DSSO, 5, TOL);
        assert_eq!(doublets.len(), 2, "{doublets:?}");

        let pairs = pair_hypotheses(&doublets, precursor, &DSSO, TOL, 500.0);
        assert_eq!(pairs.len(), 1, "{pairs:?}");
        let pair = pairs[0];
        assert!(pair.both_observed);
        assert!(TOL.contains(ALPHA, pair.alpha_mass));
        assert!(TOL.contains(BETA, pair.beta_mass));
    }

    #[test]
    fn single_doublet_infers_partner_from_precursor() {
        let spectrum = spectrum(&signature(&DSSO, ALPHA, 2, true));
        let precursor = ALPHA + BETA + DSSO.crosslink_mass;
        let doublets = find_doublets(&spectrum, &DSSO, 4, TOL);
        let pairs = pair_hypotheses(&doublets, precursor, &DSSO, TOL, 500.0);
        assert_eq!(pairs.len(), 1);
        assert!(!pairs[0].both_observed);
        assert!(TOL.contains(BETA, pairs[0].beta_mass));
    }

    #[test]
    fn dsbu_spacing_is_distinct_from_dsso() {
        let spectrum = spectrum(&signature(&DSBU, ALPHA, 2, true));
        assert_eq!(find_doublets(&spectrum, &DSBU, 4, TOL).len(), 1);
        assert!(find_doublets(&spectrum, &DSSO, 4, TOL).is_empty());
    }

    #[test]
    fn doublet_members_must_share_charge() {
        let [light, heavy] = signature(&DSSO, ALPHA, 2, true);
        // Same neutral spacing, but the heavy member is deisotoped at charge 3.
        let spectrum = spectrum(&[light, (heavy.0, 3, true, heavy.3)]);
        assert!(find_doublets(&spectrum, &DSSO, 5, TOL).is_empty());
    }

    #[test]
    fn stub_masses_are_consistent_with_linker_mass() {
        // DSSO: alkene + sulfenic acid = crosslink; thiol + water = sulfenic.
        let sulfenic = DSSO.heavy_stub + crate::mass::H2O;
        assert!((DSSO.light_stub + sulfenic - DSSO.crosslink_mass).abs() < 1e-3);
        // DSBU: Bu + BuUr = crosslink.
        assert!((DSBU.light_stub + DSBU.heavy_stub - DSBU.crosslink_mass).abs() < 1e-3);
        for linker in [DSSO, DSBU] {
            let water = linker.hydrolyzed_monolink_mass - linker.crosslink_mass;
            assert!((water - crate::mass::H2O).abs() < 1e-3);
        }
    }
}
