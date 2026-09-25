//! Fragment evidence for a glycan explanation, and its decoy twin.
//!
//! Each explanation (composition, isotope error, adducts) is checked against
//! an extended Y-ion set: every structurally plausible sub-composition that is
//! small (the core region, where most Y signal is), close to the full
//! composition (up to three residue losses from the precursor), or on the
//! chitobiose-mannose trunk (HexNAc(2)Hex(n), with any core fucose) when the
//! composition itself is oligomannose-type (at most two HexNAc, no sialic
//! acid), which carries the long Hex ladders of high-mannose and yeast
//! glycans. Ions are
//! grouped into [`CLASSES`], and only per-class counts (generated, matched) are
//! kept, which is all the likelihood score in [`crate::fdr`] needs. Spectra
//! can therefore be released right after the search and the score model
//! refitted globally.
//!
//! The decoy twin of a composition has the same precursor mass and the same
//! core ions (which every candidate composition shares), but its non-core Y
//! ions and its sialic-acid oxonium ions are moved by a deterministic
//! pseudo-random 4.5-12.5 Da. A wrong composition therefore scores like its decoy,
//! and target-decoy competition between the two estimates glycan FDR.

use sage_core::mass::{Tolerance, PROTON};
use sage_core::spectrum::ProcessedSpectrum;

use crate::composition::{
    find_at_charge, find_singly_charged, GlycanComposition, Monosaccharide, HEXNAC_CROSS_RING,
};

/// Evidence classes.
pub const CLASSES: usize = 8;
/// Core Y ions shared by every N-glycan: Y0, Y0+0,2X, HexNAc(1..2), HexNAc(2)Hex(1..3).
pub const CORE: usize = 0;
/// Core Y ions carrying fucose.
pub const CORE_FUC: usize = 1;
/// Larger Y ions without sialic acid.
pub const EXTENDED: usize = 2;
/// Y ions carrying sialic acid.
pub const SIALYL: usize = 3;
/// NeuAc oxonium ions (274.092, 292.103) when the composition has NeuAc.
pub const NEUAC_PRESENT: usize = 4;
/// The same ions when the composition has none: a match is evidence against it.
pub const NEUAC_ABSENT: usize = 5;
/// NeuGc oxonium ions (290.087, 308.098) when the composition has NeuGc.
pub const NEUGC_PRESENT: usize = 6;
/// The same ions when the composition has none.
pub const NEUGC_ABSENT: usize = 7;

pub const CLASS_NAMES: [&str; CLASSES] = [
    "core",
    "core_fuc",
    "extended",
    "sialyl",
    "neuac_present",
    "neuac_absent",
    "neugc_present",
    "neugc_absent",
];

const NEUAC_OXONIUM: [f32; 2] = [274.092_1, 292.102_7];
const NEUGC_OXONIUM: [f32; 2] = [290.087, 308.097_6];

/// Largest charge at which Y ions are looked up.
const MAX_Y_CHARGE: u8 = 4;
/// Sub-compositions with at most this many residues are always generated.
const SMALL_RESIDUES: u8 = 6;
/// Sub-compositions this many residues or fewer below the full composition
/// are always generated.
const LARGE_LOSS: u8 = 3;
/// Probes used to estimate the random match rate of one spectrum.
const Y_PROBES: usize = 64;
const OXONIUM_PROBES: usize = 32;

/// One Y ion of a composition: mass added to the bare peptide, its class,
/// and the decoy's mass.
#[derive(Copy, Clone, Debug, PartialEq)]
pub struct YIonDef {
    pub delta: f32,
    pub decoy_delta: f32,
    pub class: u8,
}

/// Y ions and oxonium ions to check for one composition.
#[derive(Clone, Debug, Default)]
pub struct IonSet {
    pub y: Vec<YIonDef>,
    /// Sialic-acid oxonium m/z for the target and decoy, with their class.
    pub oxonium: Vec<(f32, f32, u8)>,
}

/// Per-class generated and matched ion counts.
#[derive(Copy, Clone, Debug, Default, PartialEq, Eq)]
pub struct ClassCounts {
    pub generated: [u16; CLASSES],
    pub matched: [u16; CLASSES],
}

impl ClassCounts {
    pub fn y_generated(&self) -> u32 {
        self.generated[..4].iter().map(|&n| n as u32).sum()
    }

    pub fn y_matched(&self) -> u32 {
        self.matched[..4].iter().map(|&n| n as u32).sum()
    }
}

fn splitmix(mut x: u64) -> u64 {
    x = x.wrapping_add(0x9e37_79b9_7f4a_7c15);
    x = (x ^ (x >> 30)).wrapping_mul(0xbf58_476d_1ce4_e5b9);
    x = (x ^ (x >> 27)).wrapping_mul(0x94d0_49bb_1331_11eb);
    x ^ (x >> 31)
}

fn key(composition: &GlycanComposition) -> u64 {
    composition
        .0
        .iter()
        .fold(0u64, |acc, &count| (acc << 8) | count as u64)
}

/// Deterministic decoy shift of 4.5-12.5 Da either way. The range avoids
/// the mass differences that relate real glycan ions to each other: isotope
/// spacing, NH3 and H2O losses, and the 16 Da O difference between Fuc and
/// Hex or NeuAc and NeuGc.
fn decoy_shift(composition: &GlycanComposition, ion: u64) -> f32 {
    let hash = splitmix(key(composition) ^ splitmix(ion.wrapping_add(0x5a17)));
    let magnitude = 4.5 + (hash >> 11) as f32 / (1u64 << 53) as f32 * 8.0;
    if hash & 1 == 0 {
        magnitude
    } else {
        -magnitude
    }
}

/// Whether `sub` could be a Y-ion glycan on the peptide: a connected
/// fragment of an N-glycan rooted at the sequon asparagine.
fn plausible(sub: &GlycanComposition) -> bool {
    let [hexnac, hex, _fuc, neuac, neugc] = sub.0;
    if sub.0.iter().all(|&count| count == 0) {
        return true;
    }
    if hexnac == 0 {
        return false;
    }
    if hex > 0 && hexnac < 2 {
        return false;
    }
    if hexnac > 2 && hex < 2 {
        return false;
    }
    if neuac + neugc > 0 && (hexnac < 3 || hex < 4) {
        return false;
    }
    true
}

fn class_of(sub: &GlycanComposition) -> u8 {
    let [hexnac, hex, fuc, neuac, neugc] = sub.0;
    if neuac + neugc > 0 {
        SIALYL as u8
    } else if hexnac <= 2 && hex <= 3 {
        if fuc > 0 {
            CORE_FUC as u8
        } else {
            CORE as u8
        }
    } else {
        EXTENDED as u8
    }
}

/// The Y-ion and oxonium set of `composition`.
pub fn ion_set(composition: &GlycanComposition) -> IonSet {
    let full = composition.0;
    let total: u8 = full.iter().sum();
    let mut y = vec![
        YIonDef {
            delta: 0.0,
            decoy_delta: 0.0,
            class: CORE as u8,
        },
        YIonDef {
            delta: HEXNAC_CROSS_RING as f32,
            decoy_delta: HEXNAC_CROSS_RING as f32,
            class: CORE as u8,
        },
    ];
    for hexnac in 1..=full[0] {
        for hex in 0..=full[1] {
            for fuc in 0..=full[2] {
                for neuac in 0..=full[3] {
                    for neugc in 0..=full[4] {
                        let sub = GlycanComposition([hexnac, hex, fuc, neuac, neugc]);
                        let size: u8 = sub.0.iter().sum();
                        // Only oligomannose-type parents have a Hex ladder on
                        // the chitobiose trunk: in complex and hybrid glycans
                        // extra Hex sits on the HexNAc antennae.
                        let trunk = full[0] <= 2 && full[3] + full[4] == 0;
                        if sub == *composition
                            || !(size <= SMALL_RESIDUES || total - size <= LARGE_LOSS || trunk)
                            || !plausible(&sub)
                        {
                            continue;
                        }
                        let class = class_of(&sub);
                        let delta = sub.mass() as f32;
                        let decoy_delta = if class == CORE as u8 {
                            delta
                        } else {
                            delta + decoy_shift(composition, key(&sub))
                        };
                        y.push(YIonDef {
                            delta,
                            decoy_delta,
                            class,
                        });
                    }
                }
            }
        }
    }
    let mut oxonium = Vec::new();
    for (ions, present, absent, residue) in [
        (
            NEUAC_OXONIUM,
            NEUAC_PRESENT,
            NEUAC_ABSENT,
            Monosaccharide::NeuAc,
        ),
        (
            NEUGC_OXONIUM,
            NEUGC_PRESENT,
            NEUGC_ABSENT,
            Monosaccharide::NeuGc,
        ),
    ] {
        let class = if composition.count(residue) > 0 {
            present
        } else {
            absent
        };
        for (index, mz) in ions.into_iter().enumerate() {
            let shift = decoy_shift(composition, 0xff00 + residue as u64 * 8 + index as u64);
            oxonium.push((mz, mz + shift, class as u8));
        }
    }
    IonSet { y, oxonium }
}

/// Largest Hex count on the trunk ladder, see [`SpectrumEvidence::ladder`].
pub const LADDER_MAX_HEX: u8 = 14;

/// Offset of the control ladder. It avoids the isotope, NH3, H2O and
/// Fuc/Hex mass differences between real glycan ions.
const LADDER_CONTROL: f32 = 7.37;

/// Hex-ladder evidence on the chitobiose trunk: the Y ions peptide +
/// HexNAc(2)Hex(k) for k = 0..`steps`, and the same rungs at a fixed
/// offset as a control for how dense the spectrum is.
#[derive(Copy, Clone, Debug, Default, PartialEq, Eq)]
pub struct Ladder {
    /// Rungs checked (Hex counts 0..=n, so n + 1).
    pub steps: u8,
    /// Rungs matched.
    pub matched: u8,
    /// Longest run of consecutive matched rungs.
    pub run: u8,
    pub control_matched: u8,
    pub control_run: u8,
}

/// Precomputed spectrum lookups for Y ions at one peptide mass.
pub struct SpectrumEvidence<'a> {
    pub query: &'a ProcessedSpectrum,
    pub tolerance: Tolerance,
    pub max_charge: u8,
}

impl<'a> SpectrumEvidence<'a> {
    pub fn new(query: &'a ProcessedSpectrum, tolerance: Tolerance, precursor_charge: u8) -> Self {
        SpectrumEvidence {
            query,
            tolerance,
            max_charge: precursor_charge.clamp(1, MAX_Y_CHARGE),
        }
    }

    /// Most intense peak matching peptide-plus-`delta` at any charge.
    pub fn y_peak(&self, peptide_mass: f32, delta: f32) -> Option<usize> {
        let neutral = peptide_mass + delta;
        (1..=self.max_charge)
            .filter_map(|charge| {
                let mz = neutral / charge as f32 + PROTON;
                find_at_charge(self.query, mz, charge, self.tolerance)
            })
            .max_by(|&a, &b| self.query.intensities[a].total_cmp(&self.query.intensities[b]))
    }

    fn oxonium_peak(&self, mz: f32) -> Option<usize> {
        find_singly_charged(self.query, mz, self.tolerance)
    }

    /// Random match rates for Y ions (probes at peptide plus a random
    /// 60 Da..`max_delta`) and for sialic-acid oxonium ions (probes at a
    /// random 250-350 m/z). Laplace-smoothed.
    pub fn random_rates(&self, peptide_mass: f32, max_delta: f32) -> (f32, f32) {
        let span = (max_delta - 60.0).max(100.0);
        let mut y_hits = 0;
        for probe in 0..Y_PROBES {
            // Golden-ratio low-discrepancy sequence with a fixed jitter, so
            // probes avoid landing on monosaccharide multiples by design.
            let u = ((probe as f32 + 0.5) * 0.618_034).fract();
            let delta = 60.0 + u * span + 0.37;
            y_hits += usize::from(self.y_peak(peptide_mass, delta).is_some());
        }
        let mut ox_hits = 0;
        for probe in 0..OXONIUM_PROBES {
            let u = ((probe as f32 + 0.5) * 0.618_034).fract();
            ox_hits += usize::from(self.oxonium_peak(250.0 + u * 100.0 + 0.41).is_some());
        }
        (
            (y_hits as f32 + 1.0) / (Y_PROBES as f32 + 2.0),
            (ox_hits as f32 + 1.0) / (OXONIUM_PROBES as f32 + 2.0),
        )
    }

    /// The trunk Hex ladder up to `max_hex` Hex (capped at
    /// [`LADDER_MAX_HEX`]), and its control at a fixed offset.
    pub fn ladder(&self, peptide_mass: f32, max_hex: u8) -> Ladder {
        let steps = max_hex.min(LADDER_MAX_HEX) + 1;
        let trunk = 2.0 * crate::config::HEXNAC as f32;
        let hex = Monosaccharide::Hex.mass() as f32;
        let walk = |offset: f32| {
            let (mut matched, mut run, mut best) = (0u8, 0u8, 0u8);
            for k in 0..steps {
                let delta = trunk + k as f32 * hex + offset;
                if self.y_peak(peptide_mass, delta).is_some() {
                    matched += 1;
                    run += 1;
                    best = best.max(run);
                } else {
                    run = 0;
                }
            }
            (matched, best)
        };
        let (matched, run) = walk(0.0);
        let (control_matched, control_run) = walk(LADDER_CONTROL);
        Ladder {
            steps,
            matched,
            run,
            control_matched,
            control_run,
        }
    }

    /// Class counts for the target and decoy versions of `ions`, plus the
    /// target's matched Y-ion intensity fraction and whether Y0 or Y1 matched.
    pub fn count(&self, peptide_mass: f32, ions: &IonSet) -> (ClassCounts, ClassCounts, f32, bool) {
        let mut target = ClassCounts::default();
        let mut decoy = ClassCounts::default();
        let mut intensity = 0.0;
        let mut anchored = false;
        for (index, ion) in ions.y.iter().enumerate() {
            let class = ion.class as usize;
            target.generated[class] += 1;
            decoy.generated[class] += 1;
            let hit = self.y_peak(peptide_mass, ion.delta);
            if let Some(peak) = hit {
                target.matched[class] += 1;
                intensity += self.query.intensities[peak];
                // Y0 is index 0; Y1 (one HexNAc) is the first sub-composition.
                if index == 0 || (ion.delta - crate::config::HEXNAC as f32).abs() < 1e-3 {
                    anchored = true;
                }
            }
            if ion.decoy_delta == ion.delta {
                decoy.matched[class] += u16::from(hit.is_some());
            } else if self.y_peak(peptide_mass, ion.decoy_delta).is_some() {
                decoy.matched[class] += 1;
            }
        }
        for &(mz, decoy_mz, class) in &ions.oxonium {
            let class = class as usize;
            target.generated[class] += 1;
            decoy.generated[class] += 1;
            target.matched[class] += u16::from(self.oxonium_peak(mz).is_some());
            decoy.matched[class] += u16::from(self.oxonium_peak(decoy_mz).is_some());
        }
        let fraction = if self.query.total_ion_current > 0.0 {
            intensity / self.query.total_ion_current
        } else {
            0.0
        };
        (target, decoy, fraction, anchored)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn spectrum(peaks: &[(f32, u8, f32)]) -> ProcessedSpectrum {
        let mut peaks = peaks.to_vec();
        peaks.sort_by(|a, b| a.0.total_cmp(&b.0));
        ProcessedSpectrum {
            level: 2,
            masses: peaks.iter().map(|p| p.0).collect(),
            charges: peaks.iter().map(|p| p.1).collect(),
            charge_is_known: peaks.iter().map(|_| true).collect(),
            intensities: peaks.iter().map(|p| p.2).collect(),
            total_ion_current: peaks.iter().map(|p| p.2).sum(),
            ..Default::default()
        }
    }

    #[test]
    fn decoy_twin_keeps_core_ions_and_moves_the_rest() {
        let composition = GlycanComposition::parse("HexNAc(4)Hex(5)Fuc(1)NeuAc(2)").unwrap();
        let ions = ion_set(&composition);
        assert!(ions.y.len() > 20);
        for ion in &ions.y {
            let shift = (ion.decoy_delta - ion.delta).abs();
            if ion.class == CORE as u8 {
                assert_eq!(shift, 0.0);
            } else {
                assert!((4.5..=12.5).contains(&shift), "{shift}");
            }
            // The full composition is the precursor, not a Y ion.
            assert!((ion.delta as f64 - composition.mass()).abs() > 1.0);
        }
        assert_eq!(ions.oxonium.len(), 4);
        // Deterministic.
        assert_eq!(ion_set(&composition).y, ions.y);
    }

    #[test]
    fn true_y_ions_favour_the_target() {
        let composition = GlycanComposition::parse("HexNAc(4)Hex(5)Fuc(1)").unwrap();
        let ions = ion_set(&composition);
        let peptide = 1500.0f32;
        // Every target Y ion present, as neutral masses at charge 1.
        let peaks: Vec<_> = ions
            .y
            .iter()
            .map(|ion| (peptide + ion.delta, 1, 100.0))
            .collect();
        let query = spectrum(&peaks);
        let evidence = SpectrumEvidence::new(&query, Tolerance::Ppm(-10.0, 10.0), 3);
        let (target, decoy, fraction, anchored) = evidence.count(peptide, &ions);
        assert_eq!(target.y_matched(), ions.y.len() as u32);
        assert_eq!(target.y_generated(), decoy.y_generated());
        assert_eq!(decoy.matched[CORE], target.matched[CORE]);
        assert!(decoy.y_matched() < target.y_matched() / 2);
        assert!(anchored);
        assert!(fraction > 0.99);
        // Without NeuAc, NeuAc oxonium ions are generated as counter-evidence.
        assert_eq!(target.generated[NEUAC_ABSENT], 2);
        assert_eq!(target.matched[NEUAC_ABSENT], 0);
    }

    #[test]
    fn ladder_counts_rungs_runs_and_control() {
        let peptide = 1500.0f32;
        let rung = |k: u8| {
            peptide
                + 2.0 * crate::config::HEXNAC as f32
                + k as f32 * Monosaccharide::Hex.mass() as f32
        };
        // Rungs Hex 0, 1, 2, 4 and 5 present, 3 missing; one control peak.
        let mut peaks: Vec<_> = [0, 1, 2, 4, 5].map(|k| (rung(k), 1, 100.0)).to_vec();
        peaks.push((rung(1) + LADDER_CONTROL, 1, 50.0));
        let query = spectrum(&peaks);
        let evidence = SpectrumEvidence::new(&query, Tolerance::Ppm(-10.0, 10.0), 2);
        let ladder = evidence.ladder(peptide, 5);
        assert_eq!(
            ladder,
            Ladder {
                steps: 6,
                matched: 5,
                run: 3,
                control_matched: 1,
                control_run: 1,
            }
        );
        assert_eq!(evidence.ladder(peptide, 200).steps, LADDER_MAX_HEX + 1);
    }

    #[test]
    fn random_rates_are_smoothed_probabilities() {
        let query = spectrum(&[(1000.0, 1, 1.0)]);
        let evidence = SpectrumEvidence::new(&query, Tolerance::Ppm(-10.0, 10.0), 2);
        let (y, ox) = evidence.random_rates(1500.0, 3000.0);
        assert!(y > 0.0 && y < 0.05);
        assert!(ox > 0.0 && ox < 0.05);
    }
}
