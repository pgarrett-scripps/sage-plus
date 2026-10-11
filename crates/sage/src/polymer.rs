//! Polymer contamination in MS1 spectra.
//!
//! Detergents and plastics ionize as ladders of peaks spaced by one repeat
//! unit. A spectrum's polymer signal is the intensity of peaks that belong to
//! a run of consecutive ladder members, for any charge state and adduct.
//! Masses are monoisotopic, from the elemental formulas of the repeat units.

use serde::{Deserialize, Serialize};

/// Monoisotopic mass of water, the end groups of linear PEG and PPG.
const WATER: f64 = 18.010_565;

/// A polymer series: `end_group + n * repeat` for `n >= 1`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Polymer {
    pub name: &'static str,
    /// Neutral mass of the non-repeating part (zero for cyclic polymers).
    pub end_group: f64,
    /// Neutral mass of one repeat unit.
    pub repeat: f64,
}

/// Polyethylene glycol, H-(C2H4O)n-OH.
pub const PEG: Polymer = Polymer {
    name: "peg",
    end_group: WATER,
    repeat: 44.026_215,
};

/// Polypropylene glycol, H-(C3H6O)n-OH.
pub const PPG: Polymer = Polymer {
    name: "ppg",
    end_group: WATER,
    repeat: 58.041_865,
};

/// Cyclic polydimethylsiloxane, (C2H6OSi)n, e.g. the 445.12 lock mass.
pub const POLYSILOXANE: Polymer = Polymer {
    name: "polysiloxane",
    end_group: 0.0,
    repeat: 74.018_792,
};

pub const DEFAULT_POLYMERS: [Polymer; 3] = [PEG, PPG, POLYSILOXANE];

/// Charge carriers: proton, sodium, and ammonium cations.
const ADDUCTS: [f64; 3] = [1.007_276, 22.989_218, 18.033_823];

/// Default MS1 matching tolerance.
pub const DEFAULT_TOLERANCE_PPM: f64 = 10.0;

/// Consecutive ladder members required before any of them count.
pub const DEFAULT_MIN_LADDER: usize = 4;

/// Highest charge state searched.
pub const MAX_CHARGE: u8 = 3;

#[derive(Debug, Clone)]
pub struct PolymerScanner {
    pub polymers: Vec<Polymer>,
    pub tolerance_ppm: f64,
    pub min_ladder: usize,
}

impl Default for PolymerScanner {
    fn default() -> Self {
        Self {
            polymers: DEFAULT_POLYMERS.to_vec(),
            tolerance_ppm: DEFAULT_TOLERANCE_PPM,
            min_ladder: DEFAULT_MIN_LADDER,
        }
    }
}

/// Index of the most intense peak within `tolerance_ppm` of `target`.
/// `mz` must be sorted ascending.
fn best_peak(mz: &[f32], intensity: &[f32], target: f64, tolerance_ppm: f64) -> Option<usize> {
    let delta = target * tolerance_ppm * 1e-6;
    let (low, high) = (target - delta, target + delta);
    let start = mz.partition_point(|&value| (value as f64) < low);
    (start..mz.len())
        .take_while(|&index| mz[index] as f64 <= high)
        .max_by(|&left, &right| intensity[left].total_cmp(&intensity[right]))
}

impl PolymerScanner {
    /// Intensity of peaks on each polymer's ladders, in the order of
    /// `self.polymers`. `mz` must be sorted ascending and parallel to
    /// `intensity`. A peak counts once per polymer, even when several charge
    /// states or adducts claim it.
    pub fn scan(&self, mz: &[f32], intensity: &[f32]) -> Vec<f64> {
        let (Some(&first), Some(&last)) = (mz.first(), mz.last()) else {
            return vec![0.0; self.polymers.len()];
        };
        // Widen the peak range by the tolerance so ladder ends are kept.
        let slack = self.tolerance_ppm * 1e-6;
        let (first, last) = (first as f64 * (1.0 - slack), last as f64 * (1.0 + slack));
        let mut claimed = vec![false; mz.len()];
        let mut run = Vec::new();
        self.polymers
            .iter()
            .map(|polymer| {
                claimed.iter_mut().for_each(|flag| *flag = false);
                for charge in 1..=MAX_CHARGE {
                    let z = charge as f64;
                    for adduct in ADDUCTS {
                        let offset = polymer.end_group + z * adduct;
                        let start = (((first * z - offset) / polymer.repeat).ceil() as i64).max(1);
                        let mut n = start as f64;
                        run.clear();
                        loop {
                            let target = (offset + n * polymer.repeat) / z;
                            let hit = if target > last {
                                None
                            } else {
                                best_peak(mz, intensity, target, self.tolerance_ppm)
                            };
                            match hit {
                                Some(index) => run.push(index),
                                None => {
                                    if run.len() >= self.min_ladder {
                                        run.iter().for_each(|&index| claimed[index] = true);
                                    }
                                    run.clear();
                                }
                            }
                            if target > last {
                                break;
                            }
                            n += 1.0;
                        }
                    }
                }
                claimed
                    .iter()
                    .zip(intensity)
                    .filter(|(flag, _)| **flag)
                    .map(|(_, &value)| value as f64)
                    .sum()
            })
            .collect()
    }
}

/// Polymer signal of one file's MS1 spectra.
#[derive(Debug, Clone, Default, PartialEq, Serialize, Deserialize)]
pub struct PolymerFileStats {
    /// Centroided MS1 spectra scanned.
    pub ms1_spectra: usize,
    /// Profile MS1 spectra, which are skipped.
    pub skipped_profile_spectra: usize,
    /// Summed intensity of all scanned MS1 peaks.
    pub total_ion_current: f64,
    /// Percent of `total_ion_current` on each polymer's ladders.
    pub polymers: Vec<PolymerShare>,
}

#[derive(Debug, Clone, Default, PartialEq, Serialize, Deserialize)]
pub struct PolymerShare {
    pub name: String,
    pub intensity: f64,
    pub tic_pct: f64,
}

impl PolymerFileStats {
    pub fn new(polymers: &[Polymer]) -> Self {
        Self {
            polymers: polymers
                .iter()
                .map(|polymer| PolymerShare {
                    name: polymer.name.into(),
                    ..Default::default()
                })
                .collect(),
            ..Default::default()
        }
    }

    /// Add one spectrum's TIC and per-polymer intensities from
    /// [`PolymerScanner::scan`].
    pub fn add(&mut self, total_ion_current: f64, intensities: &[f64]) {
        self.ms1_spectra += 1;
        self.total_ion_current += total_ion_current;
        for (share, intensity) in self.polymers.iter_mut().zip(intensities) {
            share.intensity += intensity;
        }
    }

    pub fn merge(&mut self, other: &Self) {
        self.ms1_spectra += other.ms1_spectra;
        self.skipped_profile_spectra += other.skipped_profile_spectra;
        self.total_ion_current += other.total_ion_current;
        for (share, other) in self.polymers.iter_mut().zip(&other.polymers) {
            share.intensity += other.intensity;
        }
    }

    /// Fill in `tic_pct` from the accumulated intensities.
    pub fn finish(&mut self) {
        for share in &mut self.polymers {
            share.tic_pct = if self.total_ion_current > 0.0 {
                100.0 * share.intensity / self.total_ion_current
            } else {
                0.0
            };
        }
    }
}

#[cfg(test)]
#[path = "../tests/unit/polymer.rs"]
mod test;
