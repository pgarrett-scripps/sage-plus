//! DIA-Umpire-style pseudo-MS2 spectra.
//!
//! Each charged MS1 isotope feature is a precursor. In every isolation window
//! that contains its monoisotopic m/z, fragment hills whose apex lies within
//! `apex_tolerance` cycles of the precursor apex and whose elution profile
//! correlates with the precursor's at `min_corr` or better become the peaks
//! of one pseudo-spectrum. Grouping is deliberately loose: a fragment hill may
//! join any number of pseudo-spectra, and the search and FDR sort out the rest.

use crate::coelution::pearson;
use crate::hills::{Channel, CompactHill};

#[derive(Clone, Copy, Debug)]
pub struct PseudoSettings {
    /// Fragment apex must be within this many window cycles of the precursor apex.
    pub apex_tolerance: i64,
    /// Minimum Pearson correlation between fragment and precursor profiles.
    pub min_corr: f32,
    /// Profiles are compared within this many cycles of the precursor apex.
    pub half_window: i64,
    /// Minimum overlapping cycles for a correlation to count.
    pub min_overlap: usize,
    /// Keep at most this many peaks (by intensity) per pseudo-spectrum.
    pub max_peaks: usize,
    /// Drop pseudo-spectra with fewer peaks than this.
    pub min_peaks: usize,
    /// Fragment and precursor ion mobility (1/K0) must agree within this
    /// (ignored when either has no ion mobility).
    pub im_tolerance: f32,
}

/// Do two ion mobilities agree within `tol`? Always true when either is 0
/// (no ion mobility).
pub fn im_match(a: f32, b: f32, tol: f32) -> bool {
    a == 0.0 || b == 0.0 || (a - b).abs() <= tol
}

impl Default for PseudoSettings {
    fn default() -> Self {
        PseudoSettings {
            apex_tolerance: 2,
            min_corr: 0.5,
            half_window: 6,
            min_overlap: 3,
            max_peaks: 150,
            min_peaks: 6,
            im_tolerance: 0.03,
        }
    }
}

/// An MS1 precursor: monoisotopic m/z, charge and summed isotope profile on
/// MS1 cycles.
pub struct PrecursorTrace {
    pub mz: f32,
    pub charge: u8,
    pub intensity: f32,
    /// Ion mobility (1/K0) at the feature apex; 0 without ion mobility.
    pub im: f32,
    /// MS1 cycle of the first profile value.
    pub start: u32,
    pub apex: u32,
    pub profile: Vec<f32>,
}

impl PrecursorTrace {
    pub fn from_koth(feature: &koth_core::Feature) -> Option<Self> {
        if feature.charge == 0 {
            return None;
        }
        let (start, _, profile) = feature.elution_profile();
        let profile: Vec<f32> = profile.into_iter().map(|v| v as f32).collect();
        let (apex_off, &max) = profile
            .iter()
            .enumerate()
            .max_by(|a, b| a.1.total_cmp(b.1))?;
        Some(PrecursorTrace {
            mz: feature.monoisotopic_mz() as f32,
            charge: feature.charge,
            intensity: max,
            im: feature.im_apex() as f32,
            start: start as u32,
            apex: (start + apex_off) as u32,
            profile,
        })
    }

    fn at(&self, ms1_cycle: i64) -> f32 {
        let i = ms1_cycle - self.start as i64;
        if i < 0 || i >= self.profile.len() as i64 {
            0.0
        } else {
            self.profile[i as usize]
        }
    }
}

pub struct PseudoSpectrum {
    pub precursor_mz: f32,
    pub charge: u8,
    pub rt: f32,
    pub precursor_intensity: f32,
    /// Precursor ion mobility (1/K0); 0 when unknown.
    pub im: f32,
    /// Index of the isolation window the fragments came from.
    pub window: usize,
    /// (m/z, intensity), sorted by m/z.
    pub peaks: Vec<(f32, f32)>,
    /// Channel indices of the fragment hills used as peaks.
    pub hills: Vec<u32>,
}

/// Build the pseudo-spectra of one precursor in one window.
pub fn build(
    precursor: &PrecursorTrace,
    ms1: &Channel,
    window: &Channel,
    window_index: usize,
    settings: &PseudoSettings,
) -> Option<PseudoSpectrum> {
    let apex_rt = ms1.rts[precursor.apex as usize];
    let c = window.cycle_at(apex_rt) as i64;
    let lo = (c - settings.half_window).max(0);
    let hi = (c + settings.half_window).min(window.rts.len() as i64 - 1);
    // Precursor profile resampled on this window's cycles.
    let prec: Vec<f32> = (lo..=hi)
        .map(|j| precursor.at(ms1.cycle_at(window.rts[j as usize]) as i64))
        .collect();

    let mut peaks: Vec<(f32, f32, u32)> = Vec::new();
    let mut a = Vec::with_capacity(prec.len());
    let mut b = Vec::with_capacity(prec.len());
    for (idx, hill) in
        window.apex_between_indexed(c - settings.apex_tolerance, c + settings.apex_tolerance)
    {
        if hill.mz > precursor.mz * precursor.charge as f32
            || !im_match(hill.im, precursor.im, settings.im_tolerance)
        {
            continue;
        }
        let from = lo.max(hill.start as i64);
        let to = hi.min(hill.end() as i64);
        if to - from + 1 < settings.min_overlap as i64 {
            continue;
        }
        a.clear();
        b.clear();
        for j in from..=to {
            a.push(window.intensity(hill, j));
            b.push(prec[(j - lo) as usize]);
        }
        if pearson(&a, &b) >= settings.min_corr {
            peaks.push((hill.mz, apex_intensity(window, hill), idx as u32));
        }
    }
    finish(peaks, settings).map(|(peaks, hills)| PseudoSpectrum {
        precursor_mz: precursor.mz,
        charge: precursor.charge,
        rt: apex_rt,
        precursor_intensity: precursor.intensity,
        im: precursor.im,
        window: window_index,
        peaks,
        hills,
    })
}

type Peaks = (Vec<(f32, f32)>, Vec<u32>);

fn finish(mut peaks: Vec<(f32, f32, u32)>, settings: &PseudoSettings) -> Option<Peaks> {
    if peaks.len() < settings.min_peaks {
        return None;
    }
    if peaks.len() > settings.max_peaks {
        peaks.sort_by(|x, y| y.1.total_cmp(&x.1));
        peaks.truncate(settings.max_peaks);
    }
    peaks.sort_by(|x, y| x.0.total_cmp(&y.0));
    Some((
        peaks.iter().map(|p| (p.0, p.1)).collect(),
        peaks.iter().map(|p| p.2).collect(),
    ))
}

/// DIA-Umpire Q3-style fallback: fragment groups with no MS1 feature.
///
/// Hills not used by any MS1-anchored pseudo-spectrum (`used`) are seeded in
/// decreasing intensity; each seed collects other unused hills whose apex is
/// within `apex_tolerance` and whose profile correlates with the seed's. A
/// grouped hill is not reused as a seed. The result has no precursor mass or
/// charge (`precursor_mz` is the window centre, `charge` 0): search it with
/// the isolation window as the precursor tolerance.
pub fn build_orphans(
    window: &Channel,
    window_index: usize,
    used: &[bool],
    settings: &PseudoSettings,
) -> Vec<PseudoSpectrum> {
    let hills = window.hills();
    let mut order: Vec<usize> = (0..hills.len()).filter(|&i| !used[i]).collect();
    order.sort_by(|&x, &y| hills[y].intensity_max.total_cmp(&hills[x].intensity_max));
    let mut taken = vec![false; hills.len()];
    let mut out = Vec::new();
    let (mut a, mut b) = (Vec::new(), Vec::new());
    for seed_idx in order {
        if taken[seed_idx] {
            continue;
        }
        let seed = &hills[seed_idx];
        let c = seed.apex as i64;
        let lo = (c - settings.half_window).max(0);
        let hi = (c + settings.half_window).min(window.rts.len() as i64 - 1);
        let mut peaks = Vec::new();
        let mut members = Vec::new();
        for (idx, hill) in
            window.apex_between_indexed(c - settings.apex_tolerance, c + settings.apex_tolerance)
        {
            if used[idx] || taken[idx] || !im_match(hill.im, seed.im, settings.im_tolerance) {
                continue;
            }
            let from = lo.max(hill.start as i64).max(seed.start as i64);
            let to = hi.min(hill.end() as i64).min(seed.end() as i64);
            if to - from + 1 < settings.min_overlap as i64 {
                continue;
            }
            a.clear();
            b.clear();
            for j in from..=to {
                a.push(window.intensity(hill, j));
                b.push(window.intensity(seed, j));
            }
            if idx == seed_idx || pearson(&a, &b) >= settings.min_corr {
                peaks.push((hill.mz, apex_intensity(window, hill), idx as u32));
                members.push(idx);
            }
        }
        if let Some((peaks, hill_ids)) = finish(peaks, settings) {
            for m in members {
                taken[m] = true;
            }
            out.push(PseudoSpectrum {
                precursor_mz: ((window.lower + window.upper) / 2.0) as f32,
                charge: 0,
                rt: window.rts[c.clamp(0, window.rts.len() as i64 - 1) as usize],
                precursor_intensity: 0.0,
                im: seed.im,
                window: window_index,
                peaks,
                hills: hill_ids,
            });
        }
    }
    out
}

fn apex_intensity(window: &Channel, hill: &CompactHill) -> f32 {
    window.intensity(hill, hill.apex as i64)
}
