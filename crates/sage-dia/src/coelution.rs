//! Co-elution features for one candidate peptide in one isolation window.
//!
//! Given the candidate's theoretical fragment m/z values, the cycle of the
//! spectrum that proposed it, and its precursor m/z, look up fragment hills in
//! the window and the precursor hill in MS1, then measure agreement:
//! how many fragments have a hill at that time, how well the strongest
//! fragment profiles correlate, how tightly their apexes agree, and whether
//! the precursor hill has the same elution profile. No dense XIC matrix is
//! built: each candidate reads a few short hill slices.

use crate::hills::{Channel, CompactHill};

#[derive(Clone, Copy, Debug)]
pub struct CoelutionSettings {
    pub fragment_ppm: f32,
    pub precursor_ppm: f32,
    /// Fragments whose profiles are summed into the reference elution profile.
    pub top_k: usize,
    /// Half-width, in cycles, of the profile slice compared around the spectrum.
    pub half_window: i64,
    /// Correlation to the reference profile that counts a fragment as co-eluting.
    pub min_corr: f32,
    /// A fragment hill must reach within this many cycles of the spectrum.
    pub reach: i64,
}

impl Default for CoelutionSettings {
    fn default() -> Self {
        CoelutionSettings {
            fragment_ppm: 15.0,
            precursor_ppm: 10.0,
            top_k: 6,
            half_window: 8,
            min_corr: 0.7,
            reach: 1,
        }
    }
}

#[derive(Clone, Copy, Debug, Default)]
pub struct CoelutionFeatures {
    /// Theoretical fragments looked up.
    pub n_fragments: u32,
    /// Fragments with a hill spanning the spectrum's cycle (within `reach`).
    pub n_with_hill: u32,
    /// Fragment hills correlating with the reference profile at `min_corr`.
    pub n_coeluting: u32,
    /// Fragment hills whose apex is within 2 cycles of the median top-k apex.
    pub n_coapex: u32,
    /// Mean pairwise Pearson correlation among the top-k fragment profiles.
    pub frag_corr: f32,
    /// Mean |apex - median apex| of the top-k fragments, in cycles.
    pub apex_spread: f32,
    /// |median top-k apex - spectrum cycle|, in cycles.
    pub apex_offset: f32,
    /// Share of summed top-k intensity at the spectrum cycle vs the apex cycle.
    pub apex_fraction: f32,
    pub ms1_present: bool,
    /// Correlation between the precursor hill and the reference fragment profile.
    pub ms1_corr: f32,
    /// |precursor apex RT - fragment apex RT|, minutes.
    pub ms1_apex_delta: f32,
}

pub(crate) fn pearson(a: &[f32], b: &[f32]) -> f32 {
    let n = a.len() as f32;
    let ma = a.iter().sum::<f32>() / n;
    let mb = b.iter().sum::<f32>() / n;
    let (mut sab, mut saa, mut sbb) = (0.0f32, 0.0f32, 0.0f32);
    for (x, y) in a.iter().zip(b) {
        let (dx, dy) = (x - ma, y - mb);
        sab += dx * dy;
        saa += dx * dx;
        sbb += dy * dy;
    }
    if saa <= 0.0 || sbb <= 0.0 {
        0.0
    } else {
        sab / (saa * sbb).sqrt()
    }
}

/// Pick, among hills at a fragment's m/z, the one strongest near `cycle`.
fn best_hill<'a>(
    channel: &Channel,
    candidates: &'a [CompactHill],
    cycle: i64,
    reach: i64,
) -> Option<(&'a CompactHill, f32)> {
    candidates
        .iter()
        .filter(|h| h.start as i64 <= cycle + reach && h.end() as i64 >= cycle - reach)
        .map(|h| {
            let near = (cycle - reach..=cycle + reach)
                .map(|c| channel.intensity(h, c))
                .fold(0.0f32, f32::max);
            (h, near)
        })
        .filter(|(_, near)| *near > 0.0)
        .max_by(|a, b| a.1.total_cmp(&b.1))
}

pub fn score(
    fragment_mzs: &[f32],
    window: &Channel,
    cycle: usize,
    precursor_mz: f32,
    ms1: &Channel,
    settings: &CoelutionSettings,
) -> CoelutionFeatures {
    let cycle = cycle as i64;
    let mut out = CoelutionFeatures {
        n_fragments: fragment_mzs.len() as u32,
        ..Default::default()
    };

    // One hill per fragment; a hill claimed by two fragments counts once.
    let mut matched: Vec<(&CompactHill, f32)> = Vec::with_capacity(fragment_mzs.len());
    for &mz in fragment_mzs {
        if let Some(hit) = best_hill(
            window,
            window.find(mz, settings.fragment_ppm),
            cycle,
            settings.reach,
        ) {
            if !matched.iter().any(|(h, _)| std::ptr::eq(*h, hit.0)) {
                matched.push(hit);
            }
        }
    }
    out.n_with_hill = matched.len() as u32;
    if matched.is_empty() {
        return out;
    }
    matched.sort_by(|a, b| b.1.total_cmp(&a.1));

    let lo = (cycle - settings.half_window).max(0);
    let hi = (cycle + settings.half_window).min(window.rts.len() as i64 - 1);
    let slice = |ch: &Channel, h: &CompactHill| -> Vec<f32> {
        (lo..=hi).map(|c| ch.intensity(h, c)).collect()
    };

    let top: Vec<Vec<f32>> = matched
        .iter()
        .take(settings.top_k)
        .map(|(h, _)| slice(window, h))
        .collect();
    let mut reference = vec![0.0f32; top[0].len()];
    for p in &top {
        for (r, v) in reference.iter_mut().zip(p) {
            *r += v;
        }
    }

    if top.len() >= 2 {
        let mut sum = 0.0;
        let mut pairs = 0;
        for i in 0..top.len() {
            for j in i + 1..top.len() {
                sum += pearson(&top[i], &top[j]);
                pairs += 1;
            }
        }
        out.frag_corr = sum / pairs as f32;
    }

    // Leave-one-out reference for each top-k fragment, so a fragment is not
    // credited for correlating with itself; the rest compare to the full sum.
    out.n_coeluting = matched
        .iter()
        .enumerate()
        .filter(|(i, (h, _))| {
            let p = slice(window, h);
            let r: Vec<f32> = if *i < top.len() {
                reference.iter().zip(&p).map(|(r, v)| r - v).collect()
            } else {
                reference.clone()
            };
            pearson(&p, &r) >= settings.min_corr
        })
        .count() as u32;

    let mut apexes: Vec<i64> = matched
        .iter()
        .take(settings.top_k)
        .map(|(h, _)| h.apex as i64)
        .collect();
    apexes.sort_unstable();
    let median = apexes[apexes.len() / 2];
    out.apex_spread =
        apexes.iter().map(|a| (a - median).abs()).sum::<i64>() as f32 / apexes.len() as f32;
    out.apex_offset = (median - cycle).abs() as f32;
    out.n_coapex = matched
        .iter()
        .filter(|(h, _)| (h.apex as i64 - median).abs() <= 2)
        .count() as u32;
    let at = |c: i64| -> f32 {
        matched
            .iter()
            .take(settings.top_k)
            .map(|(h, _)| window.intensity(h, c))
            .sum()
    };
    let apex_total = at(median);
    if apex_total > 0.0 {
        out.apex_fraction = (at(cycle) / apex_total).min(1.0);
    }

    // Precursor hill: the strongest MS1 hill at the precursor m/z around the
    // spectrum's time, resampled onto this window's cycles.
    let rt = window.rts[cycle as usize];
    let ms1_cycle = ms1.cycle_at(rt) as i64;
    if let Some((hill, _)) = best_hill(
        ms1,
        ms1.find(precursor_mz, settings.precursor_ppm),
        ms1_cycle,
        settings.reach,
    ) {
        out.ms1_present = true;
        let profile: Vec<f32> = (lo..=hi)
            .map(|c| ms1.intensity(hill, ms1.cycle_at(window.rts[c as usize]) as i64))
            .collect();
        out.ms1_corr = pearson(&profile, &reference);
        let frag_apex_rt = window.rts[median.clamp(0, window.rts.len() as i64 - 1) as usize];
        out.ms1_apex_delta = (ms1.rts[hill.apex as usize] - frag_apex_rt).abs();
    }
    out
}
