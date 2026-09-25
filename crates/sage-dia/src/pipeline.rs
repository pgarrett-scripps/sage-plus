//! Turn one DIA run into pseudo-MS2 spectra that Sage can search as DDA.
//!
//! [`pseudo_spectra`] reads MS1 and MS2 spectra of a single file, detects
//! chromatographic hills with koth-core, anchors on every charged MS1
//! isotope feature, and emits one centroided MS2 [`RawSpectrum`] per
//! (precursor, isolation window) pair whose fragments co-elute with it.
//! The output carries a real precursor m/z and charge, so the normal
//! closed search applies (no wide-window mode).

use crate::hills::Channel;
use crate::pseudo::{self, PrecursorTrace, PseudoSettings, PseudoSpectrum};
use rayon::prelude::*;
use sage_core::mass::Tolerance;
use sage_core::spectrum::{Precursor, RawSpectrum, Representation};
use schemars::JsonSchema;
use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;

/// How DIA spectra are searched.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Serialize, Deserialize, JsonSchema)]
#[serde(rename_all = "lowercase")]
pub enum DiaMode {
    /// Search spectra as they are read (DDA, or DIA with `wide_window`).
    #[default]
    Off,
    /// Build MS1-anchored pseudo-MS2 spectra and search them closed.
    Pseudo,
}

/// DIA pseudo-spectrum settings (`"dia"` in the search configuration).
#[derive(Clone, Copy, Debug, PartialEq, Serialize, Deserialize, JsonSchema)]
#[serde(default, deny_unknown_fields)]
pub struct DiaSettings {
    /// `"off"` (default) or `"pseudo"`.
    pub mode: DiaMode,
    /// Minimum Pearson correlation between a fragment hill and its precursor.
    #[schemars(range(min = -1.0, max = 1.0))]
    pub min_corr: f32,
    /// Fragment apex must lie within this many cycles of the precursor apex.
    pub apex_tolerance: u32,
    /// Minimum consecutive MS2 scans for a fragment hill.
    #[schemars(range(min = 1))]
    pub ms2_min_scans: u32,
    /// Drop pseudo-spectra with fewer fragment peaks than this.
    pub min_peaks: u32,
    /// Keep at most this many fragment peaks (by intensity) per pseudo-spectrum.
    #[schemars(range(min = 1))]
    pub max_peaks: u32,
}

impl Default for DiaSettings {
    fn default() -> Self {
        let p = PseudoSettings::default();
        DiaSettings {
            mode: DiaMode::Off,
            min_corr: p.min_corr,
            apex_tolerance: p.apex_tolerance as u32,
            ms2_min_scans: 3,
            min_peaks: p.min_peaks as u32,
            max_peaks: p.max_peaks as u32,
        }
    }
}

impl DiaSettings {
    pub fn is_off(&self) -> bool {
        self.mode == DiaMode::Off
    }

    pub fn validate(&self) -> Result<(), String> {
        if !(self.min_corr.is_finite() && (-1.0..=1.0).contains(&self.min_corr)) {
            return Err("dia.min_corr must be between -1 and 1".into());
        }
        if self.ms2_min_scans == 0 || self.max_peaks == 0 {
            return Err("dia.ms2_min_scans and dia.max_peaks must be at least 1".into());
        }
        Ok(())
    }

    fn pseudo(&self) -> PseudoSettings {
        PseudoSettings {
            apex_tolerance: self.apex_tolerance as i64,
            min_corr: self.min_corr,
            min_peaks: self.min_peaks as usize,
            max_peaks: self.max_peaks as usize,
            ..Default::default()
        }
    }
}

type WindowKey = (i64, i64, i64);

/// Hills of one DIA run: MS1 channel, one channel per isolation window
/// (sorted by window), and the charged MS1 isotope features.
pub struct RunHills {
    pub ms1: Channel,
    pub windows: Vec<Channel>,
    pub precursors: Vec<PrecursorTrace>,
}

/// Isolation window `(target, lower, upper)` of an MS2 spectrum.
pub fn window_of(s: &RawSpectrum) -> Option<(f64, f64, f64)> {
    let p = s.precursors.first()?;
    let (lo, hi) = match p.isolation_window? {
        Tolerance::Da(lo, hi) => (lo, hi),
        _ => return None,
    };
    Some((p.mz as f64, (p.mz + lo) as f64, (p.mz + hi) as f64))
}

fn koth_window(w: (f64, f64, f64)) -> koth_core::IsolationWindow {
    koth_core::IsolationWindow {
        target: w.0,
        lower: w.1,
        upper: w.2,
    }
}

/// koth-core view of a spectrum (positive peaks, sorted by m/z).
pub fn to_koth(s: &RawSpectrum, index: usize) -> koth_core::Spectrum {
    let mut peaks: Vec<koth_core::Peak> =
        s.mz.iter()
            .zip(&s.intensity)
            .filter(|(_, i)| **i > 0.0)
            .map(|(&mz, &intensity)| koth_core::Peak {
                mz,
                intensity,
                ion_mobility: 0.0,
            })
            .collect();
    peaks.sort_by(|a, b| a.mz.total_cmp(&b.mz));
    koth_core::Spectrum {
        scan_index: index,
        retention_time: s.scan_start_time as f64,
        peaks,
        ms_level: s.ms_level,
        isolation_window: if s.ms_level == 2 {
            window_of(s).map(koth_window)
        } else {
            None
        },
        faims_cv: None,
    }
}

/// Detect MS1 hills + isotope features and per-window MS2 hills. `raw` must
/// be sorted by retention time.
pub fn detect_hills(raw: &[RawSpectrum], ms2_min_scans: u32) -> anyhow::Result<RunHills> {
    let kcfg = koth_core::KothConfig {
        hills_ms2: Some(koth_core::config::HillsMs2Overrides {
            min_scans: Some(ms2_min_scans as usize),
            ..Default::default()
        }),
        ..Default::default()
    };

    let ms1_rts: Vec<f32> = raw
        .iter()
        .filter(|s| s.ms_level == 1)
        .map(|s| s.scan_start_time)
        .collect();
    let ms1_hills = koth_core::hills::detect_hills_from_iter(
        raw.iter()
            .filter(|s| s.ms_level == 1)
            .enumerate()
            .map(|(i, s)| to_koth(s, i)),
        &kcfg.hills,
        &kcfg.file,
    );
    let features = koth_core::run_features(&ms1_hills, &kcfg.features, &kcfg.file)?;
    let precursors: Vec<PrecursorTrace> = features
        .iter()
        .filter_map(PrecursorTrace::from_koth)
        .collect();
    drop(features);
    let ms1 = Channel::from_hills(0.0, f64::INFINITY, ms1_rts, ms1_hills);

    // One koth MS2 detector call per isolation window, windows in parallel
    // (koth's own multi-window call runs them serially). BTreeMap keeps the
    // window order, and so pseudo-spectrum order, deterministic.
    let mut by_window: BTreeMap<WindowKey, (koth_core::IsolationWindow, Vec<&RawSpectrum>)> =
        BTreeMap::new();
    for s in raw.iter().filter(|s| s.ms_level == 2) {
        if let Some(w) = window_of(s) {
            let iw = koth_window(w);
            by_window
                .entry(iw.key())
                .or_insert_with(|| (iw, Vec::new()))
                .1
                .push(s);
        }
    }
    let ms2_cfg = kcfg.ms2_hills();
    let windows: Vec<Channel> = by_window
        .into_values()
        .collect::<Vec<_>>()
        .into_par_iter()
        .map(|(iw, spectra)| {
            let rts = spectra.iter().map(|s| s.scan_start_time).collect();
            let hills = koth_core::hills::detect_ms2_hills_from_iter(
                spectra.into_iter().enumerate().map(|(i, s)| to_koth(s, i)),
                &ms2_cfg,
                &kcfg.file,
            );
            Channel::from_hills(iw.lower, iw.upper, rts, hills)
        })
        .collect();
    Ok(RunHills {
        ms1,
        windows,
        precursors,
    })
}

/// MS1-anchored pseudo-spectra for every precursor and window containing it.
pub fn build_pseudo(hills: &RunHills, settings: &PseudoSettings) -> Vec<PseudoSpectrum> {
    hills
        .precursors
        .par_iter()
        .flat_map_iter(|p| {
            let mz = p.mz as f64;
            hills
                .windows
                .iter()
                .enumerate()
                .filter(move |(_, w)| w.lower <= mz && mz <= w.upper)
                .filter_map(move |(i, w)| pseudo::build(p, &hills.ms1, w, i, settings))
        })
        .collect()
}

/// Tier-1 subtraction: mark the fragment hills that explain confidently
/// matched fragment ions. `matches` yields, per identified pseudo-spectrum,
/// its index in `pseudo` and the experimental m/z of its matched fragments
/// (callers drop short, easily shared ions such as b1/b2/y1/y2 first). Returns
/// one `used` flag per hill, per window.
pub fn subtract<'a>(
    hills: &RunHills,
    pseudo: &[PseudoSpectrum],
    matches: impl IntoIterator<Item = (usize, &'a [f32])>,
    ppm: f32,
) -> Vec<Vec<bool>> {
    let mut used: Vec<Vec<bool>> = hills.windows.iter().map(|w| vec![false; w.len()]).collect();
    for (i, mzs) in matches {
        let p = &pseudo[i];
        for &mz in mzs {
            let tol = mz * ppm * 1e-6;
            let lo = p.peaks.partition_point(|x| x.0 < mz - tol);
            let hi = p.peaks.partition_point(|x| x.0 <= mz + tol);
            for k in lo..hi {
                used[p.window][p.hills[k] as usize] = true;
            }
        }
    }
    used
}

/// Tier 2: in every window, group the hills not subtracted by tier 1 by
/// co-elution with an intensity-ordered seed hill. Each group is searched
/// with its isolation window `(lower, upper)` as the precursor tolerance.
pub fn build_tier2(
    hills: &RunHills,
    used: &[Vec<bool>],
    settings: &PseudoSettings,
) -> Vec<(PseudoSpectrum, (f64, f64))> {
    hills
        .windows
        .par_iter()
        .enumerate()
        .flat_map_iter(|(i, w)| {
            pseudo::build_orphans(w, i, &used[i], settings)
                .into_iter()
                .map(move |p| (p, (w.lower, w.upper)))
        })
        .collect()
}

/// Keep only the MS2 peaks that belong to a fragment hill (within `ppm` of a
/// hill spanning that scan). `raw` must be sorted by retention time and be
/// the spectra `hills` was detected from.
pub fn hill_filter(raw: &mut [RawSpectrum], hills: &RunHills, ppm: f32) {
    let mut keys: Vec<WindowKey> = raw
        .iter()
        .filter(|s| s.ms_level == 2)
        .filter_map(|s| window_of(s).map(|w| koth_window(w).key()))
        .collect();
    keys.sort_unstable();
    keys.dedup();
    let mut cycle = vec![0u32; keys.len()];
    let mut jobs = Vec::new();
    for s in raw.iter_mut().filter(|s| s.ms_level == 2) {
        let Some(w) = window_of(s) else { continue };
        let k = keys
            .binary_search(&koth_window(w).key())
            .expect("known window");
        jobs.push((s, k, cycle[k]));
        cycle[k] += 1;
    }
    jobs.into_par_iter().for_each(|(s, k, c)| {
        let Some(channel) = hills.windows.get(k) else {
            s.mz.clear();
            s.intensity.clear();
            return;
        };
        let keep: Vec<bool> =
            s.mz.iter()
                .map(|&mz| {
                    channel
                        .find(mz, ppm)
                        .iter()
                        .any(|h| h.start <= c && c <= h.end())
                })
                .collect();
        let mut it = keep.iter();
        s.mz.retain(|_| *it.next().unwrap());
        let mut it = keep.iter();
        s.intensity.retain(|_| *it.next().unwrap());
        s.total_ion_current = s.intensity.iter().sum();
    });
}

/// A searchable centroided MS2 spectrum for a pseudo-spectrum. With
/// `window`, the isolation window is recorded on the precursor.
pub fn to_raw(
    file_id: usize,
    i: usize,
    p: &PseudoSpectrum,
    window: Option<(f64, f64)>,
) -> RawSpectrum {
    let mut raw = RawSpectrum::default_with_file_id(file_id);
    raw.ms_level = 2;
    raw.id = format!("pseudo={i} window={}", p.window);
    raw.scan_start_time = p.rt;
    raw.representation = Representation::Centroid;
    raw.precursors = vec![Precursor {
        mz: p.precursor_mz,
        intensity: Some(p.precursor_intensity),
        charge: (p.charge > 0).then_some(p.charge),
        isolation_window: window
            .map(|(lo, hi)| Tolerance::Da(lo as f32 - p.precursor_mz, hi as f32 - p.precursor_mz)),
        ..Default::default()
    }];
    raw.mz = p.peaks.iter().map(|x| x.0).collect();
    raw.intensity = p.peaks.iter().map(|x| x.1).collect();
    raw.total_ion_current = raw.intensity.iter().sum();
    raw
}

/// Prepare the spectra of one DIA file for searching according to
/// `settings.mode` (unchanged when off). MS1 spectra are kept when
/// `keep_ms1` is set.
pub fn prepare(
    raw: Vec<RawSpectrum>,
    file_id: usize,
    keep_ms1: bool,
    settings: &DiaSettings,
) -> anyhow::Result<Vec<RawSpectrum>> {
    match settings.mode {
        DiaMode::Off => Ok(raw),
        DiaMode::Pseudo => pseudo_spectra(raw, file_id, keep_ms1, settings),
    }
}

/// Replace the MS2 spectra of one DIA file with MS1-anchored pseudo-MS2
/// spectra. MS1 spectra are returned too when `keep_ms1` is set. Spectra of
/// other MS levels are dropped.
pub fn pseudo_spectra(
    mut raw: Vec<RawSpectrum>,
    file_id: usize,
    keep_ms1: bool,
    settings: &DiaSettings,
) -> anyhow::Result<Vec<RawSpectrum>> {
    raw.sort_by(|a, b| a.scan_start_time.total_cmp(&b.scan_start_time));
    let hills = detect_hills(&raw, settings.ms2_min_scans)?;
    let pseudo = build_pseudo(&hills, &settings.pseudo());
    log::info!(
        "DIA pseudo-spectra: {} from {} MS1 precursor features in {} isolation windows",
        pseudo.len(),
        hills.precursors.len(),
        hills.windows.len()
    );
    if hills.windows.is_empty() {
        log::warn!("DIA pseudo mode: no MS2 isolation windows found; is this a DIA file?");
    }
    drop(hills);
    raw.retain(|s| keep_ms1 && s.ms_level == 1);
    raw.extend(
        pseudo
            .iter()
            .enumerate()
            .map(|(i, p)| to_raw(file_id, i, p, None)),
    );
    Ok(raw)
}

#[cfg(test)]
mod tests {
    use super::*;

    const CYCLES: usize = 30;
    const APEX_A: f32 = 10.0;
    const APEX_B: f32 = 20.0;

    fn gauss(cycle: usize, apex: f32) -> f32 {
        let d = cycle as f32 - apex;
        1e6 * (-d * d / 4.0).exp()
    }

    fn spectrum(ms_level: u8, cycle: usize, peaks: Vec<(f32, f32)>) -> RawSpectrum {
        let mut s = RawSpectrum::default_with_file_id(3);
        s.ms_level = ms_level;
        s.id = format!("scan={}", cycle * 2 + ms_level as usize);
        s.scan_start_time = cycle as f32 * 0.05 + if ms_level == 2 { 0.01 } else { 0.0 };
        s.representation = Representation::Centroid;
        if ms_level == 2 {
            s.precursors = vec![Precursor {
                mz: 500.0,
                isolation_window: Some(Tolerance::Da(-12.5, 12.5)),
                ..Default::default()
            }];
        }
        let peaks: Vec<_> = peaks.into_iter().filter(|p| p.1 > 1.0).collect();
        s.mz = peaks.iter().map(|p| p.0).collect();
        s.intensity = peaks.iter().map(|p| p.1).collect();
        s.total_ion_current = s.intensity.iter().sum();
        s
    }

    /// Two co-isolated 2+ precursors eluting 10 cycles apart, each with its
    /// own fragments, separate into two pseudo-spectra.
    fn synthetic_run() -> Vec<RawSpectrum> {
        let iso = 1.003355 / 2.0;
        let frags_a = [300.1f32, 401.2, 512.3, 613.3, 720.4, 833.5, 901.5];
        let frags_b = [350.2f32, 455.2, 560.3, 671.4, 780.4, 866.5, 950.5];
        let mut raw = Vec::new();
        for c in 0..CYCLES {
            let (a, b) = (gauss(c, APEX_A), gauss(c, APEX_B));
            let ms1 = vec![
                (495.25, a),
                (495.25 + iso, a * 0.6),
                (495.25 + 2.0 * iso, a * 0.25),
                (505.75, b),
                (505.75 + iso, b * 0.6),
                (505.75 + 2.0 * iso, b * 0.25),
            ];
            raw.push(spectrum(1, c, ms1));
            let mut ms2: Vec<(f32, f32)> = frags_a.iter().map(|&mz| (mz, a * 0.3)).collect();
            ms2.extend(frags_b.iter().map(|&mz| (mz, b * 0.3)));
            ms2.sort_by(|x, y| x.0.total_cmp(&y.0));
            raw.push(spectrum(2, c, ms2));
        }
        raw
    }

    #[test]
    fn co_isolated_precursors_separate() {
        let settings = DiaSettings {
            mode: DiaMode::Pseudo,
            ..Default::default()
        };
        let out = pseudo_spectra(synthetic_run(), 3, false, &settings).unwrap();
        assert!(out.iter().all(|s| s.ms_level == 2 && s.file_id == 3));
        let near = |s: &RawSpectrum, mz: f32| (s.precursors[0].mz - mz).abs() < 0.01;
        let a = out.iter().find(|s| near(s, 495.25)).expect("precursor A");
        let b = out.iter().find(|s| near(s, 505.75)).expect("precursor B");
        assert_eq!(a.precursors[0].charge, Some(2));
        assert!(a.mz.iter().any(|&m| (m - 512.3).abs() < 0.01));
        assert!(!a.mz.iter().any(|&m| (m - 560.3).abs() < 0.01));
        assert!(b.mz.iter().any(|&m| (m - 560.3).abs() < 0.01));
        assert!(!b.mz.iter().any(|&m| (m - 512.3).abs() < 0.01));
    }

    #[test]
    fn keeps_ms1_on_request_and_is_deterministic() {
        let settings = DiaSettings {
            mode: DiaMode::Pseudo,
            ..Default::default()
        };
        let one = pseudo_spectra(synthetic_run(), 0, true, &settings).unwrap();
        let two = pseudo_spectra(synthetic_run(), 0, true, &settings).unwrap();
        assert_eq!(one.iter().filter(|s| s.ms_level == 1).count(), CYCLES);
        let ids = |v: &[RawSpectrum]| {
            v.iter()
                .map(|s| (s.id.clone(), s.mz.clone()))
                .collect::<Vec<_>>()
        };
        assert_eq!(ids(&one), ids(&two));
    }

    #[test]
    fn settings_parse_and_default_off() {
        let s: DiaSettings =
            serde_json::from_str(r#"{"mode": "pseudo", "min_corr": 0.6}"#).unwrap();
        assert_eq!(s.mode, DiaMode::Pseudo);
        assert_eq!(s.min_corr, 0.6);
        assert_eq!(s.apex_tolerance, 2);
        assert!(DiaSettings::default().is_off());
        assert!(serde_json::from_str::<DiaSettings>(r#"{"mode": "wide"}"#).is_err());
        assert!(serde_json::from_str::<DiaSettings>(r#"{"q3": true}"#).is_err());
    }
}
