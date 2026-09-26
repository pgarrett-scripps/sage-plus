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
    /// timsTOF only: fragment hills must lie within this ion mobility (1/K0)
    /// of the precursor feature's apex mobility.
    #[schemars(range(min = 0.0))]
    pub im_tolerance: f32,
    /// Orbitrap only: also seed charge-2 pseudo-spectra from MS1 hills that no
    /// isotope feature claimed (at least 5 scans, apex above the median of
    /// such hills). timsTOF ignores it.
    pub hill_precursors: bool,
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
            im_tolerance: p.im_tolerance,
            hill_precursors: true,
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
        if !(self.im_tolerance.is_finite() && self.im_tolerance >= 0.0) {
            return Err("dia.im_tolerance must be a non-negative number".into());
        }
        Ok(())
    }

    pub fn pseudo(&self) -> PseudoSettings {
        PseudoSettings {
            apex_tolerance: self.apex_tolerance as i64,
            min_corr: self.min_corr,
            min_peaks: self.min_peaks as usize,
            max_peaks: self.max_peaks as usize,
            im_tolerance: self.im_tolerance,
            ..Default::default()
        }
    }
}

type WindowKey = (i64, i64, i64);

/// Single-hill precursors: charge, minimum MS1 scans, and minimum apex
/// intensity as a quantile of the unclaimed hills. Tuned on the Orbitrap AIF
/// E. coli run (DIA_SEARCH.md section 12): charge 3 and weak or short hills
/// add more decoy competition than peptides.
const HILL_CHARGE: u8 = 2;
const HILL_MIN_SCANS: usize = 5;
const HILL_MIN_QUANTILE: f64 = 0.5;

/// MS1 hills that no charged isotope feature claimed, taken as charge-2
/// monoisotopic precursors. Only hills with at least [`HILL_MIN_SCANS`]
/// scans and an apex at or above the median of the unclaimed hills qualify.
pub fn single_hill_precursors(
    hills: &[koth_core::Hill],
    features: &[koth_core::Feature],
) -> Vec<PrecursorTrace> {
    let used: std::collections::HashSet<u64> = features
        .iter()
        .filter(|f| f.charge != 0)
        .flat_map(|f| f.hills.iter().map(|h| h.hill_id))
        .collect();
    let free: Vec<&koth_core::Hill> = hills
        .iter()
        .filter(|h| !used.contains(&h.hill_id))
        .collect();
    let mut apex: Vec<f64> = free.iter().map(|h| h.intensity_max).collect();
    apex.sort_by(|a, b| a.total_cmp(b));
    let Some(&min_int) = apex.get(((apex.len().max(1) - 1) as f64 * HILL_MIN_QUANTILE) as usize)
    else {
        return Vec::new();
    };
    free.into_iter()
        .filter(|h| h.n_scans >= HILL_MIN_SCANS && h.intensity_max >= min_int)
        .filter_map(|h| PrecursorTrace::from_hill(h, HILL_CHARGE))
        .collect()
}

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
/// With `hill_precursors`, strong MS1 hills that are not part of any isotope
/// feature become precursors too (see [`single_hill_precursors`]).
pub fn detect_hills(
    raw: &[RawSpectrum],
    ms2_min_scans: u32,
    hill_precursors: bool,
) -> anyhow::Result<RunHills> {
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
    let mut precursors: Vec<PrecursorTrace> = features
        .iter()
        .filter_map(PrecursorTrace::from_koth)
        .collect();
    if hill_precursors {
        let extra = single_hill_precursors(&ms1_hills, &features);
        log::info!("DIA: {} single-hill precursors", extra.len());
        precursors.extend(extra);
    }
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
            let (mz, im) = (p.mz as f64, p.im as f64);
            hills
                .windows
                .iter()
                .enumerate()
                .filter(move |(_, w)| w.contains(mz, im))
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
        inverse_ion_mobility: (p.im > 0.0).then_some(p.im),
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
    let hills = detect_hills(&raw, settings.ms2_min_scans, settings.hill_precursors)?;
    raw.retain(|s| keep_ms1 && s.ms_level == 1);
    raw.extend(from_hills(hills, file_id, settings));
    Ok(raw)
}

/// Pseudo-MS2 spectra of a timsTOF diaPASEF `.d` directory (see
/// [`crate::tims`]). MS1 spectra are not included.
pub fn pseudo_spectra_tdf(
    path: &std::path::Path,
    file_id: usize,
    settings: &DiaSettings,
) -> anyhow::Result<Vec<RawSpectrum>> {
    let hills = crate::tims::detect_hills(path, settings.ms2_min_scans)?;
    Ok(from_hills(hills, file_id, settings))
}

fn from_hills(hills: RunHills, file_id: usize, settings: &DiaSettings) -> Vec<RawSpectrum> {
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
    pseudo
        .iter()
        .enumerate()
        .map(|(i, p)| to_raw(file_id, i, p, None))
        .collect()
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
        let bad = DiaSettings {
            im_tolerance: -0.1,
            ..Default::default()
        };
        assert!(bad.validate().is_err());
    }

    #[test]
    fn ion_mobility_gates_boxes_and_fragments() {
        use crate::pseudo::im_match;
        let b = Channel::from_hills(400.0, 425.0, vec![], vec![]).with_im(0.9, 1.1);
        assert!(b.contains(410.0, 1.0));
        assert!(!b.contains(410.0, 1.2), "outside the box's 1/K0 range");
        assert!(!b.contains(430.0, 1.0), "outside the box's m/z range");
        assert!(b.contains(410.0, 0.0), "no ion mobility: m/z only");
        assert!(im_match(1.00, 1.02, 0.03));
        assert!(!im_match(1.00, 1.05, 0.03));
        assert!(
            im_match(0.0, 1.05, 0.03),
            "IM-less hills are never rejected"
        );
    }

    fn koth_hill(id: u64, mz: f64, n_scans: usize, apex: f32) -> koth_core::Hill {
        let profile: Vec<f32> = (0..n_scans)
            .map(|i| if i == n_scans / 2 { apex } else { apex / 2.0 })
            .collect();
        koth_core::Hill {
            hill_id: id,
            mz,
            mz_std: 0.0,
            mz_se: 0.0,
            rt: 1.0,
            rt_start: 0.9,
            rt_end: 1.1,
            rt_width: 0.2,
            im: 0.0,
            im_std: 0.0,
            scan_start: 4,
            scan_apex: 4 + n_scans / 2,
            scan_end: 4 + n_scans - 1,
            n_scans,
            skipped_scans: 0,
            intensity_sum: profile.iter().map(|&v| v as f64).sum(),
            intensity_max: apex as f64,
            hill_score: 1.0,
            intensity_profile: profile.into(),
            isolation_window: None,
            faims_cv: None,
        }
    }

    #[test]
    fn single_hill_precursors_skip_claimed_short_and_weak_hills() {
        let hills = vec![
            koth_hill(0, 500.0, 8, 1e6), // claimed by the feature
            koth_hill(1, 501.0, 8, 1e6), // claimed by the feature
            koth_hill(2, 600.0, 8, 4e5), // kept
            koth_hill(3, 700.0, 3, 9e5), // too short
            koth_hill(4, 800.0, 8, 1e3), // below the median apex
        ];
        let feature = koth_core::Feature {
            hills: hills[..2].to_vec(),
            charge: 2,
            cosine_score: 1.0,
            ppm_error: 0.0,
        };
        let p = single_hill_precursors(&hills, &[feature]);
        assert_eq!(p.len(), 1);
        assert_eq!(
            (p[0].mz, p[0].charge, p[0].start, p[0].apex),
            (600.0, 2, 4, 8)
        );
        assert_eq!(p[0].intensity, 4e5);
        assert!(single_hill_precursors(&[], &[]).is_empty());
    }
}
