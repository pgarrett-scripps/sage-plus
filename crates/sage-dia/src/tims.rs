//! timsTOF diaPASEF input for pseudo-spectrum mode.
//!
//! Frames are centroided with dnoise-core's watershed centroider in integer
//! `(scan, TOF index)` space, so every centroid keeps its ion mobility. MS2
//! frames are centroided per diaPASEF box (one sub-window of a window group:
//! an m/z × 1/K0 rectangle, from `DiaFrameMsMsWindows`), as dnoise does, so no
//! centroid straddles two boxes.
//!
//! koth then links MS1 centroids into hills by m/z and ion mobility, groups
//! them into isotope features with an ion-mobility apex, and detects MS2 hills
//! per box. MS1 frames are streamed in chunks; MS2 frames are read one window
//! group at a time, so only one group's centroids are in memory at once.

use crate::hills::Channel;
use crate::pipeline::RunHills;
use crate::pseudo::PrecursorTrace;
use anyhow::{anyhow, Context};
use dnoise_core::{
    filter::filter_iterated, halo::horizontal_halo_keep_mask, watershed::watershed_centroid,
    FilterParams, FlatFrame, HaloParams, MsmsFilterParams, WatershedParams,
};
use rayon::prelude::*;
use std::collections::BTreeMap;
use std::path::Path;
use std::sync::Mutex;
use std::time::Instant;
use timsrust::core::{AcquisitionType, Converter, Frame, MSLevel, ScanIndex, TofIndex};
use timsrust::{ImConverter, MzConverter, TimsTofPath};

/// Guard against pathological frames (same cap as dnoise).
const MAX_CENTROIDS: usize = 100_000;
/// MS1 frames centroided in parallel per streamed chunk.
const CHUNK: usize = 64;

/// Is `path` a Bruker timsTOF `.d` directory with TDF data?
pub fn is_tdf(path: &Path) -> bool {
    path.join("analysis.tdf").is_file()
}

struct Centroider {
    mz: MzConverter,
    im: ImConverter,
    ms1_filter: FilterParams,
    ms2_filter: FilterParams,
    halo: HaloParams,
    params: WatershedParams,
}

impl Centroider {
    fn im_of(&self, scan: u32) -> f32 {
        let scan = ScanIndex::try_from(scan).expect("scan index fits u32");
        f64::from(self.im.convert(scan)) as f32
    }

    /// Watershed-centroid the scans `lo..hi` of `frame` into koth peaks
    /// (sorted by m/z, intensity corrected for the frame's accumulation time).
    fn centroid(&self, frame: &Frame, lo: usize, hi: usize, ms2: bool) -> Vec<koth_core::Peak> {
        let ions = frame.ions();
        let offsets = ions.scan_offsets();
        if offsets.len() < 2 {
            return Vec::new();
        }
        let n = offsets.len() - 1;
        let (lo, hi) = (lo.min(n), hi.min(n));
        if hi <= lo {
            return Vec::new();
        }
        let (tofs, ints) = (ions.tof_indices(), ions.intensities());
        // The slice in a local scan origin, so the filters cannot link across
        // a box edge (as dnoise's `filter_per_window`).
        let len = offsets[hi] - offsets[lo];
        let mut flat = FlatFrame {
            frame_id: 0,
            num_scans: hi - lo,
            scan: Vec::with_capacity(len),
            tof: Vec::with_capacity(len),
            intensity: Vec::with_capacity(len),
        };
        for scan in lo..hi {
            for k in offsets[scan]..offsets[scan + 1] {
                flat.scan.push((scan - lo) as u32);
                flat.tof.push(u32::from(tofs[k]));
                flat.intensity.push(u32::from(ints[k]));
            }
        }
        // dnoise stage 1 (vertical-IM filter; its MS/MS knobs for MS2) and
        // stage 2 (horizontal halo), then watershed, all with dnoise defaults.
        let filter = if ms2 {
            &self.ms2_filter
        } else {
            &self.ms1_filter
        };
        let mut keep = filter_iterated(&flat, filter);
        let idx: Vec<usize> = (0..flat.len()).filter(|&i| keep[i]).collect();
        let pick = |v: &[u32]| idx.iter().map(|&i| v[i]).collect::<Vec<u32>>();
        let (s, t, n) = (pick(&flat.scan), pick(&flat.tof), pick(&flat.intensity));
        let halo = horizontal_halo_keep_mask(&s, &t, &n, flat.num_scans, &self.halo);
        for (k, &i) in idx.iter().enumerate() {
            keep[i] &= halo[k];
        }
        let points: Vec<(u32, u32, u32)> = (0..flat.len())
            .filter(|&i| keep[i])
            .map(|i| (flat.scan[i] + lo as u32, flat.tof[i], flat.intensity[i]))
            .collect();
        drop(flat);
        let factor = frame.info().intensity_correction_factor() as f32;
        let mut peaks: Vec<koth_core::Peak> =
            watershed_centroid(&points, &self.params, MAX_CENTROIDS)
                .into_iter()
                .filter(|p| p.2 > 0)
                .map(|(scan, tof, intensity)| koth_core::Peak {
                    mz: f64::from(
                        self.mz
                            .convert(TofIndex::try_from(tof).expect("TOF index fits u32")),
                    ) as f32,
                    intensity: intensity as f32 * factor,
                    ion_mobility: self.im_of(scan),
                })
                .collect();
        peaks.sort_by(|a, b| a.mz.total_cmp(&b.mz));
        peaks
    }
}

fn spectrum(
    index: usize,
    rt: f64,
    ms_level: u8,
    peaks: Vec<koth_core::Peak>,
    window: Option<koth_core::IsolationWindow>,
) -> koth_core::Spectrum {
    koth_core::Spectrum {
        scan_index: index,
        retention_time: rt,
        peaks,
        ms_level,
        isolation_window: window,
        faims_cv: None,
    }
}

/// One diaPASEF box: scan range and isolation window within a window group.
struct DiaBox {
    scan_lo: usize,
    scan_hi: usize,
    lower: f64,
    upper: f64,
}

/// Detect MS1 hills + isotope features and per-box MS2 hills of a diaPASEF
/// run. Channels carry their box's m/z and 1/K0 bounds, hills and precursors
/// their ion mobility.
pub fn detect_hills(path: &Path, ms2_min_scans: u32) -> anyhow::Result<RunHills> {
    let tims = TimsTofPath::new(path.to_string_lossy())
        .map_err(|e| anyhow!("{}: not a timsTOF .d: {e:?}", path.display()))?;
    let reader = tims
        .frame_reader()
        .map_err(|e| anyhow!("{}: cannot read frames: {e:?}", path.display()))?;
    let ctx = Centroider {
        mz: tims.mz_converter().context("no m/z converter")?,
        im: tims.im_converter().context("no ion mobility converter")?,
        ms1_filter: FilterParams::default(),
        ms2_filter: MsmsFilterParams::default().as_filter_params(),
        halo: HaloParams::default(),
        params: WatershedParams::default(),
    };
    let kcfg = koth_core::KothConfig {
        hills_ms2: Some(koth_core::config::HillsMs2Overrides {
            min_scans: Some(ms2_min_scans as usize),
            ..Default::default()
        }),
        ..Default::default()
    };

    // Frame table without ion data.
    let mut ms1: Vec<(usize, f64)> = Vec::new();
    let mut groups: BTreeMap<u8, Vec<(usize, f64)>> = BTreeMap::new();
    let mut quads = BTreeMap::new();
    for index in reader.iter_indices() {
        let info = reader
            .get_info(index)
            .map_err(|e| anyhow!("frame {index}: {e:?}"))?;
        let rt = info.rt_in_seconds() / 60.0;
        match info.ms_level() {
            MSLevel::MS1 => ms1.push((index, rt)),
            MSLevel::MS2 if info.acquisition_type() == AcquisitionType::DIAPASEF => {
                let g = info.window_group();
                groups.entry(g).or_default().push((index, rt));
                quads
                    .entry(g)
                    .or_insert_with(|| info.quadrupole_settings().clone());
            }
            _ => {}
        }
    }
    ms1.sort_by(|a, b| a.1.total_cmp(&b.1));
    let error: Mutex<Option<String>> = Mutex::new(None);
    let fail = |e: String| {
        error.lock().unwrap().get_or_insert(e);
    };

    // MS1: stream chunks of frames, centroided in parallel, into koth.
    let start = Instant::now();
    let ms1_rts: Vec<f32> = ms1.iter().map(|x| x.1 as f32).collect();
    let ms1_spectra = ms1.chunks(CHUNK).flat_map(|chunk| {
        chunk
            .par_iter()
            .map(|&(index, rt)| {
                let peaks = match reader.get_frame(index) {
                    Ok(frame) => ctx.centroid(&frame, 0, usize::MAX, false),
                    Err(e) => {
                        fail(format!("MS1 frame {index}: {e:?}"));
                        Vec::new()
                    }
                };
                (rt, peaks)
            })
            .collect::<Vec<_>>()
    });
    let ms1_hills = koth_core::hills::detect_hills_from_iter(
        ms1_spectra
            .enumerate()
            .map(|(i, (rt, peaks))| spectrum(i, rt, 1, peaks, None)),
        &kcfg.hills,
        &kcfg.file,
    );
    if let Some(e) = error.lock().unwrap().take() {
        return Err(anyhow!(e));
    }
    let hill_time = start.elapsed();
    let features = koth_core::run_features(&ms1_hills, &kcfg.features, &kcfg.file)?;
    log::info!(
        "timsTOF MS1: {} frames, {} hills in {:.1?}, {} isotope features in {:.1?}",
        ms1_rts.len(),
        ms1_hills.len(),
        hill_time,
        features.len(),
        start.elapsed() - hill_time
    );
    let precursors: Vec<PrecursorTrace> = features
        .iter()
        .filter_map(PrecursorTrace::from_koth)
        .collect();
    drop(features);
    let ms1 = Channel::from_hills(0.0, f64::INFINITY, ms1_rts, ms1_hills);

    // MS2: one window group at a time; its boxes' hills in parallel.
    let ms2_cfg = kcfg.ms2_hills();
    let start = Instant::now();
    let mut windows = Vec::new();
    for (g, mut frames) in groups {
        frames.sort_by(|a, b| a.1.total_cmp(&b.1));
        let quad = &quads[&g];
        let boxes: Vec<DiaBox> = (0..quad.len())
            .map(|k| DiaBox {
                scan_lo: quad.scan_starts[k],
                scan_hi: quad.scan_ends[k],
                lower: f64::from(quad.isolation_windows[k].lower()),
                upper: f64::from(quad.isolation_windows[k].upper()),
            })
            .collect();
        // [frame][box] centroids.
        let per_frame: Vec<Vec<Vec<koth_core::Peak>>> = frames
            .par_iter()
            .map(|&(index, _)| match reader.get_frame(index) {
                Ok(frame) => boxes
                    .iter()
                    .map(|b| ctx.centroid(&frame, b.scan_lo, b.scan_hi, true))
                    .collect(),
                Err(e) => {
                    fail(format!("MS2 frame {index}: {e:?}"));
                    boxes.iter().map(|_| Vec::new()).collect()
                }
            })
            .collect();
        if let Some(e) = error.lock().unwrap().take() {
            return Err(anyhow!(e));
        }
        let rts: Vec<f32> = frames.iter().map(|x| x.1 as f32).collect();
        let mut per_box: Vec<Vec<Vec<koth_core::Peak>>> = boxes
            .iter()
            .map(|_| Vec::with_capacity(frames.len()))
            .collect();
        for frame in per_frame {
            for (k, peaks) in frame.into_iter().enumerate() {
                per_box[k].push(peaks);
            }
        }
        let channels: Vec<Channel> = per_box
            .into_par_iter()
            .zip(boxes.par_iter())
            .map(|(spectra, b)| {
                let iw = koth_core::IsolationWindow {
                    target: (b.lower + b.upper) / 2.0,
                    lower: b.lower,
                    upper: b.upper,
                };
                let hills = koth_core::hills::detect_ms2_hills_from_iter(
                    spectra
                        .into_iter()
                        .enumerate()
                        .map(|(i, peaks)| spectrum(i, rts[i] as f64, 2, peaks, Some(iw))),
                    &ms2_cfg,
                    &kcfg.file,
                );
                let (a, z) = (ctx.im_of(b.scan_lo as u32), ctx.im_of(b.scan_hi as u32));
                Channel::from_hills(b.lower, b.upper, rts.clone(), hills)
                    .with_im(a.min(z) as f64, a.max(z) as f64)
            })
            .collect();
        windows.extend(channels);
    }
    log::info!(
        "timsTOF MS2: {} boxes, {} hills in {:.1?}",
        windows.len(),
        windows.iter().map(|w| w.len()).sum::<usize>(),
        start.elapsed()
    );
    Ok(RunHills {
        ms1,
        windows,
        precursors,
    })
}
