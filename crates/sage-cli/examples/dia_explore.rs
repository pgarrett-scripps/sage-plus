//! DIA spike (docs/explore/DIA_SEARCH.md): does fragment-hill co-elution
//! separate targets from decoys among wide-window search candidates?
//!
//! 1. Read a DIA Thermo `.raw` file with Sage's reader.
//! 2. Detect MS1 hills and per-window MS2 hills with koth-core (library, no files).
//! 3. Run Sage's wide-window chimeric search on the MS2 spectra; LDA + q-values.
//! 4. For every candidate, compute co-elution features from the hills, then
//!    compare targets and decoys and refit LDA with the new features.
//!
//! Modes (run each as its own process so peak RSS is per mode):
//! * `wide`: Sage wide-window chimeric search only (the baseline).
//! * `pseudo`: koth hills -> DIA-Umpire-style pseudo-MS2 spectra -> ordinary
//!   closed DDA-style Sage search. Each pseudo-spectrum is anchored on a koth
//!   MS1 isotope feature (monoisotopic m/z, charge, apex RT) and searched at
//!   +-10 ppm on that mass with isotope errors -1..=1, not with the isolation
//!   window as the precursor tolerance. `--q3` adds a DIA-Umpire Q3-style
//!   tier: fragment groups with no MS1 feature, searched wide-window.
//! * `tiered`: `pseudo` (tier 1), then subtract the fragment hills matched by
//!   tier-1 PSMs at 1% FDR (ions b/y 3+ only), group the remaining hills per
//!   window by loose co-elution (tier 2) and search those wide-window. Reports
//!   tier 1 alone, tiers FDR-controlled separately, and one LDA with the tier
//!   as a feature.
//! * `hillfilter`: `wide` on raw scans that keep only peaks on a fragment hill.
//!
//! `--raw` may also be a timsTOF diaPASEF `.d` directory (modes `wide`,
//! `tiered`, `hillfilter`). Hills then come from `sage_dia::tims` (dnoise
//! watershed centroids, ion-mobility-aware hills per m/z x 1/K0 box), and the
//! wide-window scans from Sage's Bruker reader. Tolerances are 15 ppm
//! precursor / 20 ppm fragment instead of 10 / 15.
//! * `rescore`: `wide` plus fragment-hill co-elution features for every
//!   candidate, and LDA refits with those features.
//!
//! ```text
//! cargo run --release -p sage-cli --example dia_explore -- \
//!     --raw FILE.raw --fasta FILE.fasta --out DIR [--mode wide|pseudo|tiered|hillfilter|rescore]
//! ```

use rayon::prelude::*;
use sage_core::database::IndexedDatabase;
use sage_core::database::Parameters;
use sage_core::ion_series::{IonSeries, Kind};
use sage_core::mass::{Tolerance, PROTON};
use sage_core::ml::linear_discriminant::LinearDiscriminantAnalysis;
use sage_core::scoring::{Feature, ScoreType, Scorer};
use sage_core::spectrum::ProcessedSpectrum;
use sage_core::spectrum::{RawSpectrum, SpectrumProcessor};
use sage_dia::coelution::{self, CoelutionFeatures, CoelutionSettings};
use sage_dia::hills::Channel;
use sage_dia::pipeline::{self as dia, to_koth, window_of};
use sage_dia::pseudo::{self, PrecursorTrace, PseudoSettings};
use std::collections::HashMap;
use std::io::Write;
use std::sync::atomic::{AtomicBool, Ordering};
use std::time::Instant;

/// Is the input a timsTOF `.d` (set once in `main`)?
static TIMS: AtomicBool = AtomicBool::new(false);

struct Args {
    raw: String,
    fasta: String,
    out: String,
    report_psms: usize,
    ms2_min_scans: usize,
    mode: String,
    min_corr: f32,
    apex_tolerance: i64,
    /// Also search DIA-Umpire Q3-style groups with no MS1 feature (wide window).
    q3: bool,
    t2_corr: f32,
    t2_apex: i64,
    t2_psms: usize,
    /// Diagnose mode: TSV of target precursors.
    targets: Option<String>,
}

fn parse_args() -> Args {
    let mut args = Args {
        raw: String::new(),
        fasta: String::new(),
        out: String::new(),
        report_psms: 5,
        ms2_min_scans: 3,
        mode: "rescore".into(),
        min_corr: 0.5,
        apex_tolerance: 2,
        q3: false,
        t2_corr: 0.3,
        t2_apex: 2,
        t2_psms: 1,
        targets: None,
    };
    let mut it = std::env::args().skip(1);
    while let Some(flag) = it.next() {
        let mut value = || it.next().expect("missing flag value");
        match flag.as_str() {
            "--raw" => args.raw = value(),
            "--fasta" => args.fasta = value(),
            "--out" => args.out = value(),
            "--report-psms" => args.report_psms = value().parse().unwrap(),
            "--ms2-min-scans" => args.ms2_min_scans = value().parse().unwrap(),
            "--mode" => args.mode = value(),
            "--min-corr" => args.min_corr = value().parse().unwrap(),
            "--apex-tolerance" => args.apex_tolerance = value().parse().unwrap(),
            "--q3" => args.q3 = true,
            "--t2-corr" => args.t2_corr = value().parse().unwrap(),
            "--t2-apex" => args.t2_apex = value().parse().unwrap(),
            "--t2-psms" => args.t2_psms = value().parse().unwrap(),
            "--targets" => args.targets = Some(value()),
            other => panic!("unknown flag {other}"),
        }
    }
    assert!(!args.raw.is_empty() && !args.fasta.is_empty() && !args.out.is_empty());
    args
}

/// Peak resident memory so far, in MB.
fn peak_rss_mb() -> f64 {
    std::fs::read_to_string("/proc/self/status")
        .ok()
        .and_then(|s| {
            s.lines()
                .find(|l| l.starts_with("VmHWM:"))
                .and_then(|l| l.split_whitespace().nth(1))
                .and_then(|kb| kb.parse::<f64>().ok())
        })
        .unwrap_or(0.0)
        / 1024.0
}

fn window_key(w: (f64, f64, f64)) -> (i64, i64, i64) {
    sage_dia::koth_core::IsolationWindow {
        target: w.0,
        lower: w.1,
        upper: w.2,
    }
    .key()
}

/// Target-decoy q-values for `(score, is_decoy)`; returns targets at q <= 1%.
fn q_values(rows: &mut [(f64, bool, f32)]) -> usize {
    rows.sort_by(|a, b| b.0.total_cmp(&a.0));
    let (mut t, mut d) = (0usize, 0usize);
    for r in rows.iter_mut() {
        if r.1 {
            d += 1
        } else {
            t += 1
        }
        r.2 = (d as f32 + 1.0) / t.max(1) as f32;
    }
    let mut q = 1.0f32;
    let mut pass = 0;
    for r in rows.iter_mut().rev() {
        q = q.min(r.2);
        r.2 = q;
        if q <= 0.01 && !r.1 {
            pass += 1;
        }
    }
    pass
}

struct Psm {
    feature: Feature,
    tier: u8,
    peptide: String,
    coel: CoelutionFeatures,
}

fn hill_row(c: &CoelutionFeatures) -> [f64; 10] {
    [
        (c.n_coeluting as f64).ln_1p(),
        (c.n_coapex as f64).ln_1p(),
        c.n_with_hill as f64 / c.n_fragments.max(1) as f64,
        c.frag_corr as f64,
        (c.apex_spread as f64).ln_1p(),
        (c.apex_offset as f64).ln_1p(),
        c.apex_fraction as f64,
        c.ms1_present as u8 as f64,
        c.ms1_corr as f64,
        (c.ms1_apex_delta as f64 * 60.0).ln_1p(),
    ]
}

/// Counts at 1% FDR for a score: (PSMs, precursors, peptides).
fn fdr_counts(psms: &[Psm], score: &[f64]) -> (usize, usize, usize) {
    let mut rows: Vec<(f64, bool, f32)> = psms
        .iter()
        .zip(score)
        .map(|(p, &s)| (s, p.feature.label == -1, 1.0))
        .collect();
    let n_psm = q_values(&mut rows);
    let best = |key: &dyn Fn(&Psm) -> String| {
        let mut best: HashMap<String, (f64, bool)> = HashMap::new();
        for (p, &s) in psms.iter().zip(score) {
            let e = best
                .entry(key(p))
                .or_insert((f64::NEG_INFINITY, p.feature.label == -1));
            if s > e.0 {
                *e = (s, p.feature.label == -1);
            }
        }
        let mut rows: Vec<(f64, bool, f32)> =
            best.into_values().map(|(s, d)| (s, d, 1.0)).collect();
        q_values(&mut rows)
    };
    let n_prec = best(&|p| format!("{}/{}", p.peptide, p.feature.charge));
    let n_pep = best(&|p| p.peptide.clone());
    (n_psm, n_prec, n_pep)
}

/// Percolator-style semi-supervised LDA: positives are targets passing 1% under
/// the current score, negatives all decoys; three rounds.
fn semi_supervised<const D: usize>(
    psms: &[Psm],
    init: &[f64],
    row: impl Fn(&Psm) -> [f64; D] + Sync,
) -> Vec<f64> {
    let mut score = init.to_vec();
    for _ in 0..3 {
        let mut rows: Vec<(f64, bool, f32, usize)> = psms
            .iter()
            .zip(&score)
            .enumerate()
            .map(|(i, (p, &s))| (s, p.feature.label == -1, 1.0, i))
            .collect();
        let mut tmp: Vec<(f64, bool, f32)> = rows.iter().map(|r| (r.0, r.1, 1.0)).collect();
        q_values(&mut tmp);
        rows.sort_by(|a, b| b.0.total_cmp(&a.0));
        let train: Vec<usize> = rows
            .iter()
            .zip(&tmp)
            .filter(|(r, t)| r.1 || t.2 <= 0.01)
            .map(|(r, _)| r.3)
            .collect();
        let items: Vec<&Psm> = train.iter().map(|&i| &psms[i]).collect();
        let decoy: Vec<bool> = items.iter().map(|p| p.feature.label == -1).collect();
        let Ok(lda) =
            LinearDiscriminantAnalysis::train_regularized(&items, &decoy, |p| row(p), 1e-3)
        else {
            break;
        };
        score = psms.par_iter().map(|p| lda.score(&row(p))).collect();
    }
    score
}

fn median(mut v: Vec<f64>) -> f64 {
    if v.is_empty() {
        return f64::NAN;
    }
    v.sort_by(f64::total_cmp);
    v[v.len() / 2]
}

type WindowKey = (i64, i64, i64);

struct Hills {
    ms1: Channel,
    windows: Vec<Channel>,
    /// spec_id -> (window index, cycle in that window)
    spec_cycle: HashMap<String, (usize, usize)>,
    precursors: Vec<PrecursorTrace>,
}

/// koth-core hill detection: MS1 hills + isotope features, and per-window MS2 hills.
fn detect_hills(raw: &[RawSpectrum], args: &Args, lap: &dyn Fn(&str)) -> anyhow::Result<Hills> {
    let mut ms1_rts = Vec::new();
    let mut cycles: HashMap<WindowKey, Vec<f32>> = HashMap::new();
    let mut spec_key: Vec<(String, WindowKey, usize)> = Vec::new();
    for s in raw {
        match s.ms_level {
            1 => ms1_rts.push(s.scan_start_time),
            2 => {
                if let Some(w) = window_of(s) {
                    let rts = cycles.entry(window_key(w)).or_default();
                    spec_key.push((s.id.clone(), window_key(w), rts.len()));
                    rts.push(s.scan_start_time);
                }
            }
            _ => {}
        }
    }

    let kcfg = sage_dia::koth_core::KothConfig {
        hills_ms2: Some(sage_dia::koth_core::config::HillsMs2Overrides {
            min_scans: Some(args.ms2_min_scans),
            ..Default::default()
        }),
        ..Default::default()
    };
    let t = Instant::now();
    let ms1_hills = sage_dia::koth_core::hills::detect_hills_from_iter(
        raw.iter()
            .filter(|s| s.ms_level == 1)
            .enumerate()
            .map(|(i, s)| to_koth(s, i)),
        &kcfg.hills,
        &kcfg.file,
    );
    let t_ms1 = t.elapsed().as_secs_f64();
    let n_ms1_hills = ms1_hills.len();
    let t = Instant::now();
    let features = sage_dia::koth_core::run_features(&ms1_hills, &kcfg.features, &kcfg.file)?;
    let precursors: Vec<PrecursorTrace> = features
        .iter()
        .filter_map(PrecursorTrace::from_koth)
        .collect();
    drop(features);
    lap(&format!(
        "koth MS1: {n_ms1_hills} hills in {t_ms1:.1}s; {} charged isotope features in {:.1}s",
        precursors.len(),
        t.elapsed().as_secs_f64()
    ));
    let ms1 = Channel::from_hills(0.0, f64::INFINITY, ms1_rts, ms1_hills);

    let t = Instant::now();
    // One koth MS2 detector call per isolation window, windows in parallel
    // (koth's own multi-window call runs the windows serially).
    let mut by_window: HashMap<WindowKey, Vec<&RawSpectrum>> = HashMap::new();
    for s in raw.iter().filter(|s| s.ms_level == 2) {
        if let Some(w) = window_of(s) {
            by_window.entry(window_key(w)).or_default().push(s);
        }
    }
    let ms2_cfg = kcfg.ms2_hills();
    let ms2_hills: Vec<sage_dia::koth_core::Hill> = by_window
        .into_par_iter()
        .flat_map_iter(|(_, spectra)| {
            sage_dia::koth_core::hills::detect_ms2_hills_from_iter(
                spectra.into_iter().enumerate().map(|(i, s)| to_koth(s, i)),
                &ms2_cfg,
                &kcfg.file,
            )
        })
        .collect();
    let t_ms2 = t.elapsed().as_secs_f64();
    let n_ms2_hills = ms2_hills.len();
    let koth_bytes: usize = ms2_hills
        .iter()
        .map(|h| std::mem::size_of_val(h) + h.intensity_profile.len() * 4)
        .sum();
    let mut windows = Vec::new();
    let mut index: HashMap<WindowKey, usize> = HashMap::new();
    for (iw, hills) in sage_dia::koth_core::group_ms2_hills_by_window(ms2_hills) {
        let rts = cycles.remove(&iw.key()).unwrap_or_default();
        index.insert(iw.key(), windows.len());
        windows.push(Channel::from_hills(iw.lower, iw.upper, rts, hills));
    }
    let spec_cycle = spec_key
        .into_iter()
        .filter_map(|(id, key, cycle)| Some((id, (*index.get(&key)?, cycle))))
        .collect();
    let compact: usize = windows.iter().map(Channel::heap_bytes).sum::<usize>() + ms1.heap_bytes();
    let mut widths: Vec<f64> = windows.iter().map(|w| w.upper - w.lower).collect();
    widths.sort_by(f64::total_cmp);
    let cycle_s: Vec<f64> = ms1
        .rts
        .windows(2)
        .map(|w| (w[1] - w[0]) as f64 * 60.0)
        .collect();
    lap(&format!(
        "koth MS2: {n_ms2_hills} hills in {} windows ({:.1}-{:.1} m/z wide, MS1 cycle {:.2} s) in {t_ms2:.1}s; koth MS2 Hill structs {:.0} MB -> compact MS1+MS2 channels {:.0} MB",
        windows.len(),
        widths.first().unwrap_or(&0.0),
        widths.last().unwrap_or(&0.0),
        median(cycle_s),
        koth_bytes as f64 / 1e6,
        compact as f64 / 1e6
    ));
    Ok(Hills {
        ms1,
        windows,
        spec_cycle,
        precursors,
    })
}

fn build_db(args: &Args) -> anyhow::Result<IndexedDatabase> {
    let config = serde_json::json!({
        "database": {
            "bucket_size": 32768,
            "fasta": args.fasta,
            "decoy_tag": "rev_",
            "generate_decoys": true,
            "enzyme": {"missed_cleavages": 1, "min_len": 7, "max_len": 30},
            "peptide_min_mass": 600.0,
            "peptide_max_mass": 4000.0,
            "static_mods": {"C": 57.021464},
            "variable_mods": {"M": [15.9949]},
            "max_variable_mods": 1
        },
        "precursor_tol": {"ppm": [-10.0, 10.0]},
        "fragment_tol": {"ppm": [-15.0, 15.0]},
        "mzml_paths": [args.raw]
    });
    let config_path = format!("{}/config.json", args.out);
    std::fs::write(&config_path, serde_json::to_string_pretty(&config)?)?;
    let parameters: Parameters = sage_cli::input::Input::load(&config_path)?
        .build()?
        .database;
    let fasta = sage_cloudpath::util::read_fasta(
        &sage_cloudpath::to_url(&args.fasta)?,
        &parameters.decoy_tag,
        parameters.generate_decoys,
    )?;
    Ok(parameters.build(fasta))
}

/// Search, then Sage's LDA; returns candidates and search seconds.
fn search(
    db: &IndexedDatabase,
    spectra: &[ProcessedSpectrum],
    wide: bool,
    report_psms: usize,
) -> (Vec<Feature>, f64) {
    let tims = TIMS.load(Ordering::Relaxed);
    let (prec_ppm, frag_ppm) = if tims { (15.0, 20.0) } else { (10.0, 15.0) };
    let precursor_tol = if wide {
        Tolerance::Da(-2.0, 2.0)
    } else {
        Tolerance::Ppm(-prec_ppm, prec_ppm)
    };
    let scorer = Scorer {
        db,
        precursor_tol,
        fragment_tol: Tolerance::Ppm(-frag_ppm, frag_ppm),
        min_matched_peaks: 4,
        min_isotope_err: if wide { 0 } else { -1 },
        max_isotope_err: if wide { 0 } else { 1 },
        min_precursor_charge: 2,
        max_precursor_charge: 4,
        override_precursor_charge: wide,
        max_fragment_charge: Some(2),
        chimera: wide,
        report_psms,
        wide_window: wide,
        annotate_matches: !wide,
        mass_shift_ppm: 20.0,
        score_type: ScoreType::SageHyperScore,
        mass_recalibration: None,
    };
    let t = Instant::now();
    let mut features: Vec<Feature> = spectra.par_iter().flat_map(|q| scorer.score(q)).collect();
    let secs = t.elapsed().as_secs_f64();
    if sage_core::ml::linear_discriminant::score_psms(&mut features, precursor_tol).is_err() {
        eprintln!("WARNING: Sage LDA failed; falling back to hyperscore");
        let bad = |name: &str, f: &dyn Fn(&Feature) -> f64| {
            let n = features.iter().filter(|x| !f(x).is_finite()).count();
            let lo = features.iter().map(f).fold(f64::INFINITY, f64::min);
            let hi = features.iter().map(f).fold(f64::NEG_INFINITY, f64::max);
            eprintln!("  {name}: {n} non-finite, range {lo}..{hi}");
        };
        bad("hyperscore", &|x| x.hyperscore);
        bad("delta_next", &|x| x.delta_next);
        bad("delta_best", &|x| x.delta_best);
        bad("poisson", &|x| x.poisson);
        bad("aligned_delta_mass", &|x| x.aligned_delta_mass as f64);
        bad("aligned_average_ppm", &|x| x.aligned_average_ppm as f64);
        bad("matched_intensity_pct", &|x| x.matched_intensity_pct as f64);
        bad("aligned_rt", &|x| x.aligned_rt as f64);
        bad("ims", &|x| x.ims as f64);
        bad("delta_rt_model", &|x| x.delta_rt_model as f64);
        bad("isotope_error", &|x| x.isotope_error as f64);
        for f in &mut features {
            f.discriminant_score = f.hyperscore as f32;
        }
    }
    (features, secs)
}

fn spectrum_stats(spectra: &[ProcessedSpectrum]) -> String {
    let peaks: usize = spectra.iter().map(|s| s.masses.len()).sum();
    format!(
        "{} spectra searched, {:.1} peaks/spectrum after processing",
        spectra.len(),
        peaks as f64 / spectra.len().max(1) as f64
    )
}

fn to_psms(db: &IndexedDatabase, features: Vec<Feature>) -> Vec<Psm> {
    features
        .into_par_iter()
        .map(|feature| Psm {
            peptide: db.resolve_peptide(&feature).to_string(),
            feature,
            coel: CoelutionFeatures::default(),
            tier: 1,
        })
        .collect()
}

/// Peptides passing 1% peptide-level FDR under `score`.
fn passing_peptides(psms: &[Psm], score: &[f64]) -> std::collections::HashSet<String> {
    let mut best: HashMap<&str, (f64, bool)> = HashMap::new();
    for (p, &s) in psms.iter().zip(score) {
        let e = best
            .entry(p.peptide.as_str())
            .or_insert((f64::NEG_INFINITY, p.feature.label == -1));
        if s > e.0 {
            *e = (s, p.feature.label == -1);
        }
    }
    let mut rows: Vec<(f64, bool, f32, &str)> =
        best.into_iter().map(|(k, (s, d))| (s, d, 1.0, k)).collect();
    rows.sort_by(|a, b| b.0.total_cmp(&a.0));
    let mut q: Vec<(f64, bool, f32)> = rows.iter().map(|r| (r.0, r.1, 1.0)).collect();
    q_values(&mut q);
    rows.iter()
        .zip(&q)
        .filter(|(r, q)| !r.1 && q.2 <= 0.01)
        .map(|(r, _)| r.3.to_string())
        .collect()
}

/// Sage-like PSM features for a regularized LDA (Sage's `score_psms` is
/// unregularized and fails when columns are constant, as on pseudo-spectra).
fn sage_row(p: &Psm) -> [f64; 13] {
    let f = &p.feature;
    let poisson = match (-f.poisson).ln_1p() {
        x if x.is_finite() => x,
        _ => 3.5,
    };
    [
        f.hyperscore.ln_1p(),
        f.delta_next.ln_1p(),
        poisson,
        (f.delta_mass as f64).abs(),
        f.average_ppm as f64,
        (f.matched_intensity_pct as f64).ln_1p(),
        f.matched_peaks as f64,
        (f.longest_b as f64).ln_1p(),
        (f.longest_y as f64).ln_1p(),
        (f.peptide_len as f64).ln_1p(),
        f.missed_cleavages as f64,
        f.charge as f64,
        f.isotope_error as f64,
    ]
}

fn regularized_scores(psms: &[Psm]) -> Vec<f64> {
    let init: Vec<f64> = psms.iter().map(|p| p.feature.hyperscore).collect();
    semi_supervised::<13>(psms, &init, sage_row)
}

fn lda_scores(psms: &[Psm]) -> Vec<f64> {
    psms.iter()
        .map(|p| p.feature.discriminant_score as f64)
        .collect()
}

fn fdr_line(name: &str, psms: &[Psm], score: &[f64]) -> String {
    let (a, b, c) = fdr_counts(psms, score);
    format!("{name}\tPSMs@1% {a}\tprecursors@1% {b}\tpeptides@1% {c}\n")
}

/// Per-PSM q-values under `score`, aligned with `psms`.
fn psm_q_values(psms: &[Psm], score: &[f64]) -> Vec<f32> {
    let mut order: Vec<usize> = (0..psms.len()).collect();
    order.sort_by(|&a, &b| score[b].total_cmp(&score[a]));
    let mut rows: Vec<(f64, bool, f32)> = order
        .iter()
        .map(|&i| (score[i], psms[i].feature.label == -1, 1.0))
        .collect();
    q_values(&mut rows);
    let mut q = vec![1.0; psms.len()];
    for (&i, r) in order.iter().zip(&rows) {
        q[i] = r.2;
    }
    q
}

/// Tier-aware LDA row: Sage-like features plus a tier-2 indicator.
fn tier_row(p: &Psm) -> [f64; 14] {
    let r = sage_row(p);
    let mut out = [0.0; 14];
    out[..13].copy_from_slice(&r);
    out[13] = (p.tier == 2) as u8 as f64;
    out
}

/// Hills of the input: from the `.raw` scans, or straight from `.d` frames.
fn run_hills(args: &Args, raw: &[RawSpectrum]) -> anyhow::Result<dia::RunHills> {
    if TIMS.load(Ordering::Relaxed) {
        sage_dia::tims::detect_hills(std::path::Path::new(&args.raw), args.ms2_min_scans as u32)
    } else {
        dia::detect_hills(raw, args.ms2_min_scans as u32)
    }
}

fn tiered(args: &Args, raw: Vec<RawSpectrum>, lap: &dyn Fn(&str)) -> anyhow::Result<String> {
    let hills = run_hills(args, &raw)?;
    drop(raw);
    lap(&format!(
        "hills: {} MS1 precursor features, {} windows",
        hills.precursors.len(),
        hills.windows.len()
    ));
    let n_hills: usize = hills.windows.iter().map(Channel::len).sum();
    let t1 = PseudoSettings {
        min_corr: args.min_corr,
        apex_tolerance: args.apex_tolerance,
        ..Default::default()
    };
    let t = Instant::now();
    let tier1 = dia::build_pseudo(&hills, &t1);
    let t_build1 = t.elapsed().as_secs_f64();
    lap(&format!(
        "tier 1: {} pseudo-spectra from {} precursors, {n_hills} fragment hills",
        tier1.len(),
        hills.precursors.len()
    ));
    let processor = SpectrumProcessor::new(150, false, 0.0);
    let spectra: Vec<ProcessedSpectrum> = tier1
        .par_iter()
        .enumerate()
        .map(|(i, p)| processor.process(dia::to_raw(0, i, p, None)))
        .collect();
    let stats1 = spectrum_stats(&spectra);
    let db = build_db(args)?;
    let (features, secs1) = search(&db, &spectra, false, 1);
    drop(spectra);
    let mut psms1 = to_psms(&db, features);
    let score1 = regularized_scores(&psms1);
    let q1 = psm_q_values(&psms1, &score1);
    // Subtract fragment hills matched by confident tier-1 PSMs. Short ions
    // (b1/b2/y1/y2) are shared by many peptides, so they stay.
    let matches: Vec<(usize, Vec<f32>)> = psms1
        .iter()
        .zip(&q1)
        .filter(|(p, &q)| p.feature.label == 1 && p.feature.rank == 1 && q <= 0.01)
        .filter_map(|(p, _)| {
            let i: usize = p
                .feature
                .spec_id
                .strip_prefix("pseudo=")?
                .split(' ')
                .next()?
                .parse()
                .ok()?;
            let f = p.feature.fragments.as_ref()?;
            let mzs = f
                .fragment_ordinals
                .iter()
                .zip(&f.mz_experimental)
                .filter(|(&o, _)| o >= 3)
                .map(|(_, &mz)| mz)
                .collect();
            Some((i, mzs))
        })
        .collect();
    let used = dia::subtract(
        &hills,
        &tier1,
        matches.iter().map(|(i, m)| (*i, m.as_slice())),
        15.0,
    );
    let n_used: usize = used.iter().map(|u| u.iter().filter(|&&x| x).count()).sum();
    let t2 = PseudoSettings {
        min_corr: args.t2_corr,
        apex_tolerance: args.t2_apex,
        max_peaks: 100,
        ..Default::default()
    };
    let t = Instant::now();
    let tier2 = dia::build_tier2(&hills, &used, &t2);
    let t_build2 = t.elapsed().as_secs_f64();
    drop(hills);
    lap(&format!(
        "subtracted {n_used} hills using {} tier-1 PSMs; tier 2: {} groups",
        matches.len(),
        tier2.len()
    ));
    let spectra: Vec<ProcessedSpectrum> = tier2
        .par_iter()
        .enumerate()
        .map(|(i, (p, w))| processor.process(dia::to_raw(0, i, p, Some(*w))))
        .collect();
    let stats2 = spectrum_stats(&spectra);
    let (features, secs2) = search(&db, &spectra, true, args.t2_psms);
    drop(spectra);
    let mut psms2 = to_psms(&db, features);
    for p in &mut psms2 {
        p.tier = 2;
    }
    let score2 = regularized_scores(&psms2);
    lap("tier 2 searched");

    let mut report = format!(
        "tier 1: min_corr {} apex_tolerance {}; {stats1}; build {t_build1:.1}s, search {secs1:.1}s\n",
        args.min_corr, args.apex_tolerance
    );
    report += &fdr_line("tier 1 alone, regularized LDA", &psms1, &score1);
    report += &format!(
        "subtracted {n_used} of {n_hills} fragment hills ({} tier-1 target PSMs at 1%)\n",
        matches.len()
    );
    report += &format!(
        "tier 2: min_corr {} apex_tolerance {} max_peaks 100 report_psms {}; {stats2}; build {t_build2:.1}s, search {secs2:.1}s\n",
        args.t2_corr, args.t2_apex, args.t2_psms
    );
    report += &fdr_line("tier 2 alone, regularized LDA", &psms2, &score2);
    let lda2 = lda_scores(&psms2);
    report += &fdr_line("tier 2 alone, Sage LDA", &psms2, &lda2);
    let a = passing_peptides(&psms1, &score1);
    for (name, s2) in [("regularized", &score2), ("Sage", &lda2)] {
        let b = passing_peptides(&psms2, s2);
        let added = b.difference(&a).count();
        report += &format!(
            "separate FDR per tier (tier 2 {name} LDA): peptides tier1 {} + tier2 new {added} = {}\n",
            a.len(),
            a.len() + added
        );
    }
    psms1.append(&mut psms2);
    let init: Vec<f64> = psms1.iter().map(|p| p.feature.hyperscore).collect();
    let score = semi_supervised::<14>(&psms1, &init, tier_row);
    report += &fdr_line("tiers pooled, tier as LDA feature", &psms1, &score);
    Ok(report)
}

/// Gap diagnosis: for each target precursor (`--targets` TSV: key, charge,
/// precursor m/z, RT in minutes, comma-separated theoretical fragment m/z),
/// find its MS1 isotope feature and count the theoretical fragments present in
/// its pseudo-spectrum under the shipped settings and with each grouping
/// threshold relaxed. Writes `diagnose.tsv` in `--out`.
fn diagnose(args: &Args, raw: Vec<RawSpectrum>, lap: &dyn Fn(&str)) -> anyhow::Result<String> {
    let tims = TIMS.load(Ordering::Relaxed);
    let (prec_ppm, frag_ppm) = if tims {
        (15.0f32, 20.0f32)
    } else {
        (10.0, 15.0)
    };
    let hills = run_hills(args, &raw)?;
    drop(raw);
    lap(&format!(
        "hills: {} MS1 precursor features, {} windows",
        hills.precursors.len(),
        hills.windows.len()
    ));
    let ms1 = &hills.ms1;
    let precs = &hills.precursors;
    let mut by_mz: Vec<u32> = (0..precs.len() as u32).collect();
    by_mz.sort_by(|&a, &b| precs[a as usize].mz.total_cmp(&precs[b as usize].mz));
    let sorted_mz: Vec<f32> = by_mz.iter().map(|&i| precs[i as usize].mz).collect();
    let rt_span = |p: &PrecursorTrace| {
        let last = (p.start as usize + p.profile.len()).saturating_sub(1);
        (
            ms1.rts[p.start as usize],
            ms1.rts[last.min(ms1.rts.len() - 1)],
        )
    };
    // Features at `mz` (± prec_ppm) whose RT span covers `rt` (± 0.05 min).
    let near = |mz: f32, rt: f32| -> Vec<usize> {
        let tol = mz * prec_ppm * 1e-6;
        let lo = sorted_mz.partition_point(|&m| m < mz - tol);
        let hi = sorted_mz.partition_point(|&m| m <= mz + tol);
        by_mz[lo..hi]
            .iter()
            .map(|&i| i as usize)
            .filter(|&i| {
                let (a, b) = rt_span(&precs[i]);
                rt >= a - 0.05 && rt <= b + 0.05
            })
            .collect()
    };
    let strict = PseudoSettings::default();
    use pseudo::PeakRank;
    let variants: Vec<(&str, PseudoSettings)> = vec![
        ("strict", strict),
        (
            "corr0",
            PseudoSettings {
                min_corr: 0.0,
                ..strict
            },
        ),
        (
            "corr07",
            PseudoSettings {
                min_corr: 0.7,
                ..strict
            },
        ),
        (
            "apex4",
            PseudoSettings {
                apex_tolerance: 4,
                ..strict
            },
        ),
        (
            "im015",
            PseudoSettings {
                im_tolerance: 0.015,
                ..strict
            },
        ),
        (
            "im_any",
            PseudoSettings {
                im_tolerance: 100.0,
                ..strict
            },
        ),
        (
            "rk_corr",
            PseudoSettings {
                rank: PeakRank::Correlation,
                ..strict
            },
        ),
        (
            "rk_w2",
            PseudoSettings {
                rank: PeakRank::Weighted(2.0),
                ..strict
            },
        ),
        (
            "rk_w4",
            PseudoSettings {
                rank: PeakRank::Weighted(4.0),
                ..strict
            },
        ),
        (
            "cap300",
            PseudoSettings {
                max_peaks: 300,
                ..strict
            },
        ),
        (
            "cap1k",
            PseudoSettings {
                max_peaks: 1000,
                ..strict
            },
        ),
        (
            "all",
            PseudoSettings {
                min_corr: -2.0,
                apex_tolerance: 6,
                im_tolerance: 100.0,
                min_peaks: 0,
                max_peaks: 100_000,
                min_overlap: 1,
                ..strict
            },
        ),
    ];
    let count = |peaks: &[(f32, f32)], frags: &[f32]| -> usize {
        frags
            .iter()
            .filter(|&&f| {
                let tol = f * frag_ppm * 1e-6;
                let i = peaks.partition_point(|p| p.0 < f - tol);
                i < peaks.len() && peaks[i].0 <= f + tol
            })
            .count()
    };
    let text = std::fs::read_to_string(args.targets.as_deref().expect("--targets"))?;
    let targets: Vec<(String, u8, f32, f32, Vec<f32>)> = text
        .lines()
        .filter(|l| !l.is_empty())
        .map(|l| {
            let c: Vec<&str> = l.split('\t').collect();
            (
                c[0].to_string(),
                c[1].parse().unwrap(),
                c[2].parse().unwrap(),
                c[3].parse().unwrap(),
                c[4].split(',').map(|x| x.parse().unwrap()).collect(),
            )
        })
        .collect();
    let rows: Vec<String> = targets
        .par_iter()
        .map(|(key, z, mz, rt, frags)| {
            let zf = *z as f32;
            // Classify the best feature: exact, isotope offset k, other charge.
            let mut kind = "none".to_string();
            let mut pick: Option<usize> = None;
            let exact: Vec<usize> = near(*mz, *rt)
                .into_iter()
                .filter(|&i| precs[i].charge == *z)
                .collect();
            let closest = |v: &[usize]| {
                v.iter().copied().min_by(|&a, &b| {
                    let da = (ms1.rts[precs[a].apex as usize] - rt).abs();
                    let db = (ms1.rts[precs[b].apex as usize] - rt).abs();
                    da.total_cmp(&db)
                })
            };
            if let Some(i) = closest(&exact) {
                kind = "exact".into();
                pick = Some(i);
            } else {
                for k in [-1i32, 1, -2, 2] {
                    let m = mz + k as f32 * 1.003_355 / zf;
                    let v: Vec<usize> = near(m, *rt)
                        .into_iter()
                        .filter(|&i| precs[i].charge == *z)
                        .collect();
                    if let Some(i) = closest(&v) {
                        kind = format!("iso{k:+}");
                        pick = Some(i);
                        break;
                    }
                }
                if pick.is_none() {
                    let mut other: Vec<usize> = near(*mz, *rt);
                    for k in [-1i32, 1] {
                        other.extend(near(mz + k as f32 * 1.003_355 / zf, *rt));
                    }
                    if let Some(i) = closest(&other) {
                        kind = format!("charge{}", precs[i].charge);
                        pick = Some(i);
                    }
                }
            }
            // An MS1 hill at the precursor m/z covering this RT?
            let c1 = ms1.cycle_at(*rt) as u32;
            let hill = ms1
                .find(*mz, prec_ppm)
                .iter()
                .any(|h| h.start <= c1 + 2 && c1 <= h.end() + 2);
            let mut cols = vec![key.clone(), z.to_string(), kind, (hill as u8).to_string()];
            match pick {
                Some(i) => {
                    let p = &precs[i];
                    let boxes: Vec<usize> = hills
                        .windows
                        .iter()
                        .enumerate()
                        .filter(|(_, w)| w.contains(p.mz as f64, p.im as f64))
                        .map(|(k, _)| k)
                        .collect();
                    cols.push(format!("{:.2}", (p.mz - mz) / mz * 1e6));
                    cols.push(format!("{:.3}", ms1.rts[p.apex as usize] - rt));
                    cols.push(format!("{:.3}", p.im));
                    cols.push(boxes.len().to_string());
                    for (_, s) in &variants {
                        // Best over boxes: (matched theoretical, peaks); -1 if not built.
                        let best = boxes
                            .iter()
                            .filter_map(|&k| pseudo::build(p, ms1, &hills.windows[k], k, s))
                            .map(|ps| {
                                // Null: the same fragments shifted by +11 Th.
                                let null: Vec<f32> = frags.iter().map(|f| f + 11.0).collect();
                                (
                                    count(&ps.peaks, frags),
                                    ps.peaks.len(),
                                    count(&ps.peaks, &null),
                                )
                            })
                            .max();
                        match best {
                            Some((n, np, nn)) => cols.push(format!("{n}/{np}/{nn}")),
                            None => cols.push("-1/0/0".into()),
                        }
                    }
                }
                None => {
                    cols.extend(["", "", "", "0"].iter().map(|s| s.to_string()));
                    cols.extend(variants.iter().map(|_| "-1/0/0".to_string()));
                }
            }
            cols.push(frags.len().to_string());
            cols.join("\t")
        })
        .collect();
    let mut header = vec![
        "key", "charge", "kind", "ms1_hill", "dppm", "drt", "im", "boxes",
    ];
    header.extend(variants.iter().map(|v| v.0));
    header.push("n_frags");
    let mut out = header.join("\t") + "\n";
    for r in rows {
        out += &r;
        out.push('\n');
    }
    std::fs::write(format!("{}/diagnose.tsv", args.out), out)?;
    Ok(format!("diagnosed {} targets\n", targets.len()))
}

fn hill_filtered(
    args: &Args,
    mut raw: Vec<RawSpectrum>,
    lap: &dyn Fn(&str),
) -> anyhow::Result<String> {
    let hills = run_hills(args, &raw)?;
    let before: usize = raw
        .iter()
        .filter(|s| s.ms_level == 2)
        .map(|s| s.mz.len())
        .sum();
    let t = Instant::now();
    if TIMS.load(Ordering::Relaxed) {
        hill_filter_by_rt(&mut raw, &hills, 20.0);
    } else {
        dia::hill_filter(&mut raw, &hills, 15.0);
    }
    drop(hills);
    let after: usize = raw
        .iter()
        .filter(|s| s.ms_level == 2)
        .map(|s| s.mz.len())
        .sum();
    lap(&format!(
        "hill filter kept {after} of {before} MS2 peaks in {:.1}s",
        t.elapsed().as_secs_f64()
    ));
    let db = build_db(args)?;
    let processor = SpectrumProcessor::new(150, false, 0.0);
    let spectra: Vec<ProcessedSpectrum> = raw
        .into_par_iter()
        .filter(|s| s.ms_level == 2 && !s.precursors.is_empty() && !s.mz.is_empty())
        .map(|s| processor.process(s))
        .collect();
    let stats = spectrum_stats(&spectra);
    let (features, secs) = search(&db, &spectra, true, args.report_psms);
    drop(spectra);
    let psms = to_psms(&db, features);
    let mut report =
        format!("hill-filtered MS2: kept {after} of {before} peaks; {stats}; search {secs:.1}s\n");
    report += &fdr_line(
        "hill-filtered wide-window chimeric, Sage LDA",
        &psms,
        &lda_scores(&psms),
    );
    report += &fdr_line(
        "hill-filtered wide-window chimeric, regularized LDA",
        &psms,
        &regularized_scores(&psms),
    );
    Ok(report)
}

/// timsTOF variant of `dia::hill_filter`: each wide-window scan is one
/// quadrupole window of one frame, so its channel is the box with the same
/// m/z bounds and a cycle at the scan's retention time.
fn hill_filter_by_rt(raw: &mut [RawSpectrum], hills: &dia::RunHills, ppm: f32) {
    raw.par_iter_mut()
        .filter(|s| s.ms_level == 2)
        .for_each(|s| {
            let found = window_of(s).and_then(|(_, lo, hi)| {
                hills
                    .windows
                    .iter()
                    .filter(|w| (w.lower - lo).abs() < 0.5 && (w.upper - hi).abs() < 0.5)
                    .map(|w| {
                        let c = w.cycle_at(s.scan_start_time);
                        ((w.rts[c] - s.scan_start_time).abs(), w, c as u32)
                    })
                    .min_by(|a, b| a.0.total_cmp(&b.0))
            });
            let keep: Vec<bool> = match found {
                Some((dt, w, c)) if dt < 1e-3 => {
                    s.mz.iter()
                        .map(|&mz| w.find(mz, ppm).iter().any(|h| h.start <= c && c <= h.end()))
                        .collect()
                }
                _ => vec![false; s.mz.len()],
            };
            let mut it = keep.iter();
            s.mz.retain(|_| *it.next().unwrap());
            let mut it = keep.iter();
            s.intensity.retain(|_| *it.next().unwrap());
            s.total_ion_current = s.intensity.iter().sum();
        });
}

fn main() -> anyhow::Result<()> {
    let args = parse_args();
    let tims = sage_dia::tims::is_tdf(std::path::Path::new(&args.raw));
    TIMS.store(tims, Ordering::Relaxed);
    std::fs::create_dir_all(&args.out)?;
    let started = Instant::now();
    let lap = |label: &str| {
        eprintln!(
            "[{:7.1}s, peak RSS {:6.0} MB] {label}",
            started.elapsed().as_secs_f64(),
            peak_rss_mb()
        )
    };
    let url = sage_cloudpath::to_url(&args.raw)?;
    let mut raw = if !tims {
        sage_cloudpath::util::read_thermoraw(&url, 0, None).map_err(|e| anyhow::anyhow!("{e:?}"))?
    } else if matches!(args.mode.as_str(), "wide" | "hillfilter") {
        sage_cloudpath::util::read_spectra(&url, 0, None, Default::default(), false)
            .map_err(|e| anyhow::anyhow!("{e:?}"))?
    } else if matches!(args.mode.as_str(), "tiered" | "diagnose") {
        Vec::new()
    } else {
        anyhow::bail!("mode {} does not support timsTOF input", args.mode);
    };
    raw.sort_by(|a, b| a.scan_start_time.total_cmp(&b.scan_start_time));
    let n_ms1 = raw.iter().filter(|s| s.ms_level == 1).count();
    let n_ms2 = raw.iter().filter(|s| s.ms_level == 2).count();
    let peaks_ms2: usize = raw
        .iter()
        .filter(|s| s.ms_level == 2)
        .map(|s| s.mz.len())
        .sum();
    lap(&format!(
        "read {n_ms1} MS1 + {n_ms2} MS2 spectra ({:.1}M MS2 centroids, {:.0} raw peaks/MS2)",
        peaks_ms2 as f64 / 1e6,
        peaks_ms2 as f64 / n_ms2.max(1) as f64
    ));

    let mut report = format!("mode {}\n", args.mode);
    match args.mode.as_str() {
        "wide" => {
            let db = build_db(&args)?;
            lap(&format!("index built: {} peptides", db.peptides.len()));
            let processor = SpectrumProcessor::new(150, false, 0.0);
            let spectra: Vec<ProcessedSpectrum> = raw
                .into_par_iter()
                .filter(|s| s.ms_level == 2 && !s.precursors.is_empty())
                .map(|s| processor.process(s))
                .collect();
            let stats = spectrum_stats(&spectra);
            let (features, secs) = search(&db, &spectra, true, args.report_psms);
            drop(spectra);
            lap(&format!(
                "wide-window search: {} candidates in {secs:.1}s",
                features.len()
            ));
            let psms = to_psms(&db, features);
            report += &format!("{stats}; search {secs:.1}s\n");
            report += &fdr_line("wide-window chimeric, Sage LDA", &psms, &lda_scores(&psms));
            report += &fdr_line(
                "wide-window chimeric, regularized LDA",
                &psms,
                &regularized_scores(&psms),
            );
        }
        "pseudo" => {
            let hills = detect_hills(&raw, &args, &lap)?;
            drop(raw);
            let settings = PseudoSettings {
                min_corr: args.min_corr,
                apex_tolerance: args.apex_tolerance,
                ..Default::default()
            };
            let t = Instant::now();
            let hills_ref = &hills;
            let pseudo_spectra: Vec<pseudo::PseudoSpectrum> = hills
                .precursors
                .par_iter()
                .flat_map_iter(|p| {
                    let mz = p.mz as f64;
                    hills_ref
                        .windows
                        .iter()
                        .enumerate()
                        .filter(move |(_, w)| w.lower <= mz && mz <= w.upper)
                        .filter_map(move |(i, w)| pseudo::build(p, &hills_ref.ms1, w, i, &settings))
                })
                .collect();
            let t_build = t.elapsed().as_secs_f64();
            let n_prec = hills.precursors.len();
            // Q3 tier: groups of fragment hills not used by any MS1-anchored
            // pseudo-spectrum, searched with the isolation window as tolerance.
            let orphans: Vec<(pseudo::PseudoSpectrum, (f64, f64))> = if args.q3 {
                let mut used: Vec<Vec<bool>> =
                    hills.windows.iter().map(|w| vec![false; w.len()]).collect();
                for p in &pseudo_spectra {
                    for &h in &p.hills {
                        used[p.window][h as usize] = true;
                    }
                }
                hills
                    .windows
                    .par_iter()
                    .enumerate()
                    .flat_map_iter(|(i, w)| {
                        pseudo::build_orphans(w, i, &used[i], &settings)
                            .into_iter()
                            .map(move |p| (p, (w.lower, w.upper)))
                    })
                    .collect()
            } else {
                Vec::new()
            };
            drop(hills);
            lap(&format!(
                "{} pseudo-spectra from {n_prec} precursors in {t_build:.1}s; {} Q3 groups without an MS1 feature",
                pseudo_spectra.len(),
                orphans.len()
            ));
            let processor = SpectrumProcessor::new(150, false, 0.0);
            let spectra: Vec<ProcessedSpectrum> = pseudo_spectra
                .into_par_iter()
                .enumerate()
                .map(|(i, p)| processor.process(dia::to_raw(0, i, &p, None)))
                .collect();
            let stats = spectrum_stats(&spectra);
            let db = build_db(&args)?;
            lap(&format!("index built: {} peptides", db.peptides.len()));
            let (features, secs) = search(&db, &spectra, false, 1);
            drop(spectra);
            lap(&format!(
                "pseudo-spectrum search: {} candidates in {secs:.1}s",
                features.len()
            ));
            let psms = to_psms(&db, features);
            report += &format!(
                "min_corr {} apex_tolerance {} ms2_min_scans {}; {n_prec} precursors; {stats}; build {t_build:.1}s, search {secs:.1}s\n",
                args.min_corr, args.apex_tolerance, args.ms2_min_scans
            );
            report += &fdr_line(
                "pseudo-spectrum (MS1-anchored, closed), Sage LDA",
                &psms,
                &lda_scores(&psms),
            );
            let score = regularized_scores(&psms);
            report += &fdr_line(
                "pseudo-spectrum (MS1-anchored, closed), regularized LDA",
                &psms,
                &score,
            );
            if args.q3 {
                let q3: Vec<ProcessedSpectrum> = orphans
                    .into_par_iter()
                    .enumerate()
                    .map(|(i, (p, w))| processor.process(dia::to_raw(0, i, &p, Some(w))))
                    .collect();
                let q3_stats = spectrum_stats(&q3);
                let (features, q3_secs) = search(&db, &q3, true, 1);
                drop(q3);
                let q3_psms = to_psms(&db, features);
                let q3_score = regularized_scores(&q3_psms);
                report += &format!("Q3 tier: {q3_stats}; search {q3_secs:.1}s\n");
                report += &fdr_line("Q3 tier (no MS1 feature, wide-window)", &q3_psms, &q3_score);
                let main = passing_peptides(&psms, &score);
                let extra = passing_peptides(&q3_psms, &q3_score);
                let added = extra.difference(&main).count();
                report += &format!(
                    "peptides@1%: MS1-anchored {}, Q3 {}, Q3 adds {added} new (union {}, tiers FDR-controlled separately)\n",
                    main.len(),
                    extra.len(),
                    main.len() + added
                );
            }
        }
        "tiered" => report += &tiered(&args, raw, &lap)?,
        "diagnose" => report += &diagnose(&args, raw, &lap)?,
        "hillfilter" => report += &hill_filtered(&args, raw, &lap)?,
        "rescore" => {
            let hills = detect_hills(&raw, &args, &lap)?;
            let db = build_db(&args)?;
            lap(&format!("index built: {} peptides", db.peptides.len()));
            let processor = SpectrumProcessor::new(150, false, 0.0);
            let spectra: Vec<ProcessedSpectrum> = raw
                .into_par_iter()
                .filter(|s| s.ms_level == 2 && !s.precursors.is_empty())
                .map(|s| processor.process(s))
                .collect();
            let (features, secs) = search(&db, &spectra, true, args.report_psms);
            drop(spectra);
            lap(&format!(
                "wide-window search: {} candidates in {secs:.1}s",
                features.len()
            ));
            report += &rescore(&args, &db, &hills, features, &lap)?;
        }
        other => anyhow::bail!("unknown mode {other}"),
    }
    report += &format!(
        "wall {:.1}s, peak RSS {:.0} MB, rayon threads {}\n",
        started.elapsed().as_secs_f64(),
        peak_rss_mb(),
        rayon::current_num_threads()
    );
    print!("{report}");
    std::fs::write(format!("{}/report-{}.txt", args.out, args.mode), &report)?;
    Ok(())
}

fn rescore(
    args: &Args,
    db: &IndexedDatabase,
    hills: &Hills,
    features: Vec<Feature>,
    lap: &dyn Fn(&str),
) -> anyhow::Result<String> {
    let settings = CoelutionSettings::default();
    let t = Instant::now();
    let psms: Vec<Psm> = features
        .into_par_iter()
        .filter_map(|feature| {
            let &(w, cycle) = hills.spec_cycle.get(&feature.spec_id)?;
            let window = &hills.windows[w];
            let peptide = db.resolve_peptide(&feature);
            let max_z = feature.charge.saturating_sub(1).clamp(1, 2);
            let mut frags = Vec::new();
            for kind in [Kind::B, Kind::Y] {
                for ion in IonSeries::new(&peptide, kind) {
                    for z in 1..=max_z {
                        let mz = (ion.monoisotopic_mass + z as f32 * PROTON) / z as f32;
                        if (200.0..=1800.0).contains(&mz) {
                            frags.push(mz);
                        }
                    }
                }
            }
            let precursor_mz =
                (feature.calcmass + feature.charge as f32 * PROTON) / feature.charge as f32;
            let coel = coelution::score(&frags, window, cycle, precursor_mz, &hills.ms1, &settings);
            let peptide = peptide.to_string();
            Some(Psm {
                tier: 1,
                peptide,
                feature,
                coel,
            })
        })
        .collect();
    lap(&format!(
        "co-elution features for {} candidates in {:.1}s",
        psms.len(),
        t.elapsed().as_secs_f64()
    ));

    let mut out = String::new();
    // Target vs decoy medians of each co-elution feature, rank-1 candidates only.
    let names = [
        "ln n_coeluting",
        "ln n_coapex",
        "frac_with_hill",
        "frag_corr",
        "ln apex_spread",
        "ln apex_offset",
        "apex_fraction",
        "ms1_present",
        "ms1_corr",
        "ln ms1_apex_delta_s",
    ];
    let base = lda_scores(&psms);
    let mut q: Vec<(f64, bool, f32)> = psms
        .iter()
        .zip(&base)
        .map(|(p, &s)| (s, p.feature.label == -1, 1.0))
        .collect();
    q_values(&mut q);
    let mut sorted: Vec<(f64, usize)> = base.iter().copied().zip(0..).collect();
    sorted.sort_by(|a, b| b.0.total_cmp(&a.0));
    let mut qv = vec![1.0f32; psms.len()];
    for ((_, i), r) in sorted.iter().zip(&q) {
        qv[*i] = r.2;
    }
    out += "feature\ttarget q<=1%\tall targets\tdecoys\n";
    for (k, name) in names.iter().enumerate() {
        let pick = |f: &dyn Fn(usize) -> bool| -> f64 {
            median(
                (0..psms.len())
                    .filter(|&i| f(i))
                    .map(|i| hill_row(&psms[i].coel)[k])
                    .collect(),
            )
        };
        out += &format!(
            "{name}\t{:.3}\t{:.3}\t{:.3}\n",
            pick(&|i| psms[i].feature.label != -1 && qv[i] <= 0.01),
            pick(&|i| psms[i].feature.label != -1),
            pick(&|i| psms[i].feature.label == -1),
        );
    }

    out += &fdr_line("baseline (Sage LDA)", &psms, &base);
    let hills_only = semi_supervised::<10>(&psms, &base, |p| hill_row(&p.coel));
    out += &fdr_line("hills-only LDA", &psms, &hills_only);
    let combined = semi_supervised::<11>(&psms, &base, |p| {
        let h = hill_row(&p.coel);
        let mut r = [0.0; 11];
        r[0] = p.feature.discriminant_score as f64;
        r[1..].copy_from_slice(&h);
        r
    });
    out += &fdr_line("combined LDA (Sage score + hills)", &psms, &combined);

    let mut f = std::io::BufWriter::new(std::fs::File::create(format!(
        "{}/candidates.tsv",
        args.out
    ))?);
    writeln!(
        f,
        "spec_id\tpeptide\tcharge\tlabel\trank\tsage_score\tcombined\t{}",
        names.join("\t")
    )?;
    for (i, p) in psms.iter().enumerate() {
        let row = hill_row(&p.coel);
        writeln!(
            f,
            "{}\t{}\t{}\t{}\t{}\t{:.4}\t{:.4}\t{}",
            p.feature.spec_id,
            p.peptide,
            p.feature.charge,
            p.feature.label,
            p.feature.rank,
            base[i],
            combined[i],
            row.iter()
                .map(|v| format!("{v:.3}"))
                .collect::<Vec<_>>()
                .join("\t")
        )?;
    }
    Ok(out)
}
