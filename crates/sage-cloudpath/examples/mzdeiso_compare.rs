//! SCOPING PROTOTYPE: compare Sage's deisotope() against mzdeisotope on real MS2.
//! cargo run --release --features mzdeisotope --example mzdeiso_compare -- <file.mzML> [max_ms2]
use sage_cloudpath::deconv;
use sage_core::spectrum::{deisotope, DeisotopeSettings};
use std::collections::BTreeMap;
use std::time::Instant;

#[derive(Default, Debug)]
struct Acc {
    n: usize,
    input_peaks: usize,
    sage_kept: usize,
    sage_charged: usize,
    sage_z2plus: usize,
    sage_ns: u128,
    mzd_kept: usize,
    mzd_z2plus: usize,
    mzd_ns: u128,
    mzd_matched_sage: usize,
    mzd_matched_sage_same_z: usize,
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let path = &args[1];
    let limit: usize = args.get(2).and_then(|s| s.parse().ok()).unwrap_or(usize::MAX);
    let url = url::Url::from_file_path(std::fs::canonicalize(path).unwrap()).unwrap();
    let t = Instant::now();
    let mut spectra = sage_cloudpath::util::read_mzml(&url, 0, None).unwrap();
    eprintln!("read {} spectra in {:?}", spectra.len(), t.elapsed());

    let mut engine = deconv::ms2_engine();
    let settings = DeisotopeSettings::default();
    let mut by_group: BTreeMap<String, Acc> = BTreeMap::new();
    let mut zhist_sage: BTreeMap<u8, usize> = BTreeMap::new();
    let mut zhist_mzd: BTreeMap<u8, usize> = BTreeMap::new();
    let mut seen = 0usize;
    for s in spectra.iter().filter(|s| s.ms_level == 2) {
        if seen >= limit {
            break;
        }
        seen += 1;
        let max_z = s.precursors.first().and_then(|p| p.charge).unwrap_or(3);
        let a = by_group.entry(s.acquisition.label()).or_default();
        a.n += 1;
        a.input_peaks += s.mz.len();

        let t0 = Instant::now();
        let sage = deisotope(&s.mz, &s.intensity, None, max_z, settings, 0.0);
        a.sage_ns += t0.elapsed().as_nanos();
        let mut sage_roots: Vec<(f32, Option<u8>)> = sage
            .iter()
            .filter(|p| p.envelope.is_none())
            .map(|p| {
                let z = p.charge.unwrap_or(1);
                ((p.mz - 1.007_276_4) * f32::from(z), p.charge)
            })
            .collect();
        sage_roots.sort_by(|x, y| x.0.total_cmp(&y.0));
        a.sage_kept += sage_roots.len();
        for r in &sage_roots {
            if let Some(z) = r.1 {
                a.sage_charged += 1;
                if z >= 2 {
                    a.sage_z2plus += 1;
                }
                *zhist_sage.entry(z).or_default() += 1;
            }
        }

        let t1 = Instant::now();
        let mzd = deconv::deconvolve_ms2(&mut engine, s, i32::from(max_z));
        a.mzd_ns += t1.elapsed().as_nanos();
        a.mzd_kept += mzd.len();
        for &(mz, _, z) in &mzd {
            *zhist_mzd.entry(z).or_default() += 1;
            if z >= 2 {
                a.mzd_z2plus += 1;
            }
            let mass = (mz - 1.007_276_4) * f32::from(z.max(1));
            let tol = mass * 20e-6;
            let lo = sage_roots.partition_point(|r| r.0 < mass - tol);
            let hits = sage_roots[lo..]
                .iter()
                .take_while(|r| r.0 <= mass + tol)
                .collect::<Vec<_>>();
            if !hits.is_empty() {
                a.mzd_matched_sage += 1;
                if hits.iter().any(|r| r.1.unwrap_or(1) == z) {
                    a.mzd_matched_sage_same_z += 1;
                }
            }
        }
    }

    for (label, a) in &by_group {
        let n = a.n as f64;
        println!(
            "{label}: n={} input/spec={:.1} | sage kept/spec={:.1} charged={:.1} z>=2={:.1} us/spec={:.1} | mzd kept/spec={:.1} z>=2={:.1} us/spec={:.1} | mzd peaks matching a sage root (20ppm neutral mass)={:.1}% same z={:.1}%",
            a.n,
            a.input_peaks as f64 / n,
            a.sage_kept as f64 / n,
            a.sage_charged as f64 / n,
            a.sage_z2plus as f64 / n,
            a.sage_ns as f64 / n / 1e3,
            a.mzd_kept as f64 / n,
            a.mzd_z2plus as f64 / n,
            a.mzd_ns as f64 / n / 1e3,
            100.0 * a.mzd_matched_sage as f64 / a.mzd_kept.max(1) as f64,
            100.0 * a.mzd_matched_sage_same_z as f64 / a.mzd_kept.max(1) as f64,
        );
    }
    println!("sage charge hist (assigned roots): {zhist_sage:?}");
    println!("mzd charge hist: {zhist_mzd:?}");

    let t = Instant::now();
    let stats = deconv::refit_precursors(&mut spectra, false);
    println!("precursor refit (MS1, parallel) {:?}: {stats:?}", t.elapsed());
}
