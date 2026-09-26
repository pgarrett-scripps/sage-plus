use super::*;
use crate::hills::Channel;
use koth_core::Hill;
use std::sync::Arc;

fn hill(mz: f64, start: usize, apex: usize, profile: &[f32]) -> Hill {
    Hill {
        hill_id: 0,
        mz,
        mz_std: 0.0,
        mz_se: 0.0,
        rt: 0.0,
        rt_start: 0.0,
        rt_end: 0.0,
        rt_width: 0.0,
        im: 0.0,
        im_std: 0.0,
        scan_start: start,
        scan_apex: apex,
        scan_end: start + profile.len() - 1,
        n_scans: profile.len(),
        skipped_scans: 0,
        intensity_sum: profile.iter().map(|&x| x as f64).sum(),
        intensity_max: profile.iter().copied().fold(0.0f32, f32::max) as f64,
        hill_score: 1.0,
        intensity_profile: Arc::from(profile),
        isolation_window: None,
        faims_cv: None,
    }
}

fn with_im(mut hill: Hill, im: f64) -> Hill {
    hill.im = im;
    hill
}

/// Ten cycles, 0.1 min apart, shared by MS1 and the window.
fn rts() -> Vec<f32> {
    (0..10).map(|c| c as f32 * 0.1).collect()
}

const A: [f32; 5] = [1., 4., 9., 4., 1.];
const B: [f32; 5] = [2., 8., 18., 8., 2.];
const C: [f32; 5] = [9., 4., 1., 4., 9.];

/// Hills in m/z order: 300 (0), 320 (1), 350 (2), 400 (3), 450 (4), 1200 (5).
fn window() -> Channel {
    Channel::from_hills(
        490.0,
        510.0,
        rts(),
        vec![
            hill(300.0, 2, 4, &A),
            hill(400.0, 2, 4, &B),
            // Anti-correlated with the precursor.
            hill(450.0, 2, 2, &C),
            // Heavier than a 2+ precursor at 500 m/z can produce.
            hill(1200.0, 2, 4, &A),
            // Apex outside the tolerance.
            hill(350.0, 5, 7, &A),
            // Too short to overlap three cycles.
            hill(320.0, 3, 4, &[5., 6.]),
        ],
    )
}

fn ms1() -> Channel {
    Channel::from_hills(0.0, f64::INFINITY, rts(), Vec::new())
}

fn precursor(im: f32) -> PrecursorTrace {
    PrecursorTrace {
        mz: 500.0,
        charge: 2,
        intensity: 9.0,
        im,
        start: 2,
        apex: 4,
        profile: A.to_vec(),
    }
}

fn settings(min_peaks: usize, max_peaks: usize) -> PseudoSettings {
    PseudoSettings {
        min_peaks,
        max_peaks,
        ..Default::default()
    }
}

#[test]
fn ion_mobility_matches_within_tolerance_or_when_missing() {
    assert!(im_match(0.0, 1.0, 0.03));
    assert!(im_match(1.0, 0.0, 0.03));
    assert!(im_match(1.0, 1.02, 0.03));
    assert!(!im_match(1.0, 1.05, 0.03));
}

#[test]
fn precursor_profile_is_zero_outside_its_cycles() {
    let trace = precursor(0.0);
    assert_eq!(trace.at(1), 0.0);
    assert_eq!(trace.at(4), 9.0);
    assert_eq!(trace.at(7), 0.0);
}

#[test]
fn pseudo_spectrum_keeps_coeluting_fragments_below_the_precursor_mass() {
    let spectrum = build(&precursor(0.0), &ms1(), &window(), 7, &settings(2, 150)).unwrap();
    assert_eq!(spectrum.peaks, vec![(300.0, 9.0), (400.0, 18.0)]);
    assert_eq!(spectrum.hills, vec![0, 3]);
    assert_eq!(spectrum.precursor_mz, 500.0);
    assert_eq!(spectrum.charge, 2);
    assert_eq!(spectrum.window, 7);
    assert_eq!(spectrum.precursor_intensity, 9.0);
    assert!((spectrum.rt - 0.4).abs() < 1e-6);

    // max_peaks keeps the most intense peaks.
    let capped = build(&precursor(0.0), &ms1(), &window(), 7, &settings(1, 1)).unwrap();
    assert_eq!(capped.peaks, vec![(400.0, 18.0)]);
    assert_eq!(capped.hills, vec![3]);

    // Too few peaks yields no spectrum.
    assert!(build(&precursor(0.0), &ms1(), &window(), 7, &settings(3, 150)).is_none());
}

#[test]
fn pseudo_spectrum_requires_matching_ion_mobility() {
    let window = Channel::from_hills(
        490.0,
        510.0,
        rts(),
        vec![
            with_im(hill(300.0, 2, 4, &A), 1.0),
            with_im(hill(310.0, 2, 4, &A), 1.2),
            hill(400.0, 2, 4, &A),
        ],
    );
    let spectrum = build(&precursor(1.0), &ms1(), &window, 0, &settings(1, 150)).unwrap();
    let mzs: Vec<f32> = spectrum.peaks.iter().map(|p| p.0).collect();
    assert_eq!(mzs, vec![300.0, 400.0]);
    assert_eq!(spectrum.im, 1.0);
}

#[test]
fn orphan_groups_seed_on_the_most_intense_unused_hill() {
    let window = window();
    let spectra = build_orphans(&window, 3, &[false; 6], &settings(2, 150));
    // Seed 400 gathers 300 and 1200 (same shape); 450 is anti-correlated,
    // 320 overlaps too few cycles and 350 peaks too late. No other seed
    // finds a partner.
    assert_eq!(spectra.len(), 1);
    let orphan = &spectra[0];
    assert_eq!(orphan.hills, vec![0, 3, 5]);
    assert_eq!(
        orphan.peaks,
        vec![(300.0, 9.0), (400.0, 18.0), (1200.0, 9.0)]
    );
    assert_eq!(orphan.precursor_mz, 500.0);
    assert_eq!(orphan.charge, 0);
    assert_eq!(orphan.precursor_intensity, 0.0);
    assert_eq!(orphan.window, 3);
    assert!((orphan.rt - 0.4).abs() < 1e-6);

    // Hills used by MS1-anchored spectra are neither seeds nor members.
    let mut used = [false; 6];
    used[3] = true;
    let spectra = build_orphans(&window, 3, &used, &settings(2, 150));
    assert_eq!(spectra.len(), 1);
    assert_eq!(spectra[0].hills, vec![0, 5]);

    assert!(build_orphans(&window, 3, &[true; 6], &settings(1, 150)).is_empty());
}
