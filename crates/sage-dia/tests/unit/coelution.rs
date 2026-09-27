use super::*;
use crate::hills::Channel;
use koth_core::Hill;
use std::sync::Arc;

fn hill(mz: f64, start: usize, apex: usize, profile: &[f32]) -> Hill {
    let max = profile.iter().copied().fold(0.0f32, f32::max);
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
        scan_end: start + profile.len().saturating_sub(1),
        n_scans: profile.len(),
        skipped_scans: 0,
        intensity_sum: profile.iter().map(|&x| x as f64).sum(),
        intensity_max: max as f64,
        hill_score: 1.0,
        intensity_profile: Arc::from(profile),
        isolation_window: None,
        faims_cv: None,
    }
}

/// Ten cycles, 0.1 min apart.
fn rts() -> Vec<f32> {
    (0..10).map(|c| c as f32 * 0.1).collect()
}

const A: [f32; 5] = [1., 4., 9., 4., 1.];
const B: [f32; 5] = [2., 8., 18., 8., 2.];
const C: [f32; 5] = [9., 4., 1., 4., 9.];

fn window() -> Channel {
    Channel::from_hills(
        490.0,
        510.0,
        rts(),
        vec![
            hill(500.0, 2, 2, &C),
            hill(300.0, 2, 4, &A),
            hill(400.0, 2, 4, &B),
            // Out of reach of cycle 4.
            hill(700.0, 8, 9, &[3., 5.]),
            // Hills without a profile are dropped.
            hill(800.0, 0, 0, &[]),
        ],
    )
}

fn ms1() -> Channel {
    Channel::from_hills(0.0, f64::INFINITY, rts(), vec![hill(600.0, 1, 3, &A)])
}

#[test]
fn pearson_is_scale_free_and_zero_for_constant_input() {
    assert!((pearson(&A, &B) - 1.0).abs() < 1e-6);
    // sum(A*C) - n*mean(A)*mean(C) = 59 - 5*3.8*5.4, over sqrt(saa * scc).
    let expected = (59.0 - 5.0 * 3.8 * 5.4) / ((115.0f32 - 72.2) * (195.0 - 145.8)).sqrt();
    assert!((pearson(&A, &C) - expected).abs() < 1e-5);
    assert_eq!(pearson(&[2., 2., 2.], &A[..3]), 0.0);
}

#[test]
fn channel_lookups_follow_mz_apex_and_time() {
    let channel = window();
    assert_eq!(channel.len(), 4);
    assert!(!channel.is_empty());
    assert!(channel.heap_bytes() > 0);
    // Hills are stored in m/z order.
    let mzs: Vec<f32> = channel.hills().iter().map(|h| h.mz).collect();
    assert_eq!(mzs, vec![300.0, 400.0, 500.0, 700.0]);

    assert_eq!(channel.find(300.002, 15.0).len(), 1);
    assert!(channel.find(300.01, 15.0).is_empty());

    let by_apex: Vec<f32> = channel.apex_between(3, 9).map(|h| h.mz).collect();
    assert_eq!(by_apex.len(), 3);
    assert!(by_apex.contains(&300.0) && by_apex.contains(&400.0) && by_apex.contains(&700.0));
    let indexed: Vec<usize> = channel.apex_between_indexed(0, 2).map(|(i, _)| i).collect();
    assert_eq!(indexed, vec![2]);

    let a = &channel.find(300.0, 1.0)[0];
    assert_eq!(a.end(), 6);
    assert_eq!(channel.profile(a), &A);
    assert_eq!(channel.intensity(a, 4), 9.0);
    assert_eq!(channel.intensity(a, 1), 0.0);
    assert_eq!(channel.intensity(a, 7), 0.0);

    assert_eq!(channel.cycle_at(-1.0), 0);
    assert_eq!(channel.cycle_at(0.34), 3);
    assert_eq!(channel.cycle_at(0.36), 4);
    assert_eq!(channel.cycle_at(5.0), 9);

    assert!(channel.contains(500.0, 0.0));
    assert!(!channel.contains(511.0, 0.0));
    let boxed = window().with_im(0.8, 1.0);
    assert!(boxed.contains(500.0, 0.9));
    assert!(!boxed.contains(500.0, 1.1));
    assert!(boxed.contains(500.0, 0.0));
}

#[test]
fn score_measures_fragment_and_precursor_coelution() {
    let settings = CoelutionSettings {
        min_corr: 0.6,
        ..Default::default()
    };
    // 300.002 matches the same hill as 300.0; 700.0 is out of reach; 999 has no hill.
    let fragments = [300.0, 300.002, 400.0, 500.0, 700.0, 999.0];
    let out = score(&fragments, &window(), 4, 600.0, &ms1(), &settings);

    assert_eq!(out.n_fragments, 6);
    assert_eq!(out.n_with_hill, 3);

    // Profiles over cycles 0..=9: B = 2A, and C is anti-shaped.
    let slice = |profile: &[f32; 5]| {
        let mut v = [0.0f32; 10];
        v[2..7].copy_from_slice(profile);
        v
    };
    let (a, c) = (slice(&A), slice(&C));
    let r_ac = pearson(&a, &c);
    assert!((out.frag_corr - (1.0 + 2.0 * r_ac) / 3.0).abs() < 1e-5);
    assert!((r_ac - 0.078_45).abs() < 1e-4);

    // Leave-one-out: A vs 2A + C = 0.861 and B vs A + C = 0.663 clear 0.6;
    // C vs 3A does not.
    assert_eq!(out.n_coeluting, 2);

    // Top-k apexes [2, 4, 4]: median 4.
    assert!((out.apex_spread - 2.0 / 3.0).abs() < 1e-6);
    assert_eq!(out.apex_offset, 0.0);
    assert_eq!(out.n_coapex, 3);
    assert_eq!(out.apex_fraction, 1.0);

    // The precursor hill is A shifted one cycle earlier.
    assert!(out.ms1_present);
    let precursor: Vec<f32> = (0..10)
        .map(|i| if (1..6).contains(&i) { A[i - 1] } else { 0.0 })
        .collect();
    let reference: Vec<f32> = (0..10).map(|i| 3.0 * a[i] + c[i]).collect();
    assert!((out.ms1_corr - pearson(&precursor, &reference)).abs() < 1e-6);
    assert!((out.ms1_corr - 0.6093).abs() < 1e-3);
    assert!((out.ms1_apex_delta - 0.1).abs() < 1e-6);
}

#[test]
fn score_without_matching_hills_reports_only_the_fragment_count() {
    let out = score(
        &[123.4, 700.0],
        &window(),
        4,
        600.0,
        &ms1(),
        &CoelutionSettings::default(),
    );
    assert_eq!(out.n_fragments, 2);
    assert_eq!(out.n_with_hill, 0);
    assert_eq!(out.n_coeluting, 0);
    assert!(!out.ms1_present);
    assert_eq!(out.frag_corr, 0.0);
}

#[test]
fn score_with_one_fragment_and_no_precursor_hill() {
    let out = score(
        &[300.0],
        &window(),
        4,
        650.0,
        &ms1(),
        &CoelutionSettings::default(),
    );
    assert_eq!(out.n_with_hill, 1);
    // A single fragment has no pairs and no leave-one-out reference.
    assert_eq!(out.frag_corr, 0.0);
    assert_eq!(out.n_coeluting, 0);
    assert_eq!(out.n_coapex, 1);
    assert!(!out.ms1_present);
    assert_eq!(out.ms1_corr, 0.0);
}
