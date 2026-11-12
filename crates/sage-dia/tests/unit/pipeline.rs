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
    let s: DiaSettings = serde_json::from_str(r#"{"mode": "pseudo", "min_corr": 0.6}"#).unwrap();
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
