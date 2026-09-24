use super::*;

/// Deterministic noise in roughly [-1, 1].
fn noise(i: u64) -> f32 {
    let mut x = i.wrapping_mul(0x9e3779b97f4a7c15);
    x ^= x >> 31;
    x = x.wrapping_mul(0xbf58476d1ce4e5b9);
    x ^= x >> 29;
    ((x % 20001) as f32 / 10000.0) - 1.0
}

fn points(n: u64, truth: impl Fn(f32, f32) -> f32, noise_ppm: f32) -> Vec<MassErrorPoint> {
    (0..n)
        .map(|i| {
            let rt = 5.0 + 110.0 * (i as f32 / n as f32);
            let mz = 400.0 + 800.0 * ((noise(i + 7_000_000) + 1.0) / 2.0);
            MassErrorPoint {
                rt_minutes: rt,
                mz,
                error_ppm: truth(rt, mz) + noise_ppm * noise(i),
                group: i,
            }
        })
        .collect()
}

#[test]
fn centered_noise_selects_no_model() {
    let data = points(3000, |_, _| 0.0, 2.0);
    let selection = select_model(&data, RecalibrationOptions::default());
    assert_eq!(selection.kind(), MassModelKind::None);
    assert!(selection.model.is_none());
    assert!(selection.validation_psms > 600);
}

#[test]
fn constant_offset_selects_static() {
    let data = points(3000, |_, _| 3.0, 1.0);
    let selection = select_model(&data, RecalibrationOptions::default());
    let model = selection.model.expect("static model");
    assert_eq!(model.kind, MassModelKind::Static);
    assert!((model.intercept_ppm - 3.0).abs() < 0.1);
}

#[test]
fn linear_drift_selects_rt_linear() {
    let data = points(3000, |rt, _| 1.0 + 0.05 * (rt - 60.0), 0.5);
    let selection = select_model(&data, RecalibrationOptions::default());
    let model = selection.model.expect("linear model");
    assert_eq!(model.kind, MassModelKind::Linear);
    assert_eq!(model.axes, MassModelAxes::Rt);
    let slope = model.rt.as_ref().unwrap().slope;
    assert!((slope - 0.05).abs() < 0.005, "slope {slope}");
}

#[test]
fn curved_drift_selects_smooth_and_linear_cap_holds() {
    let truth = |rt: f32, _| 4.0 * ((rt - 60.0) / 55.0).powi(2) - 1.0;
    let data = points(4000, truth, 0.4);
    let selection = select_model(&data, RecalibrationOptions::default());
    let model = selection.model.clone().expect("smooth model");
    assert_eq!(model.kind, MassModelKind::Smooth);
    assert!(model.rt.as_ref().unwrap().knots.len() <= 6);
    for rt in [10.0, 40.0, 60.0, 90.0, 110.0] {
        let err = (model.predict_ppm(rt, 700.0) - truth(rt, 700.0)).abs();
        assert!(err < 0.5, "rt {rt}: {err}");
    }

    let linear_only = select_model(
        &data,
        RecalibrationOptions {
            max_kind: MassModelKind::Linear,
            ..Default::default()
        },
    );
    assert!(linear_only.kind() <= MassModelKind::Linear);
}

#[test]
fn thin_data_falls_back() {
    let data = points(150, |rt, _| 0.05 * (rt - 60.0) + 2.0, 0.5);
    let selection = select_model(&data, RecalibrationOptions::default());
    // Too few PSMs for linear or smooth fits.
    assert!(selection.kind() <= MassModelKind::Static);
    let data = points(20, |_, _| 5.0, 0.5);
    let selection = select_model(&data, RecalibrationOptions::default());
    assert_eq!(selection.skipped.as_deref(), Some("too_few_psms"));
}

#[test]
fn predictions_are_clamped() {
    let data = points(3000, |rt, _| 0.2 * (rt - 60.0), 0.2);
    let selection = select_model(
        &data,
        RecalibrationOptions {
            max_abs_ppm: 5.0,
            ..Default::default()
        },
    );
    let model = selection.model.unwrap();
    assert!(model.predict_ppm(500.0, 700.0) <= 5.0);
    assert!(model.predict_ppm(-500.0, 700.0) >= -5.0);
    // Clamped inputs: beyond the fitted range the prediction is flat.
    let edge = model.predict_ppm(10_000.0, 700.0);
    assert_eq!(edge, model.predict_ppm(20_000.0, 700.0));
}

#[test]
fn fragment_groups_split_together() {
    // Many points per PSM; the split must be by group.
    let mut data = Vec::new();
    for group in 0..400u64 {
        for ion in 0..10u64 {
            data.push(MassErrorPoint {
                rt_minutes: 1.0 + group as f32 * 0.2,
                mz: 200.0 + ion as f32 * 100.0,
                error_ppm: -2.0 + 0.3 * noise(group * 100 + ion),
                group,
            });
        }
    }
    let selection = select_model(&data, RecalibrationOptions::default());
    assert_eq!(selection.psms, 400);
    assert_eq!(selection.fit_psms + selection.validation_psms, 400);
    assert_eq!(selection.kind(), MassModelKind::Static);
}

#[test]
fn correction_removes_predicted_error() {
    let model = MassErrorModel {
        kind: MassModelKind::Static,
        axes: MassModelAxes::None,
        intercept_ppm: 10.0,
        rt: None,
        mz: None,
        max_abs_ppm: 20.0,
    };
    let theoretical = 1000.0f32;
    let observed = theoretical * (1.0 + 10e-6);
    assert!((model.correct_mz(observed, 0.0) - theoretical).abs() < 1e-3);
    assert_ne!(stable_hash("a"), stable_hash("b"));
}
