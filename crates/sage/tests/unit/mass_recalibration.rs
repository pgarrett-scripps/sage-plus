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

#[test]
fn acquisition_groups_with_opposite_biases_get_separate_models() {
    use crate::spectrum::{AcquisitionGroup, Activation, MassAnalyzer};
    let hcd = AcquisitionGroup {
        analyzer: MassAnalyzer::Orbitrap,
        activation: Activation::Hcd,
    };
    let tof = AcquisitionGroup {
        analyzer: MassAnalyzer::Tof,
        activation: Activation::Hcd,
    };
    let trap = AcquisitionGroup {
        analyzer: MassAnalyzer::IonTrap,
        activation: Activation::Cid,
    };
    let thin = AcquisitionGroup {
        analyzer: MassAnalyzer::Orbitrap,
        activation: Activation::Etd,
    };
    // Interleaved PSMs from two analyzers with +4 and -4 ppm biases. Pooled,
    // no single offset fits both, so each group must be modelled separately.
    let mut data = Vec::new();
    for (i, point) in points(3000, |_, _| 0.0, 0.5).into_iter().enumerate() {
        let (group, bias) = if i % 2 == 0 { (hcd, 4.0) } else { (tof, -4.0) };
        let point = MassErrorPoint {
            error_ppm: point.error_ppm + bias,
            ..point
        };
        data.push((group, point));
    }
    for point in points(3000, |_, _| 30.0, 50.0).into_iter().take(600) {
        data.push((
            trap,
            MassErrorPoint {
                group: point.group + 10_000,
                ..point
            },
        ));
    }
    for point in points(20, |_, _| 5.0, 0.5) {
        data.push((
            thin,
            MassErrorPoint {
                group: point.group + 20_000,
                ..point
            },
        ));
    }

    let groups = [(tof, 1500), (hcd, 1500), (trap, 600), (thin, 20)];
    let selected = select_group_models(&data, &groups, RecalibrationOptions::default());
    assert_eq!(selected.len(), 4);
    let find = |group| selected.iter().find(|s| s.group == group).unwrap();
    let hcd_model = find(hcd).selection.model.clone().expect("orbitrap model");
    let tof_model = find(tof).selection.model.clone().expect("tof model");
    assert!(
        (hcd_model.intercept_ppm - 4.0).abs() < 0.2,
        "{}",
        hcd_model.intercept_ppm
    );
    assert!(
        (tof_model.intercept_ppm + 4.0).abs() < 0.2,
        "{}",
        tof_model.intercept_ppm
    );
    assert_eq!(
        find(trap).selection.skipped.as_deref(),
        Some("low_accuracy_analyzer")
    );
    assert!(find(trap).selection.model.is_none());
    assert_eq!(
        find(thin).selection.skipped.as_deref(),
        Some("too_few_psms")
    );

    // Corrections are looked up by group; unknown groups are left alone.
    let correction = FileMassCorrection {
        precursor: None,
        fragment: selected
            .iter()
            .map(|s| GroupMassCorrection {
                group: s.group,
                model: s.selection.model.clone(),
                fragment_tol: None,
            })
            .collect(),
        precursor_tol: None,
    };
    assert!((correction.fragment_ppm(hcd, 60.0, 700.0) - 4.0).abs() < 0.2);
    assert!((correction.fragment_ppm(tof, 60.0, 700.0) + 4.0).abs() < 0.2);
    assert_eq!(correction.fragment_ppm(trap, 60.0, 700.0), 0.0);
    assert_eq!(
        correction.fragment_ppm(AcquisitionGroup::default(), 60.0, 700.0),
        0.0
    );
}

#[test]
fn thermo_filters_map_to_acquisition_groups() {
    use crate::spectrum::{AcquisitionGroup, Activation, MassAnalyzer};
    let parse = AcquisitionGroup::from_thermo_filter;
    let group = parse("FTMS + p NSI d Full ms2 445.12@hcd28.00 [110.00-1000.00]");
    assert_eq!(
        (group.analyzer, group.activation),
        (MassAnalyzer::Orbitrap, Activation::Hcd)
    );
    let group = parse("ITMS + c NSI r d Full ms2 652.33@cid35.00 [165.00-1315.00]");
    assert_eq!(
        (group.analyzer, group.activation),
        (MassAnalyzer::IonTrap, Activation::Cid)
    );
    let group = parse("FTMS + p NSI d Full ms2 700.00@etd25.00@hcd20.00 [120.00-2000.00]");
    assert_eq!(group.activation, Activation::Ethcd);
    let group = parse("FTMS + p NSI Full ms [350.00-1500.00]");
    assert_eq!(
        (group.analyzer, group.activation),
        (MassAnalyzer::Orbitrap, Activation::Unknown)
    );
    let group = parse("ASTMS + c NSI d Full ms2 500.00@hcd25.00 [150.00-2000.00]");
    assert_eq!(group.analyzer, MassAnalyzer::Astral);
    assert_eq!(group.label(), "astral/hcd");
}

#[test]
fn auto_tolerance_narrows_to_corrected_spread() {
    // A 3 ppm offset with +/-5 ppm uniform noise.
    let data = points(3000, |_, _| 3.0, 5.0);
    let selection = select_model(&data, RecalibrationOptions::default());
    assert_eq!(selection.kind(), MassModelKind::Static);
    let spread = selection.spread.as_ref().expect("residual spread");
    // After the fit-set offset, residuals are centered.
    assert!(
        spread.fit_center_ppm.abs() < 0.3,
        "{}",
        spread.fit_center_ppm
    );
    assert_eq!(spread.fit_residuals.len(), spread.fit_points);

    let configured = Tolerance::Ppm(-20.0, 20.0);
    let estimate = auto_tolerance(&selection, configured, AutoToleranceOptions::precursor());
    assert!(estimate.narrowed);
    assert!(estimate.skipped.is_none());
    let Tolerance::Ppm(lo, hi) = estimate.tolerance else {
        panic!("ppm tolerance expected");
    };
    // 1.2 x the 99th percentile of |residual| (~4.95 ppm).
    assert!((estimate.fit_quantile_ppm - 4.95).abs() < 0.2);
    assert!((estimate.validation_quantile_ppm - estimate.fit_quantile_ppm).abs() < 0.3);
    assert!((estimate.half_width_ppm - 1.2 * estimate.fit_quantile_ppm).abs() < 1e-4);
    assert!((hi - lo - 2.0 * estimate.half_width_ppm).abs() < 1e-3);
    assert!(hi > 5.5 && hi < 6.5 && lo < -5.5 && lo > -6.5, "{lo} {hi}");
    assert_eq!(estimate.validation_coverage, 1.0);
    assert_eq!(estimate.configured, configured);

    // A configured side inside the estimate is kept.
    let estimate = auto_tolerance(
        &selection,
        Tolerance::Ppm(-5.5, 20.0),
        AutoToleranceOptions::precursor(),
    );
    let Tolerance::Ppm(lo, hi) = estimate.tolerance else {
        panic!("ppm tolerance expected");
    };
    assert_eq!(lo, -5.5);
    assert!(hi < 6.5);

    // A window that would drop held-out residuals is not applied.
    let estimate = auto_tolerance(
        &selection,
        Tolerance::Ppm(-4.0, 20.0),
        AutoToleranceOptions::precursor(),
    );
    assert!(!estimate.narrowed);
    assert_eq!(estimate.skipped.as_deref(), Some("low_coverage"));
    assert!(estimate.validation_coverage < 0.98);
    assert_eq!(estimate.tolerance, Tolerance::Ppm(-4.0, 20.0));
}

#[test]
fn auto_tolerance_follows_heavy_tails() {
    // A +/-1 ppm core with 5% of PSMs spread over +/-8 ppm: four robust
    // sigmas (~3 ppm) would cut the tail, the 99th percentile does not.
    let mut data = points(3000, |_, _| 0.0, 1.0);
    for (i, point) in data.iter_mut().enumerate() {
        if i % 20 == 0 {
            point.error_ppm = 8.0 * noise(i as u64 + 99);
        }
    }
    let selection = select_model(&data, RecalibrationOptions::default());
    let estimate = auto_tolerance(
        &selection,
        Tolerance::Ppm(-20.0, 20.0),
        AutoToleranceOptions::precursor(),
    );
    assert!(estimate.narrowed);
    assert!(4.0 * estimate.fit_sigma_ppm < 4.0);
    assert!(estimate.half_width_ppm > 6.0, "{}", estimate.half_width_ppm);
    assert!(estimate.validation_coverage >= 0.98);

    // Isotope-error matches are searched with the same window, so their
    // (wider) residuals widen it, without changing the model.
    let mut selection = select_model(
        &points(3000, |_, _| 0.0, 1.0),
        RecalibrationOptions::default(),
    );
    let before = auto_tolerance(
        &selection,
        Tolerance::Ppm(-20.0, 20.0),
        AutoToleranceOptions::precursor(),
    );
    assert_eq!(before.half_width_ppm, 3.0);
    let isotope = points(1000, |_, _| 0.0, 6.0)
        .into_iter()
        .map(|p| MassErrorPoint {
            group: p.group + 1_000_000,
            ..p
        })
        .collect::<Vec<_>>();
    let model = selection.model.clone();
    let spread = selection.spread.as_mut().expect("residual spread");
    let validation_psms = spread.validation_psms;
    spread.add_isotope_error_points(model.as_ref(), &isotope, 3);
    assert_eq!(spread.isotope_error_points, 1000);
    assert!(spread.validation_psms > validation_psms + 200);
    let after = auto_tolerance(
        &selection,
        Tolerance::Ppm(-20.0, 20.0),
        AutoToleranceOptions::precursor(),
    );
    assert_eq!(after.isotope_error_points, 1000);
    assert!(after.half_width_ppm > 6.0, "{}", after.half_width_ppm);
    assert!(after.validation_coverage >= 0.98);
}

#[test]
fn auto_tolerance_respects_floor_limits_and_units() {
    // Very tight residuals hit the floor.
    let tight = select_model(
        &points(3000, |_, _| 0.0, 0.2),
        RecalibrationOptions::default(),
    );
    let estimate = auto_tolerance(
        &tight,
        Tolerance::Ppm(-20.0, 20.0),
        AutoToleranceOptions::fragment(),
    );
    let Tolerance::Ppm(lo, hi) = estimate.tolerance else {
        panic!("ppm tolerance expected");
    };
    assert!((hi - lo - 10.0).abs() < 1e-3, "{lo} {hi}");

    // Wide residuals never widen the configured window.
    let wide = select_model(
        &points(3000, |_, _| 0.0, 9.0),
        RecalibrationOptions::default(),
    );
    let estimate = auto_tolerance(
        &wide,
        Tolerance::Ppm(-10.0, 10.0),
        AutoToleranceOptions::precursor(),
    );
    assert!(!estimate.narrowed);
    assert_eq!(estimate.tolerance, Tolerance::Ppm(-10.0, 10.0));

    // Da tolerances are left alone.
    let estimate = auto_tolerance(
        &tight,
        Tolerance::Da(-0.02, 0.02),
        AutoToleranceOptions::fragment(),
    );
    assert!(!estimate.narrowed);
    assert_eq!(estimate.skipped.as_deref(), Some("not_ppm"));

    // Too few held-out PSMs keep the configured window.
    let few = select_model(
        &points(200, |_, _| 0.0, 0.5),
        RecalibrationOptions::default(),
    );
    let estimate = auto_tolerance(
        &few,
        Tolerance::Ppm(-20.0, 20.0),
        AutoToleranceOptions::precursor(),
    );
    assert!(!estimate.narrowed);
    assert_eq!(estimate.skipped.as_deref(), Some("too_few_psms"));

    // Model selection disabled (mass_recalibration off) still reports the
    // raw spread, so tolerances can be narrowed without a correction.
    let disabled = select_model(
        &points(3000, |_, _| 2.0, 1.0),
        RecalibrationOptions {
            max_kind: MassModelKind::None,
            ..RecalibrationOptions::default()
        },
    );
    assert_eq!(disabled.skipped.as_deref(), Some("disabled"));
    assert!(disabled.model.is_none());
    let estimate = auto_tolerance(
        &disabled,
        Tolerance::Ppm(-20.0, 20.0),
        AutoToleranceOptions::precursor(),
    );
    assert!(estimate.narrowed);
    // Uncorrected: the window follows the raw 2 ppm offset.
    assert!((estimate.center_ppm - 2.0).abs() < 0.2);

    // Low-accuracy groups get no tolerance either.
    let skipped = ModelSelection::skipped(500, 5000, "low_accuracy_analyzer");
    let estimate = auto_tolerance(
        &skipped,
        Tolerance::Ppm(-20.0, 20.0),
        AutoToleranceOptions::fragment(),
    );
    assert_eq!(estimate.skipped.as_deref(), Some("low_accuracy_analyzer"));
    assert!(!estimate.narrowed);
}

#[test]
fn tolerances_are_looked_up_by_file_and_group() {
    use crate::spectrum::{AcquisitionGroup, Activation, MassAnalyzer};
    let hcd = AcquisitionGroup {
        analyzer: MassAnalyzer::Orbitrap,
        activation: Activation::Hcd,
    };
    let recalibration = MassRecalibration {
        files: vec![
            FileMassCorrection::default(),
            FileMassCorrection {
                precursor: None,
                fragment: vec![GroupMassCorrection {
                    group: hcd,
                    model: None,
                    fragment_tol: Some(Tolerance::Ppm(-6.0, 6.0)),
                }],
                precursor_tol: Some(Tolerance::Ppm(-4.0, 4.0)),
            },
        ],
    };
    assert_eq!(recalibration.tolerances(0, hcd), (None, None));
    assert_eq!(
        recalibration.tolerances(1, hcd),
        (
            Some(Tolerance::Ppm(-4.0, 4.0)),
            Some(Tolerance::Ppm(-6.0, 6.0))
        )
    );
    assert_eq!(
        recalibration.tolerances(1, AcquisitionGroup::default()),
        (Some(Tolerance::Ppm(-4.0, 4.0)), None)
    );
    // Tolerances alone do not change masses.
    assert!(recalibration.files[1].is_identity());
    assert!(recalibration.files[1].tunes_tolerances());
    assert!(recalibration.file(1).is_none());
}
