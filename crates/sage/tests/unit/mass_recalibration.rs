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
            })
            .collect(),
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
    let group = parse("FTMS + c NSI d sa Full ms2 700.00@etd25.00@cid20.00 [120.00-2000.00]");
    assert_eq!(group.activation, Activation::Etcid);
    assert_eq!(group.label(), "orbitrap/etcid");
    let group = parse("FTMS + c NSI d Full ms2 700.00@etd25.00 [120.00-2000.00]");
    assert_eq!(group.activation, Activation::Etd);
    let group = parse("FTMS + p NSI Full ms [350.00-1500.00]");
    assert_eq!(
        (group.analyzer, group.activation),
        (MassAnalyzer::Orbitrap, Activation::Unknown)
    );
    let group = parse("ASTMS + c NSI d Full ms2 500.00@hcd25.00 [150.00-2000.00]");
    assert_eq!(group.analyzer, MassAnalyzer::Astral);
    assert_eq!(group.label(), "astral/hcd");

    // Other Thermo analyzers, and filters that name none.
    assert_eq!(
        parse("TQMS + c NSI SRM ms2 500.00@cid25.00 [150.00-160.00]").analyzer,
        MassAnalyzer::Other
    );
    let group = parse("+ c ESI Full ms2 500.00@hcd25.00 [150.00-2000.00]");
    assert_eq!(
        (group.analyzer, group.activation),
        (MassAnalyzer::Unknown, Activation::Hcd)
    );
    // Text that is not a Thermo filter never becomes `other`.
    assert_eq!(
        parse("scan from vendor@site notes"),
        AcquisitionGroup::default()
    );
}

#[test]
fn filter_overrides_configured_sources() {
    use crate::spectrum::{AcquisitionGroup, Activation, MassAnalyzer};
    // The filter keeps the ETD a converter dropped from the activation terms.
    let group = AcquisitionGroup::resolve(
        Some("ITMS + c NSI r d sa Full ms3 700.00@etd25.00@hcd20.00 [120.00-2000.00]"),
        MassAnalyzer::Orbitrap,
        Activation::Hcd,
    );
    assert_eq!(group.label(), "ion_trap/ethcd");
    // Without a usable filter, the configured analyzer and terms are used.
    for filter in [None, Some("not a filter")] {
        let group = AcquisitionGroup::resolve(filter, MassAnalyzer::Tof, Activation::Cid);
        assert_eq!(group.label(), "tof/cid");
    }
    assert_eq!(Activation::Etd.combine(Activation::Hcd), Activation::Ethcd);
    assert_eq!(
        MassAnalyzer::Other.combine(MassAnalyzer::Orbitrap),
        MassAnalyzer::Orbitrap
    );
    assert_eq!(
        MassAnalyzer::Orbitrap.combine(MassAnalyzer::Other),
        MassAnalyzer::Orbitrap
    );
    assert_eq!(
        MassAnalyzer::from_psi_ms("MS:1003379"),
        Some(MassAnalyzer::Astral)
    );
}

#[test]
fn discovery_sample_is_stratified_by_acquisition_group() {
    // A strict 4-scan cycle (two Orbitrap scans, then two ion-trap scans),
    // where a stride of 4 would sample only the first scan type.
    let cycle = [0u8, 1, 2, 0];
    let keys = (0..100_000).map(|i| cycle[i % 4]).collect::<Vec<_>>();
    let sample = stratified_sample(&keys, 25_000);
    assert!(sample.len() <= 25_000);
    assert!(sample.windows(2).all(|pair| pair[0] < pair[1]));
    let count = |key: u8| sample.iter().filter(|&&i| keys[i] == key).count();
    assert_eq!((count(0), count(1), count(2)), (12_500, 6_250, 6_250));
    // Every group spans the gradient.
    for key in 0..3 {
        let members = sample
            .iter()
            .filter(|&&i| keys[i] == key)
            .collect::<Vec<_>>();
        assert!(*members[0] < 100 && *members[members.len() - 1] > 99_900);
    }
    assert_eq!(sample, stratified_sample(&keys, 25_000));

    // Small inputs are searched whole; rare groups keep at least one spectrum.
    assert_eq!(stratified_sample(&[1, 2, 3], 5), vec![0, 1, 2]);
    let mut keys = vec![0u8; 1_000];
    keys[500] = 1;
    let sample = stratified_sample(&keys, 100);
    assert!(sample.contains(&500));
    assert_eq!(sample.len(), 100);
}

#[test]
fn discovery_q_values_are_computed_per_acquisition_group() {
    use crate::scoring::Feature;
    let psm = |poisson: f64, label: i32| Feature {
        poisson,
        label,
        rank: 1,
        ..Default::default()
    };
    let mut psms = Vec::new();
    // Accurate group: 200 targets, no decoys.
    for i in 0..200 {
        psms.push(("orbitrap", i, psm(-20.0 + i as f64 * 0.05, 1)));
    }
    // Wide-window group: decoys score better than every accurate target.
    for i in 0..400 {
        let label = if i % 2 == 0 { -1 } else { 1 };
        psms.push(("ion_trap", 200 + i, psm(-40.0 + i as f64 * 0.05, label)));
    }
    let pooled = {
        let mut features = psms.iter().map(|(_, _, f)| f.clone()).collect::<Vec<_>>();
        features.sort_by(|a, b| a.poisson.total_cmp(&b.poisson));
        crate::ml::qvalue::spectrum_q_value_by(&mut features, |f| f.poisson)
    };
    assert_eq!(pooled, 0);
    let confident = confident_per_group(psms, 0.01);
    // All accurate targets pass; the wide-window group has none.
    assert_eq!(confident.len(), 200);
    assert!(confident.iter().all(|(index, _)| *index < 200));
    assert!(confident.windows(2).all(|pair| pair[0].0 < pair[1].0));
}
