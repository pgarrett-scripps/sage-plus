use super::*;

fn options() -> FitOptions {
    FitOptions {
        min_points: 10,
        min_rt_span: 5.0,
        min_linear_improvement: 0.15,
        outlier_mads: 4.0,
    }
}

#[test]
fn fits_static_offset() {
    let points = (0..40)
        .map(|i| CalibrationPoint {
            rt_minutes: i as f32,
            error_ppm: 3.0 + (i % 3) as f32 * 0.05,
        })
        .collect::<Vec<_>>();
    let fit = fit(&points, options()).unwrap();
    assert_eq!(fit.model.kind, CalibrationKind::Static);
    assert!((fit.model.intercept_ppm - 3.05).abs() < 0.06);
}

#[test]
fn selects_material_rt_drift() {
    let points = (0..60)
        .map(|i| {
            let rt = i as f32;
            CalibrationPoint {
                rt_minutes: rt,
                error_ppm: -2.0 + 0.12 * (rt - 30.0) + (i % 2) as f32 * 0.02,
            }
        })
        .collect::<Vec<_>>();
    let fit = fit(&points, options()).unwrap();
    assert_eq!(fit.model.kind, CalibrationKind::RetentionTimeLinear);
    assert!((fit.model.slope_ppm_per_min - 0.12).abs() < 0.002);
    assert!(fit.model_median_abs_residual < fit.static_median_abs_residual);
}

#[test]
fn ignores_gross_outliers() {
    let mut points = (0..40)
        .map(|i| CalibrationPoint {
            rt_minutes: i as f32,
            error_ppm: 2.0 + (i % 2) as f32 * 0.1,
        })
        .collect::<Vec<_>>();
    points.extend((0..5).map(|i| CalibrationPoint {
        rt_minutes: i as f32,
        error_ppm: 100.0,
    }));
    let fit = fit(&points, options()).unwrap();
    assert_eq!(fit.inliers, 40);
    assert!((fit.model.intercept_ppm - 2.05).abs() < 0.06);
}

#[test]
fn correction_removes_predicted_error() {
    let model = CalibrationModel {
        kind: CalibrationKind::RetentionTimeLinear,
        intercept_ppm: 2.0,
        slope_ppm_per_min: 0.1,
        rt_center: 20.0,
    };
    let theoretical = 1000.0;
    let observed = theoretical * (1.0 + model.predict_ppm(30.0) * 1e-6);
    assert!((model.correct_observed(observed, 30.0) - theoretical).abs() < 1e-4);
}

#[test]
fn fragment_alignment_preserves_spread() {
    assert_eq!(align_fragment_error(5.0, 5.0, 5.0), 0.0);
    assert_eq!(align_fragment_error(3.0, 2.0, 2.0), 1.0);
    assert_eq!(align_fragment_error(3.0, 2.0, 0.0), 3.0);
}
