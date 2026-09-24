use super::*;

fn fixture() -> serde_json::Value {
    let path =
        Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/bruker/tims_calibration_sdk.json");
    serde_json::from_str(&std::fs::read_to_string(path).unwrap()).unwrap()
}

fn coefficients(case: &serde_json::Value) -> [f64; 10] {
    let values = case["coefficients"].as_array().unwrap();
    std::array::from_fn(|index| values[index].as_f64().unwrap())
}

/// Write a minimal `analysis.tdf` with the given calibration rows and
/// `(frame, calibration)` assignments.
fn write_tdf(rows: &[(i64, i64, [f64; 10])], frames: &[(i64, i64)]) -> tempfile::TempDir {
    let directory = tempfile::tempdir().unwrap();
    let connection = rusqlite::Connection::open(directory.path().join("analysis.tdf")).unwrap();
    connection
        .execute_batch(
            "CREATE TABLE TimsCalibration (Id INTEGER PRIMARY KEY, ModelType INTEGER, \
             C0 REAL, C1 REAL, C2 REAL, C3 REAL, C4 REAL, C5 REAL, C6 REAL, C7 REAL, C8 REAL, C9 REAL);
             CREATE TABLE Frames (Id INTEGER PRIMARY KEY, TimsCalibration INTEGER);
             CREATE TABLE Precursors (Id INTEGER PRIMARY KEY, Parent INTEGER, ScanNumber REAL);",
        )
        .unwrap();
    for (id, model, c) in rows {
        connection
            .execute(
                "INSERT INTO TimsCalibration VALUES (?1, ?2, ?3, ?4, ?5, ?6, ?7, ?8, ?9, ?10, ?11, ?12)",
                rusqlite::params![id, model, c[0], c[1], c[2], c[3], c[4], c[5], c[6], c[7], c[8], c[9]],
            )
            .unwrap();
    }
    for (frame, calibration) in frames {
        connection
            .execute("INSERT INTO Frames VALUES (?1, ?2)", [frame, calibration])
            .unwrap();
    }
    directory
}

#[test]
fn model_reproduces_bruker_sdk_values() {
    let fixture = fixture();
    let cases = fixture["cases"].as_array().unwrap();
    assert_eq!(cases.len(), 4);
    let mut compared = 0;
    for case in cases {
        let model = MobilityModel::new(coefficients(case));
        let scans = case["scans"].as_array().unwrap();
        let expected = case["one_over_k0"].as_array().unwrap();
        for (scan, expected) in scans.iter().zip(expected) {
            let (scan, expected) = (scan.as_f64().unwrap(), expected.as_f64().unwrap());
            let actual = model.one_over_k0(scan);
            assert!(
                (actual - expected).abs() <= 1e-12,
                "{} calibration {}: scan {scan} gave {actual}, SDK {expected}",
                case["source"],
                case["calibration_id"]
            );
            compared += 1;
        }
    }
    assert!(compared > 60);
}

#[test]
fn frames_use_their_own_calibration_row() {
    let fixture = fixture();
    let prm = fixture["cases"]
        .as_array()
        .unwrap()
        .iter()
        .filter(|case| case["source"] == "example_prm.d")
        .collect::<Vec<_>>();
    assert_eq!(prm.len(), 2);
    let (first, second) = (coefficients(prm[0]), coefficients(prm[1]));
    let directory = write_tdf(
        &[(1, 2, first), (2, 2, second)],
        &[(1, 1), (2, 1), (3, 2), (4, 2), (5, 2)],
    );
    let calibration = MobilityCalibration::from_path(directory.path()).unwrap();
    let scan = 700.25;
    assert_eq!(
        calibration.frame(1).one_over_k0(scan),
        MobilityModel::new(first).one_over_k0(scan)
    );
    assert_eq!(
        calibration.frame(4).one_over_k0(scan),
        MobilityModel::new(second).one_over_k0(scan)
    );
    assert_ne!(
        calibration.frame(1).one_over_k0(scan),
        calibration.frame(4).one_over_k0(scan)
    );
    // Row 2 is used by more frames, and unknown frames fall back to it.
    assert_eq!(calibration.dominant(), &MobilityModel::new(second));
    assert_eq!(calibration.frame(99), &MobilityModel::new(second));
}

#[test]
fn unsupported_models_and_missing_rows_are_errors() {
    let c = coefficients(&fixture()["cases"][0]);
    let unsupported = write_tdf(&[(1, 1, c)], &[(1, 1)]);
    let error = MobilityCalibration::from_path(unsupported.path()).unwrap_err();
    assert!(matches!(
        error,
        MobilityCalibrationError::UnsupportedModel {
            id: 1,
            model: 1,
            ..
        }
    ));
    assert!(error.to_string().contains("ion_mobility_scale"));

    let missing = write_tdf(&[(1, 2, c)], &[(1, 1), (2, 7)]);
    assert!(matches!(
        MobilityCalibration::from_path(missing.path()).unwrap_err(),
        MobilityCalibrationError::MissingCalibration {
            frame: 2,
            id: 7,
            ..
        }
    ));

    let empty = write_tdf(&[], &[]);
    assert!(matches!(
        MobilityCalibration::from_path(empty.path()).unwrap_err(),
        MobilityCalibrationError::Empty { .. }
    ));
}

#[test]
fn dda_precursors_keep_fractional_average_scans() {
    let c = coefficients(&fixture()["cases"][0]);
    let directory = write_tdf(&[(1, 2, c)], &[(10, 1)]);
    let connection = rusqlite::Connection::open(directory.path().join("analysis.tdf")).unwrap();
    connection
        .execute("INSERT INTO Precursors VALUES (3, 10, 412.75)", [])
        .unwrap();
    let scans = MobilityCalibration::dda_precursor_scans(directory.path()).unwrap();
    assert_eq!(scans.get(&3), Some(&(10, 412.75)));
}

#[test]
fn scale_names_round_trip() {
    for scale in [BrukerMobilityScale::Calibrated, BrukerMobilityScale::Linear] {
        let json = serde_json::to_string(&scale).unwrap();
        assert_eq!(json, format!("\"{}\"", scale.as_str()));
        assert_eq!(
            serde_json::from_str::<BrukerMobilityScale>(&json).unwrap(),
            scale
        );
    }
    assert_eq!(
        BrukerMobilityScale::default(),
        BrukerMobilityScale::Calibrated
    );
}

#[test]
fn analysis_tdf_is_resolved_like_timsrust() {
    let directory = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/bruker/example_dia.d");
    let tdf = directory.join("analysis.tdf");
    assert_eq!(analysis_tdf(&directory), Some(tdf.clone()));
    assert_eq!(analysis_tdf(&tdf), Some(tdf.clone()));
    assert_eq!(analysis_tdf(directory.join("analysis.tdf_bin")), Some(tdf));
    // Only analysis.tdf, no analysis.tdf_bin: timsrust does not read it as TDF.
    let partial = write_tdf(&[], &[]);
    assert_eq!(analysis_tdf(partial.path()), None);

    let from_bin = MobilityCalibration::from_path(directory.join("analysis.tdf_bin")).unwrap();
    let from_directory = MobilityCalibration::from_path(&directory).unwrap();
    assert_eq!(from_bin.dominant(), from_directory.dominant());
    assert_eq!(
        MobilityCalibration::dda_precursor_scans(directory.join("analysis.tdf_bin")).unwrap(),
        MobilityCalibration::dda_precursor_scans(&directory).unwrap()
    );
}
