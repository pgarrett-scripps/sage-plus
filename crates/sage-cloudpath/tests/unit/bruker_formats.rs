use super::*;
use crate::tdf::{BrukerProcessingConfig, TdfReader};
use crate::FileFormat;

fn data(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data/bruker")
        .join(name)
}

fn parse(path: &Path, requires_ms1: bool) -> Result<Vec<RawSpectrum>, crate::Error> {
    TdfReader.parse(path, 3, BrukerProcessingConfig::default(), requires_ms1)
}

#[test]
fn tsf_ms2_spectra_take_precursors_from_frame_msms_info() {
    let path = data("example_tsf.d");
    let spectra = parse(&path, false).unwrap();
    // Frame 1 is MS1; frames 2 and 3 are MS/MS, frame 3 without peaks.
    assert_eq!(spectra.len(), 2);
    let first = &spectra[0];
    assert_eq!(first.id, "2");
    assert_eq!(first.ms_level, 2);
    assert_eq!(first.file_id, 3);
    assert_eq!(first.intensity, [7.0, 8.0]);
    assert_eq!(first.mz.len(), 2);
    assert!(first.mz.windows(2).all(|w| w[0] < w[1]));
    let precursor = &first.precursors[0];
    assert_eq!(precursor.mz, 44.99850455681751_f64 as f32);
    assert_eq!(precursor.charge, None);
    assert_eq!(precursor.spectrum_ref.as_deref(), Some("1"));
    assert_eq!(precursor.isolation_window, Some(Tolerance::Da(-1.0, 1.0)));
    assert_eq!(precursor.inverse_ion_mobility, None);
    assert_eq!(first.scan_start_time, 17.444204_f64 as f32 / 60.0);
    let second = &spectra[1];
    assert_eq!(second.id, "3");
    assert!(second.mz.is_empty());
    assert_eq!(second.precursors[0].charge, Some(2));
    assert_eq!(
        second.precursors[0].isolation_window,
        Some(Tolerance::Da(-1.5, 1.5))
    );
}

#[test]
fn tsf_ms1_uses_the_calibrated_model_and_files_inside_resolve() {
    let spectra = parse(&data("example_tsf.d/analysis.tsf"), true).unwrap();
    assert_eq!(spectra.len(), 3);
    let ms1 = spectra.iter().find(|s| s.ms_level == 1).unwrap();
    assert_eq!(ms1.id, "1");
    assert!(ms1.precursors.is_empty());
    assert_eq!(ms1.intensity, [10.0, 250.5, 3.25]);
    assert_eq!(ms1.total_ion_current, 263.75);
    // Bruker SDK m/z for these TOF indices; the linear timsrust scale is off by ~10 ppm.
    for (mz, sdk) in ms1
        .mz
        .iter()
        .zip([20.155193851563446, 111.04787914506481, 499.86573762004514])
    {
        assert_eq!(*mz, sdk as f32);
    }
}

#[test]
fn minitdf_spectra_match_timsrust() {
    let path = data("example_minitdf.ms2");
    let spectra = parse(&path, false).unwrap();
    assert_eq!(spectra.len(), 2);
    let first = &spectra[0];
    assert_eq!(first.id, "0");
    assert_eq!(first.mz, [100.0, 200.002, 300.03, 400.4]);
    assert_eq!(first.intensity, [1.0, 2.0, 3.0, 4.0]);
    assert_eq!(first.scan_start_time, 12.345_f64 as f32 / 60.0);
    assert_eq!(first.ion_injection_time, 12.345_f64 as f32);
    let precursor = &first.precursors[0];
    assert_eq!(precursor.mz, 123.4567_f64 as f32);
    assert_eq!(precursor.charge, Some(1));
    assert_eq!(precursor.inverse_ion_mobility, Some(1.234_f64 as f32));
    assert_eq!(precursor.isolation_window, Some(Tolerance::Da(-1.0, 1.0)));
    let second = &spectra[1];
    assert_eq!(second.id, "1");
    assert_eq!(second.precursors[0].charge, Some(2));
    assert_eq!(
        second.precursors[0].isolation_window,
        Some(Tolerance::Da(-1.5, 1.5))
    );
    // A spectrum file inside the directory resolves to the same run.
    let file = parse(&path.join("converter.ms2spectrum.bin"), false).unwrap();
    assert_eq!(file.len(), 2);
}

#[test]
fn minitdf_with_lfq_is_unsupported() {
    let error = parse(&data("example_minitdf.ms2"), true).unwrap_err();
    assert!(
        matches!(&error, crate::Error::Unsupported(message) if message.contains("MS1")),
        "{error}"
    );
}

#[test]
fn timsrust_parquet_spectra_stay_unsupported() {
    let directory = tempfile::tempdir().unwrap();
    let path = directory.path().join("run.d");
    std::fs::create_dir(&path).unwrap();
    std::fs::write(path.join("spectra.parquet"), b"PAR1").unwrap();
    let error = parse(&path, false).unwrap_err();
    assert!(
        matches!(&error, crate::Error::Unsupported(message) if message.contains("parquet")),
        "{error}"
    );
}

#[test]
fn tsf_and_minitdf_paths_are_bruker() {
    for path in [
        "run.d/analysis.tsf",
        "run.d/analysis.tsf_bin",
        "run/ms2writer.ms2spectrum.bin",
        "run/ms2writer.ms2spectrum.parquet",
    ] {
        assert_eq!(FileFormat::from(path), FileFormat::TDF, "{path}");
    }
    assert_eq!(
        FileFormat::from("results.parquet"),
        FileFormat::Unidentified
    );
}
