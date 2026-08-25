use super::*;

#[test]
fn test_identify_format() {
    assert_eq!(FileFormat::from("foo.mzml"), FileFormat::MzML);
    assert_eq!(FileFormat::from("foo.mzML"), FileFormat::MzML);
    assert_eq!(FileFormat::from("foo.mgf"), FileFormat::MGF);
    assert_eq!(FileFormat::from("foo.mgf.gz"), FileFormat::MGF);
    assert_eq!(FileFormat::from("foo.tdf"), FileFormat::TDF);
    assert_eq!(FileFormat::from("./tomato/foo.d"), FileFormat::TDF);
    assert_eq!(FileFormat::from("./tomato/foo.d/"), FileFormat::TDF);
    assert_eq!(FileFormat::from("foo.raw"), FileFormat::ThermoRaw);
    assert_eq!(FileFormat::from("foo.RAW"), FileFormat::ThermoRaw);
}

#[test]
fn thermoraw_rejects_signal_to_noise_mode() {
    let url = Url::parse("file:///tmp/example.raw").unwrap();
    assert!(matches!(
        read_thermoraw(&url, 0, Some(2)),
        Err(Error::Unsupported(_))
    ));
}

#[test]
fn unidentified_spectra_format_returns_an_error() {
    let url = Url::parse("file:///tmp/example.unknown").unwrap();
    assert!(matches!(
        read_spectra(&url, 0, None, BrukerProcessingConfig::default(), false),
        Err(Error::Unsupported(message)) if message.contains("determine the spectra format")
    ));
}
