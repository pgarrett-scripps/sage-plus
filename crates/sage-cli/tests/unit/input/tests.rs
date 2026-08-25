use super::LfqOptions;
use sage_core::lfq::LfqSettings;

#[test]
fn parses_lfq_rt_percent_tolerance() {
    let options: LfqOptions = serde_json::from_str(r#"{"rt_pct_tolerance": 1.25}"#).unwrap();
    let settings: LfqSettings = options.into();

    assert_eq!(settings.rt_pct_tolerance, 1.25);
}

#[test]
fn defaults_lfq_rt_percent_tolerance() {
    let options: LfqOptions = serde_json::from_str("{}").unwrap();
    let settings: LfqSettings = options.into();

    assert_eq!(settings.rt_pct_tolerance, 0.5);
}
