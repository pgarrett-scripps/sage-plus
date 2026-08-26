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
    assert!(settings.mbr);
}

#[test]
fn parses_disabled_match_between_runs() {
    let options: LfqOptions = serde_json::from_str(r#"{"mbr": false}"#).unwrap();
    let settings: LfqSettings = options.into();

    assert!(!settings.mbr);
}
