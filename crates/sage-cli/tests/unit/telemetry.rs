use super::*;

#[test]
fn serialization_contains_only_the_documented_aggregate_fields() {
    let telemetry = Telemetry {
        version: "1.2.3".into(),
        peptides: 10,
        fragments: 20,
        files: 2,
        runtime_secs: 30,
        lfq: true,
        tmt: Some(Isobaric::Tmt6),
        os_name: "test-os".into(),
        total_memory: 4096,
        cpus: 8,
    };

    let value = serde_json::to_value(telemetry).unwrap();
    let object = value.as_object().unwrap();

    assert_eq!(object.len(), 10);
    assert_eq!(value["version"], "1.2.3");
    assert_eq!(value["peptides"], 10);
    assert_eq!(value["fragments"], 20);
    assert_eq!(value["files"], 2);
    assert_eq!(value["runtime_secs"], 30);
    assert_eq!(value["lfq"], true);
    assert_eq!(value["tmt"], "Tmt6");
    assert_eq!(value["os_name"], "test-os");
    assert_eq!(value["total_memory"], 4096);
    assert_eq!(value["cpus"], 8);
}

#[test]
fn serialization_represents_missing_tmt_as_null() {
    let telemetry = Telemetry {
        version: String::new(),
        peptides: 0,
        fragments: 0,
        files: 0,
        runtime_secs: 0,
        lfq: false,
        tmt: None,
        os_name: String::new(),
        total_memory: 0,
        cpus: 1,
    };

    assert!(serde_json::to_value(telemetry).unwrap()["tmt"].is_null());
}
