use super::*;

fn temporary_root() -> PathBuf {
    let root = std::env::temp_dir().join(format!("sage-mcp-test-{}", Uuid::new_v4()));
    fs::create_dir_all(&root).unwrap();
    root
}

fn test_job(state: &State, job_id: &str, status: JobStatus, created_at_unix: u64) -> JobRecord {
    let job_directory = state.jobs_dir.join(job_id);
    let output_directory = job_directory.join("output");
    fs::create_dir_all(&output_directory).unwrap();
    let events_path = job_directory.join("events.jsonl");
    fs::write(&events_path, "").unwrap();
    JobRecord {
        job_id: job_id.into(),
        status,
        config_path: state
            .root
            .join("config.json")
            .to_string_lossy()
            .into_owned(),
        job_directory: job_directory.to_string_lossy().into_owned(),
        events_path: events_path.to_string_lossy().into_owned(),
        output_directory: output_directory.to_string_lossy().into_owned(),
        created_at_unix,
        updated_at_unix: created_at_unix,
        summary: None,
        error: None,
    }
}

#[test]
fn rejects_paths_outside_root() {
    let temp = temporary_root();
    let state = State::new(temp.clone(), None).unwrap();
    let error = state
        .resolve_existing("../outside")
        .unwrap_err()
        .to_string();
    assert!(error.contains("cannot access") || error.contains("outside"));
    fs::remove_dir_all(temp).unwrap();
}

#[test]
fn explicit_approval_is_required_before_input_access() {
    let temp = temporary_root();
    let state = State::new(temp.clone(), None).unwrap();
    let error = state
        .start_search(StartSearchArgs {
            config_path: "missing.json".into(),
            approved: false,
            batch_size: None,
        })
        .unwrap_err()
        .to_string();
    assert!(error.contains("explicit approval"));
    fs::remove_dir_all(temp).unwrap();
}

#[test]
fn rejects_remote_urls_before_attempting_io() {
    let temp = temporary_root();
    let state = State::new(temp.clone(), None).unwrap();

    let error = state
        .resolve_existing("s3://private-bucket/config.json")
        .unwrap_err()
        .to_string();
    assert!(error.contains("remote URLs are disabled"));
    fs::remove_dir_all(temp).unwrap();
}

#[test]
fn restart_marks_incomplete_jobs_as_interrupted_and_persists_status() {
    let temp = temporary_root();
    let state = State::new(temp.clone(), None).unwrap();
    let record = test_job(&state, "interrupted-job", JobStatus::Running, 10);
    write_record(&record).unwrap();
    drop(state);

    let restored = State::new(temp.clone(), None).unwrap();
    let record = restored.job("interrupted-job").unwrap();
    assert_eq!(record.status, JobStatus::Interrupted);
    let persisted: JobRecord = serde_json::from_slice(
        &fs::read(Path::new(&record.job_directory).join("job.json")).unwrap(),
    )
    .unwrap();
    assert_eq!(persisted.status, JobStatus::Interrupted);
    assert!(record.updated_at_unix >= 10);
    fs::remove_dir_all(temp).unwrap();
}

#[test]
fn corrupt_historical_job_does_not_block_valid_restoration() {
    let temp = temporary_root();
    let state = State::new(temp.clone(), None).unwrap();
    let valid = test_job(&state, "valid-job", JobStatus::Completed, 10);
    write_record(&valid).unwrap();
    let corrupt = state.jobs_dir.join("corrupt-job");
    fs::create_dir_all(&corrupt).unwrap();
    fs::write(corrupt.join("job.json"), "{definitely not json").unwrap();
    drop(state);

    let restored = State::new(temp.clone(), None).unwrap();
    assert_eq!(
        restored.job("valid-job").unwrap().status,
        JobStatus::Completed
    );
    assert!(restored.job("corrupt-job").is_err());
    fs::remove_dir_all(temp).unwrap();
}

#[test]
fn worker_panics_are_captured_for_terminal_job_updates() {
    let result = catch_job_panic(|| panic!("test worker failure"));
    assert_eq!(result.unwrap_err(), "test worker failure");
}

#[test]
fn cancellation_updates_manifest_and_signals_worker() {
    let temp = temporary_root();
    let state = State::new(temp.clone(), None).unwrap();
    let record = test_job(&state, "running-job", JobStatus::Running, 10);
    let cancellation = CancellationToken::default();
    state.jobs.write().unwrap().insert(
        record.job_id.clone(),
        JobEntry {
            record,
            cancellation: Some(cancellation.clone()),
            memory_limited: false,
        },
    );

    let cancelled = state.cancel("running-job").unwrap();
    assert_eq!(cancelled.status, JobStatus::Cancelling);
    assert!(cancellation.is_cancelled());
    let persisted: JobRecord = serde_json::from_slice(
        &fs::read(Path::new(&cancelled.job_directory).join("job.json")).unwrap(),
    )
    .unwrap();
    assert_eq!(persisted.status, JobStatus::Cancelling);
    assert!(state
        .cancel("running-job")
        .unwrap_err()
        .to_string()
        .contains("not running"));
    fs::remove_dir_all(temp).unwrap();
}

#[test]
fn memory_limited_jobs_are_exclusive() {
    let temp = temporary_root();
    let state = State::new(temp.clone(), None).unwrap();
    let record = test_job(&state, "limited-job", JobStatus::Running, 10);
    state.jobs.write().unwrap().insert(
        record.job_id.clone(),
        JobEntry {
            record,
            cancellation: Some(CancellationToken::default()),
            memory_limited: true,
        },
    );

    assert!(state.ensure_memory_execution_available(false).is_err());
    assert!(state.ensure_memory_execution_available(true).is_err());
    fs::remove_dir_all(temp).unwrap();
}

#[test]
fn job_events_are_incremental_bounded_and_available_as_resources() {
    let temp = temporary_root();
    let state = State::new(temp.clone(), None).unwrap();
    let record = test_job(&state, "event-job", JobStatus::Running, 10);
    let event_text = (0..25)
        .map(|sequence| serde_json::json!({ "sequence": sequence }).to_string())
        .collect::<Vec<_>>()
        .join("\n");
    fs::write(&record.events_path, format!("{event_text}\n")).unwrap();
    write_record(&record).unwrap();
    state.jobs.write().unwrap().insert(
        record.job_id.clone(),
        JobEntry {
            record,
            cancellation: None,
            memory_limited: false,
        },
    );

    let events = state
        .events(JobEventsArgs {
            job_id: "event-job".into(),
            after_sequence: Some(22),
            limit: Some(2),
        })
        .unwrap();
    assert_eq!(events.len(), 2);
    assert_eq!(events[0]["sequence"], 23);
    assert_eq!(events[1]["sequence"], 24);

    let summary = state.summarize("event-job").unwrap();
    assert_eq!(summary["recent_events"].as_array().unwrap().len(), 20);
    assert_eq!(summary["recent_events"][0]["sequence"], 5);
    assert!(state
        .resource("sage://jobs/event-job/events")
        .unwrap()
        .contains("\"sequence\":24"));
    assert!(state
        .resource("sage://jobs/event-job/manifest")
        .unwrap()
        .contains("event-job"));
    fs::remove_dir_all(temp).unwrap();
}

#[test]
fn jobs_are_sorted_newest_first_and_unknown_ids_are_rejected() {
    let temp = temporary_root();
    let state = State::new(temp.clone(), None).unwrap();
    for (job_id, created) in [("older", 1), ("newer", 2)] {
        let record = test_job(&state, job_id, JobStatus::Completed, created);
        state.jobs.write().unwrap().insert(
            record.job_id.clone(),
            JobEntry {
                record,
                cancellation: None,
                memory_limited: false,
            },
        );
    }

    let jobs = state.list_jobs();
    assert_eq!(jobs[0].job_id, "newer");
    assert_eq!(jobs[1].job_id, "older");
    assert!(state
        .job("missing")
        .unwrap_err()
        .to_string()
        .contains("unknown Sage job"));
    fs::remove_dir_all(temp).unwrap();
}

#[test]
fn ratios_avoid_division_by_zero() {
    assert_eq!(ratio(3, 2), Some(1.5));
    assert_eq!(percent(1, 4), Some(25.0));
    assert_eq!(ratio(1, 0), None);
    assert_eq!(percent(1, 0), None);
}
