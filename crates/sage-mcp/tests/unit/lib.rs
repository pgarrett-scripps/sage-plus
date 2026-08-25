use super::*;

fn temporary_root() -> PathBuf {
    let root = std::env::temp_dir().join(format!("sage-mcp-test-{}", Uuid::new_v4()));
    fs::create_dir_all(&root).unwrap();
    root
}

fn fixture_root() -> PathBuf {
    let root = temporary_root();
    let tests = root.join("tests");
    fs::create_dir_all(&tests).unwrap();
    let source = Path::new(env!("CARGO_MANIFEST_DIR")).join("../../tests");
    for name in ["config.json", "Q99536.fasta", "LQSRPAAPPAPGPGQLTLR.mzML"] {
        fs::copy(source.join(name), tests.join(name)).unwrap();
    }
    root
}

#[cfg(unix)]
fn worker_script(root: &Path, name: &str, body: &str) -> PathBuf {
    use std::os::unix::fs::PermissionsExt;

    let path = root.join(name);
    fs::write(&path, format!("#!/bin/sh\n{body}\n")).unwrap();
    let mut permissions = fs::metadata(&path).unwrap().permissions();
    permissions.set_mode(0o700);
    fs::set_permissions(&path, permissions).unwrap();
    path
}

fn wait_for_terminal_job(state: &State, job_id: &str) -> JobRecord {
    for _ in 0..100 {
        let record = state.job(job_id).unwrap();
        if matches!(
            record.status,
            JobStatus::Completed | JobStatus::Failed | JobStatus::Cancelled
        ) {
            return record;
        }
        thread::sleep(Duration::from_millis(20));
    }
    panic!("worker job did not reach a terminal state");
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
        worker_pid: None,
        worker_exit_code: None,
    }
}

#[test]
fn rejects_paths_outside_root() {
    let temp = temporary_root();
    let state = State::new(temp.clone(), None).unwrap();
    let error = resolve_existing_under(&state.root, "../outside")
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

    let error = resolve_existing_under(&state.root, "s3://private-bucket/config.json")
        .unwrap_err()
        .to_string();
    assert!(error.contains("remote URLs are disabled"));
    fs::remove_dir_all(temp).unwrap();
}

#[test]
fn nested_database_inputs_are_confined_to_the_mcp_root() {
    let root = fixture_root();
    let outside = std::env::temp_dir().join(format!("sage-mcp-outside-{}.tsv", Uuid::new_v4()));
    fs::write(&outside, "protein\tposition\n").unwrap();
    let config_path = root.join("tests/config.json");
    let mut config: serde_json::Value =
        serde_json::from_slice(&fs::read(&config_path).unwrap()).unwrap();
    config["database"]["custom_cleavage_sites"] =
        serde_json::Value::String(outside.to_string_lossy().into_owned());
    fs::write(&config_path, serde_json::to_vec_pretty(&config).unwrap()).unwrap();

    let error = match load_config_under(&root, &config_path.to_string_lossy()) {
        Ok(_) => panic!("outside nested path was accepted"),
        Err(error) => error.to_string(),
    };
    assert!(error.contains("outside the configured MCP root"));
    fs::remove_file(outside).unwrap();
    fs::remove_dir_all(root).unwrap();
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

#[cfg(unix)]
#[test]
fn cancellation_updates_manifest_and_signals_worker() {
    let temp = temporary_root();
    let state = State::new(temp.clone(), None).unwrap();
    let record = test_job(&state, "running-job", JobStatus::Running, 10);
    let child = Command::new("sleep").arg("60").spawn().unwrap();
    let cancel_path = Path::new(&record.job_directory).join("cancel.requested");
    let worker = WorkerHandle {
        child: Arc::new(Mutex::new(child)),
        cancel_path: cancel_path.clone(),
        cancellation_requested: Arc::new(AtomicBool::new(false)),
    };
    state.jobs.write().unwrap().insert(
        record.job_id.clone(),
        JobEntry {
            record,
            worker: Some(worker.clone()),
            memory_limited: false,
        },
    );

    let cancelled = state.cancel("running-job").unwrap();
    assert_eq!(cancelled.status, JobStatus::Cancelling);
    assert!(worker.cancellation_requested.load(Ordering::Acquire));
    assert!(cancel_path.is_file());
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
    worker.force_kill();
    fs::remove_dir_all(temp).unwrap();
}

#[cfg(unix)]
#[test]
fn isolated_worker_signal_is_reported_without_losing_server_state() {
    let root = fixture_root();
    let executable = worker_script(&root, "crash-worker.sh", "kill -9 $$");
    let state = State::new_with_worker_executable(root.clone(), None, executable).unwrap();

    let started = state
        .start_search(StartSearchArgs {
            config_path: "tests/config.json".into(),
            approved: true,
            batch_size: Some(1),
        })
        .unwrap();
    assert_ne!(started.worker_pid, Some(std::process::id()));

    let failed = wait_for_terminal_job(&state, &started.job_id);
    assert_eq!(failed.status, JobStatus::Failed);
    assert_eq!(failed.worker_exit_code, None);
    assert!(failed
        .error
        .as_deref()
        .unwrap()
        .contains("terminated by a signal"));
    assert_eq!(state.list_jobs().len(), 1);
    assert!(state.ensure_memory_execution_available(false).is_ok());
    fs::remove_dir_all(root).unwrap();
}

#[cfg(unix)]
#[test]
fn uncooperative_worker_is_force_killed_after_cancellation_grace_period() {
    let root = fixture_root();
    let executable = worker_script(&root, "sleep-worker.sh", "exec sleep 60");
    let state = State::new_with_worker_executable(root.clone(), None, executable).unwrap();

    let started = state
        .start_search(StartSearchArgs {
            config_path: "tests/config.json".into(),
            approved: true,
            batch_size: Some(1),
        })
        .unwrap();
    let cancelling = state.cancel(&started.job_id).unwrap();
    assert_eq!(cancelling.status, JobStatus::Cancelling);

    let cancelled = wait_for_terminal_job(&state, &started.job_id);
    assert_eq!(cancelled.status, JobStatus::Cancelled);
    assert!(cancelled.error.as_deref().unwrap().contains("cancelled"));
    assert!(state.ensure_memory_execution_available(false).is_ok());
    fs::remove_dir_all(root).unwrap();
}

#[test]
fn worker_spawn_failure_is_persisted_as_a_terminal_job() {
    let root = fixture_root();
    let missing = root.join("missing-worker-executable");
    let state = State::new_with_worker_executable(root.clone(), None, missing).unwrap();

    let error = state
        .start_search(StartSearchArgs {
            config_path: "tests/config.json".into(),
            approved: true,
            batch_size: Some(1),
        })
        .unwrap_err()
        .to_string();
    assert!(error.contains("failed to start isolated Sage worker"));

    let jobs = state.list_jobs();
    assert_eq!(jobs.len(), 1);
    assert_eq!(jobs[0].status, JobStatus::Failed);
    assert!(Path::new(&jobs[0].job_directory).join("job.json").is_file());
    assert!(state.ensure_memory_execution_available(false).is_ok());
    fs::remove_dir_all(root).unwrap();
}

#[cfg(unix)]
#[test]
fn server_shutdown_marks_active_worker_interrupted_before_killing_it() {
    let root = fixture_root();
    let executable = worker_script(&root, "shutdown-worker.sh", "exec sleep 60");
    let state = State::new_with_worker_executable(root.clone(), None, executable).unwrap();
    let started = state
        .start_search(StartSearchArgs {
            config_path: "tests/config.json".into(),
            approved: true,
            batch_size: Some(1),
        })
        .unwrap();

    state.shutdown_workers();
    let interrupted = state.job(&started.job_id).unwrap();
    assert_eq!(interrupted.status, JobStatus::Interrupted);
    assert!(interrupted
        .error
        .as_deref()
        .unwrap()
        .contains("server stopped"));
    thread::sleep(Duration::from_millis(100));
    assert_eq!(
        state.job(&started.job_id).unwrap().status,
        JobStatus::Interrupted
    );
    fs::remove_dir_all(root).unwrap();
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
            worker: None,
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
            worker: None,
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
                worker: None,
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
