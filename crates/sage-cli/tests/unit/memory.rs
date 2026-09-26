use super::*;

#[test]
fn converts_memory_limits() {
    let limits = MemoryLimits::from_gib(Some(8.5)).unwrap();
    assert_eq!(limits.max_bytes, Some((8.5 * GIB) as u64));
    assert_eq!(limits.max_gib(), Some(8.5));
}

#[test]
fn zero_and_missing_limits_are_disabled() {
    let limits = MemoryLimits::from_gib(Some(0.0)).unwrap();
    assert_eq!(limits, MemoryLimits::default());
    assert!(!limits.is_enabled());
}

#[test]
fn rejects_invalid_limits() {
    assert!(MemoryLimits::from_gib(Some(-1.0)).is_err());
    assert!(MemoryLimits::from_gib(Some(f64::NAN)).is_err());
    assert!(MemoryLimits::from_gib(Some(f64::INFINITY)).is_err());
}

#[test]
fn detects_configured_threshold() {
    let limits = MemoryLimits::from_gib(Some(8.0)).unwrap();
    assert!(!process_limit_reached(limits, 7 * GIB as u64));
    assert!(process_limit_reached(limits, 8 * GIB as u64));
    assert!(!process_limit_reached(MemoryLimits::default(), u64::MAX));
}

#[test]
fn scoped_guard_cancels_embedded_job_at_limit() {
    let cancellation = CancellationToken::default();
    let limits = MemoryLimits { max_bytes: Some(1) };
    let guard =
        spawn_memory_guard(limits, MemoryLimitBehavior::CancelJob(cancellation.clone())).unwrap();

    for _ in 0..100 {
        if cancellation.is_cancelled() {
            break;
        }
        std::thread::sleep(std::time::Duration::from_millis(10));
    }

    assert!(cancellation.is_memory_limit());
    assert!(guard
        .failure()
        .is_some_and(|message| message.contains("configured memory limit")));
}

#[test]
fn guard_under_generous_limits_polls_and_stops_on_drop() {
    let cancellation = CancellationToken::default();
    let limits = MemoryLimits::from_gib(Some(1_000_000.0)).unwrap();
    assert!(limits.is_enabled());
    let guard =
        spawn_memory_guard(limits, MemoryLimitBehavior::CancelJob(cancellation.clone())).unwrap();
    // Let the monitor take a few samples (poll interval is 250 ms).
    std::thread::sleep(std::time::Duration::from_millis(600));
    assert!(guard.failure().is_none());
    drop(guard);
    assert!(!cancellation.is_cancelled());
}

#[test]
fn disabled_guard_does_not_cancel_job() {
    let cancellation = CancellationToken::default();
    let guard = spawn_memory_guard(
        MemoryLimits::default(),
        MemoryLimitBehavior::CancelJob(cancellation.clone()),
    )
    .unwrap();
    drop(guard);
    assert!(!cancellation.is_cancelled());
}

#[test]
fn estimate_fits_reads_the_live_process_but_never_errors() {
    assert!(MemoryLimits::default().estimate_fits(u64::MAX));
    assert!(MemoryLimits::from_gib(Some(1_000_000.0))
        .unwrap()
        .estimate_fits(1));
    assert!(!MemoryLimits { max_bytes: Some(1) }.estimate_fits(1));
}

#[test]
fn display_limit_formats_enabled_and_disabled_values() {
    assert_eq!(display_limit(None), "disabled");
    assert_eq!(display_limit(Some(GIB as u64)), "1.00 GiB");
}

#[test]
fn memory_cancellation_reason_is_distinct_from_user_cancellation() {
    let token = CancellationToken::default();
    token.cancel_for_memory_limit();
    assert!(token.is_memory_limit());
    assert!(token
        .check()
        .unwrap_err()
        .to_string()
        .contains("memory limit"));
}

#[test]
fn allocator_trim_dispatch_matches_target_support() {
    let result = trim_allocator();

    #[cfg(all(target_os = "linux", target_env = "gnu"))]
    assert_ne!(result, AllocatorTrimResult::Unsupported);

    #[cfg(not(all(target_os = "linux", target_env = "gnu")))]
    assert_eq!(result, AllocatorTrimResult::Unsupported);
}
