use super::*;

#[test]
fn converts_memory_limits() {
    let limits = MemoryLimits::from_gib(Some(8.5), Some(2.0)).unwrap();
    assert_eq!(limits.max_bytes, Some((8.5 * GIB) as u64));
    assert_eq!(limits.min_free_bytes, Some((2.0 * GIB) as u64));
    assert_eq!(limits.max_gib(), Some(8.5));
    assert_eq!(limits.min_free_gib(), Some(2.0));
}

#[test]
fn zero_and_missing_limits_are_disabled() {
    let limits = MemoryLimits::from_gib(Some(0.0), None).unwrap();
    assert_eq!(limits, MemoryLimits::default());
    assert!(!limits.is_enabled());
}

#[test]
fn rejects_invalid_limits() {
    assert!(MemoryLimits::from_gib(Some(-1.0), None).is_err());
    assert!(MemoryLimits::from_gib(Some(f64::NAN), None).is_err());
    assert!(MemoryLimits::from_gib(None, Some(f64::INFINITY)).is_err());
}

#[test]
fn detects_configured_thresholds() {
    let limits = MemoryLimits::from_gib(Some(8.0), Some(2.0)).unwrap();
    assert!(!process_limit_reached(limits, 7 * GIB as u64));
    assert!(process_limit_reached(limits, 8 * GIB as u64));
    assert!(!reserve_limit_reached(
        limits,
        3 * GIB as u64,
        16 * GIB as u64
    ));
    assert!(reserve_limit_reached(
        limits,
        2 * GIB as u64,
        16 * GIB as u64
    ));
    assert!(!reserve_limit_reached(limits, 0, 0));
}

#[test]
fn rejects_estimates_that_cross_limits() {
    let limits = MemoryLimits::from_gib(Some(8.0), Some(2.0)).unwrap();
    assert!(!estimate_exceeds_process_limit(
        limits,
        2 * GIB as u64,
        5 * GIB as u64
    ));
    assert!(estimate_exceeds_process_limit(
        limits,
        2 * GIB as u64,
        6 * GIB as u64
    ));
    assert!(!estimate_exceeds_reserve(
        limits,
        10 * GIB as u64,
        16 * GIB as u64,
        7 * GIB as u64
    ));
    assert!(estimate_exceeds_reserve(
        limits,
        9 * GIB as u64,
        16 * GIB as u64,
        7 * GIB as u64
    ));
}

#[test]
fn scoped_guard_cancels_embedded_job_at_limit() {
    let cancellation = CancellationToken::default();
    let limits = MemoryLimits {
        max_bytes: Some(1),
        min_free_bytes: None,
    };
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
fn live_memory_preflight_accepts_and_rejects_real_process_estimates() {
    let generous = MemoryLimits::from_gib(Some(1_000_000.0), None).unwrap();
    generous.check_estimate("test", 1).unwrap();

    let tiny_process_limit = MemoryLimits {
        max_bytes: Some(1),
        min_free_bytes: None,
    };
    assert!(tiny_process_limit
        .check_estimate("test", 1)
        .unwrap_err()
        .to_string()
        .contains("max_memory_gb"));

    let impossible_reserve = MemoryLimits {
        max_bytes: None,
        min_free_bytes: Some(u64::MAX),
    };
    assert!(impossible_reserve
        .check_estimate("test", 1)
        .unwrap_err()
        .to_string()
        .contains("min_free_memory_gb"));
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
