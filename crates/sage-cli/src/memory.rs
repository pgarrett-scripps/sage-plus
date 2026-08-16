//! Lightweight process and system memory safety limits.
//!
//! Sage performs large, highly parallel allocations. Instrumenting each one would add
//! complexity and contention to hot paths, so these limits are enforced by periodically
//! sampling the process resident set and the memory available to the whole system.

use anyhow::{bail, ensure, Context, Result};
use std::thread;
use std::time::Duration;
use sysinfo::{ProcessExt, System, SystemExt};

const GIB: f64 = 1024.0 * 1024.0 * 1024.0;
const POLL_INTERVAL: Duration = Duration::from_millis(250);
const MEMORY_LIMIT_EXIT_CODE: i32 = 137;

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct MemoryLimits {
    max_bytes: Option<u64>,
    min_free_bytes: Option<u64>,
}

impl MemoryLimits {
    pub fn from_gib(max_gib: Option<f64>, min_free_gib: Option<f64>) -> Result<Self> {
        Ok(Self {
            max_bytes: gib_to_bytes("max_memory_gb", max_gib)?,
            min_free_bytes: gib_to_bytes("min_free_memory_gb", min_free_gib)?,
        })
    }

    pub fn is_enabled(self) -> bool {
        self.max_bytes.is_some() || self.min_free_bytes.is_some()
    }

    pub fn max_gib(self) -> Option<f64> {
        self.max_bytes.map(bytes_to_gib)
    }

    pub fn min_free_gib(self) -> Option<f64> {
        self.min_free_bytes.map(bytes_to_gib)
    }

    /// Reject an estimated allocation before it begins if it would cross either limit.
    pub fn check_estimate(self, stage: &str, additional_bytes: u64) -> Result<()> {
        if !self.is_enabled() {
            return Ok(());
        }

        let pid = sysinfo::get_current_pid()
            .map_err(|error| anyhow::anyhow!("could not determine the Sage process ID: {error}"))?;
        let mut system = System::new();
        system.refresh_memory();
        system.refresh_process(pid);
        let rss = system
            .process(pid)
            .context("could not inspect the Sage process")?
            .memory();

        if estimate_exceeds_process_limit(self, rss, additional_bytes) {
            bail!(
                "estimated {stage} database peak is {:.2} GiB in addition to Sage's current {:.2} GiB; this would exceed `max_memory_gb` ({:.2} GiB). Reduce variable modifications, `max_variable_mods`, or database size",
                bytes_to_gib(additional_bytes),
                bytes_to_gib(rss),
                bytes_to_gib(self.max_bytes.unwrap_or_default()),
            );
        }

        let available = system.available_memory();
        if estimate_exceeds_reserve(self, available, system.total_memory(), additional_bytes) {
            bail!(
                "estimated {stage} database peak requires {:.2} GiB, but only {:.2} GiB is currently available while preserving `min_free_memory_gb` ({:.2} GiB). Reduce variable modifications, `max_variable_mods`, or database size",
                bytes_to_gib(additional_bytes),
                bytes_to_gib(available),
                bytes_to_gib(self.min_free_bytes.unwrap_or_default()),
            );
        }

        Ok(())
    }
}

fn estimate_exceeds_process_limit(limits: MemoryLimits, rss: u64, additional_bytes: u64) -> bool {
    limits
        .max_bytes
        .map(|limit| rss.saturating_add(additional_bytes) >= limit)
        .unwrap_or(false)
}

fn estimate_exceeds_reserve(
    limits: MemoryLimits,
    available: u64,
    total: u64,
    additional_bytes: u64,
) -> bool {
    total != 0
        && limits
            .min_free_bytes
            .map(|reserve| additional_bytes.saturating_add(reserve) >= available)
            .unwrap_or(false)
}

fn gib_to_bytes(name: &str, value: Option<f64>) -> Result<Option<u64>> {
    let value = match value {
        Some(value) => value,
        None => return Ok(None),
    };

    ensure!(value.is_finite(), "`{name}` must be a finite number");
    ensure!(value >= 0.0, "`{name}` must not be negative");
    if value == 0.0 {
        return Ok(None);
    }

    ensure!(value <= u64::MAX as f64 / GIB, "`{name}` is too large");
    Ok(Some((value * GIB) as u64))
}

fn bytes_to_gib(bytes: u64) -> f64 {
    bytes as f64 / GIB
}

/// Start a detached monitor before Sage performs its large allocations.
pub fn spawn_memory_guard(limits: MemoryLimits) -> std::io::Result<()> {
    if !limits.is_enabled() {
        return Ok(());
    }

    thread::Builder::new()
        .name("sage-memory-guard".into())
        .spawn(move || guard_loop(limits))?;
    Ok(())
}

fn guard_loop(limits: MemoryLimits) {
    let pid = match sysinfo::get_current_pid() {
        Ok(pid) => pid,
        Err(error) => {
            log::error!("memory guard could not determine the Sage process ID: {error}");
            std::process::exit(MEMORY_LIMIT_EXIT_CODE);
        }
    };

    log::info!(
        "memory limits active: max Sage memory = {}, minimum free system memory = {}",
        display_limit(limits.max_bytes),
        display_limit(limits.min_free_bytes),
    );

    let mut system = System::new();
    loop {
        system.refresh_memory();
        system.refresh_process(pid);

        let rss = match system.process(pid) {
            Some(process) => process.memory(),
            None => {
                log::error!("memory guard could not inspect the Sage process");
                std::process::exit(MEMORY_LIMIT_EXIT_CODE);
            }
        };
        let available = system.available_memory();

        if process_limit_reached(limits, rss) {
            log::error!(
                "Sage reached its configured memory limit: {:.2} GiB used, {:.2} GiB allowed. Aborting to keep the system responsive. Reduce `batch_size`, reduce database complexity, or increase `max_memory_gb`.",
                bytes_to_gib(rss),
                bytes_to_gib(limits.max_bytes.unwrap_or_default()),
            );
            std::process::exit(MEMORY_LIMIT_EXIT_CODE);
        }

        if reserve_limit_reached(limits, available, system.total_memory()) {
            log::error!(
                "System available memory reached Sage's configured reserve: {:.2} GiB available, {:.2} GiB required. Sage is using {:.2} GiB. Aborting to keep the system responsive.",
                bytes_to_gib(available),
                bytes_to_gib(limits.min_free_bytes.unwrap_or_default()),
                bytes_to_gib(rss),
            );
            std::process::exit(MEMORY_LIMIT_EXIT_CODE);
        }

        thread::sleep(POLL_INTERVAL);
    }
}

fn process_limit_reached(limits: MemoryLimits, rss: u64) -> bool {
    limits.max_bytes.map(|limit| rss >= limit).unwrap_or(false)
}

fn reserve_limit_reached(limits: MemoryLimits, available: u64, total: u64) -> bool {
    total != 0
        && limits
            .min_free_bytes
            .map(|minimum| available <= minimum)
            .unwrap_or(false)
}

fn display_limit(bytes: Option<u64>) -> String {
    bytes
        .map(|bytes| format!("{:.2} GiB", bytes_to_gib(bytes)))
        .unwrap_or_else(|| "disabled".into())
}

#[cfg(test)]
mod tests {
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
}
