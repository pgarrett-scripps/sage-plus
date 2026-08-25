//! Lightweight process and system memory safety limits.
//!
//! Sage performs large, highly parallel allocations. Instrumenting each one would add
//! complexity and contention to hot paths, so these limits are enforced by periodically
//! sampling the process resident set and the memory available to the whole system.

use crate::events::CancellationToken;
use anyhow::{bail, ensure, Context, Result};
use std::sync::atomic::{AtomicBool, Ordering};
use std::sync::{Arc, Mutex};
use std::thread::{self, JoinHandle};
use std::time::Duration;
use sysinfo::{ProcessExt, System, SystemExt};

const GIB: f64 = 1024.0 * 1024.0 * 1024.0;
const POLL_INTERVAL: Duration = Duration::from_millis(250);
const MEMORY_LIMIT_EXIT_CODE: i32 = 137;

#[derive(Clone)]
pub enum MemoryLimitBehavior {
    TerminateProcess,
    CancelJob(CancellationToken),
}

pub struct MemoryGuard {
    stop: Arc<AtomicBool>,
    failure: Arc<Mutex<Option<String>>>,
    thread: Option<JoinHandle<()>>,
}

impl MemoryGuard {
    pub fn failure(&self) -> Option<String> {
        self.failure.lock().ok().and_then(|failure| failure.clone())
    }
}

impl Drop for MemoryGuard {
    fn drop(&mut self) {
        self.stop.store(true, Ordering::Release);
        if let Some(thread) = self.thread.take() {
            thread.thread().unpark();
            if thread.join().is_err() {
                log::error!("memory guard thread panicked while stopping");
            }
        }
    }
}

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

/// Start a monitor before Sage performs its large allocations.
pub fn spawn_memory_guard(
    limits: MemoryLimits,
    behavior: MemoryLimitBehavior,
) -> std::io::Result<MemoryGuard> {
    let stop = Arc::new(AtomicBool::new(false));
    let failure = Arc::new(Mutex::new(None));
    if !limits.is_enabled() {
        return Ok(MemoryGuard {
            stop,
            failure,
            thread: None,
        });
    }

    let guard_stop = stop.clone();
    let guard_failure = failure.clone();
    let thread = thread::Builder::new()
        .name("sage-memory-guard".into())
        .spawn(move || guard_loop(limits, behavior, guard_stop, guard_failure))?;
    Ok(MemoryGuard {
        stop,
        failure,
        thread: Some(thread),
    })
}

fn guard_loop(
    limits: MemoryLimits,
    behavior: MemoryLimitBehavior,
    stop: Arc<AtomicBool>,
    failure: Arc<Mutex<Option<String>>>,
) {
    let pid = match sysinfo::get_current_pid() {
        Ok(pid) => pid,
        Err(error) => {
            trigger(
                &behavior,
                &failure,
                format!("memory guard could not determine the Sage process ID: {error}"),
            );
            return;
        }
    };

    log::info!(
        "memory limits active: max Sage memory = {}, minimum free system memory = {}",
        display_limit(limits.max_bytes),
        display_limit(limits.min_free_bytes),
    );

    let mut system = System::new();
    while !stop.load(Ordering::Acquire) {
        system.refresh_memory();
        system.refresh_process(pid);

        let rss = match system.process(pid) {
            Some(process) => process.memory(),
            None => {
                trigger(
                    &behavior,
                    &failure,
                    "memory guard could not inspect the Sage process".into(),
                );
                return;
            }
        };
        let available = system.available_memory();

        if process_limit_reached(limits, rss) {
            let message = format!(
                "Sage reached its configured memory limit: {:.2} GiB used, {:.2} GiB allowed. Aborting to keep the system responsive. Reduce `batch_size`, reduce database complexity, or increase `max_memory_gb`.",
                bytes_to_gib(rss),
                bytes_to_gib(limits.max_bytes.unwrap_or_default()),
            );
            trigger(&behavior, &failure, message);
            return;
        }

        if reserve_limit_reached(limits, available, system.total_memory()) {
            let message = format!(
                "System available memory reached Sage's configured reserve: {:.2} GiB available, {:.2} GiB required. Sage is using {:.2} GiB. Aborting to keep the system responsive.",
                bytes_to_gib(available),
                bytes_to_gib(limits.min_free_bytes.unwrap_or_default()),
                bytes_to_gib(rss),
            );
            trigger(&behavior, &failure, message);
            return;
        }

        thread::park_timeout(POLL_INTERVAL);
    }
}

fn trigger(behavior: &MemoryLimitBehavior, failure: &Mutex<Option<String>>, message: String) {
    log::error!("{message}");
    match behavior {
        MemoryLimitBehavior::TerminateProcess => std::process::exit(MEMORY_LIMIT_EXIT_CODE),
        MemoryLimitBehavior::CancelJob(cancellation) => {
            if let Ok(mut stored) = failure.lock() {
                *stored = Some(message);
            }
            cancellation.cancel_for_memory_limit();
        }
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
#[path = "../tests/unit/memory.rs"]
mod tests;
