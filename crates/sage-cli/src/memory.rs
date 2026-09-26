//! Lightweight process memory safety limit.
//!
//! Sage performs large, highly parallel allocations. Instrumenting each one would add
//! complexity and contention to hot paths, so the limit is enforced by periodically
//! sampling the process resident set. Memory estimates never stop a run.

use crate::events::CancellationToken;
use anyhow::{ensure, Result};
use std::sync::atomic::{AtomicBool, Ordering};
use std::sync::{Arc, Mutex};
use std::thread::{self, JoinHandle};
use std::time::Duration;
use sysinfo::{ProcessExt, System, SystemExt};

const GIB: f64 = 1024.0 * 1024.0 * 1024.0;
const POLL_INTERVAL: Duration = Duration::from_millis(250);
const MEMORY_LIMIT_EXIT_CODE: i32 = 137;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
// Each target family constructs only the variants supported by its allocator.
#[allow(dead_code)]
pub(crate) enum AllocatorTrimResult {
    Released,
    NoRelease,
    Unsupported,
}

/// Ask the process allocator to return unused pages to the operating system.
///
/// GNU libc is currently the only supported allocator. Other targets use a no-op
/// because they do not expose an equivalent stable system API.
pub(crate) fn trim_allocator() -> AllocatorTrimResult {
    trim_allocator_impl()
}

#[cfg(all(target_os = "linux", target_env = "gnu"))]
fn trim_allocator_impl() -> AllocatorTrimResult {
    // SAFETY: `malloc_trim` accepts any padding value and is safe to call while
    // other threads exist. Sage calls this with zero once its database build has
    // completed and temporary build allocations have been dropped.
    if unsafe { libc::malloc_trim(0) } == 0 {
        AllocatorTrimResult::NoRelease
    } else {
        AllocatorTrimResult::Released
    }
}

#[cfg(not(all(target_os = "linux", target_env = "gnu")))]
fn trim_allocator_impl() -> AllocatorTrimResult {
    AllocatorTrimResult::Unsupported
}

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
}

impl MemoryLimits {
    pub fn from_gib(max_gib: Option<f64>) -> Result<Self> {
        Ok(Self {
            max_bytes: gib_to_bytes("max_memory_gb", max_gib)?,
        })
    }

    pub fn is_enabled(self) -> bool {
        self.max_bytes.is_some()
    }

    pub fn max_gib(self) -> Option<f64> {
        self.max_bytes.map(bytes_to_gib)
    }

    /// Whether Sage's current memory plus `additional_bytes` stays under
    /// `max_memory_gb`. Only used to choose how to proceed, never to stop a
    /// run: estimates can be far off, and the guard enforces the real limit.
    pub fn estimate_fits(self, additional_bytes: u64) -> bool {
        let Some(limit) = self.max_bytes else {
            return true;
        };
        let Ok(pid) = sysinfo::get_current_pid() else {
            return true;
        };
        let mut system = System::new();
        system.refresh_process(pid);
        system
            .process(pid)
            .map(|process| process.memory().saturating_add(additional_bytes) < limit)
            .unwrap_or(true)
    }
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
        "memory limit active: max Sage memory = {}",
        display_limit(limits.max_bytes),
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
        if process_limit_reached(limits, rss) {
            let message = format!(
                "Sage reached its configured memory limit: {:.2} GiB used, {:.2} GiB allowed. Aborting to keep the system responsive. Reduce `batch_size`, reduce database complexity, or increase `max_memory_gb`.",
                bytes_to_gib(rss),
                bytes_to_gib(limits.max_bytes.unwrap_or_default()),
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

fn display_limit(bytes: Option<u64>) -> String {
    bytes
        .map(|bytes| format!("{:.2} GiB", bytes_to_gib(bytes)))
        .unwrap_or_else(|| "disabled".into())
}

#[cfg(test)]
#[path = "../tests/unit/memory.rs"]
mod tests;
