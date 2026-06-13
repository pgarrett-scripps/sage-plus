//! Memory guard: a lightweight background watchdog that terminates the process
//! *before* it can exhaust system RAM and freeze the host.
//!
//! Proteomics searches can balloon in memory — most often during database
//! generation, where the number of modified peptide variants grows
//! combinatorially with `max_variable_mods` / `max_peff_variable_mods`, the
//! FASTA size, and the enzyme settings. Rather than instrument every allocation
//! (which would add atomic contention to Sage's hot, highly-parallel paths), we
//! poll resident memory a few times per second from a single background thread
//! and exit cleanly if either:
//!   - this process's resident set exceeds the configured ceiling, or
//!   - total system available memory drops below a small safety floor.
//!
//! Exiting deliberately (code 137, the conventional OOM-kill code) is far
//! kinder than letting the machine thrash into swap or invoke the OS OOM killer.

use std::thread;
use std::time::Duration;

use sysinfo::{ProcessExt, System, SystemExt};

const GIB: f64 = 1024.0 * 1024.0 * 1024.0;

/// How often the watchdog samples memory.
const POLL_INTERVAL: Duration = Duration::from_millis(250);

/// Exit code used when the guard trips (128 + SIGKILL(9), the OOM convention).
const OOM_EXIT_CODE: i32 = 137;

/// Resolve the resident-memory ceiling in bytes.
///
/// Precedence: explicit `user_gib` (CLI flag) > `SAGE_MAX_MEMORY_GB` env var >
/// default of 90% of total system RAM. A value of `0` (or negative) disables
/// the guard entirely and returns `None`.
pub fn resolve_limit_bytes(user_gib: Option<f64>) -> Option<u64> {
    let from_gib = |g: f64| (g * GIB) as u64;

    if let Some(g) = user_gib {
        return (g > 0.0).then(|| from_gib(g));
    }

    if let Ok(raw) = std::env::var("SAGE_MAX_MEMORY_GB") {
        if let Ok(g) = raw.trim().parse::<f64>() {
            return (g > 0.0).then(|| from_gib(g));
        }
    }

    let mut sys = System::new();
    sys.refresh_memory();
    let total = sys.total_memory();
    (total > 0).then(|| total / 10 * 9)
}

/// Spawn the background memory watchdog. `limit_bytes` is this process's
/// resident-set ceiling; the watchdog also independently aborts if system-wide
/// available memory falls below a safety floor.
pub fn spawn_memory_guard(limit_bytes: u64) {
    let result = thread::Builder::new()
        .name("sage-memory-guard".into())
        .spawn(move || guard_loop(limit_bytes));

    if result.is_err() {
        log::warn!("could not spawn memory guard thread; continuing without memory protection");
    }
}

fn guard_loop(limit_bytes: u64) {
    let pid = match sysinfo::get_current_pid() {
        Ok(pid) => pid,
        Err(e) => {
            log::warn!("memory guard disabled: could not determine current pid ({e})");
            return;
        }
    };

    let mut sys = System::new();
    sys.refresh_memory();
    // Safety floor: never let total system available memory drop below this, so
    // other processes (and the OS) keep breathing room. max(1 GiB, 2% of RAM).
    let total = sys.total_memory();
    let available_floor = (total / 50).max(GIB as u64);

    log::info!(
        "memory guard active: aborting if Sage exceeds {:.1} GiB resident, \
         or if system available memory falls below {:.1} GiB",
        limit_bytes as f64 / GIB,
        available_floor as f64 / GIB,
    );

    loop {
        sys.refresh_memory();
        sys.refresh_process(pid);

        let rss = sys.process(pid).map(|p| p.memory()).unwrap_or(0);
        let available = sys.available_memory();

        if rss > limit_bytes {
            log::error!(
                "Sage is using {:.1} GiB of memory, which exceeds the {:.1} GiB limit. \
                 Aborting now to keep your system responsive.\n\
                 To reduce memory: lower `database.max_variable_mods` / \
                 `database.max_peff_variable_mods`, use a smaller FASTA, narrow your \
                 precursor/fragment tolerances, or enable `database.prefilter`.\n\
                 To allow more memory: pass `--max-memory <GiB>` (or set \
                 `SAGE_MAX_MEMORY_GB`); `--max-memory 0` disables this guard.",
                rss as f64 / GIB,
                limit_bytes as f64 / GIB,
            );
            std::process::exit(OOM_EXIT_CODE);
        }

        if available != 0 && available < available_floor {
            log::error!(
                "System available memory has dropped to {:.1} GiB (Sage is using {:.1} GiB). \
                 Aborting now to prevent your machine from freezing.\n\
                 Reduce the search size (see `max_variable_mods` / `max_peff_variable_mods` / \
                 `prefilter`), or pass `--max-memory 0` to disable this guard.",
                available as f64 / GIB,
                rss as f64 / GIB,
            );
            std::process::exit(OOM_EXIT_CODE);
        }

        thread::sleep(POLL_INTERVAL);
    }
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn explicit_limit_wins_and_zero_disables() {
        assert_eq!(resolve_limit_bytes(Some(4.0)), Some((4.0 * GIB) as u64));
        assert_eq!(resolve_limit_bytes(Some(0.0)), None);
        assert_eq!(resolve_limit_bytes(Some(-1.0)), None);
    }

    #[test]
    fn default_is_ninety_percent_of_total() {
        // With no flag and no env override, the limit is a positive fraction of
        // total RAM (and strictly less than total).
        if std::env::var("SAGE_MAX_MEMORY_GB").is_err() {
            let mut sys = System::new();
            sys.refresh_memory();
            let total = sys.total_memory();
            if total > 0 {
                let limit = resolve_limit_bytes(None).unwrap();
                assert!(limit > 0 && limit < total, "limit={limit} total={total}");
            }
        }
    }
}
