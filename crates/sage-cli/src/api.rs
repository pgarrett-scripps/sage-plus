use crate::events::{CancellationToken, EventEmitter, EventKind};
use crate::input::Input;
use crate::memory;
use crate::runner::{RunSummary, Runner};
use crate::telemetry::Telemetry;

/// Execution controls shared by CLI, GUI, TUI, and protocol adapters.
pub struct JobOptions {
    pub parallel: usize,
    pub events: EventEmitter,
    pub cancellation: CancellationToken,
    pub terminate_on_memory_limit: bool,
}

impl Default for JobOptions {
    fn default() -> Self {
        Self {
            parallel: (num_cpus::get() / 2).max(1),
            events: EventEmitter::disabled(),
            cancellation: CancellationToken::default(),
            terminate_on_memory_limit: false,
        }
    }
}

pub struct JobResult {
    pub summary: RunSummary,
    pub telemetry: Telemetry,
}

/// Stable application-layer entry point for validating and running Sage jobs.
pub struct SageRunner {
    input: Input,
    options: JobOptions,
}

impl SageRunner {
    pub fn new(input: Input, options: JobOptions) -> Self {
        Self { input, options }
    }

    pub fn from_path(path: impl AsRef<str>, options: JobOptions) -> anyhow::Result<Self> {
        match Input::load(path) {
            Ok(input) => Ok(Self::new(input, options)),
            Err(error) => {
                options.events.emit(EventKind::JobFailed {
                    message: error.to_string(),
                });
                Err(error)
            }
        }
    }

    pub fn validate(&self) -> anyhow::Result<()> {
        self.input.validate()?;
        self.options.events.emit(EventKind::ConfigurationValidated {
            spectra_files: self
                .input
                .mzml_paths
                .as_ref()
                .map(Vec::len)
                .unwrap_or_default(),
        });
        self.options.events.check()
    }

    pub fn run(self) -> anyhow::Result<JobResult> {
        if let Err(error) = self.validate() {
            self.options.events.emit(EventKind::JobFailed {
                message: error.to_string(),
            });
            return Err(error);
        }

        let JobOptions {
            parallel,
            events,
            cancellation,
            terminate_on_memory_limit,
        } = self.options;
        let behavior = if terminate_on_memory_limit {
            memory::MemoryLimitBehavior::TerminateProcess
        } else {
            memory::MemoryLimitBehavior::CancelJob(cancellation.clone())
        };
        let memory_guard = memory::spawn_memory_guard(self.input.memory_limits()?, behavior)?;
        let result = (|| -> anyhow::Result<JobResult> {
            cancellation.check()?;
            let search = self.input.build()?;
            let runner =
                Runner::new_with_control(search, parallel, events.clone(), cancellation.clone())?;
            let (telemetry, summary) = runner.run_with_summary(parallel)?;
            Ok(JobResult { telemetry, summary })
        })();
        if let Some(message) = memory_guard.failure() {
            events.emit(EventKind::JobFailed {
                message: message.clone(),
            });
            return Err(anyhow::anyhow!(message));
        }
        match result {
            Ok(result) => Ok(result),
            Err(error) => {
                events.emit(
                    if cancellation.is_cancelled() && !cancellation.is_memory_limit() {
                        EventKind::JobCancelled
                    } else {
                        EventKind::JobFailed {
                            message: error.to_string(),
                        }
                    },
                );
                Err(error)
            }
        }
    }
}

#[cfg(test)]
#[path = "../tests/unit/api.rs"]
mod tests;
