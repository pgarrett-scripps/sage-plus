use crate::events::{CancellationToken, EventEmitter, EventKind};
use crate::input::Input;
use crate::memory;
use crate::runner::{RunSummary, Runner};
use crate::telemetry::Telemetry;

/// Execution controls shared by CLI, GUI, TUI, and protocol adapters.
pub struct JobOptions {
    pub parallel: usize,
    pub parquet: bool,
    pub events: EventEmitter,
    pub cancellation: CancellationToken,
}

impl Default for JobOptions {
    fn default() -> Self {
        Self {
            parallel: (num_cpus::get() / 2).max(1),
            parquet: false,
            events: EventEmitter::disabled(),
            cancellation: CancellationToken::default(),
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
            parquet,
            events,
            cancellation,
        } = self.options;
        memory::spawn_memory_guard(self.input.memory_limits()?)?;
        if let Err(error) = cancellation.check() {
            events.emit(EventKind::JobCancelled);
            return Err(error);
        }
        let search = match self.input.build() {
            Ok(search) => search,
            Err(error) => {
                events.emit(EventKind::JobFailed {
                    message: error.to_string(),
                });
                return Err(error);
            }
        };
        let runner = match Runner::new_with_control(
            search,
            parallel,
            events.clone(),
            cancellation.clone(),
        ) {
            Ok(runner) => runner,
            Err(error) => {
                events.emit(if cancellation.is_cancelled() {
                    EventKind::JobCancelled
                } else {
                    EventKind::JobFailed {
                        message: error.to_string(),
                    }
                });
                return Err(error);
            }
        };
        match runner.run_with_summary(parallel, parquet) {
            Ok((telemetry, summary)) => Ok(JobResult { telemetry, summary }),
            Err(error) => {
                events.emit(if cancellation.is_cancelled() {
                    EventKind::JobCancelled
                } else {
                    EventKind::JobFailed {
                        message: error.to_string(),
                    }
                });
                Err(error)
            }
        }
    }
}
