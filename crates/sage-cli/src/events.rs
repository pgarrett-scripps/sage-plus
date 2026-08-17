use serde::Serialize;
use std::io::Write;
use std::sync::atomic::{AtomicBool, AtomicU64, Ordering};
use std::sync::{Arc, Mutex};
use std::time::Instant;

/// A cooperative cancellation handle for a running Sage job.
#[derive(Clone, Default)]
pub struct CancellationToken(Arc<AtomicBool>);

impl CancellationToken {
    pub fn cancel(&self) {
        self.0.store(true, Ordering::Release);
    }

    pub fn is_cancelled(&self) -> bool {
        self.0.load(Ordering::Acquire)
    }

    pub(crate) fn check(&self) -> anyhow::Result<()> {
        anyhow::ensure!(!self.is_cancelled(), "Sage job cancelled");
        Ok(())
    }
}

/// Stable, machine-readable events emitted during a Sage job.
#[derive(Debug, Clone, Serialize)]
#[serde(tag = "event", rename_all = "snake_case")]
pub enum EventKind {
    ConfigurationValidated {
        spectra_files: usize,
    },
    DatabaseStarted,
    DatabaseBuilt {
        peptides: usize,
        fragments: usize,
    },
    FileStarted {
        file_id: usize,
        path: String,
    },
    FileCompleted {
        file_id: usize,
        path: String,
        spectra: usize,
    },
    SpectraProcessed {
        ms1_spectra: usize,
        msn_spectra: usize,
    },
    SearchProgress {
        files_completed: usize,
        files_total: usize,
    },
    RtModelFitted,
    RtModelSkipped {
        reason: String,
    },
    MobilityModelFitted,
    MobilityModelSkipped {
        reason: String,
    },
    FdrCompleted {
        psms: usize,
        peptides: usize,
        proteins: usize,
        protein_groups: usize,
    },
    OutputWritten {
        path: String,
    },
    JobCompleted {
        runtime_secs: u64,
        outputs: usize,
    },
    Warning {
        code: String,
        message: String,
    },
    JobCancelled,
    JobFailed {
        message: String,
    },
}

#[derive(Debug, Serialize)]
struct EventEnvelope<'a> {
    schema_version: u8,
    sequence: u64,
    elapsed_ms: u128,
    #[serde(flatten)]
    kind: &'a EventKind,
}

struct EventWriter {
    writer: Mutex<Box<dyn Write + Send>>,
    error: Mutex<Option<String>>,
}

/// Cloneable event destination that serializes one JSON object per line.
#[derive(Clone)]
pub struct EventEmitter {
    writer: Option<Arc<EventWriter>>,
    sequence: Arc<AtomicU64>,
    started: Instant,
}

impl Default for EventEmitter {
    fn default() -> Self {
        Self::disabled()
    }
}

impl EventEmitter {
    pub fn disabled() -> Self {
        Self {
            writer: None,
            sequence: Arc::new(AtomicU64::new(0)),
            started: Instant::now(),
        }
    }

    pub fn from_writer(writer: impl Write + Send + 'static) -> Self {
        Self {
            writer: Some(Arc::new(EventWriter {
                writer: Mutex::new(Box::new(writer)),
                error: Mutex::new(None),
            })),
            sequence: Arc::new(AtomicU64::new(0)),
            started: Instant::now(),
        }
    }

    pub fn is_enabled(&self) -> bool {
        self.writer.is_some()
    }

    pub fn emit(&self, kind: EventKind) {
        let destination = match &self.writer {
            Some(destination) => destination,
            None => return,
        };
        let envelope = EventEnvelope {
            schema_version: 1,
            sequence: self.sequence.fetch_add(1, Ordering::Relaxed),
            elapsed_ms: self.started.elapsed().as_millis(),
            kind: &kind,
        };
        let result = (|| -> anyhow::Result<()> {
            let mut writer = destination
                .writer
                .lock()
                .map_err(|_| anyhow::anyhow!("event writer lock poisoned"))?;
            serde_json::to_writer(&mut *writer, &envelope)?;
            writer.write_all(b"\n")?;
            writer.flush()?;
            Ok(())
        })();
        if let Err(error) = result {
            if let Ok(mut stored) = destination.error.lock() {
                if stored.is_none() {
                    *stored = Some(error.to_string());
                }
            }
        }
    }

    pub(crate) fn check(&self) -> anyhow::Result<()> {
        if let Some(destination) = &self.writer {
            let error = destination
                .error
                .lock()
                .map_err(|_| anyhow::anyhow!("event error lock poisoned"))?
                .clone();
            if let Some(error) = error {
                anyhow::bail!("failed to write structured event: {error}");
            }
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::sync::{Arc, Mutex};

    #[derive(Clone, Default)]
    struct SharedWriter(Arc<Mutex<Vec<u8>>>);

    impl Write for SharedWriter {
        fn write(&mut self, buf: &[u8]) -> std::io::Result<usize> {
            self.0.lock().unwrap().extend_from_slice(buf);
            Ok(buf.len())
        }

        fn flush(&mut self) -> std::io::Result<()> {
            Ok(())
        }
    }

    #[test]
    fn emits_versioned_json_lines() {
        let writer = SharedWriter::default();
        let output = writer.0.clone();
        let emitter = EventEmitter::from_writer(writer);
        emitter.emit(EventKind::ConfigurationValidated { spectra_files: 2 });
        emitter.emit(EventKind::DatabaseStarted);
        emitter.check().unwrap();

        let bytes = output.lock().unwrap().clone();
        let lines = String::from_utf8(bytes).unwrap();
        let values = lines
            .lines()
            .map(|line| serde_json::from_str::<serde_json::Value>(line).unwrap())
            .collect::<Vec<_>>();
        assert_eq!(values.len(), 2);
        assert_eq!(values[0]["schema_version"], 1);
        assert_eq!(values[0]["sequence"], 0);
        assert_eq!(values[0]["event"], "configuration_validated");
        assert_eq!(values[1]["sequence"], 1);
    }
}
