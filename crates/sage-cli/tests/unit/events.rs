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

#[test]
fn cancellation_token_is_shared_across_clones() {
    let token = CancellationToken::default();
    let clone = token.clone();

    assert!(token.check().is_ok());
    clone.cancel();
    assert!(token.is_cancelled());
    assert!(token.check().unwrap_err().to_string().contains("cancelled"));
}

struct FailingWriter;

impl Write for FailingWriter {
    fn write(&mut self, _buffer: &[u8]) -> std::io::Result<usize> {
        Err(std::io::Error::other("intentional test failure"))
    }

    fn flush(&mut self) -> std::io::Result<()> {
        Ok(())
    }
}

#[test]
fn writer_failures_are_reported_by_check() {
    let emitter = EventEmitter::from_writer(FailingWriter);

    emitter.emit(EventKind::DatabaseStarted);

    assert!(emitter
        .check()
        .unwrap_err()
        .to_string()
        .contains("intentional test failure"));
}
