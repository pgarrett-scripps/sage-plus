use super::*;
use std::io::Write;
use std::sync::{Arc, Mutex};

#[derive(Clone, Default)]
struct SharedWriter(Arc<Mutex<Vec<u8>>>);

impl Write for SharedWriter {
    fn write(&mut self, buffer: &[u8]) -> std::io::Result<usize> {
        self.0.lock().unwrap().extend_from_slice(buffer);
        Ok(buffer.len())
    }

    fn flush(&mut self) -> std::io::Result<()> {
        Ok(())
    }
}

#[test]
fn default_job_options_are_runnable_and_quiet() {
    let options = JobOptions::default();

    assert!(options.parallel >= 1);
    assert!(!options.events.is_enabled());
    assert!(!options.cancellation.is_cancelled());
    assert!(!options.terminate_on_memory_limit);
}

#[test]
fn loading_failure_emits_a_structured_job_failure() {
    let missing = std::env::temp_dir().join(format!(
        "sage-api-missing-{}-{}.json",
        std::process::id(),
        std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .unwrap()
            .as_nanos()
    ));
    let writer = SharedWriter::default();
    let output = writer.0.clone();
    let options = JobOptions {
        events: EventEmitter::from_writer(writer),
        ..Default::default()
    };

    assert!(SageRunner::from_path(missing.to_string_lossy(), options).is_err());

    let lines = String::from_utf8(output.lock().unwrap().clone()).unwrap();
    let event: serde_json::Value = serde_json::from_str(lines.trim()).unwrap();
    assert_eq!(event["event"], "job_failed");
    assert!(event["message"]
        .as_str()
        .is_some_and(|message| !message.is_empty()));
}
