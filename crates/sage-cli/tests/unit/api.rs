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

#[test]
fn configured_batching_overrides_legacy_api_parallelism() -> anyhow::Result<()> {
    let root = std::env::temp_dir().join(format!(
        "sage-api-batches-{}-{}",
        std::process::id(),
        std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)?
            .as_nanos()
    ));
    let workspace = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("../..");
    let mut config: serde_json::Value =
        serde_json::from_slice(&std::fs::read(workspace.join("tests/config.json"))?)?;
    config["database"]["fasta"] = workspace
        .join("tests/Q99536.fasta")
        .to_string_lossy()
        .into_owned()
        .into();
    config["mzml_paths"] =
        serde_json::json!(vec![workspace.join("tests/LQSRPAAPPAPGPGQLTLR.mzML"); 3]);
    let mut result_bytes = None;
    for batch in [1, 2] {
        config["batch_size"] = batch.into();
        config["output_directory"] = root
            .join(batch.to_string())
            .to_string_lossy()
            .into_owned()
            .into();
        let writer = SharedWriter::default();
        let output = writer.0.clone();
        let result = SageRunner::new(
            serde_json::from_value(config.clone())?,
            JobOptions {
                parallel: 3,
                events: EventEmitter::from_writer(writer),
                ..Default::default()
            },
        )
        .run()?;
        let events: Vec<serde_json::Value> = String::from_utf8(output.lock().unwrap().clone())?
            .lines()
            .map(serde_json::from_str)
            .collect::<Result<_, _>>()?;
        let progress: Vec<_> = events
            .iter()
            .filter(|e| e["event"] == "search_progress")
            .map(|e| e["files_completed"].as_u64().unwrap())
            .collect();
        assert_eq!(
            progress,
            if batch == 1 {
                vec![1, 2, 3]
            } else {
                vec![2, 3]
            }
        );
        assert_eq!(result.summary.execution.batch_size, batch as usize);
        assert_eq!(
            result.summary.execution.rayon_threads,
            rayon::current_num_threads()
        );
        assert!(!result.summary.models.mass_alignment_applied);
        assert_eq!(result.summary.models.mass_alignment_files.len(), 3);
        assert!(result
            .summary
            .models
            .mass_alignment_files
            .iter()
            .all(|file| file.precursor_skip_reason.is_some()));
        assert_eq!(
            result.summary.provenance.input_identity_mode,
            "path_size_mtime"
        );
        assert!(result
            .summary
            .provenance
            .inputs
            .iter()
            .all(|input| input.size_bytes.is_some()));
        let bytes = std::fs::read(root.join(batch.to_string()).join("results.sage.parquet"))?;
        if let Some(previous) = &result_bytes {
            assert_eq!(&bytes, previous);
        }
        result_bytes = Some(bytes);
    }
    std::fs::remove_dir_all(root)?;
    Ok(())
}

struct CancellingWriter {
    output: SharedWriter,
    cancellation: CancellationToken,
    trigger: &'static str,
}

impl Write for CancellingWriter {
    fn write(&mut self, buffer: &[u8]) -> std::io::Result<usize> {
        self.output.write(buffer)
    }

    fn flush(&mut self) -> std::io::Result<()> {
        let bytes = self.output.0.lock().unwrap();
        if bytes
            .windows(self.trigger.len())
            .any(|window| window == self.trigger.as_bytes())
        {
            self.cancellation.cancel();
        }
        Ok(())
    }
}

#[test]
fn cancellation_during_search_or_annotation_never_completes_the_run() -> anyhow::Result<()> {
    let workspace = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("../..");
    let root = std::env::temp_dir().join(format!(
        "sage-cancel-{}-{}",
        std::process::id(),
        std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)?
            .as_nanos()
    ));
    for trigger in ["file_started", "fragment_annotation_completed"] {
        let mut config: serde_json::Value =
            serde_json::from_slice(&std::fs::read(workspace.join("tests/config.json"))?)?;
        config["database"]["fasta"] = workspace
            .join("tests/Q99536.fasta")
            .to_string_lossy()
            .into_owned()
            .into();
        config["mzml_paths"] =
            serde_json::json!([workspace.join("tests/LQSRPAAPPAPGPGQLTLR.mzML")]);
        config["output_directory"] = root.join(trigger).to_string_lossy().into_owned().into();
        let output = SharedWriter::default();
        let cancellation = CancellationToken::default();
        let events = EventEmitter::from_writer(CancellingWriter {
            output: output.clone(),
            cancellation: cancellation.clone(),
            trigger,
        });
        let result = SageRunner::new(
            serde_json::from_value(config)?,
            JobOptions {
                events,
                cancellation,
                ..Default::default()
            },
        )
        .run();
        assert!(result.is_err());
        assert!(!root.join(trigger).join("run-summary.json").exists());
        let bytes = output.0.lock().unwrap();
        let text = std::str::from_utf8(&bytes)?;
        assert!(text.contains("job_cancelled"));
        assert!(!text.contains("job_completed"));
    }
    std::fs::remove_dir_all(root)?;
    Ok(())
}
