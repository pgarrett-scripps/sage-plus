use anyhow::{ensure, Context};
use rmcp::{
    handler::server::wrapper::Parameters,
    model::{
        CallToolResult, ListResourcesResult, PaginatedRequestParams, ReadResourceRequestParams,
        ReadResourceResponse, ReadResourceResult, Resource, ResourceContents, ServerCapabilities,
        ServerInfo,
    },
    service::RequestContext,
    tool, tool_handler, tool_router, ErrorData as McpError, RoleServer, ServerHandler,
};
use sage_cli::{
    api::{JobOptions, SageRunner},
    events::{CancellationToken, EventEmitter},
    input::Input,
    runner::RunSummary,
};
use schemars::JsonSchema;
use serde::{Deserialize, Serialize};
use std::{
    collections::HashMap,
    fs::{self, File},
    io::{BufRead, BufReader, BufWriter},
    path::{Path, PathBuf},
    sync::{Arc, RwLock},
    time::{SystemTime, UNIX_EPOCH},
};
use uuid::Uuid;

#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq)]
#[serde(rename_all = "snake_case")]
pub enum JobStatus {
    Queued,
    Running,
    Cancelling,
    Completed,
    Failed,
    Cancelled,
    Interrupted,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct JobRecord {
    pub job_id: String,
    pub status: JobStatus,
    pub config_path: String,
    pub job_directory: String,
    pub events_path: String,
    pub output_directory: String,
    pub created_at_unix: u64,
    pub updated_at_unix: u64,
    pub summary: Option<RunSummary>,
    pub error: Option<String>,
}

#[derive(Clone)]
struct JobEntry {
    record: JobRecord,
    cancellation: Option<CancellationToken>,
}

#[derive(Clone)]
struct State {
    root: PathBuf,
    jobs_dir: PathBuf,
    jobs: Arc<RwLock<HashMap<String, JobEntry>>>,
}

#[derive(Debug, Serialize)]
struct ConfigValidation {
    valid: bool,
    config_path: String,
    fasta_path: Option<String>,
    peptide_database_path: Option<String>,
    spectra_paths: Vec<String>,
    note: &'static str,
}

impl State {
    fn new(root: PathBuf, jobs_dir: Option<PathBuf>) -> anyhow::Result<Self> {
        let root = root
            .canonicalize()
            .with_context(|| format!("cannot access MCP root `{}`", root.display()))?;
        ensure!(root.is_dir(), "MCP root must be a directory");

        let requested_jobs = jobs_dir.unwrap_or_else(|| root.join(".sage/jobs"));
        let requested_jobs = if requested_jobs.is_absolute() {
            requested_jobs
        } else {
            root.join(requested_jobs)
        };
        fs::create_dir_all(&requested_jobs)?;
        let jobs_dir = requested_jobs.canonicalize()?;
        ensure!(
            jobs_dir.starts_with(&root),
            "jobs directory must be contained by the MCP root"
        );

        let state = Self {
            root,
            jobs_dir,
            jobs: Arc::new(RwLock::new(HashMap::new())),
        };
        state.restore_jobs()?;
        Ok(state)
    }

    fn restore_jobs(&self) -> anyhow::Result<()> {
        let mut restored = HashMap::new();
        for entry in fs::read_dir(&self.jobs_dir)? {
            let job_dir = entry?.path().canonicalize()?;
            if !job_dir.starts_with(&self.jobs_dir) {
                continue;
            }
            let manifest = job_dir.join("job.json");
            if !manifest.is_file() {
                continue;
            }
            let bytes = fs::read(&manifest)?;
            let mut record: JobRecord = serde_json::from_slice(&bytes)
                .with_context(|| format!("invalid job manifest `{}`", manifest.display()))?;
            record.job_directory = job_dir.to_string_lossy().into_owned();
            record.events_path = job_dir.join("events.jsonl").to_string_lossy().into_owned();
            record.output_directory = job_dir.join("output").to_string_lossy().into_owned();
            if matches!(
                record.status,
                JobStatus::Queued | JobStatus::Running | JobStatus::Cancelling
            ) {
                record.status = JobStatus::Interrupted;
                record.updated_at_unix = now();
                write_record(&record)?;
            }
            restored.insert(
                record.job_id.clone(),
                JobEntry {
                    record,
                    cancellation: None,
                },
            );
        }
        *self.jobs.write().expect("job registry poisoned") = restored;
        Ok(())
    }

    fn resolve_existing(&self, value: &str) -> anyhow::Result<PathBuf> {
        ensure!(
            !value.contains("://"),
            "remote URLs are disabled by sage-mcp"
        );
        let requested = Path::new(value);
        let path = if requested.is_absolute() {
            requested.to_path_buf()
        } else {
            self.root.join(requested)
        };
        let path = path
            .canonicalize()
            .with_context(|| format!("cannot access `{}`", path.display()))?;
        ensure!(
            path.starts_with(&self.root),
            "path `{}` is outside the configured MCP root",
            path.display()
        );
        Ok(path)
    }

    fn load_config(&self, config_path: &str) -> anyhow::Result<(PathBuf, Input)> {
        let config_path = self.resolve_existing(config_path)?;
        ensure!(config_path.is_file(), "configuration path must be a file");
        let mut input = Input::load(config_path.to_string_lossy())?;
        input.validate()?;

        if let Some(fasta) = input.database.fasta.as_mut() {
            *fasta = self.resolve_existing(fasta)?.to_string_lossy().into_owned();
        }
        if let Some(peptides) = input.database.peptides.as_mut() {
            *peptides = self
                .resolve_existing(peptides)?
                .to_string_lossy()
                .into_owned();
        }
        for spectra in input.mzml_paths.iter_mut().flatten() {
            *spectra = self
                .resolve_existing(spectra)?
                .to_string_lossy()
                .into_owned();
        }
        Ok((config_path, input))
    }

    fn validate_config(&self, config_path: &str) -> anyhow::Result<ConfigValidation> {
        let (config_path, input) = self.load_config(config_path)?;
        Ok(ConfigValidation {
            valid: true,
            config_path: config_path.to_string_lossy().into_owned(),
            fasta_path: input.database.fasta,
            peptide_database_path: input.database.peptides,
            spectra_paths: input.mzml_paths.unwrap_or_default(),
            note: "Local paths are normalized and constrained to the server root; remote URLs are disabled.",
        })
    }

    fn start_search(&self, args: StartSearchArgs) -> anyhow::Result<JobRecord> {
        ensure!(
            args.approved,
            "search execution requires explicit approval (`approved: true`)"
        );
        ensure!(
            args.batch_size.unwrap_or(1) > 0,
            "batch_size must be greater than zero"
        );
        let (config_path, mut input) = self.load_config(&args.config_path)?;
        let job_id = Uuid::new_v4().to_string();
        let job_dir = self.jobs_dir.join(&job_id);
        let output_dir = job_dir.join("output");
        fs::create_dir_all(&output_dir)?;
        input.output_directory = Some(output_dir.to_string_lossy().into_owned());
        let events_path = job_dir.join("events.jsonl");
        let events = EventEmitter::from_writer(BufWriter::new(File::create(&events_path)?));
        let cancellation = CancellationToken::default();
        let timestamp = now();
        let record = JobRecord {
            job_id: job_id.clone(),
            status: JobStatus::Queued,
            config_path: config_path.to_string_lossy().into_owned(),
            job_directory: job_dir.to_string_lossy().into_owned(),
            events_path: events_path.to_string_lossy().into_owned(),
            output_directory: output_dir.to_string_lossy().into_owned(),
            created_at_unix: timestamp,
            updated_at_unix: timestamp,
            summary: None,
            error: None,
        };
        write_record(&record)?;
        self.jobs.write().expect("job registry poisoned").insert(
            job_id.clone(),
            JobEntry {
                record: record.clone(),
                cancellation: Some(cancellation.clone()),
            },
        );

        let jobs = self.jobs.clone();
        let parallel = args.batch_size.unwrap_or_else(|| (num_cpus() / 2).max(1));
        let spawn = std::thread::Builder::new()
            .name(format!("sage-job-{job_id}"))
            .spawn(move || {
                update_job(&jobs, &job_id, |entry| {
                    entry.record.status = JobStatus::Running;
                });
                let result = SageRunner::new(
                    input,
                    JobOptions {
                        parallel,
                        parquet: args.parquet.unwrap_or(false),
                        events,
                        cancellation: cancellation.clone(),
                    },
                )
                .run();
                update_job(&jobs, &job_id, |entry| {
                    entry.cancellation = None;
                    match result {
                        Ok(result) => {
                            entry.record.status = JobStatus::Completed;
                            entry.record.summary = Some(result.summary);
                        }
                        Err(error) if cancellation.is_cancelled() => {
                            entry.record.status = JobStatus::Cancelled;
                            entry.record.error = Some(error.to_string());
                        }
                        Err(error) => {
                            entry.record.status = JobStatus::Failed;
                            entry.record.error = Some(error.to_string());
                        }
                    }
                });
            });
        if let Err(error) = spawn {
            update_job(&self.jobs, &record.job_id, |entry| {
                entry.cancellation = None;
                entry.record.status = JobStatus::Failed;
                entry.record.error = Some(format!("failed to start search thread: {error}"));
            });
            return Err(error.into());
        }
        Ok(record)
    }

    fn job(&self, job_id: &str) -> anyhow::Result<JobRecord> {
        self.jobs
            .read()
            .expect("job registry poisoned")
            .get(job_id)
            .map(|entry| entry.record.clone())
            .with_context(|| format!("unknown Sage job `{job_id}`"))
    }

    fn list_jobs(&self) -> Vec<JobRecord> {
        let mut jobs = self
            .jobs
            .read()
            .expect("job registry poisoned")
            .values()
            .map(|entry| entry.record.clone())
            .collect::<Vec<_>>();
        jobs.sort_by_key(|job| std::cmp::Reverse(job.created_at_unix));
        jobs
    }

    fn cancel(&self, job_id: &str) -> anyhow::Result<JobRecord> {
        let mut jobs = self.jobs.write().expect("job registry poisoned");
        let entry = jobs
            .get_mut(job_id)
            .with_context(|| format!("unknown Sage job `{job_id}`"))?;
        ensure!(
            matches!(entry.record.status, JobStatus::Queued | JobStatus::Running),
            "job is not running"
        );
        entry
            .cancellation
            .as_ref()
            .context("job cannot be cancelled after server restart")?
            .cancel();
        entry.record.status = JobStatus::Cancelling;
        entry.record.updated_at_unix = now();
        write_record(&entry.record)?;
        Ok(entry.record.clone())
    }

    fn events(&self, args: JobEventsArgs) -> anyhow::Result<Vec<serde_json::Value>> {
        let record = self.job(&args.job_id)?;
        let reader = BufReader::new(File::open(&record.events_path)?);
        let after = args.after_sequence;
        let limit = args.limit.unwrap_or(200).min(1_000);
        reader
            .lines()
            .map(|line| Ok(serde_json::from_str::<serde_json::Value>(&line?)?))
            .filter(|value| {
                value
                    .as_ref()
                    .ok()
                    .and_then(|value| value["sequence"].as_u64())
                    .map(|sequence| after.map(|after| sequence > after).unwrap_or(true))
                    .unwrap_or(true)
            })
            .take(limit)
            .collect()
    }

    fn resource(&self, uri: &str) -> anyhow::Result<String> {
        let suffix = uri
            .strip_prefix("sage://jobs/")
            .context("unsupported Sage resource URI")?;
        let (job_id, kind) = suffix
            .split_once('/')
            .context("resource URI must include a job id and resource name")?;
        let record = self.job(job_id)?;
        match kind {
            "manifest" | "summary" => Ok(serde_json::to_string_pretty(&record)?),
            "events" => Ok(fs::read_to_string(record.events_path)?),
            _ => anyhow::bail!("unknown job resource `{kind}`"),
        }
    }
}

fn now() -> u64 {
    SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .unwrap_or_default()
        .as_secs()
}

fn write_record(record: &JobRecord) -> anyhow::Result<()> {
    let path = Path::new(&record.job_directory).join("job.json");
    fs::write(path, serde_json::to_vec_pretty(record)?)?;
    Ok(())
}

fn update_job(
    jobs: &RwLock<HashMap<String, JobEntry>>,
    job_id: &str,
    update: impl FnOnce(&mut JobEntry),
) {
    let mut jobs = jobs.write().expect("job registry poisoned");
    if let Some(entry) = jobs.get_mut(job_id) {
        update(entry);
        entry.record.updated_at_unix = now();
        if let Err(error) = write_record(&entry.record) {
            entry.record.error = Some(format!("failed to persist job state: {error}"));
        }
    }
}

fn num_cpus() -> usize {
    std::thread::available_parallelism()
        .map(usize::from)
        .unwrap_or(1)
}

#[derive(Debug, Deserialize, JsonSchema)]
pub struct ConfigPathArgs {
    /// JSON configuration path, relative to the server root or absolute within it.
    pub config_path: String,
}

#[derive(Debug, Deserialize, JsonSchema)]
pub struct StartSearchArgs {
    /// JSON configuration path, relative to the server root or absolute within it.
    pub config_path: String,
    /// Must be true only after the user approves the potentially expensive search.
    pub approved: bool,
    /// Write Parquet output instead of TSV.
    pub parquet: Option<bool>,
    /// Number of spectra files to process in each batch.
    pub batch_size: Option<usize>,
}

#[derive(Debug, Deserialize, JsonSchema)]
pub struct JobIdArgs {
    pub job_id: String,
}

#[derive(Debug, Deserialize, JsonSchema)]
pub struct JobEventsArgs {
    pub job_id: String,
    /// Return only events whose sequence is greater than this value.
    pub after_sequence: Option<u64>,
    /// Maximum number of events to return (default 200, maximum 1000).
    pub limit: Option<usize>,
}

fn tool_result<T: Serialize>(result: anyhow::Result<T>) -> Result<CallToolResult, McpError> {
    Ok(match result {
        Ok(value) => match serde_json::to_value(value) {
            Ok(value) => CallToolResult::structured(value),
            Err(error) => CallToolResult::structured_error(serde_json::json!({
                "error": "serialization_failed",
                "message": error.to_string(),
            })),
        },
        Err(error) => CallToolResult::structured_error(serde_json::json!({
            "error": "sage_operation_failed",
            "message": error.to_string(),
        })),
    })
}

async fn blocking<T, F>(operation: F) -> anyhow::Result<T>
where
    T: Send + 'static,
    F: FnOnce() -> anyhow::Result<T> + Send + 'static,
{
    tokio::task::spawn_blocking(operation).await?
}

#[derive(Clone)]
pub struct SageMcp {
    state: State,
}

impl SageMcp {
    pub fn new(root: PathBuf, jobs_dir: Option<PathBuf>) -> anyhow::Result<Self> {
        Ok(Self {
            state: State::new(root, jobs_dir)?,
        })
    }
}

#[tool_router]
impl SageMcp {
    #[tool(
        description = "Validate a Sage JSON configuration and all local input paths without running a search"
    )]
    async fn validate_config(
        &self,
        Parameters(args): Parameters<ConfigPathArgs>,
    ) -> Result<CallToolResult, McpError> {
        let state = self.state.clone();
        tool_result(blocking(move || state.validate_config(&args.config_path)).await)
    }

    #[tool(
        description = "Start an approved Sage search in the background and return a persistent job id immediately"
    )]
    async fn start_search(
        &self,
        Parameters(args): Parameters<StartSearchArgs>,
    ) -> Result<CallToolResult, McpError> {
        let state = self.state.clone();
        tool_result(blocking(move || state.start_search(args)).await)
    }

    #[tool(
        description = "Get structured status, output paths, errors, and summary statistics for a Sage job"
    )]
    fn get_job_status(
        &self,
        Parameters(args): Parameters<JobIdArgs>,
    ) -> Result<CallToolResult, McpError> {
        tool_result(self.state.job(&args.job_id))
    }

    #[tool(
        description = "List known Sage jobs, including jobs restored after an MCP server restart"
    )]
    fn list_jobs(&self) -> Result<CallToolResult, McpError> {
        tool_result(Ok(self.state.list_jobs()))
    }

    #[tool(description = "Cooperatively cancel a queued or running Sage search")]
    fn cancel_search(
        &self,
        Parameters(args): Parameters<JobIdArgs>,
    ) -> Result<CallToolResult, McpError> {
        tool_result(self.state.cancel(&args.job_id))
    }

    #[tool(
        description = "Read bounded structured progress events for a Sage job; use after_sequence for incremental polling"
    )]
    async fn get_job_events(
        &self,
        Parameters(args): Parameters<JobEventsArgs>,
    ) -> Result<CallToolResult, McpError> {
        let state = self.state.clone();
        tool_result(blocking(move || state.events(args)).await)
    }
}

#[tool_handler(
    name = "sage-mcp",
    version = "0.1.0",
    instructions = "Operate Sage through validated, root-bounded local paths. Validate configurations before searches. Never call start_search with approved=true until the user has approved the compute and output location. Poll status/events rather than starting duplicate jobs."
)]
impl ServerHandler for SageMcp {
    fn get_info(&self) -> ServerInfo {
        ServerInfo::new(
            ServerCapabilities::builder()
                .enable_tools()
                .enable_resources()
                .build(),
        )
    }

    async fn list_resources(
        &self,
        _request: Option<PaginatedRequestParams>,
        _context: RequestContext<RoleServer>,
    ) -> Result<ListResourcesResult, McpError> {
        let resources = self
            .state
            .list_jobs()
            .into_iter()
            .flat_map(|job| {
                let base = format!("sage://jobs/{}", job.job_id);
                [
                    Resource::new(
                        format!("{base}/manifest"),
                        format!("{} manifest", job.job_id),
                    ),
                    Resource::new(format!("{base}/events"), format!("{} events", job.job_id)),
                ]
            })
            .collect();
        Ok(ListResourcesResult::with_all_items(resources))
    }

    async fn read_resource(
        &self,
        request: ReadResourceRequestParams,
        _context: RequestContext<RoleServer>,
    ) -> Result<ReadResourceResponse, McpError> {
        let state = self.state.clone();
        let uri = request.uri.clone();
        blocking(move || state.resource(&uri))
            .await
            .map(|text| {
                ReadResourceResult::new(vec![ResourceContents::text(text, request.uri.as_str())])
                    .into()
            })
            .map_err(|error| {
                McpError::resource_not_found(
                    error.to_string(),
                    Some(serde_json::json!({ "uri": request.uri })),
                )
            })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn rejects_paths_outside_root() {
        let temp = std::env::temp_dir().join(format!("sage-mcp-test-{}", Uuid::new_v4()));
        fs::create_dir_all(&temp).unwrap();
        let state = State::new(temp.clone(), None).unwrap();
        let error = state
            .resolve_existing("../outside")
            .unwrap_err()
            .to_string();
        assert!(error.contains("cannot access") || error.contains("outside"));
        fs::remove_dir_all(temp).unwrap();
    }

    #[test]
    fn explicit_approval_is_required_before_input_access() {
        let temp = std::env::temp_dir().join(format!("sage-mcp-test-{}", Uuid::new_v4()));
        fs::create_dir_all(&temp).unwrap();
        let state = State::new(temp.clone(), None).unwrap();
        let error = state
            .start_search(StartSearchArgs {
                config_path: "missing.json".into(),
                approved: false,
                parquet: None,
                batch_size: None,
            })
            .unwrap_err()
            .to_string();
        assert!(error.contains("explicit approval"));
        fs::remove_dir_all(temp).unwrap();
    }
}
