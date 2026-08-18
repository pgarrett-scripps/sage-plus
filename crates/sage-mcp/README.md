# sage-mcp

`sage-mcp` is a local stdio Model Context Protocol server that lets AI clients validate,
launch, monitor, and cancel Sage searches through typed tools.

The server deliberately separates inspection from execution. `start_search` requires
`approved: true`, local files must remain under `--root`, remote URLs are disabled, and all
outputs are forced into a server-managed job directory. Job manifests and JSONL events are
persistent and are restored when the server restarts.

## Build

Sage and the official MCP Rust SDK require Rust 1.88 or newer. The repository pins Rust
1.97.1 for reproducible development and CI also verifies the 1.88 minimum.

```shell
cd crates/sage-mcp
cargo build --release
```

## Run

```shell
sage-mcp --root /path/to/allowed/data
```

Example MCP client configuration:

```json
{
  "mcpServers": {
    "sage": {
      "command": "/absolute/path/to/sage-mcp",
      "args": ["--root", "/absolute/path/to/allowed/data"]
    }
  }
}
```

Tools:

- `get_capabilities`
- `validate_config`
- `inspect_config`
- `estimate_search`
- `start_search`
- `get_job_status`
- `list_jobs`
- `cancel_search`
- `get_job_events`
- `summarize_run`
- `analyze_run`
- `query_results`

Every successful Sage run writes `run-summary.json` into its output directory, whether or not
MCP is used. In addition to identification and throughput statistics, it records localized PTMs,
model/alignment outcomes, quantification, memory controls, input formats, and modification limits.
`analyze_run` reads that portable artifact and highlights basic outcomes that need attention.

`query_results` reads at most 200 matching rows and scans at most one million rows from Parquet
outputs. It supports PSMs, localized PTM sites, collapsed protein sites, and empirical spectral
library transitions, with optional q-value, protein, peptide, and modification filters. Use the
`spectral_library` dataset to query `spectral_library.sage.parquet`.

Each job also exposes `sage://jobs/{job_id}/manifest`, `sage://jobs/{job_id}/summary`, and
`sage://jobs/{job_id}/events` resources.
