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

- `validate_config`
- `start_search`
- `get_job_status`
- `list_jobs`
- `cancel_search`
- `get_job_events`

Each job also exposes `sage://jobs/{job_id}/manifest` and
`sage://jobs/{job_id}/events` resources.
