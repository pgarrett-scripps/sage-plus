use anyhow::{bail, Context};
use rmcp::{transport::stdio, ServiceExt};
use sage_mcp::{run_worker, SageMcp};
use std::path::PathBuf;

enum Mode {
    Server {
        root: PathBuf,
        jobs_dir: Option<PathBuf>,
    },
    Worker(PathBuf),
}

fn usage() -> &'static str {
    "Usage: sage-mcp [--root PATH] [--jobs-dir PATH]\n\
     \n\
     --root PATH      Restrict all configuration and input files to PATH (default: cwd)\n\
     --jobs-dir PATH  Store job manifests, events, and outputs here (default: ROOT/.sage/jobs)\n\
     -V, --version    Print version information"
}

fn arguments() -> anyhow::Result<Mode> {
    let mut root = std::env::current_dir()?;
    let mut jobs_dir = None;
    let mut args = std::env::args().skip(1).peekable();
    if args.peek().map(String::as_str) == Some("--worker") {
        args.next();
        let request = PathBuf::from(args.next().context("--worker requires a request path")?);
        if let Some(unknown) = args.next() {
            bail!("unexpected worker argument `{unknown}`");
        }
        return Ok(Mode::Worker(request));
    }
    while let Some(arg) = args.next() {
        match arg.as_str() {
            "--root" => {
                root = PathBuf::from(args.next().context("--root requires a path")?);
            }
            "--jobs-dir" => {
                jobs_dir = Some(PathBuf::from(
                    args.next().context("--jobs-dir requires a path")?,
                ));
            }
            "-h" | "--help" => {
                eprintln!("{}", usage());
                std::process::exit(0);
            }
            "-V" | "--version" => {
                println!("sage-mcp {}", env!("CARGO_PKG_VERSION"));
                std::process::exit(0);
            }
            unknown => bail!("unknown argument `{unknown}`\n{}", usage()),
        }
    }
    Ok(Mode::Server { root, jobs_dir })
}

fn main() -> anyhow::Result<()> {
    match arguments()? {
        Mode::Worker(request) => run_worker(&request),
        Mode::Server { root, jobs_dir } => tokio::runtime::Builder::new_multi_thread()
            .enable_all()
            .build()?
            .block_on(async move {
                let server = SageMcp::new(root, jobs_dir)?;
                let shutdown = server.clone();
                let service = server.serve(stdio()).await?;
                service.waiting().await?;
                shutdown.shutdown_workers();
                Ok(())
            }),
    }
}
