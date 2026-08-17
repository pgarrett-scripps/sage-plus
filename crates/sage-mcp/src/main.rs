use anyhow::{bail, Context};
use rmcp::{transport::stdio, ServiceExt};
use sage_mcp::SageMcp;
use std::path::PathBuf;

fn usage() -> &'static str {
    "Usage: sage-mcp [--root PATH] [--jobs-dir PATH]\n\
     \n\
     --root PATH      Restrict all configuration and input files to PATH (default: cwd)\n\
     --jobs-dir PATH  Store job manifests, events, and outputs here (default: ROOT/.sage/jobs)\n\
     -V, --version    Print version information"
}

fn arguments() -> anyhow::Result<(PathBuf, Option<PathBuf>)> {
    let mut root = std::env::current_dir()?;
    let mut jobs_dir = None;
    let mut args = std::env::args().skip(1);
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
    Ok((root, jobs_dir))
}

#[tokio::main]
async fn main() -> anyhow::Result<()> {
    let (root, jobs_dir) = arguments()?;
    let server = SageMcp::new(root, jobs_dir)?;
    let service = server.serve(stdio()).await?;
    service.waiting().await?;
    Ok(())
}
