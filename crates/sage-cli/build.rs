//! Record the git commit Sage is built from as `SAGE_GIT_COMMIT`, for the
//! provenance written to outputs. A `SAGE_GIT_COMMIT` environment variable
//! wins; builds outside a git checkout (or without git) record none.

use std::process::Command;

fn git(args: &[&str]) -> Option<String> {
    let output = Command::new("git").args(args).output().ok()?;
    output
        .status
        .success()
        .then(|| String::from_utf8(output.stdout).ok())
        .flatten()
        .map(|text| text.trim().to_string())
        .filter(|text| !text.is_empty())
}

fn main() {
    println!("cargo:rerun-if-env-changed=SAGE_GIT_COMMIT");
    let commit = match std::env::var("SAGE_GIT_COMMIT") {
        Ok(commit) if !commit.trim().is_empty() => Some(commit.trim().to_string()),
        _ => {
            // Rebuild when HEAD moves: HEAD itself, the branch it names, and
            // packed refs. Paths come from git so worktrees resolve correctly.
            for path in [
                git(&["rev-parse", "--path-format=absolute", "--git-path", "HEAD"]),
                git(&["symbolic-ref", "-q", "HEAD"]).and_then(|reference| {
                    git(&[
                        "rev-parse",
                        "--path-format=absolute",
                        "--git-path",
                        &reference,
                    ])
                }),
                git(&[
                    "rev-parse",
                    "--path-format=absolute",
                    "--git-path",
                    "packed-refs",
                ]),
            ]
            .into_iter()
            .flatten()
            // A missing path would rerun this script on every build.
            .filter(|path| std::path::Path::new(path).exists())
            {
                println!("cargo:rerun-if-changed={path}");
            }
            git(&["rev-parse", "HEAD"])
        }
    };
    println!(
        "cargo:rustc-env=SAGE_GIT_COMMIT={}",
        commit.unwrap_or_default()
    );
}
