# Working in this repository

Sage Plus is an independently maintained downstream distribution of
[Sage](https://github.com/lazear/sage). Contributor policy is in `CONTRIBUTING.md`, the feature
list in `FEATURES.md`, and the upstream relationship in `UPSTREAM.md`. This file covers how to
build, test and change the code.

## Building and testing

`.cargo/config.toml` sets `CMAKE_POLICY_VERSION_MINIMUM=3.5` for the whole workspace. It is
required: the default `mzmlb` feature pulls in `blosc-src`, whose vendored c-blosc declares a
CMake policy floor that CMake 4 removed. Do not drop it, and do not export it by hand — the
config file already applies to `build`, `test` and `clippy`.

Run before every commit:

```
cargo fmt
cargo clippy --workspace --all-targets --all-features --locked -- -D warnings
```

That clippy line is exactly what CI runs (`.github/workflows/rust.yml`); a shorter invocation
can pass locally and still fail CI. CI also requires at least 80% workspace line coverage
(`cargo llvm-cov --workspace --summary-only` locally).

`cargo test --workspace` is a heavy job, not a quick check: a cold build takes roughly 45
minutes and the test run peaks around 3.7 GB of resident memory. On a shared machine submit it
through the scheduler (`ajs`) with 4 cores and a 5 GB reservation rather than running it in the
foreground. A single crate (`cargo test -p sage-core`) is fine to run directly.

Use `--profile fast-release` for local experiments and work-in-progress binaries: thin LTO and
incremental compilation make rebuilds seconds instead of minutes. Use `--release` for shipped
binaries and for any timing number you report.

## Where tests go

No test implementations in production source files. A unit test file lives in
`crates/<crate>/tests/unit/<module>.rs` and is attached to its module with three lines at the
bottom of the source file:

```rust
#[cfg(test)]
#[path = "../tests/unit/ambiguity.rs"]
mod test;
```

Black-box integration tests go directly in `crates/<crate>/tests/`. This is the convention most
often broken by accident; check it before adding a test.

## Changing behaviour

Only the settings that define the search space are required in a configuration: enzyme, static
modifications, ion kinds. Every other option, existing or new, must have a default that
reproduces the previous behaviour, so a configuration written against an older version keeps
producing the same search. Adding a parameter whose default changes results is a breaking
change, not a feature.

Public API, CLI flags, configuration keys and output column names are a compatibility surface.
Renaming or removing one is a deliberate decision, not a cleanup.

Use the field's standard term in code, documentation and output columns, one term per concept:
intensity (not height), modifications or PTMs (not marks). Do not invent terminology.

## This is a fork

Keep diffs against `upstream/master` small and explainable. `git log upstream/master..main` is
the real feature list, so a change that could live upstream should go upstream instead. When a
change cannot, say why in the commit message, and record the divergence in `UPSTREAM.md` if it
is a place where Sage Plus and upstream solve the same problem differently.

## The reader crates

`sage-plus-tdf` (Bruker timsTOF `.d`/TDF) and `sage-plus-raw` (Thermo `.raw`) are separate
repositories, consumed here as published crates.io dependencies. A reader change reaches Sage
Plus only through a crates.io release.

To try an unreleased reader locally, add a `[patch.crates-io]` section or a path override —
and never commit it. A committed override breaks the build for everyone who does not have that
checkout and silently hides which reader version the results came from.
