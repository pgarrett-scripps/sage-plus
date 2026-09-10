# Beta.3 release checklist

Prepared September 10, 2026 for `v0.1.0-beta.3` on `codex/release-beta3`.
This is a hardening release. It does not change scoring defaults or establish
broader scientific calibration.

## Dependency consolidation

The release branch now incorporates the updates proposed in PRs #7, #10, #11,
#12, and #19: anyhow 1.0.104, clap 4.6.1 with clap_builder 4.6.0, itoa 1.0.18,
schemars 1.2.2, and pinned upload-artifact v7.0.1. Schemars also updates
schemars_derive to 1.2.2 and serde_derive_internals to 0.30.0, using the existing
syn 3.0.3 lockfile entry. Unrelated Windows dependency resolution changes were
excluded. The existing tempfile/getrandom fix is retained.

PR #10's branch name mentions 4.6.6, but its actual diff specifies clap 4.6.1
and clap_builder 4.6.0. This consolidation follows the reviewed diff.
DashMap remains at 5.5.3. Its major upgrade in PR #8 is deferred for focused
concurrency, performance, and scientific-output validation.

The consolidated dependency source is commit
`e405067dba6ad04b7422ef9eb3c98c3651dde938`. Fresh local evidence is retained in
`benchmarks/results/beta3-deps-20260910/`, including the committed source archive,
lockfile and binary hashes, commands, logs, and `validation-summary.json`.
The search binary SHA-256 is
`12f7574e573dc09befe3130b7a690b8044f7319a00664199f166cf987fc4a37e`.
The MCP binary SHA-256 is
`a572384757c967c096a739cc7bf19ceeaa96c86ada2e70564662214b9d90c845`.
Subsequent changes to this evidence document do not change compiled sources,
dependencies, workflows, or files shipped in release archives.

Fresh validation passed:

- 399 default and 399 minimal-feature workspace tests, 49 storage tests, and
  16 Python tests. The first storage test attempt could not bind its loopback
  server inside the sandbox. An unchanged retry with networking enabled passed.
- Strict Clippy, Rust 1.88 compatibility, formatting, version inheritance, and
  the optimized workspace build.
- The online audit found zero vulnerabilities, unsoundness findings, or yanked
  dependencies. The existing instant and paste maintenance warnings remain.
- CLI help, version, missing arguments, unknown arguments, and invalid batch-size
  responses match the preceding beta.3 binary. Generated configuration schema
  bytes and complete MCP tool definitions and schemas also match.
- Four paired workloads, each with one warmup and three measured trials per
  binary, passed without timing, memory, or identification review flags.
  Accepted PSM counts remain 2203, 2256, 3266, and 2203 for standard, modified,
  feature, and prefilter searches. Compared scores agree exactly. Every compared
  PSM, LFQ, matched-fragment, and spectral-library Parquet file is byte-identical.
- A fresh paired 20-seed entrapment experiment completed without evaluator
  failures. All stored and accepted identities, threshold counts, and FDP
  estimates agree. Independent reconstruction passed all 200 count checks.
  Spectrum q-values and hyperscores agree exactly. Maximum differences are
  8e-8 for peptide q-values and 7e-7 for discriminant scores, without acceptance
  changes. This remains a bounded regression check, not broad calibration.

The comparator for these fresh experiments is the preceding validated beta.3
binary identified below. The earlier evidence remains tied to its original
binaries and is retained separately.

[Hosted Rust checks](https://github.com/pgarrett-scripps/sage-plus/actions/runs/34536457819)
and [dependency security](https://github.com/pgarrett-scripps/sage-plus/actions/runs/34536457813)
passed on the consolidated source, including Windows and macOS smoke tests and
the coverage gate. The
[release packaging run](https://github.com/pgarrett-scripps/sage-plus/actions/runs/34536503116)
validates that same source. Successful platform archives, assembly, and container
validation are required before merge. Final statuses are available through these
run links and PR #21. The five bot PRs should be closed as superseded only after
their equivalent updates are merged through PR #21.

## Earlier local gates

| Gate | Status |
|---|---|
| Workspace version and changelog | beta.3, matching tag check passed |
| Online dependency audit | Passed, zero vulnerabilities, unsoundness findings, or yanked dependencies |
| Maintenance warnings | `instant` and `paste` remain tracked |
| Default and minimal workspace tests | 399 passed in each configuration |
| Patched filemanager compatibility | 2 new tests and 47 upstream documentation tests passed |
| Strict Clippy and Rust 1.88 compatibility check | Passed |
| Python benchmark tooling | 16 tests passed |
| Formatting and schema synchronization | Passed |
| Workspace line coverage | 81.70%, above the 80% gate |
| Optimized native archive and smoke search | Passed, both binaries report beta.3, checksums verified |
| Four paired representative workloads | Passed, no timing, memory, or identification review flags |
| Bounded entrapment comparison | 20 seeds completed, identities and FDP estimates agree, 200 independent count checks passed |

Workspace tests include the committed real Bruker TDF and Thermo RAW fixtures,
mzMLb, MCP worker isolation, and the new hardening regressions. The storage patch
also exercises S3 XML and HTTP ranges through a loopback fixture. These checks do
not certify every cloud provider's authentication or acquisition workflow.

## Earlier binary and workload evidence

The tested code and packaging source is commit `5bc965f62fa68250c1b1a8b6fa74da2fd6cc40f1`.
The Windows schema checkout fix in `379b73a` and subsequent release-evidence
updates do not change the compiled sources.
The frozen search binary SHA-256 is
`c9fa8e086d68c13eae2c68d0008afb58b4b6ae722bbefae572ac1fcfc38f3972`.
The MCP binary SHA-256 is
`a03999c21226e152ba63ef340af39921f5712cc2c6fae548649df64119688158`.

Four paired conditions used the frozen pre-hardening baseline, one warmup and
three measured trials per build, eight workers, batch size one, and alternating
engine order. All stored PSM identities and compared q-values, hyperscores, and
discriminant scores agree exactly in the first-trial comparisons.

| Condition | PSMs at 1%, both builds | Baseline median wall | Beta.3 median wall | RSS change |
|---|---:|---:|---:|---:|
| standard | 2203 | 6.77 s | 6.79 s | +0.28% |
| modified | 2256 | 14.23 s | 14.23 s | -3.87% |
| feature | 3266 | 15.56 s | 15.50 s | +0.48% |
| prefilter | 2203 | 11.26 s | 11.43 s | +0.40% |

Main results, LFQ, matched fragments, and spectral-library Parquet outputs are
byte-identical. Text mzSpecLib output differs only in three software-version
metadata lines. Exact prefiltering produces byte-identical candidate PSM output.
These timings are a bounded regression check, not a precise performance claim.

The local evidence directory `benchmarks/results/beta3-20260910/` contains the
committed source archive, build and input manifests, engineering logs, audit JSON,
paired workloads, entrapment outputs, and native archive with checksum and smoke
results. The directory is intentionally excluded from Git.

## Earlier entrapment check

The final beta.3 binary and frozen pre-hardening baseline completed the same
20-seed experiment, seeds 20260902 through 20260921. All stored identities,
accepted identities, target and entrapment counts, and FDP estimates agree at the
five reported cutoffs. Independent reconstruction verified all 200 engine, seed,
and cutoff count pairs. Spectrum q-values and hyperscores agree exactly.
Maximum paired differences are 8e-8 for peptide q-values and 5e-7 for discriminant
scores, with no identity or cutoff changes. The earlier unchanged-baseline repeat
also exhibited small floating-point differences. This is not a claim of bitwise
reproducibility.

FDRBench's JVM crashed once in its C2 compiler during the baseline FDP evaluation
for seed 20260911. The harness stopped and left no complete aggregate manifest.
The failure log and JVM crash report were retained. Retrying the incomplete stage
with unchanged executable, inputs, and arguments succeeded, and the complete
20-seed matrix then passed analysis and independent count verification. No failed
observation was converted into a zero-discovery result.

The experiment retains the limitations in the September 4 report. It uses one
HEK file and a reduced FASTA, with low power near 1% and above-nominal estimates
at higher cutoffs in both builds.

## Hosted gates

This checklist records preparation status before merge. The linked PR and
workflow runs show subsequent completion.

- [x] Release preparation is in `pgarrett-scripps/sage-plus`, branch `codex/release-beta3`, with [PR #21](https://github.com/pgarrett-scripps/sage-plus/pull/21).
- [x] Prior hosted Rust CI passed on `379b73a`, including both toolchains, Clippy, coverage, and Windows/macOS smoke tests. Logs are retained with the local evidence.
- [x] Prior packaging validation passed on `5bc965f` for all seven archives and the AMD64 container.
- [x] All seven downloaded archive checksums, required files, executable modes, and schema contents were verified.
- [x] Hosted GNU and musl x86_64 Linux archives passed smoke searches.
- [x] Current repository PR checks pass on consolidated dependency source `e405067`. Any subsequent evidence-only commit must also pass required PR checks before merge.
- [ ] Current repository release dry run passes on the final preparation.
- [ ] Reviewed preparation is merged into `main`.
- [ ] Annotated `v0.1.0-beta.3` tag points to the verified release source.
- [ ] Tag workflow publishes the prerelease, checksums, and versioned container.

The first Windows run found CRLF conversion of the generated schema. Commit
`379b73a` pins `schemas/*.json` to LF without weakening schema comparison or
changing search code. Hosted Windows tests passed after that fix. The prior
packaging validation used `5bc965f`, before the checkout fix. Compiled sources
and packaging commands are unchanged by the fix. Current repository CI and a
fresh release dry run validate the final preparation before publication.

Both downloaded Linux archive variants report beta.3, enable mzMLb, and find the
expected peptide in the synthetic smoke search. The GNU archive SHA-256 is
`b75d7ca4b87075590e71f970421d988e53df813e5b91f33bbca97b741b63de2d`.
The musl archive SHA-256 is
`fc55a359306eb56fadfe3a1d4c213ab791fb809614bea70bbfe4fb48d3f6f37e`.
Earlier CI logs, build identities, archives, and checksum verification records are
retained in `benchmarks/results/beta3-20260910/`.

Release operations use the configured `pgarrett-scripps` owner account and the
existing `pgarrett-scripps/sage-plus` repository. Check the active account and
remote before making GitHub changes. An inactive authorized account should be
selected before treating a permission error as a release blocker.

## Release scope and evidence

The release fixes configured batching, gzip completion, event ordering, input
validation, cancellation boundaries, and MCP persistence. It introduces explicit
local overwrite, run-summary schema 9, verified benchmark manifests, and CI gates.
The pinned storage patch removes the remaining vulnerable transitive XML parser.

The [September 4 results](HARDENING_RESULTS.md) retain the original hashes,
timing investigations, and 20-seed entrapment evidence. The experiment uses one
HEK file and an incomplete 2,000-protein target FASTA. It has low power near 1%,
and estimates at higher cutoffs remain above nominal in both builds. Independent
study, protein-group, localization, transfer, and library calibration remain in
the [scientific protocol](SCIENTIFIC_PROTOCOL.md).

Reproduce the local experiments with [the hardening commands](HARDENING_USAGE.md).

See the [security review](SECURITY_REVIEW.md), [changelog](../CHANGELOG.md), and
[maintainer release procedure](../RELEASING.md) for the dependency disposition,
compatibility changes, and publication commands.
