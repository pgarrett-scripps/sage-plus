# Beta.3 release checklist

Prepared September 10, 2026 for `v0.1.0-beta.3` on `codex/release-beta3`.
This is a hardening release. It does not change scoring defaults or establish
broader scientific calibration.

## Local gates

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

## Final binary and workload evidence

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

## Final entrapment check

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

- [x] Release preparation committed and available in [PR #20](https://github.com/pgarrett-scripps/sage-plus/pull/20).
- [x] Fork Rust CI passes on `379b73a`, including both toolchains, Clippy, coverage, and platform smoke tests.
- [ ] Upstream PR checks pass after maintainer approval of fork workflows.
- [x] Windows and macOS smoke tests pass on `379b73a` in the fork.
- [x] Fork manual release workflow builds all seven archives and the AMD64 container on `5bc965f`.
- [ ] Final upstream release dry run passes on the reviewed, merged preparation.
- [x] All seven release-dist checksums, required files, executable modes, and schema contents verified.
- [x] Hosted GNU and musl x86_64 Linux archives downloaded and smoke-tested.
- [ ] Maintainer merges the reviewed preparation into `main`.
- [ ] Annotated `v0.1.0-beta.3` tag points to the verified release source.
- [ ] Tag workflow publishes the prerelease, checksums, and versioned container.

The first fork Windows run found CRLF conversion of the generated schema.
Commit `379b73a` pins `schemas/*.json` to LF without weakening schema comparison
or changing search code. The corrected commit has a fresh [Rust CI run](https://github.com/afk-sapien/sage-plus/actions/runs/34532870240).
All jobs passed, including Rust 1.88 and 1.97.1 tests, strict Clippy, coverage,
optimized builds, and Windows and macOS smoke tests. The fork push run correctly
skips PR-only dependency review. The upstream PR review jobs require maintainer
approval and have not run.

The [fork release dry run](https://github.com/afk-sapien/sage-plus/actions/runs/34532258817)
uses `5bc965f`, before the schema checkout fix. The compiled sources and packaging
commands are unchanged by that fix. All seven builds, checksum assembly, and the
AMD64 container build passed. Publication was correctly skipped for manual dispatch.
All seven downloaded archive hashes match `SHA256SUMS`. Required documentation,
licenses, executable modes, and parsed schema contents were verified.

Both the GNU and static musl x86_64 archives were unpacked and smoke-tested. Both
binaries report beta.3, mzMLb is enabled, and each synthetic search finds the
expected peptide. The GNU archive SHA-256 is
`b75d7ca4b87075590e71f970421d988e53df813e5b91f33bbca97b741b63de2d`.
The musl archive SHA-256 is
`fc55a359306eb56fadfe3a1d4c213ab791fb809614bea70bbfe4fb48d3f6f37e`.
The final upstream dry run must still validate the merged release preparation.

GitHub currently authenticates as `afk-sapien` with read-only access to
`pgarrett-scripps/sage-plus`. Publishing requires a maintainer session or a
maintainer to complete merge, workflow dispatch, and tagging. A local archive
does not substitute for the hosted platform and publication gates.

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
