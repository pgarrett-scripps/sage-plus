# Sage Plus benchmark pipeline

This directory contains a small repeatable benchmark harness for Sage Plus. It is intended for
development checks and README-sized performance summaries. It is not a scientific validation
suite.

The [September 2026 hardening results](HARDENING_RESULTS.md) record fresh paired workloads,
runtime follow-up, entrapment checks, and the remaining release gates. The
[scientific protocol](SCIENTIFIC_PROTOCOL.md) defines the next validation milestone.
The [public scientific pilot](SCIENTIFIC_PILOT.md) documents the separate beta.3
validation runs, data provenance, analysis tools and interpretation limits.
The [pilot report and figures](scientific-results/20260914/SCIENTIFIC_REPORT.md)
retain the measured results and unresolved PTM and MBR confidence questions.

The harness compares the current working tree with a pinned baseline. By default the baseline is
the `v0.1.0-beta.1` release. The candidate build includes uncommitted working-tree changes. Set
`BASELINE_REF` when evaluating a different release or development boundary.

## Mass offset evaluation

`run_mass_offset.py` runs the mass offset matrix with recorded executable,
configuration, input, and output identities, and `summarize_mass_offset.py`
reduces it to the retained summary the chapter reads:

```shell
python3 benchmarks/run_mass_offset.py --root /data/sage-plus-scientific/mass-offset-20260919 \
    --sage target/release/sage
python3 benchmarks/summarize_mass_offset.py --root /data/sage-plus-scientific/mass-offset-20260919
```

Findings are in [MASS_OFFSET.md](MASS_OFFSET.md); retained evidence is under
[`scientific-results/mass-offset-20260919/`](scientific-results/mass-offset-20260919/).

## Requirements

- Linux with GNU `/usr/bin/time`
- Python 3.10 or newer
- Rust and Cargo
- `just`
- Enough disk space for two release builds

No Python packages are required.

For the hardening candidate, use `run_hardening.py` with explicit frozen baseline and candidate
binaries, the standard, modified, and feature configuration paths, an executable DuckDB path, and
a fresh output directory. It runs one warmup and three measured trials per engine and condition,
alternates engine order, records content hashes, and compares identities and exact prefilter output.

```shell
python3 benchmarks/run_hardening.py \
  --baseline /absolute/path/to/baseline/sage \
  --candidate /absolute/path/to/candidate/sage \
  --standard benchmarks/configs/local-standard.json \
  --modified benchmarks/configs/local-modifications.json \
  --feature /absolute/path/to/feature-config.json \
  --duckdb /absolute/path/to/duckdb \
  --output /absolute/path/to/fresh-hardening-results
```

Compare two existing Sage Plus runs with
`python3 benchmarks/compare_runs.py BASELINE_DIRECTORY CANDIDATE_DIRECTORY --duckdb DUCKDB_PATH --output REPORT_STEM`.
The JSON report retains exact additions, losses, and shared-identity score differences. Comparison
only covers stored rows, so choose matching output cutoffs. Peptidoform strings retain modification
mass annotations. This command does not estimate calibration.

If runtime review thresholds are exceeded, use `investigate_runtime.py` with explicit binary paths,
one or more `--config` paths from the frozen run, and a fresh `--output` directory. It repeats the
paired cases with one warmup and five measured trials by default, recording process CPU time and
host load. Preserve the initial flagged measurements alongside the follow-up. A variable shared
host cannot establish a reliable speedup from a small timing sample.

The repaired FDRBench driver accepts `--baseline-format` and `--candidate-format` explicitly,
plus `--subset`, `--mzml`, `--jar`, `--java`, `--duckdb`, `--time`, and `--prlimit`.
Use `--preflight-only` to check those dependencies before generation. For a small smoke run, set
`--seed-count 1` and provide a small FASTA subset. For paired Sage Plus builds, set both formats to
`parquet`. A complete run writes a schema-2 manifest. Analyze it using
`python3 benchmarks/analyze_fdrbench_validation.py --input RESULT_DIRECTORY`.
Paper data is written only when `--paper-output` is explicit.
Use distinct `--baseline-label` and `--candidate-label` values when comparing two Sage Plus builds.
The analyzer writes `analysis-method.json` to state its zero-discovery convention and the meaning
of its percentile ranges. These ranges describe seed variability and are not confidence intervals
for a study-level mean. Independently check threshold counts using
`python3 benchmarks/verify_entrapment_counts.py --input RESULT_DIRECTORY`.

Stage reuse requires matching content hashes for inputs and outputs. Failed stages never acquire
a completion manifest. Run only one process in each experiment directory. Legacy existence-only
results are historical evidence and are not accepted by the new analyzer. The bounded subset
experiment is a stress test and cannot establish calibration on a complete reference proteome.

## Recommended dataset

Use one representative DDA configuration with one local mzML file and the FASTA normally used for
that experiment. Keep input files on a local SSD. The configuration must be understood by both the
baseline and candidate when running the comparison benchmark.

This workspace has a suitable local workload. It uses a 219 MB HEK SILAC mzML and a 13 MB
reviewed-human FASTA. The dataset is ignored by Git and is not part of a fresh clone.

`benchmarks/configs/local-standard.json` provides the conventional database-search configuration
used for the baseline comparison. The two `upstream-compatible` configurations add common
modifications and a broad four-class PTM search using options shared by upstream Sage and Sage Plus.
`data/silac-k6r6/config.json` enables SILAC channels, LFQ, matched fragments, and spectral-library
export for a separate candidate feature check.

The one-spectrum test fixture must not be used for benchmark timing.

## Integrated memory design records

Two benchmark-backed design records document the compact database representations now integrated
into Sage Plus:

- [Protein-backed peptide sequences](PEPTIDE_INDEX_EXPERIMENT.md)
- [Lossless packed fragment index](FRAGMENT_INDEX_EXPERIMENT.md)

These records describe the alternatives considered, exactness checks, code tradeoffs, and the
benchmarks used before integration.

## Commands

List the available recipes:

```shell
just
```

Run the core suite:

```shell
just bench /absolute/path/to/config.json
```

Run the existing local HEK SILAC workload:

```shell
just bench-local
```

Run the primary upstream Sage versus Sage Plus comparison:

```shell
just bench-sage-comparison
```

This runs conventional, common-modification, and broad-PTM searches with identical inputs and
parameters for upstream Sage `v0.15.0-beta.2` and the current Sage Plus tree. It reports matched
whole-search wall time and peak resident memory.

Run the feature-heavy SILAC workload separately:

```shell
just bench-local-feature
```

Run the bounded variable-modification comparison:

```shell
just bench-local-mods
```

This workload enables methionine oxidation and peptide N-terminal acetylation. It allows at most
two variable modifications and four total variants per peptide, which expands the search space
without allowing an unbounded combinatorial search.

Run the conventional, exhaustive, and PTM site-library comparison:

```shell
just bench-ptm-library
```

The generator samples exact FASTA locations for phospho, oxidation, acetyl, and deamidation rules.
It creates nested libraries with 0, 50,000, 200,000, and 500,000 unique sites. Each configuration
allows up to three library-supported modifications per peptide and at most four peptide variants.
The recipe also measures a conventional search with no variable PTMs and an exhaustive search that
applies the same four PTM definitions at every compatible residue. The shared four-variant cap keeps
the exhaustive search inside the configured 24 GiB memory guard. The benchmark measures the complete
HEK search, including library parsing, database construction, scoring, modeling, and output writing.

Run the focused charge-aware preprocessing benchmark:

```shell
just bench-charge
```

This uses the existing deterministic `charge_matching_benchmark` example. It compares scored
deisotoping when fragment charges must be inferred with the path where charge arrays are supplied
by the input file. It reports preprocessing time, search time, PSM counts, matched peaks, and
deterministic checksums.

The core suite performs two comparisons. Both record complete-process wall time and peak RSS:

1. A normal search on the baseline and candidate
2. Candidate exact prefiltering with the feature off and on

Run an individual part:

```shell
just bench-search /absolute/path/to/config.json
just bench-prefilter /absolute/path/to/config.json
```

An optional generated pre-digested database benchmark isolates database construction. It defaults
to one million target peptides and their generated decoys:

```shell
just bench-memory
```

Run a candidate-only feature configuration:

```shell
just bench-feature /absolute/path/to/feature-config.json
```

A useful feature configuration can enable nonlinear RT alignment, RT prediction, LFQ, matched
fragment output, or any other feature under evaluation. The report copies the important fields
from `run-summary.json`, including model, quantification, and spectral-library counters when they
are present.

Check only the Python harness and Just recipes:

```shell
just bench-check
```

## Reproducibility controls

Environment variables can override the defaults:

```shell
BASELINE_REF=v0.1.0-beta.1 \
REPEATS=5 \
WARMUPS=1 \
THREADS=16 \
just bench /absolute/path/to/config.json
```

`PREFILTER_CHUNK_SIZE` controls the fixed prefilter chunk size. Defaults can also be placed in a
repository-root `.env` file because the Justfile enables dotenv loading.

Use `MEMORY_PEPTIDES=2000000 just bench-memory` to resize the optional generated database workload.

Every search uses the requested Rayon thread count, `--batch-size 1`, disabled telemetry, a fresh
output directory, and a release build made with `--locked`. Builds are completed before timing.

## Results

Reports are written under `benchmarks/results/<timestamp>-<suite>/`. Each report directory contains:

- `report.md`, the human-readable summary
- `records.json`, one record per measured trial
- `metadata.json`, machine, toolchain, commit, and dirty-tree details
- `configs/`, copies of generated benchmark configurations
- `runs/`, logs, timing data, summaries, and artifact hashes for each trial

Large analytical outputs are removed after hashing by default. Set `SAGE_BENCH_KEEP_OUTPUTS=1` to
retain them for inspection.

Peak memory is GNU Time's maximum resident set size. Search timing is external wall time, so it
includes database construction, input parsing, scoring, modeling, quantification, and output
writing. The report uses medians across measured trials.

The report labels threshold findings as `PASS` or `REVIEW`. A review is not automatically a bug.
Changes to preprocessing or scoring can intentionally alter identifications. The thresholds are
simple prompts for investigation:

- More than 1 percent fewer PSMs or peptides
- More than 10 percent slower wall time
- More than 10 percent higher peak RSS
- Different prefilter-off and prefilter-on result hashes or one-percent FDR counts

## Cleanup

Generated files are ignored by Git. Cleanup requires an explicit confirmation:

```shell
just bench-clean yes
```
