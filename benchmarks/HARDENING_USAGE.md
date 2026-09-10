# Hardening validation commands

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
