# Sage Plus benchmark pipeline

This directory contains a small repeatable benchmark harness for Sage Plus. It is intended for
development checks and README-sized performance summaries. It is not a scientific validation
suite.

The harness compares the current working tree with a pinned baseline. By default the baseline is
commit `a38387b`, the `origin/main` commit that preceded the current performance work. The candidate
build includes uncommitted working-tree changes.

## Requirements

- Linux with GNU `/usr/bin/time`
- Python 3.10 or newer
- Rust and Cargo
- `just`
- Enough disk space for two release builds

No Python packages are required.

## Recommended dataset

Use one representative DDA configuration with one local mzML file and the FASTA normally used for
that experiment. Keep input files on a local SSD. The configuration must be understood by both the
baseline and candidate when running the comparison benchmark.

This workspace has a suitable local workload. It uses a 219 MB HEK SILAC mzML and a 13 MB
reviewed-human FASTA. The dataset is ignored by Git and is not part of a fresh clone.

`benchmarks/configs/local-standard.json` provides the conventional database-search configuration
used for the baseline comparison. `data/silac-k6r6/config.json` enables SILAC channels, LFQ, matched
fragments, and spectral-library export for a separate candidate feature check.

The one-spectrum test fixture must not be used for benchmark timing.

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
