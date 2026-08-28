# Sage Plus developer benchmark results

These results are a local release check, not a published or comprehensive performance study. They
use one HEK SILAC mzML file and a reviewed-human FASTA from the untracked local `data/` directory.
See [README.md](README.md) for the repeatable commands and report format.

## Beta.2 release candidate

- CPU: Intel Core i7-10700K at 3.80 GHz
- Logical CPUs: 16
- Benchmark threads: 8
- Rust: 1.97.1
- Baseline: `v0.1.0-beta.1` at `e85a193`
- Candidate: `95fbbab`, the final integrated code before the release metadata commit
- Trials: one warmup followed by three measured runs

### Conventional database search

The conventional workload searches a 219 MB mzML against 2,803,578 target and decoy peptides and
80,833,042 theoretical fragments.

| Version | Median wall time | Median peak RSS | PSMs at 1% FDR | Peptides at 1% FDR |
|---|---:|---:|---:|---:|
| Beta.1 | 6.60 s | 2,123.9 MiB | 2,209 | 1,409 |
| Beta.2 candidate | 6.01 s | 1,497.5 MiB | 2,203 | 1,409 |

The candidate was 8.9 percent faster and used 29.5 percent less peak RSS. Peptide counts were
identical. The six-PSM difference is 0.3 percent and reflects the intentional scoring and
preprocessing changes described in the changelog.

### Exact prefiltering

Both modes used the beta.2 candidate and produced byte-identical `results.sage.parquet` files in
every measured trial.

| Mode | Median wall time | Median peak RSS | Database peptides | PSMs at 1% FDR | Peptides at 1% FDR |
|---|---:|---:|---:|---:|---:|
| Off | 6.12 s | 1,501.3 MiB | 2,803,578 | 2,203 | 1,409 |
| On | 11.11 s | 613.4 MiB | 1,031,633 | 2,203 | 1,409 |

Exact prefiltering reduced peak RSS by 59.1 percent and retained identical results. Its wall-time
cost was 81.5 percent on this workload.

### Bounded variable-modification search

This workload adds variable methionine oxidation and peptide N-terminal acetylation. It generated
7,211,871 target and decoy peptides and 219,325,108 theoretical fragments.

| Version | Median wall time | Median peak RSS | PSMs at 1% FDR | Peptides at 1% FDR |
|---|---:|---:|---:|---:|
| Beta.1 | 16.87 s | 5,463.5 MiB | 2,257 | 1,456 |
| Beta.2 candidate | 13.84 s | 3,339.5 MiB | 2,256 | 1,458 |

The candidate was 18.0 percent faster and used 38.9 percent less peak RSS. It returned one fewer
PSM and two more peptides at one-percent FDR.

### Feature-heavy SILAC search

The candidate-only workload enabled two SILAC channels, LFQ, matched fragments, retention-time
modeling, and spectral-library export. Its database contained 5,591,716 peptides and 161,210,114
theoretical fragments.

| Median wall time | Median peak RSS | PSMs at 1% FDR | Peptides at 1% FDR | LFQ features | Library entries | Library transitions |
|---:|---:|---:|---:|---:|---:|---:|
| 15.62 s | 2,593.8 MiB | 3,266 | 1,524 | 3,434 | 2,536 | 35,131 |

All three trials produced identical hashes for results, LFQ, matched fragments, Sage Parquet
library output, and PSI mzSpecLib output.

### Synthetic database memory

The generated pre-digested workload contained 1,999,980 target and decoy peptides and 47,999,520
theoretical fragments.

| Version | Median wall time | Median peak RSS |
|---|---:|---:|
| Beta.1 | 3.38 s | 1,049.0 MiB |
| Beta.2 candidate | 3.03 s | 790.1 MiB |

The candidate was 10.4 percent faster and used 24.7 percent less peak RSS.

### Charge-aware preprocessing

The deterministic synthetic scorer benchmark contains 25,000 peptides and repeatedly searches 160
spectra.

| Fragment charge source | Median wall time | Peak RSS | Preprocessing per spectrum | Search per spectrum | PSMs | Deterministic |
|---|---:|---:|---:|---:|---:|---|
| Inferred | 0.95 s | 24.2 MiB | 10.77 us | 181.14 us | 3,300 | yes |
| Supplied array | 0.91 s | 24.2 MiB | 7.28 us | 179.96 us | 3,300 | yes |

Supplied fragment charges reduced preprocessing time by 32.4 percent. Both modes produced stable
checksums and PSM counts across all three trials.

## Integrated compact database indexes

Target peptides share immutable source-protein storage where possible, as documented in
[PEPTIDE_INDEX_EXPERIMENT.md](PEPTIDE_INDEX_EXPERIMENT.md). The preliminary fragment index stores
lossless six-byte records and caps every search bucket at `database.bucket_size`, as documented in
[FRAGMENT_INDEX_EXPERIMENT.md](FRAGMENT_INDEX_EXPERIMENT.md).

## Raw local reports

When present in this workspace, the generated release-candidate reports are:

- `benchmarks/results/20260828-111601-all/report.md`
- `benchmarks/results/20260828-111814-search/report.md`
- `benchmarks/results/20260828-112043-feature/report.md`
- `benchmarks/results/20260828-112221-charge/report.md`
- `benchmarks/results/20260828-112241-memory/report.md`

The `benchmarks/results/` directory is ignored by Git because it contains machine-specific logs and
generated artifacts.
