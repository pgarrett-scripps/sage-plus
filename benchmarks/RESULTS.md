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

### Upstream Sage versus Sage Plus

Three HEK workloads were run on September 1, 2026 with upstream Sage
`v0.15.0-beta.2` at commit `df921995` and Sage Plus beta.2 release code at commit
`5aeb568`. Each workload used identical input files and search parameters for
both binaries, eight threads, one warmup, and three measured trials.

| Workload | Sage wall | Sage Plus wall | Wall change | Sage peak RSS | Sage Plus peak RSS | RSS change |
|---|---:|---:|---:|---:|---:|---:|
| Conventional | 5.38 s | 6.45 s | +19.9% | 1,911.4 MiB | 1,494.8 MiB | -21.8% |
| Common modifications | 11.94 s | 12.14 s | +1.7% | 5,495.6 MiB | 3,698.5 MiB | -32.7% |
| Broad PTMs | 25.74 s | 29.07 s | +12.9% | 12,278.0 MiB | 7,792.3 MiB | -36.5% |

Sage Plus used less peak memory in every matched workload. Upstream Sage was
faster in every workload. Upstream Sage emitted its legacy TSV output and did
not provide the run-summary fields used for direct accepted-identification
comparisons. The raw reports are `benchmarks/results/20260901-142107-search/`,
`benchmarks/results/20260901-142158-search/`, and
`benchmarks/results/20260901-142340-search/`.

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

### Modification-specific PTM search comparison

This candidate-only workload used the conventional HEK mzML and reviewed-human FASTA. It includes
a conventional search with no variable PTMs and an exhaustive search that applies phospho,
oxidation, acetyl, and deamidation definitions at every compatible residue. A seeded reservoir
sample selected exact sites for nested libraries containing 0, 50,000, 200,000, or 500,000 of the
3,720,985 eligible FASTA residues. Exhaustive and site-library searches allowed up to three PTMs
and four total variants per peptide. A 64-variant exhaustive configuration failed preflight with an
estimated 117.83 GiB additional modified-peptide peak. An eight-variant configuration also exceeded
available memory. Every reported PTM condition was rerun with the shared four-variant cap.

| Condition | Median wall time | Median peak RSS | Database peptides | Fragments | PSMs at 1% FDR | Peptides at 1% FDR |
|---|---:|---:|---:|---:|---:|---:|
| Conventional | 6.07 s | 1,499.7 MiB | 2,803,578 | 80,833,042 | 2,203 | 1,409 |
| Exhaustive four-PTM | 24.96 s | 4,605.5 MiB | 10,787,909 | 317,149,388 | 2,267 | 1,497 |
| Site library, 0 | 12.25 s | 1,584.1 MiB | 2,803,578 | 80,833,042 | 2,203 | 1,409 |
| Site library, 50,000 | 12.72 s | 1,681.8 MiB | 3,030,517 | 89,341,256 | 2,206 | 1,412 |
| Site library, 200,000 | 14.53 s | 2,085.5 MiB | 3,775,367 | 117,638,808 | 2,203 | 1,414 |
| Site library, 500,000 | 17.72 s | 2,826.1 MiB | 5,271,114 | 172,112,992 | 2,199 | 1,420 |

Against exhaustive enumeration, the 500,000-site library reduced wall time by 29.0 percent, peak
RSS by 38.6 percent, database peptides by 51.1 percent, and theoretical fragments by 45.7 percent.
From the empty library to 500,000 sites, wall time increased by 44.7 percent and peak RSS increased
by 78.4 percent. Every condition produced one stable result hash across its three measured trials.
Identification counts are descriptive because each condition searches a different modified-peptide
space. The raw reports are `20260901-140801-feature`, `20260901-140156-feature`,
`20260901-140350-feature`, `20260901-140443-feature`, `20260901-140538-feature`, and
`20260901-140640-feature`.

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
