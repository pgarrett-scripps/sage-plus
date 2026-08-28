# Sage Plus developer benchmark results

These results are a local development check, not a published or comprehensive performance study.
They use one HEK SILAC mzML file and a reviewed-human FASTA from the untracked local `data/`
directory. See [README.md](README.md) for the repeatable commands and report format.

## Environment

- CPU: Intel Core i7-10700K at 3.80 GHz
- Logical CPUs: 16
- Benchmark threads: 8
- Rust: 1.97.1
- Baseline: `a38387bef254b19b2c4de89cf0847f4c16521221`
- Candidate: `84d7f50ffd95eb767cec5b445529e3c5b74ca107`
- Trials: one warmup followed by three measured runs

## Conventional database search

The conventional workload searches a 219 MB mzML against a database containing 2,803,578 target
and decoy peptides and 80,833,042 theoretical fragments.

| Version | Median wall time | Median peak RSS | PSMs at 1% FDR | Peptides at 1% FDR | Proteins at 1% FDR |
|---|---:|---:|---:|---:|---:|
| Baseline | 7.18 s | 2,087.6 MiB | 2,209 | 1,409 | 645 |
| Candidate | 6.99 s | 2,021.8 MiB | 2,203 | 1,409 | 645 |

The candidate was 2.6 percent faster and used 3.2 percent less peak RSS. It retained the same
peptide and protein counts, with six fewer spectrum-level PSMs at 1 percent FDR.

## Exact prefiltering

Both prefilter modes used the candidate binary and produced byte-identical
`results.sage.parquet` files in every measured trial.

| Mode | Median wall time | Median peak RSS | Database peptides | PSMs at 1% FDR | Peptides at 1% FDR |
|---|---:|---:|---:|---:|---:|
| Off | 6.94 s | 2,012.9 MiB | 2,803,578 | 2,203 | 1,409 |
| On | 13.38 s | 955.3 MiB | 1,031,633 | 2,203 | 1,409 |

Exact prefiltering reduced peak RSS by 52.5 percent and retained identical results. The tradeoff
was a 92.8 percent wall-time increase on this workload.

## Bounded variable-modification search

This workload adds variable methionine oxidation and peptide N-terminal acetylation. Each
modification has a maximum count of one, with at most two modifications and four generated variants
per peptide. The resulting database contains 7,211,871 target and decoy peptides and 219,325,108
theoretical fragments.

| Version | Median wall time | Median peak RSS | PSMs at 1% FDR | Peptides at 1% FDR | Proteins at 1% FDR |
|---|---:|---:|---:|---:|---:|
| Baseline | 16.33 s | 5,498.3 MiB | 2,257 | 1,456 | 646 |
| Candidate | 13.37 s | 4,784.2 MiB | 2,256 | 1,458 | 643 |

The candidate was 18.1 percent faster and used 13.0 percent less peak RSS. It returned one fewer
PSM, two more peptides, and three fewer proteins at 1 percent FDR. The database and fragment counts
were identical between versions.

The candidate's accepted spectrum-level results included 60 oxidation-bearing PSMs, 59
acetylation-bearing PSMs, and 7 PSMs containing both modifications. Oxidation and acetylation each
appeared on 49 distinct accepted peptidoforms.

## Feature-heavy SILAC search

The candidate-only feature workload enabled two SILAC channels, LFQ, matched fragments, RT
modeling, and spectral-library export. Its database contained 5,591,716 peptides and 161,210,114
theoretical fragments.

| Median wall time | Median peak RSS | PSMs at 1% FDR | Peptides at 1% FDR | LFQ features | Library entries | Library transitions |
|---:|---:|---:|---:|---:|---:|---:|
| 15.71 s | 3,715.8 MiB | 3,266 | 1,524 | 3,434 | 2,536 | 35,131 |

All three measured feature trials produced identical counts and hashes. One earlier exploratory
attempt terminated with signal 11 during database construction after logging an impossible memory
estimate. A later standalone run and the subsequent warmup plus three measured runs all completed.
That isolated failure remains unexplained and should be investigated before treating this workload
as fully stable.

## Charge-aware preprocessing

The deterministic synthetic scorer benchmark contains 25,000 peptides and repeatedly searches 160
spectra. It compares scored deisotoping with inferred fragment charges against supplied mzML-style
fragment charge arrays.

| Fragment charge source | Median wall time | Preprocessing per spectrum | Search per spectrum | Known charges | PSMs | Deterministic |
|---|---:|---:|---:|---:|---:|---|
| Inferred | 4.41 s | 10.20 us | 1,016.43 us | 400,000 | 3,300 | yes |
| Supplied array | 4.19 s | 7.33 us | 974.67 us | 495,900 | 3,300 | yes |

Supplied fragment charges reduced preprocessing time by 28.1 percent and search time by 4.1
percent. Each mode produced stable checksums across all three measured trials.

## Raw local reports

When present in this workspace, the generated reports are:

- `benchmarks/results/20260825-235007-all/report.md`
- `benchmarks/results/20260825-235328-feature/report.md`
- `benchmarks/results/20260825-235904-search/report.md`
- `benchmarks/results/20260826-000346-charge/report.md`

The `benchmarks/results/` directory is ignored by Git because it contains machine-specific logs and
generated artifacts.
