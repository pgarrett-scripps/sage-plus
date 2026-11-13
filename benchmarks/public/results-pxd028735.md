# Receipts: five PXD028735 files, 2026-11-12

Every number in the README's results table, with the evidence behind it. The
three arms are described there; this file is the raw record.

## Machine and protocol

- 20 cores, 62 GB RAM, Linux 7.0, NVMe storage, `performance` CPU profile.
- A scheduler reserved the whole machine for the timing runs: it waited for the
  machine to go quiet, took every core, ran the three arms back to back inside
  that one reservation, and recorded what else touched the CPU meanwhile.
- 20 threads per arm (`RAYON_NUM_THREADS`), one run per arm, no repeats.
- The five conversions ran first, one file per job, outside the timing runs, so
  no conversion shares a measurement with a search.
- Both binaries were built with the same toolchain in release mode: Sage Plus
  0.1.0-beta.19 (commit 32a5ca5, reading `.raw` through `sage-plus-raw` 0.2.0),
  upstream Sage v0.15.0-beta.2 from its tag.
- FASTA: human, yeast and *E. coli* reference proteomes plus iRT peptides, the
  same file for all three arms.

## Conversion, per file

ThermoRawFileParser 1.4.5, digest-pinned container, `--format=1`, one job per
file, nothing else of ours on the machine.

| File | Seconds |
| --- | ---: |
| `LFQ_Orbitrap_DDA_Condition_A_Sample_Alpha_01` | 251.17 |
| `LFQ_Orbitrap_DDA_Condition_A_Sample_Beta_01` | 264.21 |
| `LFQ_Orbitrap_DDA_Condition_A_Sample_Gamma_01` | 259.46 |
| `LFQ_Orbitrap_DDA_Condition_B_Sample_Alpha_01` | 250.21 |
| `LFQ_Orbitrap_DDA_Condition_B_Sample_Beta_01` | 263.52 |
| **Total** | **1,288.57** |

Each conversion job peaked at 1.7–1.9 GB of resident memory and took 4m10s to
4m24s of wall time, including the container start.

## Searches

| | `plus-raw` | `plus-mzml` | `upstream-mzml` |
| --- | ---: | ---: | ---: |
| Binary version | sage 0.1.0-beta.19 | sage 0.1.0-beta.19 | sage 0.15.0-beta.2 |
| Input | 5 `.raw` | 5 mzML | 5 mzML |
| Conversion seconds | 0 | 1,288.57 (recorded) | 1,288.57 (recorded) |
| Search seconds | 61.87 | 64.23 | 101.69 |
| Total seconds | 61.87 | 1,352.80 | 1,390.26 |
| Peak RSS (GiB) | 7.810 | 7.413 | 11.159 |
| PSMs at 1% FDR | 318,226 | 318,214 | 311,016 |
| Peptides at 1% FDR | 60,429 | 60,441 | 59,051 |
| Proteins at 1% FDR | 7,034 | 7,034 | 6,977 |
| Exit code | 0 | 0 | 0 |

Search seconds are `/usr/bin/time -v` wall clock around the search process;
peak RSS is its `Maximum resident set size`. The three identification counts
are computed from each arm's own `results.sage.parquet` by one function:
targets with `spectrum_q` ≤ 0.01, and distinct `peptide` and `proteins` values
among rows whose `peptide_q` or `protein_q` is ≤ 0.01.

## What the scheduler recorded during the timing runs

The three arms ran inside one exclusive reservation, so there is one record for
all of them.

| | Recorded |
| --- | --- |
| Verdict | clean |
| Wall time of the reservation | 231.9 s for the three searches |
| Foreign load | 1.33 cores (307.9 CPU-s not attributable to the job; threshold 4.00) |
| CPU clocks | busy cores averaged 2,716 MHz (lowest 400) against a usual 2,853 MHz, 92,416 core / 32,729 package throttle events, peak 100 °C, `performance` profile |

The run was not judged contended, and the three arms ran back to back under the
same load, so they are comparable with each other. They are not comparable with
numbers from another machine: this one throttled at 100 °C throughout, which
costs every arm some clock speed, and the difference of 2.4 s between the two
Sage Plus arms is smaller than that effect.

## Reproducing it

```shell
for raw in /data/PXD028735/*.raw; do
  ./run_public_benchmark.sh --arm convert --out runs/convert --mzml-cache runs/mzml "$raw"
done
for arm in plus-raw plus-mzml upstream-mzml; do
  ./run_public_benchmark.sh --arm "$arm" --out "runs/$arm" --mzml-cache runs/mzml \
    --fasta "$SAGE_BENCH_FASTA" /data/PXD028735/*.raw
done
```

`SAGE_BENCH_PLUS`, `SAGE_BENCH_UPSTREAM` and `SAGE_BENCH_THREADS` select the
binaries and the thread count; see the README.
