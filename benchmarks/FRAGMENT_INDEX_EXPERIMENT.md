# Lossless fragment index experiment

Status: integrated into `main` by commits `4b376b9` and `66f3b45`.

## Objective

Reduce theoretical-fragment memory without changing fragment masses, peptide
identities, candidate selection, scores, or reported results.

The comparison baseline is commit `6a686d0` on `main`.

## Baseline representation

Each fragment used an 8-byte record containing a 32-bit peptide index and a
32-bit floating-point mass. The standard HEK benchmark generated 80,833,042
records, giving a 646,664,336-byte raw fragment array.

## Prototypes rejected

The first analysis measured exact block entropy. With 128-fragment blocks, a
fully bit-packed representation projected 3.629 bytes per fragment. It was not
practical in the query path. The focused scorer benchmark was 227 percent
slower because precursor searches had to decode peptide IDs before rejecting
irrelevant fragments.

Direct 16-bit lanes reduced the decoding work but remained approximately 157
percent slower when both fields were compressed. Keeping peptide IDs direct
and bit-packing only masses remained approximately 124 percent slower. Direct
16-bit mass lanes reduced the regression to approximately 49 percent.

A bucket-level 16-bit lane format removed the search regression, but whole
buckets often crossed the 16-bit suffix boundary. Its growing output vectors
also coexisted with the original raw index. The resulting allocation used
11.624 bytes per fragment and increased median peak RSS by 22.4 percent. This
demonstrated that compression after raw-index construction cannot solve peak
memory.

## Integrated design

The selected design constructs 6-byte packed fragments directly. It never
creates the original 8-byte fragment vector.

Each record stores:

- A 32-bit peptide index
- The low 12 bits of the exact `f32` mass representation, stored in a `u16`

Each mass-prefix group stores the shared high mass bits once. Groups containing
more records than the configured `bucket_size` are split into multiple search
buckets with the same prefix. This restores `bucket_size` as a hard record
limit while retaining the benchmarked 12-bit suffix layout. No floating-point
quantization occurs. Joining the shared prefix and stored suffix recreates the
original `f32` bit pattern.

Construction uses two passes over theoretical fragments. The first pass counts
records for each mass prefix. The second pass writes into exact, disjoint output
ranges in parallel. Ranges are assigned by ascending peptide chunks, so each
bucket is already ordered by peptide index and requires no global fragment
sort.

## HEK benchmark

The standard benchmark used eight Rayon threads, one warmup, and three measured
trials.

| Metric | Baseline | Packed index | Change |
|---|---:|---:|---:|
| Fragment records | 80,833,042 | 80,833,042 | Identical |
| Fragment allocation | 646.7 MB raw | 485.1 MB total | 25.0% smaller |
| Median peak RSS | 1,910.3 MiB | 1,501.1 MiB | 21.4% lower |
| Median wall time | 6.50 s | 6.59 s | 1.4% slower |
| PSMs at 1% FDR | 2,203 | 2,203 | Identical |
| Peptides at 1% FDR | 1,409 | 1,409 | Identical |

All measured Parquet result hashes were
`819f8983bb951f74be55690517434023ebcbd17f0fd1d8a2a42dfa8890f32c33`.

The focused deterministic scorer benchmark also produced identical peptide,
matched-peak, and score checksums.

## Additional validation

The modification-heavy search generated 219,325,108 fragments. Peak RSS fell
from 4,528.3 MiB to 3,345.9 MiB, a 26.1 percent reduction, with identical wall
time and output hash.

An open precursor search reduced peak RSS by 22.2 percent. Runtime increased
by 3.0 percent, within the benchmark harness's five-percent noise threshold,
and the output hash remained identical. Small semi-enzymatic and non-enzymatic
searches also produced identical hashes and fragment counts. Enabling exact
prefiltering produced the same final output as leaving prefiltering disabled.

Workspace tests include exhaustive query visitation across generated bucket
sizes from 1 through 8,192. The standard benchmark exercises a bucket size of
16,384. Dedicated tests force repeated mass prefixes, verify that no bucket
exceeds the configured limit, and check queries beginning inside a repeated
prefix run.

The remaining performance workload is:

- Wide DIA isolation windows

The implementation also uses a small audited unsafe section to fill pre-counted
parallel output ranges. Unit tests verify the 6-byte record size, exact mass
round trips, ordering, range selection, and exhaustive query visitation.
