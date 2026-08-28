# Lossless fragment index experiment

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

## Selected prototype

The selected design constructs 6-byte packed fragments directly. It never
creates the original 8-byte fragment vector.

Each record stores:

- A 32-bit peptide index
- The low 8 to 16 bits of the exact `f32` mass representation

Each mass bucket stores the shared high mass bits once. The configured bucket
size selects the suffix width. A value of 16,384 selects a 12-bit suffix, which
matches the benchmarked layout. No floating-point quantization occurs. Joining
the shared prefix and stored suffix recreates the original `f32` bit pattern.

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
| Median peak RSS | 1,887.2 MiB | 1,506.2 MiB | 20.2% lower |
| Median wall time | 6.49 s | 6.12 s | 5.7% faster |
| PSMs at 1% FDR | 2,203 | 2,203 | Identical |
| Peptides at 1% FDR | 1,409 | 1,409 | Identical |

All measured Parquet result hashes were
`819f8983bb951f74be55690517434023ebcbd17f0fd1d8a2a42dfa8890f32c33`.

The focused deterministic scorer benchmark also produced identical peptide,
matched-peak, and score checksums. Its median search time was effectively
unchanged for the bucket-level prototype. The final direct-build layout was
substantially faster in that synthetic workload because its finer mass-prefix
buckets rejected irrelevant fragments earlier.

## Remaining validation

Before integration, benchmark the packed index with:

- Open precursor searches
- Semi-enzymatic and non-enzymatic databases
- Variable modifications and larger peptide counts
- Prefilter enabled and disabled
- Narrow DDA and wide DIA isolation windows
- Different configured bucket sizes

The implementation also uses a small audited unsafe section to fill pre-counted
parallel output ranges. Unit tests verify the 6-byte record size, exact mass
round trips, ordering, range selection, and exhaustive query visitation.
