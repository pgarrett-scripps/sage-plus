# Protein-backed peptide sequence experiment

## Scope

This experiment replaces independently allocated target peptide sequences with immutable spans into shared protein storage. A peptide sequence is represented by a thin `Arc` plus `u32` start and end offsets. The value remains 16 bytes, the same size as `Arc<[u8]>`, so `Peptide` remains 144 bytes.

Sequence equality, ordering, and hashing use residue contents. Storage identity and offsets are never part of peptide identity. Identical sequences at repeated positions or in different proteins therefore continue to group and deduplicate together. Protein occurrence metadata still records every source location.

Generated decoys use owned sequence storage because they do not occur in the source protein. All modification variants produced from one unmodified peptide share the same target storage and the same reversed decoy storage. Pre-digested TSV sequences also use owned storage.

## Benchmark method

The baseline executable was captured from the current working state before the experiment. Both versions were built in release mode. Each reported result is the median of three runs with 16 Rayon threads. GNU Time measured whole-process wall time and peak resident memory.

The digestion and peptide-expansion benchmark used the reviewed human FASTA from the local HEK workload. The tryptic case used all 20,416 proteins. The semi-enzymatic case used a deterministic 1,000-protein prefix. The non-specific case used a deterministic 300-protein prefix with peptide lengths from 7 through 30.

### Unmodified digestion and grouping

| Mode | Groups | Baseline RSS | Span RSS | RSS change | Baseline wall | Span wall |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| Tryptic | 1,426,203 | 427.8 MiB | 338.7 MiB | -20.8% | 1.17 s | 0.89 s |
| Semi-enzymatic | 1,539,536 | 439.4 MiB | 320.2 MiB | -27.1% | 0.78 s | 0.51 s |
| Non-specific | 4,299,035 | 1,212.6 MiB | 716.9 MiB | -40.9% | 3.15 s | 1.38 s |

Group counts, origin counts, residue counts, and checksums were identical for every pair.

### Expanded target and decoy peptides

| Mode | Peptides | Baseline RSS | Span RSS | RSS change | Baseline wall | Span wall |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| Tryptic | 2,810,086 | 1,163.4 MiB | 1,086.5 MiB | -6.6% | 2.88 s | 2.68 s |
| Semi-enzymatic | 3,028,031 | 1,227.0 MiB | 1,119.3 MiB | -8.8% | 2.58 s | 2.21 s |
| Non-specific | 8,595,484 | 3,541.6 MiB | 3,082.5 MiB | -13.0% | 9.41 s | 6.45 s |

Peptide counts, decoy counts, residue counts, protein links, and checksums were identical for every pair.

### Full HEK search

The conventional full search used the 20,416-protein FASTA and 219 MB HEK SILAC mzML. It generated 2,803,578 peptides and 80,833,042 fragments.

| Search | Baseline RSS | Span RSS | RSS change | Baseline wall | Span wall | Runtime change |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| Conventional | 2,036.7 MiB | 1,943.1 MiB | -4.6% | 6.01 s | 5.67 s | -5.7% |
| Exact prefilter | 966.6 MiB | 890.9 MiB | -7.8% | 10.54 s | 9.80 s | -7.0% |

The baseline and span implementations produced byte-identical result Parquet files in both search modes. The conventional search reported the same 2,203 target PSMs, 1,409 target peptides, 645 target proteins, and 686 target protein groups at one percent FDR.

## Code impact

Relative to the current working state, the experiment changes 14 files with 454 insertions and 46 deletions. The main cost is one 314-line sequence module that owns conversion, byte-slice access, content-based comparison, hashing, formatting, and storage-sharing tests. Production call-site changes are small and concentrated in FASTA parsing, digestion, peptide conversion, decoy generation, and prefilter target collision tracking.

The main compatibility cost is public API churn. `Fasta::targets`, `Digest::sequence`, and `Peptide::sequence` now expose wrapper types instead of `String` or `Arc<[u8]>`. The peptide wrapper dereferences to `[u8]` and provides common conversions, but downstream code that constructs these structs directly may need `.into()`.

## Recommendation

Integrate the approach if memory-efficient semi-enzymatic and non-specific search is an important project direction. The end-to-end tryptic saving is modest but measurable, while the digestion and expanded-search improvements are large enough to change practical memory limits. The implementation keeps peptide identity correct, does not enlarge `Peptide`, improves runtime, and has no observed output changes.

Before merging, treat the public sequence type change as the main review item. If library API stability outweighs the measured memory gains, keep this branch as a proven design and defer integration until a breaking release. The internal design itself does not need a global protein index or database lifetime coupling, which keeps it substantially simpler than storing numeric protein IDs in every peptide.
