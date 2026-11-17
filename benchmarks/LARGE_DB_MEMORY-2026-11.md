# Human + 100x gut catalog: peak memory (2026-11-16)

Peak resident memory of a prefiltered search of a human reference plus a 100x gut-microbiome (IGC)
catalog, on Sage Plus main at c8219ab (reports 0.1.0-beta.21) and upstream Sage 4a051b2. Same
FASTA (sha256 `b2ff6a66…`), same mzML and configuration as the
[head-to-head](HEADTOHEAD-2026-10.md). Peak RSS is from `/usr/bin/time -v`. These were memory
runs: not exclusive, so no wall time is quoted.

| Run | Peak RSS | Kept peptides | PSMs at 1% |
|---|---:|---:|---:|
| Sage Plus, `prefilter_min_matched_peaks` 4 (default), 8 threads | 22.21 GiB | 10,548,959 | 15,063 |
| Sage Plus, default, 16 threads | 21.59 GiB | 10,548,959 | 15,063 |
| Sage Plus, `prefilter_min_matched_peaks` 6, 16 threads | 20.27 GiB | 497,364 | 14,813 |
| Sage Plus, default, 16 threads, `MALLOC_ARENA_MAX=2` | 19.07 GiB | 10,548,959 | 15,063 |
| upstream, `prefilter_low_memory: true`, 16 threads | 15.20 GiB | 2,011,050 | 14,822 |
| upstream, `prefilter_low_memory: true`, 8 threads | 15.06 GiB | 2,011,050 | 14,822 |

Sage Plus's peak comes at the end of the protein stream, when the shared digests of every protein
are gathered and regrouped; until then its resident memory stays below upstream's. The 19 GiB
quoted from Beta 11 was a different run (an MGF input, a pre-release binary reporting beta.10,
8 threads, 19.19 GiB), not a measurement of the current code.
