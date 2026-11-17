# Human + 100x gut catalog: peak memory (2026-11-16, updated 2026-11-17)

## After: shared digests tested as they stream (35a8495, 2026-11-17)

The prefilter now tests each copy of a shared digest against the spectrum index while streaming
and gathers only the shared sequences some copy could match (748,054 of 16,155,649 keys;
2,715,153 digests instead of 55,235,156). Same FASTA (sha256 `b2ff6a66…`, rebuilt from
`GUT_DB1IGC.faa` and UniProt human reviewed 2026-11), same configuration, 293T centroid mzML.
Release builds of main 2cc84ae (before) and 35a8495 (after), 16 threads, peak RSS from
`/usr/bin/time -v` with the RSS sampler (ajs 3406):

| Run | Peak RSS | Kept peptides | PSMs at 1% | Results rows |
|---|---:|---:|---:|---:|
| before (2cc84ae), default, 16 threads | 22.21 GiB | 10,548,959 | 15,063 | 20,669 |
| after (35a8495), default, 16 threads | 6.85 GiB | 10,548,959 | 15,063 | 20,669 |

Every output (results, run summary, QC and digestion tables) is identical between the two.
Exclusive timing, 20 threads, three rotated rounds (ajs 3407; not contended, not throttled; mean
core clock 2847 MHz against a usual 2640, peak 100 °C, platform profile performance):

| Round | before | after |
|---|---:|---:|
| 1 | 5:08.5 (ran first) | 4:23.7 |
| 2 | 4:31.5 | 5:07.3 (ran first) |
| 3 | 4:30.7 (ran first) | 4:24.6 |

The first run of rounds 1 and 2 was slow for whichever build ran it (cold page cache). Medians:
4:31.5 before, 4:24.6 after. A normal human prefilter search (ajs 3404, five rounds) took 29.30 s
before and 29.70 s after on average, within its 28.8-30.5 s spread, at the same 1.2 GiB.

## Before (c8219ab, 2026-11-16)

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

Before the change, Sage Plus's peak came at the end of the protein stream, when the shared digests of every protein
were gathered and regrouped; until then its resident memory stayed below upstream's. The 19 GiB
quoted from Beta 11 was a different run (an MGF input, a pre-release binary reporting beta.10,
8 threads, 19.19 GiB), not a measurement of the current code.
