# Spectrum-indexed prefilter benchmark

Beta 7 replaces the chunked prefilter. Beta 6 built a fragment index for every database chunk
and searched every spectrum against each one, rereading spectra for every chunk when they did not
fit in memory together. Beta 7 indexes the spectra once and streams each generated peptide
through that index. A peptide is kept when one of its preliminary fragments matches a peak of a
spectrum whose precursor window contains the peptide mass. This is the same rule as Beta 6, and
every window is computed with the same floating-point operations, so both keep the same peptides.

## Method

`benchmarks/run_prefilter.py` ran each configuration three times with the Beta 6 executable
(prefilter on), the candidate (prefilter on), and the candidate with the prefilter off. Engine
order alternated between repeats. Runs used one 16-thread workstation with `max_memory_gb` 28.
Evidence is retained under `/data/sage-plus-scientific/prefilter-20260923/`.

| Executable | SHA-256 prefix |
|---|---|
| Beta 6 release build | `5be3af33c4c8` |
| Candidate | `f3cc86354009` |

The candidate was built from the Beta 7 branch before the version bump and the timsTOF mobility
calibration. Neither change affects these non-Bruker searches.

| Workload | Configuration | Database peptides without prefilter |
|---|---|---:|
| HEK, oxidation and acetylation | `prefilter-closed-mods.json` | 7,211,871 |
| HEK, seven variable modifications | `prefilter-broad-ptm.json` | 18,386,093 |
| HEK, seven modifications, two per peptide | `prefilter-broad-ptm-2.json` | 70,525,438 generated |
| HEK, open search from -150 to +500 Da | `prefilter-open.json` | 7,211,871 |
| Five PXD028735 LFQ files, 567,401 spectra | `prefilter-lfq-5file.json` | 10,029,967 |
| PXD000138 HCD_1, phosphorylation mass offset | `prefilter-offset-hcd1.json` | 203,040 |
| PXD000138 HCD_2, phosphorylation mass offset | `prefilter-offset-hcd2.json` | 203,040 |

## Results

Median of three runs.

| Workload | Beta 6 prefilter | Beta 7 prefilter | Speedup | No prefilter | Retained peptides |
|---|---:|---:|---:|---:|---:|
| HEK, oxidation and acetylation | 15.6 s | 11.6 s | 1.3x | 11.5 s | 2,691,078 |
| HEK, seven variable modifications | 33.0 s | 24.1 s | 1.4x | 27.1 s | 6,545,770 |
| HEK, seven modifications, two per peptide | 105.0 s | 67.4 s | 1.6x | rejected | 20,850,246 |
| HEK, open search | 65.8 s | 54.6 s | 1.2x | 47.7 s | 7,195,546 |
| Five LFQ files | 613.4 s | 135.2 s | 4.5x | 102.9 s | 8,988,649 |
| HCD_1 mass offset | 35.7 s | 0.9 s | 40x | 0.8 s | 15,906 |
| HCD_2 mass offset | 41.0 s | 1.0 s | 41x | 0.8 s | 29,852 |

| Workload | Beta 6 prefilter | Beta 7 prefilter | No prefilter |
|---|---:|---:|---:|
| HEK, oxidation and acetylation | 1.20 GiB | 1.36 GiB | 3.41 GiB |
| HEK, seven variable modifications | 2.74 GiB | 2.93 GiB | 8.53 GiB |
| HEK, seven modifications, two per peptide | 9.11 GiB | 9.26 GiB | rejected |
| HEK, open search | 3.09 GiB | 3.24 GiB | 3.40 GiB |
| Five LFQ files | 6.37 GiB | 6.27 GiB | 6.81 GiB |
| HCD_1 mass offset | 0.09 GiB | 0.08 GiB | 0.16 GiB |
| HCD_2 mass offset | 0.08 GiB | 0.10 GiB | 0.16 GiB |

The prefilter stage itself, from database chunking to the retained set, fell from 10 s to 6 s,
21 s to 12 s, 78 s to 40 s, 19 s to 9 s, 517 s to 40 s, and about 35 to 40 s to under 1 s.

Without a prefilter, the two-modification search was rejected by the memory preflight, which
estimated 110.8 GiB for the modified database.

## Output agreement

- Beta 6 and Beta 7 retained the same number of peptides in every workload.
- `results.sage.parquet` was byte-identical between Beta 6, Beta 7, and the unfiltered search in
  every repeat of the closed, seven-modification, open, and HCD_2 workloads, and between Beta 6
  and Beta 7 in every repeat of the two-modification workload.
- HCD_1 differed in one spectrum, scan 8446, whose top hit is a tied mass-offset placement of
  `AIY[Phospho]FLSLLK`. The tie resolves differently between repeats of the Beta 6 executable
  too. Every other row matched, and all runs reported 3,141 PSMs at 1% FDR.
- The five-file LFQ search never produces identical files, even between two Beta 6 runs. Every
  search column (hyperscore, matched peaks, candidate counts, Poisson score, and the rest) agreed
  on every shared row. Only the LDA discriminant and the q-values derived from it drift, and 10 to
  21 PSMs near the 0.05 output cutoff move in or out. PSMs at 1% FDR ranged from 312,533 to
  312,538 across Beta 6 repeats and from 312,533 to 312,536 across Beta 7 repeats.

## Spectrum batches

With `SAGE_PREFILTER_INDEX_GB=0.1`, the five-file search built three spectrum indexes (243, 188,
and 101 MiB), streamed the database once per batch, and retained the same 8,988,649 peptides.
It took 140.2 s, compared with 135.2 s for a single index.

## Interpretation

- Beta 7 prefiltering is faster than Beta 6 prefiltering in every workload. The gain is largest
  when Beta 6 reread spectra for every chunk (five files) or searched many spectra against many
  small chunks (mass offsets with 200-protein chunks).
- Prefiltering still helps only when it removes much of the database. It kept 37% of the closed
  database and 36% of the seven-modification database, and matched or beat the unfiltered search
  with 60 to 66% less memory. It kept 99.8% of the open-search database and 90% of the five-file
  database, where prefiltering made the search 14 to 31% slower.
- Peak memory rose by 2 to 13% relative to Beta 6 prefiltering on the HEK searches and fell by
  2% on the five-file search. The mass-offset searches stayed at or below 0.1 GiB.
- These are single-machine measurements on public benchmark data, not a general performance
  claim.
