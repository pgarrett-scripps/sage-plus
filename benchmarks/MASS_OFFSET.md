# Mass offset search evaluation

Evaluation of `search_mode: "mass_offset"` (searching a modification as a precursor and
fragment offset) against the equivalent database expansion. Run on 2026-09-19 with
`claude/mass-offset-search` and the unmodified `0bfd5be` binary as the reference.

Artifacts, configurations, logs, and analysis scripts:
`/mnt/data1/sage-plus-mass-offset-eval` (run_eval.py, analyze.py, runs/, cascade/).
Spectra and truth come from the existing scientific corpus under
`/data/sage-plus-scientific/20260914`; identity of those inputs is recorded there and is
not re-established here.

## What was compared

| Suite | Data | Question |
| --- | --- | --- |
| Phospho | PXD000138 synthetic phosphopeptide libraries, HCD_1 and HCD_2, synthesis-defined sites | Does offset search reproduce an expanded phospho search, including localization? |
| HEK | HEK SILAC-K6R6 with reviewed human FASTA, oxidation variable, LFQ and localization enabled | Behavior on a real proteome search with quantification |
| Scale | Same, oxidation and phospho together | Feasibility when expansion is too large |
| Cascade | PXD000138 HCD_1 through Cascade discovery, native library, guided search | Compatibility with the Cascade workflow |

## Regression: configurations without offsets are unchanged

`results.sage.parquet` is byte-identical (checksum over all reported rows) between the
pre-change binary and this branch for every non-offset configuration tested:
HCD_1 (3142 rows), HCD_2 (1959 rows), and HEK (3957 rows).

## Phospho: offset vs expanded database

| Run | Indexed peptides | Fragments | PSMs @1% | Peak RSS | Search wall | Correct sites | Incorrect sites |
| --- | --- | --- | --- | --- | --- | --- | --- |
| HCD_1 database | 885,456 | 28,091,052 | 3141 | 0.45 GB | 2.1 s | 1172 | 118 |
| HCD_1 offset | 203,040 | 6,006,480 | 3141 | 0.16 GB | 1.3 s | 1171 | 118 |
| HCD_1 offset + prefilter | 15,906 | 231,236 | 3141 | 0.09 GB | 45.1 s | 1171 | 118 |
| HCD_2 database | 885,456 | 28,091,052 | 1958 | 0.45 GB | 1.9 s | 880 | 15 |
| HCD_2 offset | 203,040 | 6,006,480 | 1958 | 0.17 GB | 1.2 s | 880 | 15 |
| HCD_2 offset + prefilter | 29,852 | 470,868 | 1958 | 0.08 GB | 52.4 s | 880 | 15 |

Site counts are events at 1% spectrum q and 1% localization q, scored against the
synthesis-defined sites by `benchmarks/scientific_metrics.localization_metrics`. The
empirical site error fraction is identical in both modes (0.0915 for HCD_1, 0.0168 for
HCD_2); it includes identification and localization error and is not an arrangement-level
FLR estimate. The single differing HCD_1 site event comes from the one differing PSM
below.

Rank-1 agreement, offset against database:

- HCD_1: 3141 of 3142 spectra report the same peptidoform, all with identical hyperscore.
  The exception is an exact score tie between two isobaric sequences
  (`ALLSLIT[Phospho]FK` and `ALLSLMYFK`), where either answer is equally supported.
- HCD_2: 1959 of 1959 identical.
- `results.sage.ptm-library.tsv`, the file Cascade consumes, is byte-identical between
  modes for both files.

Prefilter on versus off agrees on every spectrum in both files, so offset-aware
retrieval survives chunked prefiltering. Prefiltering is much slower here only because
the synthetic FASTA was chunked at 200 proteins.

## HEK: real proteome search with quantification

| Run | Indexed peptides | PSMs @1% | Peptides @1% | Proteins @1% | Peak RSS | Wall |
| --- | --- | --- | --- | --- | --- | --- |
| database (oxidation) | 7,211,871 | 2256 | 1458 | 643 | 3.43 GB | 13.2 s |
| offset (oxidation) | 5,602,921 | 2267 | 1465 | 644 | 2.67 GB | 13.0 s |
| offset + prefilter | 5,602,921 | 2267 | 1465 | 644 | 2.70 GB | 12.8 s |

3909 of 3986 rank-1 spectra agree. The two searches are not equivalent here by
construction: this configuration caps expansion with `max_combinations: 4` and allows
`max_variable_mods: 2`, so the database search generates up to two oxidations but drops
variants beyond the cap, while the offset search places exactly one oxidation on every
peptide including those whose oxidized variant the cap removed. The offset search reports
slightly more identifications as a result.

Offset peptidoforms reach quantification: `lfq.parquet` holds 58 oxidized precursor rows
under offset search against 53 under expansion, with distinct precursor identities rather
than merged ones.

## Scale: two modifications at once

Oxidation and phospho together, human FASTA, `max_variable_mods: 2`, no combination cap:

- Database expansion is refused before searching: estimated modified-peptide peak of
  65.15 GiB against the configured 14 GiB ceiling.
- Offset search completed in 16 s with a 2.69 GB peak, 5,602,921 indexed peptides, and
  2296 PSMs / 1479 peptides / 652 proteins at 1% FDR, with 1177 distinct offset
  peptidoforms.

## Cascade compatibility

Cascade drives Sage Plus through its CLI and consumes the native PTM site library, so the
full workflow was run end to end with phospho declared as a mass offset:

- `cascade run` discovery produced 3029 target PSMs at 1% and a 406-site native library
  with protein coordinates, which is what iteration and guided search consume.
- `cascade guided`, which sets `site_mode: "library"` for every definition, produced 2960
  PSMs at 1% (1484 phospho). Against the same guided search with phospho as a database
  modification: 2960 of 2961 rank-1 spectra identical, every hyperscore identical, and a
  smaller index (203,040 vs 203,852 peptides).

Cascade needs two source changes, kept in
`/mnt/data1/sage-plus-mass-offset-eval/cascade-search-mode.patch`:

1. Pass `search_mode` through `mods.json` into generated Sage configurations
   (`crates/cascade/src/mods.rs`). Without it the field is silently dropped.
2. Accept Sage Plus `0.1.0-beta.3` and run-summary schema 9
   (`crates/cascade/src/sage_plus.rs`). This is the pinned-version upgrade Cascade needs
   for beta.3 regardless of this feature; the pinned revision `f9d29b4` still expects
   beta.2 and schema 8.

## FDR calibration (entrapment)

Paired 1:1 target/entrapment peptide reference
(`benchmarks/results/fdrbench-validation/seeds/20260902/paired.fasta`, FDRBench-generated),
searched against PXD001468 HEK293 DDA with oxidation as the tested modification and
N-terminal acetyl indexed in both arms. Peptide-level, non-decoy, best representative per
sequence; every identified peptide was present in the pairing file.

| Mode | Nominal q | Targets | Entrapments | Combined FDP | Lower-bound FDP |
| --- | --- | --- | --- | --- | --- |
| database | 0.01 | 1087 | 6 | 0.0110 | 0.0055 |
| mass offset | 0.01 | 1112 | 8 | 0.0143 | 0.0071 |
| database | 0.05 | 1176 | 47 | 0.0769 | 0.0384 |
| mass offset | 0.05 | 1180 | 46 | 0.0750 | 0.0375 |

Offset search tracks the database search at both thresholds. The 1% difference is two
entrapment peptides out of roughly 1100, which is within counting noise at this scale, so
these runs show no degradation rather than establishing equivalence. The 5% behavior
(both arms above nominal under the combined convention) is a property of this reference
and dataset, not of offset search.

## Performance at realistic scale

PXD001468 HEK293 DDA, full reviewed human FASTA, one 241 MB file, 20k PSMs at 1%:

| Configuration | Indexed peptides | Index build | Search | Peak RSS | PSMs @1% | Peptides @1% |
| --- | --- | --- | --- | --- | --- | --- |
| Oxidation indexed | 7,633,197 | 10.6 s | 4.4 s | 3.62 GB | 20342 | 11244 |
| Oxidation as offset | 5,602,995 | 8.0 s | 9.2 s | 2.69 GB | 20385 | 11252 |
| Oxidation + phospho as offsets | 5,602,995 | 8.9 s | 15.7 s | 2.69 GB | 20448 | 11321 |

The trade is explicit: index size, build time, and memory follow the indexed
modifications only, while search time grows with the offsets and their tested placements
(roughly 2.1x search time for the first offset, 3.6x for two, on this data). Identifications
increase slightly because offsets are applied to every indexed peptidoform.

## Other search paths

- **Determinism**: repeating an offset search reproduces byte-identical results.
- **Chimeric search**: completes with offsets and reports rank-2 PSMs (2328 PSMs, 57
  rank-2, against 2307 / 54 for the database arm on HEK).
- **Wide-window search**: completes with offsets (3324 PSMs on HEK).

## Not established here

- The entrapment arm above accepts about 1100 peptides, so it detects gross
  miscalibration only. A larger entrapment reference, or several spectrum files, would be
  needed to resolve differences of a few tenths of a percent.
- Offsets have not been evaluated with TMT, SILAC channels (they are rejected with
  channel-aware labels by design), or DIA beyond the wide-window smoke test above.
- Labile modifications that leave the fragment entirely are out of scope; offsets assume
  the modification is retained, aside from configured neutral losses.
- Combinations of offsets, or two copies of one offset, on a single peptide are not
  searched.
