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

## Not established here

- FDR calibration under offset search was not measured with entrapment. Decoy placement
  is symmetric by construction (decoys receive offsets on the same rules, and library
  placements are mirrored onto reversed sequences), and rank-1 decoy counts are
  comparable between modes (806 vs 799 on HEK), but an entrapment arm in
  `benchmarks/run_scientific.py` remains the way to establish calibration.
- Labile modifications that leave the fragment entirely are out of scope; offsets assume
  the modification is retained, aside from configured neutral losses.
- Combinations of offsets, or two copies of one offset, on a single peptide are not
  searched.
