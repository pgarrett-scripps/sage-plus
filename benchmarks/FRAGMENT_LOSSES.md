# Generic fragment losses (Beta 14)

`database.fragment_losses` adds residue-dependent loss ions (b-H2O, y-NH3 and
so on) to full candidate scoring. This page records how the integration was
chosen and why the setting is optional rather than on by default.

## What is published practice and what is ours

| Item | Source |
| --- | --- |
| Considering -H2O and -NH3 forms of b and y ions during scoring | Published practice: Comet `use_NL_ions`, Mascot's `b-H2O`/`y-NH3` ion series, X!Tandem and others. |
| Water from S, T, E, D and ammonia from R, K, N, Q | Published practice: the residue rules of Mascot's fragmentation rules, used widely since. |
| Water and ammonia losses are more prominent in low-energy ion-trap CID than in beam-type HCD | Literature claim (for example the ion-trap fragmentation statistics of Tabb et al., Anal. Chem. 2003). Tested below. |
| Losses never enter the preliminary fragment index | Our design choice: keeps first-pass speed and memory unchanged. |
| Losses as two separate rescoring features rather than hyperscore alternatives | Our design choice, selected by the benchmark below. |
| Default off | Our decision from the benchmark below. |

## Variants

- **base**: the `b14/per-mod-flr` binary (ad1b86f), no `fragment_losses` key.
- **A**: the candidate binary without the key. Must be identical to base.
- **B**: loss forms are alternatives of their cleavage inside the hyperscore
  and `matched_peaks` (at most one form per cleavage and charge, as for
  modification `neutral_losses`).
- **C**: the hyperscore and `matched_peaks` ignore loss forms. Each PSM gets
  `matched_loss_peaks` (cleavage and charge pairs whose best loss form
  matched) and `loss_intensity_pct` as two extra linear-discriminant features.
  Chimeric peak removal removes the loss peaks of the preceding rank.
- **D**: as C, but chimeric peak removal keeps peaks matched only by loss
  forms. Only meaningful with `chimera: true`, so it is compared against base
  and C with `chimera: true, report_psms: 2`.

Every variant used the approved loss configuration:

```json
"fragment_losses": {
  "Water":   {"mass": 18.010565, "sites": ["S","T","E","D"], "ion_kinds": ["b","y"], "allow_modified": false},
  "Ammonia": {"mass": 17.026549, "sites": ["R","K","N","Q"], "ion_kinds": ["y"]}
},
"max_fragment_losses": 1
```

B, C and D were run from experiment commit `6aa89ff`, which selected the
variant with the environment variable `SAGE_PLUS_FRAGMENT_LOSS_SCORING`. The
final commit hard-codes C and removes the variable; its binary reproduces the
experiment's C results row for row (checked on PXD028735 Human_01, PXD011070
and PXD004447).

## Datasets

All searches: tryptic, 1 missed cleavage, 7-50 residues, Cys carbamidomethyl,
Met oxidation, 10 ppm precursor, isotope errors -1 to 3, generated reversed
decoys, `predict_rt` and `protein_grouping` on, `report_psms: 1`.

| Dataset | Instrument / activation | Fragment tolerance | FASTA |
| --- | --- | --- | --- |
| HEK SILAC K6R6 (repository sample) | Orbitrap HCD | 20 ppm | human reviewed |
| PXD028735 `LFQ_Orbitrap_DDA_Human_01` | Orbitrap HCD | 20 ppm | human + yeast + E. coli + iRT |
| PXD001468 HEK293T (01A and 01B) | Q Exactive HCD | 20 ppm | human + yeast + E. coli + iRT |
| PXD011070 `ch_23Aug2018_HeLa_Std_1` | ion-trap CID (ITMS) | 0.5 Da, no deisotoping | human + yeast + E. coli + iRT |
| PXD004447 `F1_..._JurTryp_ETciD30_2_2` | Orbitrap ETciD, b/c/y/z-dot | 20 ppm | human + yeast + E. coli + iRT |

The four human-only samples are searched against the multi-species FASTA, so
yeast and E. coli matches are false. FDP is the combined entrapment estimate
N_e(1 + 1/r)/(N_t + N_e) from `picked_fdr_entrapment.py`, at PSM, peptide
and protein level at 1% q-value.

## Results at 1% q-value

PSMs are target PSMs at `spectrum_q <= 0.01`; peptides and proteins are
distinct targets at `peptide_q`/`protein_q <= 0.01` (single-protein rows for
proteins). Decoys are decoy PSMs at `spectrum_q <= 0.01`.

| Dataset | Variant | PSMs | Peptides | Proteins | Decoy PSMs | FDP PSM | FDP pep | FDP prot | Wall (s) | Peak RSS (GB) |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| HEK SILAC | base | 2,217 | 1,429 | 631 | 21 | - | - | - | 7.8 | 1.70 |
| | A | 2,217 | 1,429 | 631 | 21 | - | - | - | 7.8 | 1.70 |
| | B | 2,180 | 1,387 | 608 | 20 | - | - | - | 8.1 | 1.70 |
| | C | 2,219 | 1,427 | 632 | 21 | - | - | - | 8.0 | 1.70 |
| PXD028735 Human_01 | base | 65,663 | 37,904 | 4,286 | 655 | 1.29% | 1.47% | 1.79% | 30.6 | 3.14 |
| | A | 65,663 | 37,904 | 4,286 | 655 | 1.29% | 1.47% | 1.79% | 30.7 | 3.14 |
| | B | 62,971 | 36,371 | 4,172 | 628 | 1.30% | 1.58% | 2.05% | 31.1 | 3.15 |
| | C | 65,748 | 37,913 | 4,314 | 656 | 1.29% | 1.40% | 1.92% | 31.1 | 3.14 |
| PXD001468 HEK293T | base | 42,753 | 23,463 | 6,202 | 426 | 2.57% | 2.87% | 3.43% | 19.8 | 2.65 |
| | A | 42,753 | 23,463 | 6,202 | 426 | 2.57% | 2.87% | 3.43% | 20.0 | 2.65 |
| | B | 42,387 | 23,302 | 6,125 | 422 | 2.50% | 2.95% | 3.13% | 21.1 | 2.65 |
| | C | 42,758 | 23,464 | 6,202 | 426 | 2.56% | 2.87% | 3.43% | 21.3 | 2.64 |
| PXD011070 ion-trap CID | base | 7,129 | 5,927 | 1,115 | 70 | 1.36% | 1.57% | 3.44% | 13.3 | 2.34 |
| | A | 7,129 | 5,927 | 1,115 | 70 | 1.36% | 1.57% | 3.44% | 14.7 | 2.34 |
| | B | 6,304 | 5,226 | 995 | 62 | 1.48% | 1.58% | 0.89% | 14.4 | 2.34 |
| | C | 7,127 | 5,920 | 1,111 | 70 | 1.36% | 1.51% | 2.92% | 14.4 | 2.35 |
| PXD004447 ETciD | base | 13,101 | 11,750 | 2,212 | 130 | 0.82% | 0.85% | 1.87% | 16.8 | 3.40 |
| | A | 13,101 | 11,750 | 2,212 | 130 | 0.82% | 0.85% | 1.87% | 17.9 | 3.40 |
| | B | 12,397 | 11,180 | 2,086 | 122 | 0.78% | 0.87% | 1.27% | 17.6 | 3.41 |
| | C | 13,114 | 11,774 | 2,192 | 130 | 0.82% | 0.85% | 1.62% | 17.4 | 3.41 |

Chimeric search (`chimera: true`, `report_psms: 2`):

| Dataset | Variant | PSMs | Peptides | Proteins | Decoy PSMs | FDP PSM | FDP pep | FDP prot | Wall (s) | Peak RSS (GB) |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| PXD028735 Human_01 | base | 68,779 | 39,267 | 4,276 | 686 | 1.34% | 1.53% | 1.59% | 50.8* | 3.25 |
| | C | 68,809 | 39,243 | 4,308 | 687 | 1.30% | 1.42% | 1.71% | 35.5 | 3.24 |
| | D | 68,821 | 39,240 | 4,300 | 687 | 1.32% | 1.42% | 1.72% | 36.0 | 3.25 |
| PXD011070 ion-trap CID | base | 7,226 | 5,971 | 1,103 | 71 | 1.24% | 1.26% | 3.21% | 17.4 | 2.35 |
| | C | 7,226 | 5,994 | 1,094 | 71 | 1.29% | 1.55% | 2.43% | 17.8 | 2.35 |
| | D | 7,226 | 5,981 | 1,099 | 71 | 1.24% | 1.38% | 3.22% | 20.3 | 2.36 |

Wall times were measured with other agents' searches running on the same
16-core machine and vary by about 10% between repeats; the starred base run
was slowed by that contention. Peak RSS is `/usr/bin/time` maximum resident
set size. Losses add no measurable memory and 0-8% wall time (C against
base), within the contention noise.

**Off is identical to base.** A matches base in every row of
`results.sage.parquet` (duckdb `EXCEPT ALL` both ways is empty) on all five
datasets; `results.json` differs only in output paths. The Parquet files'
byte hashes differ only through provenance metadata (binary identity).

## Reading the results

- **B loses identifications everywhere**: 1-12% fewer PSMs, with the largest
  loss on ion-trap CID. Loss forms let decoys fill cleavages they would
  otherwise miss. Our reading (not separately tested) is that this lifts the
  decoy score distribution at least as much as the target one, most of all
  at the wide 0.5 Da CID tolerance. It is rejected.
- **C is neutral**: PSM, peptide and protein counts change by -0.7% to +0.7%,
  decoy PSMs at 1% change by at most one, and entrapment FDP never rises by
  more than about 0.1 percentage point at PSM or peptide level (protein FDP is
  within the noise of a few entrapment proteins: +0.13 points on Human_01,
  -0.5 on CID, -0.25 on ETciD).
- **D is not better than C**: keeping loss-only peaks for the next chimeric
  rank changes nothing measurable. C is kept because loss peaks belong to the
  peptide that explained them.

## The ion-trap CID claim

Loss ions are more prominent in ion-trap CID, as the literature says:
confident CID target PSMs carry 6.9% of MS2 intensity in matched loss peaks,
against 2.5% on Orbitrap HCD (Human_01) and 0.3% on ETciD. But at 0.5 Da
fragment tolerance loss peaks also match by chance: decoy PSMs carry 5.3%
(7.6 matched loss peaks against 11.3 for confident targets, about the same
ratio to `matched_peaks` as targets). On HCD the loss matches track the
ordinary matches equally for targets and decoys (0.32 against 0.34 loss peaks
per matched peak). So loss evidence mostly restates what the b/y matches
already say, and does not separate targets from decoys well enough to add
identifications. We did not find the larger CID gain the literature leads one
to expect; the effect here is none on CID, HCD or ETciD.

## Decision

- Integration: **C**, loss ions as two separate rescoring features, never in
  the hyperscore, `matched_peaks` or the preliminary index; chimeric peak
  removal removes them.
- Default: **off**. There is no clear gain on any dataset, so the key stays
  out of the default configuration. It remains available for users who want
  the loss evidence in their PSM tables or Percolator/mokapot features.

## Reproducing

```bash
python benchmarks/run_fragment_losses.py \
  --base <b14/per-mod-flr sage> --candidate <6aa89ff sage> --output runs \
  --variants base A B C
python benchmarks/run_fragment_losses.py ... \
  --datasets pxd028735-human01 pxd011070-itcid \
  --variants base-chimera C-chimera D-chimera
```

With the final binary, B and D are no longer selectable; `C` is simply the
approved configuration.
