# Immonium ions (Beta 14)

Immonium ions are single-residue internal fragments, `H2N=CH-R+`, at the residue mass minus CO
(27.9949) plus a proton. This note records what the literature supports, what Sage Plus
implements, which parts are our own design, and the benchmark behind the default.

## Literature: what is published

| Claim | Source | Status |
|---|---|---|
| Immonium intensity tracks the presence of F, Y, W, P, H, V and L/I in a peptide; stronger when the residue is near the N-terminus | Hohmann et al. 2008, *Anal. Chem.* 80:5596 (QTOF CID) | Published |
| Phosphotyrosine immonium ion at m/z 216.04 is a specific marker of pY (precursor-ion scanning) | Steen et al. 2001, *Anal. Chem.* 73:1440 | Published |
| In HCD, the pY immonium ion pinpoints pY peptides; about 90% of pY peptides show it at normal HCD energy (after pY enrichment), detection depends on abundance | Olsen et al. 2007, *Nat. Methods* (HCD introduction) | Published |
| Acetyl-lysine: the 126.091 ion (143.118 minus NH3) is 98.1% specific; the 143.118 ion is not specific | Trelle & Jensen 2008, *Anal. Chem.* 80:3422 | Published |
| Modification-specific immonium ions can be mined by spectral binning | Kelstrup et al. 2014, SPIID, *MCP*; Geiszler et al. 2023, PTM-Shepherd, *Nat. Commun.* | Published |
| Search engines use configurable per-modification diagnostic peaks at PSM level; site reports give presence of the diagnostic peak | MaxQuant/Andromeda (Cox et al. 2011); MSFragger-Labile (diagnostic ions) | Published |
| Mascot has an immonium ion series among its fragmentation rules | Mascot fragmentation rules | Published; how it enters the score was not verified |
| Comet has no immonium ion series (a/b/c/x/y/z/z+1 and neutral losses only) | Comet parameters | Published |
| A rescoring feature: summed diagnostic plus immonium annotation intensity over total intensity (`diagnostic_ion_ratio`) | MS2Rescore `ms2` feature generator | Published (software) |
| Percolator | no standard immonium feature | n/a |

No publication we found uses immonium ions to choose between phosphosites (pY vs pS/pT on the
same peptide) inside a localization score. The pY ion says a pY is present in the isolated
peptides, not where.

## What Sage Plus implements

Opt-in `immonium` (default off). Per PSM, on the processed spectrum the search scored:

1. Residue ions of unmodified P, V, L/I, H, F, Y, W: counts of `explained` (ion observed, residue
   in the peptide), `missing` (residue in the peptide, ion absent) and `unexplained` (ion
   observed, residue absent).
2. Modified-residue ions from a list of `{name, residue, modification, mz}` (mass plus label, no
   formulas); defaults pY 216.0420 and acK 126.0913: `modified_explained` / `modified_unexplained`.
3. Output to `immonium.tsv` and, with `write_pin`, five `.pin` columns for external rescoring.
4. `rescore: true` (default false) appends the five counts to the LDA feature row.
5. Never in the fragment index, hyperscore or localization.

Our own design, not published:

- The three-way explained / missing / unexplained split and using the counts as LDA features.
  MS2Rescore's published feature is an intensity ratio, not counts.
- `missing` is only counted when the ion is at or above the lowest retained peak, so a scan range
  starting above the ion is not read as absence.
- A peak counts only if singly charged (or of unknown charge) after deisotoping.
- A residue carrying any modification does not count as that unmodified residue; a modified ion
  matches when the modification mass is within 0.01 Da.
- The offline pY-ion test against PXD000138 known sites below. It is an evaluation only; nothing
  in localization reads immonium ions.

## Benchmark

Base is `b14/per-mod-flr` (ad1b86f); candidate is `b14/immonium`, release builds. Variants: `off`
(no `immonium` key), `on` (`immonium: true`), `rescore` (`{"rescore": true}`). Two alternating
repetitions each, one search at a time under the shared memory gate, on a machine shared with
other agents, so wall times are noisy. Counts are `run-summary.json` 1% FDR values. Runs,
configs and `analyze.py` are under `/mnt/data1/sage-plus-scientific/immonium-b14/`.

Datasets:

1. HEK SILAC K6R6, one file (`data/silac-k6r6/`, named static mods, no LFQ or library).
2. PXD028735 `LFQ_Orbitrap_DDA_Human_01`, searched against `hye-irt-defined.fasta`
   (human + yeast + E. coli); a yeast- or E. coli-only peptide is an entrapment hit.
3. PXD007058 phospho (two files, variable pS/pT/pY, prefilter, localization on).
4. PXD000138 synthetic phosphopeptide library (three files), known sites.

| Dataset | Variant | PSMs | Peptides | Proteins | Groups | Wall (s) | Peak RSS (GiB) |
|---|---|---:|---:|---:|---:|---|---:|
| SILAC | base | 3,273 | 1,530 | 662 | 711 | 10.8, 10.7 | 2.15 |
| SILAC | on | 3,273 | 1,530 | 662 | 711 | 10.7, 13.0 | 2.15 |
| SILAC | rescore | 3,285 | 1,526 | 670 | 716 | 10.6, 13.3 | 2.16 |
| HYE Human_01 | base | 65,739 | 37,912 | 4,273 | 4,370 | 28.3, 40.8 | 3.13 |
| HYE Human_01 | on | 65,739 | 37,912 | 4,273 | 4,370 | 35.1, 34.7 | 3.14 |
| HYE Human_01 | rescore | 65,982 | 38,091 | 4,278 | 4,379 | 33.5, 36.3 | 3.14 |
| PXD007058 | base | 17,292 | 7,211 | 2,449 | 2,487 | 53.6, 47.0 | 0.74 |
| PXD007058 | on | 17,292 | 7,211 | 2,449 | 2,487 | 60.6, 47.1 | 0.69 |
| PXD007058 | rescore | 17,374 | 7,257 | 2,456 | 2,493 | 47.4, 47.4 | 0.68 |
| PXD000138 | base | 10,838 | 3,745 | - | - | 20.5, 19.1 | 1.70 |
| PXD000138 | on | 10,838 | 3,745 | - | - | 18.8, 17.4 | 1.71 |
| PXD000138 | rescore | 10,980 | 3,793 | - | - | 19.8, 18.8 | 1.71 |

`off` matched base exactly on every count. With `off` and with `on`, `results.sage.parquet` and
`results.sage.ptm-sites.parquet` are identical to base (pandas frame equality) on all four
datasets: the option adds `immonium.tsv` and `.pin` columns and nothing else. Runtime and peak
RSS differences are within the run-to-run noise of this machine.

**Entrapment (HYE Human_01, rank-1 target peptides at peptide q <= 0.01).** Entrapment to human
sequence ratio r = 0.376; combined FDP = N_entrapment (1 + 1/r) / N.

| Variant | Peptides | Entrapment | Lower-bound FDP | Combined FDP |
|---|---:|---:|---:|---:|
| base / off / on | 35,310 | 161 | 0.46% | 1.67% |
| rescore | 35,468 | 168 | 0.47% | 1.73% |

The 158 extra peptides from `rescore` include 7 entrapment hits, which projects to about 16% of
them being false: the gain is not clean at this size, although 7 hits is within counting noise.

**Known sites (PXD000138, localization q <= 0.01, single-seed peptides).** Base, off and on are
identical (3,640 single-site localizations, 35 false, true FLR 0.96%, estimated 0.96%). With
`rescore`, 3,670 pass with 36 false (0.98%, estimated 0.98%).

**Feature separation (on, rank 1, mean per PSM).**

| Dataset | Set | Explained | Missing | Unexplained | Modified explained | Modified unexplained |
|---|---|---:|---:|---:|---:|---:|
| SILAC | targets at 1% | 0.59 | 0.04 | 0.90 | 0.000 | 0.002 |
| SILAC | decoys | 0.23 | 0.07 | 0.58 | 0.000 | 0.000 |
| HYE | targets at 1% | 1.04 | 0.06 | 2.07 | 0.000 | 0.003 |
| HYE | decoys | 0.77 | 0.11 | 2.56 | 0.000 | 0.003 |
| PXD007058 | targets at 1% | 0.64 | 0.14 | 0.70 | 0.000 | 0.005 |
| PXD007058 | decoys | 0.48 | 0.61 | 0.73 | 0.000 | 0.003 |
| PXD000138 | targets at 1% | 0.98 | 0.06 | 0.56 | 0.203 | 0.061 |
| PXD000138 | decoys | 0.28 | 0.40 | 1.24 | 0.043 | 0.215 |

The counts separate targets from decoys, most in the pY-rich synthetic library, but most of that
information is already in the fragment features: rescoring adds 0.4-1.3% PSMs.

**Offline pY test (PXD000138, our evaluation, not used by the engine).** Single-site phospho
localizations at q <= 0.01 on peptides containing both Y and S/T, joined to `immonium.tsv`:

| True residue | pY ion seen | pY ion absent |
|---|---:|---:|
| pY (1,000) | 750 (75%) | 250 |
| pS/pT (455) | 51 (11%) | 404 |

The ion is informative (likelihood ratio about 6.7 when present, 0.28 when absent), in line with
Olsen et al. 2007. But only 20 of 1,455 localizations are wrong across residue classes: 16 true
pS/pT called pY (14 without the ion) and 4 true pY called pS/pT (none with the ion). A rule that
demoted pY calls lacking the ion would fix at most 14 and put 250 correct pY calls at risk. The
headroom is too small, and no published method uses the ion for site choice, so localization is
unchanged.

## Recommendation

- Keep `immonium` off by default. When on, it is an output (`immonium.tsv`, `.pin` columns) for
  QC and external rescoring and changes nothing else.
- Keep `rescore` off by default. The 0.4-1.3% PSM gain is small, and on HYE the added peptides
  are enriched for entrapment hits. Revisit with more entrapment data before turning it on.
- Do not use immonium ions in localization.

## Limits

- Mass-offset PSMs are evaluated against the unmodified peptide.
- Peaks are matched at observed m/z, without `mass_recalibration` corrections; low-mass ions sit
  well inside a 10-20 ppm window.
- Only peaks kept by `max_peaks` (top 150 by default) after deisotoping are seen.
- The acetyl-lysine ion was not benchmarked: none of these searches allowed acetyl-K.
- The built-in `diagnostic_ions` QC ion `phospho_Y_immonium` is 216.0426, 0.6 mDa (2.8 ppm)
  above Y immonium + HPO3 (216.0420); it is inside that QC's 20 ppm window and left unchanged.
