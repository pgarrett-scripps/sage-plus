# Generic fragment losses (Beta 14)

`database.fragment_losses` adds residue-dependent loss ions (b-H2O, y-NH3 and
so on) to full candidate scoring. This page records how the integration was
chosen and why the setting is optional rather than on by default.

## What is published practice and what is ours

| Item | Source |
| --- | --- |
| Considering -H2O and -NH3 forms of b and y ions during scoring | Published practice, handled differently by each engine (next table). |
| Water from S, T, E, D and ammonia from R, K, N, Q | Published practice: Mascot's fragmentation rules ("b-H2O if b significant and fragment includes STED", "b-NH3 ... includes RKNQ"). |
| Water and ammonia losses are more prominent in low-energy ion-trap CID than in beam-type HCD | Common expectation, not a published comparison we could find. Tabb et al. (Anal. Chem. 2003, 75:1155) quantify losses in ion-trap CID only (b-NH3 about half of b; y about 5 times y-NH3; ammonia loss enriched with N/Q, water loss only mildly with S/T/E). Tested below. |
| Losses never enter the preliminary fragment index | Our design choice: keeps first-pass speed and memory unchanged. |
| Losses as two separate rescoring features rather than hyperscore alternatives | Our design choice, selected by the benchmark below. |
| Parent-supported losses (a loss peak counts only when its intact fragment matched) | Our design here, tested below. Andromeda reportedly offers loss peaks only when the main b/y ion is present (Cox et al., J. Proteome Res. 2011, 10:1794; wording not verified, the full text was not accessible). |
| Default off | Our decision from the benchmark below. |

How published engines score water and ammonia losses:

| Engine | Default | Residue rule | How loss peaks enter the score | Source |
| --- | --- | --- | --- | --- |
| Comet | off (`use_NL_ions = 0`) | none | Extra theoretical peaks at weight 0.2 (main peaks 1.0), 1+ fragments only, independent of the parent peak | `CometPreprocess.cpp`, `CometSearch.cpp` (github.com/UWPR/Comet); `use_NL_ions` parameter docs |
| Mascot | by instrument type (ESI-TRAP, ESI-QUAD-TOF, ESI-FTICR and others; none for ETD-TRAP) | NH3: R, K, N, Q; H2O: S, T, E, D | A separate ion series, used only "if b significant" (or y); a series at the random level is dropped | matrixscience.com/help/search_field_help.html#INSTRUMENT, fragmentation_help.html |
| X!Tandem | no loss series | - | Hyperscore uses a/b/c/x/y/z only | thegpm.org/TANDEM api docs |
| MS-GF+ | learned per parameter file | none | Ion types (-H2O, -NH3 and stacks) kept if seen in at least 15% of training spectra, scored by a learned rank log-likelihood | `ScoringParameterGeneratorWithErrors.java`, `NewScoredSpectrum` (github.com/MSGFPlus/msgfplus); Kim & Pevzner, Nat. Commun. 2014 |
| Andromeda | on (`IncludeWater`, `IncludeAmmonia`, `DependentLosses`) | amine, amide and hydroxyl side chains (reported) | Reportedly 1+ only and only when the main b/y ion is present (not verified) | MaxQuant mqpar.xml defaults; Cox et al. 2011 |

So no engine in this list adds loss forms as full-weight alternatives that
fill a cleavage in a matched-ion count. Each one either down-weights them
(Comet), gates them on the parent series or ion (Mascot, reportedly
Andromeda), learns their weight (MS-GF+), or leaves them out (X!Tandem).

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
  loss on ion-trap CID. The cause is diagnosed in "Why B loses PSMs" below:
  loss-only cleavages add to the matched-ion count, and the best wrong
  candidate collects more of them than the correct peptide.
- **C is neutral**: PSM, peptide and protein counts change by -0.7% to +0.7%,
  decoy PSMs at 1% change by at most one, and entrapment FDP never rises by
  more than about 0.1 percentage point at PSM or peptide level (protein FDP is
  within the noise of a few entrapment proteins: +0.13 points on Human_01,
  -0.5 on CID, -0.25 on ETciD).
- **D is not better than C**: keeping loss-only peaks for the next chimeric
  rank changes nothing measurable. C is kept because loss peaks belong to the
  peptide that explained them.

## Why B loses PSMs

B was rerun with the experiment commit `8761c3d`, which reproduces the first
round's B and C row for row (duckdb `EXCEPT ALL` both ways is empty on the
shared columns of PXD011070).

**It is not an implementation bug.** For a lost PSM (PXD011070 scan 4802,
target FPGQLNADLRK replaced by the decoy TKVSTATGPAPTK) an independent Python
recomputation from the mzML, with the approved site rules, 0.5 Da and charges
1-2, gives the target 25 intact matches plus 2 loss-only cleavages and the
decoy 19 plus 10, exactly the 27 and 29 `matched_peaks` Sage reported. The
unit tests cover the site rules (b fragments hold residues 0..=i, y the rest),
the neutral-mass arithmetic (a loss is subtracted from the neutral fragment
mass, so it is correct at every fragment charge) and stacking.

**Cause: a loss-only cleavage adds a matched-ion count, and wrong candidates
collect more of them than the right one.** With the approved rules 81% of b
fragments and 99.9% of y fragments carry a loss form (every tryptic y ion
contains its C-terminal K or R), so B roughly doubles the masses tried at
every cleavage. A cleavage whose intact ion already matched gains nothing;
an empty cleavage gains a new chance of a random match. The correct peptide
has few empty cleavages; each wrong candidate has many, and the runner-up is
the best of many wrong candidates, so it collects the most random fills.
Hyperscore change for confident targets that kept their peptide, against the
best other candidate for the same spectrum:

| Dataset | Correct peptide | Runner-up | Median `delta_next` base -> B | Confident targets that lost rank 1 (to a decoy) |
| --- | ---: | ---: | --- | ---: |
| PXD028735 Human_01 (HCD) | +1.8 | +6.2 | 18.9 -> 13.6 | 2,165 of 65,663 (841) |
| PXD011070 (ion-trap CID) | +4.4 | +9.7 | 17.4 -> 12.0 | 717 of 7,129 (360) |
| PXD004447 (ETciD) | +0.8 | +4.8 | 19.4 -> 13.2 | 364 of 13,101 (196) |

At the rank-1 level decoys gain as much as or more than confident targets
(CID: +4.3 matched peaks for decoys, +1.6 for confident targets; HCD: +0.86
against +0.77). The hyperscore alone then separates far worse: target-decoy
competition on the hyperscore alone passes 1,225 CID PSMs at 1% instead of
3,111. The linear discriminant recovers most of that, but not the lost rank-1
assignments.

Each candidate cause, tested with a variant that changes only that factor:

| Candidate cause | Test | Verdict |
| --- | --- | --- |
| Loss-only matches raise the matched-ion count (log-factorial term) | **B-int**: as B, but a loss-only cleavage adds its intensity and no count | **The cause.** B-int is within 0.7% of base everywhere (B: -1% to -12%). |
| Loss peaks replace the intact peak's intensity | **B-intact**: as B, but the intact form wins whenever it matched | Not the cause: B-intact equals B. |
| Decoys gain as much as targets | per-PSM comparison above | True, and worse for the best wrong candidate; this is how the count term hurts. |
| Ion-trap CID with 0.5 Da tolerance | B across datasets | Amplifies (CID -12%), not required (HCD -1% to -4%, ETciD -5%). |
| `max_fragment_losses` | **B-max2** (2 instead of 1) | Makes it worse (CID -24%, HCD -2% to -8%). |
| Over-broad site rules | **B-narrow**: water D/E, ammonia N/Q, on b and y | Halves the loss but still -1% to -8%: any residue rule still covers 76% of b and 71% of y fragments of tryptic peptides. |
| Implementation bug | independent recomputation, unit tests | None found. |

## Principled fixes

Every fix keeps the approved configuration shape; variants are selected with
the experiment variables of `8761c3d`.

| Variant | What it does | Published or ours |
| --- | --- | --- |
| B | loss forms are alternatives of their cleavage (first round) | ours (no engine above does this at full weight) |
| B-int | as B, a loss-only cleavage adds intensity but no count | ours |
| B-intact, B-narrow, B-max2 | diagnostics above | ours; B-narrow's ammonia N/Q follows Tabb et al. 2003, water D/E is ours |
| C | two LDA features, `matched_loss_peaks` and `loss_intensity_pct` (shipped) | ours |
| Cs | as C, counting only loss peaks whose intact parent (same cleavage and charge) matched: `loss_supported_fragments` | ours |
| Cf | as Cs, the count as a fraction of `matched_peaks` | ours |
| P0 | parent-supported: a supported loss adds its intensity to the parent cleavage, no count; unsupported loss peaks add nothing | ours; the parent gate resembles Andromeda's reported rule and Mascot's series gate |
| P05, P1 | as P0, and a supported loss adds 0.5 or 1 to the hyperscore's count (`matched_peaks` unchanged) | ours |
| P0-narrow | P0 with the B-narrow rules | ours |
| P0Cs | P0 in the score plus the Cs features | ours |
| W02 | every matched loss peak adds 0.2 to the count and 0.2 of its intensity | published weight: Comet `use_NL_ions` (0.2) |

"Restrict losses to CID" is answered by the per-dataset columns: no variant
does better on the ion-trap CID file than on HCD.

PSM / peptide change against base at 1% q-value (percent):

| Variant | HEK SILAC | PXD028735 Human_01 | PXD001468 HEK293T | PXD011070 ion-trap CID | PXD004447 ETciD |
| --- | ---: | ---: | ---: | ---: | ---: |
| A | +0.0 / +0.0 | +0.0 / +0.0 | +0.0 / +0.0 | +0.0 / +0.0 | +0.0 / +0.0 |
| B | -1.7 / -2.9 | -4.1 / -4.0 | -0.9 / -0.7 | -11.6 / -11.8 | -5.4 / -4.9 |
| B-int | +0.7 / -1.0 | +0.4 / +0.5 | +0.1 / +0.1 | -0.2 / -0.6 | -0.1 / -0.0 |
| B-intact | -1.5 / -3.1 | -4.2 / -4.1 | -0.9 / -0.8 | -12.0 / -11.5 | -5.4 / -4.9 |
| B-narrow | +0.0 / +0.3 | -2.1 / -2.0 | -0.7 / -0.7 | -7.5 / -7.8 | -3.7 / -3.4 |
| B-max2 | -4.3 / -5.1 | -7.7 / -7.4 | -2.1 / -2.1 | -24.1 / -22.9 | -6.3 / -5.8 |
| C | +0.1 / -0.1 | +0.1 / +0.0 | +0.0 / +0.0 | -0.0 / -0.1 | +0.1 / +0.2 |
| Cs | +0.4 / +0.6 | +0.1 / +0.3 | -0.1 / +0.1 | -0.1 / -0.1 | +0.0 / +0.2 |
| Cf | +0.0 / +0.1 | +0.1 / +0.2 | -0.0 / +0.0 | -0.1 / +0.0 | -0.0 / +0.2 |
| P0 | +0.0 / -0.1 | +0.2 / +0.0 | +0.0 / +0.0 | +0.0 / -0.3 | -0.0 / +0.0 |
| P05 | +0.5 / +0.8 | -0.1 / -0.1 | +0.1 / +0.3 | -0.2 / -1.0 | -0.0 / +0.1 |
| P1 | +0.6 / +0.4 | -0.7 / -0.3 | +0.0 / +0.2 | -1.9 / -2.6 | -0.1 / +0.1 |
| P0-narrow | -0.1 / +0.1 | +0.2 / +0.1 | +0.0 / +0.0 | +0.1 / +0.0 | +0.0 / +0.0 |
| P0Cs | +0.7 / +0.3 | +0.1 / +0.3 | -0.0 / +0.1 | -0.1 / -0.4 | +0.0 / +0.2 |
| W02 | +1.2 / -0.4 | +0.6 / +0.5 | +0.1 / +0.2 | -0.1 / -0.4 | -0.1 / +0.0 |

Entrapment FDP (full table below) stays within 0.1 point of base at PSM and
peptide level for C, Cs, Cf, P05, P1 and W02; Cs, Cf, P05 and P1 are slightly
lower. B-int, P0, P0-narrow and P0Cs raise PSM FDP on CID by 0.10-0.15 point,
and B-narrow by 0.22. Protein FDP moves by a few entrapment proteins either
way. Counted instead at a fixed 1% entrapment FDP (targets ranked by the
discriminant, cut where the combined FDP estimate exceeds 1%):

| Variant | Human_01 PSMs / peptides | CID PSMs / peptides |
| --- | ---: | ---: |
| base | 64,470 / 36,821 | 6,945 / 5,642 |
| B | 60,708 / 34,863 | 6,027 / 5,058 |
| C | 64,636 / 36,798 | 6,933 / 5,654 |
| Cs | 64,782 / 36,865 | - |
| P1 | 64,505 / 37,008 | - |
| W02 | 64,546 / 36,910 | 7,006 / 5,765 |

Only B moves beyond 1%. HEK293T's entrapment FDP is above 1% from the top of
the list for every variant, so it has no fixed-FDP count.

Parent-supported counting does fix the rank-1 problem: under P1 the correct
peptide gains 3 times as much hyperscore as the runner-up (Human_01 +8.4
against +2.8) and the median `delta_next` widens from 18.4 to 25.7. But
this does not become identifications after the discriminant: at 0.5 Da about
half of matched intact fragments also find a random peak at their loss mass,
so the extra count mostly rescales the matched-ion count for targets and
decoys alike. P1 is -1.9% on CID and -0.7% on Human_01 at 1% q-value, with
lower FDP.

<details>
<summary>All variants, full numbers at 1% q-value</summary>

| Dataset | Variant | PSMs | Peptides | Proteins | Decoy PSMs | FDP PSM | FDP pep | FDP prot |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| HEK SILAC | base | 2,217 | 1,429 | 631 | 21 | - | - | - |
|  | A | 2,217 | 1,429 | 631 | 21 | - | - | - |
|  | B | 2,180 | 1,387 | 608 | 20 | - | - | - |
|  | B-int | 2,232 | 1,415 | 635 | 21 | - | - | - |
|  | B-intact | 2,184 | 1,384 | 609 | 20 | - | - | - |
|  | B-narrow | 2,217 | 1,433 | 630 | 21 | - | - | - |
|  | B-max2 | 2,122 | 1,356 | 600 | 20 | - | - | - |
|  | C | 2,219 | 1,427 | 632 | 21 | - | - | - |
|  | Cs | 2,225 | 1,438 | 634 | 21 | - | - | - |
|  | Cf | 2,217 | 1,430 | 633 | 21 | - | - | - |
|  | P0 | 2,217 | 1,428 | 634 | 21 | - | - | - |
|  | P05 | 2,228 | 1,441 | 638 | 21 | - | - | - |
|  | P1 | 2,231 | 1,435 | 637 | 21 | - | - | - |
|  | P0-narrow | 2,214 | 1,430 | 633 | 21 | - | - | - |
|  | P0Cs | 2,233 | 1,434 | 636 | 21 | - | - | - |
|  | W02 | 2,244 | 1,423 | 632 | 21 | - | - | - |
| PXD028735 Human_01 | base | 65,663 | 37,904 | 4,286 | 655 | 1.29% | 1.47% | 1.79% |
|  | A | 65,663 | 37,904 | 4,286 | 655 | 1.29% | 1.47% | 1.79% |
|  | B | 62,971 | 36,371 | 4,172 | 628 | 1.30% | 1.58% | 2.05% |
|  | B-int | 65,938 | 38,078 | 4,279 | 658 | 1.31% | 1.52% | 1.79% |
|  | B-intact | 62,934 | 36,349 | 4,178 | 628 | 1.28% | 1.59% | 2.05% |
|  | B-narrow | 64,288 | 37,132 | 4,241 | 641 | 1.38% | 1.63% | 1.81% |
|  | B-max2 | 60,612 | 35,115 | 4,117 | 605 | 1.31% | 1.41% | 2.01% |
|  | C | 65,748 | 37,913 | 4,314 | 656 | 1.29% | 1.40% | 1.92% |
|  | Cs | 65,738 | 38,000 | 4,325 | 656 | 1.25% | 1.38% | 1.91% |
|  | Cf | 65,721 | 37,985 | 4,319 | 656 | 1.21% | 1.41% | 1.91% |
|  | P0 | 65,768 | 37,915 | 4,294 | 656 | 1.30% | 1.44% | 1.79% |
|  | P05 | 65,627 | 37,865 | 4,295 | 655 | 1.23% | 1.29% | 1.72% |
|  | P1 | 65,214 | 37,775 | 4,319 | 651 | 1.14% | 1.28% | 1.78% |
|  | P0-narrow | 65,774 | 37,949 | 4,293 | 656 | 1.29% | 1.45% | 1.79% |
|  | P0Cs | 65,748 | 38,023 | 4,322 | 656 | 1.23% | 1.40% | 1.91% |
|  | W02 | 66,055 | 38,095 | 4,301 | 659 | 1.29% | 1.47% | 1.78% |
| PXD001468 HEK293T | base | 42,753 | 23,463 | 6,202 | 426 | 2.57% | 2.87% | 3.43% |
|  | A | 42,753 | 23,463 | 6,202 | 426 | 2.57% | 2.87% | 3.43% |
|  | B | 42,387 | 23,302 | 6,125 | 422 | 2.50% | 2.95% | 3.13% |
|  | B-int | 42,776 | 23,479 | 6,204 | 426 | 2.61% | 2.85% | 3.47% |
|  | B-intact | 42,366 | 23,269 | 6,112 | 422 | 2.49% | 2.94% | 3.09% |
|  | B-narrow | 42,448 | 23,301 | 6,118 | 423 | 2.53% | 2.91% | 3.18% |
|  | B-max2 | 41,845 | 22,967 | 6,069 | 417 | 2.58% | 2.95% | 3.16% |
|  | C | 42,758 | 23,464 | 6,202 | 426 | 2.56% | 2.87% | 3.43% |
|  | Cs | 42,727 | 23,488 | 6,199 | 426 | 2.51% | 2.87% | 3.48% |
|  | Cf | 42,738 | 23,471 | 6,200 | 426 | 2.52% | 2.87% | 3.48% |
|  | P0 | 42,765 | 23,467 | 6,202 | 426 | 2.59% | 2.89% | 3.47% |
|  | P05 | 42,815 | 23,540 | 6,207 | 427 | 2.57% | 2.94% | 3.42% |
|  | P1 | 42,756 | 23,510 | 6,217 | 426 | 2.53% | 2.97% | 3.56% |
|  | P0-narrow | 42,763 | 23,466 | 6,203 | 426 | 2.59% | 2.89% | 3.47% |
|  | P0Cs | 42,737 | 23,483 | 6,198 | 426 | 2.54% | 2.87% | 3.48% |
|  | W02 | 42,794 | 23,502 | 6,206 | 426 | 2.55% | 2.93% | 3.47% |
| PXD011070 ion-trap CID | base | 7,129 | 5,927 | 1,115 | 70 | 1.36% | 1.57% | 3.44% |
|  | A | 7,129 | 5,927 | 1,115 | 70 | 1.36% | 1.57% | 3.44% |
|  | B | 6,304 | 5,226 | 995 | 62 | 1.48% | 1.58% | 0.89% |
|  | B-int | 7,112 | 5,894 | 1,115 | 70 | 1.46% | 1.46% | 3.44% |
|  | B-intact | 6,277 | 5,246 | 1,008 | 61 | 1.43% | 1.57% | 1.46% |
|  | B-narrow | 6,592 | 5,467 | 993 | 64 | 1.58% | 1.64% | 1.19% |
|  | B-max2 | 5,414 | 4,572 | 906 | 53 | 1.46% | 1.72% | 1.63% |
|  | C | 7,127 | 5,920 | 1,111 | 70 | 1.36% | 1.51% | 2.92% |
|  | Cs | 7,124 | 5,923 | 1,109 | 70 | 1.31% | 1.51% | 2.66% |
|  | Cf | 7,123 | 5,927 | 1,110 | 70 | 1.31% | 1.51% | 2.93% |
|  | P0 | 7,130 | 5,910 | 1,107 | 70 | 1.51% | 1.64% | 2.93% |
|  | P05 | 7,114 | 5,866 | 1,100 | 70 | 1.26% | 1.34% | 2.42% |
|  | P1 | 6,993 | 5,773 | 1,098 | 68 | 1.23% | 1.24% | 3.23% |
|  | P0-narrow | 7,133 | 5,928 | 1,120 | 70 | 1.46% | 1.69% | 3.43% |
|  | P0Cs | 7,123 | 5,904 | 1,113 | 70 | 1.46% | 1.70% | 2.92% |
|  | W02 | 7,120 | 5,901 | 1,107 | 70 | 1.36% | 1.52% | 2.93% |
| PXD004447 ETciD | base | 13,101 | 11,750 | 2,212 | 130 | 0.82% | 0.85% | 1.87% |
|  | A | 13,101 | 11,750 | 2,212 | 130 | 0.82% | 0.85% | 1.87% |
|  | B | 12,397 | 11,180 | 2,086 | 122 | 0.78% | 0.87% | 1.27% |
|  | B-int | 13,083 | 11,747 | 2,205 | 129 | 0.79% | 0.82% | 1.87% |
|  | B-intact | 12,398 | 11,177 | 2,082 | 122 | 0.78% | 0.87% | 1.13% |
|  | B-narrow | 12,613 | 11,352 | 2,081 | 125 | 0.82% | 0.85% | 1.13% |
|  | B-max2 | 12,271 | 11,073 | 2,055 | 121 | 0.82% | 0.91% | 1.15% |
|  | C | 13,114 | 11,774 | 2,192 | 130 | 0.82% | 0.85% | 1.62% |
|  | Cs | 13,104 | 11,774 | 2,194 | 130 | 0.79% | 0.88% | 1.48% |
|  | Cf | 13,100 | 11,771 | 2,201 | 130 | 0.79% | 0.88% | 1.61% |
|  | P0 | 13,098 | 11,750 | 2,212 | 129 | 0.82% | 0.85% | 1.87% |
|  | P05 | 13,098 | 11,767 | 2,213 | 129 | 0.79% | 0.88% | 2.00% |
|  | P1 | 13,083 | 11,756 | 2,212 | 129 | 0.79% | 0.85% | 2.00% |
|  | P0-narrow | 13,103 | 11,750 | 2,213 | 130 | 0.82% | 0.88% | 2.00% |
|  | P0Cs | 13,106 | 11,769 | 2,192 | 130 | 0.79% | 0.85% | 1.48% |
|  | W02 | 13,083 | 11,751 | 2,211 | 129 | 0.79% | 0.85% | 2.14% |

</details>

Wall time and memory on a quiet machine (PXD028735 Human_01 / PXD011070):
base 27.3 s / 12.8 s, B 29.2 / 13.1, C 29.1 / 12.8, Cs 29.4 / 12.9, P0 29.1 /
13.0, P1 29.2 / 13.0, W02 28.9 / 13.0. Losses cost 6-8% wall time on HCD and
0-2% on CID; peak RSS is unchanged (3.13-3.17 GB and 2.34-2.36 GB).

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
  removal removes them. After the second round no in-score integration gains
  beyond about 1%: the count-free fixes (B-int, P0, P0-narrow, W02) remove
  B's loss but land in the same ±1% band as C, and none of them beats base by
  more than 1% at a fixed 1% entrapment FDP. Counting supported losses (P05,
  P1) costs up to 2.6% on CID. So enabling `fragment_losses` still means C.
  Cs (parent-supported loss count as a feature) equals C within noise, with
  slightly lower PSM FDP; it was not adopted.
- Default: **off** for every instrument type. Neither HCD, ion-trap CID nor
  ETciD shows a gain that would justify turning losses on by default.
- It remains available for users who want the loss evidence in their PSM
  tables or Percolator/mokapot features.

## Reproducing

```bash
python benchmarks/run_fragment_losses.py \
  --base <b14/per-mod-flr sage> --candidate <6aa89ff sage> --output runs \
  --variants base A B C
python benchmarks/run_fragment_losses.py ... \
  --datasets pxd028735-human01 pxd011070-itcid \
  --variants base-chimera C-chimera
# second round, candidate built from 8761c3d
python benchmarks/run_fragment_losses.py --base <base sage> --candidate <8761c3d sage> \
  --output exp --variants base A B B-int B-intact B-narrow B-max2 C Cs Cf \
  P0 P05 P1 P0-narrow P0Cs W02
```

D (first round, chimeric keep-peaks) was selectable only in `6aa89ff`. The
second-round variants need the experiment variables of `8761c3d`; the final
binary ignores them and runs C whenever `fragment_losses` is set.
