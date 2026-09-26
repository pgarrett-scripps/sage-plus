# DIA search for Sage Plus

Status: M1 (tier-1 pseudo-spectrum mode) is wired into the CLI as the opt-in
`"dia": {"mode": "pseudo"}` / `--dia pseudo` (crate `crates/sage-dia`; user docs in
`DOCS.md`, "DIA pseudo-spectrum search"). The tier-2 and hill-filter experiments live in
the spike `crates/sage-cli/examples/dia_explore.rs` (see section 9) and are not shipped.
M2 (timsTOF diaPASEF in pseudo mode, ion-mobility-aware grouping) is also shipped; see
section 10.

Goal: a fast, memory-light DIA mode. It sits between spectrum-centric and
peptide-centric search, in the spirit of DIA-Umpire and MSFragger-DIA. It uses
koth hills for MS1 and MS2 and does not aim for DIA-NN depth.

Structure:

1. **Beta 10 milestone:** DIA-Umpire-style pseudo-spectrum mode. It is opt-in and ships first.
2. How other tools do DIA.
3. The full hybrid pipeline (later milestones).
4. Memory and speed design.
5. Milestones and effort.
6. Validation datasets.
7. What changes where: koth, Sage core and `sage-dia`.
8. Spike results.

---

## 1. Beta 10: pseudo-spectrum mode (DIA-Umpire style)

### Design contract (as built in the spike)

1. **Anchor.** Every pseudo-spectrum is anchored on one koth MS1 isotope
   feature. The feature gives the monoisotopic m/z, the charge, the apex RT,
   and the mobility on timsTOF (M2; the Orbitrap spike has no IM).
2. **Grouping.** Fragment hills inside that feature's isolation window are
   grouped with it by elution overlap: apex proximity plus profile
   correlation. The grouping is loose, and hills may be shared between
   pseudo-spectra.
3. **Search.** The pseudo-spectra go through a normal **closed, DDA-style**
   search: ±10 ppm on the feature mass, isotope errors -1..+1, and the
   feature's charge. They are *not* searched with the isolation window as the
   precursor tolerance.
4. **Baseline.** The existing wide-window search on the same file is the
   baseline (section 8).
5. **Optional Q3 fallback.** DIA-Umpire Q3-style groups of fragment hills with
   no MS1 feature are searched wide-window. It is implemented as
   `pseudo::build_orphans` and `--q3`. It added **0 peptides at 1%** on the
   spike file, so it stays **off by default**.

### Idea

A DIA MS2 scan holds fragments from every precursor in a 10–25 m/z window. That
is what makes wide-window search slow (a ±window×z precursor tolerance) and
chimeric. Pseudo-spectrum mode undoes the multiplexing before search:

1. **MS1 precursors.** koth builds MS1 hills and groups them into charged
   isotope features. Each feature gives a monoisotopic m/z, a charge and a
   summed isotope elution profile.
2. **MS2 fragment hills.** koth builds MS2 hills separately for each isolation
   window. It already does this: `detect_ms2_hills_from_iter` and
   `group_ms2_hills_by_window`.
3. **Grouping.** For each precursor, and each window that contains its m/z:
   - take the fragment hills whose apex is within `apex_tolerance` cycles of
     the precursor apex;
   - keep those whose profile correlates with the precursor profile at
     `min_corr` or better, over at least `min_overlap` cycles;
   - drop fragments above precursor m/z × z.

   Grouping is loose on purpose. One fragment hill may appear in many
   pseudo-spectra, and the search and FDR sort out the rest.
4. **Pseudo-MS2 spectrum.** Each kept hill becomes one peak (hill m/z, apex
   intensity), and the spectrum carries a known precursor m/z and charge.
5. **Ordinary Sage search.** Narrow ppm precursor tolerance, isotope error
   -1..1, no chimera and no wide window. Then the existing LDA, q-values,
   protein inference and LFQ, unchanged.

Why it should be faster and lighter:

- **Fewer candidates.** A search at ±10 ppm with a known charge has about 1000×
  fewer candidates per spectrum than a ±12 m/z window × 3 charges.
- **Cleaner spectra.** Each spectrum has tens of fragments rather than a
  150-peak chimeric mix.
- **Raw data is dropped early.** Only compact hills are kept once koth is done.

### Code (already in `crates/sage-dia`)

- `hills::Channel`: compact, m/z-sorted hill store for one window or for MS1.
  It uses about 24 B per hill plus 4 B per profile point, and has an
  apex-sorted index for time lookups.
- `pseudo::PrecursorTrace::from_koth(&koth_core::Feature)` builds a precursor
  from a koth feature.
- `pseudo::build(precursor, ms1, window, window_index, &PseudoSettings) -> Option<PseudoSpectrum>`
  builds one pseudo-spectrum.
- The spike turns each `PseudoSpectrum` into a `RawSpectrum` (ms_level 2,
  centroid, one `Precursor` with charge) and passes it to the normal
  `SpectrumProcessor` and `Scorer`. No Sage core change is needed.

### koth-core API used (v0.11.0, via git)

koth-core is not on crates.io; only `koth-ms` 0.10.0 is published. It is pinned
by git tag:

```toml
koth-core = { git = "https://github.com/pgarrett-scripps/koth", tag = "v0.11.0" }
```

Calls used:

- `koth_core::hills::detect_hills_from_iter(iter, &cfg.hills, &cfg.file)` for MS1 hills.
- `koth_core::run_features(&hills, &cfg.features, &cfg.file)` for charged
  isotope features. It provides `Feature::elution_profile()`,
  `monoisotopic_mz()` and `charge`.
- `koth_core::hills::detect_ms2_hills_from_iter(iter, &cfg.ms2_hills(), &cfg.file)`
  for MS2 hills, with `[hills_ms2]` overrides through `HillsMs2Overrides`.
- `koth_core::group_ms2_hills_by_window(hills)` to split MS2 hills by window.
- `koth_core::Spectrum` / `Peak` / `IsolationWindow` as the input types.
  Sage `RawSpectrum` converts to them in about 20 lines (`to_koth` in the
  spike).

**No koth change is required for Beta 10.** Nice-to-haves, none blocking:

1. **Publish koth-core 0.11 to crates.io.** A git dependency is fine for a
   beta, but a crates.io release is cleaner for `cargo install`.
2. **Parallel, push-style or per-window MS2 hill output.**
   `detect_ms2_hills_from_iter` runs the windows serially: about 40 s on the
   spike file under load. Splitting the spectra by window on the Sage side and
   calling it once per window under rayon gives the same 2.16M hills in
   **0.8 s**. That workaround is in the spike, so it is not blocking. With `emit_ms2`,
   `run_pipeline_streaming_from_spectra` buffers every spectrum (plus a clone).
   Calling the MS1 and MS2 detectors separately, as the spike does, avoids this.
   A `Ms2HillDetector::push(spectrum)` with per-window emission would let Sage
   stream a file and never hold the raw MS2 data.
3. **A slim hill type or `into_compact()`.** `Hill` is about 200 B plus an
   `Arc` profile. Sage compacts it itself today.
4. **Streaming readers for Thermo and Bruker.** Only mzML streams in koth, but
   Sage has its own readers, so this does not matter to Sage.

### Config shape (opt-in)

```json
"dia": {
  "mode": "pseudo",           // "off" (default) | "pseudo" | later "hybrid"
  "apex_tolerance": 2,        // cycles between fragment and precursor apex
  "min_corr": 0.3,            // fragment vs precursor profile Pearson
  "min_overlap": 3,           // cycles
  "half_window": 6,           // cycles compared around the apex
  "min_peaks": 6,
  "max_peaks": 150,
  "ms2_min_scans": 3,         // koth [hills_ms2].min_scans
  "koth": { }                 // optional pass-through KothConfig overrides
}
```

When `dia.mode = "pseudo"`:

- the reader output for DIA files goes through `sage-dia`;
- `wide_window` and `chimera` are ignored, and a warning is logged;
- `precursor_tol` and the isotope errors apply as for DDA.

Output: normal parquet PSMs. The `scan` field is a synthetic pseudo-spectrum id
that records the source window and precursor feature, and `rt` is the precursor
apex RT.

### Wiring (Beta 10 scope)

1. Add `sage-dia` as a dependency of `sage-cli` behind a cargo feature `dia`, on
   by default. koth-core is pure Rust.
2. `Input` gains `dia: Option<DiaOptions>`.
3. In the per-file loop, if the mode is pseudo and the file has MS2 isolation
   windows, `sage_dia::pseudo_spectra(raw, &opts) -> Vec<RawSpectrum>` replaces
   the raw MS2 list before `SpectrumProcessor`.
4. MS1 spectra are still passed through for LFQ and MS1 features as today.
   Alternatively, precursor intensity comes from the koth feature, with the
   same effect.
5. Add tests: a synthetic two-precursor DIA file whose pseudo-spectra separate,
   and an e2e smoke run on a small mzML.

Effort (see section 5 for all milestones):

| Task | Effort |
|---|---|
| Pipeline, done in the spike | done |
| CLI and config wiring, docs, tests | ~2 days |
| Parameter tuning on 2 datasets (Orbitrap and diaPASEF) | ~2 days |
| Ion-mobility dimension for diaPASEF (koth hills carry `im`; match by IM too) | ~2 days |
| **Total** | **~1–1.5 weeks** |

---

## 2. How other tools handle DIA

**DIA-Umpire** (Tsou 2015)
- Detects MS1 and MS2 features, i.e. hills: isotope clusters with elution peaks.
- Builds pseudo-MS/MS by pairing each precursor with fragments whose apex is
  close and whose elution profile correlates (Pearson). Fragments can be shared.
- Splits pseudo-spectra into quality tiers (Q1 from MS1 features, Q2 and Q3
  from unfragmented precursors in MS2). Each tier is searched with a DDA engine
  and FDR is controlled per tier.
- This is exactly the Beta 10 design, minus the tiers.

**MSFragger-DIA** (2023–24)
- Spectrum-centric: each DIA MS2 scan is searched directly against the fragment
  index with a wide precursor window, and several peptides are reported per
  scan.
- The candidates then go to peptide-centric rescoring: fragment XIC
  co-elution, MS1 support and RT/IM prediction. DIA-NN or a library quantifies
  them.
- Sage's `wide_window` + `chimera` is already the first half of this.

**DIA-NN** (used here for scoring ideas only)
- Peptide-centric: for every library precursor it extracts fragment XICs and
  scores a large feature set. The features include:
  - co-elution correlation with the best fragment;
  - apex shape;
  - MS1 correlation;
  - predicted spectrum similarity;
  - RT and IM deviation.
- An ensemble of neural networks is trained per run.
- Interference removal picks the most co-eluting fragments for quant.
- Its memory scales with the library × the XIC extraction.

**Spectronaut directDIA**
- A library-free first pass (pseudo-spectra or spectrum-centric search),
  followed by a peptide-centric targeted pass with the resulting library and
  mProphet-style semi-supervised scoring.

**Sage `wide_window`** (upstream and in Sage Plus)
- Precursor tolerance = isolation window × charge, with `chimera` subtraction.
- Fast and simple, but no chromatographic evidence is used, so it caps depth
  and FDR discrimination.

Takeaway:
- The pseudo-spectrum route (DIA-Umpire) needs no new scoring and reuses the
  whole Sage stack. That makes it the right first ship.
- Chromatographic co-elution scoring (MSFragger-DIA, DIA-NN) is what adds
  depth. It comes second, as a rescoring layer on candidates from either mode.

---

## 3. Full hybrid pipeline (after Beta 10)

- **(a) Candidate generation (spectrum-centric).** Candidates come from either
  or both of:
  - a wide-window chimeric search (`wide_window`, `chimera`, `report_psms` 5–10);
  - a pseudo-spectrum search (section 1).

  Taking the union improves recall.
- **(b) koth hills.**
  - MS1 hills, then isotope features.
  - MS2 hills per isolation window, stored as a compact `Channel`: an m/z-sorted
    array plus an apex index plus an f32 profile arena.
  - The raw spectra are then dropped.
- **(c) Peptide-centric features for each candidate** (`coelution::score`,
  implemented):
  - Look up theoretical b/y fragments (z ≤ 2) in the window's hills at ±15 ppm.
  - Features:
    - fragments with a hill;
    - fragments co-eluting with the top-k fragment sum (leave-one-out
      Pearson ≥ 0.7);
    - fragments with a co-apex;
    - mean pairwise correlation of the top-k fragments;
    - apex spread and apex offset;
    - share of signal at the apex;
    - MS1 hill present, MS1 correlation, and MS1-to-fragment apex ΔRT.
  - Later additions: predicted RT and IM deltas (Sage models) and predicted
    fragment intensities.
  - No dense XIC matrix is built: each candidate reads about 20 short hill slices.
- **(d) Semi-supervised model.**
  - Percolator-style iterative LDA: decoys are the negatives, targets at 1% are
    the positives, over 3 rounds. It reuses
    `LinearDiscriminantAnalysis::train_regularized`.
  - Cross-validation folds are needed before shipping.
  - An SVM or small GBDT is an option if LDA plateaus.
  - Keep one best candidate per precursor elution group, so co-eluting
    chimeras compete.
- **(e) FDR and quant.**
  - Existing `spectrum_q_value` / peptide / protein picked-FDR.
  - Quant from the precursor MS1 feature area (koth), or from the summed top-N
    co-eluting fragment hill areas. Fragments shared with a better-scoring
    precursor are excluded, which is a simple form of interference removal.
  - Match-between-runs through the existing LFQ RT alignment.

---

## 4. Memory and speed design

Principles:

- Never build dense XICs.
- Keep hills, not spectra.
- Store everything as f32 in one arena per window.
- Parallelise over windows (koth) and precursors or candidates (rayon).

Per-file estimate for a 2 h Orbitrap DIA (about 3k cycles × 75 windows, 230k
MS2, about 25M MS2 centroids):

| Stage | Memory | Time (16 cores) |
|---|---|---|
| Read raw (Sage reader) | ~1.5–2 GB transient (the dominant cost today) | ~10 s |
| koth MS1 + features | < 200 MB | seconds |
| koth MS2 hills (2.16M hills, windows in parallel) | koth `Hill` structs ~0.5 GB transient | 0.8 s |
| Compact channels (MS1 + 151 windows) | 145 MB (~24 B/hill + 4 B/point) | — |
| Pseudo-spectrum build (57k spectra) | ~20 MB | 0.1 s |
| Search | pseudo: 57k spectra × ~50 peaks | 0.4 s (wide-window: 40–60 s) |

The remaining large cost is holding the raw file. Two fixes:

1. Stream spectra into koth detectors as they are read (needs koth
   nice-to-have 2).
2. Drop each MS2 spectrum once its window detector has consumed it.

Both would bring peak RSS to about the index size plus the hills.

---

## 5. Milestones

| # | Milestone | Effort | Ships in |
|---|---|---|---|
| M0 | Core fix: `score_psms` uses `train_regularized` (constant feature columns make the solve fail; see section 8) | ~0.5 day | **Beta 10** |
| M1 | Pseudo-spectrum mode, opt-in (section 1): wiring, tests, tuning | 1–1.5 weeks | **Beta 10** |
| M2 | diaPASEF: IM-aware hill grouping (koth hills carry `im`) plus a timsTOF benchmark | ~1 week | Beta 10 or 11 |
| M3 | Co-elution rescoring of candidates from wide-window or pseudo (section 3 c–d), with CV folds | ~2 weeks | Beta 11 |
| M4 | Fragment-hill quant with interference exclusion, plus MBR via LFQ alignment | ~2 weeks | later |
| M5 | Streaming ingest (koth push API) to cut peak RSS to hills plus index | ~1 week (needs koth) | later |
| M6 | Optional predicted fragment intensities / spectral library scoring | 2–3 weeks | later |

---

## 6. Validation datasets (PRIDE API confirmed)

- **PXD028735** (Van Puyvelde et al. 2022, LFQ benchmark, human/yeast/E. coli)
  - Orbitrap DIA raws `LFQ_Orbitrap_AIF_*`: 0.9–2.0 GB each; single-species and
    condition A/B mixtures.
  - The spike file `LFQ_Orbitrap_AIF_Ecoli_01.raw` is windowed DIA: 75 MS2
    windows per MS1 cycle.
  - Also timsTOF diaPASEF (`.d.zip`, 12–22 GB each; over the 10 GB spike budget)
    and SCIEX SWATH.
  - The ground-truth ratios make it the primary quant benchmark.
- **PXD017703** (Meier et al. 2020, diaPASEF, HeLa and two-proteome): the IM
  benchmark for M2. Zips are 11–340 GB, so pick single runs.
- **PXD002952** (LFQbench, TripleTOF 5600 SWATH): a classic ratio benchmark
  for another vendor.

---

## 7. What changes where

**koth**
- Nothing is required.
- Proposals (none blocking):
  1. Publish koth-core 0.11.
  2. A push/per-window MS2 hill API (streaming).
  3. An optional compact hill output.

**Sage core** (`crates/sage`): nothing for M1, because pseudo-spectra are plain
`RawSpectrum`s. For M3, a small hook is needed, like the glyco/crosslink
branches:
- an optional extra-feature slot on `Feature` (e.g. `dia: Option<Box<[f32]>>`)
  that LDA can include;
- or `sage-dia` runs its own LDA after `score_psms`.

**sage-cli**
- The `dia` config block.
- The per-file branch that calls `sage_dia::pseudo_spectra`.
- Feature-gated.

**`crates/sage-dia`** (new) holds all DIA logic:
- `hills` (compact channels);
- `pseudo` (grouping);
- `coelution` (rescoring features);
- later quant and IM.

---

## 8. Spike results

Machine and file:
- File: PXD028735 `LFQ_Orbitrap_AIF_Ecoli_01.raw` (0.92 GB), searched against
  E. coli UP000000625 (4403 reviewed proteins, 422,800 peptides; tryptic,
  1 missed cleavage, Cys carbamidomethyl, Met oxidation).
- Acquisition: 3106 MS1 and 232,903 MS2 spectra, 151 isolation windows of
  8 m/z (75 per cycle, staggered), 2.9 s cycle.
- Machine: 16 cores shared with other agents' builds; load average 15–65
  during the runs, so wall times are noisy. The headline pair below was
  measured back to back at load about 15.

### Wide-window vs pseudo-spectrum (the Beta 10 comparison)

| | Wide-window (sage CLI, `wide_window` + `chimera`, report_psms 5, RT model, LDA) | Pseudo-spectrum (spike, min_corr 0.5, apex ±2 cycles, closed ±10 ppm) |
|---|---|---|
| Spectra searched | 232,903 | 57,773 |
| Mean peaks per spectrum (after processing) | 71.0 | 51.9 |
| PSMs at 1% | 82,629 | 11,037 |
| Precursors at 1% | not reported by CLI (spike, hyperscore only: 4,324) | 6,003 |
| Peptides at 1% | **6,975** | **5,425** (78%) |
| Wall time, whole file | 29.8 s | **6.3 s** (search stage 0.4 s) |
| Peak RSS | 3.0 GB | 2.35 GB (1.7 GB of it is the raw read) |

Reading the table:
- Pseudo mode is about **5× faster end to end**. It searches 4× fewer and
  cleaner spectra, and its search stage is about 100× faster.
- It is 0.65 GB lighter. Most of the remaining peak is Sage's reader holding
  the raw file.
- It reaches **78% of the wide-window peptides**.

Why the gap:
- **Ceiling.** Only 32,140 charged MS1 features anchor spectra, which caps the
  peptide count.
- **Missing models.** The spike has no RT model or mass-alignment features.
  The CLI integration would add them.
- **PSM counts are not comparable.** The wide-window PSM count is inflated by
  up to 5 PSMs per chimeric scan, while pseudo has one per spectrum. Compare
  peptides.

### Pseudo-spectrum parameter sweep (regularized LDA, peptides at 1%)

| min_corr | apex tol | MS2 min_scans | spectra | peaks/spectrum | peptides |
|---|---|---|---|---|---|
| 0.3 | 2 | 3 | 59,648 | 61.7 | 5,494 |
| 0.5 | 2 | 3 | 57,773 | 51.9 | 5,425 |
| 0.6 | 2 | 3 | 56,053 | 46.2 | 5,356 (Sage LDA: 5,450) |
| 0.7 | 2 | 3 | 53,295 | 40.3 | 5,152 |
| 0.8 | 2 | 3 | 48,565 | 33.9 | 4,960 |
| 0.5 | 3 | 2 | 59,670 | 59.0 | 5,434 |

What the sweep shows:
- Results are flat between 0.3 and 0.6. Loose grouping is fine, as
  DIA-Umpire found.
- The default should be about 0.5 with apex ±2.
- **Q3 fallback** (min_corr 0.6): 49,179 groups without an MS1 feature,
  searched wide-window, gave 0 peptides at 1%. Keep it off by default.

### Finding: Sage's `score_psms` LDA is unregularized and fails

On pseudo-spectra, and on this spike's wide-window run (which has no RT model),
several LDA columns are constant: rank, delta_best, ims and the model deltas.
`Gauss::solve` then fails, and scoring silently drops to hyperscore. Across the
sweep, Sage LDA gave about 3.1k peptides on the failing runs and about 5.4k on
the runs where it fitted.

The fix is small and belongs in core: use `train_regularized` (already present)
in `score_psms`. The spike works around it with a regularized semi-supervised
LDA on Sage-like features.

### Rescoring spike (for M3): fragment-hill co-elution on wide-window candidates

Wide-window candidates (report_psms 5, 513,677 candidates) were scored with the
features from `coelution::score`, which took 1.4 s for all candidates.

Medians of the features:

| feature | targets at q ≤ 1% | all targets | decoys |
|---|---|---|---|
| ln n_coeluting | 2.57 | 1.39 | 0.69 |
| frac fragments with a hill | 0.53 | 0.08 | 0.06 |
| top-k fragment correlation | 0.98 | 0.56 | 0.36 |
| ln apex spread | 0.00 | 0.51 | 0.69 |
| MS1 hill present | 1 | 0 | 0 |
| MS1 correlation | 0.98 | 0 | 0 |

Peptides at 1%, starting from the spike's hyperscore baseline:

| model | peptides at 1% |
|---|---|
| Baseline | 3,777 |
| Hill-only LDA | 4,006 |
| Combined | **4,253** (+13%) |

The separation is strong. Against the CLI's full-LDA baseline (6,975), the
gain still has to be measured once the features are wired into core LDA (M3).

### Reproduce

```text
CARGO_TARGET_DIR=/mnt/data1/build-cache/feat-dia cargo build --profile fast-release -p sage-cli --example dia_explore
dia_explore --raw LFQ_Orbitrap_AIF_Ecoli_01.raw --fasta ecoli.fasta --out DIR \
    --mode wide|pseudo|rescore|tiered|hillfilter [--min-corr 0.5 --apex-tolerance 2 --ms2-min-scans 3 --q3]
```

Outputs from the runs are in `/mnt/data1/explore-data/dia/runs/`.

## 9. M1 two-tier benchmark (Orbitrap AIF E. coli, PXD028735)

Peptides at 1% FDR. "CLI" rows are `sage` runs; "example" rows are `dia_explore`, which uses
the same search but its own FDR plumbing.

| Variant | Path | Peptides | Runtime | Peak RAM |
|---|---|---|---|---|
| (a) raw wide-window, chimeric | CLI | 6,975 | 30 s | 3.0 GB |
| (b) tier 1 (`dia.mode = "pseudo"`) | CLI | 5,567 | 6 s | 2.3 GB |
| (c) tier 1 + tier 2 (corr 0.3), separate q-values | example | 5,546 (tier 2 +121, about +2%) | 9 s | 2.1 GB |
| (c) tier 1 + tier 2, pooled, tier as an LDA feature | example | 4,749 | 9 s | 2.1 GB |
| (d) wide-window on hill-filtered scans | CLI | 5,299 | 21 s | 2.8 GB |

- **Tier 2.** Hills matched by tier-1 IDs at 1% q are subtracted, with tolerance for b1/b2/y1/y2. This
  removed 93k of 2.16M hills. The remaining hills are grouped per window by co-elution
  (at most 100 peaks) and searched wide-window.
- **Pooled FDR.** Pooling tier 2 with tier 1 behind a tier feature lets tier-1 targets hide
  tier-2 decoys, and the result has fewer peptides than tier 1 alone. Separate q-values are
  safe but only add about 2%.
- **Hill filter.** Keeping only peaks on hills (14.6M of 23.5M) loses about 24% of the
  wide-window peptides.
- **Decision.** (c) does not beat (a), so M1 ships tier 1 only. Wide-window stays the default
  for DIA. Tier 2 would need a per-spectrum wide/closed switch in the scorer and per-tier
  FDR. That is about 1–2 h of agent time, for about 2% more peptides.
- **DDA unchanged.** With `dia` off, a DDA timsTOF search gives a byte-identical
  `results.sage.parquet` and `lfq.parquet` against origin/main.

## 10. M2 timsTOF diaPASEF benchmark (E. coli, PXD070049)

Run: `LFQ_Ultra2_diaPASEF_15min_50ng_Ecoli_01.d` (50 ng, 15 min gradient, 1,673 MS1 frames,
8 window groups × 3 boxes = 24 m/z × 1/K0 boxes of 25 m/z). Precursor ±15 ppm (pseudo),
fragment ±20 ppm, same FASTA and search settings as section 9.

| Variant | Path | Peptides | Runtime | Peak RAM |
|---|---|---|---|---|
| (a) raw wide-window, chimeric | CLI | 8,087 | 12 min 50 s | 16.6 GB |
| (b) tier 1 (`dia.mode = "pseudo"`), defaults before section 11 | CLI | 6,403 | 3 min 7 s | 6.1 GB |
| (b) tier 1, section 11 defaults (min_corr 0.3, im_tolerance 0.01) | CLI | 7,412 | 1 min 36 s – 2 min | 6.0 GB |
| (c) tier 1 + tier 2 (corr 0.3), separate q-values | example | 5,847 (tier 1 alone 5,773; tier 2 +74, about +1.3%) | 3 min 4 s | 5.8 GB |
| (c) tier 1 + tier 2, pooled, tier as an LDA feature | example | 5,355 | 3 min 4 s | 5.8 GB |
| (c) section 11 defaults, separate q-values (tier 2 scored with Sage LDA) | example | 7,307 (tier 1 alone 7,065; tier 2 +242) | 1 min 40 s | 6.0 GB |
| (d) timsrust spectra, peaks kept only on MS2 hills (see below) | example | 6,874 (regularized LDA 6,608) | 7 min 58 s | 17.2 GB |

- **Centroiding.** dnoise-core 0.5.0 watershed (`watershed::watershed_centroid`, the mode
  dnoise and koth's validated Bruker path use), after dnoise's vertical ion-mobility filter
  (MS/MS knobs for MS2) and horizontal-halo filter, all with dnoise defaults (tuned on this
  dataset). MS2 frames are processed per box. Without the two filters, watershed alone
  gave 3.16M MS1 features, 1.19M pseudo-spectra, 6,475 peptides, 19 min and 11.6 GB. With
  them: 1.60M MS1 hills, 138k features, 5.5M MS2 hills, 94k pseudo-spectra, 6,403
  peptides. Hill detection takes about 165 s (MS1 100 s, MS2 65 s); the search takes 5 s.
- **Ion mobility.** koth links hills by m/z and 1/K0 and gives features a 1/K0 apex. A
  feature is paired only with boxes that contain its m/z and its 1/K0, and fragment hills
  must be within `dia.im_tolerance` (default 0.01 1/K0 since section 11; was 0.03) of it.
- **Memory.** MS1 frames are streamed in chunks of 64, and MS2 frames one window group at
  a time. The wide-window baseline holds every per-scan-split raw MS2 spectrum
  (524k spectra, 1.39 billion peaks), which is where its 16.6 GB comes from.
- **Decision.** Same as Orbitrap: (b) finds 79% of (a)'s peptides in a quarter of the
  time and at 37% of the memory. Tier 2 adds about 1%. Wide-window stays the default.
- **What (a) is on timsTOF.** timsrust 0.6 does not read diaPASEF as wide windows. Its
  `SpectrumReader` switches to `timsrust_centroid`'s narrow reader, which ignores
  `bruker_config.ms2`. That reader finds MS1 isotope pairs and emits one spectrum per
  (precursor, MS2 frame), with the MS2 peaks in that precursor's scan range: 524k spectra.
  The precursor m/z is the detected MS1 m/z. The isolation window is +-6.25 m/z (a quarter of
  the 25 m/z box) around it, and the RT is the MS1 frame's. So (a) on timsTOF is
  precursor-anchored already, just without elution-profile grouping. That is why its gap to
  (b) is smaller than on Orbitrap.
- **(d) debug.** The first (d) run kept 0 of 1.39 billion peaks. Its matcher looked for a
  box with the spectrum's isolation bounds at the spectrum's RT, but both are
  precursor-based, as above. The fix: the MS2 frame is in the high 32 bits of the spectrum
  index. The box is the one containing the precursor's m/z and 1/K0, and the cycle is that
  frame's (`sage_dia::tims::frame_rts`). 12% of spectra (63k) match no box, because their
  precursor lies outside every box, and they are dropped. The filter keeps 48% of peaks.
  (d) finds 6,874 peptides. That is below both (a) and (b), at (a)'s memory. Not pursued.
- **DDA unchanged.** Checked on the denoised ddaPASEF twin of this run
  (`LFQ_Ultra2_PASEF_15min_50ng_Ecoli_01.d`; the raw twin was not on disk). Results were
  byte-identical to origin/main, with 8,872 peptides.

## 11. Why pseudo mode misses peptides wide-window finds

Missed = peptides at 1% peptide FDR in (a), any of its 5 chimeric ranks, but not in (b).
For each one, `dia_explore --mode diagnose` takes its best (a) PSM's m/z, charge and RT and
finds the MS1 isotope feature (same charge within the precursor tolerance and RT span; else
an isotope offset of ±1 or ±2; else another charge). It then counts the theoretical b/y ions
(charge 1–2) in that feature's pseudo-spectrum under the shipped settings and under
variants, and also counts ions shifted by +11 Th as a random-match null. Each peptide goes
to the furthest stage it reached. Scripts: `runs/diag/{targets,bucket,peaks}.py`.

Missed peptides, defaults before this section (min_corr 0.5, im_tolerance 0.03):

| Bucket | Orbitrap (1,864 missed) | timsTOF (2,384 missed) |
|---|---|---|
| No MS1 feature at that m/z / RT (feature detection) | 674 (36%): 391 with no MS1 hill at all | 172 (7%) |
| Wrong charge, or monoisotopic off by 2 or more isotopes | 126 (7%) | 214 (9%) |
| Feature found, fewer than 4 true fragments in its pseudo-spectrum (grouping) | 651 (35%) | 1,220 (51%) |
| Built with ≥ 4 fragments, lost at scoring or FDR | 413 (22%) | 778 (33%) |

**Grouping, timsTOF.** 75% of the missed peptides' pseudo-spectra hit the 150-peak cap. The
0.03 1/K0 window let co-eluting fragments of other precursors in, and ranking by intensity
then dropped real fragments. The diagnose variants, counting feature-found targets with
at least 4 true fragments (random-match null in brackets):

| Variant | timsTOF | Orbitrap |
|---|---|---|
| shipped (corr 0.5, IM 0.03, cap 150) | 802 (7) | 420 (0) |
| min_corr 0 | 898 (9) | 649 (0) |
| im_tolerance 0.015 | 898 (5) | — |
| no IM gate | 666 (11) | — |
| cap 300 | 1,046 (42) | 424 (0) |
| rank by correlation instead of intensity | 463 (7) | 419 (0) |
| apex tolerance 4 | 796 (8) | 464 (0) |

**Fix and sweep (CLI, peptides at 1%).** timsTOF: im_tolerance 0.015 gives 7,133, 0.01 gives
7,271, 0.007 gives 7,212, and 0.01 with min_corr 0.3 gives 7,412. min_corr 0 gives 7,343, and
cap 300 at 0.01 gives 7,277, so the cap is no longer binding. Orbitrap: min_corr 0.3 gives
5,763, 0 gives 5,807 and −1 gives 5,833. Apex tolerance 3 or 4 on top of min_corr 0 lowers it
(5,779, 5,765). The defaults are now min_corr 0.3 and im_tolerance 0.01. 0.3 is within 1% of
the best Orbitrap value and is the best on timsTOF.

| | (a) wide-window | (b) before | (b) now | Share of (a) |
|---|---|---|---|---|
| Orbitrap | 6,975, 30 s, 3.0 GB | 5,567 | 5,763, 6 s, 2.3 GB | 80% → 83% |
| timsTOF | 8,087, 12 min 50 s, 16.6 GB | 6,403 | 7,412, 2 min, 6.0 GB | 79% → 92% |

**What is left** (same diagnosis, new defaults):

| Bucket | Orbitrap (1,715 missed) | timsTOF (1,892 missed; (b) also finds 1,217 that (a) does not) |
|---|---|---|
| No MS1 feature | 674 (39%) | 170 (9%) |
| Wrong charge / mono | 126 (7%) | 211 (11%) |
| Too few fragments | 497 (29%) | 1,045 (55%) |
| Lost at scoring or FDR | 418 (24%) | 466 (25%) |

On Orbitrap the largest remaining bucket is feature detection. Of those misses, 391 have no
MS1 hill: the precursor is below MS1 detection, so no MS1-anchored method can recover them.
Another 283 have a hill but no charged isotope feature. On timsTOF it is still grouping,
but no single threshold variant recovers it. Where "all relaxed" recovers it, that is at
thousands of peaks with a null that matches just as often, which is noise, not signal.

**Example vs CLI tier 1.** The example's tier 1 (5,773 before, 7,065 now) is below the CLI
(6,403, 7,412) on the same pseudo-spectra, because the example rescores with its own
regularized LDA and peptide-level FDR. Turning off the CLI's RT and ion-mobility model
features alone drops it from 6,403 to 6,136. The rest comes from the example's LDA and
FDR, which do not use Sage's mass alignment or picked-peptide FDR. So tier 2 gains are
only comparable with the example's own tier 1. With the new defaults, tier 2 scored with
Sage's LDA adds 242 peptides (+3.4%) to the example's tier 1, and adds none when scored
with the regularized LDA.

## 12. Single-hill precursors (Orbitrap)

Section 11 left 283 Orbitrap misses that have an MS1 hill but no charged isotope feature.
`dia.hill_precursors` searches every MS1 hill that no charged koth feature claimed as a
precursor on its own: m/z is the hill's m/z, the elution profile is the hill's profile, and
the charge is guessed. There are 162,773 MS1 hills; 70,889 are unclaimed.

Sweep on `LFQ_Orbitrap_AIF_Ecoli_01` (baseline: 5,763 peptides at 1% FDR). "Gained" and
"lost" are peptides relative to the baseline; the last column counts recovered misses from
the section 11 "hill, no feature" bucket:

| Charges | Min scans | Min apex | Precursors added | Peptides | Gained / lost | Hill bucket recovered |
|---|---|---|---|---|---|---|
| 2, 3 | - | - | 141,778 | 5,866 | 375 / 277 | 95 |
| 2, 3, 4 | - | - | 212,667 | 5,849 | 377 / 296 | 91 |
| 2 | - | - | 70,889 | 5,927 | 340 / 179 | 93 |
| 2, 3 | 5 | - | 92,388 | 5,909 | 350 / 208 | 94 |
| 2 | 5 | - | 46,194 | 5,937 | 321 / 150 | 89 |
| 2 | - | median | 35,445 | 5,945 | 340 / 161 | - |
| 2 | - | 75th pct | 17,723 | 5,954 | 311 / 123 | - |
| 2, 3 | - | median | 70,890 | 5,895 | 366 / 237 | - |
| **2** | **5** | **median** | **24,907** | **5,959** | 326 / 133 | - |
| 2, 3 | 5 | median | 49,814 | 5,928 | 353 / 191 | - |
| 2 | 5 | 3e5 (absolute) | 14,274 | 5,963 | 300 / 103 | 84 |
| 2, 3 | - | 1e6 (absolute) | 4,896 | 5,922 | 185 / 26 | 41 |

Every setting recovers roughly a third of the hill bucket, plus about 100 of the wrong-charge
misses. Every extra precursor also adds decoy competition, which costs 100 to 300 peptides
that the baseline found. Charge 3 and weak or short hills cost more than they recover.
Results plateau around 5,950. The default is charge 2, at least 5 scans, and apex at or above
the median of the unclaimed hills. That gives **5,959 peptides (+196, +3.4%)**, in the same
6 s and 2.3 GB. The median threshold is relative, so it does not depend on the instrument's
intensity scale; an absolute 3e5 is no better.

| | (a) wide-window | (b) pseudo | Share of (a) |
|---|---|---|---|
| Orbitrap | 6,975 | 5,959 (was 5,763) | 85% (was 83%) |
| timsTOF | 8,087 | 7,412 (unchanged) | 92% |

timsTOF ignores the setting: there are only 170 no-feature misses there, and the hills carry
1/K0, which would need its own tuning. Its default-config run still gives 7,412.
