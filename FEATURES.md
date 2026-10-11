# Sage Plus features

What Sage Plus adds to [Sage](https://github.com/lazear/sage), grouped by what it is for. Settings
and outputs are documented in [DOCS.md](DOCS.md).

## Highlights

**Search very large databases in bounded memory.** With `prefilter` on, proteins stream through an
index of the spectra, and only peptides that could match a spectrum enter the search index. A human
reference plus a 100x gut-microbiome catalog, which ran out of memory before Beta 11, now finishes
in 19 GiB (14,459 PSMs in 451 s). See [settings](DOCS.md#fasta) and the
[Beta 11 measurements](CHANGELOG.md#v010-beta11---2026-09-26).

**Localize PTM sites with a false-localization rate.** `ptm_localization` rescores each arrangement
of a PSM's modifications on site-determining ions against impossible-site decoys, with a separate
false-localization rate per modification type (7,231 phospho localizations at 1% on PXD007058
with the Beta 15 defaults; 303 when all types were pooled before Beta 14). It writes site
probabilities and localization q-values to PSM-level and protein-level site tables, a site-level
target-decoy q-value (`site_q_value`), and a reusable site library. See [PTM site localization](DOCS.md#ptm-site-localization).

**Search modifications without growing the index.** Mass-offset modifications are placed at
scoring time instead of expanded into the fragment index: a phosphorylation search used 0.16 GB
instead of 0.45 GB with the same PSMs ([evaluation](benchmarks/MASS_OFFSET.md)). With neutral
losses, this also covers crosslinker monolinks. See [details](DOCS.md#mass-offset-modifications).

**Quantify more precursors across runs.** Retention times are aligned across runs with a nonlinear
warp by default. On five PXD028735 LFQ runs this raised precursors at 1% q-value from 24,038 to
31,418 and tightened the human log2 ratio spread. See [quantification](DOCS.md#quantification).

**Read instrument files directly.** Thermo RAW and mzMLb are read without conversion (upstream Sage
reads mzML, MGF and Bruker `.d`), and timsTOF
m/z and 1/K0 use the instrument's full calibration models. Each spectrum carries its analyzer and activation (HCD, CID,
ETD, EThcD, ETciD), so recalibration never mixes scan types. See [inputs](DOCS.md#spectrum-paths).

## Features by outcome

### Find more, and more trustworthy, identifications

Better defaults and scoring for the spectra you already have.
- Initiator methionine clipping, on by default: 3.4% more PSMs at 1% FDR on HEK SILAC with protein N-terminal acetylation searched ([details](DOCS.md#initiator-methionine-clipping)).
- `z_dot` fragment ions for ETD and EThcD: 51% more ETD-scan PSMs than upstream Sage with `z` ions on PXD018176 ([head-to-head](benchmarks/HEADTOHEAD-2026-10.md)) ([details](DOCS.md#fragment-settings)).
- Search-time mass recalibration per file and per analyzer, kept only when it improves held-out error ([details](DOCS.md#other-settings)).
- Ambiguous residues: J scored as I/L, opt-in B/Z/X expansion with a `substitutions` column, I/L twins merged ([details](DOCS.md#ambiguous-residues)).
- Averagine-scored deisotoping, charge-aware fragment matching, and `ambiguity_sequence` marking unsupported regions ([details](DOCS.md#sequence-ambiguity-annotation)).
- Textbook picked target-decoy FDR for peptides, proteins and protein groups (upstream Sage also uses picked FDR; this version drops a target beaten by its own decoy, and changes IDs by at most 0.3%) ([details](DOCS.md#protein-inference)).
- Opt-in immonium-ion evidence and water/ammonia fragment losses as extra rescoring features, reported per PSM ([details](DOCS.md#immonium-ions)).
- Opt-in DIA pseudo-spectrum mode for Orbitrap DIA and diaPASEF ([details](DOCS.md#dia-pseudo-spectrum-search)).

### Search PTMs with less setup

Say exactly where a modification can go, and check it before searching.
- Named modifications with explicit sites such as `first_residue:K` or `peptide_n_term` ([details](DOCS.md#explicit-site-vocabulary)).
- Motif sites such as `motif:N*-{P}-[ST]`, checked against the source protein ([details](DOCS.md#motif-sites)).
- Per-modification `max_count` and `max_total_count`, variant caps, and a separate PEFF budget ([details](DOCS.md#static-and-variable-behavior)).
- Per-site `neutral_losses` and `immonium_ions`, e.g. the H3PO4 loss on pS/pT but not pY ([details](DOCS.md#static-and-variable-behavior)).
- PTM site libraries restrict placements to known sites instead of every eligible residue ([details](DOCS.md#ptm-site-libraries)).
- `--preview-modifications` shows eligible sites and variants for a peptide without a search ([details](DOCS.md#preview-modification-placement)).

### Run faster and use less memory

Most useful for large databases, many modifications, or many files.
- Streamed prefilter (above); `prefilter_min_matched_peaks` (default 4) trades a few weak PSMs for memory ([details](DOCS.md#fasta)).
- Compact peptide, fragment-index, and spectrum storage: 18% less peak memory than upstream Sage `v0.15.0-beta.2` on standard DDA and 41% less on a phospho search, at about the same speed on DDA ([head-to-head](benchmarks/HEADTOHEAD-2026-10.md)).
- Interleaved fragment-bucket search: 19-27% faster search phase with identical PSMs ([changelog](CHANGELOG.md)).
- `max_memory_gb` stops a run at a measured memory limit, and `--estimate` previews database size ([details](DOCS.md#memory-guard)).

### Quantify

Label-free and labeled MS1 quantification, and TMT reporter ions with unobserved channels written as null, not 0.
- LFQ with nonlinear alignment, configurable match-between-runs tolerance, and per-file MS2 confirmation ([details](DOCS.md#label-free-quantification-output)).
- A per-row LFQ `extraction_q_value` for MS2-backed and transferred rows, a diagnostic of wrong-peak quantification, not yet a calibrated transfer FDR ([details](DOCS.md#label-free-quantification-output)).
- Peak geometry per LFQ row (`apex_rt`, bounds, `fwhm`, `id_apex_offset`) and opt-in apex recentering ([details](DOCS.md#label-free-quantification-output)).
- SILAC, dimethyl, and custom label channels defined on modifications, with channel-aware LFQ ([details](DOCS.md#modification-channels)).
- Optional timsTOF MS1 denoising with dnoise for LFQ (`bruker_config.denoise`, [details](DOCS.md#spectrum-paths)).

### Get outputs you can inspect and reproduce

Typed, versioned files instead of logs you have to parse.
- Parquet is the canonical output, with versioned schemas in [`schemas/`](schemas/) and one-based protein coordinates ([details](DOCS.md#interpreting-sage-output)).
- Every Parquet file records its provenance (version, commit, configuration, FASTA hash) in the footer ([details](DOCS.md#output-provenance)).
- `run-summary.json` records outputs, warnings, peak memory, fitted models, and recommended tolerances; it and the TSV outputs have published schemas.
- QC on every search: `digestion.tsv` and a polymer contamination check, plus opt-in `diagnostic_ions.tsv` ([details](DOCS.md#quality-control-outputs)).
- Empirical spectral-library export in Parquet and PSI mzSpecLib ([details](DOCS.md#empirical-spectral-libraries)).

### Automate runs safely

For pipelines and agents that launch searches without a person watching.
- JSON Schema for configurations (`sage --write-config-schema`) and `--validate-only` checks ([details](DOCS.md#machine-readable-jobs)).
- `--threads` (or the `threads` key) sets the worker count ([details](DOCS.md#performance-and-complexity)).
- JSONL progress events (`--events-jsonl`), a Rust runner API with cancellation, and no output replaced without `--overwrite`.
- Pre-digested peptide lists and custom cleavage sites extend the database ([details](DOCS.md#fasta)).

