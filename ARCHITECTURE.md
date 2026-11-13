# Architecture

Sage Plus is an independently maintained downstream distribution of
[Sage](https://github.com/lazear/sage). This file is the map: what each crate owns, the order
the pipeline runs in and which module owns each stage. `FEATURES.md` lists what Sage Plus adds,
`UPSTREAM.md` explains the relationship to upstream, `CLAUDE.md` covers building and testing,
and `DOCS.md` documents every configuration key.

## The four crates

One Cargo workspace, four crates. Dependencies point downward only: nothing depends on
`sage-cli`, and `sage-core` depends on neither the readers nor the CLI.

| Crate | Directory | Owns |
|---|---|---|
| `sage-core` | `crates/sage` | The search itself: digestion, the fragment index, scoring, FDR, quantification. No I/O, no CLI, no file formats. |
| `sage-cloudpath` | `crates/sage-cloudpath` | Every byte in and out: spectrum file formats, local and S3/GCS paths, Parquet and TSV writers, JSON schemas. |
| `sage-cli` | `crates/sage-cli` | The `sage` binary and the embeddable API: configuration parsing, run orchestration, output artifacts, `run-summary.json`. |
| `sage-dia` | `crates/sage-dia` | DIA pseudo-spectrum construction (`dia.mode = "pseudo"`): hills, co-elution, pseudo-MS2 spectra. |

`sage-core`'s crate name is `sage-core` but its directory is `crates/sage`, and the binary
`sage` is built from `sage-cli`. Both surprise newcomers.

On the `feat/crosslink` branch a fifth crate, `sage-xlink` (`crates/sage-xlink`), holds the
MS-cleavable crosslink search. It depends on `sage-core` only and is linked into `sage-cli`
only by the `crosslink` Cargo feature, which is off by default; see
`docs/explore/CROSSLINK_SEARCH.md`.

## The pipeline, in order

`Runner::run` (`sage-cli/src/runner/execution.rs`) drives all of it. One chunk of
`batch_size` input files is read and searched at a time; everything from FDR onward runs once
over the merged results.

| Stage | Owner | Notes |
|---|---|---|
| Configuration | `sage-cli/src/input.rs`, `config_schema.rs`, `config_presets.rs` | Parses and validates the JSON configuration, applies defaults, generates `schemas/config.schema.json`. |
| Database build | `sage-core/src/database/` | `builder.rs` digests the FASTA (`digest.rs`, `../enzyme.rs`, `../cleavage.rs`), applies modifications (`../peptide/enumerate.rs`), and builds the fragment index in `fragment_index.rs`. `estimate.rs` backs `--estimate`, `limit.rs` the peptide limits, `mass_offset.rs` open-search offsets. |
| Peptides and modifications | `sage-core/src/peptide/`, `modification.rs`, `modification/serde_impls.rs`, `site_map.rs`, `motif.rs`, `ptm_library.rs` | `peptide/enumerate.rs` generates variable-modification variants, `peptide/lookup.rs` holds the compact per-peptide modification storage. The hand-written configuration deserializers live in `modification/serde_impls.rs`. |
| Spectrum read | `sage-cloudpath/src/` | `mzml.rs`, `mzmlb.rs`, `mgf.rs`, `tdf.rs`, `thermoraw.rs`, dispatched by `util.rs`. Peaks are then processed (deisotoping, top-N) by `sage-core/src/spectrum.rs`. |
| Search and scoring | `sage-core/src/scoring/`, `database/query.rs` | `scoring/matching.rs` plans the fragment windows and counts preliminary matches, `database/query.rs` reads the index, `scoring/standard.rs` is the ordinary search, `scoring/chimera.rs` the chimeric one, `scoring/feature.rs` holds the `Feature` written to the output. `../mass_recalibration.rs` runs the optional two-pass recalibration. |
| Chunk orchestration | `sage-cli/src/runner/search.rs`, `prefilter.rs` | Reads a chunk, searches it, keeps the MS1 spectra LFQ needs. `prefilter.rs` is the opt-in streamed first pass for very large databases. |
| FDR and rescoring | `sage-core/src/fdr.rs`, `ml/` | `ml/linear_discriminant.rs` rescores, `ml/qvalue.rs` assigns q-values, `ml/retention_model.rs` and `ml/mobility_model.rs` fit the prediction models, `ml/retention_alignment.rs` aligns observed retention times across files, `protein_grouping.rs` does IDPicker parsimony. Driven by `sage-cli/src/runner/postprocess.rs`. |
| PTM localization | `sage-core/src/ptm.rs`, `ambiguity.rs` | Site-level assignment and the site FDR. |
| Quantification | `sage-core/src/lfq.rs`, `lfq/grid.rs`, `lfq/signal.rs`, `tmt.rs` | `lfq.rs` builds the precursor feature map and traces each precursor, `lfq/grid.rs` is the per-precursor retention-time x file x isotope grid, `lfq/signal.rs` the peak scoring. `tmt.rs` is isobaric quantification. |
| Write | `sage-cloudpath/src/parquet/`, `sage-cli/src/runner/artifacts.rs` | `parquet/results.rs`, `parquet/lfq.rs`, `parquet/sites.rs`, `parquet/spectral_library.rs` build the files; `artifacts.rs` decides what is written, including the TSV outputs and the HTML report. `runner/summary.rs` writes `run-summary.json`, `runner/qc.rs` the QC summaries. |

Two smaller paths join the same run: `sage-dia` turns a DIA file into pseudo-MS2 spectra
before the search when `dia.mode = "pseudo"`, and `sage-core/src/spectral_library.rs` builds the
empirical library entries after FDR, which `sage-cloudpath/src/parquet/spectral_library.rs` writes.

## Where the reader crates plug in

The Bruker timsTOF and Thermo readers are separate repositories, consumed here as published
crates.io dependencies:

- `sage-plus-tdf` (Bruker `.d`/TDF, miniTDF, TSF) is used by `sage-cloudpath/src/tdf.rs`,
  `tdf_spectra.rs`, `bruker_formats.rs` and `tims_mobility.rs`, and by `sage-dia/src/tims.rs`.
  `sage-cloudpath/src/denoise.rs` holds the optional MS1 denoising.
- `sage-plus-raw` (Thermo `.raw`) is used by `sage-cloudpath/src/thermoraw.rs`.

Both are converted into `sage_core::spectrum::RawSpectrum` at that boundary, so nothing above
`sage-cloudpath` knows which format a spectrum came from. A reader change reaches Sage Plus only
through a crates.io release; to try an unreleased reader locally, use a `[patch.crates-io]`
section or a path override and never commit it.

## Where tests live

- A unit test file lives in `crates/<crate>/tests/unit/<module>.rs` and is attached to its
  module by three lines at the bottom of the source file:

  ```rust
  #[cfg(test)]
  #[path = "../tests/unit/ambiguity.rs"]
  mod test;
  ```

  Nested modules mirror the source tree (`tests/unit/runner/search.rs` for
  `src/runner/search.rs`). No test implementations in production source files.
- Black-box integration tests go directly in `crates/<crate>/tests/`:
  `crates/sage/tests/integration.rs`, `crates/sage-cli/tests/integration.rs` and
  `crates/sage-cli/tests/feature_combinations.rs`.
- Spectrum and FASTA fixtures live in `crates/sage-cloudpath/tests/data` and `tests/` of each crate.
- Output formats are pinned by the committed schemas in `schemas/`; `schemas/config.schema.json`
  is generated from `sage-cli/src/input.rs` and a unit test fails if the committed copy drifts.

## Conventions that are easy to miss

- Only the settings that define the search space are required. Every other option must have a
  default that reproduces the previous behaviour.
- Public API, CLI flags, configuration keys and output column names are a compatibility
  surface; renaming one is a deliberate decision.
- Use the field's standard term, one term per concept: intensity (not height), modifications or
  PTMs (not marks).
- Keep diffs against `upstream/master` small and explainable. `git log upstream/master..main` is
  the real feature list.
