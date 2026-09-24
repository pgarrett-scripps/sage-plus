# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

Sage Plus uses an independent release sequence beginning with `v0.1.0-beta.1`. Earlier Sage
entries are retained below for provenance.

## [Unreleased]

### Added
- `mass_recalibration` (`"off"` by default; `"static"`, `"linear"`, or `"auto"`) corrects
  precursor and fragment m/z during the search. A sampled discovery search of each file (up to
  25,000 MS2 spectra, stratified by acquisition group so interleaved scan cycles do not alias
  with the sampling stride) selects rank-1 target PSMs at 1% Poisson q-value computed within
  each acquisition group, fits a per-file precursor model on the union of those PSMs and a
  fragment model per group on its own PSMs, each on a deterministic fit split, and accepts a static offset, a linear RT and/or m/z term, or a
  smooth additive curve (`auto` only, at most five intervals) only when it improves held-out
  residuals without worsening any RT or m/z tercile. The file is then searched again with
  corrected masses for targets and decoys alike. `precursor_ppm`, `fragment_ppm`, and
  `expmass` stay raw, `calibrated_*_ppm` report residuals, and `run-summary.json` records
  the chosen model, candidate scores, and residual bins per file under
  `models.mass_recalibration`.
- Spectra now carry an acquisition group (MS2 mass analyzer and activation) read from Thermo
  filter strings, mzML instrument configurations and activation terms, or fixed as TOF/CID
  for Bruker TDF. The Astral analyzer is read from the `ASTMS` filter or from mzML
  `MS:1003379`, as msconvert and ThermoRawFileParser write it. Activations are HCD, CID, ETD, EThcD, and ETciD (`etcid`, ETD with
  supplemental CID, kept apart from ETD). `mass_recalibration` selects fragment models
  separately per group, never corrects ion-trap groups, and leaves groups with too few PSMs
  uncorrected; the run summary lists each group's model.
- `tolerance_mode` (`"fixed"` by default, or `"auto"`) narrows ppm search tolerances from
  the discovery pass: the precursor tolerance per file and the fragment tolerance per file and
  acquisition group. Residuals after correction are fitted as a signal (Gaussian, Student-t
  with 4 degrees of freedom, or two Gaussians with a shared center, chosen by held-out
  likelihood) plus a uniform background, and the window holds 99% of the signal (floors of
  3 ppm for precursors and 5 ppm for fragments). Precursor PSMs matched at a non-zero isotope
  error are fitted as their own subgroup and the window covers both. Windows are never wider
  than configured and are applied only when they hold 98% of the estimated held-out signal.
  Untrusted fits, Da and percent tolerances, ion-trap groups, sparse groups, mass-offset
  searches, and open precursor windows keep the configured tolerance. The run summary reports
  each fit (form, signal fraction, center, scale, 99% half-width, candidate likelihoods), the
  applied window or why it was kept, and held-out coverage.

### Changed
- Nonlinear retention-time alignment is the default. When `retention_time_alignment` is not
  set, RT prediction and LFQ align runs with a robust affine fit followed by a monotone
  piecewise-linear warp; runs with fewer than 16 shared landmarks fall back to the robust
  affine fit, and single-file searches are unchanged. Set `"retention_time_alignment":
  "linear"` for the previous ordinary least-squares alignment, which reproduces earlier
  results byte for byte. On five PXD028735 LFQ files, LFQ precursors at 1% q-value rose from
  24,038 to 31,418 and those quantified in all five files from 23,055 to 30,828; the spread of
  human log2(A/B) ratios (MAD) fell from 0.213 to 0.191, and the median cross-run aligned-RT
  residual fell from 0.128 to 0.077 min. PSMs, peptides, and proteins at 1% FDR changed by
  less than 0.1% (312,533 to 312,487 PSMs), and runtime was unchanged within noise.
  `run-summary.json` already records the method; its schema stays at version 9.
- The retention-time and ion-mobility models are fit on standardized features with a fixed,
  tiny ridge penalty and a Cholesky solve, replacing Gaussian elimination that retried with a
  growing perturbation. The default feature sets contain exactly redundant columns; fitted
  values no longer depend on how those columns are laid out. Additive PTM offsets keep their
  `ptm_regularization` penalty and use the same Cholesky solve. On benchmark data, predictions
  change by rounding only (a few rare-terminus peptides with physicochemical mobility features
  by up to 0.017), and identification counts are unchanged or within 5 PSMs.

### Fixed
- MS2 spectra without a precursor are skipped by the search, the discovery pass, and the
  prefilter spectrum index. Thermo RAW files whose scan events are out of step with their
  trailers can contain MS2 scans without a plausible precursor m/z; the reader already
  reported that they would not be searched, but they were searched and the run aborted with
  `missing MS1 precursor` (also in Beta 8).
- Thermo RAW files whose scan events are out of step with their trailers no longer take
  acquisition groups from those events. Their filter strings belong to other scans, so an
  Orbitrap HCD scan could be labelled ion trap and skipped by fragment recalibration. Such files
  now report every scan as `unknown/unknown`, with a warning, since the trailers do not name
  the analyzer or activation.

## [v0.1.0-beta.8] - 2026-09-24

### Fixed
- Labeled searches no longer panic when a channel has a non-zero base mass and a non-zero
  offset, such as dimethyl labels. The base definition is recovered by identity instead of by
  recomputing its mass in single precision.
- FASTA files containing `X`, `B`, `Z`, or `J` load again, as in upstream Sage. Peptides that
  contain residues without a defined mass are skipped, and one warning reports how many
  proteins contain them. Lowercase letters, digits, and symbols are still rejected.
- LFQ with `mbr: false` quantifies each file against its own identification. Files previously
  shared one retention-time grid whose centre and reference file depended on thread timing, so
  intensities were wrong and changed between runs.
- Protein-group FDR pairs each decoy with the group of the target protein it was reversed
  from. Decoys were never grouped, so no decoy group key matched a target group and picked
  competition reduced to plain target-decoy counting. On five PXD028735 LFQ files, protein
  groups at 1% FDR rose from 7,046 to 7,177; PSM, peptide, and protein results are unchanged.
- Linear regression for the retention-time and mobility models and the kernel density estimate
  for PEP and q-values sum in a fixed order. Repeated runs, including multi-file searches, now
  produce byte-identical results.
- Required neutral losses on mass offsets shift preliminary fragments by the smallest loss,
  matching the fragment index, regardless of the order losses are listed in.
- `mass_shift` and `ambiguity_sequence` exclude the matched isotope error. PSMs at a non-zero
  isotope error previously reported a spurious shift of about 1.003 Da per isotope.
- Ion-mobility residue-class features count the intended residues.
- Thermo RAW MS levels and precursors come from the scan trailers wherever OpenTFRaw's scan
  event contradicts them, as on some Orbitrap Fusion files where events are decoded out of step
  with the scans. Once any scan is contradicted, every scan with a trailer master takes its level
  and precursor from the trailer. Otherwise dependent scans follow their master scan chain, and
  scans without a master keep an MSn event only when it has a plausible precursor, as DIA and
  targeted scans do. The trailer chain uses only "Master Scan Number"; "Master Index", a 0/1 flag
  on QE and LTQ Orbitrap files, is ignored. An SPS-MS3 TMT file
  previously read as 104,433 MS1 scans and 1 MS2 scan; it now matches msconvert (14,239 MS1,
  45,139 MS2, 45,056 MS3), and MS3 scans keep their MS2 parent for reporter ions. MS2 scans
  without a plausible precursor m/z (50-20,000) are reported with a warning and are not
  searched.
- Non-finite peak and precursor m/z values are dropped during spectrum processing, so they no
  longer abort prefiltered searches.
- Calibrated 1/K0 finds `analysis.tdf` for Bruker inputs given as a path inside the `.d`
  directory, such as `analysis.tdf_bin`. Inputs without `analysis.tdf`, such as miniTDF `.ms2`
  directories, fall back to the linear scale with a warning, and `run-summary.json` records the
  scale applied (`mixed` when only some inputs fell back).
- mzML array type, compression, precursor, isolation-window, and noise state are reset for
  every array, precursor, and spectrum, so values no longer leak between them.
- MGF spectra without a title, a valid precursor mass, or finite peaks are skipped with one warning per file
  instead of failing the file. File-level `CHARGE`, `TOL`, and `TOLU` apply to the first
  spectrum, multi-digit charges parse correctly, and an empty `CHARGE=` means unknown charge.
- Spectra that share an ID within a file are annotated and localized against their own peaks
  after FDR, and MS2 TMT reporter ions stay with their own spectrum, instead of failing or
  borrowing another spectrum's features. `TmtQuant` gains an `occurrence` field and
  `serialize_features` takes the repeated-spectrum occurrence map.
- File progress events are emitted once per file; the prefilter and post-FDR rereads no longer
  repeat them.
- The HTML report keeps inputs with the same file name in different directories separate.
- The MCP `query_results` modification filter matches inside bracketed tags of PSM and
  spectral-library sequences (plain residue letters no longer match), and a job that finished
  before a cancellation request is reported as completed.
- `DOCS.md` states the correct `predict_rt` default (true).

### Added
- Motif modification sites. A `motif:` site in a named definition restricts a residue
  modification to a PROSITE-style pattern, with `*` marking the modified residue, for example
  `motif:N*-{P}-[ST]` for N-glycosylation sequons or `motif:R-x(2)-[ST]*`. Motifs are matched
  against each peptide's source protein, so they can extend past the peptide. Shared peptides
  are eligible if any occurrence matches, generated decoys mirror their target's sites, and
  localization only considers matching positions. Motifs work with static, variable,
  mass-offset, and library-restricted definitions. `--preview-before` and `--preview-after`
  supply protein context to `--preview-modifications`. Configurations without motifs produce
  identical results.

### Changed
- Prefiltered searches are 1.06 to 1.40 times faster, with byte-identical results. Spectra read
  by the prefilter are kept for the search instead of being read and processed again, up to the
  prefilter's spectrum index budget. The database memory estimate runs in parallel, and prefilter
  survivors are no longer sorted a second time. On five PXD028735 LFQ files the search took
  105-109 s instead of 133-140 s, and peak memory rose from 6.47 to 6.67 GiB.
- Update timsrust from 0.6.5 to 0.6.6, with the vendored `filemanager` patch rebased on 0.6.6.
  timsrust 0.6.6 changed its uncalibrated scan-to-1/K0 conversion from linear in sqrt(1/K0)
  to linear in 1/K0. The default calibrated scale does not use it. `ion_mobility_scale:
  "linear"` now computes the earlier conversion in Sage Plus and still reproduces Beta 6.
  Spans of the opt-in `UniformMobility` DIA window splitting follow timsrust's new scale.
- Update OpenTFRaw from 1.4.0 to 1.4.1. The fix affects profile spectra only, and Thermo RAW
  search output is unchanged.

## [v0.1.0-beta.7] - 2026-09-23

### Fixed
- Bruker timsTOF ion mobility now applies each frame's `TimsCalibration` model and matches
  the Bruker SDK. timsrust interpolated between the acquisition limits, which differed from
  the calibrated scale by up to 0.054 1/K0 on a PXD070049 run and 0.12 1/K0 on example runs. MS1 peaks are converted before
  centroiding, and DDA precursors use their fractional average scan instead of a truncated one.
  Unsupported calibration models stop the search.

### Added
- `bruker_config.ion_mobility_scale` selects `calibrated` (default) or `linear`, which reproduces
  earlier releases byte for byte. `bruker_config.ms1` and `ms2` may now be omitted. `run-summary.json` reports the scale as `models.ion_mobility_scale` for
  Bruker inputs. The field is optional and the run-summary schema stays at version 9.

### Changed
- Database prefiltering indexes the spectra and streams generated peptides through the
  spectrum index. It no longer builds a fragment index for every database chunk or searches
  every spectrum against each one, and spectra are read once instead of once per chunk.
  Retained peptides, and therefore search results, are unchanged.
- The spectrum index stores each peak once and computes tolerance and mass-offset windows at
  lookup. Wide precursor windows switch to a global peak index by estimated lookup cost.
- When the spectrum index would exceed a quarter of `max_memory_gb` (8 GiB when no limit is
  set), spectra are indexed in batches and the database is streamed once per batch.
  `SAGE_PREFILTER_INDEX_GB` overrides the budget.
- Decoy-pair closure and survivor collection after prefiltering run in parallel.
- The memory preflight no longer budgets a fragment index per prefilter chunk.
- Update DashMap from 5.5.3 to 6.2.1. Searches, LFQ, and cross-run retention alignment produced
  identical results with both versions.

## [v0.1.0-beta.6] - 2026-09-20

### Added
- Named static and variable modification dictionaries with one definition and shared
  limits across explicit site strings. Terminal groups, boundary residues, internal
  residues, and residue-conditioned terminal groups have distinct spellings.
- Configuration migration through `--migrate-modifications` and library-aware
  modification preview with protein accession and one-based peptide coordinates.
- Typed terminal-group localization and version 2 PTM libraries and site reports.
  TSV and Parquet outputs preserve attachment identity alongside protein coordinates.

### Fixed
- Library evidence for a residue cannot license an adjacent terminal-group attachment.
- Localization obeys library-only site restrictions and excludes indistinguishable
  attachment alternatives from reusable libraries.
- Shared site matching covers static application, indexed and mass-offset placement,
  memory estimation, localization, and retention and mobility feature counting.
- Conflicting overlapping static definitions fail validation instead of depending
  on hash-map iteration order.

### Compatibility
- Legacy symbol-keyed configuration remains readable and has an explicit migration path.
  Four-column PTM libraries are read as residue attachments. Legacy terminal evidence
  requires an attachment column, and strict four-column consumers need an update.
- The Beta 5 tag was retained without publication. Beta 6 includes those changes.
- No dependency upgrades are included. Terminal localization remains experimental.

## [v0.1.0-beta.5] - 2026-09-20

### Added
- Internal residue specificity such as `~K` for static, indexed variable, and
  search-time mass-offset modifications. Combine positional keys for the same
  named modification while preserving shared occurrence limits.
- `--preview-modifications PEPTIDE` prints eligible sites, limits, and bounded
  generated variants as JSON. Protein-boundary context and output limits are
  configurable, without loading spectra or writing search outputs.
- A synthetic positional-modification benchmark covers H3K9 at the first peptide
  residue, internal lysines, peptide-terminal and protein-terminal placements,
  and agreement between indexed and mass-offset search.

### Fixed
- Localization retains residue-position restrictions and full modification
  definitions, so equal-mass named modifications are not pooled and other
  modifications remain fixed during relocation.
- Modification placement, memory estimation, and retention and mobility feature
  counting consistently recognize internal residues.

### Changed
- Malformed modification keys now fail configuration loading instead of being
  logged and omitted. Runtime and JSON schema validation enforce the same syntax.

## [v0.1.0-beta.4] - 2026-09-19

### Added
- Variable modifications accept `search_mode`. `"mass_offset"` searches a modification as
  a search-time precursor and fragment offset instead of expanding it into the fragment
  index, placing at most one copy per peptide. Offset candidates compete with ordinary
  candidates in scoring, FDR, quantification, localization, and PTM site libraries, and
  are restricted to library sites under `site_mode: "library"`. Offsets are never combined
  with each other, so cost grows linearly with the number configured while index size and
  database-build memory follow the indexed modifications only. Run summaries report
  `mass_offset_definitions`, `mass_offset_psms`, and `mass_offset_peptidoforms`.
- Label-free quantification records strict MS2 confirmation and experimental per-file
  signal diagnostics (spectral angle, trace cosine, retention-time shift, and an
  uncalibrated ranking score) in new `lfq.v3` and `lfq.v4` Parquet schemas.
- `benchmarks/run_mass_offset.py` runs the mass offset evaluation matrix with recorded
  executable, configuration, input, and output identities.

### Security
- Update `rustls` to 0.23.45, resolving RUSTSEC-2026-0285, where TLS 1.3 handshake
  messages were accepted across encryption level boundaries.

### Changed
- Peptide and protein confidence fall back to count-based target-decoy q-values when the
  kernel density model is underdetermined or returns a non-finite posterior. Precursor
  q-values use the same tie-preserving counter.
- PTM localization reports no confidence for a modification whose best target arrangements
  score equally, instead of accepting one arbitrary site. Those arrangements remain in the
  false-localization-rate competition population.

## [v0.1.0-beta.3] - 2026-09-10

### Fixed
- API and MCP searches honor the configured file batch size and explicit request overrides.
- Gzip outputs finish their compression stream before returning success.
- Concurrent event output preserves sequence and elapsed-time order. Active MCP readers defer incomplete trailing records.
- Logical validation rejects unsupported percentage precursor tolerances and invalid numeric model settings before searching.
- Cancellation is checked at postprocessing and completion boundaries. MCP job records are replaced atomically.
- Mass-alignment summaries report actual per-file fits and skipped fits.

### Added
- Run-summary schema 9 records warnings, effective worker count, build feature selection, and explicitly labeled input metadata identity.
- Verified benchmark stage manifests, complete experiment checks, and typed run comparisons.
- Minimal-feature and benchmark tests in CI, macOS and Windows worker smoke tests, and scheduled dependency auditing.

### Changed
- Local output directories containing Sage artifacts require `--overwrite` or `"overwrite": true`. Explicit overwrite removes known Sage artifacts while preserving unrelated files. Use a separate directory for each concurrent run and a fresh prefix for remote outputs.
- Updated patched dependencies for crossbeam-epoch, h2, anyhow, memmap2, and direct quick-xml usage. Updated timsrust to 0.6.5 and chacha20 to 0.10.2. A pinned, licensed filemanager compatibility patch uses object_store 0.14.1 to remove the remaining vulnerable XML dependency.
- Consolidated maintenance updates for anyhow 1.0.104, clap 4.6.1, itoa 1.0.18, and schemars 1.2.2. Coverage artifact uploads now use the same pinned upload-artifact v7.0.1 as release packaging. DashMap's major upgrade remains deferred for separate validation.
- mzML binary arrays containing XML entity references now fail explicitly. Literal base64 arrays remain supported.

### Validation scope
- Representative HEK searches and a bounded 20-seed entrapment comparison preserve identification outcomes. This is regression evidence on one dataset and a reduced FASTA, not broad scientific calibration. Expanded independent-study validation remains planned. See `benchmarks/HARDENING_RESULTS.md` and `benchmarks/BETA3_RELEASE.md` for evidence and release status.

## [v0.1.0-beta.2] - 2026-08-28

### Added
- Standard builds and release binaries now read local mzMLb files.
- mzML parsing now resolves referenceable parameter groups used for shared spectrum and instrument metadata.
- Canonical PSM Parquet output now includes typed protein occurrences with one-based inclusive start and end positions plus flanking residues.
- LFQ can disable match-between-runs with `quant.lfq_settings.mbr`, and ion-mobility model fitting can be disabled independently with `ion_mobility_model.enabled`.
- A committed JSON configuration schema supports editor completion, static validation, and export with `sage --write-config-schema`.
- Searches now warn when no decoys are available for false-discovery-rate estimation.

### Changed
- **Breaking:** `Peptide::modifications` now uses `CompactModifications` instead of a dense `Vec<f32>`, and applied metadata is exposed through `Peptide::applied_modifications()`. Use `CompactModifications::from_dense`, `from_sparse`, or `from_applied` when constructing peptides directly.
- **Breaking:** Target peptide sequences now use immutable spans into shared protein storage. Generated decoys and pre-digested TSV peptides retain owned storage. Sequence equality, ordering, and hashing remain content-based.
- Peptide modifications now use two-byte position and definition IDs, modification definitions are shared across the database, and common single-protein mappings remain inline.
- The preliminary fragment index now stores lossless six-byte records. Exact floating-point mass prefixes are shared across bounded search buckets, and `database.bucket_size` remains the maximum number of fragment records in one bucket.
- Linux GNU libc builds return unused database-construction pages to the operating system before searching.
- Replaced MS2 deconvolution with averagine-scored isotope envelopes. Mass tolerance, charge ceiling, envelope length, isotope score, and intensity-ratio tolerance are configurable. Scored envelopes assign peaks exclusively and preserve charge confidence for charge-aware fragment matching. Boolean `deisotope` settings remain supported, with `true` selecting the scored defaults. mzML and mzMLb fragment charge arrays now constrain envelope assignment and flow into charge-aware matching.
- Clarified that Bioconda installs upstream Sage rather than Sage Plus.
- Added low-noise monthly dependency checks for Cargo and GitHub Actions.
- Removed the experimental bitmap preliminary search path and its configuration options.
- Updated network, compression, Parquet, Bruker, logging, and reporting dependencies to patched releases.
- Pinned CI actions to immutable commits and added dependency review for new pull requests.

### Removed
- Removed the experimental DDA spectral-library search mode. Empirical spectral-library export in Sage Parquet and PSI mzSpecLib formats remains supported.

### Fixed
- Release archives now build HDF5 and mzMLb support natively on each CPU architecture, including
  static musl builds inside matching official Rust Alpine containers and CMake 4 compatibility for
  bundled native libraries. macOS release builds also bypass an obsolete zlib 1.2.11 compatibility
  macro that conflicts with Xcode 16.3 and newer.
- Database prefiltering now filters targets and paired decoys together, preserves label-channel partners, and uses deterministic score ties. Prefiltered and full searches therefore use the same competition and FDR model.
- LFQ mass lookup now derives its coarse search margin from configured precursor ranges, avoiding missed matches at wider tolerances. Mobility boundaries remain inclusive and configuration-driven.
- Updated Bruker TDF parsing for the current `timsrust` API without dropping existing processing configuration.

## [v0.1.0-beta.1] - 2026-08-25

### Added
- **Library calibration and quantification:** Library search supports isotope-error windows, matched-fragment export, per-file precursor and fragment calibration, RT and mobility alignment, LFQ, and TMT.
- **Library rescoring:** A regularized target-decoy model combines spectral angle, explained intensity, isotope, mass-error, RT, and mobility evidence. Small searches fall back safely to spectral-angle scoring and emit a structured warning.
- **Consensus library generation:** Combine supporting PSMs with median RT, mobility, and normalized fragment intensities. Minimum PSM support and fragment frequency are configurable.
- **Modification-defined label channels:** Structured static and variable modifications can define coherent `channel_offsets` for SILAC, dimethyl, and custom precursor labels.
- **Channel-aware LFQ and outputs:** Identified channels seed exact-mass partner extraction and reference ratios. Version 2 results, LFQ, and spectral-library schemas preserve channel and label-group metadata through export.
- **MCP worker isolation:** Searches run in dedicated processes. Panics, signals, and memory-limit exits fail only the affected job, while cancellation escalates from a cooperative request to worker termination.
- **Library FDR validation:** Added a held-out and entrapment validation guide, a reporting script, and warnings when query files overlap library source files.

### Changed
- **Breaking:** Replaced `database.labels` with `channel_offsets` on structured static and variable modifications. Placement controls occupancy, while channel offsets define coherent precursor states.
- Pre-digested peptide TSV input now receives configured static, variable, and channel-aware modifications like FASTA-derived peptides.
- Embedded runner jobs now report memory-limit failures through cancellation and structured events. CLI and MCP worker processes retain exit-code 137 protection.
- Moved unit tests out of production modules and added an 80% workspace line-coverage requirement in CI.

### Fixed
- Assign equal-score PSMs as one FDR threshold group and use deterministic feature ordering for repeatable PSM identifiers and library selection.
- Preserve site-specific modification names in ambiguity annotations when labels share a mass, and serialize LFQ rows deterministically.
- Reject nonfinite, negative, or out-of-range LFQ tolerances and thresholds instead of silently normalizing them.
- Stop searches on empty, unreadable, or malformed spectrum files and emit a structured `file_failed` event. MGF errors now identify invalid fields, missing markers, and line numbers.
- Skip unavailable HTML QC plots and nonfinite plot values instead of failing report generation.
- Skip invalid historical MCP job records during server startup and report invalid local output URLs directly.
- Use portable release-version checks that do not require `ripgrep`.

### Earlier additions
- Additive protein-specific cleavage-site TSV or Parquet input with optional FASTA sequence-context validation.
- Optional `max_memory_gb` and `min_free_memory_gb` safeguards that preflight unmodified peptides, variable-modification expansion, and fragment indexes, then monitor the running process to stop Sage before configured limits are crossed.
- A top-level `batch_size` configuration option; the existing `--batch-size` command-line option takes precedence.
- **Pre-digested peptide file input** (`database.peptides`): supply a TSV of peptide sequences instead of (or in addition to) a FASTA file.
  - TSV must have a header row with a required `sequence` column; optional `protein` and `decoy` columns are also supported.
  - Ordinary variable and static modifications are not applied from config. Precursor label channels are applied when configured.
  - Decoys are still generated automatically by reversal when `generate_decoys: true`.
- If both `database.fasta` and `database.peptides` are provided, peptides from both sources are merged and deduplicated before building the index.
- Supports cloud paths (S3, GCS) via the same mechanism as FASTA loading.
- Per-modification occurrence limits using `{"mass": <mass>, "max_count": <limit>}` entries in `database.variable_mods`; existing bare-mass entries remain supported.
- `database.max_combinations` to cap the number of peptide variants (including the unmodified form) generated from variable modifications, preferring variants with fewer modifications.
- Optional names for structured static and variable modifications; named modifications are preserved on exact peptide sites and rendered in peptide output.
- Neutral-loss fragment variants for structured modifications, with `neutral_loss_mode` controlling whether retained fragments are optional or suppressed.
- Neutral-loss annotations in matched-fragment TSV and Parquet output.
- Deterministic empirical spectral-library generation from FDR-filtered target PSMs, with
  canonical long-form Parquet and PSI mzSpecLib text output.

### Earlier changes
- **Breaking:** Parquet is now the canonical analytical output. Sage no longer emits TSV variants of PSM, LFQ, matched-fragment, or PTM-site result tables, and the `--parquet` option has been removed. Purpose-specific `.pin`, JSON, HTML, and PTM-library interchange outputs are unchanged.
- LFQ remains a separate long-form `lfq.parquet`. Missing integrated signals are encoded as Parquet nulls instead of zero; each precursor/file row now includes `ms2_confirmed`, plus the LFQ peak `score` and `spectral_angle`.
- MCP result queries now scan canonical Parquet result files and return typed JSON values.
- PSM and matched-fragment Parquet output now share `output_filter.psm_q_value` (default `0.1`), an inclusive spectrum-level q-value cutoff. The effective cutoff is embedded in both files' Parquet metadata; setting it to `1.0` retains all scored PSMs.
- Matched-fragment details are now reconstructed only for retained PSMs in a batched post-FDR MS2 pass. This removes fragment-detail allocation from candidate scoring, preserves rank-ordered chimera peak removal, and shares spectrum rereads with PTM localization when both are enabled.

The entries below predate independent Sage Plus versioning and are retained from upstream Sage.

## [v0.15.0]
### Added
- Memory guard: a background watchdog that terminates Sage cleanly (exit 137) before it can exhaust system RAM and freeze the host, rather than letting the OS OOM-killer or swap-thrash take over. Aborts if Sage's resident memory exceeds a ceiling (default 90% of total RAM) or if system available memory falls below a safety floor. Configurable via `--max-memory <GiB>` / `SAGE_MAX_MEMORY_GB`; `--max-memory 0` disables it. Polls from a single thread with no allocation-hot-path overhead.
- `database.max_peff_variable_mods`: an optional, independent per-peptide budget for PEFF/positional modifications. When unset, PEFF mods continue to share the `max_variable_mods` budget with global variable mods (no behavior change); when set, the two pools get independent budgets, so a peptide may carry up to `max_variable_mods` global mods *and* up to `max_peff_variable_mods` PEFF-annotated mods simultaneously (e.g. cap global mods at 1 while allowing 3 PEFF mods).
- Sequence-ambiguity annotation, computed natively for every search (a port of the standalone SagePeptideAmbiguityAnnotator tool). Two new columns are added to `results.sage.tsv`/`.parquet`: `ambiguity_sequence`, the peptide string with residue runs lacking flanking fragment-ion evidence wrapped in `(?...)` and any residual precursor mass shift placed (`[+mass]` localized, `(...)[+mass]` regional, or leading `{+mass}` labile); and `mass_shift`, the residual `expmass - calcmass` (0.0 within `mass_shift_ppm`). The threshold for deciding a real shift is configurable via the top-level `mass_shift_ppm` parameter (default: 50.0), independent of `precursor_tol`.
- PTM site localization (`ptm_localization` setting / `--localize` flag). After spectrum FDR assignment, each passing target PSM carrying a variable modification is re-scored across candidate sites using site-determining ions. Balanced impossible-site decoy competition provides arrangement-level localization q-values, filtered by `ptm_localization.localization_q_value` (default: 0.01), alongside AScore-style deltas and per-site probabilities. Keeping localization out of the main scoring hot path avoids work on rejected PSMs. Site reports are written in the selected TSV or Parquet format.
- IDPicker-based protein grouping with picked group FDR control (`protein_grouping` setting, enabled by default). Proteins are grouped using a bipartite graph greedy set cover approach, and protein group-level q-values are reported via target-decoy competition. New output columns: `protein_groups`, `num_protein_groups`, `protein_group_q`.
- `protein_grouping_peptide_fdr` parameter to control the peptide FDR threshold used for confident peptides during protein grouping (default: 0.01)
- Initial support for LFQ on data with ion mobility.
- Speedup on the generation of databases when large number of peptides are redundant.
- Initial support for searching diaPASEF data
- `override_precursor_charge` setting that forces multiple charge states to be searched
- Cross-cloud storage support: Replaced AWS SDK with the `object_store` crate. Sage now natively supports reading/writing from Amazon S3, Google Cloud Storage, and Azure Blob Storage using `s3://`, `gs://`, and `az://` URL schemes.
- HTML QC report generation (`--write-report`)
- `lfq_settings.peptide_q_value` parameter for controlling which peptides are quantified
- Stack size configuration for Rayon threads
- Allow zero or multiple amino-acids as cleavage restrictions
- CITATION.cff file

### Fixed
- Handle negative mass errors for .pin files
- Extract precursor m/z from 'isolation window target m/z' if missing
- Selected ion m/z of 0.0 was overwriting precursor.mz
- C/N-term mixup in modification handling
- Bruker `.tdf` filename handling (use parent directory name)
- Performance optimizations on prefiltering
- Picked protein FDR now correctly uses only proteotypic peptides for competition

### Breaking Changes
- `precursor_ppm` field reports the non-absoluted average mass error, rather than the absoluted average mass error.
- Don't deisotope reporter ion regions if MS2-based TMT/iTRAQ is used
- Removed `fragment_min_mz` and `fragment_max_mz` parameters. These were decreasing the accuracy of preliminary scoring estimation when attempting to annotate multiply-charged, high-m/z ions.
- `sage-cloudpath` no longer exposes the `CloudPath` type. All paths are represented as URLs (`url::Url`). Local paths are converted to `file://` URLs internally.

## [v0.14.7]
### Added
- Added columns missing from parquet output: `semi_enzymatic` and `missed_cleavages`
### Changed
- Fixed ion mobility parsing from some mzMLs
- MGF paths were being lowercased prior to parsing

## [v0.14.6]
### Added
- Support for MGF files
- Support for writing ion mobility measurements to output files: `ion_mobility`, `predicted_mobility`, `delta_mobility` added to primary tsv and parquet reports. Ion mobility is predicted in a similar manner to RT, using a linear model trained on the data from the search.

## [v0.14.5]
### Added
- Support for semi-enzymatic digests (`database.enzyme.semi_enzymatic` parameter)
- Ability to directly export matched fragment ions (e.g. for spectral library or rescoring) with the `--annotate-matches` CLI option. This is compatible with the `--parquet` CLI option as well. Annotations will be written to `matched_fragments.sage.tsv` or `matched_fragments.sage.parquet`
- Sage sends basic telemetry data (version of Sage, run time, OS, # of CPU cores, # of peptides in database, whether LFQ is used) to a remote server. No information about your actual data is sent - e.g. identifications, quantities, organism, or modifications are NOT tracked or reported.  This data will be used to help focus efforts on improving Sage and figuring which features are most used. Please take a look at `crates/sage-cli/src/telemetry.rs` to see exactly what is sent! You can disable sending telemetry data  by using the `--disable-telemetry-i-dont-want-to-improve-sage` CLI flag.
### Changed
- Modified visibility on some crate internals to support the [sagepy project](https://github.com/theGreatHerrLebert/sagepy)
- Added `psm_id` field to various output files to match the new `--annotate-matches` option.
### Removed
- Removed the `ms1_intensity` field from CSV output, since it is essentially useless


## [v0.14.4]
### Added
- **Unstable feature**: Preliminary support for reading Bruker .d folders (ddaPASEF; no MS1/LFQ support yet)
### Changed
- Retention times are converted to minutes
### Fixed
- Fixed bug where charge state 1 would never be searched

## [v0.14.3]
### Fixed
- Hotfix for bug in parquet LFQ writer

## [v0.14.2]
### Added
- `quant.lfq_settings.combine_charge_state` boolean option. By default this is set to `true`, and LFQ is performed on the peptide-level, where all charge states are treated as the same precursor. Setting this to `false` performs LFQ on the peptide-charge-level, where each charge state will be treated separately.
### Changed
- Percolator output format now contains the integer-valued charge state encoded in the `z=other` column, if the charge state is outside the range 2-6 (e.g. a value of 7 will appear in the `z=other` column, rather than it being one-hot encoded)
- LFQ uses the the charge state range from the `precursor_charge` configuration option for tracing MS1 peaks

## [v0.14.1]
### Added
- Added additional output showing search progress if `SAGE_LOG=trace` environment variable is set
- Added additional warnings about precursor tolerances
- Added configuration option `precursor_charge` to make it explicit what charge states are being searched in the case where the mzML does not contain charge state information, or where `wide_window` is turned on.
### Changed
- Added a warning message if variable modifications are specified as single values (e.g. `15.9949`) instead of lists of values (e.g. `[15.9949]`). By v0.15 this will become a hard error and will not parse, to simply some of the internal logic.

## [v0.14.0]
### Added
- Support for parquet file format output. Search results and reporter ion quantification will be written to one file (`results.sage.parquet`) and label-free quant will be written to another (`lfq.parquet`). Parquet files tend to be significantly smaller than TSV files, faster to parse, and are compatible with a variety of distributed SQL engines.
### Changed
- Implement heapselect algorithm for faster sorting of candidate matches (#80). This is a backwards-incompatible change with respect to output - small changes in PSM ranks will be present between v0.13.4 and v0.14.0

## [v0.13.4]
### Fixed
- Bug in mzML parser, where some older specification-compliant mzMLs would not parse. If your mzMLs previously parsed, then there will be no change in behavior. Added a test case

## [v0.13.3]
### Fixed
- Bug in `database.enzyme.restrict` parameter, where `null` values were being overriden with "P" (causing Trypsin/P to behave like Trypsin)

## [v0.13.2]
### Changed
- Subtle change to TMT integration tolerance, and selection of which ion to quantify (most intense). As a result, TMT integration should be more in agreement (if not 100% so) with ProteomeDiscover/FragPipe/etc
- Remove `delta_mass` (precursor ppm) LDA feature - instead, build a delta mass (or ppm) profile using KDE/posterior error calculation code, and use the P(decoy) as a feature for LDA.

## [v0.13.1]
### Changed
- Internal performance and stability improvements for RT prediction & LDA

## [v0.13.0]
### Added
- Better error reporting thanks to @Elendol
- Added support for multiple variable mods for the same amino acid
- Added support for N/C-terminal modifications specific to an individual amino acid

New syntax:
```json
"variable_mods": {
    "M": [15.9949],
    "^Q": -17.026549,
    "^E": -18.010565,
    "[": 42.010565
}
```

Either a single floating point number (-18.0) or a list of floating point numbers ([-18.0, -15.2]) can be supplied as modifications. Support for single values may eventually be phased out to simplify the parser.

### Changed
- Changed "_fdr" columns to "_q" (e.g. "spectrum_q") in "results.sage.tsv" file
- Changed internal data representation of `Peptide` struct to allow for sharing of sequences (using `Arc`) among modified peptides
- Fragment index creation should now be faster

## [0.12.0]
### Added
- Add `wide_window` option to configuration file. This option turns off `precursor_tol`, instead using the isolation window written in the mzML file.
### Changed
- Changed internal calculation of precursor tolerances when searching with `isotope_errors`. The new version should be more accurate. This change also enables a significant boost to search speed for open searches.

## [0.11.2]
### Added
- Add rank & charge features to LDA
### Changed
- One-hot encode charge state information for percolator `.pin` files
- Change PSMId -> SpecId for Mokapot compatibility with `.pin` files

## [0.11.1]
### Added
- Support for additional fragment ion types, via the "database.ion_kinds" configuration option. Valid values are "a", "b", "c", "x", "y", "z"
### Changed
- Sort protein names alphanumerically for each peptide entry. This should enhance stability across runs, and fixes a bug with picked-protein group FDR
- Fix another bug where picked-FDR approaches assume internal decoy generation

### Changed
- Modify order of operations during deisotoping. Deisotoped peaks can contribute intensity to only 1 parent peak now, rather than potentially multiple parent peaks

## [0.11.0]
### Added
- Support for percolator output files (`--write-pin` CLI flag)
- Support for modifying file batch size (`--batch-size N` CLI flag)
- Add `delta_best` feature, which reports the delta hyperscore from the best match to current ranked PSM
- Add Sage version to `results.json` files

### Changed
- Breaking changes to `quant` section of the configuration file format
- Rename `delta_hyperscore` to `delta_next`
- Altered internal scoring algorithm. Rather than consider all MS2 peaks within a m/z tolerance window to be matches to a theoretical spectrum, consider only the closest peak. This should increase the accuracy of # of matched peaks, and subsequent scores
- Overhaul of chimeric scoring, `report_psms` can now be used to search for multiple chimeric spectra
- Completely overhauled the LFQ algorithm: added match-between runs, peak scoring using normalized spectral angle relative to theoretical isotopic envelope, target decoy scoring of MS1 integration
- Fixed bug in picked-peptide FDR that could lead to liberal FDR
- Fixed bug in picked-protein FDR that could lead to conservative FDR
- Fixed bug where using variable protein terminal (e.g. protein N-terminal acetylation) modifications could cause some determinism. This also improves the accuracy of peptide => protein assignment. Unfortunately this fix has performance implications, causing creation of the fragment index to take up to ~2x as long.

### Removed
- Remove `no-parallel` CLI flag, and `parallel` configuration file entry

## [0.10.0]
### Added
- Retention times are now globally aligned across files
- RT prediction is then performed on all files at once (on aligned RTs), rather than one file at a time - previously, there were many instances where some files in a search could not have RTs predicted, decreasing the effectiveness of delta_rt as a feature for LDA.

### Changed
- Peptide sequences within a protein are now deduplicated - previously, repeated peptides would be called multiple times for the same protein (e.g. num_proteins > 1 even if the peptide was unique)

## [0.9.4]
### Changed
- Fix issues with RT prediction (and occasionally LDA) that arise from 0's being present on the diagonals of the covariance matrix (small amount of regularization added)

## [0.9.3]
### Added
- Allow users to set minimum number of matched b+y ions for reporting PSMs (`min_matched_peaks`)

### Changed
- Internal code for calculating factorials

## [0.9.2]
### Added
- Added option for TMT signal/noise quantification, if noise values are present in mzML

## [0.9.1]
### Changed
- FASTA file path, JSON configuration file can now be specified as "s3://" paths, allowing Sage to run completely disk-free

## [0.9.0]
### Added
- Support for non-specific digests, N-terminal enzymatic digestion

## [0.8.1]
### Added
- `quant.tmt_level` configuration option to enable MS2 (or MSn) isobaric quantification

## [0.8.0]
### Added
- Support for protein N-terminal ('['), C-terminal (']') as well as peptide C-terminal ('$') modifications
- Support for k-combinations of variable modifications. This can be specified with the `database.max_variable_mods` parameter

## [0.7.1] - 2022-11-04
### Changed
- Fix bug with in silico digest: Logic around overwriting decoys with target sequences was incorrect peptides shared between targets/decoys were being annotated as decoy peptides but assigned to non-decoy proteins. We now make sure that they are assigned to non-decoy proteins and also annotated as target sequences.

## [0.7.0] - 2022-11-03
### Added
- Add support for user-specified enzymes to JSON file.  `database.enzyme.sites` and `database.enzyme.restrict` are limited to valid amino acids
- Sage can now search MS2 spectra without annotated precursor charge states. Default behavior is to search with z=2, z=3, z=4, and then merge the PSMs for scoring

### Changed
- Configuration file schema changed. `peptide_min_len`, `peptide_max_len`, `missed_cleavages` are now specified under `database.enzyme` in the JSON file
- Internal behavior of Sage was changed to enable deterministic searching
- Docker file changed from Alpine to Debian


## [0.6.0] - 2022-11-01
### Added
- Changelog
- `rank` column added to output file
- `database.generate_decoys` parameter, which turns off internal decoy generation. This enables the use of FASTA databases for SearchGUI/PeptideShaker

### Changed
- Base ProForma v2 notation is used for peptide modifications, i.e. "\[+304.2071\]-PEPTIDEM\[+15.9949\]AAC\[+57.0214\]H"
- `scannr` column now contains the full nativeID/spectrum title from the mzML file, i.e. "controllerType=0 controllerNumber=1 scan=30069"
- `discriminant_score` column renamed to `sage_discriminant_score` for PeptideShaker recognition
- `database.decoy_prefix` JSON option changed to `database.decoy_tag`. This allows decoy tagging to occur anywhere within the accession: "sp|P01234_REVERSED|HUMAN"
- Output file renamed:  `results.pin` to `results.sage.tsv`
- Output file renamed: `quant.csv` to `quant.tsv`
- Rename `pin_paths` to `output_paths` in results.json file


## [0.5.1] - 2022-10-31
### Added
- Support for selenocysteine and pyrrolysine amino acids

## [0.5.0] - 2022-10-28
### Added
- Ability to directly read/write files from AWS S3

### Changed
- Processing files in parallel processes them in batches of `num_cpus / 2` to avoid memory issues
- Fixed bug where `protein_fdr` was erroneously assigned to `peptide_fdr` output field
- Additional parallelization for assignment of PEP, FDR, writing output files

## [0.4.0] - 2022-10-18
### Added
- Label free quantification can be enabled by turning on `quant.lfq` JSON parameter 
- Commmand line arguments can be used to override configuration file

## [0.3.1] - 2022-10-06
### Added
- Workflow contributions from [@wfondrie](https://github.com/wfondrie).

### Changed
- Don't parse empty MS2 spectra

## [0.3.0] - 2015-09-15
### Added
- Retention time prediction
- Ability to filter low-number b/y-ions for faster preliminary scoring (`database.min_ion_index` option)
- Ability to toggle retention time prediction (`predict_rt`)
