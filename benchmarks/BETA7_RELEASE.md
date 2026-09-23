# Beta 7 release validation

Target: `v0.1.0-beta.7`, based on main after the published Beta 6 release.

## Scope

- Database prefiltering indexes spectra and streams generated peptides through the spectrum
  index. It no longer builds a fragment index for every database chunk or rereads spectra for
  every chunk.
- The retained peptides equal the Beta 6 exact prefilter exactly. The final search, scoring,
  and FDR are unchanged.
- The spectrum index is bounded by a quarter of `max_memory_gb` (8 GiB without a limit).
  Larger inputs are indexed in file batches, and the database is streamed once per batch.
  `SAGE_PREFILTER_INDEX_GB` overrides the budget.
- Decoy-pair closure and survivor collection run in parallel.
- The memory preflight no longer budgets a per-chunk fragment index during prefiltering.
- Bruker timsTOF 1/K0 applies each frame's `TimsCalibration` model (ModelType 2) instead of
  timsrust's interpolation between the acquisition limits. `bruker_config.ion_mobility_scale`
  selects `calibrated` (default) or `linear`, and `run-summary.json` records the scale. This
  changes reported mobility, mobility features, and quantification for timsTOF inputs only.
- The README replaces per-release sections with a feature comparison against upstream Sage.

`sage-cloudpath` now depends directly on `rusqlite` 0.35, which timsrust already used. DashMap
moves from 5.5.3 to 6.2.1 (superseding PR #8); the closed HEK, timsTOF LFQ, and five-file LFQ
searches gave the same results with both versions. No other dependency versions change. The configuration gains one optional key, and run summaries gain
one optional field while keeping schema version 9.

## Reproducible evidence

```shell
cargo test --workspace --locked
cargo test --workspace --no-default-features --locked
cargo clippy --workspace --all-targets --all-features --locked -- -D warnings
python3 benchmarks/run_prefilter.py --baseline /path/to/sage-beta6 --candidate target/release/sage \
    --config benchmarks/configs/prefilter-closed-mods.json \
    --config benchmarks/configs/prefilter-broad-ptm.json \
    --config benchmarks/configs/prefilter-broad-ptm-2.json \
    --config benchmarks/configs/prefilter-open.json \
    --config benchmarks/configs/prefilter-lfq-5file.json \
    --config benchmarks/configs/prefilter-offset-hcd1.json \
    --config benchmarks/configs/prefilter-offset-hcd2.json \
    --output /data/sage-plus-scientific/prefilter-20260923/matrix --repeats 3 --prefilter-off
python3 benchmarks/run_prefilter.py --baseline /path/to/sage-beta6 --candidate target/release/sage \
    --config benchmarks/configs/prefilter-lfq-5file.json \
    --output /data/sage-plus-scientific/prefilter-20260923/batched --repeats 1 \
    --env SAGE_PREFILTER_INDEX_GB=0.1
```

Mobility tests compare the calibration model with Bruker SDK values (`libtimsdata`
`tims_scannum_to_oneoverk0`) for four calibration rows from three example acquisitions,
including fractional scans and scans beyond both model limits. The fixture is
`crates/sage-cloudpath/tests/data/bruker/tims_calibration_sdk.json`, written by
`benchmarks/generate_tims_calibration_fixture.py`. Synthetic `analysis.tdf`
files cover per-frame row selection, unsupported models, and missing rows. The bundled DIA
example is read on both scales.

### timsTOF mobility on PXD070049

`LFQ_Ultra2_PASEF_15min_50ng_Ecoli_01.d` (DDA-PASEF, one calibration row, 15,054 frames) was
searched against the E. coli reference without its seven `X`-containing entries, with LFQ and
retention-time prediction enabled. Evidence is retained under
`/data/sage-plus-scientific/mobility-calibration-20260923/`.

- The model matched `libtimsdata` to 4.4e-16 1/K0 at every half scan. The uncalibrated scale
  differed from it by up to 0.054 1/K0.
- All 17,565 reported precursor 1/K0 values on the calibrated scale equal the SDK conversion of
  a fractional precursor scan. On the linear scale, 113 did.
- `ion_mobility_scale: "linear"` produced `results.sage.parquet` and `lfq.parquet` byte-identical
  to Beta 6.
- Calibrated precursor mobility shifted by a median of 0.045 and at most 0.054 1/K0.
  Identifications changed little: 16,837 against 16,841 PSMs, 8,872 against 8,874 peptides,
  and 1,297 proteins in both at 1% FDR. LFQ reported 8,602 against 8,619 precursors. The median
  absolute mobility-model residual was 0.0146 against 0.0130 1/K0.

The calibrated scale makes reported mobility comparable with Bruker software. It is not
expected to increase identifications.

Unit tests compare the spectrum-indexed survivors with the Beta 6 `exact_prefilter` for closed
ppm and Da tolerances, unknown and overridden charges, wide windows, open searches, and mass
offsets with required neutral losses. Each case runs both lookup paths and batched indexes.
Results and interpretation are in [PREFILTER.md](PREFILTER.md).

## Release gates

- [x] Survivor-equivalence unit tests for every precursor and fragment lookup path.
- [x] Paired benchmark matrix with Beta 6, including prefilter-off and forced spectrum batches.
- [x] Mobility model agreement with the Bruker SDK, including extrapolation.
- [x] Complete workspace and minimal-feature tests, formatting, and strict Clippy.
- [x] Calibrated versus linear search on a PXD070049 DDA-PASEF run.
- [ ] Release version check and optimized locked build.
- [ ] Hosted required checks and preparation pull request merge.
- [ ] Manual packaging and native archive inspection.
- [ ] Annotated Beta 7 tag and verified publication with all assets.
- [ ] Paper refresh against the published Beta 7 executable.

## Paper refresh after publication

The paper compares upstream Sage with the latest published Sage Plus release, so it is refreshed
after the Beta 7 assets are verified.

- Replace the exact-prefilter figure and table row, which report the Beta 2 time-for-memory
  tradeoff, with prefilter-on, prefilter-off, and Beta 6 measurements from the Beta 7 executable.
- Record the evidence under `paper/analysis/data/` with executable and input hashes, and
  generate every reported number through `gen_stats.py`.
- Revise the storage and execution text in `paper.typ` to describe spectrum-indexed prefiltering
  and its memory budget.
- Report timsTOF mobility on the calibrated scale wherever the paper uses Bruker inputs, and
  state that earlier releases used timsrust's uncalibrated interpolation.
- Rerun the release comparison searches with the Beta 7 executable, as the paper requires when
  the evaluated executable changes, and update the recorded versions and hashes. Search results
  are expected to match Beta 6, and any difference must be explained before publication.
