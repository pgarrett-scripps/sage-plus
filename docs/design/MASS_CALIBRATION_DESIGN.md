# PSM and fragment mass-tolerance correction

## Recommendation

The first implementation is a one-search, per-file alignment of precursor and
fragment error features before final FDR fitting. It uses provisional 1% q-value
rank-1 target PSMs, applies the fitted correction equally to targets and decoys,
and preserves the original output errors.

The calibration primitive in `crates/sage/src/mass_calibration.rs` supplies the
robust static/RT-linear model and conservative model selection. It is wired into
FDR features, but deliberately does not change candidate or fragment matching.

## Why correction must happen during search

Sage uses `precursor_tol` to retrieve peptide candidates and `fragment_tol` in
both preliminary candidate scoring and full hyperscoring. Correcting only the
reported `delta_mass` or `average_ppm` features cannot recover candidates or
fragment matches that fell outside those windows.

A future search-correction workflow would therefore be:

1. Search each raw file with discovery windows (the configured windows, or a
   documented wider calibration window).
2. Select high-confidence, rank-1 target PSMs using target-decoy q-values.
3. Split calibration PSMs deterministically into fit and validation sets.
4. Fit separate precursor and fragment models per file.
5. Accept a model only if it improves validation residuals and has enough
   observations and RT coverage.
6. Search the already-loaded spectra again using corrected experimental masses
   and the user's original tolerances.
7. Report both raw and residual mass errors plus the fitted model diagnostics.

`Runner::process_chunk` already retains processed spectra while searching, so a
second scoring pass need not repeat file I/O. It will, however, approximately
double database-search CPU for calibrated files. Searching only a representative
subset of spectra in the discovery pass can reduce that cost later.

## Error definitions and fitting units

Use one signed convention everywhere:

```text
error_ppm = (observed - theoretical) / theoretical * 1e6
corrected_observed = observed / (1 + predicted_error_ppm * 1e-6)
```

The existing fragment feature is an absolute, intensity-weighted error, so it
cannot be used to learn the direction of correction. Full scoring should also
accumulate a signed fragment error. To avoid giving long/highly fragmented PSMs
disproportionate influence, derive one robust signed fragment estimate per PSM
(for example, the median of its matched-ion ppm errors), then fit across PSMs.

The precursor point is the isotope-adjusted observed precursor mass against the
theoretical peptide mass. Exclude open-search mass shifts and ambiguous isotope
assignments from the calibration set.

## Applying the model in Sage

The least invasive integration is to add per-file precursor and fragment
`CalibrationModel`s to `Scorer` and look them up with `query.file_id`.

- Candidate retrieval: correct the experimental precursor and fragment masses
  before `IndexedDatabase::query` and `IndexedQuery::page_search`.
- Full scoring: center fragment matching on the expected observed mass, or
  search a corrected mass view. `select_most_intense_peak` already supports a
  static Da offset, but a corrected mass view gives one consistent convention.
- Chimeric removal: use the same corrected matching logic as full scoring.
- Output: preserve raw `expmass`, add corrected/residual precursor ppm, signed
  raw/residual fragment ppm, model kind, and predicted offsets. Do not silently
  change the meaning of existing output columns.

Calibration should be disabled for wide-window/DIA searches initially. Their
precursor coordinates describe isolation windows rather than a single precise
analyte, while fragment calibration can be revisited independently.

## Static versus RT-linear

A static median is the right MVP: it is robust, easy to diagnose, and addresses
the common whole-run offset. RT-linear drift is cheap and interpretable, but it
should only replace static correction when all of these hold on validation data:

- at least 200 high-confidence PSMs per file (more for short/noisy runs),
- at least 20 minutes or a substantial fraction of the run is represented,
- median absolute residual improves by at least 10-15%, and
- early-, middle-, and late-run bins all improve rather than trading one region
  for another.

Residual diagnostics should also be plotted against m/z. Published calibration
methods often model both RT and m/z; a strong remaining m/z trend is evidence
that an RT-only line is underfit. A later model could use a small 2-D RT-by-m/z
grid with shrinkage toward the static offset, but that is not necessary for the
first experiment.

## Is it worth it?

Probably, as an opt-in feature with benchmark gates—not yet as an unconditional
default.

Expected benefit is highest for long gradients with instrument drift, older or
poorly lock-mass-calibrated data, TOF data, cross-run/open-search work, and users
who want to tighten fragment windows. On well-calibrated modern Orbitrap data
searched at roughly +/-10 to 20 ppm, a static correction may mainly sharpen mass
features rather than add many identifications. Fragment correction can still add
matched ions and improve localization or discrimination even when precursor
candidate recall is already saturated.

Advance the feature only if a multi-instrument benchmark shows, at fixed 1% FDR:

- no loss of accepted PSMs on any well-calibrated control run,
- repeatable PSM/peptide gains on drifted runs,
- lower held-out precursor and fragment median absolute residuals,
- stable target-decoy calibration and modification localization, and
- acceptable runtime (initial target: less than 1.8x total wall time by sampling
  the discovery pass).

The critical ablations are: no correction, static-only, RT-linear with
validation gating, precursor-only, fragment-only, and both. Also compare fixed
user tolerances against narrower post-calibration tolerances; much of the value
may come from improved specificity rather than raw PSM count.
