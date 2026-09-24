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

## Search-time implementation (Beta 9)

`mass_recalibration` (`off` by default, `static`, `linear`, `auto`) now corrects
masses during the search, in `crates/sage/src/mass_recalibration.rs` and
`Runner::process_chunk`:

1. Discovery: search at most 25,000 evenly spaced MS2 spectra per file with the
   configured windows, and keep rank-1 target PSMs at 1% spectrum q-value.
   Precursor points are monoisotopic matches without mass offsets; fragment
   points are every matched ion of those PSMs.
2. Split PSMs 7:3 into fit and held-out sets by a hash of the spectrum id, and
   drop points more than 4 robust sigmas from the fit-set median.
3. Walk none -> static -> linear -> smooth (additive RT and/or m/z terms). A
   step is kept only if the held-out median absolute residual improves (2% for
   static, 5% afterwards) and no RT or m/z tercile gets more than 2% worse.
   Minimum fit PSMs are 50, 200, and 1000. Predictions are capped at the
   tolerance half-width.
4. Refit the chosen form on all inliers, and search the whole file again with
   corrected precursor and fragment m/z. Targets and decoys see the same
   correction; output keeps raw errors plus `calibrated_precursor_ppm`.

The precursor model is per file. Fragment models are per file and acquisition
group (mass analyzer, activation), read from mzML/mzMLb CV terms, Thermo filter
strings, and timsTOF metadata. Groups never pool; ion-trap groups and groups
with too few PSMs stay uncorrected.

### Results

Five PXD028735 Orbitrap DDA files (+/-10 ppm precursor, +/-20 ppm fragment,
isotope errors -1..3) and one SILAC file, 1% FDR:

| Run | PSMs | Peptides | Proteins | Wall |
|---|---|---|---|---|
| PXD off | 283,901 | 59,004 | 7,082 | 162 s |
| PXD static | 283,988 | 58,937 | 7,091 | 121 s |
| PXD linear | 283,903 | 58,938 | 7,082 | 137 s |
| PXD auto | 283,797 | 58,909 | 7,072 | 119 s |
| SILAC off | 3,266 | 1,524 | 663 | 19 s |
| SILAC auto | 3,258 | 1,521 | 664 | 15 s |

Wall times vary by about 20% between runs on this shared machine, so only large
differences mean anything. On PXD, raw precursor medians drifted from 0.07 to
1.25 ppm across RT and from -0.2 to 1.3 ppm across m/z. Correction flattened
every bin to about 0 ppm and cut the held-out median absolute residual from 1.48
to 1.29 ppm. Most files chose a linear RT x m/z precursor model; one chose
smooth. Four files took a static fragment offset (2.32 -> 2.21 ppm). IDs did
not change (within 0.1%). The nonlinear step adds nothing measurable here.
Entrapment FDP at 1% peptide q-value was 1.09% without correction and 1.11%
with it, so there is no sign of overfitting.

Recommendation: keep `mass_recalibration` off by default. It removes real,
structured bias, but on well-calibrated Orbitrap data at +/-10 ppm the
candidate windows are already wide enough, so IDs do not change.

## Search tolerances from discovery residuals

`tolerance_mode: auto` reuses the discovery pass to narrow ppm windows: the
precursor per file, and fragments per file and acquisition group. Residuals
after the fit-set model are fitted as `pi * signal + (1 - pi) * Uniform(window)`
by EM, and the half-width is `max(h99, floor)` around the signal center, where
`h99` holds 99% of the signal. The floor is 3 ppm for precursors and 5 ppm for
fragments. The window is clipped to the configured one and applied only if it
holds 98% of the estimated held-out signal. Da and percent windows, mass-offset
searches, and windows beyond +/-100 ppm keep the configured window.

Two rules came first; both are kept below as the record of why the mixture is
needed.

The first rule tried was `4 * 1.4826 * MAD`, and it failed. Mass errors of
confident PSMs are heavy-tailed. For SILAC monoisotopic PSMs, the 50/90/99%
quantiles of |error| are 0.7/2.7/7.4 ppm, while four robust sigmas is about
4.2 ppm. Isotope-error PSMs, 38% of confident SILAC discovery PSMs, are wider
still. That rule narrowed the SILAC precursor window to about +/-4.3 ppm and
the fragment window to about +/-13 ppm, losing 9% of PSMs (3,266 -> 2,981). It
cost 1.6% of PXD PSMs (279,321). It also raised entrapment FDP from 1.09% to
1.23%: the windows cut true matches, not random ones.

The second rule was `max(1.2 * q99, floor)` around the fit-set residual
median, with `q99` the 99th percentile of absolute fit-set residuals. With
that percentile rule:

| Run | Fixed PSMs / peptides / proteins | Auto PSMs / peptides / proteins | Wall fixed -> auto |
|---|---|---|---|
| PXD +/-10/20 ppm | 283,901 / 59,004 / 7,082 | 283,972 / 58,999 / 7,082 | 162 -> 121 s |
| SILAC +/-10/20 ppm | 3,266 / 1,524 / 663 | 3,266 / 1,524 / 663 | 19 -> 16 s |
| PXD +/-50/50 ppm | 291,319 / 56,241 / 6,892 | 294,637 / 56,791 / 6,924 | 136 -> 179 s |
| SILAC +/-50/50 ppm | 3,539 / 1,571 / 687 | 3,560 / 1,580 / 687 | 18 -> 19 s |

Entrapment (one PXD file, combined FDP at 1% / 5% peptide q-value):
+/-10/20 ppm fixed 1.09% / 5.28%, auto 1.09% / 5.28%. At +/-50 ppm, fixed
1.10% / 5.47% and auto 1.05% / 5.35%.

The percentile rule is safe, but it narrows very little. The precursor window
was never narrowed. At +/-10 ppm the precursor q99 is about 9 ppm, and at
+/-50 ppm it is about 48 ppm: even among monoisotopic PSMs 1% have errors near
45 ppm at +/-50. So confident target PSMs fill whatever window is searched,
and a percentile of them cannot shrink it. Fragment q99 also tracks the window
(17 ppm at +/-20, 39 ppm at +/-50) because noise matches fill it. Narrowing at
+/-20 ppm was below 0.2 ppm. At +/-50 ppm it came to about +/-45 ppm, for +1.1%
PSMs, +1.0% peptides and slightly lower entrapment FDP. Fit-set and held-out
q99 agreed within 0.3 ppm everywhere.

### Signal plus background mixture

The current rule fits the residuals as a signal plus a uniform background
over the searched window. EM starts from the median and 1.4826 * MAD, keeps
the scale in [0.02 ppm, half the window] and `pi` in (0, 1), and stops after
500 iterations or a relative log-likelihood change below 1e-8. Inputs are
thinned to 40,000 points. Three signal forms are fitted on the fit set:

- a Gaussian;
- a Student-t with 4 degrees of freedom (EM with weights `5 / (4 + d^2)`);
- two Gaussians with a shared center.

A heavier form replaces a lighter one only if it converged and its held-out
log-likelihood is at least 0.005 nats per point higher. On PXD the t4 beat
the Gaussian by 0.001-0.04 nats per point, and two Gaussians beat the t4 by
0.002-0.016, so precursors mostly chose t4 or two Gaussians and fragments
always chose two Gaussians. Mass errors really are heavier-tailed than a
Gaussian.

A fit is not trusted, and the configured window is kept with the reason in
`skipped`, when EM did not converge, the scale is at a bound, `pi < 0.2`, or
`h99` exceeds 90% of the configured half-width (signal and background cannot
be told apart). No benchmark fit was rejected.

Isotope-error PSMs are fitted as their own subgroup, not pooled and not
dropped. They differ from monoisotopic PSMs: on PXD from +/-10 ppm their
signal fraction is about 0.5 (0.93-0.96 for monoisotopic), their center is usually
lower (by up to 0.8 ppm), and their Gaussian sigma is 2.1-3.0 ppm against 1.2-1.3
ppm for the narrow monoisotopic component. Pooling would force one scale on
both. Dropping them would size the window from monoisotopic PSMs alone,
though every isotope hypothesis is searched with the same window. So each
subgroup gets its own fit, and the window is the union of the two 99%
windows. With fewer than 200 isotope-error points they are left out, and an
untrusted isotope fit keeps the configured window (`isotope_error_*`).

Results (one run each; wall times vary by about 20% between runs):

| Run | PSMs / peptides / proteins at 1% FDR | Wall | Precursor window (ppm) | Fragment window (ppm) |
|---|---|---|---|---|
| PXD +/-10/20 fixed | 283,901 / 59,004 / 7,082 | 162 s | +/-10 | +/-20 |
| PXD +/-10/20 auto | 273,266 / 58,001 / 6,982 | 91 s | about -7 to +8 | about -14 to +12.5 |
| PXD +/-50 fixed | 291,319 / 56,241 / 6,892 | 136 s | +/-50 | +/-50 |
| PXD +/-50 auto | 281,418 / 58,867 / 7,052 | 97 s | about -8.5 to +9.5 | about -18 to +16.5 |
| SILAC +/-10/20 fixed | 3,266 / 1,524 / 663 | 19 s | +/-10 | +/-20 |
| SILAC +/-10/20 auto | 3,111 / 1,470 / 649 | 14 s | -6.10 to 6.26 | -12.56 to 11.82 |
| SILAC +/-50 fixed | 3,539 / 1,571 / 687 | 18 s | +/-50 | +/-50 |
| SILAC +/-50 auto | 3,272 / 1,527 / 667 | 15 s | -10.75 to 10.89 | -15.96 to 15.29 |

Entrapment (one PXD file, combined FDP at 1% / 5% peptide q-value): +/-10/20
fixed 1.09% / 5.28%, auto 1.20% / 5.30%. At +/-50, fixed 1.10% / 5.47%, auto
1.11% / 5.28%.

SILAC ablations from +/-10/20: narrowing only the precursor to the auto
window gives 3,136 PSMs, and only the fragment gives 3,244. Fixing both at
the windows auto chose from +/-50 gives 3,271, the same as auto (3,272).

What this shows:

- From +/-50 ppm, auto converges to roughly +/-9 ppm (PXD) or +/-11 ppm
  (SILAC) for precursors and +/-15-18 ppm for fragments. That is close to the
  well-chosen fixed +/-10/20. IDs match it within 1%: PXD peptides -0.2% and
  proteins -0.4%, SILAC PSMs +0.2%. Entrapment FDP is unchanged. Against
  +/-50 fixed it gains 4.7% PXD peptides but loses 3.4% PSMs (PXD) and 7.5%
  (SILAC).
- From +/-10/20 ppm, auto loses 3.7% PSMs and 1.7% peptides on PXD and 4.7%
  PSMs on SILAC, and entrapment FDP rises from 1.09% to 1.20%. Most of the
  SILAC loss is the precursor window. A 99% signal window drops more than 1%
  of matches: the searched window is only about 1.3-1.5 times the signal
  width, so the uniform term absorbs part of the true tail, and `h99` shrinks
  with the searched window (about 7.5 ppm from +/-10 against 9 ppm from +/-50).
- The fitted window is used faithfully. Fixed at auto's windows gives the
  same IDs as auto.

Recommendation: keep `tolerance_mode` fixed by default. The mixture fixes the
percentile rule's failure (it now narrows, and from a too-wide window it
lands near a sensible fixed window). But it does not beat a well-chosen fixed
tolerance, and from a sensible window it costs IDs. It is useful only as a
guard against a window set much too wide.
