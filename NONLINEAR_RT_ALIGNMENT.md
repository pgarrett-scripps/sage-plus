# Nonlinear retention-time alignment exploration

Enable the prototype explicitly in the search JSON:

```json
"retention_time_alignment": "nonlinear"
```

This explicitly enables nonlinear alignment of observed retention times, independently
of `predict_rt`. Retention-time prediction and LFQ use alignment internally; if no
method is specified for those workflows, they retain Sage's existing linear alignment.

Sage currently maps each run to a global consensus with one ordinary least-squares
line. That is fast and easy to extrapolate, but a few incorrect peptide landmarks can
move the fit and one line cannot represent local changes in a chromatography gradient.

## Prototype in this branch

The prototype uses shared, high-confidence peptide IDs as landmarks and fits each run
in two stages:

1. Normalize RT by the run's maximum observed RT and use the across-run median RT as
   each peptide's consensus target. Peptides observed in only one run are excluded.
2. Fit many deterministic two-point affine hypotheses and retain the one with the
   smallest median residual. An adaptive median-absolute-deviation cutoff removes gross
   outliers without treating ordinary curve shape as an outlier.
3. Divide the remaining landmarks into RT quantile bins and take median coordinates.
4. Apply weighted isotonic regression to the bin targets, then use piecewise-linear
   interpolation between the monotone knots.
5. Fall back to a robust affine alignment when there are fewer than 16 landmarks or
   when they cover less than 25% of the normalized gradient.

The same `Alignment::transform` method is used for PSM RTs and LFQ MS1 scan times.

## Options considered

- **RANSAC affine only:** robust to bad IDs and large offsets, but does not solve local
  nonlinear drift. It is useful as the first stage and as a sparse-data fallback.
- **LOESS:** captures smooth local changes, but unconstrained LOESS can reverse local
  elution order and behaves poorly near sparse boundaries.
- **Cubic or smoothing splines:** smoother derivatives than linear interpolation, but
  ordinary splines can overshoot. A monotone cubic Hermite spline is a promising later
  refinement if piecewise-linear corners matter in real datasets.
- **Dynamic time warping:** useful when dense comparable chromatograms are available,
  but Sage's current inputs naturally provide sparse peptide landmarks. DTW also needs
  explicit regularization to avoid implausible warps.
- **Monotone piecewise linear:** simple, dependency-free, predictable outside local
  perturbations, and directly compatible with LFQ range lookup. This is the prototype.

## Evaluation before merging

The synthetic unit tests cover outliers, large shifts, monotonicity, nonlinear error,
and sparse fallback. A dataset-level comparison should additionally measure:

- held-out shared-peptide median absolute RT error per run;
- error by RT decile, especially gradient boundaries;
- the number of LFQ features recovered at a fixed FDR;
- coefficient of variation for confidently quantified peptides;
- knot count, rejected-landmark fraction, and runtime/memory; and
- behavior for disconnected run groups with few or no shared peptides.

The current global consensus is computed once. If experiments contain several weakly
connected batches, a useful next experiment is graph-based reference selection followed
by one or two consensus/refit iterations.

## Related work

- Prince and Marcotte, *ChromA: signal-based retention time alignment for
  chromatography-mass spectrometry data* (2009):
  <https://pmc.ncbi.nlm.nih.gov/articles/PMC2722998/>
- Kirchner et al., *amsrpm: Robust Point Matching for Retention Time Alignment of
  LC/MS Data with R* (2007): <https://doi.org/10.18637/jss.v018.i04>
- Fischer et al., *Retention Time Alignment Algorithms for LC/MS Data Must Consider
  Non-Linear Shifts* (2009): <https://doi.org/10.1093/bioinformatics/btp052>
- Christin et al., *Time Alignment Algorithms Based on Selected Mass Traces for
  Complex LC-MS Data* (2010): <https://doi.org/10.1021/pr9010124>
