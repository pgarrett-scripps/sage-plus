# Nonlinear retention-time alignment

Nonlinear alignment is the default from Beta 9. Alignment runs whenever `predict_rt` is
true (the default), LFQ is enabled, or `retention_time_alignment` is set. To restore the
ordinary least-squares alignment used before Beta 9, set:

```json
"retention_time_alignment": "linear"
```

With `"linear"`, results are byte-identical to the previous default. `run-summary.json` records the method
used under `models.retention_time_alignment`.

Upstream Sage maps each run to a global consensus with one ordinary least-squares
line. That is fast and easy to extrapolate, but a few incorrect peptide landmarks can
move the fit and one line cannot represent local changes in a chromatography gradient.

## Method

The method uses shared, high-confidence peptide IDs as landmarks and fits each run
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
   when they cover less than 25% of the normalized gradient. A single-file search has no
   shared landmarks and is left unwarped, which gives the same aligned times as `linear`.

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
  perturbations, and directly compatible with LFQ range lookup. This is the chosen method.

## Evaluation

Synthetic unit tests cover outliers, large shifts, monotonicity, nonlinear error, sparse
fallback, and single-run behavior. The dataset comparison below made nonlinear the
default. Both searches used five PXD028735 Orbitrap LFQ files (HYE mixture,
`hye-irt-defined.fasta`, LFQ enabled, 8 threads). Residuals are the absolute difference
between each peptide's aligned RT and its across-run median, for target peptides at 1%
peptide FDR seen in at least two runs, converted to minutes.

| Metric | Linear | Nonlinear |
| --- | ---: | ---: |
| PSMs / peptides / proteins at 1% FDR | 312,533 / 62,817 / 7,051 | 312,487 / 62,795 / 7,045 |
| Protein groups at 1% FDR | 7,177 | 7,171 |
| LFQ precursors at 1% q-value | 24,038 | 31,418 |
| LFQ precursors quantified in all 5 files | 23,055 | 30,828 |
| Human log2(A/B) MAD (expected 0) | 0.213 | 0.191 |
| Yeast log2(A/B) median / MAD | 1.111 / 0.208 | 1.114 / 0.198 |
| E. coli log2(A/B) median / MAD | -2.225 / 0.388 | -2.215 / 0.397 |
| Alignment residual median / 90th percentile (min) | 0.128 / 0.913 | 0.077 / 0.391 |

Each run received 24 knots. LFQ decoys at 1% q-value stayed near 1% of targets (239 and
313). Alignment and LFQ take the same time with either method; end-to-end runtime
differences were within run-to-run variation on a shared machine. A single-file HEK SILAC
search produced byte-identical results with both methods.

Still worth measuring on other datasets:

- held-out shared-peptide RT error, by RT decile, especially at gradient boundaries;
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
