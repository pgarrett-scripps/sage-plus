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
default. Both searches used five PXD028735 Orbitrap LFQ files (HYE mixture) and the
committed configuration
[`benchmarks/configs/prefilter-lfq-5file.json`](../../benchmarks/configs/prefilter-lfq-5file.json):
`hye-irt-defined.fasta`, M oxidation and N-terminal acetylation as variable modifications
(at most one of each and two per peptide), 10 ppm precursor and 20 ppm fragment tolerance,
LFQ enabled with match-between-runs on (`quant.lfq_settings.mbr = true`), 8 threads. The
figures below are that configuration; a search without the two variable modifications finds
about 11% fewer PSMs and is not comparable to them. The precursor counts below depend on that
setting: match-between-runs is what fills a precursor's missing files, so with it off the
same searches quantify about 41% of precursors in all five files instead of 98%.
Residuals are the absolute difference between each peptide's aligned RT and its
across-run median, for target peptides at 1% peptide FDR seen in at least two runs,
converted to minutes.

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

### The same configuration on the current code

The table above is the Beta 9 measurement that made nonlinear the default, and it is kept as
the record of that decision rather than restated. Rerunning the same configuration on the
current code gives:

| Metric | Linear | Nonlinear |
| --- | ---: | ---: |
| PSMs at 1% FDR | 312,867 | 312,744 |
| Peptides at 1% FDR | 62,204 | 62,218 |
| Protein groups at 1% FDR | 6,970 | 6,965 |
| LFQ precursors at 1% q-value | 20,714 | 31,376 |
| LFQ precursors quantified in all 5 files | 19,859 | 30,812 |

Identifications are within 0.2% of the Beta 9 figures and the nonlinear LFQ column within
0.2% as well. The linear column is 13.8% lower, for a reason that has nothing to do with
either alignment method: upstream Sage commit `75a5300` changed how LFQ breaks a tie between
candidate warp offsets. A trace whose reference carries no signal used to be pushed out of
its own integration window by the largest offset in the search and die; it is now kept at
zero shift, integrates, and competes. Such traces are almost all shifted decoys, so target
precursors and their peak scores do not move while the surviving decoy population grows,
which pushes the 1% precursor q-value cut up. Under linear alignment the target peak scores
overlap the decoy tail, so the cut has to travel far enough to cost 3,159 precursors; under
nonlinear alignment it travels about a fifth as far and costs 245.

Both columns are therefore more conservative than they were at Beta 9, and the gain from
nonlinear alignment measured on the current code (+51%) is larger than the +31% measured at
Beta 9. The Beta 9 figure is the one quoted in the README and FEATURES.md, because the
difference between the two is a correction to decoy counting rather than an improvement in
alignment.

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
