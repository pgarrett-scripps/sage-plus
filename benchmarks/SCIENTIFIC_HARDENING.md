# Scientific hardening follow-up

This development experiment preserves the archived 20260914 pilot. Outputs and compact builds live in `/data/sage-plus-scientific/hardening-20260914`. This is a development rerun of previously inspected data, not a reserved validation set.

The frozen pilot source modules are unchanged. The follow-up uses their existing explicit binary override and writes distinct jobs with input and output hashes.

## PTM confidence repair

The first synthetic HCD run has 3,141 target PSMs and one decoy PSM. A Gaussian KDE fitted to the singleton decoy class has zero bandwidth. Its undefined probabilities propagated into peptide confidence, leaving every peptide q-value at one.

When either competition class has fewer than two observations, constant scores, or nonfinite fitted posterior values, the candidate falls back to cumulative target-decoy counts with a +1 decoy correction. Complete score ties receive identical q-values. A warning makes the fallback explicit. This repairs numerical behavior. Sparse decoys and the restricted synthetic database still do not establish production calibration.

The comparison endpoint remains joint PSM, peptide and localization q-values at 1%. Synthesis-site consistency remains a diagnostic until acquisition files are independently mapped to synthesized libraries. Site-library candidate generation does not constrain the localizer's residue-compatible positions.

## LFQ evidence and experimental transfer estimates

The existing precursor q-value is retained and explicitly identified as applying across files. It is not a per-file transfer estimate. Exact score ties in count-based precursor confidence now share a common threshold.

LFQ schemas 3 and 4 add strict direct MS2 evidence, requiring both PSM and peptide q-values at the configured LFQ peptide threshold. Existing `ms2_confirmed` semantics are retained for compatibility.

Each positive file signal exposes isotope spectral angle at the selected shared apex, cosine similarity to the warped reference trace, and the local warp offset in grid bins. The experimental ranking score is:

`isotope_angle^3 * trace_cosine * max(0, 1 - abs(warp_bins) / grid_bins)`

This score is fixed before inspecting candidate outcomes. It does not change integration or remove rows. Missing signals have null diagnostics. A transfer candidate lacks strict direct target MS2 evidence in the acquisition file. The same eligibility condition applies to the corresponding shifted decoy, even if no target peak survives integration. MBR-off signals are not transfer candidates.

The separate analyzer estimates experimental q-values from cumulative shifted-decoy counts with a +1 correction within each recipient file, using only transfer candidates and keeping exact score ties together. No precursor q-value filter is applied before estimation. Results at 1% and 5% are reported alongside the legacy precursor-filtered population.

The pure-human check excludes I/L-equivalent sequences present in multiple species. Foreign yeast or E. coli sequences provide a diagnostic lower bound on incorrect transfers under the sample-purity and reference-completeness assumptions. Human-to-human mistakes remain invisible. Empty accepted sets are not evidence of calibration.

Transfer scores inherit the existing shared-apex selection, local warping and shifted-decoy construction. Their null model may fail. If foreign assignments contradict the proposed threshold, the estimate stays experimental and no production FDR filter or calibration claim is introduced. This experiment must precede any reserved-data validation.

LFQ schemas 5 and 6 add `extraction_q_value` for every target precursor/file row. Each target's shifted decoy (+11.06 m/z, earlier retention time) is evaluated at the target's selected apex, window and reference trace, with its own per-file warp search, and these paired decoy rows compete with all target rows by `file_score` in one pool. Evaluating decoy precursors at their own apex made decoy reference-file rows as strong as target rows (209 rows passed 1%); using the target's warps without a search made decoys too weak (239,874 rows passed, and foreign rows were 29.7% of passing Human_01 transfers). On the PXD028735 five-file HYE set with MBR on, 172,866 of 293,684 target rows pass 1% (139,101 MS2-backed, 33,765 transfers). In the pure-human file, 2,726 transfers pass and 617 of them (22.8% of species-assigned rows) are yeast or E. coli, a lower bound on false transfers. Ranking transfers only against paired decoy transfers gives 18,529 passing transfers and 19.3% foreign. A decoy offset that keeps the peptide mass defect (+11.0055 m/z) still gave 19.0% foreign while collapsing precursor-level discovery, so the decoy offset is not the limiting factor: MS1 extraction evidence alone does not separate these foreign transfers from real ones. The value is reported for all rows and is not claimed as a calibrated transfer FDR.

Transfer-specific confidence populations are also distinguished in the [IonQuant study](https://pmc.ncbi.nlm.nih.gov/articles/PMC8131922/). The candidate score and simple count estimator here do not reproduce IonQuant's mixture model.

## Follow-up after inspecting PTM acceptance

The first repaired run revealed confident site calls with exactly tied best target arrangements. In HCD 1, 103 of the 118 inconsistent accepted site events had that ambiguity. In HCD 2, 14 of 15 did. This diagnosis was made after inspecting the development reruns.

The second candidate retains those competitions when estimating the dataset curve, then assigns localization q-value 1 to unresolved target ties. A single possible arrangement remains eligible despite a zero delta score. The change prevents arbitrary tie breaking from becoming a confident site call. It does not establish calibration of the remaining subset. The second candidate, source patch and reruns are recorded separately.
