#import "stats.typ": lit, s
#import "assets.typ": fig, tbl

This Supporting Information records the inputs, software identities, and
additional paired results for the Sage and Sage Plus versions named in the main
text.

= Frozen software and evidence <sec:si-provenance>

Sage `v0.15.0-beta.2` was obtained from #link(
  "https://github.com/lazear/sage",
)[the Sage repository], source commit
`df9219951cc9a54cf4cd55d76541af24b687bd3d`. Sage Plus `v0.1.0-beta.3` was
obtained from #link("https://github.com/pgarrett-scripps/sage-plus")[the Sage
  Plus repository], source commit `9bfc8acbc3002d83aa3f27af5f848376b5c5c9be`.
FDRBench `1.1.1` used source commit `3b619a9acf60d7292fb651a00da55f58cb67fb79`.

The local benchmark archive is
`/data/sage-plus-scientific/20260914/artifacts/sage-plus-scientific-pilot-20260914.tar.gz`.
Executable hashes, configurations, and input identities are recorded with the
repository summary at `benchmarks/scientific-results/20260914/`. Large spectra
remain at the paths recorded in the archive manifests. The worker-scaling and
LFQ extension records are retained at
`/data/sage-plus-scientific/report-extension-20260915`. These locations are
local evidence records, not an external deposition.

= Development source map <sec:si-development>

The development inventory uses the Sage and Sage Plus source commits recorded in
@sec:si-provenance. The file `analysis/data/development-changes.json` records
the table text, source paths, and content digests of the inspected release
files. The following paths are relative to those repository snapshots.

Input readers are in `crates/sage-cloudpath/src/`. Sage Plus adds `mzmlb.rs` and
`thermoraw.rs`. Peptide and cleavage input changes are in
`crates/sage/src/database.rs` and `cleavage.rs`. Canonical result schemas and
writers are in `crates/sage-cloudpath/src/parquet/`, including `results.rs`,
`lfq.rs`, and `sites.rs`.

Peptide storage is defined in `crates/sage/src/peptide.rs` and `sequence.rs`.
The fragment index is in `database.rs`, and spectrum storage and isotope
processing are in `spectrum.rs`. Resource estimation and limits are in
`crates/sage-cli/src/memory.rs` and the runner modules. API events and summaries
are defined by `api.rs`, `events.rs`, and `runner/artifacts.rs`. Persistent
worker orchestration is in `crates/sage-mcp/src/lib.rs`.

Scientific extensions are implemented in `crates/sage/src/mass_calibration.rs`,
`scoring.rs`, and the retention and mobility modules under `ml/`. Modification
semantics and neutral losses are defined in `modification.rs` and
`ion_series.rs`. Quantification, localization, and empirical library export are
in `lfq.rs`, `ptm.rs`, `ambiguity.rs`, and `spectral_library.rs`. The released
`CHANGELOG.md`, `DOCS.md`, and `UPSTREAM.md` supply the corresponding interface
and maintenance descriptions.

This map establishes the provenance of the implementation description. Benchmark
provenance remains tied to the frozen executable hashes. A documented capability
and a measured endpoint are distinguished in @sec:development.

= Public input order and reference amendment <sec:si-inputs>

The HEK file indices are defined by the following deposited mzXML basenames.
Their derived MGF files preserve those names.

- `PXD001468-0`: `b1906_293T_proteinID_01A_QE3_122212`.
- `PXD001468-1`: `b1937_293T_proteinID_01B_QE3_122212`.

The mixture indices follow the order in the frozen comparison plan. All use the
`LFQ_Orbitrap_DDA_` prefix and the first injection.

- `PXD028735-0`: `Condition_A_Sample_Alpha_01`.
- `PXD028735-1`: `Condition_B_Sample_Alpha_01`.
- `PXD028735-2`: `Condition_A_Sample_Beta_01`.
- `PXD028735-3`: `Condition_B_Sample_Beta_01`.

Repeat timing used the index-zero file of each study. The acquisition receipts
provide complete source filenames and URLs. Public identification pairs are
indexed in the summary, so the association can be checked without interpreting
figure order.

The original mixture reference contained undefined residues in these _E. coli_
entries:

```text
P33369  P45766  P75901  P39901  P76000  P37003  P58095
```

The amended `hye-irt-defined.fasta` excluded these entries in full for both
engines. The original reference, comparison plan, and failed Sage Plus outputs
were retained alongside the amended `public-comparison-v2` results.

The community sample annotation named a different yeast species from the primary
methods. The reagent identified in those methods was Promega V7461. Its
manufacturer documentation identified _S. cerevisiae_, which determined the
pilot's reference choice. #link(
  "https://www.promega.com/-/media/files/resources/protocols/technical-manuals/101/ms-compatible-yeast-and-human-protein-extracts-protocol.pdf",
)[Manufacturer technical manual TM410] records that reagent identity. Reference
completeness and the precise _E. coli_ strain remain limitations rather than
assumptions silently resolved by the amendment.

= Additional paired resource measurements <sec:si-resource>

@tbl:si-timing retains the exact summaries behind the public timing figure.

@fig:si-enlarged gives the absolute resource measurements for the entrapment
workloads compared on a relative scale in @fig:workloads.

#figure(
  fig("fig.supplement-entrapment-resources", width: 100%),
  caption: [Resource use with entrapment-expanded references. Large markers and
    labels show median wall time (A) and peak resident memory (B). Faint points
    show two selected files with three construction seeds per study, offset
    vertically for visibility. These are different file-by-seed searches, not
    repeated timing trials of a fixed input. The mixture reference includes the
    documented exclusions.],
) <fig:si-enlarged>

The local HEK conventional and common-modification results appear in
@fig:si-local. These searches used the same frozen release pair as the main
text. Acquisition could overlap these runs, so their wall times are not part of
the controlled public timing endpoint. The common-modification workload adds
methionine oxidation and peptide N-terminal acetylation under the shared search
configuration. The local biological input has unverified redistribution
provenance and is not offered as a public replication dataset.

#figure(
  fig("fig.supplement-local", width: 100%),
  caption: [Contextual local HEK searches with standard and common-modification
    settings. Panels compare wall time (A), peak resident memory (B), and
    accepted target PSMs (C). Large markers and labels show medians of three
    measured trials after one warmup. Faint points show individual trials,
    offset vertically. PSMs pass the reported one-percent spectrum q-value
    threshold. Possible acquisition overlap limits timing interpretation.],
) <fig:si-local>

The broad PTM workload did not complete with either engine under the selected
limits. Plus applied its additional modified-database guard and Sage encountered
allocation failure under the shared address-space ceiling. Failed warmups and
attempts are retained. Their elapsed time is not summarized as the time required
to finish the search.

= Threshold-dependent entrapment results <sec:si-thresholds>

@fig:si-entrapment extends the primary-threshold comparison in @fig:entrapment
with calibration curves and individual file-by-seed differences.

#figure(
  fig("fig.scientific-entrapment-detail", width: 100%),
  caption: [Entrapment calibration for HEK (A) and mixture (B), with mean
    conservative paired FDP profiles and faint lines for individual files and
    construction seeds. Panels C and D show file-by-seed differences at the
    nominal one-percent peptide threshold and study means with conditional
    bootstrap intervals. Seed labels follow the order in this supplement.
    Nominal peptide q-values are distinct from independently estimated FDP.],
) <fig:si-entrapment>

@fig:si-threshold-yield shows target yield at every prespecified nominal
threshold, averaging the same file-by-seed combinations as the FDP curves.

#figure(
  fig("fig.supplement-threshold-yield", width: 100%),
  caption: [Target peptide yield across nominal peptide q-value thresholds for
    HEK (A) and mixture (B). Each point averages two files and three shared
    construction seeds. Lines join prespecified thresholds on a logarithmic
    horizontal scale. Panel-specific vertical scales start at zero. These are
    nominal-threshold yields. The independently estimated FDP curves appear in
    @fig:si-entrapment, so a higher count alone does not demonstrate greater
    sensitivity at matched error.],
) <fig:si-threshold-yield>

Construction seeds were `20260914`, `20260915`, and `20260916`. The bootstrap
used seed `20260914`, with bounds at the #lit("2.5")th and #lit("97.5")th
percentiles. The main Methods describe the paired resampling and audit
conventions. Search failures were retained rather than converted to empty
discovery sets.

== Peptide yield at a common estimated FDP ceiling <sec:si-matched-fdp>

The retained peptide-level outputs permitted evaluation of every distinct
observed peptide q-value threshold, without interpolation between the
prespecified operating points. Each threshold retained its complete q-value tie
group. For each engine, file, and construction seed, we selected the threshold
giving the largest target peptide count at or below a conservative paired FDP of
#lit("1") percent. Ties in target count were resolved by taking the largest
q-value threshold. Exact target-partner score ties used the same conservative
bound as the primary analysis.

The selected FDP estimates ranged from #s("matched.fdp.min") to #s(
  "matched.fdp.max",
) percent. @fig:si-matched-fdp reports the achieved ranges and target peptide
counts. Counts exclude entrapment peptides, while the paired FDP denominator
includes both target and entrapment discoveries. The main-text yield changes are
means of paired file-by-seed percentage differences, rather than a comparison of
pooled counts.

#figure(
  fig("fig.supplement-matched-fdp", width: 100%),
  caption: [Exploratory peptide yield at a common one-percent estimated paired
    FDP ceiling. Panel A shows mean target peptide counts across two files and
    three construction seeds. Panels B and C show the achieved FDP and selected
    nominal peptide q-value for each file-by-seed search, offset vertically.
    Horizontal bars span observed ranges, not confidence intervals. The dotted
    line marks the FDP ceiling. Thresholds retain complete observed q-value
    steps. FDP estimation and threshold selection reuse the same entrapment
    observations, so this is not held-out error control.],
) <fig:si-matched-fdp>

Discrete thresholds need not attain the same FDP exactly. Reusing the entrapment
observations to select thresholds can favor apparent yield, so the comparison
does not establish greater sensitivity or equivalent true error rates.

= Public identification counts and worker scaling

@tbl:si-overlap records the exact accepted PSM agreement counts plotted in
@fig:overlap. Shared identities contribute once to the union.

@fig:si-identification shows threshold-dependent changes in accepted PSMs and
peptidoforms without duplicating their absolute counts.

#figure(
  fig("fig.report-identification", width: 100%),
  caption: [Mean file-level percentage changes from Sage to Sage Plus across
    nominal spectrum and peptide q-value thresholds for HEK (A) and mixture (B).
    Shading spans the observed file range, not a confidence interval. PSMs use
    spectrum q-values and peptidoforms use peptide q-values. The same files are
    used at every threshold.],
) <fig:si-identification>

@fig:si-public-counts gives the absolute accepted PSM and peptidoform counts
behind the primary-threshold comparison, together with retained decoy PSMs.

#figure(
  fig("fig.supplement-public-counts", width: 100%),
  caption: [Accepted identifications for each public file at the one-percent
    threshold of the indicated confidence level. Panels show target PSMs (A),
    target peptidoforms (B), and decoy PSMs (C). Engine markers are offset
    vertically for visibility. PSMs use spectrum q-values and peptidoforms use
    peptide q-values. Each count axis starts at zero. File labels follow
    @sec:si-inputs. Decoy PSM counts describe retained search output and are not
    independent entrapment estimates.],
) <fig:si-public-counts>

@tbl:si-scaling reports every engine and worker-count summary from the
extension. Each range spans three measured repeats following one warmup. The PSM
range is shown alongside time so that a repeat with changed yield cannot be
hidden by a single aggregate timing value. Assignment differences are analyzed
separately in the main text.

= Quantitative comparison and control denominators

@tbl:si-lfq retains the exact engine-specific values summarized in
@fig:lfq-endpoints, alongside the ratio-pair denominators.

@fig:si-lfq-shared restricts the accuracy comparison to peptide and preparation
pairs accepted by both releases. Combined-charge precursors use the same empty
charge key in both formats. Numeric modification spelling is normalized without
changing mass values or modification positions. The shared set is intersected at
the ratio-pair level, so both conditions have positive accepted intensity in
both engines.

#figure(
  fig("fig.supplement-lfq-shared", width: 100%),
  caption: [Ratio accuracy on identical accepted peptide-preparation pairs.
    Panel A compares each engine's median absolute error against the expected
    species log-base-two B/A ratio. Panel B shows the median absolute difference
    between the engines' observed ratios on those same pairs. Labels give values
    and shared pair counts. Engine markers in A are offset vertically. These
    descriptive medians have no uncertainty intervals. The shared sets
    complement the engine-specific accepted sets in @fig:lfq-endpoints.],
) <fig:si-lfq-shared>

The control counts in @tbl:si-control use the species and direct-MS2 rules
defined in the main Methods. Retention-time standards and unmapped proteins are
excluded from the denominator. Numerators and denominators both require positive
finite accepted intensities, and direct evidence is drawn from rank-one PSM
results.

= Synthetic phosphorylation challenge <sec:si-ptm>

== Search design

The phosphorylation challenge evaluated a PTM search setting. It used
higher-energy collisional dissociation (HCD) acquisitions from PRIDE PXD000138
and synthesis-defined peptide sequences @marx2013. The restricted database
contained the synthetic sequences, with one variable phosphorylation permitted
on serine, threonine, or tyrosine. Peptides were searched as intact library
entries, with fixed cysteine carbamidomethylation and a shared target-decoy
configuration. This restricted search does not reproduce the original
full-background database experiment.

Both releases were evaluated at the spectrum-only threshold and at the joint
spectrum and peptide thresholds. The primary localization endpoint additionally
required an accepted site with a reported localization q-value. Sage Plus
exports site-level results, but the Sage output used here does not provide an
identical site-confidence field. We therefore compared the shared upstream
acceptance stages first.

The secondary Sage Plus synthesis-consistency diagnostic pools the available
library truth because the mapping between acquisition files and individual
synthetic libraries was not independently verified. It combines identification
and localization discrepancies and cannot estimate arrangement-level false
localization rate.

The synthesis reference and truth were prepared before the original searches,
and all-site configurations are retained for both releases. Comparative accuracy
claims exclude oracle site-prior runs using synthesis truth and development
repairs made after the released snapshot.

== Peptide acceptance

Evaluation of the new PTM functionality identified a limit of the released
workflow. The restricted phosphorylation challenge produced the same
spectrum-level accepted counts in both releases. Each engine accepted #s(
  "report.ptm.1.upstream.spectrum_accepted",
) target PSMs in the first HCD file and #s(
  "report.ptm.2.upstream.spectrum_accepted",
) in the second. Requiring the peptide threshold reduced the accepted sets to
#s("report.ptm.1.upstream.joint_accepted") for both releases and both files.
@fig:ptm shows this loss at the confidence-filtering stage.

All retained rank-one matches had peptide q-values equal to one. The shared
peptide acceptance stage therefore precluded site-accuracy evaluation before the
Sage Plus localization cutoff was considered. Omitting the peptide filter does
not produce a validated phosphosite set. Sparse decoy evidence further limits
interpretation of this restricted synthetic search.

#figure(
  fig("fig.report-ptm", width: 100%),
  caption: [Acceptance stages in the matched synthetic phosphorylation search.
    Panels A and B correspond to the two HCD input files. Bars compare target
    PSMs passing the spectrum threshold alone with those passing both spectrum
    and peptide thresholds at one percent. The joint accepted sets are empty in
    both engines. The figure reports successful search outputs with failed
    primary acceptance, not search failures or a zero localization-error
    estimate.],
) <fig:ptm>

== Site diagnostic

@fig:si-ptm-diagnostic records the released Sage Plus site output after spectrum
and localization filtering, showing the events withheld by the primary peptide
confidence requirement.

#figure(
  fig("fig.supplement-ptm-diagnostic", width: 100%),
  caption: [Secondary Sage Plus synthesis-consistency diagnostic. Bars partition
    site events into synthesis-consistent and inconsistent events, with counts
    and inconsistent fractions labeled. Events pass spectrum and localization
    q-values at one percent but omit the peptide filter. Restoring that filter
    accepts no sites in either file. Consistency uses pooled unambiguous
    synthesis truth with file-to-library mapping unaudited. It is not an
    arrangement-level localization error estimate.],
) <fig:si-ptm-diagnostic>

#pagebreak()
= Tables <sec:si-tables>

These numerical reference tables retain exact counts, denominators, and timing
ranges that complement the comparative figures. Resource summaries correspond to
@fig:timing and @fig:scaling. Identification counts correspond to @fig:overlap,
and quantitative endpoints correspond to @fig:lfq-endpoints and @fig:control.

#figure(
  tbl("tbl.scientific-timing"),
  caption: [Public timing medians for the specified Sage and Sage Plus releases.
    Parentheses contain observed minimum and maximum wall times, not confidence
    intervals. Peak memory is maximum resident set size in MiB. Target PSMs pass
    the reported one-percent spectrum q-value threshold. Warmups are excluded.],
) <tbl:si-timing>

#figure(
  tbl("tbl.scientific-overlap"),
  caption: [Accepted target PSM agreement at each engine's reported one-percent
    spectrum q-value threshold. File indices identify the frozen input order in
    @sec:si-inputs. Mixture rows use the amended reference for both engines.
    Jaccard is the shared count divided by the union of accepted PSM
    identities.],
) <tbl:si-overlap>

#figure(
  tbl("tbl.report-scaling"),
  caption: [Matched worker-scaling summaries from the report extension. Seconds
    and peak resident MiB are medians. Time and PSM ranges span measured
    repeats. A single PSM count indicates identical yield in every repeat.
    Warmups remain in the execution records but are excluded from these
    summaries.],
) <tbl:si-scaling>

#figure(
  tbl("tbl.report-lfq"),
  caption: [Quantification endpoints for the matched released engines at the
    primary precursor threshold. Ratio pairs combine a peptide and preparation.
    Bias and absolute error are on the log-base-two B/A scale. CV describes
    Alpha versus Beta preparation variability within condition. Missingness uses
    the engine-specific observed union, not the theoretical proteome.],
) <tbl:si-lfq>

#figure(
  tbl("tbl.report-control"),
  caption: [Exact human-only control denominators at the primary LFQ threshold.
    Foreign means a yeast or E. coli assignment without a human-reference
    equivalent. Direct MS2 requires matching peptide and file with both spectrum
    and peptide q-values at most one percent. No direct MS2 is a subset of the
    complete control set, not an additional independent sample.],
) <tbl:si-control>
