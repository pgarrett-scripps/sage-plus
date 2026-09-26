#import "stats.typ": lit, s
#import "assets.typ": fig, tbl

This Supporting Information records the inputs, software identities, and
additional paired results for the Sage and Sage Plus versions named in the main
text.

= Frozen software and evidence <sec:si-provenance>

Sage `v0.15.0-beta.2` was obtained from #link(
  "https://github.com/lazear/sage",
)[the Sage repository], source commit
`df9219951cc9a54cf4cd55d76541af24b687bd3d`. Sage Plus `v0.1.0-beta.6` was
obtained from #link("https://github.com/pgarrett-scripps/sage-plus")[the Sage
  Plus repository], source commit `3e30135fb8786ec8a12c1f62e0ff9300e57f9567`.
FDRBench `1.1.1` used source commit `3b619a9acf60d7292fb651a00da55f58cb67fb79`.

The comparison was refreshed using the published Sage Plus Linux executable,
verified against the release SHA256SUMS manifest. The evidence directory is
`runs/paper-refresh-20260920/`, relative to the Sage Plus working repository.
Its `evidence/` and `extension/` directories retain commands, configurations,
executable and input hashes, run outcomes, and analytical outputs. Large spectra
and unchanged reference inputs remain at their recorded acquisition paths. These
are local reproducibility records, not an external deposition. All comparative
results in this manuscript use the release pair specified above.

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
@tbl:execution includes all comparative attempts, including warmups and
failures.

#figure(
  tbl("tbl.execution"),
  caption: [Execution outcomes for the two compared engines. Attempts include
    warmups and any explicitly recorded repeat after a failure. Completed
    searches passed output extraction. Failed attempts are retained separately
    and excluded from successful-run resource and identification summaries.],
) <tbl:execution>

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
text. Their wall times are reported separately from the public timing endpoint.
The common-modification workload adds methionine oxidation and peptide
N-terminal acetylation under the shared search configuration. The local
biological input has unverified redistribution provenance and is not offered as
a public replication dataset.

#figure(
  fig("fig.supplement-local", width: 100%),
  caption: [Contextual local HEK searches with standard and common-modification
    settings. Panels compare wall time (A), peak resident memory (B), and
    accepted target PSMs (C). Large markers and labels show medians of completed
    measured trials after one warmup. Faint points show individual trials,
    offset vertically. PSMs pass the reported one-percent spectrum q-value
    threshold. Unverified biological input provenance limits interpretation.],
) <fig:si-local>

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

Both engines accepted #s("report.ptm.1.upstream.spectrum_accepted") target PSMs
in the first HCD file and #s("report.ptm.2.upstream.spectrum_accepted") in the
second at the spectrum threshold. Adding the peptide threshold retained #s(
  "report.ptm.1.plus.joint_accepted",
) and #s("report.ptm.2.plus.joint_accepted") PSMs in Sage Plus, respectively.
Upstream Sage retained #s("report.ptm.1.upstream.joint_accepted") in each file.
@fig:ptm distinguishes spectrum-level acceptance from the joint confidence rule.

The count-based fallback in Sage Plus permits peptide acceptance when density
modeling is underdetermined. This behavior does not by itself establish error
calibration. Sparse decoy evidence in this restricted synthetic database limits
interpretation, and the absence of accepted upstream peptides prevents a matched
site-accuracy comparison.

#figure(
  fig("fig.report-ptm", width: 100%),
  caption: [Acceptance stages in the matched synthetic phosphorylation search.
    Panels A and B show the two HCD inputs. Bars count target PSMs passing the
    spectrum threshold alone or both spectrum and peptide thresholds at one
    percent. Sage Plus retains a joint accepted set in each file. Sage retains
    none. These are confidence-filtering outcomes from completed searches.],
) <fig:ptm>

== Site diagnostic

After the joint spectrum, peptide, and localization confidence requirements,
Sage Plus retained #s("report.ptm.1.plus.sites") and #s(
  "report.ptm.2.plus.sites",
) assessable site events in the two HCD files. The synthesis-inconsistent counts
were #s("report.ptm.1.plus.inconsistent") and #s(
  "report.ptm.2.plus.inconsistent",
), respectively (@fig:si-ptm-diagnostic). Consistency pools unambiguous
synthesis truth across libraries because the acquisition-to-library assignment
was not independently verified. It includes identification and localization
discrepancies and is not an arrangement-level false-localization rate estimate.

#figure(
  fig("fig.supplement-ptm-diagnostic", width: 100%),
  caption: [Sage Plus synthesis consistency after joint confidence filtering.
    Bars partition assessable site events into synthesis-consistent and
    inconsistent events. Spectrum, peptide, and localization q-values must each
    be at most one percent. Counts and inconsistent fractions are labeled. The
    acquisition-to-library mapping remains unaudited, so this diagnostic does
    not establish localization-error calibration.],
) <fig:si-ptm-diagnostic>

= Attachment identity and modification reuse <sec:si-attachments>

The published Sage Plus executable passed #s("named.checks") synthetic CLI
checks. The fixtures test first and last residue placements, N-terminal and
C-terminal groups, and indexed and mass-offset search. The same named definition
shares its occurrence limit across listed sites. Tests cover discovery, guided
reuse, and iteration, checking that attachment identity survives each step.
Indistinguishable first-residue and terminal-group alternatives must not
generate reusable evidence, including at permissive reporting thresholds. These
tests assess implementation correctness rather than empirical
terminal-localization calibration. The schematic in @fig:modification-model
describes the tested identity distinction.

PTM-library and site-report schemas preserve attachment identity alongside
protein coordinates. Legacy four-column libraries represent residue evidence.
Terminal-group evidence requires an explicit attachment field. This distinction
prevents a residue annotation from authorizing a neighboring terminal group.
Migration accepts legacy modification syntax, and library-aware previews expose
eligible sites before spectrum searching. The fixture configuration, binary
identity, and check outcomes accompany the refreshed analysis records.

#pagebreak()
= Tables <sec:si-tables>

These numerical reference tables retain exact counts, denominators, and timing
ranges that complement the comparative figures. Resource summaries correspond to
@fig:timing and @fig:scaling. Identification counts correspond to @fig:overlap,
and quantitative endpoints correspond to @fig:lfq-endpoints and @fig:control.

#figure(
  tbl("tbl.scientific-timing"),
  caption: [Public timing medians for the specified Sage and Sage Plus releases.
    Completed trials are shown against planned measured attempts, with failures
    excluded from timing summaries. Parentheses contain observed minimum and
    maximum wall times, not confidence intervals. Peak memory is maximum
    resident set size in MiB. Target PSMs pass the reported one-percent spectrum
    q-value threshold. Warmups are excluded.],
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
  caption: [Matched worker-scaling summaries. Completed measured trial counts
    are shown for each engine and worker setting. Seconds and peak resident MiB
    are medians. Time and PSM ranges span measured repeats. A single PSM count
    indicates identical yield in every repeat. Warmups remain in the execution
    records but are excluded from these summaries.],
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

#figure(
  tbl("tbl.mass-offset"),
  caption: [Mass offset cost, index size, and accepted identifications on one
    public HEK file. Search and wall seconds and peak resident memory are
    medians of two measured searches. Database and identification counts are
    identical across those repeats. Offset PSMs count reported candidates of any
    rank whose peptidoform carries an offset, including decoys.],
) <tbl:si-mass-offset>

#figure(
  tbl("tbl.large-db"),
  caption: [Every large-database search with its outcome. Refused searches
    stopped at the preflight estimate before building, and stopped searches hit
    the runtime memory guard. Peptides searched counts the database peptides
    after prefiltering where it applied. PSMs and peptides are accepted at one
    percent q-value.],
) <tbl:si-large-db>

#figure(
  tbl("tbl.prefilter-sweep"),
  caption: [Prefilter match threshold and peak cap on the HEK file against the
    human reference with the ten-times catalog subset. Min. matches is the
    preliminary fragment matches one precursor hypothesis needs to keep a
    peptide, and Peaks limits the prefilter to each spectrum's most intense
    peaks. Kept is the percent of streamed peptides retained. PSMs are accepted
    at one percent q-value, and the entrapment FDP is combined at the peptide
    level.],
) <tbl:si-prefilter-sweep>
