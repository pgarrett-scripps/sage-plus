#import "stats.typ": lit, s
#import "assets.typ": fig, tbl

This Supporting Information records the inputs, software identities, and
additional paired results for the Sage and Sage Plus versions named in the main
text.

= Frozen software and evidence <sec:si-provenance>

The software pair is Sage `v0.15.0-beta.2` and Sage Plus `v0.1.0-beta.3`. The
tag names belong to different projects. Version identity is fixed for this
report, not updated automatically when a repository changes. Release metadata
was checked on September #lit("15"), #lit("2026").

The Sage source commit is:

```text
df9219951cc9a54cf4cd55d76541af24b687bd3d
```

The Sage executable has the following Secure Hash Algorithm 256-bit (SHA-256)
digest:

```text
065f9c7b445d2d5f4d3f23be007f123231e135b0344a586a656666585f8c402e
```

The Sage Plus release tag resolves to source commit:

```text
9bfc8acbc3002d83aa3f27af5f848376b5c5c9be
```

Sage Plus executable digest:

```text
44ed3fbde2a159b9bb4cdbfa1bf70d4e50a450aeb4a1e6d228efe81b857f43f1
```

FDRBench `1.1.1` uses source commit:

```text
3b619a9acf60d7292fb651a00da55f58cb67fb79
```

The frozen summary is
`benchmarks/scientific-results/20260914/pilot-summary.json`. Its environment
block records executable hashes and the host. The adjacent `FINALIZATION.json`,
`evidence-audit.json`, and `archive-verification.json` record the evidence
checks. Detailed results, reference receipts, conversion records, and frozen
plans remain under `/data/sage-plus-scientific/20260914`.

The local archive was reopened and #s("archive.members") members were
individually verified. Large original and converted spectra remain external to
that archive, with their identities recorded in receipts and manifests. A local
archive is not an external deposition. Rebuilding tables from the committed
summary is distinct from reacquiring inputs and repeating the searches. Original
failures and amended configurations remain separately identifiable.

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

The amended `hye-irt-defined.fasta` excluded each listed entry for both engines.
It did not replace unknown residues with invented amino acids or retain the
known subsequences of excluded proteins. The original comparison plan and failed
Sage Plus outputs were retained alongside the amended `public-comparison-v2`
results. Unamended Sage results were not substituted for the matched
amended-reference comparison.

The community sample annotation named a different yeast species from the primary
methods. The reagent identified in those methods was Promega V7461. Its
manufacturer documentation identified _S. cerevisiae_, which determined the
pilot's reference choice. #link(
  "https://www.promega.com/-/media/files/resources/protocols/technical-manuals/101/ms-compatible-yeast-and-human-protein-extracts-protocol.pdf",
)[Manufacturer technical manual TM410] records that reagent identity. Reference
completeness and the precise _E. coli_ strain remain limitations rather than
assumptions silently resolved by the amendment.

= Additional paired resource measurements <sec:si-resource>

@tbl:si-enlarged summarizes the entrapment workloads. Each engine has one search
per selected file and construction seed. These medians combine different inputs,
so they are descriptive workload summaries rather than repeat timing estimates.

#figure(
  tbl("tbl.scientific-enlarged"),
  caption: [Resource use with entrapment-expanded references. Each engine
    searched two selected files with three construction seeds per study. Medians
    span those file-by-seed searches. The mixture reference includes the
    documented exclusions. MiB describes peak resident memory, not the process
    address-space limit.],
) <tbl:si-enlarged>

The local HEK conventional and common-modification results appear in
@tbl:si-local. These searches used the same frozen release pair as the main
text. Acquisition could overlap these runs, so their wall times are not part of
the controlled public timing endpoint. The common-modification workload adds
methionine oxidation and peptide N-terminal acetylation under the shared search
configuration. The local biological input has unverified redistribution
provenance and is not offered as a public replication dataset.

#figure(
  tbl("tbl.scientific-local"),
  caption: [Contextual local HEK comparison of the specified Sage and Sage Plus
    releases. Values are medians of three measured trials after one warmup.
    Target PSMs pass the reported one-percent spectrum q-value threshold.
    Possible overlap with acquisition activity limits timing interpretation.],
) <tbl:si-local>

The broad PTM workload did not complete with either engine under the selected
limits. Plus applied its additional modified-database guard and Sage encountered
allocation failure under the shared address-space ceiling. Failed warmups and
attempts are retained. Their elapsed time is not summarized as the time required
to finish the search.

= Threshold-dependent entrapment results <sec:si-thresholds>

@tbl:si-thresholds retains every prespecified nominal threshold from the
independent audit. Target counts and FDP estimates are averaged across the same
file-by-seed combinations. A higher target count at a nominal q-value is not a
matched-error sensitivity result.

#figure(
  tbl("tbl.scientific-thresholds"),
  caption: [Threshold-specific peptide entrapment results. Each row averages two
    files and three shared construction seeds for the named engine and study.
    FDP uses the conservative paired bound over exact target-partner score ties.
    Nominal q-values and independently estimated FDP are reported in percent.],
) <tbl:si-thresholds>

Construction seeds were `20260914`, `20260915`, and `20260916`. The bootstrap
used seed `20260914` and #s("pilot.bootstrap") resamples. It resampled the file
and construction-seed dimensions while preserving engine pairing. The bounds in
the main text use the #lit("2.5")th and #lit("97.5")th percentiles. These are
conditional descriptive intervals over the selected files and seeds.

The audit checked counts against official FDRBench output and explicitly
considered exact score ties. Empty discovery denominators remain undefined. This
convention prevents absent evidence from being presented as a zero-error
success. Search failures are also retained as failures rather than converted to
empty discovery sets.

= Public identification counts and worker scaling

@tbl:si-public-counts gives the absolute counts underlying the main
identification figure. Peptidoform counts use peptide q-values and PSM counts
use spectrum q-values. These filters are distinct. Decoy PSM counts document the
retained confidence-filtered output, but they are not substituted for the
independent entrapment estimator.

#figure(
  tbl("tbl.report-public-counts"),
  caption: [Per-file accepted identifications at the one-percent threshold of
    the indicated confidence level. The common normalizer reads underlying
    result tables from both engines. Input indices follow @sec:si-inputs.
    Peptidoforms are counted separately from spectra.],
) <tbl:si-public-counts>

@tbl:si-scaling reports every engine and worker-count summary from the
extension. Each range spans three measured repeats following one warmup. The PSM
range is shown alongside time so that a repeat with changed yield cannot be
hidden by a single aggregate timing value. Assignment differences are analyzed
separately in the main text.

#figure(
  tbl("tbl.report-scaling"),
  caption: [Matched worker-scaling summaries from the report extension. Seconds
    and peak resident MiB are medians. Time and PSM ranges span measured
    repeats. Warmups remain in the execution records but are excluded from these
    summaries.],
) <tbl:si-scaling>

= Quantitative comparison denominators

@tbl:si-lfq-shared restricts the accuracy comparison to peptide and preparation
pairs accepted by both releases. Combined-charge precursors use the same empty
charge key in both formats. Numeric modification spelling is normalized without
changing mass values or modification positions. The shared set is intersected at
the ratio-pair level, so both conditions have positive accepted intensity in
both engines.

#figure(
  tbl("tbl.report-lfq-shared"),
  caption: [Ratio accuracy on identical accepted peptide-preparation pairs.
    Absolute error is measured against the expected species log-base-two B/A
    ratio. The last column is the median absolute difference between the two
    engines' observed ratios on the same pairs. These metrics assess common
    features and complement the engine-specific sets in the main text.],
) <tbl:si-lfq-shared>

The control denominator in @tbl:si-control includes positive finite accepted
precursor-file intensities with an unambiguous species assignment. A foreign
sequence also occurring in the human reference after isoleucine/leucine
normalization is excluded. Retention-time standards and unmapped proteins are
outside the species-control denominator. The numerator is a subset of the
denominator under the same rule. The absence of direct evidence is computed from
rank-one PSM results at the joint spectrum and peptide threshold. It does not
use an exported confirmation flag whose meaning may differ between releases.

#figure(
  tbl("tbl.report-control"),
  caption: [Exact human-only control denominators at the primary LFQ threshold.
    Foreign means a yeast or E. coli assignment without a human-reference
    equivalent. Direct MS2 requires matching peptide and file with both spectrum
    and peptide q-values at most one percent. No direct MS2 is a subset of the
    complete control set, not an additional independent sample.],
) <tbl:si-control>

The primary LFQ threshold applies to a precursor-level score shared across
files. A positive intensity in a file is therefore not accompanied by an
independently calibrated recipient-file transfer probability in this analysis.
The human-only control tests the intended absent-species condition. The main
Limitations section describes the boundaries of its error interpretation.

= Synthetic-site diagnostic retained separately

@tbl:si-ptm-diagnostic records the released Sage Plus site output after spectrum
and localization filtering, deliberately before peptide filtering. All primary
jointly accepted site sets were empty. The diagnostic rows therefore explain
what was withheld by the primary confidence requirement. They are not an
alternative validated result obtained by relaxing that requirement.

#figure(
  tbl("tbl.report-ptm-diagnostic"),
  caption: [Secondary Sage Plus synthesis-consistency diagnostic from the same
    released executable. Diagnostic events pass spectrum and localization
    q-values at one percent but omit the peptide filter. The final column
    restores the primary joint requirement. Consistency is assessed against
    pooled unambiguous synthesis truth, with file-to-library mapping unaudited.
  ],
) <tbl:si-ptm-diagnostic>

The synthesis reference and truth were prepared before the original searches.
The original all-site configurations are retained for both releases. Oracle
site-prior runs used information from synthesis truth and are excluded from
comparative accuracy claims. Development repairs after the released snapshot are
also excluded. A result from a modified executable must be assigned its own
software identity and independently re-evaluated before entering a later report.

= Report extension and analysis provenance

The extension plan is `paper/analysis/data/report-extension/plan.json`. Its
execution directory is `/data/sage-plus-scientific/report-extension-20260915`.
The frozen runner records executable and input hashes, the host, commands,
process limits, and output identities. Its matrix-status file retains the
completion state of every planned job. The LFQ extension has one search per
engine and does not support a repeated runtime estimate. Worker scaling has one
warmup and three measured trials per engine and worker count.

All #s("report.extension.completed_jobs") planned extension jobs completed. The
closing audit verified #s("report.extension.source_files_verified") source and
output file identities.

The new analysis does not modify the original pilot summary or its analysis
scripts. `collect_report.py` reads the retained raw results, normalizes the two
output representations, and records SHA-256 identities of the source files in
separate report snapshots. Public disagreement classifications are checked
against the frozen engine-only counts. Quantitative ratio counts are checked
against the common scientific metric function. The shared-set calculation uses
identical peptide, charge, and preparation keys across releases.

The additional analyses were motivated by gaps identified after inspection of
the pilot. They reuse previously selected public files and are explicitly post
hoc extensions. Neither their additional compute repeats nor their many
precursor measurements increase the number of independent biological studies.
The original pilot, extension measurements, and engine-specific diagnostics
remain distinguishable in the figures and captions.

= Regenerating the report <sec:si-repro>

The report reads the frozen pilot summary through
`analysis/scripts/_scientific.py`. The statistics generator supplies computed
values used in prose. Table generators read the same inputs, and the figure
generators plot the same frozen evidence and new report snapshots. Each asset
records the analysis code and data identities that produced it. These steps do
not rerun biological searches or modify the finalized evidence.

From the repository root, rebuild the report with:

```shell
cd paper
just assets
just fmt
just docx
just paper
just verify
just check-stats-deep
```

To repeat the biological analysis, first inspect
`benchmarks/SCIENTIFIC_PILOT.md` and the retained acquisition, conversion,
reference, and search plans. Supply the pinned executables and verified inputs
at the paths expected by a copied plan. Use a new output directory for a new
experiment. Do not overwrite the frozen evidence or describe a changed binary as
the originally evaluated release.

The small report snapshots can be regenerated from the retained raw results:

```shell
cd paper/analysis
uv run scripts/collect_report.py public
uv run scripts/collect_report.py ptm
uv run scripts/collect_report.py lfq
uv run scripts/collect_report.py scaling
```

These commands derive tables and distributions from existing searches. Repeating
an executable benchmark is a different operation that requires the pinned
binary, complete input files, and a newly named output directory. The report's
figures can be rebuilt from the snapshots without copying the large spectra into
the report directory.
