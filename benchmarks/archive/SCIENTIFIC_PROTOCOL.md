Scientific validation protocol, draft 1

Prepared September 4, 2026. This document specifies future evaluation. It does not certify
calibration or change any scoring default. The hardening HEK runs are engineering comparisons.

**Decision and endpoints**

The primary release question is whether a candidate preserves supported search behavior while
correcting execution defects. Compare it to the frozen pre-hardening Sage Plus binary. Keep the
upstream Sage comparison separately labeled. A result against upstream does not isolate the
effect of this release.

The later scientific question is whether reported peptide error control is consistent with
independent evidence under a matched, complete reference search. Use a nominal 1% peptide
threshold as the primary endpoint. Prespecify 0.1%, 0.5%, 2%, and 5% as secondary thresholds.
Retain all targets and decoys needed to reconstruct the full ranking. The historical subset
harness cuts stored rows at 20% and uses a different secondary grid. Preserve that experiment
as a stress case rather than silently changing its meaning.

**Data and feasibility**

The local HEK mzML and reviewed-human FASTA are available for feasibility and engineering tests.
The acquisition accession and reuse terms for this particular local file have not been verified.
Do not redistribute it or treat it as a newly verified independent study.

Candidate external studies are HEK293 DDA PXD001468 and the ISB18 Mix 7 Orbitrap collection.
Both appear in the data availability statement of the
[FDRBench paper](https://www.nature.com/articles/s41592-025-02719-x).
Before acquisition, verify each selected file's availability, format, sample composition, reuse
terms, and checksum. Record conversion software and settings for RAW conversions. Do not equate
technical injections or additional shuffled databases with independent biological studies.

The dataset manifest lists candidates with pending fields explicitly. A final evaluation manifest
must have at least two independent studies, multiple files per study, exact URLs, reference
proteome releases, content hashes, and an explanation of search-space completeness. Select files
by acquisition metadata before comparing candidate outcomes. Reserve separate development and
evaluation files. Known mixtures require their actual component sequences and any justified
contaminants, not an arbitrary reduced human proteome.

Initial acquisition budget: 20 GiB total downloads and 40 GiB total scratch. Initial compute
budget: four hours, sequential searches, eight Rayon workers, and at most 14 GiB configured
search memory. These are feasibility limits, not a power calculation. Pilot two files per study
where available. Record elapsed time, RSS, usable spectra, target and entrapment counts, and
failure reasons. If limits prevent an adequate pilot, publish that limitation and revise the
budget before freezing evaluation. Do not choose the budget from favorable final results.

**Estimator contract**

Use a pinned FDRBench binary and record its JAR and runtime hashes. Record generation seed,
enzyme, cleavage constraints, I/L normalization, target-to-entrapment ratio, pairing rules, and
collision checks. Evaluate whether the target and entrapment construction meets each estimator's
assumptions. Paired, combined, and lower-bound FDP outputs have distinct meanings. Report their
names and denominators separately. A low lower bound is not evidence of valid error control.
The [FDRBench implementation](https://github.com/Noble-Lab/FDRBench) is the executable reference
to review before freezing the estimator contract.

At peptide level, normalize sequence and modification identity consistently for both engines.
Define how multiple PSMs compete for a peptide and use one coherent ranking statistic and
threshold mapping. The historical subset extractor uses minimum peptide q-value and maximum
hyperscore across PSMs grouped by artificial FASTA accession. Audit that correspondence before
using it for final evaluation. Preserve ties as complete score blocks. Keep targets and decoys
subject to the same ranking and grouping rules.

Distinguish a completed run with zero discoveries from an absent or invalid result. Report zero
discoveries and an undefined FDP when its denominator is zero. Never fill a missing run with
zero. Verify the full study, file, seed, engine, and configuration matrix before aggregation.

**Uncertainty and decision rules**

Report each file and study separately, paired candidate-minus-baseline differences, and
variability across entrapment seeds. The local 20-seed experiment estimates construction
variability conditional on one spectrum file. It does not measure between-study uncertainty.

Use a prespecified hierarchical resampling analysis that preserves file and study dependence.
With only two studies, between-study generalization remains descriptive. Report an interval
within each study and the consistency of direction across studies. Do not bootstrap individual
PSMs as if independent. Fix resampling seeds, the number of replicates, and estimator-specific
interval construction before evaluation.

Proposed practical inflation margin at nominal 1%: 0.5 percentage points above nominal. This is
a proposed decision tolerance, not a demonstrated detectable effect. Use the pilot to determine
the number of files and entrapment replicates needed for an informative interval. Freeze the
margin and resource plan before final runs. If the resulting interval cannot distinguish the
margin, report insufficient precision. More reported discoveries alone is not a pass criterion.

**Separate feature claims**

Protein grouping needs a protein or group-level estimator with shared-peptide cases. PTM and
label analyses need subgroup denominators and partner checks. Localization needs independent
site truth. LFQ needs known-ratio mixtures, missingness and ratio error, with MBR enabled and
disabled. Transferred evidence must remain separate from observed MS2 evidence. Library export
needs structural checks and entry-level selection, with independent spectra for any claimed
identification benefit. None of these claims follows from a peptide entrapment result.

**Learned scoring experiments**

Freeze the validation data before prototyping model changes. Audit provisional selection, mass
calibration, RT alignment, base RT prediction, PTM offsets, mobility prediction, LDA, and final
error estimation for information reuse. Group repeated spectra, peptide variants, label partners,
and related decoys consistently when constructing folds. Fit preprocessing and model selection
only within training partitions. Check score comparability across folds.

Compare current scoring, held-out RT alone, held-out rescoring alone, and their combination.
Use explicit experimental settings with unchanged defaults. Evaluate calibration first, then
yield at comparable error control, stability, time, and memory. The motivation for holding out
evaluation information follows the
[proteomics cross-validation analysis](https://noble.gs.washington.edu/papers/granholm2012cross-validation.html).
It is not a finding that the current implementation is miscalibrated.

**Freeze gate**

Before final scientific runs, replace all pending manifest fields, finish the estimator audit,
record pilot results, finalize uncertainty and power calculations, and freeze source, binaries,
tools, configurations, dataset splits, and hashes in a versioned manifest. Any tuning after
examining final results creates a new experiment requiring fresh evaluation evidence.
