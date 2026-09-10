Sage Plus hardening and scientific validation plan

Prepared September 4, 2026, following review of commit `5aeb568` and the existing uncommitted benchmark and paper work. Release target: `v0.1.0-beta.3`. The version is prepared on `codex/release-beta3`. Implementation status is recorded below.

The plan has two milestones. Milestone A delivers a bounded hardening release with reproducible regression evidence. Milestone B broadens scientific validation and evaluates changes to learned scoring. Experimental status remains appropriate until the evidence supports more specific claims about supported workflows.

**Scope and starting evidence**

The review passed 390 Rust tests with default features, 390 with minimal features, strict Clippy, formatting, five Python harness tests, and release-version checks. It independently reproduced four production defects and one defect in the unfinished FDRBench harness:

- API batching processed files in groups of two despite a configured batch size of one.
- Gzip writing returned success for a stream that an independent decoder rejected as truncated.
- Concurrent event output contained descending sequence numbers.
- Percentage precursor tolerances passed validation and then caused a scoring panic.
- The FDRBench harness reused existing results despite nonexistent binary and configuration paths.

These are the first implementation targets. The baseline tests do not establish scientific calibration, complete platform coverage, or absence of dependency vulnerabilities.

**Milestone A: beta.3 hardening**

| Package | Work | Completion criteria |
|---|---|---|
| A0: Establish the baseline | Inventory existing changes and validation assets. Record exact baseline commit, binary hashes, effective configurations, and representative outputs. Preserve the existing paper and benchmark work as identifiable changes. | Baseline can be reconstructed from recorded inputs. Existing results remain identifiable and are never relabeled as fresh evidence. |
| A1: Fix execution and output defects | Correct batching, gzip completion, and concurrent event ordering. | Regressions fail on the reviewed implementation and pass after the fixes. CLI, API, and MCP behavior agree. |
| A2: Harden validation and failure handling | Reject unsupported precursor tolerances before work begins. Audit adjacent input constraints, event readers, cancellation, and completion records. | Invalid logical settings return actionable errors. Interrupted jobs cannot appear successfully completed. Valid historical manifests remain readable. |
| A3: Repair the scientific harness | Replace existence-only caches with verified manifests. Make tools and engine formats explicit. Validate complete experiment matrices. | Changed inputs invalidate reuse. Failed, missing, corrupt, and partial outputs cannot become successful observations. |
| A4: Strengthen CI | Add minimal-feature and Python tests, platform smoke tests, and scheduled vulnerability scanning. | Relevant checks run on pull requests or the documented schedule. Existing lint, MSRV, and coverage gates continue to pass. |
| A5: Add minimal provenance and run comparison | Record effective execution settings and model outcomes. Add a developer comparison report using typed outputs. | Every candidate comparison is attributable, differences are inspectable, and summaries describe actual execution. |
| A6: Validate and package | Rerun representative workloads and existing bounded FDRBench experiments through the repaired harness. Prepare release notes and exercise release packaging. | All release gates below have recorded evidence. |

Recommended merge order: A0, A1, A2, A3, A4, A5, A6. Keep each defect and its regression test in a focused change. Record test commands and material compatibility changes with each implementation package.

**A1 implementation contracts**

Batching: resolve effective file batch size once. Explicit CLI or MCP request overrides configuration, which overrides the existing default. Use the resolved value for initial search, prefiltering, and deferred postprocessing. Report it accurately. Keep worker thread count distinct. Preserve the public `JobOptions.parallel` compatibility path where possible, and document any required migration instead of silently changing its meaning.

Acceptance: three inputs with batch size one produce three batch boundaries through CLI, API, and MCP. Overrides win consistently. Changing batch size preserves analytical results on a fixed fixture, including annotation and prefilter paths.

Gzip: finish the encoder before obtaining bytes. Verify empty, small, and larger round trips using a decoder that checks complete stream termination. Exercise the exported writer API, not only an internal compression helper.

Events: allocate sequence numbers and elapsed timestamps under the write lock. Verify concurrent monotonicity and paginated retrieval without omissions. Readers may defer a trailing partial record while a worker is active. Malformed complete records must still surface as errors.

**A2 implementation contracts**

Reject percentage precursor tolerances with an explanatory error for this release. Supporting another tolerance throughout the scoring pipeline would be a separate feature. Validate finite, ordered mass bounds and strictly positive physical charges without rejecting valid asymmetric tolerances. Resolve default peak limits before comparing them. Check PTM-library prerequisites and unsupported model settings consistently.

Maintain the documented distinction between logical validation and filesystem preflight. `--validate-only` should not begin reading large spectra or building the database.

Add focused failure tests around cancellation during search and postprocessing, event-output failures, and worker failure before completion. Make persisted job records atomic where replacement is supported, with platform tests. A completed manifest must reference a complete output inventory. Audit reuse of output directories for stale optional files, then define an explicit compatible overwrite policy before changing CLI behavior.

**A3 implementation contracts**

Create one versioned manifest format for benchmark provenance where practical, reusing the existing harness helpers. Include binary checksums, code identity including dirty source state when relevant, effective configuration, spectra and database hashes, library inputs, evaluator identity, random seeds, resource settings, and expected outputs. Resolve executable paths before checking caches.

Cache keys must include every input that affects the stage. Mark a stage complete only after successful exit and output validation. Write completion metadata atomically. Unknown or mismatched old cache entries require fresh runs. Keep previous experiments in separate directories.

Replace binary-filename heuristics with explicit engine and output-format selection. Make DuckDB, Java, the evaluator JAR, and dataset paths configurable. Add dependency preflight and a documented small smoke command that works outside the original maintainer's machine.

The analyzer must distinguish a completed run with zero discoveries from a missing, failed, or corrupt run. Require the planned seed and condition matrix before aggregation. Test cache invalidation, custom paths, failed timing records, incomplete matrices, tie handling during threshold extraction, and zero discoveries.

**A4 and A5 implementation contracts**

CI: keep Rust 1.88 and the pinned toolchain checks, default tests, strict Clippy, formatting, schema synchronization, and the 80% workspace line-coverage gate. Add minimal-feature tests and all benchmark unit tests on Linux. Add bounded Windows and macOS smoke coverage for configuration, paths, input reading, and worker execution. Retain the full architecture matrix in release packaging. Add a scheduled whole-lockfile advisory scan and document narrowly justified exceptions. Avoid unrelated dependency upgrades in the defect-fix changes.

Repository hygiene: narrow broad generated-file ignore rules only after checking what they would expose. Explicitly retain generated data and local experiment exclusions. Add a small redistributable regression fixture with its source, checksum, and permitted use recorded.

Provenance: extend existing run artifacts with actual batching, Rayon thread count, software/build identity, stage outcomes, warnings, and per-file model status. Record mass alignment as fitted or skipped, with reasons, rather than unconditionally marking it applied. Follow existing schema-version and backward-reader policies. Input content hashes are mandatory for scientific benchmarks. Ordinary searches should expose hashing cost explicitly and may record a weaker identity mode when full hashing is disabled.

Comparison: begin with a developer command or script, plus JSON and Markdown reports. Compare stable peptide or peptidoform identities and spectrum keys rather than generated PSM IDs. Report additions, losses, overlap, confidence changes, model fallbacks, runtime, peak memory, and optional artifact differences. Compare matching score semantics only. Include decoys and sufficient score ranges in validation outputs.

**Milestone A release gates**

- All four reproduced production defects have regression coverage. The harness cache reproduction is closed.
- Relevant CI checks pass on the exact candidate commit, including coverage and supported feature configurations.
- No unresolved critical or high dependency advisory affects the shipped configuration without a documented applicability assessment and disposition.
- Conventional, modified-peptide, feature-heavy, and exact-prefilter workloads run with recorded inputs and fresh candidate outputs. Preserve the existing one-warmup, three-trial performance protocol initially.
- Existing harness thresholds of more than 10% runtime or RSS increase and more than 1% identification-count loss trigger investigation. They are review thresholds, not proof of scientific equivalence. Any intended analytical difference must be explained using identities and scores.
- Exact prefiltering produces equivalent canonical analytical results on the selected fixtures. Compare byte hashes when formats and metadata are identical, and typed rows when provenance metadata necessarily differs.
- The bounded entrapment experiment is rerun after cache repair against the pinned baseline and candidate. Any adverse shift is investigated. Low-power results are reported as inconclusive where appropriate.
- Release documentation states the tested workflows and limitations. It does not broaden calibration claims based on passing engineering tests.
- The manual release workflow successfully builds supported archives and the container. Inspect and smoke-test packaged binaries. Confirm version, schemas, checksums, and documentation before the later publishing action.

If an unresolved issue compromises default search correctness or reported confidence, hold the release. A missing specialized validation experiment limits the corresponding feature claim. It does not automatically require changing default algorithms or expanding this milestone indefinitely.

**Milestone B: scientific validation**

First deliverable: a committed protocol and dataset manifest before evaluating model alternatives. Choose publicly accessible data with compatible formats and documented reuse terms. Inventory available compute, then run a small feasibility pilot. Set the final replication and resource budget from the pilot without using the final evaluation results to choose favorable settings.

Target at least two independent DDA studies with multiple files and complete appropriate reference proteomes, covering different acquisition conditions where feasible. Retain the existing HEK subset experiment as a historical stress case. Additional shuffled databases measure entrapment variability and do not substitute for independent biological or acquisition datasets.

| Validation track | Evidence to collect | Required distinction |
|---|---|---|
| Core identification | PSM and peptide calibration curves, discoveries, entrapment estimates, runtime, and memory | Evaluate baseline and candidate under matched supported configurations and estimator assumptions. |
| Protein inference | Protein and group-level calibration with shared-peptide cases | Peptide-level evidence does not establish protein-group calibration. |
| PTMs and labeling | Modified/unmodified and label-group strata, search-space caps, and partner consistency | Adequate pooled calibration can hide a weak subgroup. |
| Localization | Known-site or otherwise independently grounded localization benchmarks | Identification FDR and localization error are separate endpoints. |
| LFQ and transfers | Known-ratio mixtures, missingness, ratio error, false-transfer evidence, MBR on/off | Observed MS2 and transferred evidence require separate reporting. |
| Spectral-library export | Entry-level selection, structural correctness, consensus reproducibility, and independent evaluation where claimed | Transition rows are not independent library entries. Reusing source spectra is not independent validation. |

Use nominal 1% as the primary identification threshold, with a prespecified secondary curve such as 0.1%, 0.5%, 2%, and 5%. Before final runs, specify each estimator, denominator, grouping rule, treatment of zero discoveries, uncertainty method, and a practically meaningful inflation margin. Determine sample size or replication sufficient to assess that margin. Preserve file and study dependence in uncertainty calculations. Report uncertainty or insufficient power explicitly instead of labeling a noisy result calibrated.

The protocol should distinguish within-study uncertainty, between-study consistency, and exploratory subgroup findings. Freeze the final evaluation set. Retuning after inspecting it creates a new experiment requiring fresh evaluation evidence.

**Milestone B: learned scoring experiments**

Audit the full sequence of provisional selection, mass calibration, RT alignment, base RT prediction, PTM offsets, mobility prediction, LDA, and error estimation for reuse of evaluation information. Existing mobility cross-fitting and RT offset folds are useful starting points.

Prototype grouped held-out base RT predictions first, then rescoring. Group related spectra and peptide variants as required to prevent leakage, including the relevant label and decoy relationships. Learn preprocessing and model-selection decisions within the appropriate training partition. Check comparability of scores pooled from different folds and apply transformations consistently to targets and decoys.

Compare current behavior, held-out RT alone, held-out rescoring alone, and the combined approach. Keep settings explicit and experimental until the frozen evaluation protocol supports promotion. Measure calibration first, then yield at comparable error control, stability, time, and memory. Do not select a method solely because it reports more identifications at its own q-value threshold.

Add invariance tests for ties at each FDR level, sparse decoys, label groups, PTM variants, input permutations, thread counts, and prefilter equivalence. Correct an estimator only after its intended statistical contract is written down and independently exercised.

Method references: [FDRBench entrapment assessment](https://www.nature.com/articles/s41592-025-02719-x) and [cross-validation in shotgun proteomics](https://noble.gs.washington.edu/papers/granholm2012cross-validation.html). These support the evaluation design, not a finding that Sage Plus is currently miscalibrated.

**Feature follow-through and deferred work**

After A5, extend the existing HTML QC report with per-file fits, skipped-model reasons, residual diagnostics, missingness, and transfer evidence. Make the developer run comparison available through normal CLI and MCP workflows when its schema and semantics have settled. Expose CLI resource preflight by sharing MCP estimation logic.

Defer reusable database indexes, broad module rewrites, major scoring-default changes, and streaming LFQ or Parquet redesign to separately profiled changes. These are valuable candidates, but each has its own compatibility or scientific validation burden. Use measured stage costs to select performance work.

**Execution checklist**

- [x] A0: Freeze baseline source, binaries, input hashes, and existing work in the local evidence directory.
- [x] A1: Implement batching, gzip, and event-order fixes with regressions.
- [x] A2: Implement validation, partial-event handling, atomic MCP records, cancellation boundaries, and explicit local overwrite policy.
- [x] A3: Implement verified harness reuse, explicit tools and formats, and analyzer integrity tests.
- [x] A4: Configure missing CI checks and add an original redistributable fixture. Hosted platform execution remains a release gate.
- [x] A5: Add truthful model outcomes, metadata provenance, and a typed run comparison command.
- [ ] A6: Record regression evidence and complete beta.3 packaging checks.
- [ ] B1: Commit protocol, dataset manifest, feasibility results, and final resource budget.
- [ ] B2: Run the expanded baseline and candidate validation matrix.
- [ ] B3: Evaluate held-out modeling alternatives under the frozen protocol.
- [ ] B4: Publish scoped findings in repository documentation and choose the next algorithm or feature release based on evidence.

The hardening implementation has passed local default and minimal-feature tests, strict Clippy,
Rust 1.88 checks, and the workspace line-coverage gate. Fresh representative and entrapment runs
are recorded in [the hardening results](benchmarks/HARDENING_RESULTS.md). A scientific protocol and candidate dataset inventory are
available in `benchmarks/SCIENTIFIC_PROTOCOL.md` and `benchmarks/datasets.json`. Their pending
fields deliberately keep B1 open.

The September 10 dependency patch resolves the XML findings and the fresh online audit passes.
Local validation, hosted Rust CI, all seven release archive builds, and the container build pass.
The preparation is now in the existing repository under the configured owner account.
Current repository PR checks, the final dry run, review and merge, and tagging remain release gates. Follow [the beta.3 checklist](benchmarks/BETA3_RELEASE.md).
No scoring-default change or external dataset calibration claim has been made.
