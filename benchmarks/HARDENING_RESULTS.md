Hardening release evidence

Historical evidence prepared September 4, 2026. For the September 10 dependency fix, beta.3 version, and final candidate checks, see [the beta.3 release checklist](BETA3_RELEASE.md). The following records the earlier candidate as tested.

At that snapshot, the core hardening implementation was ready for review. Publication remained on hold for dependency disposition and hosted platform validation. The workspace version was `0.1.0-beta.2`. No scoring default or paper result was changed by this hardening work.

**Identity and scope**

The baseline commit is `5aeb56852fc176e9bbc5ae82e9d931927480760a`. Baseline and candidate snapshots include the respective dirty working trees, with source and input hashes. Existing benchmark and paper work was preserved. The search binary used for every candidate comparison has SHA-256 `0fd79d98b6ee192f8cea9f95203504c89108af6e874df85f9a0c40ab3ac72acf`. The frozen baseline binary has SHA-256 `ff7d837db3bc5d65b2a50cd9ce1017807340f9615ff8ca285f9baae840cce76b`.

Local evidence is intentionally excluded from Git under `benchmarks/results/`. Keep the source archives and binary manifests with any shared result bundle. Subsequent MCP-only UTF-8 handling was separately built and tested. It did not change the frozen search binary or its dependencies. The final Linux package includes that MCP correction.

**Engineering checks**

| Check | Result |
|---|---|
| Default workspace tests | 399 passed |
| Minimal-feature workspace tests | 399 passed |
| Python benchmark tests | 16 passed |
| Workspace line coverage | 81.70%, above the 80% gate |
| Strict Clippy | Full workspace passed, followed by the final MCP-specific check |
| Rust 1.88 | Workspace compatibility check passed |
| Formatting, schema synchronization, release version inheritance | Passed |
| Workflow YAML | Parsed locally |
| Native Linux archive | Built, unpacked, hashes verified, synthetic search passed |

Regression tests reproduce and close configured batching, incomplete gzip streams, concurrent event ordering, and unsupported precursor percentage tolerances. Additional checks cover invalid logical settings, cancellation, atomic MCP records, partial event records including UTF-8 boundaries, explicit local overwrite, cache integrity, failed timing records, and incomplete experiment matrices.

Hosted Windows, macOS, architecture-specific release archives, and container packaging have not been executed from this working tree. The full native release workflow is also deliberately blocked by the unsuppressed audit findings. Local results do not substitute for those gates.

**Representative workloads**

The initial suite ran one warmup and three measured trials per engine for four conditions, with eight Rayon workers, batch size one, matched inputs, sequential execution, and alternating engine order. The 219 MB local HEK mzML and full local reviewed-human FASTA were reused. The modified condition uses bounded oxidation and N-terminal acetylation. The feature condition enables SILAC, LFQ, annotations, and library export.

| Condition | PSMs at 1%, both builds | Peptides at 1%, both builds | Analytical comparison |
|---|---:|---:|---|
| standard | 2203 | 1409 | Byte-identical analytical outputs |
| modified | 2256 | 1458 | Byte-identical analytical outputs |
| feature | 3266 | 1524 | Byte-identical analytical outputs |
| prefilter | 2203 | 1409 | Byte-identical analytical outputs |

Every stored PSM identity, compared q-value, hyperscore, and discriminant score agrees in these four first-trial paired outputs. In the feature condition, LFQ, matched fragments, both spectral-library formats, and main PSM output are byte-identical. Exact prefiltering also produces byte-identical candidate PSM output relative to the unfiltered search. JSON summaries differ intentionally because provenance and output paths differ.

The first timing set flagged the standard and feature conditions. Standard medians were 10.01 versus 13.09 seconds, and feature medians were 15.75 versus 17.49 seconds. Individual standard times ranged from 6.29 to 18.51 seconds including warmups. These measurements are retained rather than discarded.

Five additional paired trials per flagged condition, with their own warmups and CPU/load diagnostics, did not reproduce either threshold violation:

| Condition | Baseline median wall | Candidate median wall | Change | Baseline CPU | Candidate CPU | RSS change |
|---|---:|---:|---:|---:|---:|---:|
| standard | 6.22 s | 6.26 s | +0.64% | 22.61 s | 22.60 s | +0.32% |
| feature | 15.62 s | 15.29 s | -2.11% | 52.46 s | 52.27 s | -0.70% |

Interpretation: the flagged slowdown was not reproducible in the follow-up. Similar CPU costs and variable host timing support treating the initial wall-time changes cautiously. This is bounded regression evidence, not a precise speedup claim or a guarantee on other workloads.

**Entrapment stress case**

Twenty fresh paired entrapment databases and 40 searches completed through the repaired harness. Seeds run from 20260902 through 20260921. Both engine labels explicitly identify Sage Plus builds. The FDRBench JAR, Java executable, DuckDB executable, search binaries, configurations, input files, and stage outputs have recorded content hashes.

Target and entrapment counts at all 200 engine/seed/cutoff combinations were independently reconstructed from extracted peptide identities and the pair tables. Every count matched. Retained PSM identities, accepted identities, threshold target counts, threshold entrapment counts, and paired FDP estimates agree between builds at every evaluated cutoff.

The largest peptide-q difference across paired stored PSMs is 4e-8, and the largest discriminant difference is 5e-7. Spectrum q-values and hyperscores agree exactly. A repeat with the unchanged baseline also produced peptide-q differences up to 2e-8 and discriminant differences up to 3e-7, without identity or cutoff changes. Bitwise floating-point reproducibility therefore remains a separate concern from the observed hardening comparison.

| Nominal peptide q cutoff | Mean paired FDP, both builds | Seed percentile range | Mean target discoveries |
|---|---:|---:|---:|
| 1% | 0.58% | 0.00% to 3.27% | 103.10 |
| 2% | 2.05% | 0.00% to 5.88% | 126.90 |
| 5% | 7.89% | 4.16% to 14.10% | 142.55 |
| 10% | 15.87% | 10.20% to 22.94% | 153.75 |
| 20% | 26.01% | 18.38% to 34.06% | 172.75 |

These are empirical 2.5th and 97.5th percentiles across database seeds, not confidence intervals for a study mean. At 1%, 18 of 20 seeds have discoveries. The aggregation convention contributes zero for a completed empty discovery set and rejects missing runs. The raw per-seed summaries retain empty-cutoff values explicitly.

The experiment uses one spectrum file and an incomplete 2,000-protein target subset. It has low power near 1%, and the higher-cutoff estimates remain above nominal in both builds. This does not certify calibration on complete proteomes, independent studies, protein groups, PTM sites, transfers, or library entries. The new protocol keeps those claims gated on broader evidence.

A second full harness invocation reused the completed stages. All 140 stage-log timestamps remained unchanged. Cache integrity tests additionally cover changed or missing inputs, corrupted outputs, settings changes, input mutation during execution, and failed stages.

**Release disposition and compatibility**

The audit updates resolve the direct parser findings and the crossbeam-epoch, h2, anyhow, and memmap2 findings. Two high transitive XML advisories remain in the legacy cloud dependency tree. Normal search entry points reject remote Bruker URLs, but broader library applicability and the dependency path still need disposition. The scanner continues to fail. A yanked chacha20 dependency and three unmaintained dependencies also remain tracked. See [the security review](SECURITY_REVIEW.md).

Local output reuse now requires explicit overwrite and removes known Sage artifacts while preserving unrelated files. Use separate directories for concurrent jobs and fresh remote prefixes. Run-summary schema 9 records actual model outcomes, structured warnings, workers, and explicitly labeled metadata identity. Ordinary input metadata does not replace content hashes. mzML entity references inside binary arrays now fail explicitly, with literal base64 arrays supported.

Publication and a beta.3 version bump remain pending. The next science milestone starts with [the draft protocol](SCIENTIFIC_PROTOCOL.md) and [candidate dataset inventory](datasets.json). Dataset reuse terms, exact external files, final estimator review, and pilot-based power/resource planning remain open. Held-out RT and rescoring experiments have not been implemented or evaluated.

**Evidence inventory**

- `results/hardening-baseline/`: original source snapshot, input identities, binaries, and manifest.
- `results/hardening-candidate/`: frozen search build and source snapshot.
- `results/hardening-20260904/`: paired workload outputs, commands, hashes, comparisons, and prefilter equivalence.
- `results/hardening-runtime-followup/`: diagnostic paired trials with CPU and host load records.
- `results/fdrbench-hardening-20260904/`: stage manifests, 20-seed outputs, analysis conventions, independent count checks, typed comparisons, and cache-reuse check.
- `results/hardening-reproducibility/`: repeated unchanged-baseline search and score comparison.
- `results/hardening-package-review/`: native archive, unpacked checksum verification, and synthetic smoke output.
- `results/hardening-evidence/`: engineering logs and audit JSON.
