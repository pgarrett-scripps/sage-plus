# Scientific hardening development results

Two PTM bugs are repaired. MBR transfer calibration remains unresolved. These are development comparisons on previously inspected pilot data. Reserved validation and a new release have not been performed.

The existing Typst manuscript is in `paper/paper.typ`. It was not modified. The frozen 20260914 scientific pilot remains intact, with 3,862 original integrity checks passing.

## PTM changes and evidence

A singleton decoy class made the Gaussian peptide confidence model undefined. The candidate detects underdetermined fits and uses cumulative target-decoy counts with a +1 correction, issuing an explicit warning. The first synthetic sample now accepts 1,132 peptides at 1%, compared with zero in the published binary. This is numerical recovery, not a certification of calibration in a restricted synthetic database.

The localizer could assign confidence to one of several exactly tied best target arrangements. The second candidate keeps those competitions in the dataset calculation, but assigns unresolved target arrangements localization q-value 1. Single-arrangement cases retain their confidence.

At joint 1% PSM, peptide and localization thresholds, the all-sites experiments give:

| Sample | After acceptance repair | After tie repair | Remaining disagreement |
|---|---:|---:|---:|
| HCD 1 | 1,172 consistent, 118 inconsistent | 1,128 consistent, 15 inconsistent | 1.31% |
| HCD 2 | 880 consistent, 15 inconsistent | 876 consistent, 1 inconsistent | 0.11% |

The values describe consistency with synthesis sites across all libraries. Acquisition-to-library mapping remains independently unaudited. They do not establish production localization FDR. Oracle coverage runs at 0%, 25%, 50% and 100% are retained for both samples and both candidate stages. Their truth-derived annotations are development controls.

## MBR findings

LFQ schemas 3 and 4 retain the existing precursor q-value and add strict per-file MS2 evidence and file-level signal diagnostics. Strict evidence requires both PSM and peptide acceptance. Metadata identifies the precursor q-value as applying across files and marks the new file score as experimental and uncalibrated.

A separate analyzer evaluated an experimental transfer estimator using shifted-decoy counts within each recipient file. It does not filter production output.

- At 1%, it accepted zero transfers in every file, including the mixtures. This provides no useful sensitivity or calibration evidence.
- At 5%, the pure-human control had 69 accepted target transfers. Four ambiguous sequences were excluded. Of the remaining 65, 39 were foreign, or 60%.
- The legacy 1% precursor filter retained 8,717 assessable strict-unconfirmed control signals, including 5,826 foreign sequences. This is a different confidence unit and is not directly comparable to a transfer threshold.

These results reject the candidate transfer estimator as a basis for an FDR-controlled MBR claim. Shared peak selection and decoy reference traces are plausible sources of bias. That is a hypothesis for the next experiment, not a demonstrated explanation. Sample purity, carryover and reference coverage remain control assumptions.

## Verification and next work

All 24 development searches completed. The two routine HEK result Parquet files are byte-for-byte identical to the archived beta.3 outputs. Their accepted target PSM counts are 16,917 and 19,471. Scores and peptide q-values are unchanged.

All 383 relevant Rust tests and 38 Python benchmark tests passed. Strict Clippy checks passed for the core, output and CLI packages. Formatting and whitespace checks passed. New builds and evidence are under `/data`, with debug symbols and incremental compilation disabled. These concurrent development runs are not performance benchmarks.

Next, test a transfer decoy construction with matched donor evidence and symmetric recipient selection, independently establish the PTM acquisition mapping, then freeze a successful method before reserved-data evaluation. The paper must continue to distinguish numerical correctness, synthesis consistency and calibrated error control.

No commit, pull request, release or manuscript update was made in this development pass.
