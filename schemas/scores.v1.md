# Sage score definitions, version 1

## PSM output

- `hyperscore`: X!Tandem-style fragment-match score for the candidate PSM; larger is better.
- `sage_discriminant_score`: linear-discriminant score used to order PSMs for spectrum-level target-decoy competition; larger is better.
- `posterior_error`: estimated local false-identification probability for the PSM; smaller is better.
- `spectrum_q`, `peptide_q`, `protein_q`, `protein_group_q`: monotonic minimum estimated false-discovery rates at the named aggregation level; smaller is better.

## LFQ output

- `score`: score of the cross-run LFQ peak selected from the traced isotope signal; larger is better. Its exact construction is selected by `lfq_settings.peak_scoring` and is not a calibrated probability.
- `spectral_angle`: intensity-weighted normalized agreement between the observed and theoretical isotope patterns for the selected cross-run peak, in the interval `[0, 1]`; larger is better.
- `q_value`: precursor-level q-value from cumulative target and shifted-decoy counts over LFQ `score`, with a +1 correction. It is repeated across files and does not estimate individual transfer confidence.
- `intensity`: integrated MS1 signal for the precursor/file row. Null means no positive finite signal was integrated; zero is not used as a missing-value sentinel.
- `ms2_confirmed`: whether the same precursor has an accepted target PSM in that acquisition file at `lfq_settings.peptide_q_value`. This is direct-identification evidence, not a statement that a different LFQ algorithm was used. Every LFQ intensity is produced by the same cross-run feature-tracing workflow.

## LFQ file evidence in schemas 3 and 4

- `ms2_confirmed_strict`: direct target evidence passing both PSM and peptide q-values at `lfq_settings.peptide_q_value`. The original `ms2_confirmed` field continues to require only peptide acceptance.
- `file_spectral_angle`: isotope agreement in this file at the selected shared apex.
- `file_trace_cosine`: normalized similarity to the reference trace after local warping.
- `file_rt_shift_bins`: signed local warp offset in retention-time grid bins.
- `file_score`: experimental ranking score combining isotope agreement, trace similarity and warp proximity. It is not a probability or q-value and is not used to filter output.
- `transfer_candidate`: MBR is enabled and the underlying target precursor lacks strict direct evidence in this file. Shifted decoys use the underlying target's eligibility. Null indicates no integrated signal.

All file diagnostics are null when intensity is null. Metadata identifies the precursor q-value scope and marks the file score as experimental and uncalibrated. The exact score definition and development validation are documented in `benchmarks/SCIENTIFIC_HARDENING.md`.
