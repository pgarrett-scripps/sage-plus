# Sage score definitions, version 1

## PSM output

Directions: "larger is better" or "smaller is better". Columns not listed here are descriptive;
DOCS.md ("Interpreting Sage Output") gives the direction for every PSM column.

- `hyperscore`: X!Tandem-style fragment-match score for the candidate PSM; larger is better.
- `delta_next`: hyperscore minus the next-ranked candidate's hyperscore; larger is better.
- `delta_best`: best candidate's hyperscore minus this candidate's, 0 for rank 1; smaller is better.
- `matched_peaks`, `longest_b`, `longest_y`, `longest_y_pct`, `matched_intensity_pct`, `ms2_intensity`: fragment-match evidence; larger is better. `longest_y_pct` is a fraction from 0 to 1. `ms2_intensity` is the summed intensity of matched fragments.
- `precursor_ppm`, `fragment_ppm`, `calibrated_precursor_ppm`, `calibrated_fragment_ppm`: mass errors; closer to zero is better.
- `delta_rt_model`, `delta_mobility`: absolute differences from the predicted value; smaller is better.
- `poisson`: log10 of the Poisson probability mass of the PSM's matched-peak count, with the expected count set to the mean over the spectrum's scored candidates. Always 0 or negative. It is a point probability, not a tail p-value. For top-ranked PSMs, which match more peaks than the mean, smaller (more negative) is better.
- `sage_discriminant_score`: linear-discriminant score used to order PSMs for spectrum-level target-decoy competition; larger is better.
- `posterior_error`: log10 of the estimated posterior error probability (local false-identification probability) for the PSM. Always 0 or negative; -324 marks a probability that underflows to zero. Smaller is better.
- `spectrum_q`, `peptide_q`, `protein_q`, `protein_group_q`: monotonic minimum estimated false-discovery rates at the named aggregation level; smaller is better. `protein_q` and `protein_group_q` use only peptides unique to one protein or one group; shared peptides are assigned 1.
- `spectral_angle`, `explained_library_intensity`, `explained_query_intensity` in PSM output: retained for schema compatibility and always 0.

## LFQ output

- `score`: score of the cross-run LFQ peak selected from the traced isotope signal; larger is better. Its exact construction is selected by `lfq_settings.peak_scoring` and is not a calibrated probability.
- `spectral_angle`: intensity-weighted normalized agreement between the observed and theoretical isotope patterns for the selected cross-run peak, in the interval `[0, 1]`; larger is better.
- `q_value`: precursor-level q-value from cumulative target and shifted-decoy counts over LFQ `score`, with a +1 correction. It is repeated across files and does not estimate individual transfer confidence.
- `intensity`: integrated MS1 signal for the precursor/file row. Null means no positive finite signal was integrated; zero is not used as a missing-value sentinel.
- `ms2_confirmed`: whether the same precursor has an accepted target PSM in that acquisition file at `lfq_settings.peptide_q_value`. This is direct-identification evidence, not a statement that a different LFQ algorithm was used. Every LFQ intensity is produced by the same cross-run feature-tracing workflow.

## LFQ file evidence in schemas 3 to 6

- `ms2_confirmed_strict`: direct target evidence passing both PSM and peptide q-values at `lfq_settings.peptide_q_value`. The original `ms2_confirmed` field continues to require only peptide acceptance.
- `file_spectral_angle`: isotope agreement in this file at the selected shared apex.
- `file_trace_cosine`: normalized similarity to the reference trace after local warping.
- `file_rt_shift_bins`: signed local warp offset in retention-time grid bins.
- `file_score`: experimental ranking score combining isotope agreement, trace similarity and warp proximity. It is not a probability or q-value and is not used to filter output.
- `transfer_candidate`: MBR is enabled and the underlying target precursor lacks strict direct evidence in this file. Shifted decoys use the underlying target's eligibility. Null indicates no integrated signal.

- `extraction_q_value` (schemas 5 and 6): target-decoy q-value of this precursor/file extraction, computed for every target row with an integrated signal, MS2-backed rows and transfers alike. Each target precursor's shifted decoy is evaluated at the target's own peak (same apex, integration window and reference trace), with the same per-file warp search towards that reference; the resulting paired decoy rows compete with target rows by `file_score` in one pool across files, with cumulative counts, a +1 correction and monotonic q-values. Decoy precursor rows are null. It does not replace the precursor-level `q_value`, is not used to filter output, and is not calibrated for transfers (see `benchmarks/SCIENTIFIC_HARDENING.md`).

All file diagnostics are null when intensity is null. Metadata identifies the precursor q-value scope and marks the file score as experimental and uncalibrated. The exact score definition and development validation are documented in `benchmarks/SCIENTIFIC_HARDENING.md`.
