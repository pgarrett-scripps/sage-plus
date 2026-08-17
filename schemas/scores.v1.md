# Sage score definitions, version 1

## PSM output

- `hyperscore`: X!Tandem-style fragment-match score for the candidate PSM; larger is better.
- `sage_discriminant_score`: linear-discriminant score used to order PSMs for spectrum-level target-decoy competition; larger is better.
- `posterior_error`: estimated local false-identification probability for the PSM; smaller is better.
- `spectrum_q`, `peptide_q`, `protein_q`, `protein_group_q`: monotonic minimum estimated false-discovery rates at the named aggregation level; smaller is better.

## LFQ output

- `score`: score of the cross-run LFQ peak selected from the traced isotope signal; larger is better. Its exact construction is selected by `lfq_settings.peak_scoring` and is not a calibrated probability.
- `spectral_angle`: intensity-weighted normalized agreement between the observed and theoretical isotope patterns for the selected cross-run peak, in the interval `[0, 1]`; larger is better.
- `q_value`: precursor-level q-value from picked target-decoy competition over LFQ `score`; smaller is better.
- `intensity`: integrated MS1 signal for the precursor/file row. Null means no positive finite signal was integrated; zero is not used as a missing-value sentinel.
- `ms2_confirmed`: whether the same precursor has an accepted target PSM in that acquisition file at `lfq_settings.peptide_q_value`. This is direct-identification evidence, not a statement that a different LFQ algorithm was used. Every LFQ intensity is produced by the same cross-run feature-tracing workflow.
