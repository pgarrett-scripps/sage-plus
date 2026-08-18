# DDA spectrum-to-library search spike

## Status

This branch contains a tested search kernel, not a user-facing search mode. It can load Sage's
native `spectral_library.sage.parquet`, index entries by precursor neutral mass and charge, and
rank candidate library spectra against a processed DDA MS2 spectrum. It intentionally does not
emit identifications or q-values yet.

The spike lives in `sage_core::spectral_library_search`. Native Parquet loading is implemented by
`sage_cloudpath::parquet::deserialize_spectral_library`.

## Implemented scoring path

1. Convert the observed precursor m/z and charge to neutral mass.
2. Use a sorted precursor-mass index to find entries inside the configured precursor tolerance.
3. Require the library and observed precursor charges to agree.
4. Match library fragments to processed query peaks inside the fragment tolerance.
5. Greedily enforce a one-to-one assignment so one observed peak cannot satisfy two library
   fragments.
6. Square-root transform intensities and calculate normalized spectral angle.
7. Report matched peaks, explained library/query intensity, precursor ppm, retention-time delta,
   and ion-mobility delta.

This reuses Sage's existing `ProcessedSpectrum` representation, including charge-deconvoluted
fragment masses, but does not reuse theoretical hyperscore. Candidate lookup is approximately
`O(log N + C)`, where `N` is the number of library precursors and `C` is the small precursor-window
candidate set.

## Why it is not connected to the runner yet

The difficult boundary is false-discovery-rate control, not peak matching.

Sage's current FDR pipeline assumes every `Feature` points to a `PeptideIx` in an
`IndexedDatabase`. Target and reversed peptides compete using the same theoretical scoring path,
then picked-peptide and protein FDR use the target/decoy relationships stored in that database.

The current empirical library export contains target entries only. Sending those scores into the
existing target-decoy machinery would be invalid because no equivalently scored library decoys
exist. Giving targets empirical intensities while giving decoys theoretical or uniform
intensities would also create a systematic score advantage for targets.

A production mode therefore needs one explicit policy:

- Require paired target and decoy spectra in the input library; or
- Generate sequence decoys and transfer target intensities to valid decoy fragments using a
  documented, validated method; or
- Define a different, empirically validated competition strategy.

The native Parquet v1 reader marks every entry as target because the schema does not yet carry
decoy provenance. A later schema should include `is_decoy` and, preferably, a target/decoy pair
identifier before standalone library search is released.

## Recommended production sequence

1. **Define and validate paired library decoys.** Add decoy identity and pairing to the canonical
   schema and mzSpecLib metadata where representable.
2. **Map library analytes into Sage's peptide model.** Prefer matching library ProForma entries to
   peptides built from the configured FASTA so peptide/protein FDR, protein grouping, and existing
   result schemas remain meaningful.
3. **Add an alternate candidate scorer.** The runner can reuse spectrum reading, batching,
   cancellation, events, and output writing, while selecting the theoretical or library scorer
   through an explicit mode.
4. **Extend `Feature`.** Add spectral angle, explained intensities, library RT/IM deltas, and source
   library entry ID. Retrain or bypass the existing LDA deliberately rather than silently treating
   these as hyperscore features.
5. **Validate before optimizing.** Use held-out DDA runs and an entrapment database to test 1% FDR
   calibration, target/decoy score symmetry, identification yield, and cross-instrument transfer.
6. **Optimize only after calibration.** Candidate indexing is already cheap; fragment matching can
   later use compact arrays, vectorization, and batch scratch buffers if profiling justifies it.

## Scope estimate

The search kernel and native reader are small. A credible end-to-end DDA MVP is a medium-sized
change because it also needs analyte mapping, paired decoys, `Feature` integration, configuration,
output columns, and validation fixtures. A production-quality mode is primarily a scientific
validation project. It should remain separate from chromatogram-based DIA library searching,
which requires peak-group extraction and substantially different algorithms.
