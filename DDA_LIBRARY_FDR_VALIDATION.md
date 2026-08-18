# DDA library-search FDR validation

Library-search FDR must be evaluated on spectra that were not used to construct or tune the
library. A same-run round trip is useful for software regression testing, but it is not an FDR
validation because the query spectra contributed the stored library intensities.

## Recommended experiment

1. Split acquisitions by biological run before building the library. Do not split individual PSMs
   from one run between construction and evaluation.
2. Construct a library containing genuine entries from the real proteome and an unrelated
   entrapment proteome. The entrapment entries must be deliberately populated—for example from
   separate reference acquisitions or a suitable external library. Accidental false entrapment
   identifications from a pure-target construction run are not a balanced candidate set.
3. Export the library only from construction-run PSMs passing the configured spectrum and peptide
   q-value thresholds.
4. Search untouched evaluation runs with `library_search`, internally generated shuffled decoys,
   and the intended isotope-error range.
5. Run the bundled report over `results.sage.parquet`. The protein regex must identify only the
   entrapment proteome. For a human/yeast experiment, `_YEAST` is a typical regex.

```bash
scripts/validate-library-fdr.sh \
  heldout/results.sage.parquet \
  construction/spectral_library.sage.parquet \
  '_YEAST'
```

The report derives the non-entrapment/entrapment scaling factor from distinct entries in Sage's
library Parquet file. An explicit ratio can be supplied as a fourth argument when the experimental
design requires a different weighting. Inspect the full curve rather than only 1%. Material
anti-conservative deviation should block production claims and trigger review of decoy generation,
score calibration, and library/source separation.

Repeat the experiment across independent instruments, gradients, library sources, and modification
profiles. Preserve the construction/evaluation manifests with the results so source-spectrum
leakage can be audited.
