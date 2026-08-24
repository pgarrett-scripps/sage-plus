# DDA spectral-library search

Library search is a separate search mode. A configuration must contain exactly one of `database`
or `library_search`; a FASTA is not read or required for library search.

```json
{
  "library_search": {
    "path": "library.mzspeclib.txt",
    "decoy_tag": "rev_",
    "decoy_attempts": 32
  },
  "mzml_paths": ["run.mzML"],
  "precursor_tol": { "ppm": [-10, 10] },
  "fragment_tol": { "ppm": [-10, 10] }
}
```

The input may be mzSpecLib 1.0 text or Sage's `spectral_library.sage.parquet`. Each analyte must
provide a ProForma peptide, precursor charge and mass, supported fragment annotations and
intensities, and at least one protein accession. Protein mappings are taken directly from the
library.

For every target, Sage Plus generates a deterministic composition-preserving shuffled decoy. The
precursor mass, charge, terminal modifications, protein pairing, and intensity vector are retained;
fragment m/z values are recalculated after moving residue/modification tokens. Several shuffles are
evaluated and the candidate with the lowest target-fragment overlap is selected. Targets and decoys
then use the same precursor filtering and initial spectral-angle scorer. Final q-values use a
regularized target-decoy model over spectral angle, explained intensity, calibrated mass errors,
isotope error, and aligned retention-time and mobility residuals. Searches with fewer than 20
targets or decoys retain spectral-angle scoring and emit `library_rescoring_fallback`.

The mode writes the normal `results.sage.parquet`, `results.json`, and `run-summary.json` outputs.
Spectrum, peptide, protein, and protein-group q-values are calculated within the library-search
path. The run summary records library entry and transition counts separately and reports zero
database peptides/fragments.

Library search supports the same `isotope_errors` range as database search. Current restrictions
are validated explicitly: library search does not yet support chimera/DIA wide-window search,
bitmap search, PTM localization, or spectral-library re-export. Library matches support matched
fragment export, direct library RT/mobility alignment, LFQ, and TMT quantification.

Sage Parquet libraries retain each entry's source file and spectrum. When a query filename matches
a recorded library source filename, Sage emits a `library_source_overlap` warning because that run
is useful for regression testing but not independent FDR validation.

Before calling the mode production-validated, benchmark its FDR calibration with held-out and
entrapment datasets across library sources and instruments. See `DDA_LIBRARY_FDR_VALIDATION.md`
and `scripts/validate-library-fdr.sh` for the reproducible workflow and threshold report.
