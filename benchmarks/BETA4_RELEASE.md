# Beta.4 release checklist

Prepared September 19, 2026 for `v0.1.0-beta.4` on `claude/mass-offset-search`.
This release adds search-time mass offset modifications and the LFQ and
confidence-model changes carried since beta.3. It does not change scoring
defaults for existing configurations.

## Scope

- Mass offset modifications: `search_mode: "mass_offset"` searches a
  modification as a precursor and fragment offset instead of expanding it into
  the fragment index. One offset is placed per peptide and offsets are never
  combined with each other. Behavior and limits are in
  [DOCS.md](../DOCS.md#mass-offset-modifications).
- Per-file LFQ evidence in the new `lfq.v3` and `lfq.v4` schemas. The per-file
  scores are experimental and uncalibrated, and their metadata says so.
- Count-based fallbacks when the peptide and protein confidence model is
  underdetermined, and no localization confidence for equally scoring target
  arrangements.

## Compatibility

- Configurations without `search_mode` are unaffected. Result Parquet output for
  such configurations is byte-identical to the beta.3 binary on the three
  workloads checked in [the evaluation](MASS_OFFSET.md).
- `lfq.parquet` gains columns and moves to schema 3, or 4 with labels. Readers
  pinned to schema 1 or 2 must be updated. Consumers that pin the run-summary
  schema are unaffected: it remains 9, with additive `modifications` fields.
- Cascade pins Sage Plus `0.1.0-beta.2` and run-summary schema 8, so it requires
  its own upgrade before it can drive this release; the mass offset workflow
  additionally needs `search_mode` passed through its modification file.

## Validation

Local evidence, recorded September 19, 2026:

- Workspace tests, strict Clippy, and formatting pass on the merged branch.
- Mass offset evaluation: [MASS_OFFSET.md](MASS_OFFSET.md) and the retained
  [evidence summary](scientific-results/mass-offset-20260919/summary.json).
- Regression: identical results to the beta.3 binary for configurations without
  offsets.
- Entrapment error, localization against synthesis-defined sites, determinism,
  chimeric and wide-window smoke tests: recorded in the evaluation.

## Pending before tagging

- [ ] Hosted required checks on the release preparation pull request.
- [ ] A manual `Release Sage Plus` workflow run and archive inspection.
- [ ] Dependency audit refresh for this release.
- [ ] Decide whether the Cascade upgrade lands alongside this tag.
