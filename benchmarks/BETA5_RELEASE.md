# Beta 5 release checklist

Target: `v0.1.0-beta.5`. Prepared September 20, 2026 against the Beta 4 baseline.

## Scope

- `~K` and equivalent residue keys restrict static, indexed variable, and mass-offset
  modifications to internal peptide residues.
- Combined positional keys retain shared named occurrence limits and deduplicate
  overlapping placements. Protein-boundary rules remain explicit.
- Modification preview reports compatible sites, limits, and bounded generated
  variants without running a search.
- Invalid modification keys fail loading. The committed JSON schema uses the same
  accepted key syntax.
- Localization retains full definitions and positional rules. Equal-mass named
  modifications remain separate and differently modified occupied sites remain fixed.

## Reproduce positional validation

```shell
cargo run --locked -p sage-core --example positional_modifications
cargo run --locked -p sage-cli --bin sage -- benchmarks/configs/positional-modifications.json --preview-modifications KSTGGKAPR
cargo run --locked -p sage-cli --bin sage -- benchmarks/configs/positional-modifications.json --preview-modifications ASQKSTGGK --peptide-position cterm
```

The eight synthetic scoring and localization cases pass in both indexed and
offset modes. Each mode recovers the specified peptide and site with 16 matched
peaks. Retained results are in
[POSITIONAL_MODIFICATIONS_RESULTS.json](POSITIONAL_MODIFICATIONS_RESULTS.json).
These are regression fixtures, not empirical calibration or a performance claim.

Core tests additionally cover lengths one through three, shared occurrence limits,
overlapping keys, PTM-library restrictions, label channels, conservative memory
counts, retention and mobility feature counting, and equal-mass identities.

The VAT1 and original synthetic search fixtures produce byte-identical PSM Parquet
outputs to the published Beta 4 Linux GNU executable, with two Rayon threads.
VAT1 matched-fragment output is also byte-identical. Recorded hashes are in
[BETA5_REGRESSION_RESULTS.json](BETA5_REGRESSION_RESULTS.json).

## Compatibility

- Existing valid modification keys keep their placement meaning. Malformed keys
  now produce an error, including keys that previously lost trailing characters.
- Residue-qualified terminal rules participate in positional localization. Bare
  terminal-group modifications remain fixed.
- Core callers can continue passing mass-only localization rules. Full identity
  protection requires the new `IndexedDatabase::localization_mods` definitions.
- Programmatically constructed `Builder` values can call
  `validate_modification_keys` to receive a recoverable error before conversion.
  Unchecked invalid maps passed directly to conversion fail instead of dropping
  entries.
- Preview is bounded and supports exhaustive placement only. It explicitly rejects
  PTM-library configurations and does not infer protein context from input files.
- No analytical output schema version changes or dependency upgrades are required.

## Release gates

- [x] Synthetic positional benchmark.
- [x] Schema syntax checked against valid and malformed modification keys.
- [x] Final workspace and minimal-feature tests.
- [x] Strict Clippy and formatting.
- [ ] Optimized release build and storage patch tests.
- [x] Existing-fixture regression against the published Beta 4 executable.
- [ ] Hosted required checks and release preparation pull request merge.
- [ ] Manual release packaging and native archive inspection.
- [ ] Annotated Beta 5 tag and successful publication with all assets.
