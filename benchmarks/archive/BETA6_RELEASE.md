# Beta 6 release validation

Target: `v0.1.0-beta.6`, based on main after the unpublished Beta 5 tag.

## Scope

- Named static and variable definitions with explicit site strings.
- Residue-conditioned terminal-group attachments remain distinct from boundary residues.
- Typed library matching, localization, and TSV/Parquet export.
- Shared occurrence limits, legacy configuration migration, and library-aware preview.
- Conflicting fixed modifications fail validation.
- Legacy four-column libraries are interpreted as residue evidence. New outputs use
  five columns and preserve attachment identity. PTM analytical schemas are version 2.

No dependency versions are changed. The existing Beta 5 tag is retained unchanged.

## Reproducible regression evidence

```shell
cargo test --workspace --locked
cargo test --workspace --no-default-features --locked
cargo clippy --workspace --all-targets --all-features --locked -- -D warnings
cargo run --locked -p sage-core --example positional_modifications
python3 benchmarks/run_named_modifications.py --sage target/debug/sage --output benchmarks/NAMED_MODIFICATIONS_RESULTS.json
sage benchmarks/configs/named-modifications.json --preview-modifications KSTGGKAPR
```

The named-modification fixture runs first/last residue and N/C-terminal cases in
indexed and mass-offset modes. It verifies discovery output, guided reuse, and
iteration preserve the same typed library records. An indexed case carries both
terminal modifications together. Indistinguishable first-residue and N-terminal
alternatives produce no reusable library row, even at a permissive report threshold.

These are synthetic correctness regressions. They do not establish empirical
terminal-localization FLR calibration or production performance.

## Release gates

- [x] New parser, placement, shared-limit, typed-library, localization, and preview tests.
- [x] Full CLI discovery, guided reuse, and iteration regression (29 checks).
- [x] Complete workspace and minimal-feature tests.
- [x] Formatting, strict Clippy, and release version checks.
- [ ] Optimized build and existing fixture regression.
- [ ] Hosted required checks and preparation pull request merge.
- [ ] Manual packaging and native archive inspection.
- [ ] Annotated Beta 6 tag and verified publication with all assets.
