# Sage Plus dependency patch

This directory contains the published `filemanager` 0.6.5 source from
MannLabs/timsrust commit `5e09572fa3f1aace86e6ba64244a01fd0850fbf3`, subdirectory
`crates/filemanager`. The upstream Apache-2.0 license is retained in `LICENSE`.
`UPSTREAM_SHA256.json` records the original copied files before patching.

Sage Plus changes the normalized manifest to use `object_store` 0.14.1, imports
`ObjectStoreExt`, converts storage byte ranges to `u64`, and checks that object
sizes fit `usize`. Redundant development dependency entries are removed so Cargo
can run the added tests from the root lockfile without a second workspace. This removes the old `quick-xml` 0.37 dependency and its
RUSTSEC-2026-0194 and RUSTSEC-2026-0195 findings. No advisory is suppressed.

The added tests cover upload, metadata, ranges, listing, buffered file upload,
and download through the storage API. A local HTTP fixture also exercises S3
XML listing and HTTP range responses without external credentials. These tests
do not certify every provider's authentication. Workspace tests separately exercise
the committed real Bruker TDF fixture.

Run the patch tests from the repository root:

```shell
cargo test -p filemanager --locked
```

Remove this directory and the root Cargo patch when a published timsrust and
filemanager dependency chain uses a patched XML parser. Re-run these compatibility
checks, the complete workspace gates, and the dependency audit before removal.

The untouched `Cargo.toml.orig` is retained for provenance. Cargo builds from the
patched normalized `Cargo.toml`.
