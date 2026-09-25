# Sage Plus dependency patch

This directory contains the published `filemanager` 0.6.6 source from
MannLabs/timsrust commit `80e235a331057e0ef922eb3c784b7a631c9b9666`, subdirectory
`crates/filemanager`. The crate metadata declares Apache-2.0, while the pinned
repository root contains an MIT license. The upstream copyright and MIT text are
retained unchanged in `LICENSE`. `LICENSE-APACHE` contains the standard Apache 2.0
text matching the crate declaration. Both are distributed with executable artifacts.
`UPSTREAM_SHA256.json` records the original copied files before patching.

Sage Plus changes the normalized manifest to use `object_store` 0.14.1, imports
`ObjectStoreExt`, converts storage byte ranges to `u64`, and checks that object
sizes fit `usize`. Redundant development dependency entries are removed so Cargo
can run the added tests from the root lockfile without a second workspace. This removes the old `quick-xml` 0.37 dependency and its
RUSTSEC-2026-0194 and RUSTSEC-2026-0195 findings. It also moves `arrow` and `parquet` from 57
to 59 and `serde_arrow` from 0.13 to 0.15 (`arrow-59`), and replaces the deprecated
`set_max_row_group_size` with `set_max_row_group_row_count`. Parquet 58+ no longer depends on
the `thrift` crate, which removes the thrift < 0.23.0 excessive-allocation advisory. No
advisory is suppressed.

Sage Plus also removes `cloud` from the default features. `timsrust-core` depends on
`filemanager` with `default-features = true`, so the default feature set is what every Sage
Plus build compiles. Sage Plus reads Bruker `.d` directories only from local paths (see
`read_tdf` in `sage-cloudpath`), and timsrust needs only the `sql` and `parquet` features
(the binary reader and `Uri` are always built). With `cloud` off, no Sage Plus build compiles filemanager's
`object_store` client. Sage Plus's own S3/GCS/Azure support is the separate `cloud`
feature of `sage-cloudpath`. The upstream no-`cloud` path did not compile, so
`src/cloud_store.rs` re-exports `CloudObject` as `pub(crate)`. `src/cloud_store/no_cloud.rs`
also imports `CloudProvider` from its module, derives `Debug`, and recognizes S3, GCS,
and Azure schemes, so those URIs fail with `FeatureNotEnabled` instead of being read as
local paths. `parquet` stays enabled because `timsrust-minitdf` and
`timsrust-parquet-spectra` use `formats::parquet::ParquetReader` unconditionally.

The added tests cover upload, metadata, ranges, listing, buffered file upload,
and download through the storage API. A local HTTP fixture also exercises S3
XML listing and HTTP range responses without external credentials. These tests
do not certify every provider's authentication. They need the `cloud` feature, so they
compile out of the shipped feature set that the command below tests. Workspace tests separately exercise
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
