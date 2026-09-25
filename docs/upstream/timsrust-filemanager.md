# Draft upstream issue for MannLabs/timsrust

Status: draft only. Nothing has been opened upstream. Opening it is a maintainer decision.

---

**Title:** Bump filemanager dependencies (object_store, parquet/arrow) and let timsrust-core use `filemanager` without default features

**Body:**

Thanks for timsrust. We use it in a downstream search engine to read Bruker `.d`
(TDF) directories from local disk. Two dependency issues in the 0.6.6 release
chain make that harder than it needs to be.

### 1. Advisories reached through `filemanager` 0.6.6

`filemanager` 0.6.6 pins `object_store = "0.11"` and `arrow`/`parquet = "57"`,
`serde_arrow = "0.13"`. Through these, `cargo audit` reports:

- `quick-xml` 0.37 (pulled in by `object_store` 0.11's cloud features):
  RUSTSEC-2026-0194 and RUSTSEC-2026-0195.
- `thrift` < 0.23 (pulled in by `parquet` < 58): the excessive-allocation advisory.

The source changes needed to move to `object_store` 0.14 and `arrow`/`parquet` 59
are small. We carry them as a local `[patch.crates-io]` and all our tests pass:

- In `Cargo.toml`: `object_store = "0.14.1"`, `arrow = "59"`, `parquet = "59"`,
  `serde_arrow = { version = "0.15", features = ["arrow-59"] }`.
- `object_store` 0.14 moved the convenience methods (`get`, `put`, `head`,
  `get_range`, ...) to the `ObjectStoreExt` trait, so `use object_store::ObjectStoreExt;`
  is needed where the store is used.
- `get_range` and friends now take `Range<u64>`, so byte ranges are converted from
  `usize` (plus a checked `usize::try_from` on object sizes).
- `WriterProperties::set_max_row_group_size` is deprecated in parquet 59, and
  `set_max_row_group_row_count` replaces it.

We are happy to open a PR with exactly these changes if that helps.

### 2. `timsrust-core` forces all `filemanager` default features

`crates/timsrust-core/Cargo.toml` declares:

```toml
filemanager = { workspace = true, default-features = true, optional = true }
```

This overrides the workspace entry's `default-features = false`. It always turns on
`filemanager`'s `cloud` feature (object_store with AWS, GCP and Azure clients, and
tokio) and its `json` feature, even for users who only read local files. Downstream
crates cannot turn these off, because Cargo features are additive.

For reading local `.d` directories, the crates use `filemanager::Uri`,
`formats::sql::SqlReader`/`SqlError` and `formats::binary::BinaryReader`/`BinaryError`.
`timsrust-minitdf` and `timsrust-parquet-spectra` also use `formats::parquet::ParquetReader`.
So the minimum is `sql` plus `parquet`. `cloud` and `json` are not needed.

Suggested change:

```toml
# crates/timsrust-core/Cargo.toml
filemanager = { workspace = true, default-features = false, features = ["sql", "parquet"], optional = true }

[features]
io = ["filemanager"]
cloud = ["filemanager?/cloud"]
```

The `timsrust` crate could then forward a `cloud` feature (on by default if you want
to keep today's behaviour). Local-only users could build without object_store, and
optionally without arrow if `parquet` were also made a feature of the minitdf and
parquet-spectra readers.

The `filemanager` build without `cloud` does not currently compile. There are three problems:

- `cloud_store.rs` has `pub use no_cloud::CloudObject`, but the struct is `pub(crate)` (E0365).
- `no_cloud.rs` imports `crate::CloudProvider`, which is not exported at the crate root.
- `formats::binary` needs `CloudObject: Debug`.

We fixed these with `pub(crate) use`, `use crate::cloud_store::{CloudError, CloudProvider}`,
and `#[derive(Debug)]`. We also made the no-cloud `CloudProvider::parse` recognize
`s3`/`gs`/`az` schemes, so those URIs return `CloudError::FeatureNotEnabled` instead of
being treated as local paths. Happy to include this in the same PR.

### Environment

- timsrust / timsrust-core / filemanager 0.6.6 (crates.io, matching commit `80e235a`)
- Rust 1.88+
