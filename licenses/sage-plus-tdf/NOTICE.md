# Source attribution

This repository preserves the Git history of MannLabs TimsRust through commit `80e235a331057e0ef922eb3c784b7a631c9b9666` (0.6.6).

The TDF byte shuffling, scan offsets and TOF reconstruction in `src/compression.rs` derive from the upstream `crates/timsrust-tdf/src/frame_reader/compression1.rs` and `compression2.rs`. Bounds checks, resource limits and error paths were rewritten for the narrowed reader API.

Mobility equations and the reference fixture come from the local Sage Plus `crates/sage-cloudpath/src/tims_mobility.rs` and `tests/data/bruker/tims_calibration_sdk.json`. The latter records outputs of Bruker `tims_scannum_to_oneoverk0`, including fractional coordinates and boundary continuation. No Bruker SDK code or binaries are included.

Synthetic SQLite and compressed frame fixtures were generated locally for these regression tests. No private acquisition is included.

The original repository license notice and upstream package license declaration are both retained. Dependencies keep their own licenses.

The ModelType 1 mass equations were first adapted from MSConvert RS source calibration, which documents the timsrust-calibration 0.1.1 equations. The full model, including the `dC2` and `T2` drift and the `C2`, `C3` and `C4` terms, is checked against Bruker `tims_index_to_mz` and `tims_mz_to_index` outputs recorded in `tests/fixtures/mz_calibration_sdk.json`. That fixture holds only calibration coefficients, temperatures, coordinates and SDK outputs for the tdfpy v2.2.0 public test acquisitions described below. No Bruker SDK code or binaries are included.


## Source revisions and license interpretation

The MannLabs root `LICENSE` at `80e235a331057e0ef922eb3c784b7a631c9b9666` contains the MIT grant and Copyright 2026 MannLabs. At the same revision, `crates/timsrust-tdf/Cargo.toml` declares Apache-2.0. Both texts are retained. The package expression `Apache-2.0 AND MIT` preserves the obligations of the mixed source origins and does not assert that recipients can choose either license.

The Sage Plus source was inspected at revision `d7ab67d80ea01304f644b9e9335acd8d022fe7dd`. Its mobility adapter was last changed at `ec12bd8055dc5d77467d95610c93f72dd140d0c9`. The source crate declares MIT. `LICENSE-MIT-SAGE` retains the original Copyright (c) 2022 Michael Lazear notice from that repository.

The MSConvert RS source calibration adapter was inspected at `e6b9936d6137dd66b404e73430298406b95d45a5`, whose repository license is Apache-2.0. Its documented equation source is `timsrust-calibration` 0.1.1, whose manifest and bundled license both declare Apache-2.0. The applicable text is included as `LICENSE-APACHE`.

This derivative replaces upstream workspace and spectrum-processing APIs with read-only frame parsing, adds bounded decoding and metadata execution, and adds explicit source-calibration models and typed metadata tables. The package is renamed `sage-plus-tdf`. These modifications do not remove the notices of the retained source material.


## Calibration fixture provenance

The calibration reference JSON was generated through the documented Sage Plus adapter using metadata from the public `tacular-omics/tdfpy` v2.2.0 test acquisitions. The `analysis.tdf` SHA256 values recorded in the JSON match that tag's `tests/data/example_dda.d`, `example_dia.d` and `example_prm.d` inputs exactly. Source: [tdfpy v2.2.0 test acquisitions](https://github.com/tacular-omics/tdfpy/tree/v2.2.0/tests/data).

Only the calibration coefficients, scan coordinates, expected numeric outputs and source hashes are included here. The original public data source's MIT notice is retained in `LICENSE-MIT-TDFPY`. No vendor SDK code or binary is redistributed. Numerical reference values were obtained from Bruker `tims_scannum_to_oneoverk0` as documented in the fixture.
