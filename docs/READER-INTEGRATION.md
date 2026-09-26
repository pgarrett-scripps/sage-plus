# Native reader integration

Thermo input uses the strict `opentfraw-reader` API through the existing `opentfraw` import alias. MSConvert uses full-precision profile arrays. Sage propagates extraction errors instead of returning a partial acquisition. Existing consumer acquisition gates and scientifically justified processing policies remain in place. Sage retains its tested trailer-level safeguard for Fusion event alignment because the strict library API does not independently establish correct alignment for every inherited event layout.

Bruker raw frame IO and metadata use `timsrust-reader`. The consumer-owned `vendor/timsrust-tdf` adapter preserves the public upstream 0.6.6 processing algorithms, including DDA precursor merging, DIA windows, smoothing and centroiding. See its MIGRATION.md for source provenance and the narrow adapter changes. No private library source is vendored. Source calibration equations, approximation flags and Sage application denoising remain consumer policy. Sage aborts imports on frame decoding failures.

The Rust SQLite dependency raises the required compiler to Rust 1.97. The normal build graph may retain C dependencies for other formats and application processing. This is not a claim that the entire consumer is Rust-only.

## Private dependency access

Final Cargo manifests pin reviewed Git revisions. Building requires read access to both `pgarrett-scripps/opentfraw` and `pgarrett-scripps/timsrust-reader`. Configure a credential helper locally. CI needs a trusted credential restricted to reading those two repositories, and Cargo Git CLI fetching. The ordinary repository GITHUB_TOKEN does not grant cross-repository private access. Do not put tokens in manifests, URLs, source files or logs.

Sage Plus is public. Keep this integration local while the dependencies are private. Public pull requests and release jobs must not be sent private credentials or private reader source. MSConvert RS is private, but its CI still needs explicit access to the reader repositories.

## Validation

Run the workspace tests and Clippy on Rust 1.97.1. The new guard regressions prove the active processing adapter invokes the native reader and corrupted frames fail the entire import. Compare real input spectra and writer outputs against a frozen pre-migration checkout. Explain intentional corrected Thermo metadata differences separately from numerical changes. Qualification evidence is recorded by the integration run, not inferred from successful compilation.

## Expected metadata corrections

The Angiotensin qualification fixture has v66 reaction tag 1, raw reaction energy 35 and an HCD Energy trailer value of -1. The latter is an unset sentinel. The former fields do not establish an independently qualified activation method or energy unit. The previous reader inferred HCD and energy from those ambiguous fields. The strict reader leaves them unknown. Peak arrays are expected to remain exact, while activation metadata and Sage acquisition grouping may change. These omissions are reported in MSConvert warnings and Sage logs.

The trimmed HCD fixture also loses a dangling precursor link to scan 3163 outside its acquisition range of 1 through 1. A missing source parent cannot be represented as a valid internal spectrum reference.

## Recorded adoption qualification

See `READER-ADOPTION-QUALIFICATION.json` for baseline and reader revisions, exact numerical comparisons, intentional metadata corrections and qualification limits. The three Bruker conversion cases match the prior implementation at zero numerical tolerance across 13,090 exported spectra. The two Thermo acquisitions preserve 5,455 scans in both profile and centroid modes. Sage preserves complete processed spectra for the Bruker DDA and DIA cases, while the larger Thermo acquisition now reads successfully. These checks qualify conversion and search input. They do not claim peptide search result equivalence.
