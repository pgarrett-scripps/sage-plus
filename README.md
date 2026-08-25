<img src="figures/logo.png" width="300">

# Sage Plus

> [!WARNING]
> **Sage Plus is an experimental, actively developed downstream distribution.** Features are
> being integrated rapidly and have not yet been comprehensively validated together. APIs,
> configuration, output formats, and scientific behavior may change or contain unresolved
> issues. Sage Plus is independently maintained, is not an official Sage release, and should
> not be adopted accidentally as a drop-in replacement. Most users should use
> [upstream Sage](https://github.com/lazear/sage). If you evaluate Sage Plus, pin an exact commit
> or release and independently validate results before relying on them for production, clinical,
> or published analyses.

[![Rust](https://github.com/pgarrett-scripps/sage-plus/actions/workflows/rust.yml/badge.svg)](https://github.com/pgarrett-scripps/sage-plus/actions/workflows/rust.yml) [![Upstream Sage](https://img.shields.io/badge/upstream-lazear%2Fsage-blue)](https://github.com/lazear/sage)

Sage Plus is a downstream distribution of the
[Sage proteomics search engine](https://github.com/lazear/sage). It preserves Sage's core
workflow while integrating experimental PTM, modeling, performance, automation, and
agent-facing capabilities.

## Sage Plus additions

- Faster, lower-memory searching with compact peptide/spectrum storage and optional bitmap preliminary scoring.
- Search-space memory estimation, runtime memory limits, minimum-free-memory protection, and configurable file batching.
- Per-modification limits, total variant caps, named modifications, and optional or required neutral-loss fragments.
- Modification-defined SILAC, dimethyl, and custom precursor channels on required static or optional variable modifications.
- Channel-aware LFQ with exact-mass partner extraction, reference ratios, label-aware FDR, and spectral-library round trips.
- Robust per-file precursor and fragment mass-error alignment before final FDR rescoring.
- Configurable LFQ match-between-runs retention-time tolerance.
- PTM localization, ambiguity-aware sequences, site reports, and target/decoy false-localization-rate q-values.
- Enriched linear retention-time features with regularized, cross-validated variable-PTM offsets.
- Optional robust nonlinear cross-run retention-time alignment shared by prediction and LFQ.
- PTM-aware, peptide-grouped ion-mobility prediction with cross-validated enriched features.
- Direct reading of local Thermo Fisher RAW files.
- Empirical spectral-library export in canonical Parquet and PSI mzSpecLib formats, using either a deterministic best PSM or robust consensus spectra.
- Standalone DDA library search with shuffled decoys, isotope-aware matching, per-file mass calibration, RT and mobility alignment, matched-fragment export, LFQ, and TMT.
- Regularized library rescoring with spectrum, peptide, protein, and protein-group FDR reporting.
- A structured runner API with JSONL events, validation-only mode, cancellation, and an automatic `run-summary.json` artifact.
- A root-bounded MCP server with persistent jobs and isolated search workers for configuration, estimation, safe execution, cancellation, monitoring, analysis, and result queries.

Most additions are opt-in, and upstream Sage defaults are retained where practical.

## Build and run

Sage Plus currently requires Rust 1.88 or newer.

```shell
git clone https://github.com/pgarrett-scripps/sage-plus.git
cd sage-plus
cargo build --release --workspace
./target/release/sage config.json
```

The release build produces the standard `sage` executable and the optional `sage-mcp` server.
Run `sage --help` for CLI options.

Prebuilt binaries are available from [Sage Plus releases](https://github.com/pgarrett-scripps/sage-plus/releases),
and versioned Linux AMD64 container images are published as
`ghcr.io/pgarrett-scripps/sage-plus:<release-tag>`.

## Documentation

- [Sage Plus configuration and outputs](DOCS.md)
- [Sage MCP server](crates/sage-mcp/README.md)
- [Maintainer release procedure](RELEASING.md)
- [Upstream relationship and synchronization](UPSTREAM.md)
- [DDA spectral-library search](DDA_LIBRARY_SEARCH.md)
- [Upstream Sage documentation](https://sage-docs.vercel.app/docs)

## Attribution and citation

Sage Plus retains Sage's Git history, authorship, citation metadata, and MIT license. It is not
endorsed by or released on behalf of the upstream Sage maintainers. When publishing work that
uses Sage Plus, cite the original Sage paper listed in [CITATION.cff](CITATION.cff) and report the
exact Sage Plus commit or release used.
