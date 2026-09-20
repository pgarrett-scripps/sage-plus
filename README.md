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

- Lower-memory searching with protein-backed target peptide sequences, compact spectrum storage,
  and a lossless six-byte theoretical-fragment index.
- Search-space memory estimation, runtime memory limits, minimum-free-memory protection, and configurable file batching.
- Per-modification limits, total variant caps, named modifications, and optional or required neutral-loss fragments.
- Search-time mass offset modifications: one offset placed per peptide without expanding the fragment index, competing with ordinary candidates and localized like any other modification.
- Modification-defined SILAC, dimethyl, and custom precursor channels on required static or optional variable modifications.
- Channel-aware LFQ with exact-mass partner extraction, reference ratios, and label-aware FDR.
- Robust per-file precursor and fragment mass-error alignment before final FDR rescoring.
- Configurable LFQ match-between-runs retention-time tolerance.
- Optional LFQ match-between-runs and independent ion-mobility model controls.
- Local mzMLb input in standard Sage builds.
- Typed protein occurrence coordinates in the canonical PSM Parquet output.
- A machine-readable JSON configuration schema, available with `sage --write-config-schema`.
- PTM localization, ambiguity-aware sequences, site reports, and target/decoy false-localization-rate q-values.
- Enriched linear retention-time features with regularized, cross-validated variable-PTM offsets.
- Optional robust nonlinear cross-run retention-time alignment shared by prediction and LFQ.
- PTM-aware, peptide-grouped ion-mobility prediction with cross-validated enriched features.
- Direct reading of local Thermo Fisher RAW files.
- Empirical spectral-library export in canonical Parquet and PSI mzSpecLib formats, using either a deterministic best PSM or robust consensus spectra.
- A structured runner API with JSONL events, validation-only mode, cancellation, and an automatic `run-summary.json` artifact.
- A root-bounded MCP server with persistent jobs and isolated search workers for configuration, estimation, safe execution, cancellation, monitoring, analysis, and result queries.

Most additions are opt-in, and upstream Sage defaults are retained where practical.

## Beta.6 named modifications and explicit sites

Define a modification once and list its complete attachment rules:

```json
"variable_mods": {
  "Acetyl": {
    "mass": 42.010565,
    "sites": ["first_residue:K", "internal_residue:K", "protein_last:K"],
    "max_count": 2
  },
  "Phospho": {"mass": 79.966331, "sites": ["S", "T", "Y"], "max_count": 3}
}
```

`peptide_n_term:K` modifies the terminal group when the peptide starts with K.
`first_residue:K` modifies the K residue itself. Typed site libraries preserve this
distinction through search, localization, and reuse. Ambiguous attachments are not
promoted into the reusable library. Static definitions use the same site vocabulary.

Convert older configurations with `sage old.json --migrate-modifications > new.json`.
Library-aware previews accept `--preview-protein` and `--preview-start`.
See [the modification guide](DOCS.md#modifications) and
[Beta 6 validation](benchmarks/BETA6_RELEASE.md). PTM libraries and site reports now
carry an attachment column, which requires updates in consumers expecting four columns.

## Beta.5 positional modifications

Use `~K` to restrict a modification to internal peptide residues. Combine `^K`
and `~K` to include the first residue, or `~K` and `$K` to include the last.
The same syntax works for static, indexed variable, and mass-offset modifications.
Entries sharing a name retain one occurrence limit across their combined sites.

Preview compatible sites and generated variants before searching:

```shell
sage config.json --preview-modifications KSTGGKAPR
```

Localization now preserves positional restrictions and distinguishes named
modifications with equal masses. Malformed modification keys fail configuration
loading. See [the positional modification guide](DOCS.md#positional-residue-modifications)
and [release validation](benchmarks/BETA5_RELEASE.md). Beta 5 was tagged but not published.
Its changes are included in Beta 6.

## Beta.4 mass offset search

A variable modification can set `"search_mode": "mass_offset"` to be searched as a
precursor and fragment offset instead of being expanded into the fragment index:

```jsonc
"variable_mods": {
  "S": [{"mass": 79.966331, "name": "Phospho", "search_mode": "mass_offset"}],
  "T": [{"mass": 79.966331, "name": "Phospho", "search_mode": "mass_offset"}],
  "Y": [{"mass": 79.966331, "name": "Phospho", "search_mode": "mass_offset"}],
  "M": [{"mass": 15.994915, "name": "Oxidation"}]
}
```

Each spectrum is searched once per offset against a translated precursor window, with
fragment lookups at both the unshifted and shifted masses, and every compatible placement
competes as its own candidate. At most one offset is placed on a peptide and offsets are
never combined with each other, so search cost grows linearly with the number configured
while the index stays at its unmodified size. Placements become ordinary peptidoforms
before FDR, quantification, localization, and site libraries, so an offset is never
reported as precursor mass error. See [DOCS.md](DOCS.md#mass-offset-modifications) for
behavior and limits and [the evaluation](benchmarks/MASS_OFFSET.md) for measured cost and
agreement with database expansion.

## Beta.3 hardening

Beta.3 fixes batching, gzip completion, event ordering, validation, and worker persistence.
Existing local Sage outputs now require explicit `--overwrite`, and run summaries use schema 9.
The release also adds verified benchmark provenance and resolves the audited dependency findings.
See the [release checklist and validation scope](benchmarks/BETA3_RELEASE.md) and
[changelog](CHANGELOG.md) for compatibility details and publication status.

## Memory and performance

Sage Plus `v0.1.0-beta.2` reduces the largest in-memory search structures without lossy mass
rounding:

- Target peptides reference immutable source-protein sequences by coordinates where possible,
  avoiding a separate sequence allocation for every generated peptide.
- The theoretical-fragment index uses lossless six-byte records while retaining bounded search
  buckets.
- Spectrum storage uses compact representations for fragment data and optional charge arrays.
- Exact database prefiltering can trade additional preparation time for a substantially smaller
  peptide and fragment search index.
- Memory estimation, configurable batching, and runtime memory guards help keep searches within
  available system memory.

On the local beta.2 release benchmark, conventional search used 29.5% less peak memory than
beta.1 and was 8.9% faster. A variable-modification search used 38.9% less peak memory and was
18.0% faster. Exact prefiltering reduced peak memory by 59.1% with byte-identical output, at an
81.5% wall-time cost for that workload. These measurements are machine and workload specific.
See the [benchmark methodology and complete results](benchmarks/RESULTS.md).

## Build and run

Sage Plus currently requires Rust 1.88 or newer.

```shell
git clone https://github.com/pgarrett-scripps/sage-plus.git
cd sage-plus
cargo build --release --workspace
./target/release/sage config.json
```

mzMLb support is included in standard builds and release binaries. Minimal source builds can omit
the HDF5-based mzMLb reader with `cargo build --release --workspace --no-default-features`.

The release build produces the standard `sage` executable and the optional `sage-mcp` server.
Run `sage --help` for CLI options.

Prebuilt binaries are available from [Sage Plus releases](https://github.com/pgarrett-scripps/sage-plus/releases),
and versioned Linux AMD64 container images are published as
`ghcr.io/pgarrett-scripps/sage-plus:<release-tag>`.

Sage Plus uses its own version sequence, beginning with `v0.1.0-beta.1`, independently of
upstream Sage releases.

## Documentation

- [Sage Plus configuration and outputs](DOCS.md)
- [Sage MCP server](crates/sage-mcp/README.md)
- [Maintainer release procedure](RELEASING.md)
- [Upstream relationship and synchronization](UPSTREAM.md)
- [Developer benchmark pipeline and results](benchmarks/RESULTS.md)
- [Upstream Sage documentation](https://sage-docs.vercel.app/docs)

## Attribution and citation

Sage Plus retains Sage's Git history, authorship, citation metadata, and MIT license. It is not
endorsed by or released on behalf of the upstream Sage maintainers. When publishing work that
uses Sage Plus, cite the original Sage paper listed in [CITATION.cff](CITATION.cff) and report the
exact Sage Plus commit or release used.
