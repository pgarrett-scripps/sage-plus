<img src="figures/logo.png" width="300">

# Sage Plus

[![Rust](https://github.com/pgarrett-scripps/sage-plus/actions/workflows/rust.yml/badge.svg)](https://github.com/pgarrett-scripps/sage-plus/actions/workflows/rust.yml) [![Upstream Sage](https://img.shields.io/badge/upstream-lazear%2Fsage-blue)](https://github.com/lazear/sage)

Sage Plus is a fork of the [Sage](https://github.com/lazear/sage) proteomics search engine. It
keeps Sage's speed and configuration, and adds Thermo RAW input, PTM site localization with
error rates, PTM searches without index blow-up, metaproteome-scale databases, ETD scoring, and
better label-free quantification. Current release: **v0.1.0-beta.17** (prerelease).

## At a glance

Compared with upstream Sage v0.15.0-beta.2. Numbers link to their source.

| | Sage | Sage Plus |
|---|:---:|---|
| **Input** | | |
| Thermo RAW, no conversion | ✗ | ✓ RAW to results [2.7x faster](benchmarks/HEADTOHEAD-2026-10.md) than converting first |
| mzMLb | ✗ | ✓ |
| timsTOF m/z and 1/K0 from the instrument calibration | linear | ✓ up to [0.054 1/K0](DOCS.md#spectrum-paths) closer |
| **PTMs** | | |
| Site localization with a false-localization rate | ✗ | ✓ [7,231 phosphosites](DOCS.md#ptm-site-localization) at 1% FLR |
| Site-level FDR and ProForma output | ✗ | ✓ [site q-values](DOCS.md#site-level-fdr), `localized_peptide` |
| Mods scored as mass offsets, not indexed | ✗ | ✓ phospho search in [1.8 vs 8.1 GiB](benchmarks/HEADTOHEAD-2026-10.md) |
| Named mods with explicit and motif sites | ✗ | ✓ e.g. `first_residue:K`, [`motif:N*-{P}-[ST]`](DOCS.md#motif-sites) |
| Per-site neutral losses, per-mod caps, site libraries | ✗ | ✓ losses add [+7.8% localizations](DOCS.md#static-and-variable-behavior) |
| **Scale** | | |
| Huge databases | chunked prefilter | ✓ streamed: human + 100x gut catalog in [19 GiB](DOCS.md#fasta) |
| Peak memory, standard DDA | 7.6 GiB | ✓ [6.2 GiB](benchmarks/HEADTOHEAD-2026-10.md) (−18%) |
| Memory limit and size estimate before a run | ✗ | ✓ [`max_memory_gb`, `--estimate`](DOCS.md#memory-guard) |
| **Identifications** | | |
| ETD/EThcD z• ions | ✗ | ✓ [+51% ETD PSMs](benchmarks/HEADTOHEAD-2026-10.md) |
| Initiator Met clipping | ✗ | ✓ [+3.4% PSMs](DOCS.md#initiator-methionine-clipping) with N-term acetylation |
| False hits at a reported 1% (entrapment) | 1.64% PSMs | ✓ [1.25% PSMs](benchmarks/HEADTOHEAD-2026-10.md), lower at every level |
| Mass recalibration, ambiguous residues (B/J/Z/X), I/L merge | ✗ | ✓ [details](DOCS.md#ambiguous-residues) |
| **Quantification** | | |
| Nonlinear RT alignment | linear | ✓ [+31% LFQ precursors](DOCS.md#other-settings) |
| SILAC and dimethyl channels in LFQ | ✗ | ✓ [channels](DOCS.md#modification-channels) |
| Per-row LFQ q-value, peak shape, MS2 confirmation | ✗ | ✓ [details](DOCS.md#label-free-quantification-output) |
| **Outputs** | | |
| Provenance in every Parquet file, versioned schemas | ✗ | ✓ [details](DOCS.md#output-provenance) |
| QC: digestion summary, polymer check, diagnostic ions | ✗ | ✓ [details](DOCS.md#quality-control-outputs) |
| Config schema, `--validate-only`, presets, progress events | ✗ | ✓ [details](DOCS.md#machine-readable-jobs) |
| Telemetry | sent by default | ✓ none |

Every feature, with settings and caveats, is in [FEATURES.md](FEATURES.md).

## Head-to-head

Same configuration, same machine, run back to back ([full results](benchmarks/HEADTOHEAD-2026-10.md)).

| Test | Upstream Sage | Sage Plus | |
|---|---:|---:|---|
| Standard DDA, 4 Orbitrap runs (PXD028735) | 247,102 PSMs, 7.6 GiB | 250,060 PSMs, 6.2 GiB | +1.2% PSMs, −18% memory, same speed |
| ETD spectra (PXD018176) | 4,515 PSMs | 6,799 PSMs | +51% with z• ions |
| Thermo RAW to results, 2 runs (PXD007058) | 118.6 s (89.6 s converting) | 44.5 s | 2.7x faster |
| Entrapment false hits at a reported 1% | 1.64% PSMs, 2.91% groups | 1.25% PSMs, 2.02% groups | lower at every level |

Where Sage Plus is not ahead: on the phospho dataset it finds 5.4% fewer peptides at the same PSM
count and is slower unless phospho is searched as a mass offset; on standard DDA the speed is the
same within noise. Entrapment estimates are conservative for both engines.

## Install

Binaries for Linux, macOS and Windows are on the
[releases page](https://github.com/pgarrett-scripps/sage-plus/releases). Container:

```shell
docker run --rm -v "$PWD":/data ghcr.io/pgarrett-scripps/sage-plus:<tag> sage -o /data /data/config.json
```

From source (Rust 1.97+): `cargo build --release --workspace`.

## Quickstart

```shell
sage --write-config trypsin-hcd > config.json
sage -f proteins.fasta config.json run1.raw run2.mzML
```

> [!NOTE]
> Sage Plus is a beta, maintained independently of upstream Sage. Behavior can change between
> betas: pin a release, and validate results on your own data before publishing them.

## Learn more

- [FEATURES.md](FEATURES.md): every addition, grouped by what it is for
- [DOCS.md](DOCS.md): configuration and outputs; [upstream Sage docs](https://sage-docs.vercel.app/docs)
- [CHANGELOG.md](CHANGELOG.md), [benchmarks](benchmarks/README.md), [upstream relationship](UPSTREAM.md)

## Attribution and citation

Sage Plus retains Sage's Git history, authorship, citation metadata, and MIT license. It is not
endorsed by or released on behalf of the upstream Sage maintainers. When publishing work that uses
Sage Plus, cite the original Sage paper listed in [CITATION.cff](CITATION.cff) and report the exact
Sage Plus commit or release used.

Sage Plus compiles in modification data from [Unimod](https://www.unimod.org). Its notice and
license are in [THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md) and
[`crates/sage/data/LICENSE-unimod.txt`](crates/sage/data/LICENSE-unimod.txt).

Initiator methionine clipping, ambiguous-residue expansion, the digestion summary, the polymer
contamination check, diagnostic-ion reporting, and tolerance recommendations were inspired by
[sageRecon](https://github.com/usnistgov/sageRecon) from NIST. Sage Plus implements them
independently; no sageRecon code is included.