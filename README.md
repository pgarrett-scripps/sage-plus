<img src="figures/logo.png" width="300">

# Sage Plus

[![Rust](https://github.com/pgarrett-scripps/sage-plus/actions/workflows/rust.yml/badge.svg)](https://github.com/pgarrett-scripps/sage-plus/actions/workflows/rust.yml) [![Upstream Sage](https://img.shields.io/badge/upstream-lazear%2Fsage-blue)](https://github.com/lazear/sage)

Sage Plus is a fork of the [Sage](https://github.com/lazear/sage) proteomics search engine. It
keeps Sage's speed and configuration, and adds Thermo RAW input, PTM site localization with a
false-localization rate, ETD-aware scoring, large-database search in bounded memory, and better
LFQ. The current release is **v0.1.0-beta.17**, a prerelease.

## Sage Plus vs upstream Sage

Upstream Sage v0.15.0-beta.2 and Sage Plus run back to back on one machine with the same
configuration file ([full results](benchmarks/HEADTOHEAD-2026-10.md)).

| Test | Upstream Sage | Sage Plus | |
|---|---:|---:|---|
| Standard DDA, 4 Orbitrap runs (PXD028735) | 247,102 PSMs, 7.6 GiB | 250,060 PSMs, 6.2 GiB | +1.2% PSMs, 18% less memory, same speed |
| ETD spectra (PXD018176) | 4,515 PSMs | 6,799 PSMs | +51% with `z_dot` ions |
| Thermo RAW to results, 2 runs (PXD007058) | 118.6 s (89.6 s converting) | 44.5 s | 2.7x faster, no conversion step |
| False hits at a reported 1%, entrapment test | 1.64% PSMs, 2.91% protein groups | 1.25% PSMs, 2.02% protein groups | Fewer false hits at every level |

Also only in Sage Plus: a human reference plus a 100x gut-microbiome catalog searches in 19 GiB
(14,459 PSMs in 451 s, [Beta 11](CHANGELOG.md#v010-beta11---2026-09-26)), and phosphosites come
with a per-type false-localization rate (7,231 phospho localizations at 1% FLR on PXD007058,
[details](DOCS.md#ptm-site-localization)).

For a standard search the two engines are close. The gains are in the cases above. Entrapment
estimates are conservative for both engines; see the full results for the method.

## Install

Download a binary for Linux, macOS or Windows from the
[releases page](https://github.com/pgarrett-scripps/sage-plus/releases), or use the container:

```shell
docker run --rm -v "$PWD":/data ghcr.io/pgarrett-scripps/sage-plus:<tag> sage -o /data /data/config.json
```

From source (Rust 1.97+): `cargo build --release --workspace`. Standard builds include mzMLb and
S3/GCS/Azure paths; `--no-default-features` drops both.

## Quickstart

```shell
sage --write-config trypsin-hcd > config.json
sage -f proteins.fasta config.json run1.raw run2.mzML
```

> [!NOTE]
> Sage Plus is a beta, maintained independently of upstream Sage. Behavior can change between
> betas: pin a release, and validate results on your own data before publishing them.

## Learn more

- [All features](FEATURES.md), grouped by what they are for
- [Configuration and outputs](DOCS.md) and the [upstream Sage documentation](https://sage-docs.vercel.app/docs)
- [Changelog](CHANGELOG.md) and [benchmarks](benchmarks/README.md)
- [Upstream relationship](UPSTREAM.md) and [release procedure](RELEASING.md)

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