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

## Differences from upstream Sage

The tables compare Sage Plus with upstream Sage
[`v0.15.0-beta.2`](https://github.com/lazear/sage/releases/tag/v0.15.0-beta.2), the latest
published Sage release. Each row lists a feature present in the current release and the
published Sage Plus release that first shipped it in its current form. Upstream capabilities
such as protein grouping, cloud storage, HTML reports, and Percolator output are not repeated
here. Most additions are opt-in, and upstream defaults are retained where practical. Features
that were later removed or replaced, such as the DDA spectral-library search mode and Beta 5
positional keys, are recorded only in the [changelog](CHANGELOG.md).

Measured benefits come from the linked benchmarks. They are workload and machine specific.
Other benefits describe the intended effect and have not all been validated independently.

### Search performance and memory

| Feature | Since | Why it was added | Benefit |
|---|---|---|---|
| Spectrum-indexed exact prefilter | beta.7 | The chunked prefilter rebuilt a fragment index and reread spectra for every database chunk | Same retained peptides and results; 1.2 to 41 times faster than Beta 6 prefiltering ([results](benchmarks/PREFILTER.md)) |
| Protein-backed peptide sequences and compact modification records | beta.2 | Every generated peptide allocated its own sequence and a dense modification vector | 29.5% less peak memory on a conventional search and 38.9% less with variable modifications ([results](benchmarks/RESULTS.md)) |
| Lossless six-byte fragment index | beta.2 | Fragment records dominate index memory in large searches | Smaller index with exact masses and bounded search buckets |
| Compact spectrum storage | beta.2 | Loaded spectra repeat fragment and charge data | Lower resident memory for large file batches |
| Memory estimation, `max_memory_gb`, and `min_free_memory_gb` | beta.1 | Large searches could exhaust workstation memory | Oversized searches are rejected before they start, and running searches stop before the limit |
| Configurable `batch_size` | beta.1 | Only a command-line option controlled file batching | Batching can be set per configuration |

### Modifications and PTMs

| Feature | Since | Why it was added | Benefit |
|---|---|---|---|
| Motif modification sites (`motif:N*-{P}-[ST]`) | beta.8 | Residue sites could not require a sequence context such as the N-glycosylation sequon or a kinase motif | PROSITE-style patterns are evaluated against the source protein, including residues beyond the peptide, with mirrored decoys and motif-restricted localization |
| Named modifications with explicit sites | beta.6 | Residue keys could not separate a terminal group from the residue at that terminus, or exclude terminal residues | One definition and one occurrence limit across attachment rules such as `first_residue:K`, `internal_residue:K`, and `peptide_n_term`; `--migrate-modifications` converts older configurations |
| Modification preview (`--preview-modifications`) | beta.6 | Placement rules could only be checked by running a search | Eligible sites and generated variants for a peptide, optionally in protein context, without loading spectra |
| Typed terminal-group localization and version 2 PTM libraries | beta.6 | Libraries recorded residues only | Terminal and residue attachments stay distinct through search, localization, and reuse |
| Mass-offset modifications | beta.4 | Every variable modification multiplies the fragment index | The index keeps its unmodified size; phosphorylation search used 0.16 GB instead of 0.45 GB with the same PSMs ([evaluation](benchmarks/MASS_OFFSET.md)) |
| Per-modification limits, variant caps, and neutral-loss fragments | beta.1 | Combinatorial expansion was only bounded globally | Bounded search spaces and neutral-loss fragment matching |
| Separate PEFF modification budget (`max_peff_variable_mods`) | beta.1 | PEFF and global modifications shared one budget | Independent limits for annotated and global modifications |
| PTM localization with false-localization-rate q-values | beta.1 | Search reports a peptidoform without site confidence | Site probabilities, site reports, and target/decoy localization confidence |
| Ambiguity-aware sequences and residual mass shifts | beta.1 | Unsupported residue orders were reported as certain | Sequence regions without fragment evidence are marked in the output |

### Scoring and modeling

| Feature | Since | Why it was added | Benefit |
|---|---|---|---|
| Count-based confidence fallback | beta.4 | The density model can fail on small or unusual score distributions | Peptide and protein q-values remain defined |
| Averagine-scored isotope envelopes and charge-aware fragment matching | beta.2 | Deisotoping could assign peaks to several envelopes | Each peak belongs to one envelope, and fragment charges constrain matching |
| Per-file precursor and fragment mass-error alignment | beta.1 | Systematic mass error differs between files | Mass errors are corrected before final rescoring |
| Enriched retention-time model and nonlinear cross-run alignment | beta.1 | A linear model misses modification effects and nonlinear drift | Better retention-time features for rescoring and LFQ |
| PTM-aware ion-mobility prediction | beta.1 | Mobility features ignored modifications | Cross-validated mobility features for rescoring |

### Quantification and libraries

| Feature | Since | Why it was added | Benefit |
|---|---|---|---|
| LFQ confirmation and per-file signal diagnostics | beta.4 | Integrated signals had no quality evidence | MS2 confirmation, spectral angle, trace cosine, and retention shift per file |
| Match-between-runs and ion-mobility model switches | beta.2 | Neither could be disabled independently | Controlled quantification and modeling experiments |
| Modification-defined SILAC, dimethyl, and custom label channels | beta.1 | Labels were not tied to modifications | Coherent precursor channels with channel-aware LFQ, reference ratios, and label-aware FDR |
| Configurable match-between-runs retention tolerance | beta.1 | The transfer window was fixed | Tolerance can match the chromatography |
| Empirical spectral-library export | beta.1 | Search results could not be reused as libraries | Parquet and PSI mzSpecLib libraries from best or consensus spectra |

### Inputs and outputs

| Feature | Since | Why it was added | Benefit |
|---|---|---|---|
| Calibrated timsTOF ion mobility | beta.7 | timsrust interpolates 1/K0 between the acquisition limits instead of applying the instrument calibration | Reported 1/K0 equals the Bruker SDK value; the old scale was off by up to 0.054 1/K0 on a PXD070049 run, with nearly unchanged identifications ([validation](benchmarks/BETA7_RELEASE.md)). `bruker_config.ion_mobility_scale: "linear"` restores the old scale |
| mzMLb input | beta.2 | Compressed HDF5 spectra required conversion | Read directly in standard builds |
| Typed protein occurrences in Parquet output | beta.2 | Protein positions required re-mapping | One-based coordinates and flanking residues for each protein |
| JSON configuration schema | beta.2 | Configuration errors surfaced only at run time | Editor completion and static validation (`sage --write-config-schema`) |
| Thermo RAW input | beta.1 | RAW files required conversion | Local RAW files are read directly |
| Pre-digested peptide and custom cleavage-site inputs | beta.1 | Only FASTA digestion was supported | Peptide lists and protein-specific cleavage sites extend the database |
| Parquet as the canonical analytical output | beta.1 | TSV and Parquet outputs diverged | One typed output with nulls for missing signals |

### Automation and safety

| Feature | Since | Why it was added | Benefit |
|---|---|---|---|
| Overwrite protection for existing outputs | beta.3 | Reruns could silently replace results | Replacing Sage outputs requires `--overwrite` |
| Runner API, JSONL events, and `run-summary.json` | beta.1 | Runs could only be followed through logs | Validation-only runs, progress events, cancellation, and a machine-readable summary |
| MCP server with isolated search workers | beta.1 | Agents needed safe, persistent search jobs | A failed or cancelled search affects only its own worker |

## Prefilter performance

Database prefiltering keeps every peptide that can match a fragment in any spectrum, then builds
the search index from those peptides only. It gives the same results as a full search. Beta 7
indexes the spectra once and streams the generated peptides through that index.

| Workload | Beta 6 prefilter | Beta 7 prefilter | No prefilter | Peptides kept |
|---|---:|---:|---:|---:|
| HEK, oxidation and acetylation | 15.6 s | 11.6 s | 11.5 s | 37% |
| HEK, seven variable modifications | 33.0 s | 24.1 s | 27.1 s | 36% |
| HEK, seven modifications, two per peptide | 105.0 s | 67.4 s | exceeds memory | 30% |
| HEK, open search | 65.8 s | 54.6 s | 47.7 s | 99.8% |
| Five LFQ files, 567,401 spectra | 613.4 s | 135.2 s | 102.9 s | 90% |
| Phosphorylation mass offset, HCD_1 | 35.7 s | 0.9 s | 0.8 s | 8% |

Median of three runs on a 16-thread workstation. Prefiltering pays off when it removes much of
the database: the first two searches used 60 to 66% less memory than unfiltered searches, and the
two-modification search only fits with it. Open searches and multi-file runs that keep most
peptides are faster without it. See the [prefilter benchmark](benchmarks/PREFILTER.md) for
memory, output agreement, and repeats.

## Roadmap

These additions are planned. None of them is available yet.

| Addition | Goal | Notes |
|---|---|---|
| dnoise MS1 denoising for Bruker timsTOF | Keep ion-mobility streaks and drop MS1 noise before feature extraction, which evaluations show improves quantified coverage. | Optional, MS1 only, off by default. Targets dnoise 0.5.0, which shares Sage Plus's timsrust and SQLite versions. |
| koth feature finding | Extract MS1 and fragment ion chromatograms (hills) with retention time and ion mobility. | Replaces or complements the current MS1 feature extraction for LFQ, and provides the signals for the DIA mode below. |
| Chromatogram-based DIA mode | Score precursor and fragment chromatograms together, possibly in a peptide-centric search, instead of treating DIA scans as wide-window spectra. | Builds on koth feature finding, retention-time and mobility models, and spectral libraries. |
| Custom residue compositions and isotope-labeled residues | Define residues by elemental composition, and label whole residues or backbones with heavy isotopes (for example <sup>15</sup>N or <sup>13</sup>C metabolic labeling, or deuterium). | Label mass shifts then depend on the sequence rather than on a fixed modification mass. |
| Glycopeptide search | Search glycan compositions on sequon or motif-restricted sites with oxonium and Y-ion evidence. | Builds on mass offsets, neutral losses, and motif modification sites. |
| Crosslink search | Identify crosslinked peptide pairs, starting with simple or MS-cleavable linkers. | Needs pair-aware candidate generation and crosslink-specific FDR. |

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
- [Release validation records](benchmarks/BETA7_RELEASE.md) and the [changelog](CHANGELOG.md)
- [Upstream Sage documentation](https://sage-docs.vercel.app/docs)

## Attribution and citation

Sage Plus retains Sage's Git history, authorship, citation metadata, and MIT license. It is not
endorsed by or released on behalf of the upstream Sage maintainers. When publishing work that
uses Sage Plus, cite the original Sage paper listed in [CITATION.cff](CITATION.cff) and report the
exact Sage Plus commit or release used.
