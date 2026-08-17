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

This is my personal experimental project for exploring how Sage, proteomics models, structured
workflows, and related analysis tools can be integrated into automated and AI-assisted pipelines.
I am sharing it publicly for transparency, reproducibility, and collaboration while I test ideas;
it is not intended to compete with, replace, or speak for the upstream Sage project.

[![Rust](https://github.com/pgarrett-scripps/sage-plus/actions/workflows/rust.yml/badge.svg)](https://github.com/pgarrett-scripps/sage-plus/actions/workflows/rust.yml) [![Upstream Sage](https://img.shields.io/badge/upstream-lazear%2Fsage-blue)](https://github.com/lazear/sage)

Sage Plus is a downstream distribution of the [Sage proteomics search engine](https://github.com/lazear/sage). It preserves Sage's core workflow while integrating additional PTM, modeling, performance, automation, and agent-facing capabilities. It is independently maintained and is not an official Sage release.

See [DOCS.md](DOCS.md) for Sage Plus configuration and the [upstream Sage documentation](https://sage-docs.vercel.app/docs) for the core search engine.

## Sage Plus additions

This distribution adds the following features on top of upstream Sage:

- Faster, lower-memory searching with compact peptide/spectrum storage and optional bitmap preliminary scoring.
- Search-space memory estimation, runtime memory limits, minimum-free-memory protection, and configurable file batching.
- Per-modification limits, total variant caps, named modifications, and optional or required neutral-loss fragments.
- Robust per-file precursor and fragment mass-error alignment before final FDR rescoring.
- Configurable LFQ match-between-runs retention-time tolerance.
- PTM localization, ambiguity-aware sequences, site reports, and target/decoy false-localization-rate q-values.
- Enriched linear retention-time features with regularized, cross-validated variable-PTM offsets.
- Optional robust nonlinear cross-run retention-time alignment shared by prediction and LFQ.
- PTM-aware, peptide-grouped ion-mobility prediction with cross-validated enriched features.
- Direct reading of local Thermo Fisher RAW files.
- A structured runner API with JSONL events, validation-only mode, cancellation, and an automatic `run-summary.json` artifact.
- A root-bounded MCP server for AI-assisted configuration, estimation, execution, monitoring, and result queries.

Most additions are opt-in; existing Sage defaults are retained where practical. See [UPSTREAM.md](UPSTREAM.md) for the maintenance and attribution policy.


# Introduction
 
Sage is, at it's core, a proteomics database search engine - 
    a tool that transforms raw mass spectra from proteomics experiments into peptide identifications 
    via database searching & spectral matching. 

However, Sage includes a variety of advanced features that make it a one-stop shop: retention time prediction, quantification (both isobaric & LFQ), peptide-spectrum match rescoring, and FDR control. You can directly use results from Sage without needing to use other tools for these tasks.

Additionally, Sage was designed with cloud computing in mind - massively parallel processing and the ability to directly stream compressed mass spectrometry data to/from AWS S3 enables unprecedented search speeds with minimal cost. 

 Sage also runs just as well reading local files from your Mac/PC/Linux device!

## Why use Sage instead of other tools?

Sage is **simple to configure**, **powerful** and **flexible**. 
It also happens to be well-tested, **mind-boggingly fast**, open-source (MIT-licensed) and free.

## Citation

If you use Sage in a scientific publication, please cite the following paper:

[Sage: An Open-Source Tool for Fast Proteomics Searching and Quantification at Scale](https://doi.org/10.1021/acs.jproteome.3c00486)


## Features

- Incredible performance out of the box
- [Effortlessly cross-platform](https://sage-docs.vercel.app/docs/started#download-the-latest-binary-release) (Linux/MacOS/Windows), effortlessly parallel (uses all of your CPU cores)
- [Fragment indexing strategy](https://sage-docs.vercel.app/docs/how_it_works) allows for blazing fast narrow and open searches (> 500 Da precursor tolerance)
- [Isobaric quantification](https://sage-docs.vercel.app/docs/how_it_works#tmt-based) (MS2/MS3-TMT, or custom reporter ions)
- [Label-free quantification](https://sage-docs.vercel.app/docs/how_it_works#label-free): consider all charge states & isotopologues *a la* FlashLFQ
- Capable of searching for [chimeric/co-fragmenting spectra](https://sage-docs.vercel.app/docs/configuration/additional)
- Wide-window (dynamic precursor tolerance) search mode - [enables WWA/PRM/DIA searches](https://sage-docs.vercel.app/docs/configuration/tolerance#wide-window-mode)
- Retention time prediction models fit to each LC/MS run
- [PSM rescoring](https://sage-docs.vercel.app/docs/how_it_works#machine-learning-for-psm-rescoring) using built-in linear discriminant analysis (LDA)
- PEP calculation using a non-parametric model (KDE)
- FDR calculation using target-decoy competition and picked-peptide & picked-protein approaches
- Percolator/Mokapot [compatible output](https://sage-docs.vercel.app/docs/configuration#env)
- Configuration by [JSON file](https://sage-docs.vercel.app/docs/configuration#file)
- Model Context Protocol server for validated, agent-operated searches
- Built-in support for reading gzipped mzML and local Thermo Fisher RAW files
- Support for reading/writing directly from [AWS S3](https://sage-docs.vercel.app/docs/configuration/aws), Google Cloud, or Azure.

### Enriched linear retention-time model

Linear regression with the original feature set remains the default retention-time model. Enriched sequence features and additive variable-PTM offsets can be enabled with:

```json
"retention_time_model": {
  "features": "additive_ptm",
  "folds": 3,
  "seed": 42,
  "ptm_regularization": 25.0
}
```

The `additive_ptm` feature set uses exact N1/C1 residue indicators, position-binned hydrophobicity and residue-property fractions. It then learns cross-validated, regularized offsets for each configured variable modification plus residue/terminal identity. Static modifications are excluded. Use `"features": "physicochemical"` for enriched linear features with generic modification descriptors.

When nonzero precursor ion-mobility values are present, Sage can independently fit a peptide-grouped, cross-validated linear mobility model:

```json
"ion_mobility_model": {
  "enabled": true,
  "features": "additive_ptm",
  "folds": 3,
  "seed": 42,
  "ptm_regularization": 25.0,
  "min_training_psms": 200
}
```

The enriched mobility model adds exact N1/C1 residues, physicochemical composition, peptide mass, m/z and charge features. `additive_ptm` learns regularized variable-PTM offsets with more strongly shrunk charge-specific deviations. It does not fit PTM-pair interactions, ignores static modifications, and falls back to the basic linear mobility features if the enriched fit fails.

## Interoperability

Sage is well-integrated into the open-source proteomics ecosystem. The following projects support analyzing results from Sage (typically in addition to other tools), or redistribute Sage binaries for use in their pipelines. 

- [SearchGUI](http://compomics.github.io/projects/searchgui): a graphical user interface for running searches
- [PeptideShaker](http://compomics.github.io/projects/peptide-shaker): visualize peptide-spectrum matches
- [MS2Rescore](http://compomics.github.io/projects/ms2rescore): AI-assisted rescoring of results
- [Picked group FDR](https://github.com/kusterlab/picked_group_fdr): scalable protein group FDR for large-scale experiments
- [sagepy](https://github.com/theGreatHerrLebert/sagepy): Python bindings to the sage-core library
- [quantms](https://github.com/bigbio/quantms): nextflow pipeline for running searches with Sage
- [OpenMS](https://github.com/OpenMS/OpenMS): Sage is included as a "TOPP" tool in OpenMS
- [sager](https://github.com/UCLouvain-CBIO/sager): R package for analyzing results from Sage searches
- [Sage results to mzIdentML](https://github.com/magnuspalmblad/shic/blob/main/shims/Peptide_identification_in_TSV_to_Peptide_identification_in_mzIdentML.sh): Bash script to convert `results.sage.tsv` files to mzIdentML
- [i2MassChroQ](http://pappso.inrae.fr/bioinfo/i2masschroq/): a graphical user interface for proteomics analysis
- [annotator](https://github.com/snijderlab/annotator): a graphical user interface for visualizing peptide-spectrum matches
- [rustyms](https://github.com/snijderlab/rustyms): a Rust library (with Python bindings) to handle peptides and identified peptide files
- If your project supports Sage and it's not listed, please open a pull request! If you need help integrating or interfacing with Sage in some way, please reach out.

Check out the (now outdated) [blog post introducing the first version of Sage](https://lazear.github.io/sage/) for more information and full benchmarks!
