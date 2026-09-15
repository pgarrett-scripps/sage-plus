# Sage Plus Development and Benchmark Report

[Read the report](report.pdf)

This repository technical report explains the changes in Sage Plus and evaluates
their consequences using the frozen Sage `v0.15.0-beta.2` and Sage Plus
`v0.1.0-beta.3` executables. It includes resource use, identification agreement,
independent peptide entrapment, quantification, PTM acceptance, compatibility
failures, and incomplete workloads. Favorable and unfavorable results are retained.

## Attribution and status

[Sage](https://github.com/lazear/sage) was developed by Michael R. Lazear and the
upstream Sage contributors. Credit for the original engine and its continuing
maintenance belongs to that project. Please cite the
[original Sage paper](https://doi.org/10.1021/acs.jproteome.3c00486) when using Sage
or Sage Plus, and identify the exact downstream version when applicable.

Sage Plus is independently maintained and is not an official Sage release. This
report is project documentation, with no individual author byline. It does not
imply authorship, review, or endorsement by the upstream creators or maintainers.
It is not a peer-reviewed paper. Contribution history and licenses remain in the
repository. Dataset and methodology references are included in the report.

## What is included

- `paper.typ`, `config.typ`, `si-body.typ`, and `references.bib` contain the text.
- `figures/` and `si/` contain the generated figures and tables used in the report.
- `stats.json` and `assets.json` connect the text and assets to their generators.
- `analysis/data/` contains the reviewed development inventory and frozen report
  snapshots, including peptide-level ratios derived from public mixture data.
- `../benchmarks/scientific-results/20260914/` contains the pilot summaries,
  original failure records, and audit receipts used by the report.
- `analysis/scripts/` regenerates figures, tables, and prose statistics.

Large spectra, reference downloads, raw search outputs, the local evidence
archive, virtual environments, and Word exports are not included. The contextual
local HEK measurements retain their provenance limitation. The underlying local
spectrum is not distributed. Earlier exploratory assets and unrelated working-tree
repairs are outside this report package.

## Rebuild from the included snapshots

Install Typst 0.14 or newer, just, uv, Python 3.11 or newer, and typstyle.
The first build may download the Typst template and locked Python dependencies.

```sh
cd paper
just assets
just fmt
just docx
just paper
just verify
just check-stats-deep
```

This rebuild uses the small included snapshots. It does not rerun searches or
require the original `/data/` paths. The optional Word export is generated locally
as `paper.docx`. The compiled `paper.pdf` is also local build output. The reviewed
public snapshot is `report.pdf`. After changing report sources, refresh it with:

```sh
just public-report
```

The verification gate checks formatting, prose extraction, declared statistics,
asset provenance, and build freshness. The deep statistics check recomputes the
values used in the prose. These checks verify report construction, not general
scientific validity of either search engine.

## Repeating the underlying experiments

See [the pilot protocol](../benchmarks/SCIENTIFIC_PILOT.md) and the supporting
information. Original execution receipts retain historical local paths and
content hashes as provenance. Those paths are not public download links and
must be mapped to newly acquired inputs for another execution.

The `collect_report.py` and `audit_report.py` scripts describe derivation from
the retained raw results. Running them requires the pinned executables, external
tools, original input files, and full search outputs. The original frozen
snapshots should not be overwritten by new experiments. No full external evidence
deposition is claimed.

## Maintenance

Report errors and reproducibility questions belong in the
[Sage Plus issue tracker](https://github.com/pgarrett-scripps/sage-plus/issues).
The build scaffold comes from `pgarrett-scripps/paper-scaffold`. Its license is
preserved in `LICENSE.scaffold`. Editing guidance is in `AGENTS.md` and `STYLE.md`.
