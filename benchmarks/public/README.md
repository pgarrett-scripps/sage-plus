# A public benchmark anyone can rerun

Sage Plus reads Thermo `.raw` files directly. Upstream Sage reads mzML, so a
`.raw` dataset has to be converted first. This benchmark measures what that
costs end to end: the same files, the same search settings, the same machine,
with the conversion counted against the arm that needs it.

Everything here is public: the dataset is on PRIDE, the FASTA comes from
UniProt, the converter is a digest-pinned container, and both engines build
from a tagged commit.

## What is measured, and what is not

| Reported | How |
| --- | --- |
| Conversion time | wall clock around the converter, mzML arm only |
| Search time | `/usr/bin/time -v` around the search |
| Total time | conversion plus search |
| Peak memory | `Maximum resident set size` from the same `time -v` |
| PSMs, peptides, proteins at 1% FDR | counted from each arm's own `results.sage.parquet` by one function, so neither engine's summary format is trusted over the other |

Not measured here: quantification, localization, or anything only one of the
two engines has. The point of the comparison is the cost of getting from a
vendor file to identifications.

Identification counts are *not* expected to match exactly. The two engines
differ in scoring and FDR details, and a single search is not reproducible to
the last PSM in either engine: the linear-discriminant model is fitted on the
run's own PSMs, so repeats of the same search move by tens of PSMs. Treat a
difference under a few tenths of a percent as a tie, and run three repeats if
a number matters.

## Inputs

### Dataset: PXD028735

A three-species HYE mixture acquired on an Orbitrap, five DDA files:

```
LFQ_Orbitrap_DDA_Condition_A_Sample_Alpha_01.raw
LFQ_Orbitrap_DDA_Condition_A_Sample_Beta_01.raw
LFQ_Orbitrap_DDA_Condition_B_Sample_Alpha_01.raw
LFQ_Orbitrap_DDA_Condition_B_Sample_Beta_01.raw
LFQ_Orbitrap_DDA_Human_01.raw
```

They are in the project's PRIDE archive:
<https://www.ebi.ac.uk/pride/archive/projects/PXD028735>. Each `.raw` file is
about 3.5 GB, so the five together need about 17 GB of disk, and about 6 GB
more for the mzML arm: a converted file is about 1.2 GB.

### FASTA

The mixture needs the three proteomes it was made from, plus the iRT standard
peptides. Download the UniProt reference proteomes and concatenate them:

- human `UP000005640`, yeast `UP000002311`, *E. coli* `UP000000625` from
  <https://www.uniprot.org/proteomes/>, one FASTA each.

Record the UniProt release you downloaded: proteome contents change between
releases, so the absolute identification counts move with them. Both arms must
use the same file, which is what makes the comparison fair; which release it is
only matters when comparing against numbers published elsewhere.

Proteins containing ambiguous residues (`X`, `B`, `Z`, `J`) are searched by Sage
Plus and skipped by some other tools. If you want a FASTA both engines treat
identically, drop those proteins before searching.

### Engines

- **Sage Plus**: `cargo build --release` at the tag you want to test, from
  <https://github.com/pgarrett-scripps/sage-plus>. The binary reports its
  version with `sage --version`.
- **Upstream Sage**: the release binary from
  <https://github.com/lazear/sage/releases>, or `cargo build --release -p sage-cli`
  at the matching tag. Upstream asks to be told not to send telemetry; the
  script passes `--disable-telemetry-i-dont-want-to-improve-sage` for you.

Build both with the same toolchain, and note that a bundled HDF5 in the mzMLb
feature needs `CMAKE_POLICY_VERSION_MINIMUM=3.5` under CMake 4.

### Converter

ThermoRawFileParser 1.4.5, pinned by digest so everyone runs the same build:

```
quay.io/biocontainers/thermorawfileparser@sha256:465b623c81a8f300d2802196246307e4280e9135ebfcc4c4526a2beac7e7ea8f
```

`SAGE_BENCH_CONVERTER` points the script at a local ThermoRawFileParser
instead, called as `--input=FILE --output_directory=DIR --format=1`.
`msconvert` from a pinned ProteoWizard container is an equally valid converter,
with its own flags (`--mzML --outdir`); if you use it, say so, because the
converter is part of the mzML arm's cost.

## Running it

```shell
export SAGE_BENCH_PLUS=/path/to/sage-plus/target/release/sage
export SAGE_BENCH_UPSTREAM=/path/to/sage/target/release/sage
export SAGE_BENCH_FASTA=/path/to/hye-irt.fasta       # or --fasta
export SAGE_BENCH_THREADS=8                          # optional, else every core

# Sage Plus, straight from the vendor files
./run_public_benchmark.sh --arm plus-raw --out runs/plus \
  --fasta "$SAGE_BENCH_FASTA" /data/PXD028735/*.raw

# Upstream Sage, converting first
./run_public_benchmark.sh --arm upstream-mzml --out runs/upstream \
  --fasta "$SAGE_BENCH_FASTA" /data/PXD028735/*.raw
```

Each run writes into its `--out` directory:

- `receipt.json`: the arm, both binaries' versions, conversion and search
  seconds, peak RSS, and the three identification counts
- `config.used.json`: the configuration as searched, with the paths filled in
- `search.time`, `search.stdout`, `convert.*.log`: the raw evidence
- `search/`: the engine's own output, including `results.sage.parquet`

Two differences between the engines' command lines are handled for you, and
neither changes the search: upstream Sage has no `--overwrite`, so each run
gets a fresh output directory, and it writes TSV unless asked for `--parquet`,
which the script does so both arms are counted from the same file format. The
thread count comes from `RAYON_NUM_THREADS`, which both engines honour, rather
than from a flag only one of them has.

`pxd028735.json` holds the search settings, and only settings both engines
implement: tryptic with one missed cleavage, 7 to 50 residues, carbamidomethyl
on C, variable oxidation on M and N-terminal acetylation at most two per
peptide, 10 ppm precursor and 20 ppm fragment tolerance, isotope errors −1 to
3, b and y ions, deisotoping on, at least 6 matched peaks, one PSM reported per
spectrum. The script fills in `fasta` and `mzml_paths` and leaves everything
else alone, so both arms search the same space.

## A one-file pilot, to show what the output looks like

One PXD028735 file (`Condition_A_Sample_Alpha_01`, 3.5 GB) through both arms,
8 threads each, one machine with nothing else running. This is a pilot that
proves the harness works, not a published result: one file, one repeat.

| | Sage Plus 0.1.0-beta.18, `.raw` | Upstream Sage 0.15.0-beta.2, mzML |
| --- | ---: | ---: |
| Conversion | none | 250.2 s |
| Search | 21.0 s | 29.1 s |
| Total | **21.0 s** | **279.3 s** |
| Peak memory | 5.16 GiB | 7.48 GiB |
| PSMs at 1% FDR | 61,274 | 59,861 |
| Peptides at 1% FDR | 41,581 | 40,421 |
| Proteins at 1% FDR | 5,784 | 5,723 |

Conversion is 90% of the mzML arm's wall time here, and it is the part that
scales with the dataset: the same five-file run is about 1 minute of Sage Plus
against about 23 minutes, of which 21 minutes is conversion. Identifications
agree within 2.4% of each other, which is the range two engines with different
scoring and FDR details should be expected to land in; see the caveat above
before reading anything into a difference that size.

## Reporting honestly

Run the arms one at a time on an idle machine. Timing on a shared or thermally
limited machine is not comparable: record the machine, the thread count, and
whether anything else was running. If you have a scheduler, take the whole
machine for the timing runs.

State the versions of both engines, the converter, the UniProt release and the
dataset files. A benchmark without those is not reproducible, whichever
direction it points.
