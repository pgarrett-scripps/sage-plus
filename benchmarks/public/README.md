# A public benchmark anyone can rerun

Sage Plus reads Thermo `.raw` files directly. Upstream Sage reads mzML, so a
`.raw` dataset has to be converted first. This benchmark measures what that
costs end to end, and separates it from the search itself: the same files, the
same search settings, the same machine, three arms.

| Arm | Engine | Input | What it isolates |
| --- | --- | --- | --- |
| `plus-raw` | Sage Plus | the vendor `.raw` files | the whole cost of going from vendor files to identifications |
| `plus-mzml` | Sage Plus | the converted mzML | the same engine on the other input, so the reader's share is visible |
| `upstream-mzml` | upstream Sage | the same converted mzML | the other engine on byte-identical input, so the engine's share is visible |

The third arm is why the comparison is fair. `plus-raw` against
`upstream-mzml` mixes two effects — a different reader and a different engine —
and `plus-mzml` sits between them: compared with `plus-raw` it shows what the
reader costs, and compared with `upstream-mzml` it shows what the engine does
on exactly the same input.

Conversion time is always reported as its own number. There is no single
"faster by N times" figure here, because such a figure hides whether the
difference was conversion or search.

Everything here is public: the dataset is on PRIDE, the FASTA comes from
UniProt, the converter is a digest-pinned container, and both engines build
from a tagged commit.

## What is measured, and what is not

| Reported | How |
| --- | --- |
| Conversion time | wall clock around the converter, per file, recorded beside the mzML; the mzML arms report the sum for the files they read |
| Search time | `/usr/bin/time -v` around the search |
| Total time | conversion plus search, both shown |
| Peak memory | `Maximum resident set size` from the same `time -v` |
| PSMs, peptides, proteins at 1% FDR | counted from each arm's own `results.sage.parquet` by one function, so neither engine's summary format is trusted over the other |

Not measured here: quantification, localization, or anything only one of the
two engines has. The point of the comparison is the cost of getting from a
vendor file to identifications.

Identification counts are *not* expected to match exactly, and the section
["Where the two programs differ"](#where-the-two-programs-differ) lists every
reason they do not. On top of that, a single search is not reproducible to the
last PSM in either engine: the linear-discriminant model is fitted on the run's
own PSMs, so repeats of the same search move by tens of PSMs. Treat a
difference under a few tenths of a percent as a tie, and run three repeats if a
number matters.

## Inputs

### Dataset: PXD028735

A three-species HYE mixture acquired on an Orbitrap, five DDA files:

```
LFQ_Orbitrap_DDA_Condition_A_Sample_Alpha_01.raw
LFQ_Orbitrap_DDA_Condition_A_Sample_Beta_01.raw
LFQ_Orbitrap_DDA_Condition_A_Sample_Gamma_01.raw
LFQ_Orbitrap_DDA_Condition_B_Sample_Alpha_01.raw
LFQ_Orbitrap_DDA_Condition_B_Sample_Beta_01.raw
```

Any five DDA files from the project work; these are the ones the table below
was measured on, so quote the file names with any number you publish. They are
in the project's PRIDE archive:
<https://www.ebi.ac.uk/pride/archive/projects/PXD028735>. Each `.raw` file is
about 3.5 GB, so the five together need about 18 GB of disk, and about 6 GB
more for the mzML arms: a converted file is about 1.2 GB.

### FASTA

The mixture needs the three proteomes it was made from, plus the iRT standard
peptides. Download the UniProt reference proteomes and concatenate them:

- human `UP000005640`, yeast `UP000002311`, *E. coli* `UP000000625` from
  <https://www.uniprot.org/proteomes/>, one FASTA each.

Record the UniProt release you downloaded: proteome contents change between
releases, so the absolute identification counts move with them. Every arm must
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
converter is part of the mzML arms' cost.

## Running it

```shell
export SAGE_BENCH_PLUS=/path/to/sage-plus/target/release/sage
export SAGE_BENCH_UPSTREAM=/path/to/sage/target/release/sage
export SAGE_BENCH_FASTA=/path/to/hye-irt.fasta       # or --fasta
export SAGE_BENCH_THREADS=8                          # optional, else every core
CACHE=runs/mzml                                      # the converted files

# Convert once, one file at a time, each conversion timed on its own
for raw in /data/PXD028735/*.raw; do
  ./run_public_benchmark.sh --arm convert --out runs/convert --mzml-cache "$CACHE" "$raw"
done

# Arm (a): Sage Plus, straight from the vendor files
./run_public_benchmark.sh --arm plus-raw --out runs/plus-raw \
  --fasta "$SAGE_BENCH_FASTA" /data/PXD028735/*.raw

# Arm (b): Sage Plus on the converted mzML
./run_public_benchmark.sh --arm plus-mzml --out runs/plus-mzml \
  --mzml-cache "$CACHE" --fasta "$SAGE_BENCH_FASTA" /data/PXD028735/*.raw

# Arm (c): upstream Sage on the same mzML
./run_public_benchmark.sh --arm upstream-mzml --out runs/upstream-mzml \
  --mzml-cache "$CACHE" --fasta "$SAGE_BENCH_FASTA" /data/PXD028735/*.raw
```

Every arm is given the vendor files; the mzML arms look up each file's
conversion in the cache. A missing conversion is done and timed on the spot, so
the arms run correctly without the loop above, but then the conversion shares a
timing run with the search — convert first if the numbers are going to be
published. `receipt.json` says which it was: `convert_time_source` is
`recorded` for a cache hit and `measured` when the arm converted.

Each run writes into its `--out` directory:

- `receipt.json`: the arm, the binary and its version, conversion seconds
  (total and per file), search seconds, peak RSS, and the three identification
  counts
- `config.used.json`: the configuration as searched, with the paths filled in
- `search.time`, `search.stdout`, `convert.*.log`: the raw evidence
- `search/`: the engine's own output, including `results.sage.parquet`

Two differences between the engines' command lines are handled for you, and
neither changes the search: upstream Sage has no `--overwrite`, so each run
gets a fresh output directory, and it writes TSV unless asked for `--parquet`,
which the script does so every arm is counted from the same file format. The
thread count comes from `RAYON_NUM_THREADS`, which both engines honour, rather
than from a flag only one of them has.

`pxd028735.json` holds the search settings, and only settings both engines
implement: tryptic with one missed cleavage, 7 to 50 residues, carbamidomethyl
on C, variable oxidation on M and N-terminal acetylation at most two per
peptide, 10 ppm precursor and 20 ppm fragment tolerance, isotope errors −1 to
3, b and y ions, deisotoping on, at least 6 matched peaks, one PSM reported per
spectrum. The script fills in `fasta` and `mzml_paths` and leaves everything
else alone, so every arm searches the same space.

## Where the two programs differ

The configuration pins every setting both programs have and that matters for
the search space. What it cannot pin is a setting one program does not have.
These are the differences that remain, and they are why the identification
counts are close but not equal.

| Setting | Sage Plus | Upstream Sage 0.15.0-beta.2 | Why identifications can differ |
| --- | --- | --- | --- |
| `database.clip_n_term_met` | `true` by default | no such option: the initiator methionine is never removed | Sage Plus also searches the Met-clipped form of every protein N-terminal peptide, so its search space is slightly larger |
| `database.merge_isoleucine_leucine` | `true` by default | no such option | isoleucine and leucine are one residue when peptides are deduplicated, so an I/L-isobaric pair is one target rather than two; this moves peptide counts and target-decoy competition |
| `retention_time_alignment` | `nonlinear` by default: robust outlier filtering, then a monotone piecewise-linear warp | ordinary least squares, not configurable | the `delta_rt_model` column of the discriminant is computed from a different alignment, so PSMs are ranked differently |
| `retention_time_model` | cross-validated in three folds, with ridge-penalised retention-time offsets per variable modification | one linear fit over all PSMs | same column, second reason |
| `mass_recalibration` | `off` by default, but the post-hoc fragment alignment still feeds the model: the discriminant's eighth column is `aligned_fragment_ppm` | the same column is `average_ppm`, the raw mean | one of the twenty discriminant columns holds a different number |
| spectrum q-values | one FDR estimate per tied-score group | one per row in sort order | a handful of PSMs on the 1% boundary fall on the other side |
| `output_filter.psm_q_value` | `0.1`: PSMs worse than 10% spectrum q are not written | every PSM is written | the parquet files have different row counts; the counts at 1% FDR do not change, because q-values are computed before the filter |
| reader | `sage-plus-raw` reads the `.raw` file | ThermoRawFileParser writes mzML, which upstream reads | centroiding, charge and metadata decisions belong to whoever read the vendor file; arm (a) against arm (b) is the measurement of this |

Everything else both programs default to the same way, and the configuration
does not override it: `min_matched_peaks` 4 (pinned to 6 here), `min_peaks` 15,
`max_peaks` 150, `precursor_charge` 2–4, `isotope_errors` 0–0 (pinned to −1–3),
`deisotope` on, `predict_rt` on, `chimera` off, `wide_window` off,
`report_psms` 1, `bucket_size` 8192, `min_ion_index` 2, `decoy_tag` `rev_`,
decoy generation on, peptide mass 500–5000 Da, `max_variable_mods` 2, protein
grouping on at a 1% peptide FDR, and the prefilter off.

Two conclusions follow, and the arms are what make them separable:

- Arm (b) against arm (c) is the engine difference on byte-identical input. No
  reader is involved, so whatever gap remains comes from the rows above.
- Arm (a) against arm (b) is the reader difference in the same engine. A gap
  here is the vendor reader and the converter disagreeing about the spectra, not
  a scoring difference.

<!-- RESULTS -->

## Reporting honestly

Run the arms one at a time on an idle machine. Timing on a shared or thermally
limited machine is not comparable: record the machine, the thread count, and
whether anything else was running. If you have a scheduler, take the whole
machine for the timing runs, and keep the conversions out of them.

State the versions of both engines, the converter, the UniProt release and the
dataset files. A benchmark without those is not reproducible, whichever
direction it points.
