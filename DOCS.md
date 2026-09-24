# Sage Plus Documentation

This document covers Sage Plus configuration, outputs, and downstream features. The
[upstream Sage documentation](https://sage-docs.vercel.app/docs) remains the reference for
general Sage concepts.


## Features & Information

### Assign multiple peptides to complex spectra

<img src="figures/chimera_27525.png" width="800">

- When chimeric searching is enabled, multiple peptide identifications can be reported for each MS2 scan

### Sage trains machine learning models for FDR refinement and posterior error probability calculation

- Retention times are globally aligned across runs
- Boosts PSM identifications using prediction of retention times with a [linear regression](https://doi.org/10.1021/ac070262k) model
- Hand-rolled, 100% pure Rust implementations of Linear Discriminant Analysis and KDE-mixture models for refinement of false discovery rates
- Models demonstrate 1:1 results with scikit-learn, but have increased performance
- No need for a second post-search pipeline step

<img src="figures/SageLDA.png" width="600px">

## Installation

Sage Plus is distributed as source code, prebuilt release archives, and a versioned container.
It is not currently published through Conda.

### Installing upstream Sage via Conda

The [Bioconda package](https://anaconda.org/bioconda/sage-proteomics) installs upstream Sage, not
Sage Plus. It does not include the downstream features documented here. Use it only when you
specifically want the upstream distribution:

```
$ conda install -c bioconda -c conda-forge sage-proteomics
$ sage --help
```

### Compiling the development version

1. Install the [Rust programming language compiler](https://rustup.rs/)
2. Download Sage Plus source code via git: `git clone https://github.com/pgarrett-scripps/sage-plus.git` or by [zip file](https://github.com/pgarrett-scripps/sage-plus/archive/refs/heads/main.zip)
3. Compile: `cargo build --release --workspace`
4. Run: `./target/release/sage config.json`

Once you have Rust installed, you can copy and paste the following lines into your terminal to complete the above instructions, and run Sage on the example mzML provided in the repository (a single scan from PXD016766)

```sh
git clone https://github.com/pgarrett-scripps/sage-plus.git
cd sage-plus
cargo run --release tests/config.json 
```

### Downloading a Sage Plus release

1. Visit the [Sage Plus releases](https://github.com/pgarrett-scripps/sage-plus/releases) website.
2. Download the correct pre-compiled binary for your operating system.
3. Run: `sage <path/to/config.json>`

### Interfacing with AWS S3

Sage Plus can natively read and write files through AWS S3:

- S3 paths should be specified as `s3://bucket/prefix/key.mzML.gz` or `s3://bucket/prefix` for output folder
- See [AWS docs](https://docs.aws.amazon.com/sdk-for-rust/latest/dg/credentials.html) for configuring your credentials
- Using S3 may incur data transfer charges as well as multi-part upload request charges.

## Usage 

```shell
Usage: sage [OPTIONS] <parameters> [mzml_paths]...

🔮 Sage 🧙 - Proteomics searching so fast it feels like magic!

Arguments:
  <parameters>     Path to configuration parameters (JSON file)
  [mzml_paths]...  Paths to mzML, MGF, Bruker TDF, or Thermo RAW files to process. Overrides files listed in the configuration file.

Options:
  -f, --fasta <fasta>
          Path to FASTA database. Overrides the FASTA file specified in the configuration file.
  -o, --output_directory <output_directory>
          Path where search and quant results will be written. Overrides the directory specified in the configuration file.
      --batch-size <batch-size>
          Number of files to search in parallel (default = number of CPUs/2)
      --write-pin
          Write percolator-compatible `.pin` output files
      --max-memory <GiB>
          Abort if Sage's memory use exceeds this many GiB, to keep the system responsive
          (default: 90% of total RAM; 0 disables). Also settable via SAGE_MAX_MEMORY_GB.
      --events-jsonl <PATH>
          Stream versioned JSONL job events to PATH (use '-' for stdout)
      --validate-only
          Validate the configuration and overrides without running a search
  -h, --help
          Print help information
  -V, --version
          Print version information
```

Sage is called from the command line using and requires a path to a JSON-encoded parameter file as an argument (see below). 

Example usage: `sage config.json`

Some options in the parameters file can be over-written using the command line interface. These are:

1. The paths to the mzML data
2. The path to the database (fasta file)
3. The output directory

For example: 

```
# Specify fasta and output dir:
sage -f proteins.fasta -o output_directory config.json

# Specify mzML files:
sage -f proteins.fasta config.json *.mzML

# Specify mzML file located in an S3 bucket
sage config.json s3://my-bucket/YYYY-MM-DD_expt_A_fraction_1.mzML.gz
```

Running Sage will produce several output files (located in either the current directory, or `output_directory` if that option is specified):
- A record of search parameters (`results.json`) and a portable basic-statistics artifact (`run-summary.json`) are created for every successful search
- MS2 search results are stored in `results.sage.parquet`. TMT reporter-ion values, when enabled, are a nested array on each PSM row.
- Label-free quantification is stored separately in long-form `lfq.parquet`, with one precursor/file row.
- `results.json` records the effective configuration and `run-summary.json` records portable run statistics and output paths.

Local output directories must be fresh unless `--overwrite` or `"overwrite": true` is explicit.
Overwrite removes known Sage artifacts, including stale optional outputs, and preserves unrelated
files. It does not coordinate concurrent writers. Use one directory per job and a fresh object-store
prefix for remote searches. A failed run can leave partial analytical files without a completion summary.

Run-summary schema 9 adds `warnings`, `provenance`, per-file mass-alignment outcomes, and
`execution.rayon_threads`. The existing `execution.parallelism` field describes file batching.
The provenance mode `path_size_mtime` records local input metadata, not content hashes or a
cryptographically exact build identity. Benchmark manifests separately record SHA-256 hashes of
inputs and binaries. Older summaries remain readable through defaults for the new fields.

For library callers, `JobOptions.parallel` remains the fallback file batch size when configuration
does not specify `batch_size`. It does not set the Rayon worker count. CLI and MCP batch overrides
take precedence over the configuration.

Parquet is the canonical analytical output format. Sage does not emit parallel TSV copies of the PSM, LFQ, matched-fragment, or PTM-site result tables. Purpose-specific interchange artifacts such as Percolator `.pin` files and the reusable PTM-library TSV remain available.

The versioned physical schemas and score definitions are published in [`schemas/`](schemas/). Canonical Parquet files embed `sage.schema.name` and `sage.schema.version` metadata so downstream tools can select the matching contract.

#### Memory guard

A search can balloon in memory — most often during database generation, where the number of modified peptide variants grows combinatorially with `max_variable_mods` / `max_peff_variable_mods`, the FASTA size, and enzyme settings. To prevent a runaway search from exhausting RAM and freezing the host, Sage runs a lightweight background watchdog that terminates the process **cleanly** (exit code 137) if either:

- Sage's own resident memory exceeds a ceiling (default: 90% of total system RAM), or
- system-wide available memory drops below a small safety floor (max of 1 GiB or 2% of RAM).

The ceiling is set with `--max-memory <GiB>` (or the `SAGE_MAX_MEMORY_GB` environment variable); `--max-memory 0` disables the guard entirely. The watchdog polls a few times per second from a single thread and adds no overhead to the allocation hot path. When it trips it prints how to reduce the search size (e.g. lower `max_variable_mods` / `max_peff_variable_mods`, use a smaller FASTA, narrow tolerances, or enable `prefilter`).

#### Sequence-ambiguity annotation

Every PSM row carries two additional columns, `ambiguity_sequence` and `mass_shift`, that encode which residues are actually supported by fragment-ion evidence (a native port of the [SagePeptideAmbiguityAnnotator](https://github.com/pgarrett-scripps/SagePeptideAmbiguityAnnotator) tool):

- **ambiguity_sequence**: the peptide string in which any run of residues lacking *both* forward (a/b/c) and reverse (x/y/z) ion cleavage evidence is wrapped in `(?...)`. For example `(?LQ)SRPAAPPAPGPGQLTLR` means the leading `L`/`Q` could be reordered without changing the matched peaks. When the experimental precursor mass does not match the peptide's calculated mass (e.g. in an open search), the residual mass is placed using the same coverage:
  - `...T[+79.96633]...` — localized to a single residue,
  - `(...)[+mass]` — confined to a region but not a single residue,
  - a leading `{+mass}` — labile / cannot be localized (forward and reverse coverage overlap).
- **mass_shift**: the residual `expmass - calcmass` (in Da) that was placed, or `0.0` when the precursor matches within `mass_shift_ppm`.

These are computed for every search; mods are rendered in the same `[+mass]`/`[Name]` notation as the `peptide` column. The threshold used to decide whether a precursor delta mass is a real shift is configurable via the top-level **`mass_shift_ppm`** parameter (default: 50.0). It is deliberately independent of `precursor_tol`, so wide/open searches still surface and place real shifts.

## Machine-readable jobs

Use `--validate-only` to check configuration and CLI overrides without reading the FASTA,
spectra, or creating the output directory:

```shell
sage config.json --validate-only
```

The committed [JSON Schema](schemas/config.schema.json) provides editor completion and static
validation for configuration files. An installed binary can copy its matching schema to a file or
standard output:

```shell
sage --write-config-schema sage-config.schema.json
sage --write-config-schema -
```

Use `--events-jsonl <path>` to stream versioned, newline-delimited JSON events while a
search runs. `--events-jsonl -` writes events to standard output. Human-readable logs remain
on standard error, so standard output can be consumed directly by workflow engines and other
applications.

```shell
sage config.json --events-jsonl run.events.jsonl
```

Every event contains `schema_version`, a monotonically increasing `sequence`, `elapsed_ms`,
and an `event` discriminator. Events cover configuration validation, database construction,
file reads, spectra processing, search progress, model fitting or fallback, FDR, written
outputs, and terminal job state. Consumers should ignore unknown fields and event names so
that compatible events can be added to schema version 1.

Rust callers can use `sage_cli::api::SageRunner` rather than invoking the CLI. `JobOptions`
accepts an `EventEmitter` and a cloneable `CancellationToken`; `run` returns a structured
`RunSummary` alongside telemetry. This application layer is intended to be shared by future
protocol servers and user interfaces.

### MCP server for AI clients

The `sage-mcp` binary exposes the runner to MCP-compatible coding agents and assistants over
local standard input/output. Build it with `cargo build --release -p sage-mcp`, then configure
the client to launch it with a directory that contains every allowed configuration and input:

```shell
sage-mcp --root /path/to/allowed/data
```

The server can inspect and validate configurations, estimate database expansion and memory,
start approved background searches, monitor or cancel jobs, summarize completed runs, and make
basic analysis from the portable run summary, and bounded queries over TSV PSM and PTM-site results. Searches require `approved: true`, remote URLs
are disabled, local inputs cannot escape `--root`, and outputs are written beneath
`ROOT/.sage/jobs`. See `crates/sage-mcp/README.md` for client configuration and tool details.

## Configuration file schema

### Notes

- The majority of parameters are optional - only "database.fasta", "precursor_tol", and "fragment_tol" are required. Sage will try and use reasonable defaults for any parameters not supplied
- Tolerances are specified on the *experimental* m/z values. To perform a -100 to +500 Da open search (mass window applied to *theoretical*), you would use `"da": [-500, 100]`

### Decoys

Using decoy sequences is critical to controlling the false discovery rate in proteomics experiments. Sage can use decoy sequences in the supplied FASTA file, or it can generate internal sequences. Sage reverses tryptic peptides (not proteins), so that the [picked-peptide](https://pubmed.ncbi.nlm.nih.gov/36166314/) approach to FDR can be used.

If `database.generate_decoys` is set to true (or unspecified), then decoy sequences in the FASTA database matching `database.decoy_tag` will be *ignored*, and Sage will internally generate decoys. It is __critical__ that you ensure you use the proper `decoy_tag` if you are using a FASTA database containing decoys and have internal decoy generation turned on - otherwise Sage will treat the supplied decoys as hits!

Internally generated decoys will have protein accessions matching "{decoy_tag}{accession}", e.g. if `decoy_tag` is "rev_" then a protein accession like "rev_sp|P01234|HUMAN" will be listed in the output file.

### FASTA digestion

Sage will process a protein into peptides via several routes listed below. Currently, one and only one is supported.

- Enzymatic: `database.enzyme.cleave_at = "KR"` - configuration option set to a sequence of amino acids (e.g. "KR" for trypsin, "FWYL" for chymotrypsin)
- Non-enzymatic: `database.enzyme.cleave_at = ""` - All potential peptides between `min_len` and `max_len` will be generated from the sequence
- No digestion: `database.enzyme.cleave_at = "$"` - FASTA entries will be used as-is, subject to `min_len` and `max_len` options

Protein-specific cleavage sites can be added to any FASTA digest with
`database.custom_cleavage_sites = "cleavage-sites.tsv"`. TSV and Parquet files
are supported and require `protein` and `position` columns. In Parquet,
`protein` must be UTF-8 and `position` must be an integer. `position` is the
zero-based index of the residue immediately before the cut; for example,
position `0` cuts between the first and second residues. An optional UTF-8
`context` column validates a short sequence window, with `|` marking the cut:

```text
protein	position	context
P12345	86	KLGF|APQT
```

Both products adjacent to each site are generated using the configured enzyme
boundaries, missed-cleavage allowance, length, mass, and modification limits.
Normal digest peptides remain unchanged. Context mismatches and terminal or
out-of-range positions are errors; sites without context are accepted with a
warning.


### Example configuration file

For additional information about configuration options and output file formats, please see [the new documentation](https://sage-docs.vercel.app/docs)

```jsonc
// Note that json does not allow comments, they are here just as explanation
// but need to be removed in a real config.json file
{
  "database": {
    "bucket_size": 32768,           // Maximum fragments in each internal search bucket
    "enzyme": {               // Optional. Default is trypsin, using the parameters below
      "missed_cleavages": 2,  // Optional[int], Number of missed cleavages for tryptic digest
      "min_len": 5,           // Optional[int] {default=5}, Minimum AA length of peptides to search
      "max_len": 50,          // Optional[int] {default=50}, Maximum AA length of peptides to search
      "cleave_at": "KR",      // Optional[str] {default='KR'}. Amino acids to cleave at
      "restrict": "P",        // Optional[str] {default='P'}. Do not cleave if one of these AAs follows the cleavage site
      "c_terminal": false,      // Optional[bool] {default=true}. Cleave at c terminus of matching amino acid
      "semi_enzymatic": false      // Optional[bool] {default=false}. Generate semi-enzymatic peptides
    },
    "peptide_min_mass": 500.0,      // Optional[float] {default=500.0}, Minimum monoisotopic mass of peptides to fragment
    "peptide_max_mass": 5000.0,     // Optional[float] {default=5000.0}, Maximum monoisotopic mass of peptides to fragment
    "ion_kinds": ["b", "y"],        // Optional[List[str]] {default=["b","y"]} Which fragment ions to generate and search?
    "min_ion_index": 2,     // Optional[int] {default=2}, Do not generate b1/b2/y1/y2 ions for preliminary searching. Does not affect full scoring of PSMs
    "static_mods": {
      "TMT": {"mass": 304.207, "sites": ["peptide_n_term", "K"]},
      "Carbamidomethyl": {"mass": 57.0215, "sites": ["C"]}
    },
    "variable_mods": {
      "Oxidation": {"mass": 15.9949, "sites": ["M"], "max_count": 1},
      "Acetyl": {"mass": 42.0106, "sites": ["K"], "max_count": 1},
      "PyroGlu-Q": {"mass": -17.026549, "sites": ["first_residue:Q"]},
      "PyroGlu-E": {"mass": -18.010565, "sites": ["first_residue:E"]}
    },
    "max_variable_mods": 2, // Optional[int] {default=2} Limit modifications on each peptide
    "max_total_variable_mods": 2, // Exhaustive + PTM-library placements
    "max_combinations": 8,  // Optional[int] {default=null} Limit total variants per peptide
    "decoy_tag": "rev_",    // Optional[str] {default="rev_"}: See notes above
    "generate_decoys": false, // Optional[bool] {default="true"}: Ignore decoys in FASTA database matching `decoy_tag`
    "fasta": "dual.fasta",  // str: mandatory path to FASTA file
    "custom_cleavage_sites": "cleavage-sites.tsv" // Optional protein-specific sites
  },
  "quant": {                // Optional - specify only if TMT or LFQ
    "tmt": "Tmt16",         // Optional[str] {default=null}, one of "Tmt6", "Tmt10", "Tmt11", "Tmt16", or "Tmt18"
    "tmt_settings": {
      "level": 3,           // Optional[int] {default=3}, MS-level to perform TMT quantification on
      "sn": false           // Optional[bool] {default=false}, use Signal/Noise instead of intensity for TMT quant. Requires noise values in mzML
    },
    "lfq": true,            // Optional[bool] {default=null}, perform MS1 feature quantification
    "lfq_settings": {
      "peak_scoring": "Hybrid", // See DOCS.md for details - recommend that you do not change this setting
      "integration": "Sum",   // Optional["Sum" | "Apex"], use sum of MS1 traces in peak, or MS1 intensity at peak apex
      "spectral_angle": 0.7,  // Optional[float] {default = 0.7}, normalized spectral angle cutoff for calling an MS1 peak
      "ppm_tolerance": 5.0,    // Optional[float] {default = 5.0}, tolerance (in p.p.m.) for DICE window around calculated precursor mass
      "rt_pct_tolerance": 0.5, // Optional[float] {default = 0.5}, symmetric match-between-runs RT tolerance as percent of total gradient length
      "mbr": true,             // Optional[bool] {default = true}, trace precursors into runs without direct MS2 evidence
      // Optional[bool] {default = true}. Combine all charge states for quantification. Setting this to false
      // quantifies each peptide-charge precursor in `precursor_charge` range (see below) separately
      "combine_charge_states": true
    }
  },
  "precursor_tol": {        // Tolerance can be either "ppm" or "da"
    "da": [
      -500,                 // This value is substracted from the experimental precursor to match theoretical peptides
      100                   // This value is added to the experimental precursor to match theoretical peptides
    ]
  },
  "fragment_tol": {         // Tolerance can be either "ppm" or "da"
    "ppm": [
     -10,                   // This value is subtracted from the experimental fragment to match theoretical fragments 
     10                     // This value is added to the experimental fragment to match theoretical fragments 
    ]
  },
  // Optional[Tuple[int, int]] {default=[2, 4]}
  // If charge states are not annotated in the mzML, or if `wide_window` mode is turned on, then consider
  // all precursors at z=2, z=3, z=4
  "precursor_charge": [2, 4]
  "isotope_errors": [       // Optional[Tuple[int, int]] {default=[0,0]}: C13 isotopic envelope to consider for precursor
    -1,                     // Consider -1 C13 isotope
    3                       // Consider up to +3 C13 isotope (-1/0/1/2/3) 
  ],
  "deisotope": true,        // Optional[bool|object] {default=true}: perform MS2 deisotoping and charge deconvolution
  "chimera": false,         // Optional[bool] {default=false}: search for chimeric/co-fragmenting PSMS
  "wide_window": false,     // Optional[bool] {default=false}: _ignore_ `precursor_tol` and search in wide-window/DIA mode
  "predict_rt": false,    // Optional[bool] {default=true}: use retention time prediction model as a feature for LDA
  "ion_mobility_model": {
    "enabled": false       // Optional[bool] {default=true}: retain observed mobility without fitting a prediction model
  },
  "retention_time_alignment": "nonlinear", // Optional["linear" | "nonlinear"]: explicitly enable observed-RT alignment
  "min_peaks": 15,          // Optional[int] {default=15}: only process MS2 spectra with at least N peaks
  "max_peaks": 150,         // Optional[int] {default=150}: take the top N most intense MS2 peaks to search,
  "min_matched_peaks": 6,   // Optional[int] {default=4}: minimum # of matched b+y ions to use for reporting PSMs
  "max_fragment_charge": 1, // Optional[int] {default=null}: maximum fragment ion charge states to consider,
  "report_psms": 1,         // Optional[int] {default=1}: number of PSMs to report for each spectra. Higher values might disrupt PSM rescoring.
  "output_filter": {         // Optional: rows written to PSM and matched-fragment Parquet files
    "psm_q_value": 0.1       // Optional[float] {default=0.1}: maximum spectrum-level q-value, inclusive
  },
  "max_memory_gb": 16,      // Optional[float] {default=null}: stop Sage if its resident memory reaches this many GiB; 0 disables
  "min_free_memory_gb": 2,  // Optional[float] {default=null}: stop Sage if system-available memory falls to this many GiB; 0 disables
  "batch_size": 1,          // Optional[int] {default=# of CPUs/2}: number of input files to load and search at once
  "output_directory": "s3://bucket/prefix", // Optional[str] {default=`.`}: Place output files in a given directory or S3 bucket/prefix
  "mzml_paths": [           // List[str]: representing paths to mzML (or gzipped-mzML) files for search
    "local/path.mzML",
    "s3://bucket/PXD0000001/foo.mzML.gz"
  ]       
}
```

## Using the docker image

Sage can be used from a docker image!

```shell
$ docker pull ghcr.io/pgarrett-scripps/sage-plus:v0.1.0-beta.7
$ docker run -it --rm -v ${PWD}:/data ghcr.io/pgarrett-scripps/sage-plus:v0.1.0-beta.7 sage -o /data /data/config.json
# The sage executable is located in /app/sage in the image
```

Container images currently target Linux AMD64. Use the native ARM64 executable archive from
the GitHub release on ARM64 systems until a native multi-architecture image is available.

> `-v ${PWD}:/data` means it will mount your current directory as `/data`
> in the docker image. Make sure all the paths in your command and configuration
> use the location in the image and not your local directory

# Further Details

This documentation covers the parameters in the JSON configuration file for the proteomics search engine. The configuration file contains information about the search engine's settings, including database, enzyme, modifications, and other settings. For a complete example of a configuration file, please see the [online docs](https://sage-docs.vercel.app/docs)

## Database

- **bucket_size**: Integer. The maximum number of theoretical-fragment records in one internal search bucket (default: 8192). Values are normalized to a power of two. A mass-prefix group containing more records is split into multiple buckets with the same prefix, while smaller groups allocate only the records they contain. Tweaking this parameter can affect search performance for wide precursor or fragment searches.

The preliminary fragment index stores each fragment as a 32-bit peptide index and a 16-bit exact
mass suffix. The shared upper bits of the original `f32` mass are stored once per mass-prefix
group. Recombining the prefix and suffix reproduces the original floating-point bit pattern, so
this representation does not round or quantize fragment masses. Empty mass-prefix groups are not
stored.

### Enzyme

The enzyme section contains parameters related to the enzyme used for digestion. The default enzyme is trypsin, with the parameters specified below.

- **missed_cleavages**: Integer. The number of missed cleavages for tryptic digest (default: 1).
- **min_len**: Integer. The minimum amino acid (AA) length of peptides to search (default: 5).
- **max_len**: Integer. The maximum AA length of peptides to search (default: 50).
- **cleave_at**: String. Amino acids to cleave at (default: 'KR').
- **restrict**: String. Do not cleave if one of these amino acids follows the cleavage site (default: 'P').
- **c_terminal**: Boolean. Cleave at the C-terminus of matching amino acids (default:true).

Example: 
```json
"database": {
  "enzyme": {
    "missed_cleavages": 1,
    "min_len": 5,
    "max_len": 50,
    "cleave_at": "KR",
    "restrict": "P",
    "c_terminal": true
  }
}
```

### Fragment Settings

- **peptide_min_mass**: Float. The minimum monoisotopic mass of peptides to fragment *in silico* (default: 500.0).
- **peptide_max_mass**: Float. The maximum monoisotopic mass of peptides to fragment *in silico* (default: 5000.0).
- **ion_kinds**: List of strings. Which fragment ions to produce? Allowed values: "a", "b", "c", "x", "y", "z". (default: ["b", "y"])
- **min_ion_index**: Integer. Do not generate b1/bN/y1/yN ions for preliminary searching if `min_ion_index = N`. Does not affect full scoring of PSMs (default: 2).

Example:
```json
"database": {
  "peptide_min_mass": 500.0,
  "peptide_max_mass": 5000.0,
  "ion_kinds": ["b", "y"],
  "min_ion_index": 2
}
```

### Modifications

Define each modification once under its stable name in `static_mods` or
`variable_mods`. Each definition contains `mass` and a nonempty `sites` list.
The same site syntax works for static, indexed variable, and mass-offset search.

```json
{
  "database": {
    "static_mods": {
      "Carbamidomethyl": {"mass": 57.021464, "sites": ["C"]}
    },
    "variable_mods": {
      "Acetyl": {
        "mass": 42.010565,
        "sites": ["first_residue:K", "internal_residue:K", "protein_last:K"],
        "max_count": 2
      },
      "Phospho": {
        "mass": 79.966331,
        "sites": ["S", "T", "Y"],
        "max_count": 3,
        "neutral_losses": [97.976896]
      }
    },
    "max_variable_mods": 3
  }
}
```

The dictionary key is the modification identity used by localization and site
libraries. A separate `name` field is unnecessary. If present, it must equal the
key. Equal masses do not merge distinct identities. All sites in one definition
share its `max_count`, neutral losses, channel offsets, and search policy.
Use distinct names when those policies need to differ.

#### Explicit site vocabulary

| Site | Attachment |
| --- | --- |
| `K` | Any K residue |
| `first_residue:K` | First peptide residue, if K |
| `internal_residue:K` | K strictly between the first and last peptide residues |
| `last_residue:K` | Last peptide residue, if K |
| `protein_first:K` | First protein residue, if K |
| `protein_last:K` | Last protein residue, if K |
| `peptide_n_term` | Peptide N-terminal group, any boundary residue |
| `peptide_c_term` | Peptide C-terminal group, any boundary residue |
| `protein_n_term` | Protein N-terminal group |
| `protein_c_term` | Protein C-terminal group |
| `peptide_n_term:K` | Peptide N-terminal group only when the first residue is K |
| `peptide_c_term:K` | Peptide C-terminal group only when the last residue is K |
| `protein_n_term:K` | Protein N-terminal group only when the first residue is K |
| `protein_c_term:K` | Protein C-terminal group only when the last residue is K |

Substitute any supported uppercase one-letter residue for K. Each string encodes
one complete rule. Entries are alternatives. Overlapping entries count the same
physical attachment once.

A terminal-group attachment and a modification on its adjacent residue are distinct
sites and can coexist if the limits allow it. For example,
`["peptide_n_term:K", "first_residue:K"]` allows either attachment on a peptide
starting with K. A limit of one permits alternatives, while a limit of two also
permits both simultaneously. Identical fragment evidence does not establish which
attachment occurred.

Use `["peptide_n_term", "peptide_c_term"]` to permit either terminal group regardless
of boundary residue. No residue wildcard or additional field is required.

#### Static and variable behavior

Static definitions apply at matching unoccupied sites after variable modifications,
preserving the existing variable-before-static behavior. Conflicting fixed
definitions that can occupy the same site are rejected. One name cannot appear in
both static and variable sections.

For indexed variable modifications, `max_count` limits occurrences of that identity
across all its sites. `max_variable_mods` limits exhaustive placements per peptide.
`max_total_variable_mods` limits exhaustive and library-supported placements combined
and defaults to `max_variable_mods`. `max_combinations` caps generated variants,
including the unmodified form. Existing mass-offset limitations are described below.

Optional `neutral_losses` contains positive fragment-loss masses.
`neutral_loss_mode` is `optional` by default or `required` to omit the retained
fragment form. Required mode needs at least one loss. Masses and channel offsets
must be finite. Static definitions do not accept variable-only limits or policies.

#### Positional residue modifications

The Acetyl example includes first and internal lysines, plus protein-last lysines.
It therefore includes H3K9 when K is the first residue of `KSTGGKAPR`.
An internal-only rule would exclude that placement. Protein position is evaluated
from the digest context, independently of which enzyme generated the peptide.

Length-one and length-two peptides have no internal residues. On a length-one
peptide, first and last refer to the same residue and do not double-apply a mod.

#### Mass offset modifications

Set `search_mode` to `mass_offset` on a named variable definition:

```json
{
  "variable_mods": {
    "Phospho": {
      "mass": 79.966331,
      "sites": ["S", "T", "Y"],
      "search_mode": "mass_offset"
    },
    "Oxidation": {"mass": 15.994915, "sites": ["M"]}
  }
}
```

An offset is searched at scoring time instead of expanded into the fragment index.
Each spectrum is searched against a translated precursor window with both shifted
and unshifted fragment lookups. Every eligible placement competes as a candidate.
Offsets and ordinary indexed candidates are scored and target-decoy competed together.

At most one offset copy is placed on a peptide. This is independent of `max_count`,
`max_variable_mods`, `max_total_variable_mods`, and `max_combinations`, which govern
indexed modifications. Multiple offsets are not combined. Offset and indexed
modifications can coexist at different unoccupied sites.

Offsets follow explicit sites and typed library restrictions. They require nonzero
mass, cannot use `channel_offsets`, and are limited to 254 distinct definitions.
Fragments retain the offset mass, or mass minus the first required neutral loss.
The placed peptidoform supplies mass error, FDR, quantification, localization, and
output identity. The offset is not reported as precursor mass error. Prefiltering
uses the same offset-aware retrieval.

#### Preview modification placement

```shell
sage config.json --preview-modifications KSTGGKAPR
sage config.json --preview-modifications ASQKSTGGK --peptide-position cterm --preview-limit 100
sage library-config.json --preview-modifications KSTGGKAPR --preview-protein P68431 --preview-start 9
```

Preview prints explicit rules, eligible physical sites, limits, and bounded generated
variants without reading spectra. `--peptide-position` supplies protein boundary
context and defaults to `internal`. Library preview loads the configured TSV or
Parquet library and requires an accession and one-based peptide start coordinate.
It trusts the supplied sequence and boundary context rather than loading a FASTA.
The limit ranges from 1 to 10000. `truncated` reports when returned variants were
limited. Static and variable occupancy is reflected in generated variants.

#### Migration from symbol keys

Existing residue-keyed numeric and structured configurations remain readable.
New named definitions accept only explicit spellings in `sites`. Do not mix named
and legacy entries within one section.

```shell
sage old-config.json --migrate-modifications > new-config.json
```

The command prints a converted configuration and leaves its input untouched.
Repeated named entries are grouped only when their definitions agree. Unnamed
entries receive distinct deterministic `legacy_static_mods_N` or
`legacy_variable_mods_N` identities, preserving separate occurrence limits.
It does not infer chemical identity from mass.

Legacy bare `^`, `$`, `[`, and `]` become the corresponding terminal-group names.
Legacy `^K`, `$K`, `[K`, and `]K` become first/last residue rules, preserving their
meaning. `~K` becomes `internal_residue:K`.

#### Modification channels

Both static and variable definitions may contain `channel_offsets`. Effective mass
is the base mass plus the selected offset. All channel-aware modifications on one
peptide resolve to one coherent channel.

```json
{
  "static_mods": {
    "SILAC-K": {
      "mass": 0.0,
      "sites": ["K"],
      "channel_offsets": {"light": 0.0, "heavy": 8.014199}
    },
    "SILAC-R": {
      "mass": 0.0,
      "sites": ["R"],
      "channel_offsets": {"light": 0.0, "heavy": 10.008269}
    }
  }
}
```

All channel dictionaries must contain the same names and at least two chemically
distinct channels. The unique channel whose offsets are all zero is the reference.
Variable channel modifications consume normal indexed-modification budgets.
Zero-offset duplicates are searched once while channel partners remain available
for LFQ. Complete incorporation and fixed mass shifts are assumed. Partial
incorporation and isotope-purity correction are not modeled.

#### PTM site libraries

A library contains observed locations. Modification chemistry and search settings
remain in the named definitions. `site_mode` selects:

- `exhaustive`, the default, generates candidates from the configured sites.
- `library` requires a matching typed library attachment as well as a configured site.
- `both` allows exhaustive placement and recognizes matching library-supported sites.

Set `database.ptm_library` to `{"path": "sites.tsv", "strict": true}`.
TSV, TSV.gz, and Parquet are supported. When a library is configured, indexed
variable modifications require `max_count`. Names are supplied by dictionary keys.
Library-supported indexed placements bypass `max_variable_mods` but still consume
`max_total_variable_mods` and the shared per-modification limit. Library evidence
never overrides the configured site restrictions.

Beta 6 writes five columns: `protein`, `position`, `residue`, `modification`, and
`attachment`. Position is a one-based protein coordinate. Modification is the exact
configured identity. Attachment is `residue`, `peptide_n_term`, `peptide_c_term`,
`protein_n_term`, or `protein_c_term`. Additional evidence columns may be present.

For a terminal-group record, position identifies its adjacent boundary residue.
A peptide N-terminal record supports only peptides starting there, and a C-terminal
record supports only peptides ending there. Protein-terminal records additionally
require the corresponding protein-terminal digest context. Terminal observations
are exported with peptide-terminal attachment names and remain constrained by the
configured protein-terminal rules when reused.

A residue observation never supports a terminal-group attachment at the same
coordinate. Library unions must include attachment in their identity key.
Legacy four-column libraries are read as residue attachments. Manually authored
legacy terminal libraries must add explicit attachment values before reuse.
Consumers that insist on exactly four columns need an update before reading Beta 6
output. Do not silently discard attachment when merging libraries.

With localization enabled, FASTA searches emit
`results.sage.ptm-library.tsv` and `results.sage.ptm-library.parquet`.
Typed Parquet library and site reports embed schema version 2. Equal-scoring or
indistinguishable attachment alternatives are excluded from the reusable library,
even with a permissive localization threshold. Preserve the named definitions
alongside the library, because the location table does not embed chemical masses.

### Decoys

- **decoy_tag**: String. The tag used to identify decoy entries in the FASTA database (default: "rev_").
- **generate_decoys**: Boolean. If true, ignore decoys in the FASTA database matching `decoy_tag`, and generate internally reversed peptides (default: false).

### FASTA

- **fasta**: String. The path to the FASTA file, either a local path or s3 object URI.
- **prefilter**: Boolean. Retain only peptides that can contribute a preliminary fragment match
  before building the search index. The spectra are indexed once, and the database is generated
  in chunks and streamed through the spectrum index, so no fragment index is built for discarded
  peptides. Targets, paired decoys, and label-channel partners are retained together, so the final
  search uses the same FDR competition and produces the same results as a full database search.
  The spectrum index is limited to a quarter of `max_memory_gb`, or 8 GiB without a limit. Larger
  inputs are indexed in file batches, and the database is streamed once per batch. Set the
  `SAGE_PREFILTER_INDEX_GB` environment variable to override the budget.
- **prefilter_chunk_size**: Integer. Approximate number of FASTA sequences per generated chunk.
  A value of zero selects the chunk size from the estimated number of modified peptides.
- **prefilter_low_memory**: Deprecated and ignored. Exact prefiltering always uses compact survivor
  tracking.

## Quantification

The quant section is optional and should be specified only if TMT or LFQ is used. Precursor channel
chemistry remains on modification definitions through `channel_offsets`. LFQ automatically becomes
channel-aware when these offsets are configured.


- **tmt**: String. One of "Tmt6", "Tmt10", "Tmt11", "Tmt16", or "Tmt18" (default: null).
- **tmt_settings**: Object containing TMT-specific settings.
  - **level**: Integer. The MS-level to perform TMT quantification on (default: 3).
  - **sn**: Boolean. Use Signal/Noise instead of intensity for TMT quantification. Requires noise values in mzML (default: false).
- **lfq**: Boolean. Perform MS1 feature quantification. This is label-free without channel-aware
  modifications and channel-aware when `channel_offsets` are configured (default: null).
- **lfq_settings**: Object containing LFQ-specific settings.
  - **peak_scoring**: String. The method used for scoring peaks in LFQ, one of: "Hybrid", "RetentionTime", "SpectralAngle" (default: "Hybrid").
  - **integration**: String. The method used for integrating peak intensities, either "Sum" or "Max" (default: "Sum").
  - **spectral_angle**: Float. Threshold for the spectral angle similarity measure, ranging from 0 to 1 (default: 0.7).
  - **ppm_tolerance**: Float. Tolerance for matching MS1 ions in parts per million (default: 5.0).
  - **rt_pct_tolerance**: Float. Symmetric retention-time tolerance for match-between-runs, as a percentage of total gradient length (default: 0.5). For example, `0.5` searches +/-0.5% around the aligned retention time.
  - **mbr**: Boolean. Trace identified precursors into runs without direct MS2 evidence. Set this to `false` to quantify a precursor only in runs where it was identified (default: true).

Example: 
```json
 "quant": {
    "tmt": "Tmt16",
    "tmt_settings": {
      "level": 3,
      "sn": false
    },
    "lfq": true,
    "lfq_settings": {
      "peak_scoring": "Hybrid",
      "integration": "Sum",
      "spectral_angle": 0.7,
      "ppm_tolerance": 5.0,
      "rt_pct_tolerance": 0.5,
      "mbr": true
    }
  }
```


## Precursor Tolerance

- **precursor_tol**: Dictionary with either "ppm" or "da" as keys, and lists of two integers as values (default: {}).
  - Example: Tolerance of [-500, 100] in daltons.
    ```json
    "precursor_tol": {
      "da": [-500, 100]
    }
    ```

## Fragment Tolerance

- **fragment_tol**: Dictionary with either "ppm" or "da" as keys, and lists of two integers as values (default: {}).
  - Example: Tolerance of [-10, 10] in parts per million.
    ```json
    "fragment_tol": {
      "ppm": [-10, 10]
    }
    ```

## Isotope Errors

- **isotope_errors**: List of two integers. The C13 isotopic envelope to consider for precursor (default: [0, 0]).
  - Example: Consider -1 and up to +3 C13 isotopes (-1/0/1/2/3).
    ```json
    "isotope_errors": [-1, 3]
    ```

**NOTE**: Searching with isotope errors is slower than searching with a wider precursor tolerance that encompasses the isotope errors, e.g. `"da": [-3.5, 1.25]`. Using the wider precursor tolerance will generally increase the number of confidently identified PSMs as well.

## Other Settings

Note on the settings below:

Retention-time alignment and prediction are separate features. `retention_time_alignment` aligns observed times even when `predict_rt` is false. Prediction uses aligned times and therefore runs linear alignment when no method is specified. LFQ also requires alignment, but does not require retention-time prediction.

- **deisotope**: Boolean or object. Perform scored averagine deisotoping and charge state deconvolution on MS2 spectra (default: true). Use `false` to disable deconvolution or an object to tune the bounded isotope-envelope scoring. Sage excludes the reporter-ion region from MS2 deisotoping when TMT or iTRAQ quantification is configured.

  ```json
  "deisotope": {
    "enabled": true,
    "ppm_tolerance": 10.0,
    "max_charge": null,
    "min_envelope_peaks": 2,
    "max_envelope_peaks": 4,
    "min_score": 0.45,
    "max_isotope_log2_ratio": 1.5
  }
  ```

  `ppm_tolerance` controls isotope-spacing matches. `max_charge` optionally caps the precursor-derived charge search. Envelope sizes are bounded between two and four peaks. `min_score` is the minimum Bhattacharyya isotope-pattern score required to merge an envelope and treat its charge as known. `max_isotope_log2_ratio` limits the difference between observed and averagine-predicted adjacent isotope ratios. Boolean `true` uses the object defaults shown above. When mzML or mzMLb provides the fragment charge binary array `MS:1000516`, positive values constrain isotope-envelope assignment and are used directly for charge-aware fragment matching. Zero values remain unknown and use scored inference.
- **chimera**: Boolean. Search for chimeric/co-fragmenting PSMs (default: false).
- **wide_window**: Boolean. Ignore `precursor_tol` and search spectra in wide-window/dynamic precursor tolerance mode (default: false).
- **predict_rt**: Boolean. Use retention time prediction model as a feature for LDA (default: true).
- **ion_mobility_model.enabled**: Boolean. Fit and use the ion-mobility model when mobility observations are present (default: true). Set this to `false` to keep observed mobility data without fitting predictions.
  - Example:
    ```json
    "ion_mobility_model": {
      "enabled": false
    }
    ```
- **retention_time_alignment**: Explicitly align observed retention times across experiments. `"linear"` uses Sage's existing ordinary least-squares alignment. `"nonlinear"` enables robust outlier filtering followed by a monotone piecewise-linear warp. This operates independently of `predict_rt`.
- **min_peaks**: Integer. Only process MS2 spectra with at least N peaks (default: 15).
- **max_peaks**: Integer. Take the top N most intense MS2 peaks to search (default: 150).
- **min_matched_peaks**: Integer. The minimum number of matched b+y ions to use for reporting PSMs (default: 4).
- **max_fragment_charge**: Integer. The maximum fragment ion charge states to consider (default: null - use precursor z-1).
- **report_psms**: Integer. The number of PSMs to report for each spectrum. Higher values might disrupt LDA (default: 1).
- **annotate_matches**: Boolean. Write `matched_fragments.sage.parquet` for PSMs passing `output_filter.psm_q_value` (default: false). Detailed annotations are reconstructed in a batched post-FDR MS2 pass rather than allocated for every candidate during scoring. When PTM localization is also enabled, both operations share the same spectrum reread. Chimera ranks replay preceding-rank peak removal before annotation.
- **spectral_library**: Object. Build an empirical library from confident target PSMs. See [Empirical Spectral Libraries](#empirical-spectral-libraries).
- **output_filter.psm_q_value**: Float from 0 to 1. Maximum spectrum-level PSM q-value written to `results.sage.parquet` and `matched_fragments.sage.parquet` (default: 0.1). The boundary is inclusive. Set it to `1.0` to retain every scored PSM. This is an output-only filter: scoring, FDR estimation, LFQ, PTM localization, `.pin` output, and the HTML report continue to use their existing inputs and thresholds. Target and decoy PSMs that pass the threshold are retained so downstream target-decoy analyses remain possible.
- **max_memory_gb**: Number. Abort the search if Sage's resident memory reaches this many GiB. Zero disables this limit (default: disabled).
- **min_free_memory_gb**: Number. Abort the search if system-available memory falls to this many GiB, preserving capacity for the operating system and other applications. Zero disables this limit (default: disabled).
- **batch_size**: Integer. Number of input files to load and search at once. Smaller values reduce temporary spectrum memory at the cost of throughput (default: half the number of CPUs, with a minimum of one). The `--batch-size` command-line option overrides this value.

When either memory limit is enabled, Sage estimates the unmodified digest, variable-modification expansion, and fragment/index sizes before allocating them. Unsafe database searches return an error before expansion begins. Estimates are conservative and are backed by a runtime memory monitor for allocations outside database construction.

## Empirical Spectral Libraries

Set `spectral_library.enabled` to build a library directly from the spectra identified in the
current search. `--spectral-library` is a shortcut that enables the feature with the configured
values or defaults.

```json
"spectral_library": {
  "enabled": true,
  "psm_q_value": 0.01,
  "peptide_q_value": 0.01,
  "strategy": "best_psm",
  "min_matched_peaks": 6,
  "max_fragments": 20,
  "min_relative_intensity": 0.01,
  "min_consensus_psms": 1,
  "min_fragment_frequency": 0.5,
  "include_chimeric": false,
  "formats": ["sage_parquet", "mzspeclib"]
}
```

The current `best_psm` strategy groups eligible target PSMs by exact modified peptide and
precursor charge, then chooses one representative deterministically: lowest spectrum q-value,
lowest peptide q-value, highest discriminant score, highest hyperscore, and finally lowest PSM
ID. By default, only rank-one PSMs are eligible. `include_chimeric: true` also permits later
chimera ranks. The PSM and peptide cutoffs are independent of `output_filter.psm_q_value`.

The `consensus` strategy combines every eligible PSM in each peptidoform and charge group.
Retention time, aligned retention time, ion mobility, and normalized fragment intensities use
robust medians. `min_consensus_psms` controls the minimum group size. Groups of one remain valid
by default. `min_fragment_frequency` controls the fraction of supporting spectra in which a
fragment must appear before it enters the consensus spectrum.

For the selected spectrum, Sage retains matched fragments at or above
`min_relative_intensity`, keeps at most `max_fragments` by observed intensity, and reports them
in theoretical-m/z order. Intensities are normalized to the most intense retained candidate
peak. Detailed annotations are reconstructed in the same deferred MS2 pass used by matched-ion
output and PTM localization, so enabling more than one of these features does not add a separate
spectrum reread.

Available formats are:

- `sage_parquet`: `spectral_library.sage.parquet`, the canonical long-form table with one row
  per transition. It includes source-spectrum provenance, mass-delta ProForma, precursor data,
  aligned retention time, ion mobility, q-values, supporting-PSM count, fragment identity,
  theoretical fragment m/z, and relative intensity. Its versioned schema is in
  `schemas/spectral_library.sage.v1.parquet.schema`.
- `mzspeclib`: `spectral_library.mzspeclib.txt`, a PSI mzSpecLib 1.0 text library containing
  singleton or consensus spectra and mzPAF peak annotations.

This empirical export is distinct from `database.ptm_library`, which restricts which protein
modification sites are searched and is not a spectral library.

- **ptm_localization**: Object. Configure PTM site localization and site-level reports. See [PTM Site Localization](#ptm-site-localization).
  - **enabled**: Boolean. Enable localization (default: false). The `--localize` CLI flag is a shortcut that sets this to true.
  - **psm_q_value**: Float from 0 through 1. Spectrum-level identification q-value cutoff for PSMs localized and included in the site reports (default: 0.01). It is not a PTM localization probability or false-localization-rate threshold.
  - **localization_q_value**: Float from 0 through 1. Arrangement-level false localization rate cutoff for reported PTM localizations (default: 0.01).

## PTM Site Localization

When `ptm_localization.enabled` is true, sage attempts to pinpoint which physical attachment carries each variable modification on a confidently-identified peptide, analogous to MaxQuant's site tables or MSFragger/PTMProphet.

Example configuration:

```json
"ptm_localization": {
  "enabled": true,
  "psm_q_value": 0.01,
  "localization_q_value": 0.01
}
```

For each FDR-passing target PSM (spectrum q-value ≤ `ptm_localization.psm_q_value`), and for each distinct variable-modification identity it carries, sage:
1. recovers candidate residue and terminal attachments from the configured sites and library restrictions,
2. enumerates every way to distribute the modification(s) across those candidate sites, keeping all other modifications pinned,
3. re-scores each arrangement against the experimental spectrum using only *site-determining ions* (fragments whose mass differs between arrangements), and
4. scores a balanced set of impossible-site decoy arrangements alongside the valid target arrangements,
5. converts target/decoy competition scores across the dataset into monotonic localization q-values, and
6. reports target arrangements at or below `ptm_localization.localization_q_value`, together with an AScore-style delta and per-site localization probabilities.

The current implementation combines one AScore-inspired, site-determining-ion strategy with balanced impossible-site target/decoy competition. It is intentionally not presented as a configurable strategy yet: a future strategy name should select a genuinely different, validated scoring or FLR model rather than act as an alias for the same calculation.

Two Parquet site reports are written:

- **results.sage.ptm-sites.parquet**: one row per localized modification site of each PSM. Columns include `peptide`, `modification`, `position` (1-based, within the peptide), `residue`, `localization_probability`, `delta_localization_score`, `target_decoy_score`, `localization_q_value`, `candidate_sites`, site-determining-ion counts, and `site_probabilities`.
- **results.sage.protein-sites.parquet**: the best localization for each (protein, modified peptide site) aggregated across all supporting PSMs, including `best_localization_q_value`.

For example, the PSM-site report contains rows shaped like this (positions are 1-based within the peptide):

```text
psm_id  peptide            modification  position  residue  localization_probability  localization_q_value  site_probabilities
42      AAS[+79.966]AATAA  Phospho       3         S        0.982                     0.008                 S3:0.982;T6:0.018
```

Notes:
- Residue and terminal-group modifications are localized within their configured rules. Site reports include `attachment`, and terminal coordinates identify the adjacent residue.
- The current backbone-ion model cannot distinguish equal masses on a terminal group and its adjacent residue in some configurations. Such alternatives retain ambiguity and are excluded from reusable libraries. Terminal localization is experimental and has synthetic regression coverage, not an empirical FLR calibration study.
- Localization runs after spectrum FDR assignment and only for passing target PSMs. Sage re-reads MS2 spectra for this optional pass rather than retaining the full experiment in memory.
- `ptm_localization.psm_q_value` controls identification quality; `ptm_localization.localization_q_value` controls arrangement-level localization FLR. `localization_probability` remains a within-PSM marginal site probability.
- A modification without enough eligible impossible residues to construct a balanced decoy search space is not included in the FDR-controlled reports.
- FASTA searches preserve protein coordinates during indexing. The canonical PSM output attaches each protein accession to its one-based inclusive start and end positions plus the preceding and following amino acids. Pre-digested peptide TSV and spectral-library inputs omit coordinates when they are unavailable.

## Spectrum Paths

- **mzml_paths**: List of strings. Despite the legacy field name, Sage accepts mzML, mzMLb, MGF, Bruker TDF, and Thermo Fisher RAW inputs. mzML and MGF paths may be local or use a configured object-store URL. mzMLb, Thermo RAW, and Bruker TDF inputs must be local because their readers require seekable files. mzMLb support is included in standard builds and release binaries. Minimal source builds created with `--no-default-features` omit it. Files ending in ".gz" or ".gzip" are inferred to be compressed. MGF spectra that cannot be searched (no TITLE, no PEPMASS, or no peaks) are skipped with one warning per file; unparseable values and missing BEGIN IONS/END IONS markers stop the search.
  - Thermo RAW input uses centroid peak lists directly. TMT signal-to-noise mode (`quant.tmt_settings.sn: true`) still requires mzML containing a noise array.
  - Bruker TDF ion mobility (1/K0) uses each frame's `TimsCalibration` model from `analysis.tdf`, matching the Bruker SDK. Every scan is converted before MS1 centroiding, and DDA precursors convert their fractional average scan. DIA window centers use the calibration row shared by most frames. Inputs that point at a file inside the `.d` directory, such as `analysis.tdf_bin`, read the same `analysis.tdf`. Inputs without an `analysis.tdf`, such as miniTDF `.ms2` directories, have no calibration table and fall back to the linear scale with a warning. Only ModelType 2 is supported; other models stop the search. Set `"bruker_config": {"ion_mobility_scale": "linear"}` to reproduce the uncalibrated scale of Beta 6 and earlier, which interpolates between the acquisition limits. `run-summary.json` records the scale as `models.ion_mobility_scale`. Mobility tolerances are relative and apply unchanged.
  - Example:
    ```json
    "mzml_paths": [
      "local/path.mzML",
      "local/path.mzMLb",
      "local/path.raw",
      "s3://my-mass-spec-data/PXD0000001/foo.mzML.gz"
    ]
    ```
  
## Output directory:

- **output_directory**: Local directory, or S3 location where output files will be written. If the local directory does not already exist, it will be created. Write permissions are required for the directory or S3 path.
  - Possible analytical output files are `results.sage.parquet`, `lfq.parquet`, `matched_fragments.sage.parquet`, `results.sage.ptm-sites.parquet`, `results.sage.protein-sites.parquet`, and `spectral_library.sage.parquet`. Optional purpose-specific artifacts include `spectral_library.mzspeclib.txt`, `results.sage.pin`, the HTML report, and PTM-library Parquet/TSV files. `results.json` and `run-summary.json` are always written after a successful run; the summary contains runtime, database size, 1% FDR counts, localized-PTM counts and thresholds, spectral-library entries and transitions, model/alignment outcomes, quantification counts, memory and batching controls, input-format counts, modification-expansion limits, and output paths.
  - Example:
  ```json
  "output_directory": "s3://my-mass-spec-results/PXD003881/"
  ```

# Interpreting Sage Output

The `results.sage.parquet` file contains the following columns:

Rows satisfy the configured `output_filter.psm_q_value` threshold. The same PSM IDs define the rows emitted to `matched_fragments.sage.parquet`, so that file never contains fragments for a PSM omitted from the main result table. Both files record the effective threshold as `sage.output_filter.spectrum_q_max` in Parquet key-value metadata.

- `peptide`: Peptide sequence, including modifications (e.g., NC\[+57.021\]HKGSFK).
- `proteins`: Proteins containing the peptide sequence.
- `protein_sites`: Typed list of protein occurrences. Each item contains `protein`, one-based inclusive `start` and `end`, plus nullable `prev_aa` and `next_aa` flanking residues.
- `num_proteins`: Number of proteins assigned to the peptide sequence.
- `filename`: File containing this PSM
- `scannr`: Spectrum identifier from mzML file.
- `rank`: Rank of the PSM. If `report_psms > 1`, then the best match will have rank = 1, the second best match will have rank = 2, etc. 
- `label`: Target/Decoy label (-1: decoy, 1: target).
- `expmass`: Experimental mass of the peptide.
- `calcmass`: Calculated mass of the peptide.
- `charge`: Reported precursor charge.
- `pepide_len`: Length of the peptide sequence.
- `missed_cleavages`: Number of missed cleavages.
- `isotope_error`: C13 isotope error.
- `precursor_ppm`: Difference between experimental mass and calculated mass, reported in parts-per-million.
- `fragment_ppm`: Average parts-per-million (delta mass) for matched fragment ions compared to theoretical ions.
- `hyperscore`: X!Tandem hyperscore for the PSM.
- `delta_next`: Difference between the hyperscore of this candidate and the next best candidate.
- `delta_bext`: Difference between the hyperscore of the best candidate (rank=1) and this candidate.
- `rt`: Retention time.
- `aligned_rt`: Globally aligned retention time.
- `predicted_rt`: Predicted retention time, if enabled.
- `delta_rt_model`: Difference between predicted and observed retention time.
- `matched_peaks`: Number of matched theoretical fragment ions.
- `longest_b`: Longest b-ion series.
- `longest_y`: Longest y-ion series.
- `longest_y_pct`: Longest y-ion series, divided by peptide length (as a percentage).
- `matched_intensity_pct`: Fraction of MS2 intensity explained by matched b- and y-ions (as a percentage of total MS2 intensity for this spectrum).
- `scored_candidates`: Number of scored candidates for this spectrum.
- `poisson`: Probability of matching exactly N peaks across all candidates (Pr(x=k)).
- `sage_discriminant_score`: Combined score from linear discriminant analysis, used for FDR (False Discovery Rate) calculation.
- `posterior_error`: Posterior error probability for this PSM / local FDR.
- `spectrum_q`: Assigned spectrum-level q-value.
- `peptide_q`: Assigned peptide-level q-value.
- `protein_q`: Assigned protein-level q-value.
- `ms1_intensity`: Intensity of the selected MS1 precursor ion (not label-free quant)
- `ms2_intensity`: Total intensity of MS2 spectrum

These columns provide comprehensive information about each candidate peptide spectrum match (PSM) identified by the Sage search engine.

## Label-free quantification output

`lfq.parquet` is a separate long-form table with one row per quantified precursor and acquisition file. All intensities are produced by Sage's cross-run feature-tracing workflow; `ms2_confirmed` records whether that precursor also has an accepted target PSM in the specific file.

- `peptide`: Modified peptide sequence.
- `stripped_peptide`: Unmodified amino-acid sequence.
- `charge`: Precursor charge, or null when charge states were combined.
- `proteins`: Protein assignments.
- `is_decoy`: Whether the LFQ precursor is a decoy.
- `q_value`: Precursor-level q-value assigned by picked target-decoy competition.
- `score`: Cross-run LFQ peak score used for precursor-level competition.
- `spectral_angle`: Intensity-weighted normalized isotope-pattern spectral angle for the selected cross-run peak.
- `filename`: Acquisition file represented by this row.
- `intensity`: Integrated MS1 signal. A missing signal is a Parquet null, never a numeric zero sentinel.
- `ms2_confirmed`: Boolean indicating direct accepted MS2 identification evidence for this precursor in this file. `false` does not mean the intensity used a different quantification algorithm; all LFQ intensities use the same cross-run workflow.

Sage does not report a `missing_reason`: it cannot reliably distinguish biological absence from detection-limit, alignment, extraction, or scoring causes for a null intensity.
