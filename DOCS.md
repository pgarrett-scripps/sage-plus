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
    "bucket_size": 32768,           // How many fragments are in each internal mass bucket
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
    "static_mods": {        // Static modification masses or structured objects
      "^": 304.207,         // Apply static modification to N-terminus of peptide
      "K": 304.207,         // Apply static modification to lysine
      "C": {"mass": 57.0215, "name": "Carbamidomethyl"}
    },
    "variable_mods": {    // Variable modification masses or structured objects
      "M": [{             // Variable mods are applied *before* static mod
        "mass": 15.9949,
        "max_count": 1,
        "name": "Oxidation",
        "neutral_losses": [17.0265],
        "neutral_loss_mode": "optional"
      }],
      "K": [{"mass": 42.0106, "max_count": 1}, 14.0157],
      "^Q": [-17.026549],
      "^E": [-18.010565], // Applied to N-terminal glutamic acid
      "$": [49.2, 22.9],  // Applied to peptide C-terminus
      "[": [42.0],          // Applied to protein N-terminus
      "]": [111.0]          // Applied to protein C-terminus
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
  "deisotope": false,       // Optional[bool] {default=false}: perform deisotoping and charge state deconvolution
  "chimera": false,         // Optional[bool] {default=false}: search for chimeric/co-fragmenting PSMS
  "wide_window": false,     // Optional[bool] {default=false}: _ignore_ `precursor_tol` and search in wide-window/DIA mode
  "predict_rt": false,    // Optional[bool] {default=true}: use retention time prediction model as a feature for LDA
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
$ docker pull ghcr.io/pgarrett-scripps/sage-plus:v0.1.0-beta.1
$ docker run -it --rm -v ${PWD}:/data ghcr.io/pgarrett-scripps/sage-plus:v0.1.0-beta.1 sage -o /data /data/config.json
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

- **bucket_size**: Integer. The number of fragments in each internal mass bucket (default: 8192). Tweaking this parameter can increase search performance for wide precursor or fragment searches.

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

#### Static Modifications

- **static_mods**: Dictionary with characters as keys and bare masses or structured modification objects. Represents static modifications applied to amino acids or termini (default: {}). Static modifications are applied after variable modifications.
  - Example: Apply a static modification of 304.207 to the N-terminus of the peptide and lysine, and 57.0215 to cysteine.
    ```json
    "database": {
      "static_mods": {
        "^": 304.207,
        "K": 304.207,
        "C": {"mass": 57.0215, "name": "Carbamidomethyl"}
      }
    }
    ```

#### Variable Modifications

- **max_variable_mods**: Integer. Limit the total variable modifications on each peptide (default: 2).
- **max_total_variable_mods**: Integer. Limit exhaustive and PTM-library-supported variable modifications combined. Defaults to `max_variable_mods` and cannot be lower.
- **max_combinations**: Integer. Optional hard cap on the total variants generated per input peptide, including the unmodified form. Variants with fewer modifications are retained first. Values below 1 are treated as 1 (default: unlimited).
- **static_mods** and **variable_mods** accept existing bare numeric masses or structured objects. Structured objects require `mass` and may contain `name`, `neutral_losses`, and `neutral_loss_mode`; variable modifications may additionally contain `max_count` and `site_mode`.
- **name**: Optional display label. Named peptide modifications render as `[Name]`; unnamed and legacy entries retain numeric mass rendering.
- **neutral_losses**: Optional list of positive neutral-loss masses. During full scoring, retained and loss forms from the same cleavage and charge are alternatives and contribute at most one match. When multiple applicable modified sites occur in one fragment, their allowed loss choices are combined and duplicate total losses are removed.
- **neutral_loss_mode**: Either `"optional"` (default) or `"required"`. Optional mode generates the retained fragment plus configured losses. Required mode suppresses the retained form for fragments containing the modification and requires at least one configured neutral loss. Preliminary indexing retains one canonical form per cleavage to avoid favoring modifications with more configured fragment alternatives.
- **site_mode**: Either `"exhaustive"` (default), `"library"`, or `"both"`. Library-backed modes require `name` and `max_count`. Entries with the same name across residue specificities are one logical modification and must have identical mass, limit, name, and neutral-loss settings.
  - Example: Apply a variable modification of 15.9949 to methionine, 49.2022 to the C-terminus of the peptide, 42.0 to the N-terminus of the protein, and 111.0 to the C-terminus of the protein.
    ```jsonc
    "database": {
      "variable_mods": {
        "M": [{
          "mass": 15.9949,
          "name": "Oxidation",
          "neutral_losses": [17.0265],
          "neutral_loss_mode": "optional"
        }],
        "K": [{"mass": 42.0106, "max_count": 1, "name": "Acetyl"}, 14.0157],
        "^Q": [-17.026549],
        "^E": [-18.010565], // Applied to N-terminal glutamic acid
        "$": [49.2022],     // Applied to peptide C-terminus
        "[": 42.0,          // Applied to protein N-terminus
        "]": 111.0          // Applied to protein C-terminus
      }
    }
    ```

#### Modification channels

Structured static and variable modifications may define `channel_offsets`. The effective mass is
the modification's `mass` plus the selected channel offset. Static versus variable placement keeps
its normal meaning, while every channel-aware modification on a peptide resolves to one coherent
channel.

```json
"database": {
  "static_mods": {
    "C": {"mass": 57.021464, "name": "Carbamidomethyl"},
    "K": {
      "mass": 0.0,
      "name": "SILAC-K",
      "channel_offsets": {"light": 0.0, "medium": 4.025107, "heavy": 8.014199}
    },
    "R": {
      "mass": 0.0,
      "name": "SILAC-R",
      "channel_offsets": {"light": 0.0, "medium": 6.020129, "heavy": 10.008269}
    }
  },
  "variable_mods": {
    "M": [{"mass": 15.994915, "name": "Oxidation"}]
  }
}
```

The same field is valid on a variable modification. For example, an optional two-channel lysine
label can be written as `{"mass": 0.0, "name": "Optional-Lys8", "max_count": 2,
"channel_offsets": {"light": 0.0, "heavy": 8.014199}}`.

All `channel_offsets` dictionaries must contain exactly the same channel names and at least two
chemically distinct channels. A unique channel whose offsets are all zero is inferred as the
reference channel. Variable channel modifications consume the existing variable-modification and
combination limits. Chemically identical zero-offset variants are searched once, while their full
set of channel partners is retained for LFQ extraction.

When `quant.lfq` is enabled, one identified channel can seed extraction of its configured channel
partners at their exact precursor masses. `lfq.parquet` schema version 2 records `label_channel`,
`label_group`, and `ratio_to_reference`. Without channel-aware modifications, Sage retains the
existing LFQ behavior and writes schema version 1.

Current label channels assume complete incorporation with fixed site mass shifts. Partial
incorporation, whole-proteome nitrogen labeling, NeuCode resolution, and isotope-purity correction
are not yet modeled.

#### PTM site libraries

A PTM library is a Parquet or TSV table containing observed locations only. Modification
masses, names, limits, and neutral losses remain defined in `variable_mods`. The format is
selected from the `.parquet`, `.tsv`, or `.tsv.gz` filename extension.

Required columns are `protein` (UTF-8), `position` (one-based integer), `residue`
(one-letter UTF-8), and `modification` (the exact configured modification name).
Additional evidence columns are allowed and ignored during database construction.
When `ptm_library` is configured, every variable modification must specify `max_count`;
library-referenced modifications must also have a unique, non-empty `name`.

```json
"database": {
  "variable_mods": {
    "S": [{
      "name": "Phospho",
      "mass": 79.966331,
      "max_count": 3,
      "site_mode": "both",
      "neutral_losses": [97.976896]
    }]
  },
  "max_variable_mods": 1,
  "max_total_variable_mods": 3,
  "max_combinations": 1000,
  "ptm_library": {
    "path": "discovery-sites.tsv",
    "strict": true
  }
}
```

`max_variable_mods` limits placements generated exhaustively. Library-supported
placements do not consume that budget, but do consume `max_total_variable_mods`, the
named modification's `max_count`, and `max_combinations`. All candidates are enumerated
together before decoy generation and indexing. With no library configured, existing
modification behavior is unchanged.

When PTM localization is enabled for a FASTA search, Sage also writes both
`results.sage.ptm-library.parquet` and `results.sage.ptm-library.tsv`. They contain the
passing localized sites whose modification names match the configured variable
modifications; either file can be supplied to a later search through `ptm_library.path`.
  - Syntax:
    "^X": Modification to be applied to amino acid X if it appears at the N-terminus of a peptide
    "$X": Modification to be applied to amino acid X if it appears at the C-terminus of a peptide
    "[X": Modification to be applied to amino acid X if it appears at the N-terminus of a protein
    "]X": Modification to be applied to amino acid X if it appears at the C-terminus of a protein

### Decoys

- **decoy_tag**: String. The tag used to identify decoy entries in the FASTA database (default: "rev_").
- **generate_decoys**: Boolean. If true, ignore decoys in the FASTA database matching `decoy_tag`, and generate internally reversed peptides (default: false).

### FASTA

- **fasta**: String. The path to the FASTA file, either a local path or s3 object URI.

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
      "rt_pct_tolerance": 0.5
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

- **deisotope**: Boolean. Perform deisotoping and charge state deconvolution on MS2 spectra (default: false). Recommended for high-resolution MS2 scans. This setting may interfere with TMT-MS2 quantification, use at your own risk.
- **chimera**: Boolean. Search for chimeric/co-fragmenting PSMs (default: false).
- **wide_window**: Boolean. Ignore `precursor_tol` and search spectra in wide-window/dynamic precursor tolerance mode (default: false).
- **predict_rt**: Boolean. Use retention time prediction model as a feature for LDA (default: false).
- **ion_mobility_model.enabled**: Boolean. Fit and use the ion-mobility model when mobility observations are present (default: true). Set this to `false` to keep observed mobility data without fitting predictions.
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

When `ptm_localization.enabled` is true, sage attempts to pinpoint which residue carries each variable modification on a confidently-identified peptide, analogous to MaxQuant's site tables or MSFragger/PTMProphet.

Example configuration:

```json
"ptm_localization": {
  "enabled": true,
  "psm_q_value": 0.01,
  "localization_q_value": 0.01
}
```

For each FDR-passing target PSM (spectrum q-value ≤ `ptm_localization.psm_q_value`), and for each distinct variable-modification delta mass it carries, sage:
1. recovers the candidate residues from the search's modification specificity rules (e.g. all S/T/Y for Phospho),
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
- All variable modifications are localized; terminal-specificity modifications (peptide/protein N- and C-term) are not relocated.
- Localization runs after spectrum FDR assignment and only for passing target PSMs. Sage re-reads MS2 spectra for this optional pass rather than retaining the full experiment in memory.
- `ptm_localization.psm_q_value` controls identification quality; `ptm_localization.localization_q_value` controls arrangement-level localization FLR. `localization_probability` remains a within-PSM marginal site probability.
- A modification without enough eligible impossible residues to construct a balanced decoy search space is not included in the FDR-controlled reports.
- FASTA searches preserve protein coordinates during indexing. The canonical PSM output attaches each protein accession to its one-based inclusive start and end positions plus the preceding and following amino acids. Pre-digested peptide TSV and spectral-library inputs omit coordinates when they are unavailable.

## Spectrum Paths

- **mzml_paths**: List of strings. Despite the legacy field name, Sage accepts mzML, mzMLb, MGF, Bruker TDF, and Thermo Fisher RAW inputs. mzML and MGF paths may be local or use a configured object-store URL. mzMLb, Thermo RAW, and Bruker TDF inputs must be local because their readers require seekable files. mzMLb support is optional and requires building with `--features mzmlb`. Files ending in ".gz" or ".gzip" are inferred to be compressed.
  - Thermo RAW input uses centroid peak lists directly. TMT signal-to-noise mode (`quant.tmt_settings.sn: true`) still requires mzML containing a noise array.
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
