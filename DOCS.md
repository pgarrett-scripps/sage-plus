# Sage Documentation

The most up-to-date documentation now lives [here](https://sage-docs.vercel.app/docs).


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

Sage is distributed as source code, and as a standalone executable file.

### Installing via conda

Sage can be installed from [bioconda](https://anaconda.org/bioconda/sage-proteomics):

```
$ conda install -c bioconda -c conda-forge sage-proteomics
$ sage --help
```

### Compiling the development version

1. Install the [Rust programming language compiler](https://rustup.rs/)
2. Download Sage source code via git: `git clone https://github.com/lazear/sage.git` or by [zip file](https://github.com/lazear/sage/archive/refs/heads/master.zip)
3. Compile: `cargo build --release`
4. Run: `./target/release/sage config.json`

Once you have Rust installed, you can copy and paste the following lines into your terminal to complete the above instructions, and run Sage on the example mzML provided in the repository (a single scan from PXD016766)

```sh
git clone https://github.com/lazear/sage.git
cd sage
cargo run --release tests/config.json 
```

### Downloading the latest release

1. Visit the [Releases](https://github.com/lazear/sage/releases/latest) website.
2. Download the correct pre-compiled binary for your operating system.
3. Run: `sage <path/to/config.json>`

### Interfacing with AWS S3

Sage is capable of natively reading & writing files to AWS S3:

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
      --parquet
          Write parquet files instead of tab-separated files
      --write-pin
          Write percolator-compatible `.pin` output files
      --max-memory <GiB>
          Abort if Sage's memory use exceeds this many GiB, to keep the system responsive
          (default: 90% of total RAM; 0 disables). Also settable via SAGE_MAX_MEMORY_GB.
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
- Record of search parameters (`results.json`) will be created that details input/output paths and all search parameters used for the search
- MS2 search results will be stored as a tab-separated file (`results.sage.tsv`) file - this is a tab-separated file, which can be opened in Excel/Pandas/etc
- MS2 and MS3 quantitation results will be stored as a tab-separated file (`tmt.tsv`, `lfq.tsv`) if `quant.tmt` or `quant.lfq` options are used in the parameter file

If `--parquet` is passed as a command line argument, `results.sage.parquet` (and optionally, `lfq.parquet`) will be written. These have a similar set of columns, but TMT values are stored as a nested array alongside PSM features

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
    "max_combinations": 8,  // Optional[int] {default=null} Limit total variants per peptide
    "decoy_tag": "rev_",    // Optional[str] {default="rev_"}: See notes above
    "generate_decoys": false, // Optional[bool] {default="true"}: Ignore decoys in FASTA database matching `decoy_tag`
    "fasta": "dual.fasta"   // str: mandatory path to FASTA file
  },
  "quant": {                // Optional - specify only if TMT or LFQ
    "tmt": "Tmt16",         // Optional[str] {default=null}, one of "Tmt6", "Tmt10", "Tmt11", "Tmt16", or "Tmt18"
    "tmt_settings": {
      "level": 3,           // Optional[int] {default=3}, MS-level to perform TMT quantification on
      "sn": false           // Optional[bool] {default=false}, use Signal/Noise instead of intensity for TMT quant. Requires noise values in mzML
    },
    "lfq": true,            // Optional[bool] {default=null}, perform label-free quantification
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
  "predict_rt": false,    // Optional[bool] {default=true}: use retention time prediction model as an feature for LDA
  "min_peaks": 15,          // Optional[int] {default=15}: only process MS2 spectra with at least N peaks
  "max_peaks": 150,         // Optional[int] {default=150}: take the top N most intense MS2 peaks to search,
  "min_matched_peaks": 6,   // Optional[int] {default=4}: minimum # of matched b+y ions to use for reporting PSMs
  "max_fragment_charge": 1, // Optional[int] {default=null}: maximum fragment ion charge states to consider,
  "report_psms": 1,         // Optional[int] {default=1}: number of PSMs to report for each spectra. Higher values might disrupt PSM rescoring.
  "percolator": {            // Optional: cross-fitted semi-supervised rescoring; LDA remains the default
    "enabled": false,
    "model": "svm",          // "svm" (default) or "perpetual"
    "budget": 0.3,
    "svm_c": 1.0,
    "svm_epochs": 100,
    "folds": 3,
    "iterations": 3,
    "train_fdr": 0.01,
    "min_positive_psms": 50,
    "min_decoy_psms": 50,
    "seed": 42
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
$ docker pull ghcr.io/lazear/sage:master
$ docker run -it --rm -v ${PWD}:/data ghcr.io/lazear/sage:master sage -o /data /data/config.json
# The sage executable is located in /app/sage in the image
```

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

### Percolator-style rescoring

Set `percolator.enabled` to `true` to replace LDA scoring with cross-fitted,
semi-supervised rescoring. `model` can be `svm` (the default linear SVM) or
`perpetual` (boosted trees). Each training fold fits its own mass-error KDE. If
training cannot produce enough positive or decoy PSMs, Sage logs a warning and
falls back to LDA. `svm_c` controls the SVM penalty and `budget` controls tree
complexity.

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
- **max_combinations**: Integer. Optional hard cap on the total variants generated per input peptide, including the unmodified form. Variants with fewer modifications are retained first. Values below 1 are treated as 1 (default: unlimited).
- **static_mods** and **variable_mods** accept existing bare numeric masses or structured objects. Structured objects require `mass` and may contain `name`, `neutral_losses`, and `neutral_loss_mode`; variable modifications may additionally contain `max_count`.
- **name**: Optional display label. Named peptide modifications render as `[Name]`; unnamed and legacy entries retain numeric mass rendering.
- **neutral_losses**: Optional list of positive neutral-loss masses. During full scoring, retained and loss forms from the same cleavage and charge are alternatives and contribute at most one match. When multiple applicable modified sites occur in one fragment, their allowed loss choices are combined and duplicate total losses are removed.
- **neutral_loss_mode**: Either `"optional"` (default) or `"required"`. Optional mode generates the retained fragment plus configured losses. Required mode suppresses the retained form for fragments containing the modification and requires at least one configured neutral loss. Preliminary indexing retains one canonical form per cleavage to avoid favoring modifications with more configured fragment alternatives.
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

The quant section is optional and should be specified only if TMT or LFQ is used.


- **tmt**: String. One of "Tmt6", "Tmt10", "Tmt11", "Tmt16", or "Tmt18" (default: null).
- **tmt_settings**: Object containing TMT-specific settings.
  - **level**: Integer. The MS-level to perform TMT quantification on (default: 3).
  - **sn**: Boolean. Use Signal/Noise instead of intensity for TMT quantification. Requires noise values in mzML (default: false).
- **lfq**: Boolean. Perform label-free quantification (default: null).
- **lfq_settings**: Object containing LFQ-specific settings.
  - **peak_scoring**: String. The method used for scoring peaks in LFQ, one of: "Hybrid", "RetentionTime", "SpectralAngle" (default: "Hybrid").
  - **integration**: String. The method used for integrating peak intensities, either "Sum" or "Max" (default: "Sum").
  - **spectral_angle**: Float. Threshold for the spectral angle similarity measure, ranging from 0 to 1 (default: 0.7).
  - **ppm_tolerance**: Float. Tolerance for matching MS1 ions in parts per million (default: 5.0).
  - **rt_pct_tolerance**: Float. Symmetric retention-time tolerance for match-between-runs, as a percentage of total gradient length (default: 0.5). For example, `0.5` searches +/-0.5% around the aligned retention time.

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

`predict_rt` is incompatible with `quant.lfq = true`. Setting `quant.lfq = true` will automatically turn on global retention time alignment and prediction, which are crucial for accurate direct ion current extraction.

- **deisotope**: Boolean. Perform deisotoping and charge state deconvolution on MS2 spectra (default: false). Recommended for high-resolution MS2 scans. This setting may interfere with TMT-MS2 quantification, use at your own risk.
- **chimera**: Boolean. Search for chimeric/co-fragmenting PSMs (default: false).
- **wide_window**: Boolean. Ignore `precursor_tol` and search spectra in wide-window/dynamic precursor tolerance mode (default: false).
- **predict_rt**: Boolean. Use retention time prediction model as a feature for LDA (default: false).
- **min_peaks**: Integer. Only process MS2 spectra with at least N peaks (default: 15).
- **max_peaks**: Integer. Take the top N most intense MS2 peaks to search (default: 150).
- **min_matched_peaks**: Integer. The minimum number of matched b+y ions to use for reporting PSMs (default: 4).
- **max_fragment_charge**: Integer. The maximum fragment ion charge states to consider (default: null - use precursor z-1).
- **report_psms**: Integer. The number of PSMs to report for each spectrum. Higher values might disrupt LDA (default: 1).
- **max_memory_gb**: Number. Abort the search if Sage's resident memory reaches this many GiB. Zero disables this limit (default: disabled).
- **min_free_memory_gb**: Number. Abort the search if system-available memory falls to this many GiB, preserving capacity for the operating system and other applications. Zero disables this limit (default: disabled).
- **batch_size**: Integer. Number of input files to load and search at once. Smaller values reduce temporary spectrum memory at the cost of throughput (default: half the number of CPUs, with a minimum of one). The `--batch-size` command-line option overrides this value.

When either memory limit is enabled, Sage estimates the unmodified digest, variable-modification expansion, and fragment/index sizes before allocating them. Unsafe database searches return an error before expansion begins. Estimates are conservative and are backed by a runtime memory monitor for allocations outside database construction.

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

Two site reports are written. They use TSV by default and Parquet when `--parquet` is selected:

- **results.sage.ptm-sites.tsv** or **results.sage.ptm-sites.parquet**: one row per localized modification site of each PSM. Columns include `peptide`, `modification`, `position` (1-based, within the peptide), `residue`, `localization_probability`, `delta_localization_score`, `target_decoy_score`, `localization_q_value`, `candidate_sites`, site-determining-ion counts, and `site_probabilities`.
- **results.sage.protein-sites.tsv** or **results.sage.protein-sites.parquet**: the best localization for each (protein, modified peptide site) aggregated across all supporting PSMs, including `best_localization_q_value`.

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
- Protein coordinates are not resolved (the FASTA is consumed during indexing), so the protein-site report uses peptide-relative positions attributed to each mapped protein.

## Spectrum Paths

- **mzml_paths**: List of strings. Despite the legacy field name, Sage accepts mzML, MGF, Bruker TDF, and Thermo Fisher RAW inputs. mzML and MGF paths may be local or use a configured object-store URL. Thermo RAW and Bruker TDF inputs must be local because their readers require seekable files. Files ending in ".gz" or ".gzip" are inferred to be compressed.
  - Thermo RAW input uses centroid peak lists directly. TMT signal-to-noise mode (`quant.tmt_settings.sn: true`) still requires mzML containing a noise array.
  - Example:
    ```json
    "mzml_paths": [
      "local/path.mzML",
      "local/path.raw",
      "s3://my-mass-spec-data/PXD0000001/foo.mzML.gz"
    ]
    ```
  
## Output directory:

- **output_directory**: Local directory, or S3 location where output files will be written. If the local directory does not already exist, it will be created. Write permissions are required for the directory or S3 path.
  - Possible output files are: "results.json", "results.sage.tsv", "lfq.tsv", "tmt.tsv", "results.sage.ptm-sites.tsv", and "results.sage.protein-sites.tsv". With `--parquet`, the Sage and PTM result tables use the corresponding `.parquet` names instead.
  - Example:
  ```json
  "output_directory": "s3://my-mass-spec-results/PXD003881/"
  ```

# Interpreting Sage Output

The "results.sage.tsv" file contains the following columns (headers):

- `peptide`: Peptide sequence, including modifications (e.g., NC\[+57.021\]HKGSFK).
- `proteins`: Proteins containing the peptide sequence.
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
