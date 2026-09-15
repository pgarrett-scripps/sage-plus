# Scientific pilot protocol

The September 14, 2026 pilot runs all large acquisitions and search outputs under
`/data/sage-plus-scientific/20260914`. The released beta.3 executable remains
frozen. No scoring defaults or earlier result values are changed by this work.

## Experiments

| Track | Pilot evidence | Interpretation boundary |
| --- | --- | --- |
| Independent error validation | Full-reference FDRBench searches on public HEK and mixed-species DDA files, three construction seeds, both engines | Conditional pilot estimates. Two files and three seeds do not establish universal calibration |
| Fair comparisons | Common extraction of target PSMs, peptidoforms and sequences from upstream TSV and beta.3 Parquet, nominal thresholds from 0.1% to 5% | Identification counts are separate from independently estimated error |
| Scaling and ablations | One, two, four and eight threads, exact prefilter off/on, explicit virtual-memory limits | Exact-prefilter comparisons isolate that setting. Historical compact-storage experiments remain separately labeled |
| PTM and quantification | Synthesis-defined phosphosites, oracle site-coverage controls, known-ratio mixtures, MBR off/on and a pure-human control | Oracle sites do not demonstrate benefit from prior biological annotations. Site-event error includes identification error |
| Reproducibility | Input receipts, pinned binaries, conversion records, frozen plans, raw results, independent count checks and verified archive members | The local archive does not imply publication, external deposition or final scientific acceptance |

The source of public inputs is PRIDE. File selection precedes the corresponding
experimental outcomes. Remaining acquisitions are reserved for a later frozen
evaluation. The main selection is 16.44 GiB. The added pure-human control brings
spectra and reference-library acquisition to 19.89 GiB, within the initial
20 GiB budget. Supplementary design documents and reference proteomes are small
additional inputs. Search processes run sequentially with eight or fewer Rayon
threads, a 10 GiB Sage Plus memory guard and a common 16 GiB address-space ceiling.
Explicit memory-sensitivity jobs use smaller address-space ceilings.
An additional timing plan freezes one selected file from each public study
before the corresponding public comparison outcomes. It runs after acquisition,
conversion and the primary search matrices, using one warmup and three measured
trials per engine with alternating engine order. The other public-file timings
remain single observations.

## Sources and reference assumptions

- PXD001468: the first named A and B HEK fraction acquisitions, distributed as
  mzXML. These are two files from one study, not independent biological studies.
  MS2 arrays are converted to MGF without centroiding, with charge, precursor,
  retention time and array-length checks. MS1 data is deliberately absent from
  this identification-only derivative.
- PXD028735: first injections of the Alpha and Beta mixture preparations for
  conditions A and B. A pure-human acquisition provides an absence control.
  The primary methods specify B/A ratios of 1 for human, 0.5 for yeast and 4 for
  E. coli. [Dataset paper](https://doi.org/10.1038/s41597-022-01216-6).
  Species assignments are checked against all three reference proteomes after
  I/L normalization. Any occurrence in multiple species excludes that peptide
  from ratio analysis, even if only one species appears in the engine's assignment.
  This conservative check includes occurrences outside the configured digestion.
- The community SDRF labels yeast as Candida albicans. The primary methods name
  Promega V7461. The manufacturer's TM410 identifies that reagent as
  Saccharomyces cerevisiae, which determines this pilot's reference choice.
  [Manufacturer documentation](https://www.promega.com/-/media/files/resources/protocols/technical-manuals/101/ms-compatible-yeast-and-human-protein-extracts-protocol.pdf).
- References are the reviewed canonical human set, yeast and E. coli K12
  proteomes from UniProt release 2026_03, with response headers and SHA-256
  receipts. The mixed reference includes the stated iRT standards. E. coli
  strain identity, incidental contaminants and the coverage of human isoforms
  remain limitations to audit before final evaluation.
  The original combined reference has 30,897 proteins. Beta.3 rejects seven
  E. coli entries containing undefined `X` residues. The repaired mixed-data
  experiments use `hye-irt-defined.fasta`, retaining 30,890 proteins and listing
  every excluded accession in its reference manifest. Both engines use this
  same variant. Entire affected proteins, including their known subsequences,
  are excluded. The original reference, compatibility failures and plans remain
  in the evidence. This exception qualifies claims of reference completeness.
- PXD000138: the first two distributed HCD MGF files. Supplementary Table 2 and
  the deposited FASTA reconstruct 96 synthetic libraries and 101,520 sequences
  with synthesis-defined sites. Each library's seed, length and expected variant
  count are checked. The search database contains the complete synthesized
  sequence set, without the original IPI background. File-to-library assignment
  is not yet independently established. [Synthesis paper](https://doi.org/10.1038/nbt.2585).

RAW conversion uses the official self-contained ThermoRawFileParser
`v.2.0.0-dev` Linux distribution. This is a development release, selected because
its runtime is self-contained. The entire converter directory is hashed in each
conversion record. Conversion enables the vendor's default peak picking and
retains MS1 and MS2 in mzML. Converted data is identical between engines.

## Estimator audit

FDRBench 1.1.1 is paired with its tagged source commit
`3b619a9acf60d7292fb651a00da55f58cb67fb79`. Its ranking uses peptide q-value
first and the supplied score second. The new extractor uses a single actual
best-discriminant PSM per peptidoform and verifies consistent peptide q-values.
It does not combine a minimum q-value with a maximum hyperscore from unrelated
rows. Static-only entrapment searches require exactly one peptidoform per sequence.

The independent audit reproduces target counts, entrapment counts and the
combined estimator at every prespecified threshold. It also bounds the paired
estimator over exact target-partner score ties, then checks that the official
output falls within that interval. Empty denominators remain undefined. Failed
or absent runs cannot become zero-discovery runs.

The resampling analysis keeps the engines paired and resamples shared
construction seeds jointly across files. Its intervals are descriptive within
the selected study. Individual PSMs are never treated as independent replicates.

## Learned-model audit

The released pipeline selects provisional targets, fits mass-error alignment and
base retention-time predictions, then fits LDA and posterior-error models to the
search results. The base retention-time fit and LDA are not fully held out from
their scored observations. Additive PTM retention offsets have grouped folds,
but those folds do not make the whole pipeline cross-fitted. This is a documented
information-reuse boundary, not proof that reported FDR is invalid.

The synthetic pilot triggers the existing heuristic rescoring fallback. Its
peptide-level q-values and localization outputs must therefore be reported
together. Agreement with a synthesis site does not rescue an unaccepted peptide,
and a nominal localization threshold does not by itself certify empirical error.

The oracle site library constrains peptide generation. The released localizer
then considers all residue-compatible positions for the modification mass, using
`potential_mods` rather than the protein site library. Consequently, 100% oracle
coverage does not force the exported localization to the synthesis site.
The coverage experiment assesses this complete behavior and must not be described
as localization constrained to known sites.

The primary LFQ endpoint retains the 1% LFQ q-value threshold. After observing
that the MBR-off run accepts no features at that threshold, an exploratory
diagnostic measures ratios before the LFQ confidence filter, requiring direct
MS2 evidence passing both 1% PSM and peptide thresholds. This diagnostic does
not replace the primary endpoint. The engine's exported MS2 flag applies the
LFQ peptide threshold without an additional PSM q-value requirement, so the
absence-control analysis also reports an independently checked joint threshold.
The exported LFQ q-value belongs to a precursor peak and is repeated across
files. It does not quantify confidence in each transfer. The engine's internal
discovery counter uses 5%, so the report also shows yields at that threshold
while preserving the primary 1% analysis. These reporting units must remain
distinct in the benchmark report.

## Running and auditing

The scripts are standalone benchmark tooling, with standard-library Python and
the pinned DuckDB CLI for result extraction. The synthesis-design reader uses
the bundled Python runtime with openpyxl. Preserve its recorded version.

```bash
python3 benchmarks/acquire_scientific.py --root /data/sage-plus-scientific/20260914
python3 benchmarks/run_scientific.py PLAN.json --output OUTPUT_DIRECTORY
python3 benchmarks/scientific_entrapment.py --reference REFERENCE.fasta --spectra INPUT.mgf --output OUTPUT_DIRECTORY --seed 20260914
python3 benchmarks/analyze_scientific.py --root /data/sage-plus-scientific/20260914 --output /data/sage-plus-scientific/20260914/pilot-summary.json
python3 -m unittest discover -s benchmarks/tests -v
```

`continue_scientific.py` completes the selected pilot sequentially after the
initial scaling matrix. Its stage ledger retains every nonzero exit status.
It is specific to this prespecified pilot and is not a general scheduler.
`--resume-after-reference-repair` resumes after the original public matrix and
uses a separate stage ledger and comparison directory. It also runs the amended
public timing plan. The orchestration interruption and reference amendment are
recorded under `continuation/reference-repair.json`.

`bundle_scientific.py` creates a local evidence archive, includes the released
and upstream source snapshots, and reopens the archive to hash every member.
Large raw and converted spectra remain under `/data` with their receipts and
hashes. The archive retains the results needed to regenerate the figures.
`audit_scientific.py` first rechecks recorded search inputs, conversion inputs,
outputs and acquisition hashes. `report_scientific.py` and `plot_scientific.py`
derive the report and PNG/SVG figures from the same machine-readable summary.
The plotting environment uses Matplotlib 3.10.7 with an installation receipt
under `tools/plot-install.json`, including resolved dependency download hashes.
The archive retains the DuckDB executable. Analysis can replay extracted results
at a new root with `--duckdb PATH_TO_RETAINED_DUCKDB`, without rerunning searches.

Pilot configuration repairs receive new output directories. Original failed
attempts remain visible. Successful runs are reusable only when their input
identities and retained output hashes still match. Do not overwrite beta.2 paper
measurements with beta.3 values or promote these pilot results into final claims
without resolving the stated limitations and freezing a separate evaluation.
