# Sage Plus scientific pilot

The released beta.3 binary and upstream Sage commit df9219951cc9a54cf4cd55d76541af24b687bd3d are frozen. These measurements assess the selected public pilot and do not certify production calibration.

Recorded jobs: 123 complete, 21 failed or invalid, 0 running.

Host CPU: Intel(R) Core(TM) i7-10700K CPU @ 3.80GHz. Logical processors reported: 16. Search worker limits are recorded per job.

## Paired engineering comparisons

Counts come from the underlying PSM files. These local HEK timing observations can overlap acquisition work and must not be presented as controlled final timing experiments.

| Workload | Engine | Measured trials | Median seconds | Median peak MiB | Median target PSMs at 1% | Median peptidoforms at 1% |
| --- | --- | ---: | ---: | ---: | ---: | ---: |
| common-mods | plus | 3 | 15.06 | 3633.5 | 2270 | 1463 |
| common-mods | upstream | 3 | 12.47 | 5348.5 | 2266 | 1460 |
| standard | plus | 3 | 6.41 | 1502.6 | 2203 | 1409 |
| standard | upstream | 3 | 5.71 | 1903.9 | 2214 | 1409 |

## Public identification agreement

Matched files use identical search configurations. Overlap is measured among target PSMs passing each engine's reported 1% spectrum q-value threshold.

The original mixed reference contains seven E. coli proteins with unknown X residues and is rejected by beta.3. The v2 comparison, mixed-reference entrapment, quantification and repeated timing use a documented variant excluding those same seven proteins for both engines. Their known subsequences are also excluded. The original reference and failed attempts are retained.

| Input pair | Shared PSMs | Upstream only | Plus only | Jaccard | Status |
| --- | ---: | ---: | ---: | ---: | --- |
| public-comparison/PXD001468-0 | 16752 | 114 | 165 | 0.9836 | complete |
| public-comparison/PXD001468-1 | 19291 | 115 | 180 | 0.9849 | complete |
| public-comparison/PXD028735-0 | unavailable | unavailable | unavailable | unavailable | incomplete |
| public-comparison/PXD028735-1 | unavailable | unavailable | unavailable | unavailable | incomplete |
| public-comparison/PXD028735-2 | unavailable | unavailable | unavailable | unavailable | incomplete |
| public-comparison/PXD028735-3 | unavailable | unavailable | unavailable | unavailable | incomplete |
| public-comparison-v2/PXD028735-0 | 54639 | 818 | 973 | 0.9683 | complete |
| public-comparison-v2/PXD028735-1 | 55090 | 850 | 954 | 0.9683 | complete |
| public-comparison-v2/PXD028735-2 | 57416 | 921 | 1092 | 0.9661 | complete |
| public-comparison-v2/PXD028735-3 | 56282 | 892 | 1039 | 0.9668 | complete |

These are single paired public-file observations. They do not estimate timing variability.

| Public search | Seconds | Peak MiB | Target PSMs at 1% | Peptidoforms at 1% |
| --- | ---: | ---: | ---: | ---: |
| public-comparison/PXD001468-0-plus | 9.70 | 1500.6 | 16917 | 9293 |
| public-comparison/PXD001468-0-upstream | 10.02 | 2121.5 | 16866 | 9264 |
| public-comparison/PXD001468-1-plus | 9.39 | 1494.6 | 19471 | 11103 |
| public-comparison/PXD001468-1-upstream | 9.80 | 2087.3 | 19406 | 11084 |
| public-comparison/PXD028735-0-upstream | 26.46 | 3055.2 | 55453 | 39005 |
| public-comparison/PXD028735-1-upstream | 27.59 | 3056.2 | 55946 | 38917 |
| public-comparison/PXD028735-2-upstream | 30.89 | 3103.4 | 58333 | 40044 |
| public-comparison/PXD028735-3-upstream | 31.18 | 3143.9 | 57184 | 39161 |
| public-comparison-v2/PXD028735-0-plus | 24.94 | 2411.9 | 55612 | 39117 |
| public-comparison-v2/PXD028735-0-upstream | 26.15 | 3079.6 | 55457 | 39008 |
| public-comparison-v2/PXD028735-1-plus | 25.39 | 2446.3 | 56044 | 38982 |
| public-comparison-v2/PXD028735-1-upstream | 26.93 | 3056.4 | 55940 | 38925 |
| public-comparison-v2/PXD028735-2-plus | 27.29 | 2509.5 | 58508 | 40117 |
| public-comparison-v2/PXD028735-2-upstream | 29.08 | 3125.3 | 58337 | 40047 |
| public-comparison-v2/PXD028735-3-plus | 26.84 | 2468.1 | 57321 | 39271 |
| public-comparison-v2/PXD028735-3-upstream | 29.21 | 3164.1 | 57174 | 39147 |

## Repeated public timing

One selected file per study, one warmup and three measured trials per engine. These runs follow acquisition, conversion and the primary search matrices. Parentheses show the observed minimum and maximum, not a confidence interval.

| Study | Engine | Trials | Median seconds (range) | Median peak MiB | Median target PSMs at 1% |
| --- | --- | ---: | ---: | ---: | ---: |
| PXD001468 | plus | 3 | 9.48 (9.43, 9.50) | 1491.8 | 16917 |
| PXD001468 | upstream | 3 | 9.91 (9.91, 9.94) | 2134.1 | 16866 |
| PXD028735 | plus | 3 | 25.03 (24.95, 25.38) | 2419.9 | 55612 |
| PXD028735 | upstream | 3 | 26.05 (25.80, 26.79) | 3071.8 | 55457 |

## Independent entrapment pilot

The table reports each file and construction seed. FDP estimators describe realized discovery sets. A small estimate or similarity between engines alone is not proof of valid FDR control.

| Experiment | File and engine | Target peptides | Entrapments | Conservative paired FDP at nominal 1% | Seconds | Peak MiB |
| --- | --- | ---: | ---: | ---: | ---: | ---: |
| entrapment-human-20260914 | file-0-plus | 9221 | 57 | 1.23% | 24.05 | 3063.9 |
| entrapment-human-20260914 | file-0-upstream | 9191 | 54 | 1.17% | 22.09 | 3957.8 |
| entrapment-human-20260914 | file-1-plus | 10920 | 47 | 0.84% | 23.16 | 3164.3 |
| entrapment-human-20260914 | file-1-upstream | 10912 | 54 | 0.97% | 19.55 | 3772.5 |
| entrapment-human-20260915 | file-0-plus | 9208 | 48 | 1.05% | 22.81 | 3139.7 |
| entrapment-human-20260915 | file-0-upstream | 9190 | 53 | 1.16% | 20.08 | 3749.3 |
| entrapment-human-20260915 | file-1-plus | 10897 | 57 | 1.04% | 23.29 | 3157.8 |
| entrapment-human-20260915 | file-1-upstream | 10876 | 56 | 1.03% | 19.05 | 4025.8 |
| entrapment-human-20260916 | file-0-plus | 9211 | 45 | 0.95% | 21.68 | 3195.6 |
| entrapment-human-20260916 | file-0-upstream | 9182 | 43 | 0.91% | 18.37 | 4011.2 |
| entrapment-human-20260916 | file-1-plus | 10949 | 49 | 0.88% | 21.07 | 3085.9 |
| entrapment-human-20260916 | file-1-upstream | 10927 | 48 | 0.87% | 17.38 | 3775.1 |
| entrapment-hye-20260914 | file-0-plus | 37348 | 182 | 0.95% | 48.83 | 4315.2 |
| entrapment-hye-20260914 | file-0-upstream | 37206 | 171 | 0.89% | 39.30 | 5279.6 |
| entrapment-hye-20260914 | file-1-plus | 37126 | 185 | 0.97% | 47.95 | 4348.1 |
| entrapment-hye-20260914 | file-1-upstream | 37082 | 196 | 1.04% | 40.26 | 5313.8 |
| entrapment-hye-20260915 | file-0-plus | 37374 | 224 | 1.16% | 47.78 | 4335.9 |
| entrapment-hye-20260915 | file-0-upstream | 37225 | 210 | 1.10% | 39.52 | 5345.0 |
| entrapment-hye-20260915 | file-1-plus | 37170 | 213 | 1.12% | 49.79 | 4339.4 |
| entrapment-hye-20260915 | file-1-upstream | 37112 | 199 | 1.05% | 42.62 | 5363.2 |
| entrapment-hye-20260916 | file-0-plus | 37471 | 205 | 1.06% | 49.02 | 4328.1 |
| entrapment-hye-20260916 | file-0-upstream | 37264 | 198 | 1.04% | 41.83 | 5373.4 |
| entrapment-hye-20260916 | file-1-plus | 37259 | 206 | 1.08% | 47.89 | 4346.2 |
| entrapment-hye-20260916 | file-1-upstream | 37105 | 202 | 1.06% | 41.39 | 5316.6 |

Conditional resampling summaries:

- human: mean candidate-minus-upstream difference -0.019 percentage points, descriptive interval [-0.083, 0.054]. Two files and three shared seeds do not support broad generalization.
- hye: mean candidate-minus-upstream difference 0.028 percentage points, descriptive interval [-0.039, 0.066]. Two files and three shared seeds do not support broad generalization.

## Exact prefilter

Exact PSM, score and q-value equality passed in 12 of 12 completed measured pairs. Missing or failed pairs are not counted as passes.

## PTM synthesis consistency

These are site events passing reported 1% PSM and localization thresholds. The synthesis database is restricted, file-to-library mapping is not independently established, and the searches trigger heuristic rescoring. The table explicitly retains the peptide-level acceptance limitation.

| Search | Consistent sites | Inconsistent sites | Inconsistent fraction | Site rows with peptide q = 1 |
| --- | ---: | ---: | ---: | ---: |
| library-1-plus-all | 1172 | 118 | 9.15% | 100.00% |
| library-2-plus-all | 880 | 15 | 1.68% | 100.00% |
| library-1-plus-oracle-0 | 0 | 0 | undefined | undefined |
| library-1-plus-oracle-100 | 1173 | 106 | 8.29% | 100.00% |
| library-1-plus-oracle-25 | 461 | 80 | 14.79% | 100.00% |
| library-1-plus-oracle-50 | 761 | 92 | 10.79% | 100.00% |
| library-2-plus-oracle-0 | 0 | 0 | undefined | undefined |
| library-2-plus-oracle-100 | 880 | 15 | 1.68% | 100.00% |
| library-2-plus-oracle-25 | 414 | 22 | 5.05% | 100.00% |
| library-2-plus-oracle-50 | 659 | 19 | 2.80% | 100.00% |

No site events survive the joint 1% PSM, peptide and localization thresholds in this restricted synthetic pilot. Error among that empty accepted set is undefined.

Oracle coverage is derived from synthesis truth. It measures sensitivity to missing known sites and cannot demonstrate the benefit of an independently annotated biological site library. The site library constrains peptide generation, but the released localizer reconsiders all residue-compatible positions. Full oracle coverage therefore does not force synthesis-consistent localization.

## Known-ratio quantification

B/A ground truth is human 1, yeast 0.5 and E. coli 4. Ratios pair Alpha and Beta preparations without imputation. CVs measure preparation variability, not technical-repeat precision.

| Search | Species | Ratio pairs | Median log2 bias | Median absolute log2 error | Median preparation CV | Missing fraction |
| --- | --- | ---: | ---: | ---: | ---: | ---: |
| mbr-0 | ecoli | 0 | undefined | undefined | undefined | undefined |
| mbr-0 | human | 0 | undefined | undefined | undefined | undefined |
| mbr-0 | yeast | 0 | undefined | undefined | undefined | undefined |
| mbr-1 | ecoli | 1826 | 0.189 | 0.447 | 0.275 | 0.000547 |
| mbr-1 | human | 35230 | -0.045 | 0.235 | 0.297 | 0.008838 |
| mbr-1 | yeast | 10439 | -0.141 | 0.272 | 0.259 | 0.000048 |

Pure-human control, mbr-0: 0 quantified foreign-species rows. Among 0 quantified rows without direct MS2 confirmation, 0 are assigned exclusively to foreign species. I/L-indistinguishable matches to human are excluded. This measures the detectable foreign component of transfer error. Sample purity and reference completeness remain assumptions.
An independent check requires both PSM and peptide q-values at most 1% for direct MS2 evidence: 0 foreign-species rows among 0 control rows lacking that evidence. The engine's exported MS2 flag uses its LFQ peptide threshold without a separate PSM q-value requirement, so the two counts need not agree.


Pure-human control, mbr-1: 5872 quantified foreign-species rows. Among 8349 quantified rows without direct MS2 confirmation, 5767 are assigned exclusively to foreign species. I/L-indistinguishable matches to human are excluded. This measures the detectable foreign component of transfer error. Sample purity and reference completeness remain assumptions.
An independent check requires both PSM and peptide q-values at most 1% for direct MS2 evidence: 5836 foreign-species rows among 8765 control rows lacking that evidence. The engine's exported MS2 flag uses its LFQ peptide threshold without a separate PSM q-value requirement, so the two counts need not agree.


The serializer repeats one precursor-peak q-value across all files. It does not report an individual transfer q-value. Foreign assignments in the control therefore cannot be interpreted as a direct calibration test at the global precursor unit. The engine's discovery counter uses 5%, while this pilot retains its primary 1% LFQ endpoint. The following counts expose that threshold distinction before species filtering.

| Search | LFQ threshold | Target precursors | Positive target file rows |
| --- | ---: | ---: | ---: |
| mbr-0 | 1% | 0 | 0 |
| mbr-0 | 5% | 42242 | 137604 |
| mbr-1 | 1% | 24510 | 121179 |
| mbr-1 | 5% | 46846 | 231930 |

### Exploratory ratio diagnostic before LFQ filtering

This diagnostic was added after observing that MBR off accepts no features at the primary 1% LFQ threshold. It uses positive intensities supported by direct MS2 evidence passing both 1% PSM and peptide thresholds, with cross-species I/L ambiguity excluded. It deliberately omits the LFQ q-value filter and does not change the primary result.

| Search | Species | Ratio pairs | Median log2 bias | Median absolute log2 error |
| --- | --- | ---: | ---: | ---: |
| mbr-0 | ecoli | 1186 | 0.043 | 0.361 |
| mbr-0 | human | 42352 | -0.052 | 0.245 |
| mbr-0 | yeast | 10517 | -0.106 | 0.274 |
| mbr-1 | ecoli | 1229 | 0.078 | 0.345 |
| mbr-1 | human | 45122 | -0.035 | 0.253 |
| mbr-1 | yeast | 11079 | -0.106 | 0.266 |

## Failed and invalid jobs

Original rejected configurations and resource failures are retained. Configuration repairs use separate output directories.

The broad-PTM stress case exceeds the selected resource budgets. Plus estimates a 23.8 GiB additional modified-database peak against its 10 GiB guard, while upstream fails allocation under the 16 GiB address-space ceiling. Initial oracle configurations omit the required max_count field and are corrected in ptm-oracle-v2. The 2 GiB address-space stress case fails with the prefilter off. Original mixed-reference failures are the undefined-residue case described above.

- local-paired/broad-ptm-plus-0: failed.
- local-paired/broad-ptm-plus-1: failed.
- local-paired/broad-ptm-plus-2: failed.
- local-paired/broad-ptm-plus-3: failed.
- local-paired/broad-ptm-upstream-0: failed.
- local-paired/broad-ptm-upstream-1: failed.
- local-paired/broad-ptm-upstream-2: failed.
- local-paired/broad-ptm-upstream-3: failed.
- ptm/library-1-plus-oracle-0: failed.
- ptm/library-1-plus-oracle-100: failed.
- ptm/library-1-plus-oracle-25: failed.
- ptm/library-1-plus-oracle-50: failed.
- ptm/library-2-plus-oracle-0: failed.
- ptm/library-2-plus-oracle-100: failed.
- ptm/library-2-plus-oracle-25: failed.
- ptm/library-2-plus-oracle-50: failed.
- public-comparison/PXD028735-0-plus: failed.
- public-comparison/PXD028735-1-plus: failed.
- public-comparison/PXD028735-2-plus: failed.
- public-comparison/PXD028735-3-plus: failed.
- scaling/address-space-2-prefilter-0: failed.

## Follow-up priorities

- Establish and validate confidence for individual MBR transfers. Keep precursor-level confidence and transfer-level confidence distinct in outputs and manuscript claims.
- Resolve the restricted PTM pilot's peptide acceptance failure, independently establish file-to-library truth, and test localization separately from identification and oracle candidate generation.
- Freeze a separate evaluation after the pilot, including held-out learned-model comparisons and a sample-size justification. Preserve the enlarged-reference runtime costs and explicit reference exclusions in any broad performance or coverage claims.

## Reproducibility and remaining scope

- Pilot selection and finite precision do not establish production calibration.
- Local HEK origin remains unverified and local timing may overlap acquisition or conversion work.
- Timing comes from one Linux workstation and does not establish performance across platforms or hardware.
- Seven E. coli proteins with undefined X residues are excluded from the repaired mixed-reference experiments.
- PTM site-library controls use synthesis truth and are oracle sensitivity experiments.
- PTM synthesis consistency includes identification and localization error and is not automatically an arrangement-level FLR estimate.
- Preparation variability, technical injections and entrapment seeds are not biological replication.
- Contaminant coverage and pure-species sample purity need validation before a final FDR or false-transfer claim.

The base RT and final LDA fits reuse scored observations. Grouped additive PTM-offset folds do not make the entire pipeline cross-fitted. This requires a separate model-validation experiment before strong learned-model claims.

Machine-readable source SHA-256: `c2028d0ebf86a1eb3dba0b208f144cd290ba443dc63bb25aa00e4d71fb717379`.
The evidence bundle includes configs, result files, reference sequences, binary hashes, source snapshots, extraction scripts and source receipts. Raw and converted spectra remain on /data. No external deposition has been performed.

## Figures

- [entrapment pilot](figures/entrapment-pilot.svg)
- [ptm pilot](figures/ptm-pilot.svg)
- [public timing pilot](figures/public-timing-pilot.svg)
- [quantification diagnostic](figures/quantification-diagnostic.svg)
- [quantification pilot](figures/quantification-pilot.svg)
- [scaling pilot](figures/scaling-pilot.svg)

Captions and image checksums are recorded in [figures.json](figures/figures.json).
