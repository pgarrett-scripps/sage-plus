# TMT known-ratio benchmark: HYpro16 (PXD020815), preliminary

Status: **draft, preliminary**. Candidate is Sage Plus v0.1.0-beta.12, which writes missing
reporter channels as 0.0. Beta 13 writes null. The analysis treats 0, null and NaN as missing,
so a Beta 13 rerun is directly comparable.

Rerun with a new binary (about 10 minutes: four Plus searches plus about 4 minutes of analysis; upstream runs are reused):

```
benchmarks/scientific-results/20260927-tmt/run.sh /path/to/sage b13
```

## Design

- **Sample:** HYpro16 (Navarrete-Perea, Gygi, Paulo, JASMS 2021 32:247, doi:10.1021/jasms.0c00299).
  - Human HCT116 peptides are 1:1 in all 16 TMTpro channels.
  - Yeast is 0 in 126, then 1/20, 1/10, 1/4 and 1/2 in triplicate, then 1x in 133N, 133C and 134N.
  - The map was checked against the paper text and the PRIDE SDRF: [channel-map-source.md](channel-map-source.md).
    The paper's Methods list prints only 15 entries, dropping one 1/2. The Results text ("single empty channel, followed by triplicates") and the SDRF resolve this.
- **Files:** Orbitrap Fusion Lumos + FAIMS, 90 min gradient, from the deposit's `2020_Hyp.xlsx`.
  - hrMS2: a18203, a18209, a18212.
  - SPS-MS3: a18207, a18210, a18213.
- **Engines:**
  - Upstream Sage `df92199` (0.15.0-beta.2), the same binary as the 20260914 LFQ pilot. SHA-256 `065f9c7b…402e`.
  - Sage Plus v0.1.0-beta.12 release binary (x86_64-unknown-linux-gnu, tag commit `54668e1`). SHA-256 `e11976c9a614e298e91c96881a4ba5ac7827d78e078bf3641358b04febc57c78`. It was verified against the release `SHA256SUMS`.
- **Spectra:**
  - mzML was converted by the pilot's ThermoRawFileParser build (`-f=1 -m=0`, zip SHA-256 `19566762…0d52`).
  - Sage Plus was also run on the RAW files directly.
  - Every input SHA-256 is in `summary-b12.json` under `provenance.inputs`.
- **FASTA:** UniProt 2026_03 human (reviewed) + yeast (UP000002311), both from the pilot, with receipts.
  - Human SHA-256 `dffdff38…9c42`.
  - Yeast SHA-256 `1aaa8186…342a`.
  - Concatenated file SHA-256 `64b4055b…70dd`.
- **Config:** identical JSON for both engines (the plan in `tmt-known-ratio-candidates.md` §3).
  - Tryptic, 2 missed cleavages, 7-50 residues.
  - Static modifications: TMTpro +304.207146 on K and the N-terminus; carbamidomethyl C.
  - Variable modification: Met oxidation.
  - Precursor tolerance ±10 ppm; isotope errors -1..3; generated decoys.
  - hrMS2: fragment tolerance ±20 ppm, deisotoping on, `tmt_settings.level` 2.
  - SPS-MS3: fragment tolerance ±0.5 Da (ion-trap CID), no deisotoping, `level` 3.
  - Reporters are raw intensities (`sn: false`). Neither engine does impurity correction or purity filtering.
  - All other parameters are each engine's own default.

## Metrics

- **Filtering:**
  - Rank-1 target PSMs with spectrum q ≤ 0.01 and peptide q ≤ 0.01.
  - Species comes from the accession suffix. Peptides are dropped if their I/L-normalized sequence also occurs in the other proteome.
- **Normalization:** per file, each channel is scaled to the human median (human is 1:1 by design).
- **Ratio:** for yeast, log2(channel / mean of the 1x channels). Each peptide takes the median over its PSMs, then peptide × channel values are pooled per level.
- **Bias, MAE and fold:**
  - Bias is the median of observed minus expected log2.
  - MAE is the mean |error|.
  - The fold is 2^-median, which is the paper's "1x : level" ratio.
- **Slope:** the OLS slope of observed on expected log2 over the four spiked levels. A slope of 1 means no compression.
- **Leakage:** yeast 126 / 1x mean, over PSMs where 126 was observed.
- **Missing-channel rate:** fraction of channels 127N-134N that were missing (0, null or NaN) among quantified PSMs.
- **Human CV:** CV across 16 normalized channels for human PSMs with all channels observed.

## Results (Sage Plus beta.12 vs upstream df92199)

| Run | Quant PSMs human / yeast | Yeast peptides | Bias 2:1 | 4:1 | 10:1 | 20:1 | MAE 20:1 | Fold 20:1 (paper) | Slope | 126 leak | Missing | Human CV |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| upstream hrMS2 mzML | 23,849 / 8,577 | 6,657 | 0.230 | 0.483 | 1.083 | 1.702 | 1.778 | 6.15 (7.69) | 0.542 | 0.129 | 0.47% | 8.87% |
| Plus hrMS2 mzML | 24,017 / 8,637 | 6,701 | 0.231 | 0.484 | 1.084 | 1.702 | 1.777 | 6.15 (7.69) | 0.543 | 0.130 | 0.48% | 8.88% |
| Plus hrMS2 RAW | 23,944 / 8,616 | 6,697 | 0.231 | 0.485 | 1.085 | 1.702 | 1.778 | 6.15 (7.69) | 0.542 | 0.130 | 0.48% | 8.89% |
| upstream SPS-MS3 mzML | 15,530 / 5,528 | 4,496 | 0.049 | 0.104 | 0.339 | 0.586 | 0.888 | 13.32 (13.55) | 0.768 | 0.041 | 3.27% | 7.80% |
| Plus SPS-MS3 mzML | 15,707 / 5,590 | 4,539 | 0.052 | 0.107 | 0.343 | 0.591 | 0.893 | 13.28 (13.55) | 0.767 | 0.041 | 3.25% | 7.83% |
| Plus SPS-MS3 RAW | 15,654 / 5,572 | 4,536 | 0.052 | 0.107 | 0.343 | 0.592 | 0.894 | 13.27 (13.55) | 0.767 | 0.042 | 3.24% | 7.83% |

Observed fold (1x : level) against the paper:

| Method | Engine | 2:1 | 4:1 | 10:1 | 20:1 |
|---|---|---|---|---|---|
| hrMS2 | upstream / Plus | 1.71 / 1.71 | 2.86 / 2.86 | 4.72 / 4.72 | 6.15 / 6.15 |
| hrMS2 | paper | 1.81 | 3.16 | 5.67 | 7.69 |
| SPS-MS3 | upstream / Plus | 1.93 / 1.93 | 3.72 / 3.71 | 7.90 / 7.89 | 13.32 / 13.28 |
| SPS-MS3 | paper | ~1.9 | n/a (RTS 3.64) | n/a (RTS 8.01) | ~13.55 |

**High-signal sensitivity analysis.** The paper keeps only peptides with S/N ≥ 10 in every channel and summed S/N ≥ 200. Neither engine outputs S/N from these inputs. As an approximation, this analysis keeps yeast PSMs with 127N-134N all observed and summed intensity at or above the median.

| Run | PSMs | 2:1 | 4:1 | 10:1 | 20:1 | Slope |
|---|---|---|---|---|---|---|
| upstream hrMS2 | 4,095 | 1.75 | 3.06 | 5.43 | 7.60 | 0.615 |
| Plus hrMS2 mzML | 4,121 | 1.75 | 3.05 | 5.44 | 7.61 | 0.616 |
| upstream SPS-MS3 | 2,464 | 1.97 | 3.89 | 8.71 | 16.00 | 0.851 |
| Plus SPS-MS3 mzML | 2,492 | 1.96 | 3.88 | 8.68 | 15.95 | 0.850 |

**Yield at 1% q:**

| Run | Accepted PSMs | Human proteins | Yeast proteins |
|---|---|---|---|
| upstream hrMS2 | 32,981 | 5,092 | 2,080 |
| Plus hrMS2 | 33,228 | 5,089 | 2,087 |
| upstream SPS-MS3 | 22,385 | 4,122 | 1,567 |
| Plus SPS-MS3 | 22,676 | 4,134 | 1,576 |

Protein counts are distinct protein strings on quantified PSMs with protein q ≤ 0.01.

- Plus quantifies 0.7-1.4% more PSMs in every file (see `paired_upstream_vs_plus_per_file`).
- About 4% of SPS-MS3 identifications have no linked MS3 reporter scan in either engine.

**RAW vs mzML parity (Sage Plus):**

| Method | Shared quantified spectra | Same peptide | Reporter values equal | Channels missing in one input only |
|---|---|---|---|---|
| hrMS2 | 33,114 | 33,114 (100%) | 100% | 0 |
| SPS-MS3 | 21,689 | 21,689 (100%) | 100% | 0 |

- On shared spectra, reporter values are bit-identical. This includes the SPS-MS3 path, where the MS3 is linked to its MS2 through the RAW "Master Scan Number".
- The RAW runs accept about 0.3% fewer PSMs than the mzML runs.

## Interpretation

1. **Upstream and Plus quantify identically.**
   - Per-level bias differs by ≤0.005 log2.
   - Slope, leakage, missing rate and CV agree to the third decimal place.
   - The small differences come from the slightly different PSM sets, not from reporter extraction.
2. **Ratio compression is the expected physics.** hrMS2 compresses strongly (slope 0.54, 20:1 read as 6.1:1). SPS-MS3 is much closer to the truth (slope 0.77, 20:1 read as 13.3:1).
3. **Agreement with the paper:**
   - SPS-MS3 matches the published medians without any filter.
   - hrMS2 matches once low-signal PSMs are removed: 7.60 vs 7.69 at 20:1 and 5.43 vs 5.67 at 10:1.
   - So the unfiltered gap at hrMS2 is the paper's S/N filter, not an extraction problem.
   - With the high-signal filter, SPS-MS3 overshoots at 20:1 (16.0 vs 13.55). The S/N and intensity filters are not equivalent.
4. **Interference in the zero channel (126).**
   - hrMS2 yeast PSMs carry 13% of the 1x signal in 126, which is 70% of the 1/20 level.
   - SPS-MS3 carries 4.1%, which is 46% of the 1/20 level.
   - The 1/20 level is therefore dominated by co-isolated human signal in hrMS2.
5. **Missing channels.**
   - Missing channels are rare in hrMS2: 0.5% of 127N-134N values, 1.1% for yeast.
   - They are more frequent in SPS-MS3: 3.3%.
   - 126 is missing for 17% (hrMS2) and 35% (SPS) of yeast PSMs.
   - The beta.12 zeros inflate nothing here because they are treated as missing. A consumer that averages zeros would bias the 1/20 and 126 channels downward. This is the Beta 13 fix.

## Caveats

- Three files per method gives descriptive medians, not inference.
- Intensities, not S/N: the paper's S/N filter is only approximated.
- No impurity correction was applied, and none is published for this lot in the deposit.
- FAIMS multi-CV files were searched as-is.
- Wall time and RSS are recorded per job in `summary-b12.json`. They come from a shared machine with other agents building and are not a performance claim.
- The SPS-MS3 4:1 and 10:1 paper values are for RTS-MS3 only. The research write-up listed them under SPS, which is corrected here.

## Files

- `tmt_hypro16.py`: the whole pipeline (fetch, convert, search, analyze).
- `run.sh`: the uv wrapper with pinned pyarrow and numpy.
- `test_tmt_hypro16.py`: checks that 0, null and NaN are all treated as missing.
- `summary-b12.json` and `table-b12.md`: full metrics, provenance and configs.
- `channel-map-source.md`: the source for the channel map.
- Data, run outputs and logs: `/mnt/data1/sage-plus-scientific/tmt-hypro16/` (4.2 GB).
