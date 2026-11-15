# Head-to-head: upstream Sage v0.15.0-beta.2 vs Sage Plus main (9483f1c)

Measured on 2026-10-11, back to back on one machine. Engines ran strictly one at a time.
All numbers below come from files in this directory. `results.json` holds the same numbers in machine-readable form, and `make_results.py` regenerates it.

## Summary

| # | Comparison | Upstream Sage v0.15.0-beta.2 | Sage Plus main 9483f1c | Note |
|---|---|---|---|---|
| 1 | DDA HYE (PXD028735, 4 files, 3 repeats) | 247,102 PSMs / 57,009 peptides / 6,916 protein groups; median wall 107.6 s; median peak RSS 7.57 GiB | 250,060 / 57,700 / 7,028; median wall 113.7 s; median peak RSS 6.19 GiB | IDs were identical in all 3 repeats. Wall time is within noise, and the Plus median is slower. |
| 2 | ETD (PXD018176, dual HCD/ETD) | `z` ions: 14,395 PSMs (ETD scans 4,515) / 7,495 peptides / 711 groups; 7.5 s; 1.27 GiB | `z_dot` ions: 16,818 PSMs (ETD scans 6,799) / 7,631 / 725; 6.31 s; 0.82 GiB. With `z` ions: 14,686 (ETD 4,674) / 7,588 / 719 | Medians of 3 repeats, searching the mzML. ETD-scan PSMs: +51% with `z_dot`, +3.5% with `z`. Re-measured on 2026-11-14 after the c/x/z fragment-mass fix, reading the `.raw`: 16,812 PSMs (ETD 6,803) / 7,637 / 723, +45.5% against our own `z` arm and the +51% against upstream unchanged (see section 2). |
| 3 | Phospho (PXD007058, 2 files, max_variable_mods 1) | 14,825 PSMs / 6,583 peptides / 4,845 phospho peptides; 26.6 s; 8.06 GiB. Site localization: not available | 14,838 / **6,226** / **4,484**; **44.3 s**; 4.75 GiB. Localization: 2,930 phospho protein sites at 1% site FDR | Plus finds fewer peptides and is slower. The mass-offset arm gives 15,227 / 6,405 / 4,664 in 29.0 s and 1.83 GiB. |
| 4 | Huge DB (human + 100x IGC gut catalog) | not measured | not measured | Blocked: needs 24 GB free under memgate; the shared machine had about 11 GiB available. |
| 5 | Thermo RAW end to end (2 PXD007058 files) | ThermoRawFileParser 89.6 s + search 29.0 s = 118.6 s (14,825 PSMs) | Direct .raw read: 44.5 s (14,834 PSMs) | Plus is 2.7x faster end to end. Plus IDs from .raw are within 4 PSMs of its own mzML run. |
| 6 | Entrapment FDP at 1% (Human_01 vs human+yeast+E. coli) | Combined FDP: PSM 1.64%, peptide 1.77%, protein 2.80%, protein group 2.91% | PSM 1.25%, peptide 1.51%, protein 1.71%, protein group 2.02% | Plus is closer to the 1% target at every level. Both still exceed 1% by the combined estimator. |

### Where Sage Plus is worse or within noise

1. **C3 phospho peptides:** 6,226 vs 6,583 (-5.4%), and phospho peptides 4,484 vs 4,845 (-7.5%). PSMs are tied (14,838 vs 14,825), and phospho PSMs are identical (10,044).
2. **C3 phospho speed:** 44.3 s vs 26.6 s. Plus is 1.7x slower here, in a single run.
3. **C1 wall time:** Plus median 113.7 s vs upstream 107.6 s. Both ranges overlap (Plus 103–126 s, upstream 105–119 s), so this is within noise under a loaded machine.
4. **Ties:** C1 proteins 6,893 vs 6,888, and C3 PSMs 14,838 vs 14,825.

## Setup

- **Upstream binary:** the official release tarball `sage-v0.15.0-beta.2-x86_64-unknown-linux-gnu`, unpacked on a desktop workstation, local drive, not retained. It was run with `--disable-telemetry-i-dont-want-to-improve-sage`.
- **Sage Plus binary:** built for this session on a desktop workstation, local drive, not retained.
  - Built with `cargo build --release` from a clean worktree at main `9483f1c4eb33f3f2ee4b4b5c1b4475fa6d207b34` ("Use the reader crates from crates.io"). The build log is `build.log` (5m42s).
  - The binary reports version 0.1.0-beta.17. The version string has not been bumped since Beta 17.
  - sha256 `745de67125784e1ad2d151e5598071fd255fac43f39739314805731fbc87529a`
- **Hardware:** Intel i7-10700K (8 cores / 16 threads), 31 GiB RAM, local data drive.
- **Threads:** all 16 threads by default (`RAYON_NUM_THREADS` unset), except C4, which was set to 8 to match Beta 11.
- **Machine load:** the desktop was in use throughout. User apps (browser, chat apps) held about 19–20 GiB, swap was full, and MemAvailable was 3–11 GiB.
  - `loadlog.tsv` samples the 1-minute load average and MemAvailable every 30 s. Over the whole session the mean was 13.2, the max 24.7, and the minimum MemAvailable 3.2 GiB.
  - The load average includes the engine's own 16 threads.
  - Per-run load windows are in `results.json` under `load`. Runs shorter than 30 s have no sample.
- **Memory gate:** every search ran under a `memgate.sh` wrapper on the same workstation (local drive, not retained) with `MEMGATE_HEADROOM_GB=2`. The default of 3 never admitted an 8 GB job on this machine.
- **Run directories:** each one (`<comparison>/<arm>/runN/`) holds:
  - `config.used.json`, `command.txt`, `env.txt` (threads, gate, binary sha256)
  - `stdout.log`, `stderr.log`, `time.txt` (`/usr/bin/time -v`)
  - `started_at.txt`, `finished_at.txt`, `after.txt` (exit code, uptime, free)
- **Driver scripts:** `run.sh` (one run), `chain.sh`, `chain2.sh`, `chain3.sh` (run order), `watchdog.sh`, and `chain.log` (timestamps).
- **ID definitions (`analyze_runs.py`, same as `benchmarks/HEADTOHEAD.md`):**
  - PSMs: target rows with `spectrum_q` ≤ 0.01.
  - Peptides: distinct `peptide` among targets with `peptide_q` ≤ 0.01.
  - Protein groups: distinct `protein_groups` among 1% PSMs with `protein_group_q` ≤ 0.01.
  - Upstream is read from `results.sage.tsv`, and Sage Plus from `results.sage.parquet`.
  - Upstream's own log line "N PSMs at 1% FDR" includes decoys, so it is higher than the target count reported here.
- **Sage Plus defaults that differ from upstream:** Met clip, I/L merge, picked protein FDR, nonlinear RT, and averagine deisotoping. Both engines were given the same config file, so these defaults are part of what is compared.
- **Configs:** all are in `configs/`.
  - Since Beta 16, Sage Plus requires `ion_kinds` and `semi_enzymatic` to be explicit. They were set to upstream's defaults (`b,y`, `false`) wherever a config lacked them.
  - Every config passed `sage --validate-only`.

## 1. DDA HYE: PXD028735

- **Files:** `LFQ_Orbitrap_DDA_Condition_{A,B}_Sample_{Alpha,Beta}_01.mzML` (4 files)
- **FASTA:** `hye-irt-defined.fasta`
- **Config:** `configs/c1-dda.json`, the same config as the Beta 16 head-to-head. LFQ is on.
- **Run order:** upstream run1, Plus run1, upstream run2, and so on, alternating. Memgate 8.

| | Upstream | Sage Plus |
|---|---|---|
| PSMs | 247,102 | 250,060 (+1.2%) |
| Peptides | 57,009 | 57,700 (+1.2%) |
| Proteins | 6,888 | 6,893 |
| Protein groups | 6,916 | 7,028 (+1.6%) |
| Wall time, runs 1/2/3 | 105.1 / 107.6 / 118.6 s | 103.2 / 113.7 / 126.2 s |
| Wall time, median | 107.6 s | 113.7 s |
| Peak RSS, median | 7.57 GiB | 6.19 GiB (-18%) |
| Mean load during runs 1/2/3 | 11.7 / 17.0 / 19.9 | 17.1 / 18.0 / 22.1 |

- The ID counts were identical across all three repeats of each engine. They also match the Beta 16 head-to-head exactly.

**LFQ.** These are precursors with `q_value` ≤ 0.01 that were quantified in all 4 files. The ratio is log2 A/B; the expected values are human 0, yeast +1, E. coli −2.

| Species | Upstream n / median / IQR / MAD | Sage Plus n / median / IQR / MAD |
|---|---|---|
| Human | 26,077 / +0.009 / 0.367 / 0.180 | 26,676 / +0.017 / 0.351 / 0.169 |
| Yeast | 7,427 / 1.117 / 0.448 / 0.224 | 7,585 / 1.114 / 0.413 / 0.205 |
| E. coli | 1,370 / −2.219 / 0.866 / 0.435 | 1,400 / −2.215 / 0.849 / 0.421 |

- Precursors at 1% / in all 4 files: upstream 35,127 / 35,119, Plus 35,990 / 35,981.
- Details are in `c1-dda/metrics.json` and `c1-dda/lfq_ratios.parquet`. `c1-dda/analyze.py` produces them.

## 2. ETD: PXD018176

- **File:** `B191108_09_Lumos_AB_DE_165_Ec_dualHCD_ETD_1.mzML`
- **FASTA:** `ecoli_k12_sp.fasta`
- **Search settings:** missed cleavages 2, length 7–50, CAM fixed, Met-ox variable, max_variable_mods 2, precursor ±10 ppm, fragment ±20 ppm, isotope errors −1..3, deisotope on, max_fragment_charge 2.
- **Ion kinds:** `b,y,c,z` (`configs/c2-etd-z.json`) or `b,y,c,z_dot` (`configs/c2-etd-zdot.json`).
- **Runs:** 3 repeats per arm, alternating. Memgate 4.
- **Scan types:** the HCD/ETD split came from a scan-group table for `B191108_09_Lumos_AB_DE_165_Ec_dualHCD_ETD_1` (PXD018176) built on a desktop workstation, local drive, not retained. It was rebuilt from the `.raw` for the 2026-11-14 re-measurement below, which agrees with it to within 15 HCD-scan PSMs.

| Arm | PSMs | HCD-scan PSMs | ETD-scan PSMs | Peptides | Protein groups | Median wall | Median RSS |
|---|---|---|---|---|---|---|---|
| Upstream, `z` | 14,395 | 9,880 | 4,515 | 7,495 | 711 | 7.50 s (7.50 / 8.41 / 6.67) | 1.27 GiB |
| Sage Plus, `z_dot` | 16,818 | 10,019 | 6,799 | 7,631 | 725 | 6.31 s (6.47 / 6.31 / 6.04) | 0.82 GiB |
| Sage Plus, `z` | 14,686 | 10,012 | 4,674 | 7,588 | 719 | 6.15 s (6.33 / 6.15 / 6.08) | 0.81 GiB |

- The IDs were identical across the 3 repeats of every arm.
- Upstream has no `z_dot` ion kind.
- Plus with `z` measures the engine difference at equal ion kinds. Plus with `z_dot` adds the effect of the ion kind.

### Re-measured 2026-11-14, after the c/x/z fragment-mass fix

The ammonia constant in `IonSeries::new` was one electron light until main `fa072b2`, so c ions
were 0.549 mDa too light and x, z and z-dot ions 0.549 mDa too heavy. The arms above were
measured before that fix. They were re-measured afterwards, and **this is not the same input**:
the mzML of the table above and the `ecoli_k12_sp.fasta` it used are both gone, so the
re-measurement reads the deposited `.raw` of the same file,
`B191108_09_Lumos_AB_DE_165_Ec_dualHCD_ETD_1.raw` from PXD018176 (604,032,589 bytes, sha256
begins `03f00f59`), against the submitters' own FASTA deposit `UP000000625_83333.fasta`
(4,391 entries, sha256 begins `45ede201`). Binaries: main `0e2eff8` before, `fa072b2` after,
each built in its own worktree and target directory. The HCD/ETD split was rebuilt from the
`.raw`'s own scan activations (6,765 MS1, 13,846 MS2 HCD, 13,522 MS2 ETD) and agrees with the
split used above to within 15 HCD-scan PSMs. These were not timing runs, so no wall time or
peak memory from them is quoted. Every number is from each arm's `run-summary.json`
(`psms/peptides/proteins/protein_groups_at_one_percent_fdr`). Full record:
`audits/etd-mass-fix-check-2026-11.md` in the project's shared audit directory, with its own
`README.md` and `checksums.txt` beside the data.

| Arm | | PSMs | ETD-scan PSMs | HCD-scan PSMs | Peptides | Protein groups |
|---|---|---:|---:|---:|---:|---:|
| `b,y,c,z` | before (`0e2eff8`) | 14,678 | 4,673 | 10,005 | 7,589 | 718 |
| | after (`fa072b2`) | 14,676 | 4,675 | 10,001 | 7,590 | 718 |
| `b,y,c,z_dot` | before (`0e2eff8`) | 16,796 | 6,792 | 10,004 | 7,629 | 724 |
| | after (`fa072b2`) | **16,812** | **6,803** | 10,009 | **7,637** | 723 |

- **The October baseline reproduces through the format change:** every column of the `before`
  rows is within 0.13% of the published table (PSMs 14,686 to 14,678 and 16,818 to 16,796),
  across a different file format, a different FASTA deposit and a month of commits.
- **The fix costs nothing and gains a little:** +16 PSMs and +8 peptides (+0.10%) on the
  `z_dot` arm, −2 PSMs on the even-electron `z` arm.
- **The evidence that it is right is the mass error, not the counts.** Median signed fragment
  mass error over matched rank-1 target fragments: b and y identical to four decimals before
  and after (−0.268 and −0.411 ppm), c moves +0.226 to −0.526 ppm and z-dot −1.236 to
  −0.456 ppm, each by 0.549 mDa at its series' median m/z and in opposite directions. The
  z-dot series, which disagreed with b and y by about 0.8-1.0 ppm, now agrees with them inside
  0.2 ppm, and the whole-file `recommended_tolerances.fragment.sigma_ppm` falls from 1.628 to
  1.528 ppm on that arm.

## 3. Phospho: PXD007058

- **Files:** `SF_200217_U2OS_TiO2_HCD_nlETcaD_quad_OT_rep1.mzML` and `SF_200217_U2OS_TiO2_HCD_nlEThcD_OT_rep2.mzML`
- **FASTA:** `human_sp.fasta`
- **Search settings:** missed cleavages 1, CAM fixed, Met-ox and STY +79.966331 variable, `b,y`, precursor ±10 ppm, fragment ±10 ppm.
- **Runs:** one run per arm. Memgate 9, plus a 12 GiB RSS watchdog.

**Changed setting: max_variable_mods is 1, not 2.**
- At 2, the digest is about 32M peptides and 1.5–1.7B fragments, which the `--estimate` output confirmed. Upstream would need about 20 GB, more than the machine had free.
- At 1, upstream peaked at 8.06 GiB.
- Both engines used the same value.

| Arm | Config | PSMs | Peptides | Protein groups | Phospho PSMs | Phospho peptides | Index (peptides / fragments) | Wall | Peak RSS |
|---|---|---|---|---|---|---|---|---|---|
| Upstream | `c3-phospho.json` | 14,825 | 6,583 | 2,467 | 10,044 | 4,845 | 11.78M / 434M | 26.55 s | 8.06 GiB |
| Sage Plus | `c3-phospho.json` | 14,838 | 6,226 | 2,491 | 10,044 | 4,484 | 11.89M / 439M | 44.26 s | 4.75 GiB |
| Sage Plus, mass offset | `c3-phospho-offset.json` | 15,227 | 6,405 | 2,506 | 10,433 | 4,664 | 3.85M / 121M | 29.00 s | 1.83 GiB |
| Sage Plus, localization | `c3-phospho-loc.json` | 14,838 | 6,226 | 2,491 | 10,044 | 4,484 | 11.89M / 439M | 46.68 s | 4.89 GiB |

**Mass-offset arm.** Phospho uses `search_mode: mass_offset` and is not enumerated as a variable mod. The `max_variable_mods: 1` limit applies only to the enumerated mods (Met-ox), so a phospho offset can sit on top of an oxidation. Its search space is therefore not identical to the other arms.

**Localization arm (`ptm_localization: {enabled: true}`, Sage Plus only):**
- 10,346 localized PSMs.
- 6,014 phospho PSM-sites at localization q ≤ 1%, which is 2,838 unique peptide+position pairs.
- 2,930 phospho protein sites at site q ≤ 1% (2,977 target and 23 decoy protein sites were scored).
- Upstream has no site localization, so its value is "not available".

**Peptide gap.** Upstream 6,583 vs Plus 6,226. The measured facts:
- After merging I and L the counts are unchanged, so the I/L merge does not explain it.
- 5,716 peptides are shared, 867 are upstream-only, and 510 are Plus-only.
- Distinct stripped sequences are upstream 5,449 vs Plus 5,474. The gap is therefore in modified forms, not in sequences.
- No root cause was measured.

## 4. Huge DB: blocked, not measured

- **Configs:**
  - `configs/c4-plus.json`: copied from the Beta 11 run `large-db-20260925/runs-v2-mp4/hek-igc-100x-prefilter/config.json`, with the input paths changed to this workstation's local drive (not retained) and `semi_enzymatic` / `ion_kinds` made explicit.
  - `configs/c4-upstream.json`: the same settings with the Plus-only keys removed (`max_memory_gb`, `min_free_memory_gb`, `output_filter`, `score_type`, `prefilter*`).
- **FASTA:** `human-igc-100x.fasta` (human + gut IGC catalog).
- **Input:** `PXD001468 b1906_293T_proteinID_01A_QE3_122212.mgf`.
- **Run settings:** `RAYON_NUM_THREADS=8`, and memgate 24 with budget 24 for both engines. Upstream also had a watchdog that kills it by PID above 24 GiB RSS. See `chain3.sh`.
- **Why it is blocked:** memgate needs 24 GiB + 2 GiB headroom of MemAvailable. From 09:34:38 to 09:50:28 (about 16 min) it reported only 10–11 GiB available, because user apps held about 20 GiB and swap was full. The job was never admitted. The waiting chain was then stopped by PID, as logged in `chain.log`, so upstream never started.
  - The wait log is `c4-hugedb/plus/run1/stderr.log`.
  - Neither engine produced any number.
- **FASTA discrepancy:** the task named `gut-db2mg-human.fasta`. The Beta 11 result usually quoted (14,459 PSMs / 451 s / 19.19 GiB) actually used `human-igc-100x.fasta`, `RAYON_NUM_THREADS=8` and a pre-release binary that reported beta.10, so C4 was set up to match that run.

## 5. Thermo RAW end to end: PXD007058, the same 2 files as C3

- **Converter:** ThermoRawFileParser 2.0.0-dev, the official Linux release.
  - Unpacked on a desktop workstation, local drive, not retained; the release zip's sha256 begins `19566762`.
  - Arguments: `-f 2 -l 3`, the same as the original conversions. The commands are in `c5-raw/convert/*.command.txt`.
- **Upstream:** convert both files, then search the fresh mzMLs (`configs/c5-upstream-mzml.json`).
- **Sage Plus:** search the `.raw` files directly (`configs/c5-plus-raw.json`).
- **Search settings:** the same as C3.

| Step | Upstream path | Sage Plus path |
|---|---|---|
| Convert rep1 | 45.91 s (0.76 GiB RSS) | — |
| Convert rep2 | 43.66 s (0.74 GiB RSS) | — |
| Search | 29.01 s (8.36 GiB RSS) | 44.48 s (5.10 GiB RSS) |
| **End to end** | **118.58 s** | **44.48 s** |
| PSMs / peptides / protein groups | 14,825 / 6,583 / 2,467 | 14,834 / 6,226 / 2,488 |
| Phospho PSMs / phospho peptides | 10,044 / 4,845 | 10,042 / 4,485 |

- Upstream's IDs on the fresh conversions are identical to C3.
- Plus reading `.raw` differs from Plus reading mzML (C3) by 4 PSMs, 0 peptides, 3 protein groups, 2 phospho PSMs and 1 phospho peptide.
- Single run per step.

## 6. Entrapment FDP

- **Sample:** `LFQ_Orbitrap_DDA_Human_01.mzML`, a human-only sample.
- **FASTA:** `hye-irt-defined.fasta` (human + iRT targets, with yeast + E. coli as entrapment).
- **Config:** `configs/c6-entrapment.json` (C1 settings, LFQ off).
- **Method:** the estimators and filters of `benchmarks/picked_fdr_entrapment.py`, applied by `entrapment_both.py`, which reads both output formats.
  - Lower bound = N_e/(N_t+N_e).
  - Combined = N_e(1+1/r)/(N_t+N_e), from Wen et al. 2025.
  - r = 0.387 (peptides) and 0.512 (proteins).
  - Protein and protein-group rows are single-protein / single-group only.
- **Runs:** one run each. Upstream 45.9 s / 4.45 GiB, Plus 44.0 s / 3.20 GiB.

At q ≤ 1%, each cell is targets / entrapments / lower bound / combined:

| Level | Upstream | Sage Plus |
|---|---|---|
| PSM | 64,580 / 296 / 0.46% / 1.64% | 65,275 / 229 / 0.35% / 1.25% |
| Peptide | 37,225 / 185 / 0.49% / 1.77% | 37,620 / 159 / 0.42% / 1.51% |
| Protein | 4,286 / 41 / 0.95% / 2.80% | 4,284 / 25 / 0.58% / 1.71% |
| Protein group | 4,318 / 43 / 0.99% / 2.91% | 4,360 / 30 / 0.68% / 2.02% |

Combined FDP at the other thresholds:

| Threshold | Upstream PSM / peptide / protein / group | Sage Plus PSM / peptide / protein / group |
|---|---|---|
| 0.5% | 1.04 / 1.03 / 2.11 / 2.50% | 0.71 / 0.80 / 0.98 / 1.24% |
| 5% | 5.75 / 6.08 / 8.06 / 7.75% | 5.26 / 5.71 / 6.89 / 7.62% |

All thresholds, including 2%, are in `c6-entrapment/entrapment-{upstream,plus}.json`.

## Deviations from the requested plan

1. **C3:** max_variable_mods was 1 instead of 2, for the memory reason given in section 3. Both engines used the same value.
2. **C4:** not run, because memory was insufficient on the shared machine.
3. **Repeats:** C3, C5 and C6 are single runs. Only C1 and C2 have 3 repeats.
4. **Memgate:** headroom was 2 GiB instead of 3, and C1/C6 used gate 8 instead of 9/10. The higher values were never admitted under the background memory load.
5. **Run order:** an earlier attempt at the chain was restarted after a scripting error that killed the chain shell; the orphaned upstream run was killed by PID and rerun. Only completed runs with exit 0 are reported.
