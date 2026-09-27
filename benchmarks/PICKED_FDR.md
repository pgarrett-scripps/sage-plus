# Picked peptide / protein / protein-group q-values (Beta 13)

Verdict: the q-values inherited from upstream Sage were **not** textbook picked FDR and were
wrong by construction (category c): after picking, both pair members were ranked, so a target
beaten by its own decoy could still pass, and the decoy count was replaced by a model-based sum.
On two human-only entrapment searches the error is at the noise level (one protein per dataset
that lost to its own decoy passed 1%). The code now implements textbook picked FDR:
winners only, `(decoys + 1) / targets`, losers q = 1.

## 1. What the code computed (beta.12, `crates/sage/src/fdr.rs`)

The q-value loop is identical to upstream Sage (`lazear/sage` master, `Competition::assign_q_value`).
Sage Plus only added a count-based fallback for when the KDE cannot be fitted.

For every competition key (peptide: modified sequence, decoys mapped to their target;
protein: the single accession of a unique peptide; protein group: the IDPicker group, decoys
mapped to their target's group), the target score `f` and decoy score `r` are the best PSM
`sage_discriminant_score` on each side.

1. **PEP model.** A KDE mixture is fitted on the winners only: score `max(f, r)`, label decoy
   when `r >= f`. `PEP(s) = pi_d f_d(s) / (pi_d f_d(s) + pi_t f_t(s))`, that is
   P(winner is a decoy | s), made monotone.
2. **Ranking.** Both members of every pair that were observed are put in one list (winners
   *and* losers), sorted by their own score.
3. **q.** Walking down the list, `q_i = (1 + sum_{j<=i} PEP(s_j)) / #targets_{<=i}`, where the
   sum runs over targets and decoys alike; then the reverse cumulative minimum (capped at 1).
   There was no tie grouping in this path.

Spectrum q (`ml/qvalue.rs`) is plain TDC: rank all PSMs, `(decoys + 1) / targets` after each
complete tie group, reverse cumulative minimum. It is correct and was not changed. The per-PSM
`posterior_error` column is the LDA-score KDE and is separate from this code.

## 2. Comparison with the references

- **Elias & Gygi 2007 TDC, +1 (Levitsky 2017; He et al.; Käll 2008).** For a competition
  where each false target is as likely to win as its decoy, `(D + 1) / T` over the ranked
  winners is a conservative FDR estimate with a finite-sample guarantee.
- **Picked protein FDR (Savitski et al. 2015)** and **picked group FDR (The et al. 2022)**:
  a target and its decoy compete; the loser is *discarded*, and TDC runs over the winners.
  Picking is what removes the inflation of decoy counts by decoys of present proteins
  (classic protein TDC is conservative for large datasets) without double-counting.
- **Picked peptide (Lin et al. 2022)**: same double competition at peptide level.

Where Sage departed from this:

1. **Losers were ranked.** A target that lost to its own decoy stayed in the list with its own
   score and got its own q. In picked FDR it is removed (its best evidence is weaker than a
   random-match decoy). Its false-ness is only counted through `PEP(s)`, which is the decoy
   probability of a *winner* at that score, far below the real error probability of a
   target that lost to its decoy. This is anti-conservative. Losing decoys also add PEP mass
   (conservative), and losing targets add to the denominator (anti-conservative), so the
   net direction depends on the data.
2. **Sum of PEP instead of decoy count.** Over winners, `E[sum PEP] = #decoy winners`, so
   `(1 + sum PEP) / T` is a smoothed `(D + 1) / T`, valid in expectation *if* the KDE is
   calibrated. It loses the TDC finite-sample guarantee, depends on the bandwidth in the high
   score tail where decoys are sparse, and the `+ 1` no longer has its TDC meaning. Note
   `PEP` here is P(decoy | s), not the target PEP; summing target PEPs over targets would be
   the Käll "expected FDR", a different estimator.
3. **No tie grouping**, so equal scores could get different q-values.

The PEP-sum is a legitimate idea (it is approximately count-based), so on its own it would be
category (b). Ranking the losers makes the method wrong by construction (c), even if the
effect is small in practice.

## 3. New implementation

`Competition::assign_q_value` now keeps one row per pair (the winner; a tie goes to the
decoy), runs the existing count-based `assign_count_q_values` (ties grouped, `(D + 1) / T`,
reverse cummin, cap 1), and returns q = 1 for every losing member. The KDE is no longer used
for peptide/protein/group q. Keys, pairing and unique-peptide rules are unchanged.
Unit tests in `crates/sage/tests/unit/fdr.rs` cover hand-computed examples: +1 with a sparse
decoy, a decoy that beats its target (target discarded even though its score would pass), a
losing decoy that is not ranked, ties across a target/decoy boundary, a peptide-level
reversed-decoy pair, and protein q with a shared peptide.

No schema version bump: the columns, types and meaning (a q-value per level) are the same;
only the estimator changed, as with the Beta 8 group-pairing fix.

## 4. Entrapment measurement

Design: human-only samples searched against `hye-irt-defined.fasta` (20,416 human + 11 iRT
targets; 6,067 yeast + 4,396 E. coli entrapment proteins). Any yeast / E. coli
identification is false. Peptides shared with a human protein count as human.
FDP estimates: lower bound `N_e / (N_t + N_e)` and the combined estimator
`N_e (1 + 1/r) / (N_t + N_e)` (Wen et al. 2025) with `r` = 0.387 (distinct tryptic
peptides, I/L merged) for PSM and peptide level and `r` = 0.512 (protein count) for protein
and group level.

- `human01`: PXD028735 `LFQ_Orbitrap_DDA_Human_01.mzML` (Orbitrap DDA, 1 file).
- `hek293t`: PXD001468 `b1906_293T_proteinID_01A` + `b1937_..._01B` (Q Exactive, 2 files).

Search: trypsin, 1 missed cleavage, 7-50 residues, CAM fixed, Met-ox variable, ±10 ppm
precursor, ±20 ppm fragment, generated decoys, RT prediction on, `output_filter.psm_q_value`
1.0. Baseline = origin/main (beta.12, `54668e1`), picked = this branch. Only the FDR code
differs, so PSMs and scores are identical. Configs, logs and outputs:
`/mnt/data1/sage-plus-scientific/picked-fdr-20260927/`; analysis script
`benchmarks/picked_fdr_entrapment.py`.

### At q ≤ 1%

| data | level | method | targets | entrapment | lower bound | combined |
|---|---|---|---:|---:|---:|---:|
| human01 | PSM (reference) | both | 65,426 | 237 | 0.36% | 1.29% |
| human01 | peptide | baseline | 37,657 | 148 | 0.39% | 1.40% |
| human01 | peptide | picked | 37,749 | 155 | 0.41% | 1.47% |
| human01 | protein | baseline | 4,249 | 26 | 0.61% | 1.80% |
| human01 | protein | picked | 4,260 | 26 | 0.61% | 1.79% |
| human01 | protein group | baseline | 4,330 | 31 | 0.71% | 2.10% |
| human01 | protein group | picked | 4,339 | 31 | 0.71% | 2.09% |
| hek293t | PSM (reference) | both | 42,447 | 306 | 0.72% | 2.57% |
| hek293t | peptide | baseline | 23,294 | 191 | 0.81% | 2.91% |
| hek293t | peptide | picked | 23,275 | 188 | 0.80% | 2.87% |
| hek293t | protein | baseline | 6,130 | 73 | 1.18% | 3.47% |
| hek293t | protein | picked | 6,130 | 72 | 1.16% | 3.43% |
| hek293t | protein group | baseline | 6,254 | 87 | 1.37% | 4.05% |
| hek293t | protein group | picked | 6,251 | 86 | 1.36% | 4.01% |

### Other thresholds (lower bound / combined)

| data | level | q | baseline | picked |
|---|---|---:|---:|---:|
| human01 | peptide | 0.5% | 0.22% / 0.78% | 0.23% / 0.84% |
| human01 | peptide | 5% | 1.58% / 5.67% | 1.58% / 5.66% |
| human01 | protein | 0.5% | 0.31% / 0.92% | 0.31% / 0.92% |
| human01 | protein | 5% | 2.27% / 6.69% | 2.36% / 6.98% |
| hek293t | peptide | 0.5% | 0.68% / 2.43% | 0.68% / 2.43% |
| hek293t | peptide | 5% | 1.95% / 6.97% | 1.93% / 6.92% |
| hek293t | protein | 0.5% | 0.98% / 2.88% | 0.93% / 2.74% |
| hek293t | protein | 5% | 2.64% / 7.81% | 2.81% / 8.31% |

### Differences at 1%

- Targets at protein q ≤ 1% that lost to their own decoy: baseline 1 (human01) and 1
  (hek293t, an E. coli entrapment protein); picked 0 by construction.
- Peptides only in baseline / only in picked: human01 2 / 101 (7 entrapment);
  hek293t 22 (3 entrapment) / 0. Proteins: human01 1 / 12 (0 entrapment); hek293t
  1 (entrapment) / 0.
- Identification counts at 1%: peptides +0.24% (human01), -0.08% (hek293t); proteins
  +0.26%, 0; protein groups +0.21%, -0.05%.

### Reading the numbers

- Both methods give the same calibration within a few entrapment hits at every level and
  threshold. The old method's defect is real but rare on these data.
- The entrapment FDP already exceeds nominal at the PSM level, where q is plain TDC, and
  peptide / protein levels track it. The inflation is therefore not from the picked-FDR step.
  A large part is homology: of the entrapment peptides at 1% peptide q (baseline), 41 of 147
  (human01) and 26 of 185 (hek293t) are one substitution away from a human tryptic peptide,
  i.e. likely a real human spectrum explained by a conserved yeast / E. coli paralog. The
  combined estimator treats those as random false hits, so it is an upper-side estimate
  here, and the lower bound is the more reliable figure. By the lower bound, peptide and
  protein FDP at 1% are 0.4-0.8% and 0.6-1.2%.
- Protein-level FDP above the PSM level is expected for any protein FDR computed from the
  best peptide: one false peptide makes a false entrapment protein, while human proteins
  collect many true peptides.

## 5. Not changed

- FASTA-supplied decoys (`generate_decoys: false`): decoy accessions do not pair with their
  targets for `protein_q`, so every observed protein is its own winner (classic protein TDC,
  conservative). Documented in DOCS; FASTA decoys are not the endorsed path.
- LFQ precursor q (`picked_precursor`) is a count-based TDC over all target and decoy peaks,
  not a picked competition; outside the scope of this change.
