# Crosslink search: exploration and design recommendation

Status: milestone 1 (cleavable DSSO/DSBU) is implemented on the `feat/crosslink` branch.
The search lives in `crates/sage-xlink`. It is built into `sage` only with the `crosslink`
Cargo feature and runs only when the config has a `"crosslink"` block. See
[Milestone 1](#milestone-1-implementation) at the end of this document. Nothing here ships
in Beta 9.

## Recommendation in one paragraph

Start with **MS-cleavable linkers (DSSO, DSBU) searched through signature doublets**,
with monolinks as mass offsets alongside. A doublet gives the mass of one chain directly
and the precursor gives its partner, so every pair hypothesis becomes two *closed-window*
lookups in the existing fragment index. That avoids the n² pair space entirely, reuses
the mass-offset retrieval path (`page_search_shifted`) for link-site fragments, reuses
neutral-loss variants for stub fragments, and gives a testable first milestone on the two
synthetic-library benchmarks with known ground truth (PXD014337, PXD029252).
Non-cleavable linkers (DSS/BS3) need the "alpha then beta" open search. Do that
second.

## 1. How established tools do it

| Tool | Linkers | Pair candidate generation | Notes |
| --- | --- | --- | --- |
| xiSEARCH (Rappsilber lab) | any, incl. cleavable | Alpha by linear fragment match with open precursor window, then beta mass = precursor − linker − alpha, looked up closed | Scores linear and crosslinked together; FDR in the separate xiFDR (boosting, self vs between) |
| pLink 2 / 3 | any | Fragment-index open search of one chain, partner by mass | Two-stage coarse and fine scoring; separate intra/inter-protein FDR |
| MeroX | cleavable (DSBU, DSSO, CDI) | Signature doublet ("Rise" modes) prefilters chain masses | Also non-cleavable mode by brute force on small databases |
| XlinkX (Proteome Discoverer) | cleavable | Doublet detection in MS2, chain masses searched as linear peptides (MS2 or MS3) | Target-decoy on CSMs |
| MS Annika | cleavable (DSSO, DSBU, PhoX-style via cleavable mode) | Doublet detection, then MS Amanda scoring of each chain as a linear peptide with a stub modification | Residue-pair and CSM FDR |
| Kojak | any | Scores every peptide with an open "partner mass" modification at the link site, keeps top alpha list, then pairs within precursor mass | Classic alpha-then-beta, with Percolator downstream |
| Scout | cleavable (primary) | Doublet (and stub) evidence, indexed chain lookup, ML rescoring | Proteome-scale; separate FDR levels up to PPI |
| OpenPepXL (OpenMS) | non-cleavable + labeled | Enumerates pairs under precursor mass with prefiltering by linear ion matches | Slower, exhaustive-ish |

Common themes:

- **Open-modification view.** A crosslinked chain is a linear peptide carrying one
  modification of mass (partner + linker) at the link residue. Fragments N- or C-terminal
  to the link site are ordinary ions; fragments containing it are shifted by that mass.
  Linear ions alone rank alpha candidates; the shifted ions confirm them.
- **Alpha then beta.** Take the top-k alpha candidates from an open precursor window,
  then look up beta in a closed window at `precursor − linker − alpha`. Cost is
  O(spectra × k), not O(n²).
- **Cleavable linkers.** DSSO cleaves in HCD to leave an alkene (+54.0106) or an
  unsaturated thiol (+85.9826; sulfenic acid +103.9932 is the thiol plus water) on each
  chain. DSBU leaves Bu (+85.0528) or BuUr (+111.0320). Each released chain gives a
  same-charge doublet spaced 31.9721 Da (DSSO) or 25.9792 Da (DSBU), so a four-peak
  signature fixes both chain masses. Link-site fragments also appear with stubs, not the
  full partner. PhoX (+209.9637) is not MS-cleavable, but it is enrichable, so search it as
  non-cleavable.
- **Dead-ends and loops.** Monolinks (one end hydrolyzed: DSSO +176.0143, DSBU +214.0954,
  DSS/BS3 +156.0786, or +155.0946 amidated) are ordinary variable modifications. Loop-links
  (both ends on one peptide) add the linker mass once. Backbone cleavages *between* the two
  sites do not separate the peptide, so those ions are missing.
- **Two-chain annotation.** Each fragment is labeled with its chain (alpha/beta), ion
  type, ordinal, whether it contains the link site, and for cleavable linkers which stub it
  carries. A peak explained by both chains is counted once.

## 2. Crosslink FDR

- **Classes.** Each chain is a target or a decoy, so a CSM is TT, TD (or DT), or DD.
  A TD hit can arise from a false target with a random decoy, and a DD from two random
  chains. So the estimate is `FDR = (TD − DD) / TT` (Walzthoeni et al. 2012). Clamp at
  ≥ 0 and monotonize like `assign_count_q_values` in `crates/sage/src/fdr.rs`. With equal
  target and decoy databases, random TD ≈ 2 × DD; a deviation is a useful diagnostic.
- **Intra vs inter (self vs between).** Random pairs are overwhelmingly between proteins
  (n² vs n). Pooling lets easy intra-protein links hide a high between-protein error rate.
  Estimate FDR separately for the two groups (Fischer and Rappsilber 2017; Lenz et al. 2021).
  A pair is intra when the two chains share any protein accession, and inter otherwise.
  Overlapping identical sequences are evidence of homomultimers and are reported as a flag.
- **Levels.** CSM → peptide pair (unordered peptidoforms) → residue pair (protein and
  position of each link site) → PPI (protein pair). Error grows on aggregation because false
  hits spread over many unique residue pairs while true ones collapse. Each level needs its
  own TT/TD/DD count on its own best-scoring representatives, as `picked_peptide` and
  `picked_protein` already do for linear peptides. Report q at every level. Recommend
  users filter at the level they publish (usually residue pair or PPI).
- **Decoys.** Reuse Sage's reversed-internal decoys per chain (`Peptide::reverse`). They keep
  mass and the terminal residue, and they move internal link sites to mirrored positions.
  Decoy chains are retrieved by the same closed lookups, so both classes get equal
  opportunity.

## 3. Mapping onto Sage Plus

### What is reusable as-is

| Need | Existing code | How |
| --- | --- | --- |
| Monolinks | `SearchMode::MassOffset` (`crates/sage/src/modification.rs:54`), `MassOffset` (`database.rs:1652`) | Configure e.g. DSSO hydrolyzed +176.0143 on `internal_residue:K` and `protein_first:*` with `search_mode: "mass_offset"`. **Works today, configuration only.** For DSSO/DSBU add stub neutral losses (176.0143 − 54.0106, 176.0143 − 85.9826), `neutral_loss_mode: optional` |
| Linkable residues | `ModificationSpecificity::Internal` (`internal_residue:K`), protein N-term, motif sites (`motif.rs`) | Trypsin does not cut at a linked K, so the link site is internal except at protein N-term |
| Link-site shifted fragment retrieval | `IndexedQuery::page_search_shifted` (`database.rs:2257`), `offset_query` (`scoring.rs:469`) | A chain with known partner mass *is* a mass-offset hypothesis with `mass = partner + linker` and `shift = partner + linker` (non-cleavable) or `shift ∈ {light_stub, heavy_stub}` (cleavable). This is the same pair of lookups the offset search already does |
| Placing the link and scoring one chain | `Peptide::with_mass_offset` (`peptide.rs:723`), `mass_offset_sites` (`database.rs:2079`), `score_peptide` (`scoring.rs:1187`) | Build the chain with a synthetic `ModificationDefinition` of mass partner + linker. Its `neutral_losses` are `partner + linker − stub` for each stub, so `IonGroupSeries` generates intact and stub variants per cleavage. `score_peptide` already picks one variant per cleavage and charge |
| Per-site competition | `score_hypothesis` loop over sites (`scoring.rs:1140`) | Each link-site placement competes as its own candidate |
| Candidate prefilter | `spectrum_index.rs` (per-spectrum precursor windows plus shifted probes) | Chain-mass windows from doublets are just more `PrecursorWindow`s with stub shifts |
| q-value and competition code | `fdr.rs` `assign_count_q_values`, `Competition` | Pattern for TT/TD/DD counting at each level |

### What is new

1. **Dynamic per-spectrum offsets.** `PreScore.offset: u8` indexes the static
   `db.mass_offsets`. Crosslinks need a per-spectrum hypothesis table (chain mass, partner
   mass, shifts). Recommended: generalize `offset_query` and `fragment_shift` to take a
   small `OffsetHypothesis { precursor_delta, shifts: SmallVec<[f32; 2]> }` and let
   `matched_peaks_with_isotope` take a slice of them. Static mass offsets become the special
   case, so there is one code path and no fork of the scorer.
2. **Doublet detection** (prototyped): `crosslink::find_doublets` and
   `crosslink::pair_hypotheses`.
3. **Pair scoring.** Score alpha and beta separately with `score_peptide`, then combine,
   with a shared peak mask so one peak is not counted for both chains. Candidate combined
   score: hyperscore over the union of matches, plus `min(alpha, beta)` as the weak-chain
   feature. The weak chain is what separates TD from TT. Keep per-chain matched counts,
   doublet evidence (0, 1, or 2 chains observed), and stub-ion counts as rescoring features.
4. **Pair FDR** (`fdr.rs` or a new `crosslink_fdr.rs`): TT/TD/DD at CSM, peptide pair,
   residue pair, and PPI, each split into intra and inter. Start count-based. The LDA in
   `ml/linear_discriminant.rs` is hard-wired to `Feature`; a crosslink feature vector is a
   later step.
5. **Output.** Do not overload `results.sage.parquet`. Add `crosslinks.sage.parquet` (new
   schema file under `schemas/`, writer next to `crates/sage-cloudpath/src/parquet/results.rs`).
   Columns: `psm_id, filename, scannr, charge, expmass, calcmass, linker, link_type
   (crosslink|looplink|monolink), alpha_peptide, alpha_proteins, alpha_link_site,
   alpha_protein_site, beta_*` (same set), `decoy_class (TT|TD|DD), protein_relation
   (intra|inter|homomeric), alpha_score, beta_score, combined_score, doublets_observed,
   csm_q, peptide_pair_q, residue_pair_q, ppi_q`. Fragment annotations go into a
   `matched_fragments`-style table with extra `chain` and `stub` columns. A later step adds a
   CSV export in the xiVIEW input format for visualization. This bumps the run-summary
   schema version. Cascade pins it (see memory note on Cascade pins), so coordinate.
6. **Linear competition.** Keep the linear search running on the same spectra (xiSEARCH
   does this). A spectrum explained as well by a linear peptide should not be reported as a
   crosslink. Compare the best linear and best crosslink score per spectrum and keep the
   delta as a feature.

### Where it plugs into the flow

`Scorer::score` → new `Scorer::score_crosslinked(query)` when `crosslink` is configured:

1. `find_doublets` on the recalibrated spectrum (`Scorer::recalibrated`), for each precursor
   charge and isotope error.
2. `pair_hypotheses` gives up to N chain-mass pairs (cap at about 10 per spectrum).
3. For each chain mass: `db.query(chain_mass, precursor_tol, fragment_tol)`, then count
   `page_search` plus `page_search_shifted(peak, stub)` hits and keep the top-k chains
   (reuse `matched_peaks_with_isotope` via the generalized offsets).
4. Cross the top-k alpha with top-k beta (k about 5, so 25 pairs), place link sites with
   `mass_offset_sites`, score with `score_peptide`, and combine.
5. Postprocess in `crates/sage-cli/src/runner/postprocess.rs` with crosslink FDR, then the
   new writer.

## 4. Staged plan

| Stage | Scope | Effort |
| --- | --- | --- |
| M0 | Monolinks as mass offsets. Config recipe plus a documented example; validate monolink IDs on PXD014337 (the library contains monolinked peptides) | 2–3 days |
| **M1 (first real milestone)** | DSSO and DSBU crosslinks through doublets (steps 1–5 above). TT/TD/DD FDR at CSM and residue pair, intra/inter split. `crosslinks.sage.parquet`. Validate true FDR on PXD014337 and PXD029252 | 4–5 weeks: engine 1.5 wk, FDR 0.5 wk, output and CLI 1 wk, validation and tuning 1–1.5 wk |
| M2 | Non-cleavable (DSS/BS3, PhoX) with alpha-then-beta: open alpha window `[min_chain, M − linker − min_chain]` on linear ions only, top-k, closed beta | 3–4 weeks, most of it performance |
| M3 | Loop-links with correct fragments (suppress ions whose cleavage falls between the two sites, an `IonGroupSeries` filter). PPI-level FDR, LDA or ML rescoring on crosslink features, xiVIEW export | 2–3 weeks |
| M4 | Proteome-scale (whole lysate) validation and tuning | open-ended |

About 8–10 weeks of focused work from M0 through M1, M2 and M3 together. M1 alone is useful
to the DSSO and DSBU users who make up most current XL-MS studies.

### Validation data (verified on PRIDE)

- **PXD014337**: Beveridge et al. 2020, *Nat Commun* 11:742. Synthetic crosslinked peptide
  library built from 95 Cas9 tryptic peptides in groups. Crosslinks are true only within a
  group, so every identification can be classified true or false. DSS, DSBU and DSSO;
  Orbitrap Fusion Lumos and Q Exactive; HCD, MS3, and ETD methods. Contains monolinks.
  The paper reports actual false crosslink identification rates of 2.4–32% depending on
  the analysis strategy, which is the baseline to beat.
- **PXD029252**: Matzinger et al. 2022, *Nat Commun* 13:3975. Synthetic peptides from
  E. coli ribosomal proteins that mimic a protein complex, with intra- and inter-protein
  links. DSSO, DSBU, CDI, ADH, DHSO, azide-A-DSBSO; Q Exactive HF. It comes with IMP-X-FDR,
  which computes the experimentally validated FDR and compares results across search engines.

Metrics: true FDR at the residue-pair level against nominal 1% and 5%, split into intra
and inter. Also report true identifications at 1%. Compare with the published numbers for
MS Annika and MeroX on the same raw files. Run time and peak RSS go next to the existing
benchmark tables. For M4, a proteome-scale lysate set (for example Lenz et al. 2021, E. coli)
comes later. Its accession is not verified here.

## 5. Risks

- **n² search space.** Human tryptic has about 5.6 M indexed peptides (see
  `benchmarks/MASS_OFFSET.md`). Naive pairs number about 10¹³. M1 avoids this entirely:
  work per spectrum is (doublet hypotheses) × 2 closed lookups. M2 is O(spectra × open
  alpha query), comparable to Sage open search. The Beta 7 open search (−150 to +500 Da)
  already costs about 4× a closed search (`benchmarks/PREFILTER.md`), and the alpha window
  is wider (up to several kDa). Mitigations: restrict the index to linkable peptides
  (has `internal_residue:K` or protein N-term), a large cut. Rank alpha on linear ions,
  with `min_ion_index` pruning. Keep top-k small.
- **Memory.** No pair index is ever built. Memory stays at the linear index plus
  per-spectrum hypothesis vectors. The mass-offset path showed that index size stays fixed
  as search-time offsets are added.
- **Spurious doublets.** Random same-charge pairs at the spacing occur in dense spectra.
  The DSSO spacing (31.9721, one S) is 0.018 Da from O₂ (31.9898, for example two oxidations),
  which is about 12 ppm at 1500 Da. Mitigations: require deisotoped (known-charge) peaks
  or top-N intensity, use tight fragment tolerance, cap hypotheses per spectrum, and use the
  four-peak signature as a score feature, not a hard requirement.
- **Prefilter.** `spectrum_index.rs` works when windows are closed. Doublet chain masses
  are closed windows, so the M1 prefilter is an extension that adds windows per spectrum.
  For M2 the alpha window is open, so the prefilter keeps almost everything. Use the
  linkable-peptide restriction instead.
- **Precursor isotope errors.** Crosslinked precursors are large and often picked off the
  monoisotope. `isotope_errors` multiplies hypotheses, and the offset benchmarks show it is
  the largest cost lever.
- **Scoring calibration.** Hyperscore of a chain carrying a big "modification" is not on the
  same scale as linear hyperscores. The weak-chain score, not the sum, has to drive FDR,
  or TD/DD will be underestimated. Validate on ground truth before trusting nominal q.
- **Schema and downstream.** A new output file and a run-summary schema bump affect Cascade
  pins.

## Prototype

The prototype has been superseded by milestone 1. The doublet code now lives in
`crates/sage-xlink/src/doublets.rs`. The prototype had these parts:

- `CleavableLinker` with `DSSO` and `DSBU` constants. A test checks the stub masses sum to
  the linker mass.
- `find_doublets(spectrum, linker, max_charge, fragment_tol)` finds same-charge doublets.
  It handles deisotoped peaks and unknown-charge peaks the same way `FragmentMatchIndex`
  does.
- `pair_hypotheses(doublets, precursor_mass, linker, precursor_tol, min_chain_mass)` turns
  doublets into (alpha, beta) chain masses. Pairs where both chains were observed rank
  first, and symmetric duplicates are merged.

- `rank_chain_candidates(db, spectrum, chain_mass, shifts, ...)` runs a closed `db.query`
  at one chain mass. It counts the preliminary fragment matches, both unshifted and shifted
  by each stub (and by partner plus linker), keeps linkable (K-containing) peptides, and
  returns the top k.

The harness `crosslink_doublets` example has been removed. It wrote one JSON line per MS2
spectrum with:
- the doublets;
- doublet counts at four decoy spacings (real spacing −1.9, −0.7, +0.55, +1.3 Da);
- up to 10 pair hypotheses, each with its top-5 alpha and beta candidates.

## Results on PXD014337 (DSSO, Q Exactive HF-X, stepped HCD)

Setup:
- File: `XLpeplib_Beveridge_QEx-HFX_DSSO_stHCD.raw`, 36,568 MS2 spectra.
- Database: `cas9_crapome.fasta` from the paper's Supplementary Data 2 (33,058 indexed
  peptides with decoys).
- Tolerances: ±10 ppm precursor, ±20 ppm fragment, isotopes 0–2.
- Truth: 481 target CSMs from the paper's XlinkX 1% FDR list (`xlinkx_QExHFX_DSSO_1perFDR_crapome_CSMs.csv`),
  joined by scan number. Each CSM is labelled correct or incorrect from the peptide groups
  in SI Table 1: 423 correct, 58 incorrect.
- Harness cost: 12.5 s wall time, 0.74 GB peak RSS.

Scripts and outputs are in `explore-data/crosslink/runs/` on the data drive:
`analysis/eval_doublets.py` and `analysis/eval_m0.py`.

**Doublets.**
- 97.6% of MS2 spectra carry at least one DSSO doublet, averaging 14.7 per spectrum.
- 18.9% of spectra have a pair hypothesis in which both chains were seen as doublets.
- This is a crosslink library, so almost every spectrum is a crosslink or monolink, and
  both produce doublets. The prevalence is therefore not a specificity measure.

**Chain lookup, correct XlinkX CSMs (n = 423):**

| Metric | Count | Share |
| --- | ---: | ---: |
| Spectrum has any doublet | 423 | 100% |
| Doublet at the alpha chain mass | 416 | 98.3% |
| Doublet at the beta chain mass | 421 | 99.5% |
| True (alpha, beta) pair among the ≤10 hypotheses | 377 | 89.1% |
| True pair ranked first among the hypotheses | 368 | 87.0% |
| True alpha in the top 5, given the true pair | 377 | 100% |
| True alpha ranked first | 376 | 99.7% |
| True beta in the top 5, given the true pair | 377 | 100% |
| True beta ranked first | 368 | 97.6% |

- End to end, both true chains are in the top 5 for 377 of 423 spectra (89.1%).
- Most of the 46 misses trace to the precursor, not the doublets: 38 of the 481
  truth CSMs do not match our precursor mass within isotopes 0–2 (mispicked monoisotope or
  a larger isotope error).
- The 58 incorrect (cross-group) XlinkX CSMs also get doublets: 44 have both chains in the
  top 5. Doublets alone cannot catch those errors. They need chain scoring and FDR.

**False doublets.**
- *Chance doublets.* The three clean decoy spacings give 0.65–0.92 doublets per spectrum
  against 14.7 at the true spacing, so about 5% of doublets are chance. Per spectrum, 18–36%
  of spectra have at least one chance doublet.
- *Rejected decoy spacing.* The −1.9 Da decoy spacing (30.07 Da) sits on a real mass
  difference, near CH₂O (30.011 Da): 2.8 doublets per spectrum. It is excluded.
- *Doublets away from the chain masses.* Inside correct CSMs, only 13.9% of doublets
  (1,246 of 8,982; 21 per spectrum) sit at a true chain mass. The rest are mostly
  stub-carrying b/y fragments, which also show the 31.97 Da spacing, plus chance pairs.
- *Consequence.* Pair hypotheses must be capped and ranked. Ranking both-observed pairs
  first, then by intensity, puts the true pair first 87% of the time. Only 9 true pairs
  rank 2nd–7th.

**M0 monolinks** (normal `sage` binary; configs in `runs/m0-*.json`). Each form is a K
`mass_offset` with its stub neutral losses:

| Config | Target PSMs at 1% | Monolink PSMs | Monolink peptidoforms |
| --- | ---: | ---: | ---: |
| Linear only (no K offsets) | 130 | – | – |
| Hydrolyzed (+176.0143) | 791 | 666 | 76 (60 sequences) |
| Hydrolyzed + Tris (+279.0777) + amidated (+175.0303) | 969 | 845 | 145 (74 sequences) |

- Split of the 845 monolink PSMs in the three-form run: 536 hydrolyzed, 10 Tris, 299
  amidated.
- *Amidated is inflated.* The amidated and hydrolyzed forms differ by 0.984 Da. That is
  within 10 ppm of one isotope step (1.003 Da), so with `isotope_errors` [0, 2], 135 of
  the 299 amidated PSMs are hydrolyzed monolinks assigned as amidated plus an isotope
  error. The M0 recipe should ship hydrolyzed plus Tris only, or amidated with
  `isotope_errors` [0, 0].
  - *Fixed in milestone 1.* On exact hyperscore ties the scorer now prefers the smaller
    isotope error. That cuts isotope-shifted amidated PSMs from 177 to 74, and the linear
    search is unchanged at 130 PSMs.
- *Peptide-level q.* It came out as 0 peptides in the one- and zero-offset runs but as 165
  in the three-form run. The picked-peptide step looks unstable on a database this small.
  Treat the PSM level as the M0 metric here.
- *Conclusion.* The mass-offset path gives about 6× the linear IDs on this sample with no
  code changes, which confirms M0.

## Milestone 1 implementation

Branch `feat/crosslink`. The search handles cleavable DSSO and DSBU crosslinks, end to end,
from the normal `sage` binary.

### Architecture

- **Core (`crates/sage`), generic hooks only.** Normal searches are unchanged: 568 of 568
  workspace tests pass with and without the feature.
  - `Scorer::score_offset_hypotheses` scores peptides against caller-supplied precursor
    masses and mass offsets. Each offset has its own neutral losses and preliminary
    shifts.
  - Fragments carry matched peak indices.
  - On exact score ties, isotope 0 is preferred.
  - `Scorer` derives `Clone`.
- **`crates/sage-xlink`:**
  - `doublets.rs`: signature doublets.
  - `linker.rs`: linkers and settings.
  - `search.rs`: the per-spectrum search.
  - `fdr.rs`: the discriminant and q-values.
  - `output.rs`: parquet output.
- **`sage-cli`:** a `crosslink` Cargo feature (off by default). When on:
  - `process_chunk` runs the crosslink search over MS2 spectra after the linear search,
    reusing the recalibrated scorer.
  - The run writes `crosslinks.sage.parquet` next to `results.sage.parquet`.
- **Per spectrum:**
  1. Doublets and the precursor give up to `max_pairs` chain-mass pairs. Precursor
     isotopes −1..3 are tried. A pair of observed chains is also accepted inside the
     isolation window.
  2. Each chain is looked up as a peptide carrying a spectrum-specific mass offset
     (partner + linker, on a linkable residue), with the stubs as optional neutral losses.
     Preliminary matching uses the stub shifts only.
  3. Chain candidates are combined, and the hyperscore is computed over the union of
     matched peaks.
- **FDR:**
  - A 12-feature LDA scores each CSM.
  - q = (TD − DD) / TT, made monotone, computed separately for intra- and inter-protein
    links.
  - It is applied at the CSM level and at the residue-pair level (unordered
    protein:position pairs).
  - Each residue pair is scored by a second, 2-feature LDA:
    - the best CSM discriminant score of the pair;
    - ln(1 + the number of the pair's CSMs with CSM q ≤ 5%).
  - A pair seen in a single spectrum stays reachable through its best score. At 1% FDR, 5
    such pairs pass on Beveridge (all correct) and 22 on the ribosome file (21 correct).

### Config

The `"crosslink"` block is valid only in builds with the feature. `prefilter` and
`wide_window` must be off.

```json
"crosslink": {
  "linker": "DSSO",
  "residues": "K", "protein_n_term": true,
  "isotope_errors": [-1, 3], "missing_charges": [3, 6], "isolation_half_width": 1.0,
  "min_chain_mass": 400.0, "max_pairs": 12,
  "preliminary_candidates": 5, "chain_candidates": 3,
  "min_chain_matched_peaks": 2, "output_q_value": 1.0
}
```

- `linker` is `"DSSO"`, `"DSBU"`, or a custom object.
- All other fields are optional; the values shown are the defaults.
- Keep the M0 monolink `mass_offset` mods in `variable_mods`, so that monolinks are not
  forced into crosslinks.

### Results, true FDR at 1% estimated CSM FDR

"Wrong" means the two chains come from different synthetic peptide groups, using the
IMP-X-FDR substring rule. Scripts are in `explore-data/crosslink/runs/`:
`analysis/eval_xl.py` and `analysis/xi_mzid.py`.

| Dataset / tool | CSMs (true FDR) | Residue pairs (true FDR) |
| --- | ---: | ---: |
| PXD014337 Beveridge Cas9, **Sage Plus** | 1828 (2.2%) | 180 (4.4%) |
| PXD014337, XlinkX (paper) | 481 (12.1%) | 181 crosslinks (29%) |
| PXD014337, MeroX Rise / Riseup (paper) | – | 125 (0.8%) / 168 (11%) |
| PXD029252 ribosome rep1, **Sage Plus** | 2322 (2.1%) | 529 (2.6%) |
| PXD029252 rep1, xiSEARCH (deposited mzid) | 1415 (2.1%) | 805 peptide pairs (3.4%) |

- PXD029252 was searched against the 671-sequence entrapment FASTA: 171 E. coli proteins
  plus 500 abundant human proteins. 28 of the 2322 CSMs involve a human chain. Those 28 are
  included in the 48 wrong CSMs.
- xiSEARCH searched only the 171 E. coli proteins, with 1% link-level and CSM-level FDR.
- PXD014337 is intra-protein only (Cas9). PXD029252 is 92% inter-protein.
- The estimate is slightly optimistic at the CSM level: 2.1–2.2% true FDR at a nominal 1%.
- **Residue-pair support feature.** Before it, pairs were scored by their best CSM alone.
  Results at a nominal 1% pair FDR, before → after (CSM results unchanged):

  | Dataset | Pairs | Correct | Wrong | True FDR |
  | --- | --- | --- | --- | --- |
  | Beveridge | 189 → 180 | 177 → 172 | 12 → 8 | 6.3% → 4.4% |
  | Ribosome | 475 → 529 | 466 → 515 | 9 → 14 | 1.9% → 2.6% |

  - Variants tried offline (script: `analysis/pair_variants.py`) did no better on both
    sets at once:
    - counting all CSMs;
    - adding distinct charge states;
    - adding distinct peptide forms.
  - Charges plus forms reached 2.1% on Beveridge, but cost 38 correct pairs and pushed the
    ribosome file to 3.6%.
  - A hard filter (3 or more supporting CSMs) was not adopted, because it discards
    single-spectrum pairs.
  - The remaining gap is mostly the CSM-level optimism carried through to the pairs.

### Where the CSM-level excess comes from

Wrong target CSMs at 1% estimated FDR fall into two kinds. Scripts:
`analysis/random_fdr.py`, `analysis/tail.py` and `analysis/csm_variants.py`.

- **Random:** at least one chain is outside the library (entrapment, crapome, or a
  missed-cleavage extension). Target-decoy is built to estimate these.
- **Cross-group:** both chains are library peptides, but from different groups.

| Dataset | True FDR at 1% | Random | Cross-group |
| --- | ---: | ---: | ---: |
| Beveridge | 2.2% | 0.2% | 2.0% (36 CSMs, 16 peptide pairs, 27 CSMs in recurring pairs) |
| Ribosome | 2.1% | 1.4% (entrapment estimate 1.5%) | 0.6% |

- **Decoys match false targets in bulk.** In score bands, TD − DD matches the count and the
  feature distributions of wrong TT on both files.
- **The Beveridge excess is not random.** The same cross-group pairs recur with strong
  scores on both chains. For example, MIAKSEQEIGK–LVDSTDKADLR (groups 1 and 10) appears
  in 12 CSMs, with both chain hyperscores at 40 or above. This is consistent with cross-contamination
  between groups, which no target-decoy model can see.
  - Even at 0.2% estimated FDR, Beveridge is 1.8% true.
- **Ribosome:** a stricter cutoff reaches the target. At 0.5% estimated FDR the true FDR is
  1.0%, but correct CSMs fall from 2399 to 1658.
- **Ranking variants tried offline, none calibrates better:**
  - per-chain doublet flags;
  - chain length and length per charge;
  - an isotope ≥ 2 flag;
  - the weak chain as the primary term.

  Weak-chain-primary cuts Beveridge depth by about 70%.
- **One variant adds depth.** "Weak-chain matched fragments per residue":
  - At 1% estimated FDR: +7% correct CSMs on Beveridge and +23% on the ribosome file, but
    true FDR rises to 2.6% and 2.5%.
  - At 0.2% estimated FDR: 2343 correct at 1.14% true on the ribosome file (98% of the
    current depth). On Beveridge it is 1189 correct at 2.0%.
  - Whether to adopt it, and at which cutoff, is left as a decision.

### Cost (`/usr/bin/time -v`, 16 cores, other jobs running)

| File | Linear: CPU, wall, RSS | Crosslink: CPU, wall, RSS |
| --- | --- | --- |
| PXD014337 (36.6k MS2) | 47 s, 16 s, 754 MB | 200 s, 37 s, 753 MB (crosslink step 21 s) |
| PXD029252 rep1 | 151 s, 38 s, 961 MB | 395 s, 41 s, 963 MB (crosslink step 21 s) |

- Peak RSS does not change, because chain lookups reuse the fragment index.
- CPU scales with the number of pair hypotheses:
  - `max_pairs` 4 costs half the CPU but loses 10% of CSMs.
  - Candidate caps of 5/3 match 10/5 at about 5% less CPU.
  - Dropping the intact-linker preliminary shift saved about 13% CPU with no loss.
- Wall times are noisy; the load average was 12–50 during these runs.
