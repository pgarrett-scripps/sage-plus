# Intact glycopeptide search: design exploration

Status: milestone 1 is implemented on branch `feat/glyco` as the opt-in `glyco` Cargo
feature (section 7). Nothing here is released. Based on `release/beta9` at `0ebb881`.

The README roadmap already lists "Glycopeptide search ... builds on mass offsets,
neutral losses, and motif modification sites". This document checks that claim against
the code. Most of the pieces do exist. The existing mass-offset search cannot carry a
glycan library as it stands, though: one offset per composition costs about one extra
search pass per composition. The recommended design therefore adds a single open "glycan
hypothesis" to retrieval and resolves the composition afterwards.

Sections 1 to 6 are the design exploration. It used a prototype in
`crates/sage/src/glycan.rs` and a real-data harness,
`crates/sage-cli/examples/glyco_explore.rs` (both at commit `fe742af`). Milestone 1
moved the prototype to `crates/sage-glyco/src/composition.rs` and replaced the harness
with the `sage` binary itself. Section 7 has the implementation and its before/after
numbers.

## 1. How established tools do it

| Tool | Strategy | Glycan space | Fragmentation |
| --- | --- | --- | --- |
| MSFragger-Glyco (Polasky 2020) | Mass-offset search of the peptide: each glycan mass is an offset, and the fragment index is queried with shifted precursor windows. In labile mode, b/y ions are matched bare or carrying a remainder (HexNAc). Oxonium ions filter the offsets that are allowed. Y ions are extra fragments. Glycan composition is assigned afterwards in PTM-Shepherd with its own glycan FDR. | Curated composition list (N: the 182-composition human list from Byonic; O: a short mucin list) | HCD, sceHCD, EThcD |
| O-Pair (Lu 2020, MetaMorpheus) | Open HCD search finds the peptide and total glycan mass. The paired EThcD scan localizes O-glycans with a graph over site and glycan combinations (site probabilities, "level 1/1b/2/3" confidence). | O-glycan compositions, several per peptide | HCD-pd-EThcD pairs |
| pGlyco 2 / pGlyco3 | Glycan-first. Y ions and oxonium ions nominate glycan compositions and the peptide mass (Y1 = peptide + HexNAc). A narrow peptide search follows. Separate peptide, glycan, and glycopeptide FDR with decoy glycans. pGlyco3 adds modified-glycan support and EThcD site localization (pGlycoSite). | Structure database (GlycomeDB-derived), scored as compositions | sceHCD primarily, EThcD in pGlyco3 |
| Byonic | Exhaustive: peptide × glycan combinations as "rare" modifications with a glycan-aware fragment model and a two-dimensional posterior. | 182 human N-glycans, 6–78 O-glycan lists | HCD, EThcD |
| GlycoDecipher, StrucGP | Glycan-first with structure-level (not only composition) interpretation from Y-ion and B-ion topology. StrucGP targets structure-specific N-glycan evidence in HCD. | Structures | HCD |

Chemistry that every tool relies on:

- **Oxonium ions** (1+ glycan B ions) mark a spectrum as glyco: HexNAc 204.087, its
  fragments 186.076, 168.066, 144.066, 138.055, HexHexNAc 366.140, NeuAc 292.103 and
  274.092 (−H2O), NeuGc 308.098, Hex 163.060. HexNAc 204.087 plus at least one other
  is the usual gate. The sialic acid ions also decide between compositions (see risks).
- **Y ions** are the peptide carrying part of the glycan. In HCD of N-glycopeptides the
  glycosidic bonds break first. The spectrum shows Y0 (bare peptide), Y0+83.037 (0,2X
  cross-ring of the core HexNAc), Y1 (+HexNAc), the chitobiose and trimannosyl core
  ladder, and intact-minus-antenna ions, often at charge 2+ to 4+. Y1 is usually the
  most intense Y ion and pins the peptide mass independently of the glycan.
- **Peptide backbone ions**: in HCD, b/y ions have mostly lost the glycan or keep one
  HexNAc. The site is therefore weakly determined. For N-glycans the sequon settles it
  when a peptide has only one sequon. ETD and EThcD keep the glycan intact on c/z•
  ions, which makes O-glycan site localization possible.
- **Stepped-energy HCD (sceHCD)** gives oxonium, Y, and b/y ions in one scan and is the
  default for N-glycoproteomics. EThcD (or HCD-product-triggered EThcD) is used for
  O-glycans and site localization.

## 2. Glycopeptide FDR

A glycopeptide identification has three separable claims, and the tools control each:

1. **Peptide backbone**: ordinary target-decoy on reversed or shuffled peptides.
   MSFragger-Glyco uses its normal PSM FDR (with Percolator or MSBooster). pGlyco
   computes a peptide FDR from peptide decoys.
2. **Glycan composition**: whether the composition is right given a correct peptide.
   pGlyco 2/3 build decoy glycans by moving the Y-ion masses of each target composition
   by random amounts while keeping the precursor mass. The target composition competes
   with its decoys, and a glycan FDR is estimated from the fraction of decoy wins.
   MSFragger-Glyco (PTM-Shepherd glycan assignment) does the same: decoy compositions
   with randomized Y and oxonium masses, and a 1% glycan FDR applied after PSM FDR.
3. **Glycopeptide**: pGlyco2 combines them as FDR_gp ≈ FDR_p + FDR_g − FDR_p·FDR_g and
   filters on that. pGlyco3 estimates it jointly.
4. **Site** (O-glycans, or N-glycopeptides with two sequons): O-Pair site probabilities
   and pGlycoSite site probabilities. Most tools do not report site-level FDR in the
   target-decoy sense. They report probability thresholds. Sage Plus already has a
   competition-based localization q-value (`ModLocalization::set_competition_q_value`
   in `crates/sage/src/ptm.rs`), which could extend to glycosites.

Sage Plus reports spectrum, peptide, protein, and protein-group q-values
(`crates/sage/src/fdr.rs`: `picked_peptide`, `picked_protein`, ...). Peptide q-values
carry over unchanged for the backbone. A glycan q-value is new.

## 3. Mapping onto Sage Plus

### What exists and is reused

| Need | Existing code | Fit |
| --- | --- | --- |
| Restrict N-glycans to the sequon | Motif sites, `motif:N*-{P}-[ST]` (`crates/sage/src/motif.rs`, DOCS.md "Motif sites"). Evaluated against protein context with mirrored decoys. | Direct reuse. The same rule supplies placements and localization rules. |
| Test a modification at search time instead of expanding the index | `SearchMode::MassOffset` (`modification.rs`), `MassOffset` and `IndexedDatabase::mass_offset_sites` (`database.rs`), retrieval in `Scorer::matched_peaks_with_isotope` via `offset_query` and `IndexedQuery::page_search_shifted` (`scoring.rs`), placement scoring in `Scorer::score_hypothesis`, and identity assignment in `IndexedDatabase::materialize_mass_offsets`. | The machinery is right: precursor translated, fragments looked up bare and shifted, every placement scored, and offset peptidoforms materialized so FDR, LFQ, and site reports see ordinary peptides. The problem is the cardinality (below). |
| Glycan lost from b/y ions in HCD | `NeutralLossMode::Required` with a list of losses (`modification.rs`). `IonGroupSeries::losses` (`ion_series.rs`) turns them into alternative variants of one cleavage, and scoring picks at most one per cleavage. `MassOffset::fragment_shift` uses the smallest loss. | A glycan defined as `mass = G`, `neutral_losses = [G, G − 203.079]`, `required` gives exactly the MSFragger "labile" fragment model: b/y ions containing the site appear bare or +HexNAc, never intact. `fragment_shift` becomes 203.079, the same for every composition. |
| Activation-aware behaviour | `AcquisitionGroup { analyzer, activation }` on every `ProcessedSpectrum` (`spectrum.rs`), with `Activation::{Hcd, Ethcd, Etd, ...}`. Today it only selects fragment recalibration models (`mass_recalibration.rs`). | Available per spectrum. Ion kinds (`Parameters::ion_kinds`) and neutral-loss behaviour are still global, so EThcD needs new per-activation settings. |
| ETD/EThcD fragments | `Kind::C`, `Kind::ZDot` (just added). | Needed for milestone 3. |
| Site localization and site reports | `ptm::localize` (site-determining ions, AScore-style delta, probabilities, competition q-value) and the `ptm_sites` / `protein_sites` Parquet outputs (`crates/sage-cloudpath/src/parquet/sites.rs`). | Reused for multi-sequon peptides and later O-glycans. A glycan is one delta mass per composition. |
| Rescoring | 20-feature LDA (`crates/sage/src/ml/linear_discriminant.rs`, `FEATURE_NAMES`). | Glyco features are added here (fixed-size array, so a glyco feature set). |
| Isotope errors, precursor recalibration, chimeric search, LFQ | Unchanged. LFQ already quantifies materialized offset peptidoforms. | Reused. |

### Why one offset per composition does not scale

`benchmarks/MASS_OFFSET.md` measured the offset path on PXD001468: search time grows
linearly with offsets, about +5 to +7 s of search stage per offset on a 4.4 s baseline.
Each offset is another translated precursor window and another pass of shifted fragment
lookups. `MAX_MASS_OFFSETS` is 254 because `PreScore::offset` is a `u8`. A realistic
human N-glycan library is 150–350 compositions (the prototype's rule-based space has
352 without NeuGc, and the Byonic/MSFragger list has 182). As offsets, that is roughly
180–350× the unmodified search time, and above the cap. MSFragger survives this because
its index and its offset retrieval are organized differently. Porting that would change
Sage's core retrieval.

### Recommended retrieval: one open glycan hypothesis

Every glycan composition shares the same fragment model: bare b/y ions plus b/y+HexNAc.
The shifted lookup is therefore identical for all compositions, and only the precursor
window differs. That removes the need for one pass per composition:

1. **Glyco index.** Build the fragment index from peptides that carry a sequon site only.
   On the reviewed human FASTA used by the scientific corpus (20,416 proteins, 59,520
   sequons), 212,014 of 2,329,982 tryptic peptides (7–50 aa, ≤2 missed cleavages,
   ≤5000 Da) carry an in-context sequon N: **9.1%**
   (`docs/explore/sequon_fraction.py`). The glyco index is about a tenth
   of the normal one. Implementation: a peptide filter in database construction that
   keeps a peptide only if `compatible_sites` returns a site for the glycan
   specificity (the `database.rs` loop around `peptide.monoisotopic <=
   self.peptide_max_mass`).
2. **Gate.** Compute `glycan::oxonium_evidence` per spectrum before retrieval. Spectra
   that fail the gate skip the glycan hypothesis. In enriched samples the gate passes
   most spectra, and in unenriched samples it passes a small fraction.
3. **One hypothesis.** For a gated spectrum at neutral precursor mass M, add one
   retrieval pass with the precursor window `[M − G_max − tol, M − G_min + tol]` and
   fragment shift +203.079. This is `offset_query` with a window instead of a point,
   and the existing `page_search` / `page_search_shifted` loop is reused unchanged. On
   the glyco index this costs about one open search of roughly 420k peptides (targets
   and decoys) per gated spectrum. Sage already runs open searches of this kind.
4. **Resolve compositions at scoring.** For each preliminary candidate with peptide mass
   m, `GlycanLibrary::explain(M − m, tol, isotope_errors)` returns the compositions that
   fit. Each (composition, sequon site) becomes a candidate through the same path as
   `score_hypothesis`. The peptide is `with_mass_offset(site, definition)`, where the
   definition is a per-composition `ModificationDefinition` (mass G, required losses
   `[G, G − HexNAc]`). A preliminary hit with no explaining composition is dropped.
5. **Identity.** `MassOffsetAssignment` becomes an enum (`Offset { offset, site }` or
   `Glycan { composition, site }`), and `materialize_mass_offsets` and
   `resolve_peptide` gain a glycan arm. After materialization, glycopeptides are
   ordinary peptides for FDR, LFQ, and site reports.

The unshifted offset-0 hypothesis stays in place, so a gated spectrum still competes
against an unmodified explanation (in-source fragments, co-isolated peptides).

An alternative is pGlyco-style glycan-first retrieval: find Y1 in the spectrum, take
peptide mass = Y1 − HexNAc, and run a narrow closed search. It is cheaper but fails when
Y1 is weak or chimeric. It fits best as an accelerator, or as an extra feature, after
the open hypothesis works.

### What is new

| Component | Location | Notes |
| --- | --- | --- |
| Glycan composition library | `crates/sage/src/glycan.rs` (prototyped) | One composition per line, `HexNAc(4)Hex(5)Fuc(1)NeuAc(2)` (Byonic/MSFragger naming) or pGlyco `N(4)H(5)F(1)A(2)`. Config: `glycan: { library: "n_glycans.txt", sites: ["motif:N*-{P}-[ST]"], oxonium_min: 2 }`. Ship the curated human N list as a built-in, with its license checked before it is shipped. |
| Oxonium gate and features | `glycan::oxonium_evidence` (prototyped) | Count, intensity fraction, and a sialic-acid flag (292/274 seen). |
| Y-ion features | `glycan::n_glycan_core_y_ions`, `y_ion_evidence` (prototyped) | Matched count, intensity fraction, and Y0/Y1 anchor. Kept out of hyperscore at first, so b/y scoring and its calibration are unchanged, and fed to LDA. |
| Glyco retrieval hypothesis | `scoring.rs` (`matched_peaks`, `score_hypothesis`), `database.rs` | As above. |
| Glyco LDA features | `ml/linear_discriminant.rs` | Oxonium count and fraction, Y matched and fraction, anchored flag, glycan mass error, and the number of explaining compositions. |
| Glycan FDR | new `glycan_fdr` module, after peptide FDR | Decoy compositions: same precursor mass, Y-ion masses shifted by random offsets drawn once per target (seeded), scored with the same Y and oxonium evidence. Report `glycan_q` and a combined `glycopeptide_q`. |
| Output | results schema v3 (`schemas/results.sage.v*.parquet.schema`, `crates/sage-cloudpath/src/parquet/results.rs`) | Optional columns: `glycan_composition` (utf8), `glycan_mass`, `glycan_q`, `oxonium_ions` (int), `oxonium_intensity_pct`, `y_ions_matched`, `y_ion_intensity_pct`, `glycan_candidates` (number of compositions that fit). Peptide strings use ProForma `N[Glycan:HexNAc4Hex5Fuc1]`. Site reports reuse `ptm_sites` with the composition as the modification name. |

## 4. Staged plan

### Milestone 1: N-glycopeptides, HCD and sceHCD (smallest useful)

Scope: one glycan library, sequon-restricted, glyco index, oxonium gate, one open glycan
hypothesis, labile b/y model, Y-ion and oxonium features in LDA, glycan columns in the
results. Peptide-level FDR only, with composition reported as "best fit" plus
`glycan_candidates`. No glycan q-value yet. Documentation states that composition is
not FDR-controlled in this milestone.

Validate on:

- **PXD005411, PXD005412, PXD005413** (mouse brain, kidney, heart; pGlyco 2 paper, Liu
  et al. 2017, sceHCD). Also reanalyzed in the MSFragger-Glyco paper, so published
  pGlyco2 and MSFragger-Glyco results give an external comparison at the peptide and
  site level. PXD005553 and PXD005555 are the same study's further mouse sets.
- **PXD005565** (fission yeast, *S. pombe*, same study; the runs pool unlabeled,
  15N- and 13C-labeled cells, so only about a third of precursors have natural
  isotopes). Yeast N-glycans are high mannose. Searching with the mammalian library makes every sialylated or fucosylated
  complex composition a known-false glycan assignment. This gives an empirical glycan
  FDP for free and is the main tool for milestone 2.
- Peptide-level calibration: the existing FDRBench entrapment tooling
  (`benchmarks/run_fdrbench_validation.py`) applied to the glyco index.
- Regression: byte-identical `results.sage.parquet` for configurations without a glycan
  section, as in `MASS_OFFSET.md`.

Effort: about 3–4 weeks of focused work (glyco index filter 2–3 days; retrieval
hypothesis and assignment enum 5–7 days; features and LDA 3 days; output schema and
docs 2–3 days; data acquisition and conversion, runs, and comparison 5 days).

### Milestone 2: glycan FDR

Decoy compositions, `glycan_q`, combined glycopeptide q, calibrated on PXD005565 (yeast
known-false compositions) and on a NeuGc entrapment in human data. Humans cannot make
NeuGc, so NeuGc compositions in human tissue are almost all false. Dietary incorporation
is low but not zero, so treat this as an upper bound. About 1.5–2 weeks.

### Milestone 3: EThcD and O-glycans

Per-activation fragment model: `AcquisitionGroup::activation` selects ion kinds (b/y for
HCD; b/y/c/z• for EThcD) and the glycan loss mode (labile on b/y, retained on c/z•). This
needs `Parameters::ion_kinds` and the loss mode to become per activation, which also
affects the preliminary index (index both kinds and filter at scoring). O-glycans:
S/T sites, several glycans per peptide (at least a sum-composition hypothesis with O-Pair
style placement enumeration), and site probabilities through `ptm::localize` extended to
multiple glycans. Validate on:

- **PXD017646** (Riley et al. 2020, "Optimal dissociation methods differ for N- and
  O-glycopeptides"; HCD, sceHCD, EThcD on the same samples). This isolates the
  activation question.
- **PXD009476** (EXoO, Yang et al. 2018). OpeRATOR cleaves N-terminal to O-glycosylated
  S/T, so the site is the peptide's first residue. This gives ground truth for site
  localization.
- **PXD011533** (Riley et al. 2019, AI-ETD N-glycoproteome) for ETD-family N-glycan
  scoring.

About 4–6 weeks. Multi-glycan O-site enumeration is the open-ended part.

All accessions above were confirmed to exist through the PRIDE API on 2026-09-24. Titles
were checked. File lists and instrument methods were not.

## 5. Risks

1. **Compositions share masses.** From the prototype `ambiguity_report` test, over a
   rule-based N-glycan space:

   | Space | Compositions | Ambiguous at 0.02 Da, no isotope shift | Ambiguous with ±1 isotope error |
   | --- | --- | --- | --- |
   | Human (no NeuGc) | 352 | 0 | 236 (67%) |
   | Mammalian (with NeuGc) | 766 | 426 (56%) | 677 (88%) |

   Two sources drive this. First, NeuAc + Hex and NeuGc + Fuc are exactly isomeric.
   Second, Fuc2 lies 17 mDa from NeuAc plus one neutron, the classic precursor
   monoisotopic-peak error. With Sage's default `isotope_errors`, most human
   compositions have a twin. Mitigations: restrict glyco isotope errors to `[0, 1]` and
   report `glycan_candidates`. Break ties with evidence: sialic acid oxonium ions
   (292/274) versus their absence, and Y ions carrying two Fuc. That is exactly what
   glycan FDR has to control, and why it is milestone 2 rather than an afterthought.
   Fucose migration in HCD can create misleading Fuc-containing Y ions. Keep the
   core-fucose Y ions (Y1+Fuc) as the main fucose evidence.
2. **O-glycan site ambiguity.** HCD gives almost no O-site information. S/T-rich mucin
   domains produce many equivalent placements, and several glycans per peptide multiply
   them. Do not report O-sites from HCD. Require EThcD and report O-Pair-style
   confidence levels.
3. **In-source fragments.** Glycopeptides lose sialic acid and fucose in the source. The
   fragment is a real precursor that matches a real, smaller composition on the right
   peptide at the same retention time. This is a correct-peptide, wrong-glycoform
   error, and target-decoy does not see it. Mitigation: flag an identification when the
   same peptide and site is identified with a strictly larger composition within a small
   RT window (about 0.1 min), and prefer the larger form in site and glycoform summaries.
   LFQ traces can confirm co-elution.
4. **Open-window false positives.** A wide precursor window lets a non-glyco or chimeric
   spectrum match a sequon peptide plus some composition. The mitigations: the oxonium
   gate, the Y0/Y1 anchor feature, the unmodified hypothesis competing, and decoys
   drawn from the same glyco index (so target-decoy sees the same window).
5. **Spectrum preprocessing.** `max_peaks` (default 150) and deisotoping must keep
   low-m/z oxonium ions and high-charge Y ions. Check this on real sceHCD data before
   tuning features. Y ions can also be confused with b/y ions of the backbone at
   higher charge, and `max_fragment_charge` needs to allow 3+ Y ions.
6. **Cost for unenriched data.** When few spectra pass the gate, the cost is small. For
   enriched data it is one open search of the glyco index per spectrum. Measure this in
   milestone 1 before committing to the open-window design over glycan-first retrieval.

## 6. Real-data test (2026-09-25)

One run from each study, searched with the `glyco_explore` harness (commit `fe742af`).
Outputs are in
`/mnt/data1/explore-data/glyco/runs/<run>/` (`summary.tsv`, `glyco_psms.tsv`,
`config.json`).

### What the harness does

1. Reads the Thermo .raw directly, keeps the 300 most intense peaks, and deisotopes.
2. Gate: the HexNAc oxonium ion 204.087 plus at least one other oxonium ion, at 20 ppm.
3. Index: tryptic peptides with 2 missed cleavages, 5–50 aa and 500–5000 Da, keeping
   only those with an `N-{P}-[ST]` sequon (`Peptide::compatible_sites`, then
   `Parameters::build_from_peptides`). The index has target and decoy peptides.
4. Search: the stock `Scorer` searches gated spectra with an asymmetric Da precursor
   window that spans the glycan library, from −(Gmax + 1.5) to −(Gmin − HexNAc) + 0.1.
   A labile HexNAc mass offset on the sequon (neutral loss of 203.079, optional) lets
   fragments carry either nothing or the innermost HexNAc. This is milestone 1's single
   open hypothesis, built from existing parts.
5. Assignment: takes the best-ranked of 5 candidates whose precursor delta a library
   composition explains within 20 ppm (isotope 0 or +1, optionally plus one NH3
   adduct). Ties among compositions are broken in order by:
   1. sialic-acid oxonium consistency (274/292 for NeuAc, 308 for NeuGc)
   2. core Y-ion matches
   3. no adduct, then isotope 0, then mass error
6. FDR: the stock LDA (`score_psms`, Da mode) and spectrum q-values. This is
   peptide-level only. There is no glycan FDR yet.

### Gate and Y-ion ladder

"Ladder" is a peptide-free measure: the longest chain of peaks, at one charge, separated
by HexNAc, Hex or Fuc residue masses above 500 Da. Ungated spectra are the background.

| Run | MS2 | Gated | Ladder ≥4 steps, gated | Ladder ≥4 steps, ungated | Y0/Y1 seen in 1% PSMs |
| --- | --- | --- | --- | --- | --- |
| MouseBrain-Z-T-1 | 53,102 | 96.2% | 34.7% | 4.2% | 89–93% |
| cwq_mix2-1_726 (yeast) | 59,609 | 73.3% | 49.1% | 21.0% | 84–95% |
| Human HCD, not glyco-enriched (control) | 60,676 | 0.002% (1 spectrum) | — | — | — |

The gate is specific. On an ordinary human HCD run it passed 1 of 60,676 spectra. Both
glyco datasets are enriched, so nearly everything passes and the gate saves little time
here. The ladder of three or fewer steps is not specific (57–97% of ungated spectra show
two steps). Four or more steps separates well in mouse. In yeast the background is
higher: 21% of ungated spectra show a ladder of four or more steps, which suggests the
gate misses some glycopeptides. After identification, the peptide-anchored Y0 or
Y1 (peptide+HexNAc) ion is present in about 90% of glycoPSMs, which supports Y-ion
evidence as the main composition feature.

### Mouse brain (PXD005411, run 1 of 5, mouse Swiss-Prot)

| Glycan list | NH3 adduct | GlycoPSMs at 1% | Decoys | Unique peptide+glycan | Unique peptides | Several compositions fit | Isotope +1 | Time |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 182 (FragPipe mouse) | no | 4,850 | 47 | 3,007 | 1,251 | 12.9% | 27.3% | 169 s |
| 1670 (pGlyco mouse large) | no | 5,546 | 54 | 3,566 | 1,282 | 33.6% | 27.6% | 186 s |
| 1670 | yes, sialic first | **6,168** | 60 | 3,785 | 1,347 | 47.4% | 29.0% | 184 s |

The 182-list run with the adduct lost its LDA fit and fell back to ranking by
hyperscore, which gave 1,966 PSMs. The fallback needs a more robust
discriminant, such as regularization or the Poisson score, before this is productized.

Top compositions (1670 list with adduct): HexNAc(2)Hex(5) 874, HexNAc(2)Hex(8) 422,
HexNAc(2)Hex(9) 379, HexNAc(2)Hex(6) 342, HexNAc(4)Hex(3)Fuc(1) 235, HexNAc(2)Hex(7)
227, HexNAc(5)Hex(3)Fuc(1) 200, HexNAc(5)Hex(3)Fuc(1)NeuAc(1) 167. High mannose makes up
38% and fucosylated compositions 47%. Brain is known to be rich in fucosylated
complex glycans.

NeuGc works as an entrapment here, because brain carries almost none. It is assigned to
2.3% of PSMs (139). That falls to 1.7% when two or more Y ions are required.
MSFragger/PTM-Shepherd report 1.1% total entrapment for brain.

Moving to sialic-first ranking with the adduct changed the composition mix. NeuAc rose
from 23.6% to 38.4% and Fuc≥2 fell from 1,527 to 1,134. Two Fuc weigh 1.02 Da more
than one NeuAc, so with isotope +1 allowed the two compositions are indistinguishable
by precursor mass. That is the same Fuc-for-NeuAc error PTM-Shepherd attributes to
pGlyco3. Oxonium evidence is the only thing that separates them.

### Yeast (PXD005565, run 1 of 3)

The provided `yeast_sp.fasta` is *S. cerevisiae*. The sample is *S. pombe*
(Polasky et al. 2022), and the *S. cerevisiae* search gave **0** glycoPSMs at 1%. The
runs below use reviewed *S. pombe* (UniProt, 5,129 entries) plus mouse Swiss-Prot as a
peptide entrapment, following the published design.

| Glycan list | NH3 | GlycoPSMs at 1% | Non-high-mannose (glycan FDP) | Non-HM with Y ≥2 | Isotope +1 | Adduct |
| --- | --- | --- | --- | --- | --- | --- |
| 1670 + HexNAc(2)Hex(3–20) | no | 1,495 | **52.7%** | — | 23.5% | — |
| 1670 + HexNAc(2)Hex(3–20) | yes | 1,682 | **28.5%** | 17.4% | 25.3% | 36.6% |
| HexNAc(2)Hex(3–20) only (18) | yes | 1,437 | 0 by construction | — | 27.3% | 38.6% |

- **Ammonium adducts dominate the errors.** Without adducts, the top false compositions
  were HexNAc(4)Hex(3)Fuc(2)NeuAc(1) and HexNAc(4)Hex(5)Fuc(3)NeuGc(1). Each is exactly
  HexNAc(2)Hex(n) + NH3 within 7 ppm (a 17.053 versus 17.027 Da difference). An adduct
  explains 37% of assignments in yeast and 33% in mouse, so milestone 1 must include
  the adduct.
- **Y ions separate the remaining errors.** Of the assignments with 0–1 Y-ion matches,
  82% (238/290) are non-high-mannose. At ≥4 Y ions the rate is 12.9%; at ≥7 it is
  6.2%. The remaining top false assignments are HexNAc(4)Hex(5)Fuc(4) and
  HexNAc(6)Hex(3–4). Those are the kinds of compositions a decoy-composition glycan
  FDR (milestone 2) has to catch.
- Labeled precursors (15N/13C, about two thirds of the pool) cannot be explained by
  natural-isotope masses. Some of them probably survive as wrong compositions.

### Comparison with published numbers

Published counts are per dataset. The per-run figures divide by the number of runs
(5 brain, 3 yeast).

| | Mouse brain glycoPSMs/run | Yeast glycoPSMs/run | Yeast non-yeast glycans |
| --- | --- | --- | --- |
| MSFragger-Glyco + PTM-Shepherd, 1% peptide and glycan FDR, 1670 list, NH3 | ≈8,990 (44,931/5) | ≈2,745 (8,234/3) | 3.8% |
| pGlyco3, same lists, NH3 | ≈5,400–6,400 (40–66% fewer than above) | ≈1,895 (5,684/3) | 7.5% |
| pGlyco3, no NH3 | — | ≈1,135 (3,405/3) | 7.3% |
| This harness, 1% peptide FDR only, 1670 list, NH3 | 6,168 | 1,682 | 28.5% (17.4% at Y≥2) |

Sources: Polasky et al., MCP 2022 (PMC8933705), Tables 1 and 4; Polasky et al., Nat
Methods 2020 (PMC7606558); Liu et al., Nat Commun 2017 (PMC5585273).

### Conclusions

1. **Retrieval is fine.** The sequon index plus the open window and labile HexNAc
   offset is already at pGlyco3's per-run glycoPSM level in brain. It runs in about
   3 minutes per 2 GB run on 16 cores, and needs no core changes. The open-window
   design holds up. Glycan-first retrieval is not needed for milestone 1.
2. **Composition assignment is the gap.** Assignment is about 4–8× worse than the
   published tools on the yeast entrapment. It needs the three pieces the harness
   lacks:
   1. NH3 adducts, which already gave a 2× improvement
   2. Y-ion and oxonium evidence as scored features, not tie-breaks
   3. a glycan-level FDR with decoy compositions (milestone 2)
   The yeast Y-ion bins show a combined score should separate true from false
   compositions. Glycan FDR should move up to be part of milestone 1, not follow it.
3. **The candidate budget limits recall.** Only 16,147 of 50,633 searched brain spectra
   had an explainable candidate among the top 5. Composition-aware rescoring over more
   candidates is the next recall lever.

## 7. Milestone 1 implementation (2026-09-25)

### How to run it

Build with `cargo build --release -p sage-cli --features glyco` and add a block to a
normal config:

```json
"glyco": {
  "glycan_files": ["Mouse_N-glycans-1670-pGlyco.glyc"],
  "high_mannose": false,
  "ammonium_adducts": 1
}
```

Every field is optional. Without a glycan list the built-in mammalian composition space
is used. The normal search runs as before and writes its usual outputs. The glyco pass
adds `glyco.sage.parquet`, with one row per explained candidate (schema
`schemas/glyco.sage.v1.parquet.schema`, marked experimental). A build without the
feature rejects a config that has the block. Without the block nothing changes.

### Where the code lives

1. **Core (`crates/sage`), two hooks.**
   1. `IndexedDatabase::offsets_only` skips the unshifted hypothesis when every
      precursor carries a mass offset.
   2. `glycan.rs` moved out to `sage-glyco`.
   Normal searches are bit-identical: 8,674 mouse and 7,940 yeast PSMs, before and
   after.
2. **`crates/sage-glyco`.** Config, sequon index build, oxonium gate, glycan
   explanation, Y-ion and oxonium evidence with decoy twins, two-level FDR, and Parquet
   output.
3. **`sage-cli`.** Behind the `glyco` feature. The glyco index is built next to the
   main one. Each spectrum chunk is searched once more with the glyco scorer, in the
   same single pass over the files, and the spectra are released as usual. The
   candidates carry only per-class ion counts, so FDR runs at the end without keeping
   spectra.

### Glycan FDR

Each candidate explanation (composition, isotope error, NH3 adducts) gets a decoy
twin. The twin has the same precursor mass and the same core Y ions, but its other
Y ions and its sialic-acid oxonium ions are moved by a pseudo-random 4.5–12.5 Da. Each
explanation is scored by a per-class log-likelihood ratio against the spectrum's own
random match rate. The class rates are refitted once from confident PSMs. The best of
all targets and twins wins. Two choices made the FDR bite:

1. **Exact ties are a coin flip, not a target win.** When nothing beyond the shared
   core matches, a target and its twin score the same. The first version gave such
   ties to the target, so only 2 of 5,355 mouse candidates ever lost to a decoy, and
   the yeast non-yeast rate stayed at 8.8%.
2. **Competition on score plus lead.** Target-decoy competition ranks by the winner's
   score plus its lead over the best explanation of the other label. Both terms are
   symmetric when a wrong composition swaps with its twin. On yeast this passes 781
   glycoPSMs at 2.4% non-yeast glycans, against 708 at 1.3% for the score alone.
   The minimum of the two passed 875 but at 5.8%, so it was rejected.

Peptide FDR is an LDA over 20 standardized features (22 from milestone 2). Glyco features include the core,
best-Y and Y-intensity fractions, Y0/Y1 anchoring, oxonium count and intensity, the
best glycan score, precursor error and ambiguity. If the fit fails it is retried with
ridge 1e-3, 1e-2 and 1e-1, and then falls back to the hyperscore, so a run always
produces q-values.

### Before and after

Same runs as section 6. All numbers are at 1% peptide FDR and 1% glycan FDR unless noted.
Machine: 16 cores, shared with other jobs, so wall times vary by about ±50% between
repeats. Glyco search time is the time spent in the glyco scorer.

| | Mouse glycoPSMs | Yeast glycoPSMs | Yeast non-yeast glycans | Mouse wall / peak RSS | Mouse glyco search |
| --- | --- | --- | --- | --- | --- |
| Plain Sage, no glyco block | — | — | — | 18–25 s / 3.9 GB | — |
| Section 6 harness (peptide FDR only, 300 peaks) | 6,168 | 1,682 | 28.5% | about 184 s, separate program | — |
| First integration (ties to target, Met-ox in glyco index) | 5,353 | 1,160 | 8.8% | 580 s / 5.07 GB | 489 s |
| Milestone 1 (defaults, 150 peaks) | **5,458** | **781** | **2.4%** | 97–154 s / 4.4 GB | 88–129 s |
| Milestone 1, `max_peaks: 300` | 6,024 | 900 | 4.2% | 239 s / 4.5 GB | 191 s |
| MSFragger-Glyco + PTM-Shepherd (per run) | ≈8,990 | ≈2,745 | 3.8% | — | — |
| pGlyco3 (per run) | ≈5,400–6,400 | ≈1,895 | 7.5% | — | — |

Yeast wall time is 106–205 s against 20–30 s plain, with 5.6 GB against 4.9 GB peak
RSS.

What moved the numbers:

1. **Speed: glyco index without user variable mods.** The open glyco window spans
   about 3,000 Da, so every peptide variant is a candidate for every spectrum. Met-ox
   with `max_variable_mods: 2` gave 634,667 sequon peptides. Dropping it gives 393,009,
   and glyco search time fell from 489 s to under 130 s. Mouse glycoPSMs held (5,346
   against 5,353) even though the stricter tie rule landed in the same run. `glyco.variable_mods: true` turns them back on.
2. **Error: the two FDR choices above.** They took yeast from 8.8% to 2.4% non-yeast
   glycans, and mouse went up by 105 glycoPSMs.
3. **Recall: peaks.** `max_peaks: 300` recovers most of the harness count in mouse
   (6,024). In yeast it gives 900 glycoPSMs, but at 4.2% non-yeast glycans. The
   default stays at Sage's 150.

### Gaps

1. **Yeast recall is well below MSFragger-Glyco** (781 against about 2,745 per run).
   Only 1,553 yeast candidates pass peptide FDR, and the glycan FDR then drops half of
   them. Most decoy winners (256 of 314) are coin-flip ties, where no ion beyond the
   core tells target from twin. Composition-aware rescoring of more
   than the first explained candidate is the next lever. So is an EM refit of the Hex
   ladder for long yeast mannans.
2. **Mouse is at pGlyco3's level, not MSFragger's.** The candidate budget from section 6
   still applies.
3. **Wall-time noise.** The machine was shared during these runs. Repeat on an idle
   machine before quoting a single number.

## 8. Milestone 2: more candidates per spectrum (2026-09-25)

Goal: recover yeast recall without raising the wrong-glycan rate above about 4%, and
keep mouse at or above 5,458.

### What changed

1. **Top-N rescoring.** Each spectrum keeps up to `explain_candidates` (default 5)
   explained peptide candidates, not just the first. Every candidate gets its glycan
   assignment from the initial score model. The one with the best glycan score
   represents the spectrum, and ties keep the better-ranked peptide. Two peptide
   features were added (22 in total): the log peptide rank, and the lead of the
   chosen candidate over the next one.
2. **Deeper peptide list.** `report_candidates` went from 5 to 50. Most top-ranked
   sequon peptides have no glycan that explains the precursor delta, so a deep list
   is what lets more spectra reach an explained candidate. Yeast explained spectra
   rose from 9,114 to 22,050 and mouse from 16,440 to 36,024, with no measurable time
   or memory cost.
3. **Trunk Y ladder for oligomannose.** For compositions with at most two HexNAc and
   no sialic acid, every HexNAc(2)Hex(n) Y ion is generated, not just the core and the
   ions near the precursor. This breaks target/twin ties on long yeast mannans. Applying
   the ladder to every composition gave 920 yeast glycoPSMs but at 4.8% non-yeast
   glycans, because complex compositions then shared the Hex ladder. So it is limited
   to oligomannose parents.

### Before and after

Same runs, settings and machine as section 7. Each row adds to the one above it.

| | Mouse glycoPSMs | Yeast glycoPSMs | Yeast non-yeast glycans | Yeast decoy winners | Yeast wall / peak RSS |
| --- | --- | --- | --- | --- | --- |
| Milestone 1 | 5,458 | 781 | 2.4% | 314 | 106–205 s / 5.6 GB |
| Top-5 rescoring | 5,491 | 791 | 2.7% | 267 | 188 s / 5.6 GB |
| Trunk ladder on all compositions (rejected) | 5,474 | 920 | 4.8% | 228 | 236 s / 5.5 GB |
| Trunk ladder on oligomannose only | 5,490 | 881 | 2.5% | 246 | 160 s / 5.6 GB |
| `report_candidates: 50` (**new default**) | **5,909** | **992** | **2.6%** | 219 | 139 s / 5.6 GB |

Mouse wall time is 130 s at 4.5 GB peak RSS with the new defaults. Plain searches
(no glyco block) are unchanged: 8,674 mouse and 7,940 yeast PSMs.

### What is left

1. **Yeast is still at about a third of MSFragger-Glyco** (992 against about 2,745).
   Of 41,401 oxonium-gated yeast spectra, 22,050 get an explained candidate but only
   1,579 pass peptide FDR.
   The likely limit is retrieval (not yet measured): the right sequon peptide is often not in the top 50 by
   peptide fragments alone. Glycan-aware candidate ranking, for example Y1-shifted
   fragments in the index, is the next lever.
2. **Ties remain.** 219 decoy winners pass peptide FDR in yeast. Most of them match
   only core ions. Oxonium intensity ratios (Hex versus HexNAc) could separate
   mannans from complex glycans without relying on Y ions.

## 9. Milestone 3: Y ions in the fragment index (2026-09-25)

Milestone 2 left 22,050 explained yeast spectra, of which only 1,579 passed peptide FDR.
That pointed at retrieval. The glyco window spans about 3,000 Da, so the peptide mass is
free, and preliminary retrieval ranked every sequon peptide in it on b and y ions alone.
In sceHCD spectra, the strongest peaks are often the Y ions, which carry the whole peptide
mass.

### What changed

1. **Core hook.** `Parameters::build_from_peptides_with_extra_fragments` appends
   caller-supplied masses to each peptide's preliminary fragments. They only count toward
   preliminary retrieval, and full scoring is unchanged. `build_from_peptides` passes no
   extra masses, so normal searches are identical: 8,674 mouse and 7,940 yeast PSMs,
   re-checked after the change.
2. **Six Y ions per sequon peptide** (`glyco.index_y_ions`, on by default): Y0, Y1,
   Y1+Fuc and HexNAc(2)Hex(1..=3). The HexNAc hypothesis also looks up every indexed
   fragment at +HexNAc, so Y2 and Y2+Fuc count too. Decoys are reversed peptides with
   the same mass, so they get the same Y ions, and target-decoy competition stays fair.

### Before and after

Release builds, same runs and machine as section 8, with `report_candidates: 50`.

| | Mouse glycoPSMs | Yeast glycoPSMs | Yeast non-yeast glycans | Yeast peptide-FDR passes | Yeast wall / peak RSS | Mouse wall / peak RSS |
| --- | --- | --- | --- | --- | --- | --- |
| Milestone 2 | 5,909 | 992 | 2.6% | 1,579 | 139 s / 5.6 GB | 130 s / 4.5 GB |
| Milestone 3, Y-ion index (**default**) | **6,860** | **1,229** | **2.4%** | 1,620 | 109 s / 5.6 GB | 97 s / 4.5 GB |

Index cost: glyco fragments grow by 14%, from 22.7 to 25.9 million in yeast and from 17.0
to 19.3 million in mouse. That is about 19 MB and 14 MB at 6 bytes per fragment. Peak
RSS did not change measurably. Wall time went down, within machine noise.

Tried and not adopted (fast-release builds, yeast):

1. **Y ladder to HexNAc(2)Hex(6)**, nine ions: 1,227 yeast, 6,835 mouse. No gain for 6%
   more fragments.
2. **`report_candidates: 100`**: 1,217. **`explain_candidates: 20`**: 1,237. Retrieval
   depth is no longer the limit.

### What is left

1. **Yeast is at 45% of MSFragger-Glyco** (1,229 against about 2,745). 23,195 yeast
   spectra are explained, but only 1,620 pass peptide FDR, and deeper candidate lists no
   longer help. So the limit is now how well the peptide model separates correct
   candidates from wrong ones. Full scoring still ranks on b and y ions only. Y-ion counts
   in the hyperscore, or a glyco-specific peptide score, is the next lever.
2. **Ties.** 143 yeast decoy winners pass peptide FDR, down from 219.

## Appendix: code and measurements

- `crates/sage-glyco/src/composition.rs` (the prototype's `glycan.rs`):
  `GlycanComposition` (parse, mass, containment), `GlycanLibrary` (`explain`,
  `indistinguishable_pairs`), `OXONIUM_IONS`, `oxonium_evidence`,
  `n_glycan_core_y_ions`, `y_ion_evidence`, and `n_glycan_composition_space`. Tests
  check the residue and oxonium masses, both notations, the isomeric and near-isobaric
  pairs, the gate, and the Y ladder. Run the ambiguity numbers with
  `cargo test -p sage-glyco --lib composition::tests::ambiguity_report -- --ignored --nocapture`.
- `crates/sage-glyco/examples/glyco_summary.rs`: counts, non-high-mannose share and top
  compositions of a `glyco.sage.parquet`, optionally dumping every row as TSV.
  `cargo run --release -p sage-glyco --example glyco_summary -- glyco.sage.parquet [rows.tsv]`.
- The section 6 harness, `glyco_explore`, is at commit `fe742af`. Glycan lists come
  from FragPipe `tools/Glycan_Databases` (not vendored).
- Sequon fraction: `docs/explore/sequon_fraction.py`, a tryptic digest of
  `/mnt/data1/sage-plus-scientific/20260914/references/human.fasta` counting peptides
  that overlap an `N[^P][ST]` match in protein context (K/R not before P, 7–50 aa, ≤2
  missed cleavages, ≤5000 Da, carbamidomethyl C).
