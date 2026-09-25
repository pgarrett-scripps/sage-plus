# Intact glycopeptide search: design exploration

Status: exploration for a future beta. Nothing here is scheduled or released. Branch
`worktree-agent-a4498022d70b5b3e7`, based on `release/beta9` at `0ebb881`.

The README roadmap already lists "Glycopeptide search ... builds on mass offsets,
neutral losses, and motif modification sites". This document checks that claim against
the code. Most of the pieces do exist. The existing mass-offset search cannot carry a
glycan library as it stands, though: one offset per composition costs about one extra
search pass per composition. The recommended design therefore adds a single open "glycan
hypothesis" to retrieval and resolves the composition afterwards.

A prototype of the parts that do not depend on retrieval is in
`crates/sage/src/glycan.rs`: composition parsing, the composition library, oxonium
gating, and the Y-ion ladder, with tests. It is not wired into the search.

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
- **PXD005565** (yeast glycoproteome, same study). Yeast N-glycans are high
  mannose. Searching with the mammalian library makes every sialylated or fucosylated
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

## Appendix: prototype and measurements

- `crates/sage/src/glycan.rs`: `GlycanComposition` (parse, mass, containment),
  `GlycanLibrary` (`explain`, `indistinguishable_pairs`), `OXONIUM_IONS`,
  `oxonium_evidence`, `n_glycan_core_y_ions`, `y_ion_evidence`, and
  `n_glycan_composition_space`. Tests check the residue and oxonium masses, both
  notations, the isomeric and near-isobaric pairs, the gate, and the Y ladder. Run the
  ambiguity numbers with
  `cargo test -p sage-core --lib glycan::tests::ambiguity_report -- --ignored --nocapture`.
- Sequon fraction: `docs/explore/sequon_fraction.py`, a tryptic digest of
  `/mnt/data1/sage-plus-scientific/20260914/references/human.fasta` counting peptides
  that overlap an `N[^P][ST]` match in protein context (K/R not before P, 7–50 aa, ≤2
  missed cleavages, ≤5000 Da, carbamidomethyl C).
