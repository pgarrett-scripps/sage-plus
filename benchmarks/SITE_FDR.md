# Site-level FDR benchmark

Beta 13 adds `site_q_value` to `results.sage.protein-sites.parquet` and
`results.sage.ptm-sites.parquet`. Decoy PSMs that pass the same PSM and localization cutoffs
now form decoy sites. These decoy sites compete with target sites on the best supporting PSM
discriminant score. This benchmark compares two filters on the protein-site table: the new
`site_q_value <= 0.01`, and the old proxy `best_spectrum_q <= 0.01`, which is the minimum PSM
q-value of the site.

## Method

- Data: PXD007058 (Ferries et al. 2017), U2OS TiO2 phosphopeptide enrichment, Orbitrap Fusion.
  Two files: `SF_200217_U2OS_TiO2_HCD_nlETcaD_quad_OT_rep1` (rep1) and
  `SF_200217_U2OS_TiO2_HCD_nlEThcD_OT_rep2` (rep2).
- Database: human Swiss-Prot (20,416 proteins) plus the same number of entrapment proteins
  (`entrap_` prefix). Each entrapment protein shuffles the residues between K/R cleavage sites
  of its source and keeps an initiator Met. Sage generated reversed decoys for both halves.
- Search:
  - Trypsin with one missed cleavage, peptide length 7 to 50.
  - Fixed carbamidomethyl C. Variable phospho STY (at most 2) and oxidation M, at most 2
    variable modifications per peptide.
  - 10 ppm precursor and fragment tolerance, `prefilter: true`, `max_memory_gb: 7`.
  - The ion-trap EThcD scans (about 10%) were searched at the same 10 ppm, so they contribute
    little.
- Entrapment FDP is 2 x entrapment sites / all target sites (1:1 entrapment, combined
  estimator).
- Build: `b13/site-fdr` at commit 11aa304. Each search took under 1 minute and under 0.6 GB RSS.

Three `ptm_localization` settings:

| Setting | `psm_q_value` | `localization_q_value` |
|---|---:|---:|
| default | 0.01 | 0.01 |
| wide | 1.0 | 0.01 |
| loose | 1.0 | 1.0 |

## Results

Phospho PSMs at spectrum q <= 0.01: about 6,200 per file.

| Run | Target sites | Decoy sites competed | Old: min PSM q <= 0.01 | New: site q <= 0.01 | New only | Old only | Entrapment old / new | FDP old / new |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| rep1 default | 156 | 0 | 156 | 156 | 0 | 0 | 0 / 0 | 0 / 0 |
| rep1 wide | 156 | 0 | 156 | 156 | 0 | 0 | 0 / 0 | 0 / 0 |
| rep1 loose | 5,038 | 803 | 4,178 | 4,192 | 14 | 0 | 15 / 16 | 0.72% / 0.76% |
| rep2 default | 156 | 0 | 156 | 156 | 0 | 0 | 0 / 0 | 0 / 0 |
| rep2 wide | 0 | 0 | 0 | 0 | 0 | 0 | - | - |
| rep2 loose | 4,806 | 748 | 4,036 | 4,066 | 30 | 0 | 20 / 22 | 0.99% / 1.08% |

In rep1 and rep2 default, the 156-site count is a coincidence: the two sets share 105 sites.

## Findings

1. **Localization cutoff:** with the default 1% localization cutoff, site FDR changes nothing.
   - The localization FLR admits only the sites ranked above the first decoy-winning
     arrangement. That is 156 sites (localization q 1/125 and 1/126) out of about 6,200 phospho
     PSMs, or none in rep2 wide, where the minimum localization q is 1/63.
   - No decoy PSM clears that score, so no decoy site competes.
   - This bottleneck is in the localization FLR, which is unchanged from Beta 12, not in the
     site FDR.
2. **Without a localization cutoff:** site FDR is the only gate, and it changes little.
   - Site q <= 0.01 keeps every site the old proxy keeps, and adds 14 (rep1) and 30 (rep2).
   - Sites and decoy sites both collapse repeated PSMs, so the site-level decoy fraction at a
     given score is slightly lower than the PSM-level fraction.
   - Entrapment FDP stays near 1% under both filters: 0.72 to 1.08%, from 15 to 22 entrapment
     sites.
3. **What the numbers support:**
   - The old proxy was not badly anticonservative on this data.
   - The new column is the controlled quantity, and it now exists together with the decoy
     counts behind it (`decoy_protein_sites` in `run-summary.json`).
   - Entrapment checks site identity, not localization correctness.

## Caveats

- Sites are keyed by peptide position, not protein position. One residue seen in two peptides
  (missed cleavage, oxidation) counts as two sites, which inflates both target and decoy counts.
- The entrapment counts are small (15 to 22), so the FDP differences between filters are
  within noise.
- The localization FLR result (finding 1) needs its own fix before site FDR matters at default
  settings. Candidate causes:
  - Ranking by absolute target-decoy separation puts confident decoy wins at the top.
  - The +1-corrected estimate cannot fall below 1 / (targets before the first decoy).
