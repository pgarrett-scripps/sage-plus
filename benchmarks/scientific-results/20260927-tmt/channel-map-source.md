# HYpro16 channel -> yeast ratio map (frozen 2026-09-27)

| Channel | Yeast : 1x | Human |
|---|---|---|
| 126 | 0 | 1 |
| 127N 127C 128N | 1/20 | 1 |
| 128C 129N 129C | 1/10 | 1 |
| 130N 130C 131N | 1/4 | 1 |
| 131C 132N 132C | 1/2 | 1 |
| 133N 133C 134N | 1 | 1 |

Sources checked:
1. Navarrete-Perea J, Gygi SP, Paulo JA. JASMS 2021 32(1):247-254, doi:10.1021/jasms.0c00299,
   PMC8210950 (full text via NCBI efetch db=pmc id=8210950, 2026-09-27).
   - Results: "Using the TMTpro16 format, the arrangement of yeast consisted of a single empty
     channel, followed by triplicates channels representing a 20:1,10:1, 4:1, and 2:1 decrease in
     the amount of yeast peptides from the final three channels (Figure 1A)."
   - Abstract: "TMTpro16-labeled yeast peptides in triplicate at 20:1, 10:1, 4:1, and 2:1 ratios."
   - Methods: "Yeast peptides were combined at a 0:(1/20):(1/20):(1/20):(1/10):(1/10):(1/10):(1/4):
     (1/4):(1/4):(1/2):(1/2):1:1:1 ratio across the 16 channels. Human peptides were combined at a
     1:1 ratio across all 16 channels." The printed list has 15 entries (one 1/2 is missing); the
     Results/abstract triplicate description fixes 16 = 1 empty + 5 x 3.
   - Channel order assumed to be the standard TMTpro order (126, 127N, 127C, ... 134N).
2. PXD020815_community_annotated.sdrf.tsv (PRIDE; annotates a18182/85/88 only, same pooled sample):
   TMT126 = 0; 127N/127C/128N = 1/20; 128C/129N/129C = 1/10; 130N/130C/131N = 1/4;
   131C/132N/132C = 1/2; 133N/133C/134N = 1. Agrees with (1).
3. 2020_Hyp.xlsx (deposit): hrMS2 = a18203, a18209, a18212; SPS-MS3 = a18207, a18210, a18213;
   RTS-MS3 = a18239, a18240, a18241; LumosFAIMS / LumosNoFAIMS / HFX isolation-width series.

Paper observed medians (1x : level), Lumos + FAIMS, Sequest, S/N >= 10 per channel and summed S/N >= 200:
- 2:1: hrMS2 ~1.81; SPS-MS3 and RTS-MS3 ~1.9
- 4:1: hrMS2 ~3.16; RTS-MS3 ~3.64 (SPS-MS3 not stated separately)
- 10:1: hrMS2 ~5.67; RTS-MS3 ~8.01 (SPS-MS3 not stated separately)
- 20:1: hrMS2 ~7.69; SPS-MS3 and RTS-MS3 ~13.55
The research write-up's SPS column for 4:1/10:1 (3.64, 8.01) is the RTS-MS3 value, not SPS.
