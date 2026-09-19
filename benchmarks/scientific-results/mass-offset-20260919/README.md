# Mass offset evaluation evidence

Prepared September 19, 2026 for `v0.1.0-beta.4`. This compares searching a
modification as a search-time mass offset against the equivalent database
expansion. It is engineering evidence on public spectra, not an independent
biological study, and it does not certify production calibration.

- 14 search jobs completed, none failed.
- Cost and index size use one public HEK file with the reviewed human reference,
  two measured searches per configuration.
- Localization uses the synthetic phosphopeptide libraries and their
  synthesis-defined sites.
- Entrapment uses the retained FDRBench paired target and entrapment peptide
  reference at the peptide level.

The searched executable is `4046e0291e36c8cd1a0db83b5b1cae04c7e9a64025d56b602b502de064252ba9`.
Every job records the executable, configuration, input, and output identities,
the command, thread budget, timings, and peak resident memory. The complete
records are under `/data/sage-plus-scientific/mass-offset-20260919`, with
`matrix.json` at its root and `job.json` in each job directory. The retained
[summary](summary.json) is the reduction the chapter reads; regenerate it with
`benchmarks/summarize_mass_offset.py`.

Spectra and references are external to this repository and identified by hash:

- `/data/sage-plus-scientific/20260914/converted/PXD001468/b1906_293T_proteinID_01A_QE3_122212.mgf`
- `/data/sage-plus-scientific/20260914/inputs/PXD000138/HCD_1.raw.-1.mgf`
- `/data/sage-plus-scientific/20260914/inputs/PXD000138/HCD_2.raw.-1.mgf`
- `/data/sage-plus-scientific/20260914/ptm-truth/synthetic.fasta`
- `/data/sage-plus-scientific/20260914/references/human.fasta`
- `/mnt/data1/sage-plus-scientific/mass-offset-20260919/inputs/paired.fasta`

Measured cost depends on the thread budget (8
threads here) and the isotope-error range, which multiplies the searched
windows. Accepted-identification differences between the two modes reflect
different searched hypothesis sets, not a correction of one by the other.
