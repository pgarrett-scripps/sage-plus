"""Count tryptic peptides carrying an N-glycosylation sequon (N-{P}-[ST]).

Sequon is evaluated against protein context, as Sage Plus motif sites are.
"""
import re
import sys

MASS = dict(zip("ACDEFGHIKLMNPQRSTVWY", [71.03711, 103.00919, 115.02694, 129.04259, 147.0684, 57.02146,
    137.05891, 113.08406, 128.09496, 113.08406, 131.0405, 114.04293, 97.05276, 128.05858, 156.1011,
    87.03203, 101.04768, 99.06841, 186.07932, 163.06332]))
H2O = 18.010565


def proteins(path):
    name, seq = None, []
    for line in open(path):
        if line.startswith(">"):
            if name:
                yield "".join(seq)
            name, seq = line, []
        else:
            seq.append(line.strip())
    if name:
        yield "".join(seq)


total = set()
sequon = set()
sites = 0
for protein in proteins(sys.argv[1]):
    positions = {m.start() for m in re.finditer(r"(?=N[^P][ST])", protein)}
    sites += len(positions)
    cuts = [0] + [i + 1 for i, aa in enumerate(protein[:-1]) if aa in "KR" and protein[i + 1] != "P"] + [len(protein)]
    for i in range(len(cuts) - 1):
        for mc in range(3):
            j = i + 1 + mc
            if j >= len(cuts):
                break
            start, end = cuts[i], cuts[j]
            peptide = protein[start:end]
            if not 7 <= len(peptide) <= 50 or any(aa not in MASS for aa in peptide):
                continue
            mass = sum(MASS[aa] for aa in peptide) + H2O + 57.02146 * peptide.count("C")
            if mass > 5000:
                continue
            total.add(peptide)
            if any(start <= p < end for p in positions):
                sequon.add(peptide)
print(f"sequon sites in proteome: {sites}")
print(f"tryptic peptides (7-50 aa, <=2 missed, <=5000 Da): {len(total)}")
print(f"with an in-peptide sequon N: {len(sequon)} ({100 * len(sequon) / len(total):.1f}%)")
