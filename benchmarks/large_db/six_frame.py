"""Translate genome FASTA records in all six reading frames.

Each frame is split at stop codons and at codons containing an ambiguous base,
and every open segment of at least MIN_LENGTH residues is written as one
protein entry named by its genome coordinates.
"""

import argparse
import itertools
import sys

BASES = "TCAG"
AMINO = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
CODONS = {
    "".join(codon): aa
    for codon, aa in zip(itertools.product(BASES, repeat=3), AMINO)
}
COMPLEMENT = str.maketrans("ACGTN", "TGCAN")
MIN_LENGTH = 7


def read_fasta(path):
    name, chunks = None, []
    with open(path) as handle:
        for line in handle:
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(chunks)
                name, chunks = line[1:].split()[0], []
            else:
                chunks.append(line.strip().upper())
    if name is not None:
        yield name, "".join(chunks)


def segments(sequence):
    """Yield (start, end, peptide) for open segments of one frame, in codon units."""
    start, residues = 0, []
    for ix in range(0, len(sequence) - 2, 3):
        aa = CODONS.get(sequence[ix : ix + 3], "*")
        if aa == "*":
            if len(residues) >= MIN_LENGTH:
                yield start, ix, "".join(residues)
            start, residues = ix + 3, []
        else:
            residues.append(aa)
    if len(residues) >= MIN_LENGTH:
        yield start, start + 3 * len(residues), "".join(residues)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("genome", nargs="+")
    parser.add_argument("--prefix", required=True)
    parser.add_argument("--exclude", nargs="*", default=[])
    args = parser.parse_args()
    out = sys.stdout
    for path in args.genome:
        for name, sequence in read_fasta(path):
            if name in args.exclude:
                continue
            length = len(sequence)
            reverse = sequence.translate(COMPLEMENT)[::-1]
            for strand, seq in (("+", sequence), ("-", reverse)):
                for frame in range(3):
                    for start, end, peptide in segments(seq[frame:]):
                        lo, hi = frame + start, frame + end
                        if strand == "-":
                            lo, hi = length - hi, length - lo
                        out.write(
                            f">{args.prefix}|{name}_{strand}{frame + 1}_{lo + 1}-{hi}\n"
                        )
                        for ix in range(0, len(peptide), 60):
                            out.write(peptide[ix : ix + 60] + "\n")


if __name__ == "__main__":
    main()
