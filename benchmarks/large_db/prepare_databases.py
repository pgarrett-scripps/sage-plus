"""Build the large-database search references.

Scaling series: the human reference plus nested random subsets of the
Integrated Gene Catalog (IGC) of the human gut microbiome, sized as multiples
of the human residue count. A protein enters every subset whose fraction
exceeds a hash of its name, so each subset contains the smaller ones.

Catalog entries lose one terminal stop. Entries that still contain a stop or a
residue outside the twenty standard amino acids are dropped in their entirety,
because Sage Plus rejects undefined residues. Counts are written to a JSON
receipt next to the outputs.
"""

import argparse
import hashlib
import json
from pathlib import Path

STANDARD = set("ACDEFGHIKLMNPQRSTVWY")


def records(path):
    name, chunks = None, []
    with open(path) as handle:
        for line in handle:
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(chunks)
                name, chunks = line[1:].rstrip("\n"), []
            else:
                chunks.append(line.strip())
    if name is not None:
        yield name, "".join(chunks)


def clean(sequence):
    sequence = sequence.upper()
    if sequence.endswith("*"):
        sequence = sequence[:-1]
    return sequence if sequence and set(sequence) <= STANDARD else None


def key(name):
    digest = hashlib.sha256(name.split()[0].encode()).digest()
    return int.from_bytes(digest[:8], "big") / 2**64


def write(handle, name, sequence):
    handle.write(f">{name}\n")
    for ix in range(0, len(sequence), 60):
        handle.write(sequence[ix : ix + 60] + "\n")


def residues(path):
    return sum(len(seq) for _, seq in records(path))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--human", type=Path, required=True)
    parser.add_argument("--catalog", type=Path, required=True)
    parser.add_argument("--catalog-residues", type=int, required=True,
                        help="cleaned catalog residues, from a counting pass")
    parser.add_argument("--multiples", type=float, nargs="+", required=True,
                        help="a multiple at or above the whole catalog writes human-igc-full")
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)

    human_text = args.human.read_text()
    human_residues = residues(args.human)
    fractions = {m: m * human_residues / args.catalog_residues for m in args.multiples}
    handles = {}
    for multiple in args.multiples:
        label = "full" if fractions[multiple] >= 1 else f"{multiple:g}x"
        path = args.out / f"human-igc-{label}.fasta"
        handles[multiple] = path.open("w")
        handles[multiple].write(human_text)
    stats = {m: {"proteins": 0, "residues": 0} for m in args.multiples}
    dropped = kept = 0
    for name, sequence in records(args.catalog):
        sequence = clean(sequence)
        if sequence is None:
            dropped += 1
            continue
        kept += 1
        k = key(name)
        for multiple, fraction in fractions.items():
            if k < fraction:
                write(handles[multiple], name, sequence)
                stats[multiple]["proteins"] += 1
                stats[multiple]["residues"] += len(sequence)
    for handle in handles.values():
        handle.close()
    receipt = {
        "human": str(args.human),
        "human_residues": human_residues,
        "catalog": str(args.catalog),
        "catalog_proteins_kept": kept,
        "catalog_proteins_dropped": dropped,
        "subsets": {f"{m:g}": {"fraction": fractions[m], **stats[m]} for m in args.multiples},
    }
    (args.out / "scaling-databases.json").write_text(json.dumps(receipt, indent=1))


if __name__ == "__main__":
    main()
