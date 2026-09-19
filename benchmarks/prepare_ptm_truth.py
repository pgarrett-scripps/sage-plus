#!/usr/bin/env python3
"""Recover synthesis-defined peptides and sites from the published design."""

import argparse
import csv
import json
import random
from collections import defaultdict
from pathlib import Path

import openpyxl

from generate_ptm_library_benchmark import fasta_entries
from provenance import atomic_json, file_identities


def recover_truth(fasta, design_rows):
    source = dict(fasta_entries(fasta))
    truth = defaultdict(set)
    library_for = defaultdict(set)
    checks = []
    for library, seed, position, length, _, expected in design_rows:
        library, position, length = int(library), int(position), int(length)
        header = f"IPI:{chr(65 + (library - 1) // 12)}{(library - 1) % 12 + 1}"
        sequence = source[header]
        if len(sequence) % length:
            raise ValueError(f"Library sequence is not divisible by seed length: {header}")
        peptides = [sequence[i:i + length] for i in range(0, len(sequence), length)]
        if len(peptides) * 2 != int(expected):
            raise ValueError(f"Expected modified plus unmodified count disagrees: {header}")
        seed_plain = seed.replace("p", "")
        if len(seed_plain) != length or seed_plain not in peptides:
            raise ValueError(f"Seed does not match published FASTA library: {header}")
        if seed[position - 1:position + 1] != "p" + seed_plain[position - 1]:
            raise ValueError(f"Phosphosite numbering disagrees: {header}")
        for peptide in peptides:
            if peptide[position - 1] not in "STY":
                raise ValueError(f"Nonphosphorylatable truth position: {header}")
            truth[peptide].add(position)
            library_for[peptide].add(library)
        checks.append({"library": library, "header": header, "seed": seed,
                       "position": position, "variants": len(peptides)})
    return truth, library_for, checks


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--root", type=Path, required=True)
    args = p.parse_args()
    inputs = args.root / "inputs/PXD000138"
    fasta = inputs / "009606_IPI_v3_72_with_P_Lib_062111.fasta"
    design = inputs / "synthesis-design.xlsx"
    book = openpyxl.load_workbook(design, read_only=True, data_only=True)
    rows = list(book["Seed_Peptides"].values)[2:]
    truth, libraries, checks = recover_truth(fasta, rows)
    output = args.root / "ptm-truth"
    output.mkdir(exist_ok=True)
    entries = []
    with (output / "synthetic.fasta").open("w") as f, (output / "truth.tsv").open("w") as t:
        writer = csv.writer(t, delimiter="\t")
        writer.writerow(["protein", "sequence", "positions", "libraries", "unambiguous_site"])
        for index, (sequence, positions) in enumerate(sorted(truth.items())):
            accession = f"synth_{index:06d}"
            f.write(f">{accession}\n{sequence}\n")
            writer.writerow([accession, sequence, ",".join(map(str, sorted(positions))),
                             ",".join(map(str, sorted(libraries[sequence]))), len(positions) == 1])
            entries.extend((accession, position, sequence[position - 1], "Phospho") for position in positions)
    rng = random.Random(20260914)
    rng.shuffle(entries)
    for fraction in (0, 25, 50, 100):
        with (output / f"oracle-sites-{fraction}.tsv").open("w") as f:
            writer = csv.writer(f, delimiter="\t")
            writer.writerow(["protein", "position", "residue", "modification"])
            writer.writerows(entries[:len(entries) * fraction // 100])
    atomic_json(output / "manifest.json", {
        "status": "synthesis_design_recovered", "inputs": file_identities([fasta, design, Path(__file__).resolve()]),
        "openpyxl_version": openpyxl.__version__, "libraries": checks,
        "unique_synthetic_sequences": len(truth),
        "ambiguous_sequences": sum(len(p) != 1 for p in truth.values()),
        "truth_scope": "Known synthesis sites across all 96 libraries. File-to-library assignment is not established. Site accuracy excludes sequences with ambiguous synthesis positions.",
        "search_space": "Complete synthesized sequence library without the IPI background. Restricted synthetic benchmark, not a human proteome production search.",
        "oracle_warning": "Site subsets are truth-derived controls for coverage sensitivity. They cannot establish the benefit of a previously known biological annotation library.",
        "outputs": file_identities([p for p in output.iterdir() if p.suffix in (".fasta", ".tsv")])})
    print(json.dumps({"libraries": len(checks), "sequences": len(truth), "sites": len(entries)}))


if __name__ == "__main__":
    main()
