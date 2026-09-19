#!/usr/bin/env python3
"""Retain the original reference and record every excluded ambiguous protein."""

import argparse
from pathlib import Path

from generate_ptm_library_benchmark import fasta_entries
from provenance import atomic_json, file_identities


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("source", type=Path)
    p.add_argument("--output", type=Path, required=True)
    args = p.parse_args()
    if args.output.exists():
        raise ValueError("Preserve reference variants with distinct output paths")
    excluded = []
    retained = 0
    with args.output.open("w") as stream:
        for accession, sequence in fasta_entries(args.source):
            invalid = set(sequence) - set("ACDEFGHIKLMNPQRSTVWYUO")
            if invalid:
                excluded.append({"accession": accession, "length": len(sequence),
                                 "invalid_residues": sorted(invalid)})
            else:
                stream.write(f">{accession}\n{sequence}\n")
                retained += 1
    atomic_json(args.output.with_suffix(".reference.json"), {
        "status": "reference_variant", "retained_proteins": retained,
        "excluded_proteins": excluded,
        "inputs": file_identities([args.source, Path(__file__).resolve()]),
        "outputs": file_identities([args.output]),
        "scope": "Entire proteins containing undefined residue masses are excluded identically for both engines. Known subsequences from those proteins are also lost. No unknown residues are replaced with invented amino acids. The original reference is retained."})
    print(retained, "retained,", len(excluded), "excluded")


if __name__ == "__main__":
    main()
