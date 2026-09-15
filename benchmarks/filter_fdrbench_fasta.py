#!/usr/bin/env python3
"""Filter unsupported peptide residues from FDRBench FASTA and pair files."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


ALLOWED = frozenset("ACDEFGHIKLMNPQRSTVWYU")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fasta", type=Path, required=True)
    parser.add_argument("--pairs", type=Path, required=True)
    parser.add_argument("--output-fasta", type=Path, required=True)
    parser.add_argument("--output-pairs", type=Path, required=True)
    return parser.parse_args()


def read_fasta(path: Path) -> list[tuple[str, str]]:
    entries: list[tuple[str, str]] = []
    header = ""
    sequence: list[str] = []
    with path.open() as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if line.startswith(">"):
                if header:
                    entries.append((header, "".join(sequence)))
                header = line
                sequence = []
            elif line:
                sequence.append(line)
    if header:
        entries.append((header, "".join(sequence)))
    return entries


def main() -> int:
    args = parse_args()
    entries = read_fasta(args.fasta)
    kept = [(header, sequence) for header, sequence in entries if set(sequence) <= ALLOWED]
    args.output_fasta.parent.mkdir(parents=True, exist_ok=True)
    with args.output_fasta.open("w") as handle:
        for header, sequence in kept:
            handle.write(f"{header}\n{sequence}\n")
    with args.pairs.open(newline="") as source:
        reader = csv.DictReader(source, delimiter="\t")
        fields = reader.fieldnames
        if fields is None or "sequence" not in fields:
            raise RuntimeError("pair file does not contain a sequence column")
        rows = [row for row in reader if set(row["sequence"]) <= ALLOWED]
    with args.output_pairs.open("w", newline="") as target:
        writer = csv.DictWriter(target, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    if len(kept) != len(rows):
        raise RuntimeError(
            f"filtered FASTA has {len(kept):,} entries but table has {len(rows):,} rows"
        )
    print(
        f"kept {len(kept):,} of {len(entries):,} entries, "
        f"removed {len(entries) - len(kept):,}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
