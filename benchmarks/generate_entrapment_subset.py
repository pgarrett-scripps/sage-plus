#!/usr/bin/env python3
"""Create a deterministic protein subset for a bounded entrapment benchmark."""
from __future__ import annotations

import argparse
import random
from pathlib import Path


REPO = Path(__file__).resolve().parents[1]
DEFAULT_INPUT = REPO / "data/silac-k6r6/human-reviewed.fasta"
DEFAULT_OUTPUT = REPO / "benchmarks/.work/entrapment-small/human-subset.fasta"


def fasta_entries(path: Path) -> list[tuple[str, str]]:
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


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--proteins", type=int, default=2_000)
    parser.add_argument("--seed", type=int, default=20260901)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    entries = fasta_entries(args.input)
    if args.proteins < 1 or args.proteins > len(entries):
        raise RuntimeError(
            f"requested {args.proteins:,} proteins from {len(entries):,} entries"
        )
    rng = random.Random(args.seed)
    selected = sorted(rng.sample(range(len(entries)), args.proteins))
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w") as handle:
        for index in selected:
            header, sequence = entries[index]
            handle.write(f"{header}\n")
            for start in range(0, len(sequence), 60):
                handle.write(sequence[start:start + 60] + "\n")
    print(
        f"wrote {args.output.relative_to(REPO)} with "
        f"{len(selected):,} of {len(entries):,} proteins"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
