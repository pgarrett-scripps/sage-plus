#!/usr/bin/env python3
"""Generate deterministic, dense PTM site libraries and search configs."""

from __future__ import annotations

import argparse
import json
import random
from pathlib import Path


REPO = Path(__file__).resolve().parents[1]
DEFAULT_FASTA = REPO / "data/silac-k6r6/human-reviewed.fasta"
DEFAULT_MZML = REPO / "data/silac-k6r6/HEK_SILAC-K6R6.mzML"
DEFAULT_OUTPUT = REPO / "benchmarks/.work/ptm-library"
MAX_COMBINATIONS = 4
MODIFICATIONS = {
    "S": "Phospho",
    "T": "Phospho",
    "Y": "Phospho",
    "M": "Oxidation",
    "K": "Acetyl",
    "N": "Deamidated",
    "Q": "Deamidated",
}


def fasta_entries(path: Path):
    accession = None
    sequence: list[str] = []
    with path.open() as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if accession is not None:
                    yield accession, "".join(sequence)
                accession = line[1:].split()[0]
                sequence = []
            else:
                sequence.append(line)
    if accession is not None:
        yield accession, "".join(sequence)


def sample_sites(fasta: Path, count: int, seed: int) -> tuple[list[tuple[str, int, str, str]], int]:
    rng = random.Random(seed)
    selected: list[tuple[str, int, str, str]] = []
    eligible = 0
    for protein, sequence in fasta_entries(fasta):
        for position, residue in enumerate(sequence, start=1):
            modification = MODIFICATIONS.get(residue)
            if modification is None:
                continue
            site = (protein, position, residue, modification)
            eligible += 1
            if len(selected) < count:
                selected.append(site)
                continue
            replacement = rng.randrange(eligible)
            if replacement < count:
                selected[replacement] = site
    if eligible < count:
        raise RuntimeError(f"requested {count:,} sites but FASTA has {eligible:,} eligible sites")
    rng.shuffle(selected)
    return selected, eligible


def variable_modifications(site_mode: str = "library") -> dict[str, dict[str, object]]:
    return {
        "Phospho": {
            "mass": 79.966331,
            "max_count": 3,
            "site_mode": site_mode,
            "neutral_losses": [97.976896],
            "sites": ["S", "T", "Y"],
        },
        "Oxidation": {
            "mass": 15.994915,
            "max_count": 2,
            "site_mode": site_mode,
            "sites": ["M"],
        },
        "Acetyl": {
            "mass": 42.010565,
            "max_count": 2,
            "site_mode": site_mode,
            "sites": ["K"],
        },
        "Deamidated": {
            "mass": 0.984016,
            "max_count": 2,
            "site_mode": site_mode,
            "sites": ["N", "Q"],
        },
    }


def search_config(fasta: Path, mzml: Path, library: Path) -> dict[str, object]:
    return {
        "database": {
            "bucket_size": 16384,
            "enzyme": {
                "missed_cleavages": 1,
                "cleave_at": "KR",
                "restrict": "P",
                "min_len": 7,
                "max_len": 50,
            },
            "static_mods": {
                "Carbamidomethyl": {"mass": 57.021464, "sites": ["C"]},
            },
            "variable_mods": variable_modifications(),
            "max_variable_mods": 1,
            "max_total_variable_mods": 3,
            "max_combinations": MAX_COMBINATIONS,
            "ptm_library": {"path": str(library), "strict": True},
            "decoy_tag": "rev_",
            "generate_decoys": True,
            "fasta": str(fasta),
        },
        "precursor_tol": {"ppm": [-10.0, 10.0]},
        "fragment_tol": {"ppm": [-20.0, 20.0]},
        "isotope_errors": [-1, 3],
        "deisotope": True,
        "chimera": False,
        "max_fragment_charge": 2,
        "min_matched_peaks": 6,
        "report_psms": 1,
        "output_filter": {"psm_q_value": 0.05},
        "max_memory_gb": 24.0,
        "min_free_memory_gb": 2.0,
        "batch_size": 1,
        "mzml_paths": [str(mzml)],
        "score_type": "SageHyperScore",
    }


def exhaustive_search_config(fasta: Path, mzml: Path) -> dict[str, object]:
    config = search_config(fasta, mzml, Path("unused"))
    database = config["database"]
    assert isinstance(database, dict)
    database["variable_mods"] = variable_modifications(site_mode="exhaustive")
    database["max_variable_mods"] = 3
    database.pop("ptm_library")
    return config


def write_library(path: Path, sites: list[tuple[str, int, str, str]]) -> None:
    with path.open("w") as handle:
        handle.write("protein\tposition\tresidue\tmodification\n")
        for protein, position, residue, modification in sites:
            handle.write(f"{protein}\t{position}\t{residue}\t{modification}\n")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fasta", type=Path, default=DEFAULT_FASTA)
    parser.add_argument("--mzml", type=Path, default=DEFAULT_MZML)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--sizes", type=int, nargs="+", default=[0, 50_000, 200_000, 500_000])
    parser.add_argument("--seed", type=int, default=20260901)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    sizes = sorted(set(args.sizes))
    if not sizes or sizes[0] < 0:
        raise RuntimeError("library sizes must be nonnegative")
    args.output.mkdir(parents=True, exist_ok=True)
    selected, eligible = sample_sites(args.fasta.resolve(), sizes[-1], args.seed)

    manifest = {
        "seed": args.seed,
        "eligible_sites": eligible,
        "sizes": sizes,
        "modifications": sorted(set(MODIFICATIONS.values())),
        "max_combinations": MAX_COMBINATIONS,
        "fasta": str(args.fasta.resolve()),
        "mzml": str(args.mzml.resolve()),
    }
    for size in sizes:
        library = args.output / f"sites-{size}.tsv"
        config = args.output / f"config-{size}.json"
        write_library(library, selected[:size])
        config.write_text(json.dumps(
            search_config(args.fasta.resolve(), args.mzml.resolve(), library.resolve()),
            indent=2,
            sort_keys=True,
        ) + "\n")
        print(f"wrote {library.relative_to(REPO)} and {config.relative_to(REPO)}")
    exhaustive_config = args.output / "config-exhaustive.json"
    exhaustive_config.write_text(json.dumps(
        exhaustive_search_config(args.fasta.resolve(), args.mzml.resolve()),
        indent=2,
        sort_keys=True,
    ) + "\n")
    print(f"wrote {exhaustive_config.relative_to(REPO)}")
    (args.output / "manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    print(f"eligible modification sites in FASTA: {eligible:,}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
