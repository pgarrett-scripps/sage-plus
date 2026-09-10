#!/usr/bin/env python3
"""Compare canonical Sage result identities and confidence between two runs."""

from __future__ import annotations

import argparse
import json
import math
import shutil
import statistics
import subprocess
from pathlib import Path

from provenance import atomic_json, sha256


def load_rows(path: Path, duckdb: Path) -> list[dict]:
    quoted = str(path.resolve(strict=True)).replace("'", "''")
    query = f"SELECT * EXCLUDE (protein_sites, reporter_ion_intensity) FROM read_parquet('{quoted}')"
    result = subprocess.run([str(duckdb), "-json", "-c", query], check=True, capture_output=True, text=True)
    return json.loads(result.stdout or "[]")


def identity(row: dict) -> tuple:
    return (row["filename"], row["scannr"], row["rank"], row["peptide"],
            row["charge"], row["is_decoy"], row.get("label_channel"))


def indexed(rows: list[dict]) -> dict[tuple, dict]:
    result = {}
    for row in rows:
        key = identity(row)
        if key in result:
            raise ValueError(f"ambiguous duplicate spectrum identity: {key}")
        for column in ("spectrum_q", "peptide_q"):
            value = row[column]
            if value is None or not math.isfinite(value) or not 0 <= value <= 1:
                raise ValueError(f"invalid {column} for {key}")
        result[key] = row
    return result


def set_change(left: set, right: set) -> dict:
    return {
        "baseline": len(left), "candidate": len(right), "shared": len(left & right),
        "added": sorted(right - left, key=repr), "lost": sorted(left - right, key=repr),
        "jaccard": len(left & right) / len(left | right) if left | right else None,
    }


def compare(left_rows: list[dict], right_rows: list[dict], threshold: float) -> dict:
    if not math.isfinite(threshold) or not 0 <= threshold <= 1:
        raise ValueError("q threshold must be finite and between zero and one")
    left, right = indexed(left_rows), indexed(right_rows)
    passing = lambda rows: {key for key, row in rows.items() if not row["is_decoy"] and row["spectrum_q"] <= threshold}
    peptides = lambda rows: {row["peptide"] for row in rows.values() if not row["is_decoy"] and row["peptide_q"] <= threshold}
    shifts = {}
    for column in ("spectrum_q", "peptide_q", "hyperscore", "sage_discriminant_score"):
        differences = [right[key][column] - left[key][column] for key in left.keys() & right.keys()
                       if left[key].get(column) is not None and right[key].get(column) is not None
                       and math.isfinite(left[key][column]) and math.isfinite(right[key][column])]
        shifts[column] = {
            "compared": len(differences), "changed": sum(value != 0 for value in differences),
            "mean_candidate_minus_baseline": statistics.mean(differences) if differences else None,
            "max_absolute_change": max(map(abs, differences), default=None),
        }
    return {
        "q_threshold": threshold,
        "stored_psms": set_change(set(left), set(right)),
        "accepted_target_psms": set_change(passing(left), passing(right)),
        "accepted_peptidoforms": set_change(peptides(left), peptides(right)),
        "confidence_shifts": shifts,
        "scope": "stored result rows, including decoys for score comparisons",
    }


def markdown(report: dict) -> str:
    lines = ["Sage run comparison", "", f"Q-value threshold: {report['q_threshold']}", "",
             "| Population | Baseline | Candidate | Shared | Added | Lost |", "|---|---:|---:|---:|---:|---:|"]
    for name in ("stored_psms", "accepted_target_psms", "accepted_peptidoforms"):
        values = report[name]
        lines.append(f"| {name} | {values['baseline']} | {values['candidate']} | {values['shared']} | {len(values['added'])} | {len(values['lost'])} |")
    lines.extend(["", "Exact identities and confidence shifts are recorded in the JSON report.",
                  "Comparisons cover stored rows. Differences in output cutoffs can change coverage.",
                  "Discriminant differences are descriptive and require matching scoring semantics.", ""])
    return "\n".join(lines)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("baseline", type=Path)
    parser.add_argument("candidate", type=Path)
    parser.add_argument("--q", type=float, default=0.01)
    parser.add_argument("--duckdb", type=Path, default=Path(shutil.which("duckdb") or "duckdb"))
    parser.add_argument("--output", type=Path, required=True, help="Output stem for .json and .md reports")
    args = parser.parse_args()
    paths = [path / "results.sage.parquet" if path.is_dir() else path for path in (args.baseline, args.candidate)]
    report = compare(*(load_rows(path, args.duckdb) for path in paths), args.q)
    report["runs"] = {}
    for name, path in zip(("baseline", "candidate"), paths):
        summary = path.parent / "run-summary.json"
        report["runs"][name] = {
            "path": str(path.resolve()), "sha256": sha256(path),
            "summary": json.loads(summary.read_text()) if summary.exists() else None,
        }
    atomic_json(args.output.with_suffix(".json"), report)
    args.output.with_suffix(".md").write_text(markdown(report))
    print(markdown(report))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
