#!/usr/bin/env python3
"""Independently check reported threshold counts against extracted peptide identities."""

import argparse
import csv
import json
from pathlib import Path

from analyze_fdrbench_validation import completed_summaries, LIMITS
from provenance import atomic_json, file_identities


def rows(path):
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True)
    args = parser.parse_args()
    root = args.input.resolve()
    summaries = completed_summaries(root)
    manifest = json.loads((root / "manifest.json").read_text())
    labels = manifest.get("engine_labels", ["Sage", "Sage Plus"])
    checked, inputs = [], [root / "manifest.json"]
    for summary in summaries:
        directory = root / "seeds" / str(summary["seed"])
        pair = directory / "paired.txt"
        types = {}
        for row in rows(pair):
            key, kind = row["sequence"], row["peptide_type"]
            if key in types and types[key] != kind:
                raise RuntimeError(f"ambiguous peptide type in {pair}")
            types[key] = kind
        inputs.append(pair)
        for engine, slug in zip(labels, ["sage", "sage-plus"]):
            path = directory / slug / "peptides.tsv"
            peptides = rows(path)
            if len({row["peptide"] for row in peptides}) != len(peptides):
                raise RuntimeError(f"duplicate extracted peptide in {path}")
            if any(types.get(row["peptide"]) not in ("target", "p_target") for row in peptides):
                raise RuntimeError(f"unknown peptide identity in {path}")
            inputs.append(path)
            for threshold in LIMITS:
                selected = [row for row in peptides if float(row["q_value"]) <= threshold]
                targets = sum(types[row["peptide"]] == "target" for row in selected)
                entrapments = len(selected) - targets
                point = summary["engines"][engine]["q_points"][str(threshold)]
                if (point is None and selected) or (point is not None and
                        (point["targets"] != targets or point["entrapments"] != entrapments)):
                    raise RuntimeError(f"threshold counts disagree for {engine}, seed {summary['seed']}, q {threshold}")
                checked.append({"seed": summary["seed"], "engine": engine, "q": threshold,
                                "targets": targets, "entrapments": entrapments})
    atomic_json(root / "independent-count-check.json", {
        "status": "passed", "scope": "threshold counts, not an independent FDP estimator validation",
        "inputs": file_identities(inputs), "checks": checked})
    print(f"Verified {len(checked)} threshold count pairs")


if __name__ == "__main__":
    main()
