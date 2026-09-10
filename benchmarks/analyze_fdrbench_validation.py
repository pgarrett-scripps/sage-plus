#!/usr/bin/env python3
"""Aggregate shuffled, modified, and foreign FDRBench validation outputs."""

from __future__ import annotations

import argparse
import csv
import math
import json
import statistics
from pathlib import Path

from provenance import atomic_json, sha256


REPO = Path(__file__).resolve().parents[1]
ROOT = REPO / "benchmarks/results/fdrbench-validation"
PAPER_DATA = REPO / "paper/analysis/scripts/fdrbench_validation_results.csv"
LIMITS = (0.01, 0.02, 0.05, 0.10, 0.20)
ENGINES = ("Sage", "Sage Plus")


def percentile(values: list[float], probability: float) -> float:
    ordered = sorted(values)
    if not ordered:
        raise RuntimeError("cannot take a percentile of an empty list")
    position = (len(ordered) - 1) * probability
    lower = int(position)
    upper = min(lower + 1, len(ordered) - 1)
    fraction = position - lower
    return ordered[lower] * (1 - fraction) + ordered[upper] * fraction


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def point_at_q(rows: list[dict[str, str]], limit: float) -> dict[str, str] | None:
    eligible = [row for row in rows if float(row["q_value"]) <= limit]
    if not eligible:
        return None
    return max(
        eligible,
        key=lambda row: (float(row["q_value"]), -float(row["score"])),
    )


def timing(path: Path) -> tuple[float, float]:
    wall, rss, status = path.read_text().strip().split("\t")
    if status != "0":
        raise RuntimeError(f"nonzero search status recorded in {path}")
    return float(wall), int(rss) / 1024


def completed_summaries(root: Path) -> list[dict]:
    manifest = json.loads((root / "manifest.json").read_text())
    if manifest.get("schema_version") != 2:
        raise RuntimeError("validated schema-2 experiment manifest required, rerun the repaired harness")
    seeds = manifest["seeds"]
    if not seeds or len(set(seeds)) != len(seeds) or manifest["seed_count"] != len(seeds):
        raise RuntimeError("experiment seed matrix is empty, duplicated, or incomplete")
    summaries = []
    engines = manifest.get("engine_labels", ENGINES)
    if len(engines) != 2 or len(set(engines)) != 2:
        raise RuntimeError("two distinct engine labels are required")
    for seed in seeds:
        path = root / "seeds" / str(int(seed)) / "summary.json"
        if manifest["summary_hashes"].get(str(path.resolve())) != sha256(path):
            raise RuntimeError(f"summary checksum mismatch: {path}")
        summary = json.loads(path.read_text())
        if summary.get("status") != "completed" or summary["seed"] != seed or set(summary["engines"]) != set(engines):
            raise RuntimeError(f"incomplete engine matrix: {path}")
        for values in summary["engines"].values():
            for kind in ("q_points", "matched_yield"):
                points = values[kind]
                if set(points) != {str(limit) for limit in LIMITS}:
                    raise RuntimeError(f"incomplete threshold matrix: {path}")
                for point in points.values():
                    if point is not None and (
                        not all(math.isfinite(float(v)) for v in point.values())
                        or point["targets"] < 0 or point["entrapments"] < 0
                        or point["paired_fdp"] < 0
                    ):
                        raise RuntimeError(f"invalid threshold values: {path}")
        summaries.append(summary)
    return summaries


def main() -> int:
    global ROOT, PAPER_DATA, ENGINES
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=REPO / "benchmarks/results/fdrbench-hardening")
    parser.add_argument("--paper-output", type=Path)
    parser.add_argument("--include-exploratory", action="store_true")
    args = parser.parse_args()
    ROOT = args.input.resolve()
    PAPER_DATA = args.paper_output
    summaries = completed_summaries(ROOT)
    ENGINES = tuple(json.loads((ROOT / "manifest.json").read_text()).get("engine_labels", ENGINES))
    if args.include_exploratory and ENGINES != ("Sage", "Sage Plus"):
        raise RuntimeError("historical exploratory directories require the original engine labels")
    seeds = [summary["seed"] for summary in summaries]
    threshold_lookup = {}
    matched_lookup = {}
    for summary in summaries:
        for engine, values in summary["engines"].items():
            for limit in LIMITS:
                threshold_lookup[(summary["seed"], engine, limit)] = values["q_points"][str(limit)]
                matched_lookup[(summary["seed"], engine, limit)] = values["matched_yield"][str(limit)]
    output_rows: list[dict[str, object]] = []
    for engine in ENGINES:
        for limit in LIMITS:
            fdp_values = []
            target_values = []
            entrapment_values = []
            observed = 0
            for seed in seeds:
                row = threshold_lookup.get((seed, engine, limit))
                if row is None:
                    fdp_values.append(0.0)
                    target_values.append(0.0)
                    entrapment_values.append(0.0)
                else:
                    observed += 1
                    fdp_values.append(float(row["paired_fdp"]))
                    target_values.append(float(row["targets"]))
                    entrapment_values.append(float(row["entrapments"]))
            output_rows.append(
                {
                    "analysis": "shuffled_calibration",
                    "engine": engine,
                    "threshold": limit,
                    "mean": statistics.mean(fdp_values),
                    "lower": percentile(fdp_values, 0.025),
                    "upper": percentile(fdp_values, 0.975),
                    "mean_targets": statistics.mean(target_values),
                    "mean_entrapments": statistics.mean(entrapment_values),
                    "replicates": len(seeds),
                    "replicates_with_discoveries": observed,
                }
            )
            yields = []
            for seed in seeds:
                row = matched_lookup.get((seed, engine, limit))
                yields.append(float(row["targets"]) if row is not None else 0.0)
            output_rows.append(
                {
                    "analysis": "matched_yield",
                    "engine": engine,
                    "threshold": limit,
                    "mean": statistics.mean(yields),
                    "lower": percentile(yields, 0.025),
                    "upper": percentile(yields, 0.975),
                    "mean_targets": statistics.mean(yields),
                    "mean_entrapments": "",
                    "replicates": len(seeds),
                    "replicates_with_discoveries": sum(value > 0 for value in yields),
                }
            )

    if args.include_exploratory:
        foreign_root = ROOT / "foreign"
        for engine in ENGINES:
            slug = "sage-plus" if engine == "Sage Plus" else "sage"
            rows = read_rows(foreign_root / slug / "fdp.csv")
            for limit in LIMITS:
                row = point_at_q(rows, limit)
                if row is None:
                    continue
                for name, column in (
                    ("foreign_combined", "combined_fdp"),
                    ("foreign_lower_bound", "lower_bound_fdp"),
                ):
                    output_rows.append(
                        {
                            "analysis": name,
                            "engine": engine,
                            "threshold": limit,
                            "mean": float(row[column]),
                            "lower": float(row[column]),
                            "upper": float(row[column]),
                            "mean_targets": int(row["n_t"]),
                            "mean_entrapments": int(row["n_p"]),
                            "replicates": 1,
                            "replicates_with_discoveries": 1,
                        }
                    )

        modified_root = ROOT / "modified/seed-20260902"
        for engine in ENGINES:
            slug = "sage-plus" if engine == "Sage Plus" else "sage"
            for analysis, filename in (
                ("modified_overall", "fdp.csv"),
                ("modified_variable", "fdp-variable.csv"),
                ("modified_unmodified", "fdp-unmodified.csv"),
            ):
                rows = read_rows(modified_root / slug / filename)
                for limit in LIMITS:
                    row = point_at_q(rows, limit)
                    if row is None:
                        continue
                    output_rows.append(
                        {
                            "analysis": analysis,
                            "engine": engine,
                            "threshold": limit,
                            "mean": float(row["paired_fdp"]),
                            "lower": float(row["paired_fdp"]),
                            "upper": float(row["paired_fdp"]),
                            "mean_targets": int(row["n_t"]),
                            "mean_entrapments": int(row["n_p"]),
                            "replicates": 1,
                            "replicates_with_discoveries": 1,
                        }
                    )

    fields = [
        "analysis",
        "engine",
        "threshold",
        "mean",
        "lower",
        "upper",
        "mean_targets",
        "mean_entrapments",
        "replicates",
        "replicates_with_discoveries",
    ]
    for path in [ROOT / "analysis-summary.csv"] + ([PAPER_DATA] if PAPER_DATA else []):
        path.parent.mkdir(parents=True, exist_ok=True)
        with path.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fields)
            writer.writeheader()
            writer.writerows(output_rows)

    resources = []
    for seed in seeds:
        summary = json.loads((ROOT / "seeds" / str(seed) / "summary.json").read_text())
        for engine in ENGINES:
            values = summary["engines"][engine]["timing"]
            resources.append(
                {
                    "analysis": "shuffled",
                    "replicate": seed,
                    "engine": engine,
                    "wall_seconds": values["wall_seconds"],
                    "peak_rss_mib": values["peak_rss_mib"],
                }
            )
    if args.include_exploratory:
        for analysis, directory in (
            ("modified", ROOT / "modified/seed-20260902"),
            ("foreign", ROOT / "foreign"),
        ):
            for engine in ENGINES:
                slug = "sage-plus" if engine == "Sage Plus" else "sage"
                wall, rss = timing(directory / slug / "timing.tsv")
                resources.append(
                    {
                        "analysis": analysis,
                        "replicate": 20260902 if analysis == "modified" else 1,
                        "engine": engine,
                        "wall_seconds": wall,
                        "peak_rss_mib": rss,
                    }
                )
    with (ROOT / "resource-summary.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(resources[0]))
        writer.writeheader()
        writer.writerows(resources)

    atomic_json(ROOT / "analysis-method.json", {
        "source_manifest_sha256": sha256(ROOT / "manifest.json"),
        "engine_labels": list(ENGINES),
        "zero_discovery_convention": "completed empty discovery sets contribute zero to the mean, missing runs are rejected",
        "interval": "2.5th and 97.5th empirical percentiles across entrapment seeds, not a confidence interval for the mean",
        "replication_scope": "entrapment construction variability conditional on the supplied spectrum file",
    })

    print(f"wrote {ROOT / 'analysis-summary.csv'}")
    if PAPER_DATA:
        print(f"wrote {PAPER_DATA}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
