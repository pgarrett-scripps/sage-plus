"""Read the complete-threshold entrapment yield analysis."""
import json
from pathlib import Path

PAPER = Path(__file__).resolve().parents[2]
MATCHED_INPUTS = ["analysis/data/matched-fdp.json"]


def load():
    return json.loads((PAPER / MATCHED_INPUTS[0]).read_text())


def add_matched_stats(stats):
    data = load()
    for row in data["summaries"]:
        stats.add(f"matched.{row['study']}.change", row["mean_paired_percent_change"],
                  fmt="+.2f", desc="Mean paired target peptide yield change at the common estimated FDP ceiling in percent",
                  between=(-100, 100))
    points = [row["selected"] for row in data["runs"]]
    for bound, value in (("min", min(p["paired_fdp"] for p in points)),
                         ("max", max(p["paired_fdp"] for p in points))):
        stats.add(f"matched.fdp.{bound}", 100 * value, fmt=".3f",
                  desc="Achieved conservative paired FDP bound across selected thresholds in percent",
                  between=(0, 1))

