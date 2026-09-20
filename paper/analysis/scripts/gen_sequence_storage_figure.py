#!/usr/bin/env python3
"""Plot memory savings from protein-backed target peptide sequences."""
from __future__ import annotations

import csv
from pathlib import Path

import matplotlib

from _assets import record

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

HERE = Path(__file__).resolve().parent
PAPER = HERE.parent.parent
SRC = HERE / "sequence_storage_results.csv"
OUT = PAPER / "figures" / "sequence_storage_memory.png"


def main() -> int:
    with SRC.open(newline="") as fh:
        rows = list(csv.DictReader(fh))

    stages = ("Unmodified digestion", "Expanded targets and decoys")
    colors = {"Tryptic": "#2563eb", "Semi-enzymatic": "#0f766e", "Non-specific": "#9333ea"}
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.0), dpi=300, sharey=True)
    for ax, stage in zip(axes, stages, strict=True):
        selected = [row for row in rows if row["stage"] == stage]
        labels = [row["cleavage"] for row in selected]
        values = [float(row["rss_reduction_percent"]) for row in selected]
        bars = ax.barh(labels, values, color=[colors[label] for label in labels])
        ax.set_title(stage)
        ax.set_xlabel("Peak RSS reduction (%)")
        ax.set_xlim(0, 46)
        ax.grid(axis="x", alpha=0.2)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.bar_label(bars, labels=[f"{value:.1f}%" for value in values], padding=3, fontsize=8)
    fig.suptitle("Protein-backed peptide sequences reduce peak memory", fontsize=11, weight="bold")
    fig.tight_layout(pad=0.8)

    OUT.parent.mkdir(exist_ok=True)
    fig.savefig(OUT, metadata={"Software": None})
    plt.close(fig)
    record(
        "fig.sequence-storage",
        str(OUT.relative_to(PAPER)),
        kind="figure",
        inputs=[str(SRC.relative_to(PAPER))],
        desc="Peak-memory reductions from shared protein-backed sequences",
    )
    print(f"wrote {OUT.relative_to(PAPER)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
