#!/usr/bin/env python3
"""Plot the exact-prefilter resource and database-size tradeoff."""
from __future__ import annotations

import csv
from pathlib import Path

import matplotlib

from _assets import record

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

HERE = Path(__file__).resolve().parent
PAPER = HERE.parent.parent
SRC = HERE / "benchmark_results.csv"
OUT = PAPER / "figures" / "exact_prefilter.png"


def main() -> int:
    with SRC.open(newline="") as fh:
        row = next(row for row in csv.DictReader(fh) if row["suite"] == "exact prefilter")

    comparisons = (
        ("Wall time", "Wall time (s)", [float(row["baseline_wall_seconds"]), float(row["candidate_wall_seconds"])], "{:.2f}"),
        ("Peak resident memory", "Peak RSS (MiB)", [float(row["baseline_rss_mib"]), float(row["candidate_rss_mib"])], "{:,.0f}"),
        ("Final search database", "Database peptides (millions)", [float(row["baseline_database_peptides"]) / 1e6, float(row["candidate_database_peptides"]) / 1e6], "{:.2f}"),
    )
    fig, axes = plt.subplots(1, 3, figsize=(7.2, 2.8), dpi=300)
    fig.suptitle("Exact prefilter trades time for a smaller search database", fontsize=11, weight="bold")
    for ax, (title, ylabel, values, value_format) in zip(axes, comparisons, strict=True):
        bars = ax.bar(["Off", "On"], values, color=["#94a3b8", "#0f766e"], width=0.60)
        ax.bar_label(bars, labels=[value_format.format(value) for value in values], padding=3, fontsize=8)
        ax.set_title(title, fontsize=9.2)
        ax.set_ylabel(ylabel)
        ax.margins(y=0.22)
        ax.grid(axis="y", alpha=0.2)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

    fig.tight_layout(rect=(0, 0, 1, 0.88), pad=0.9)
    OUT.parent.mkdir(exist_ok=True)
    fig.savefig(OUT, metadata={"Software": None})
    plt.close(fig)
    record(
        "fig.prefilter",
        str(OUT.relative_to(PAPER)),
        kind="figure",
        inputs=[str(SRC.relative_to(PAPER))],
        desc="Exact-prefilter wall time, memory, and database size",
    )
    print(f"wrote {OUT.relative_to(PAPER)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
