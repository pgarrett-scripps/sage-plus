#!/usr/bin/env python3
"""Plot dense PTM site-library scaling results."""
from __future__ import annotations

import csv
from pathlib import Path

import matplotlib

from _assets import record

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

HERE = Path(__file__).resolve().parent
PAPER = HERE.parent.parent
SRC = HERE / "ptm_library_results.csv"
BASELINES = HERE / "ptm_comparison_results.csv"
OUT = PAPER / "figures" / "ptm_library_scaling.png"


def main() -> int:
    with SRC.open(newline="") as fh:
        rows = list(csv.DictReader(fh))
    with BASELINES.open(newline="") as fh:
        comparisons = {row["condition"]: row for row in csv.DictReader(fh)}

    labels = ["0", "50k", "200k", "500k"]
    panels = (
        ("Whole-search wall time", "Seconds", "wall_seconds", 1.0, "{:.2f}"),
        ("Peak resident memory", "Peak RSS (GiB)", "peak_rss_mib", 1024.0, "{:.2f}"),
        ("Database peptides", "Peptides (millions)", "database_peptides", 1e6, "{:.2f}"),
        ("Theoretical fragments", "Fragments (millions)", "database_fragments", 1e6, "{:.1f}"),
    )
    fig, axes = plt.subplots(2, 2, figsize=(7.2, 5.0), dpi=300)
    fig.suptitle("PTM site-library scaling with conventional and exhaustive baselines", fontsize=11, weight="bold")
    for ax, (title, ylabel, key, divisor, value_format) in zip(axes.flat, panels, strict=True):
        values = [float(row[key]) / divisor for row in rows]
        conventional = float(comparisons["Conventional"][key]) / divisor
        exhaustive = float(comparisons["Exhaustive"][key]) / divisor
        ax.plot(labels, values, color="#7c3aed", marker="o", linewidth=2.0, label="Site library")
        ax.axhline(conventional, color="#64748b", linestyle=":", linewidth=1.5, label="Conventional")
        ax.axhline(exhaustive, color="#dc2626", linestyle="--", linewidth=1.5, label="Exhaustive")
        for index, value in enumerate(values):
            ax.annotate(value_format.format(value), (index, value), xytext=(0, 5), textcoords="offset points", ha="center", fontsize=7.0)
        ax.set_title(title, fontsize=9.5)
        ax.set_xlabel("PTM library sites")
        ax.set_ylabel(ylabel)
        ax.margins(y=0.14)
        ax.grid(axis="y", alpha=0.2)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

    handles, legend_labels = axes.flat[0].get_legend_handles_labels()
    fig.legend(handles, legend_labels, loc="upper center", bbox_to_anchor=(0.5, 0.92), ncol=3, frameon=False, fontsize=8)
    fig.tight_layout(rect=(0, 0, 1, 0.84), pad=1.0)
    OUT.parent.mkdir(exist_ok=True)
    fig.savefig(OUT, metadata={"Software": None})
    plt.close(fig)
    record(
        "fig.ptm-library",
        str(OUT.relative_to(PAPER)),
        kind="figure",
        inputs=[str(SRC.relative_to(PAPER)), str(BASELINES.relative_to(PAPER))],
        desc="Dense PTM site-library scaling with conventional and exhaustive baselines",
    )
    print(f"wrote {OUT.relative_to(PAPER)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
