#!/usr/bin/env python3
"""Compare conventional, exhaustive, and site-library PTM searches."""
from __future__ import annotations

import csv
from pathlib import Path

import matplotlib

from _assets import record

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

HERE = Path(__file__).resolve().parent
PAPER = HERE.parent.parent
SRC = HERE / "ptm_comparison_results.csv"
OUT = PAPER / "figures" / "ptm_strategy_comparison.png"


def main() -> int:
    with SRC.open(newline="") as fh:
        rows = list(csv.DictReader(fh))

    labels = ["Conventional", "Exhaustive", "500k-site\nlibrary"]
    colors = ["#94a3b8", "#dc2626", "#7c3aed"]
    panels = (
        ("Whole-search wall time", "Seconds", [float(row["wall_seconds"]) for row in rows], "{:.2f}"),
        ("Peak resident memory", "Peak RSS (GiB)", [float(row["peak_rss_mib"]) / 1024 for row in rows], "{:.2f}"),
        ("Database peptides", "Peptides (millions)", [float(row["database_peptides"]) / 1e6 for row in rows], "{:.2f}"),
        ("Theoretical fragments", "Fragments (millions)", [float(row["database_fragments"]) / 1e6 for row in rows], "{:.1f}"),
    )
    fig, axes = plt.subplots(2, 2, figsize=(7.2, 5.2), dpi=300)
    fig.suptitle("Site-specific PTM search against explicit baselines", fontsize=11, weight="bold")
    for ax, (title, ylabel, values, value_format) in zip(axes.flat, panels, strict=True):
        bars = ax.bar(labels, values, color=colors, width=0.62)
        ax.bar_label(bars, labels=[value_format.format(value) for value in values], padding=3, fontsize=7.1)
        ax.set_title(title, fontsize=9.5)
        ax.set_ylabel(ylabel)
        ax.margins(y=0.20)
        ax.grid(axis="y", alpha=0.2)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

    fig.tight_layout(rect=(0, 0, 1, 0.93), pad=1.0)
    OUT.parent.mkdir(exist_ok=True)
    fig.savefig(OUT, metadata={"Software": None})
    plt.close(fig)
    record(
        "fig.ptm-comparison",
        str(OUT.relative_to(PAPER)),
        kind="figure",
        inputs=[str(SRC.relative_to(PAPER))],
        desc="Conventional, exhaustive, and site-specific PTM search comparison",
    )
    print(f"wrote {OUT.relative_to(PAPER)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
