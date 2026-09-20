#!/usr/bin/env python3
"""Plot exploratory FDRBench calibration by modification stratum."""

from __future__ import annotations

import csv
from pathlib import Path

import matplotlib

from _assets import record

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402


HERE = Path(__file__).resolve().parent
PAPER = HERE.parent.parent
SRC = HERE / "fdrbench_validation_results.csv"
OUT = PAPER / "figures" / "fdrbench_modification_strata.png"


def main() -> int:
    with SRC.open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    engines = ("Sage", "Sage Plus")
    colors = {"Sage": "#64748b", "Sage Plus": "#2563eb"}
    panels = (
        ("modified_unmodified", "Unmodified best evidence", 0.24),
        ("modified_variable", "Variable-modified best evidence", 0.95),
    )
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.15), dpi=300)
    fig.suptitle("Exploratory modification-stratified calibration", fontsize=11, weight="bold")
    for ax, (analysis, title, ceiling) in zip(axes, panels, strict=True):
        ax.plot([0, 0.20], [0, 0.20], color="#0f172a", linewidth=1, linestyle=":", label="Ideal")
        for engine in engines:
            selected = sorted(
                (
                    row
                    for row in rows
                    if row["analysis"] == analysis and row["engine"] == engine
                ),
                key=lambda row: float(row["threshold"]),
            )
            ax.plot(
                [float(row["threshold"]) for row in selected],
                [float(row["mean"]) for row in selected],
                marker="o",
                markersize=3.5,
                linewidth=1.6,
                color=colors[engine],
                label=engine,
            )
        ax.set_title(title, fontsize=9.5)
        ax.set_xlabel("Reported peptide q-value")
        ax.set_ylabel("Paired FDP")
        ax.set_xlim(0, 0.205)
        ax.set_ylim(0, ceiling)
        ax.grid(alpha=0.2)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
    axes[0].legend(frameon=False, fontsize=7.5, loc="upper left")
    fig.tight_layout(rect=(0, 0, 1, 0.91), pad=1.0)
    OUT.parent.mkdir(exist_ok=True)
    fig.savefig(OUT, metadata={"Software": None})
    plt.close(fig)
    record(
        "fig.fdrbench-modification-strata",
        str(OUT.relative_to(PAPER)),
        kind="figure",
        inputs=[str(SRC.relative_to(PAPER))],
        desc="Exploratory paired FDP calibration by modification stratum",
    )
    print(f"wrote {OUT.relative_to(PAPER)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
