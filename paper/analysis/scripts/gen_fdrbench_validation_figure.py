#!/usr/bin/env python3
"""Plot repeated shuffled entrapment calibration and matched-FDP yield."""

from __future__ import annotations

import csv
from pathlib import Path

import matplotlib

from _assets import record

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402


HERE = Path(__file__).resolve().parent
PAPER = HERE.parent.parent
SRC = HERE / "fdrbench_validation_results.csv"
OUT = PAPER / "figures" / "fdrbench_validation.png"


def main() -> int:
    with SRC.open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    engines = ("Sage", "Sage Plus")
    colors = {"Sage": "#64748b", "Sage Plus": "#2563eb"}
    labels = {"Sage": "Sage", "Sage Plus": "Sage Plus"}
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.35), dpi=300)
    fig.suptitle(
        "Repeated paired entrapment comparison",
        fontsize=11,
        weight="bold",
    )

    ax = axes[0]
    ax.plot([0, 0.20], [0, 0.20], color="#0f172a", linewidth=1, linestyle=":", label="Ideal")
    for engine in engines:
        selected = sorted(
            (
                row
                for row in rows
                if row["analysis"] == "shuffled_calibration" and row["engine"] == engine
            ),
            key=lambda row: float(row["threshold"]),
        )
        x = np.array([float(row["threshold"]) for row in selected])
        mean = np.array([float(row["mean"]) for row in selected])
        lower = np.array([float(row["lower"]) for row in selected])
        upper = np.array([float(row["upper"]) for row in selected])
        ax.fill_between(x, lower, upper, color=colors[engine], alpha=0.13)
        ax.plot(x, mean, marker="o", markersize=3.5, linewidth=1.6, color=colors[engine], label=labels[engine])
    ax.set_title("Reported q-value calibration", fontsize=9.5)
    ax.set_xlabel("Reported peptide q-value")
    ax.set_ylabel("Mean paired FDP")
    ax.set_xlim(0, 0.205)
    ax.set_ylim(0, 0.39)
    ax.grid(alpha=0.2)
    ax.legend(frameon=False, fontsize=7.5, loc="upper left")

    ax = axes[1]
    limits = [0.01, 0.05, 0.10, 0.20]
    x = np.arange(len(limits))
    width = 0.34
    for index, engine in enumerate(engines):
        selected = {
            float(row["threshold"]): row
            for row in rows
            if row["analysis"] == "matched_yield" and row["engine"] == engine
        }
        means = np.array([float(selected[limit]["mean"]) for limit in limits])
        lower = np.array([float(selected[limit]["lower"]) for limit in limits])
        upper = np.array([float(selected[limit]["upper"]) for limit in limits])
        position = x + (index - 0.5) * width
        ax.bar(position, means, width, color=colors[engine], label=labels[engine])
        ax.errorbar(
            position,
            means,
            yerr=np.vstack([means - lower, upper - means]),
            fmt="none",
            ecolor="#0f172a",
            elinewidth=0.8,
            capsize=2,
        )
        for bar_x, value in zip(position, means, strict=True):
            ax.text(
                bar_x,
                value - 6,
                f"{value:.0f}",
                ha="center",
                va="top",
                color="white",
                fontsize=6.8,
                weight="bold",
            )
    ax.set_title("Targets at matched paired FDP", fontsize=9.5)
    ax.set_xlabel("FDRBench paired FDP limit")
    ax.set_ylabel("Mean target peptides")
    ax.set_xticks(x, [f"{limit:.0%}" for limit in limits])
    ax.margins(y=0.20)
    ax.grid(axis="y", alpha=0.2)
    ax.legend(frameon=False, fontsize=7.5, loc="upper left")

    for ax in axes:
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
    fig.tight_layout(rect=(0, 0, 1, 0.92), pad=1.0)
    OUT.parent.mkdir(exist_ok=True)
    fig.savefig(OUT, metadata={"Software": None})
    plt.close(fig)
    record(
        "fig.fdrbench-validation",
        str(OUT.relative_to(PAPER)),
        kind="figure",
        inputs=[str(SRC.relative_to(PAPER))],
        desc="Repeated shuffled entrapment calibration and matched-FDP yield",
    )
    print(f"wrote {OUT.relative_to(PAPER)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
