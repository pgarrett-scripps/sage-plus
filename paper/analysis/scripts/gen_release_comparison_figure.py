#!/usr/bin/env python3
"""Plot matched search workloads for upstream Sage and Sage Plus."""
from __future__ import annotations

import csv
from pathlib import Path
from statistics import median

import matplotlib

from _assets import record

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

HERE = Path(__file__).resolve().parent
PAPER = HERE.parent.parent
SRC = HERE / "sage_vs_plus_results.csv"
OUT = PAPER / "figures" / "release_comparison.png"


def main() -> int:
    with SRC.open(newline="") as fh:
        rows = list(csv.DictReader(fh))

    workloads = ["conventional", "common_mods", "broad_ptm"]
    labels = ["Conventional", "Common\nmodifications", "Broad PTMs"]
    releases = ["Sage v0.15.0-beta.2", "Sage Plus beta.2"]
    colors = ["#64748b", "#2563eb"]
    grouped = {
        (workload, release): [
            row for row in rows
            if row["workload"] == workload and row["release"] == release
        ]
        for workload in workloads
        for release in releases
    }

    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.15), dpi=300)
    fig.suptitle("Sage Plus beta.2 compared with upstream Sage", fontsize=11, weight="bold")
    comparisons = (
        (axes[0], "Wall time", "Wall time (s)", "wall_seconds", "{:.2f}"),
        (axes[1], "Peak resident memory", "Peak RSS (MiB)", "peak_rss_mib", "{:,.0f}"),
    )
    x = np.arange(len(workloads))
    width = 0.34
    for ax, title, ylabel, key, value_format in comparisons:
        for release_index, release in enumerate(releases):
            values = [
                median(float(row[key]) for row in grouped[workload, release])
                for workload in workloads
            ]
            positions = x + (release_index - 0.5) * width
            bars = ax.bar(
                positions,
                values,
                width=width,
                color=colors[release_index],
                label="Sage" if release_index == 0 else "Sage Plus",
            )
            value_labels = []
            for workload_index, value in enumerate(values):
                if release_index == 0:
                    value_labels.append(value_format.format(value))
                else:
                    baseline = median(
                        float(row[key])
                        for row in grouped[workloads[workload_index], releases[0]]
                    )
                    delta = (value / baseline - 1) * 100
                    value_labels.append(f"{value_format.format(value)}\n{delta:+.1f}%")
            ax.bar_label(bars, labels=value_labels, padding=3, fontsize=6.8)
            for workload_index, workload in enumerate(workloads):
                trials = [float(row[key]) for row in grouped[workload, release]]
                offsets = np.linspace(-0.055, 0.055, len(trials))
                ax.scatter(
                    positions[workload_index] + offsets,
                    trials,
                    s=8,
                    color="#0f172a",
                    alpha=0.6,
                    zorder=3,
                )
        ax.set_title(title, fontsize=9.5)
        ax.set_ylabel(ylabel)
        ax.set_xticks(x, labels)
        ax.margins(y=0.22)
        ax.grid(axis="y", alpha=0.2)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

    handles, legend_labels = axes[0].get_legend_handles_labels()
    fig.legend(
        handles,
        legend_labels,
        loc="upper center",
        bbox_to_anchor=(0.5, 0.90),
        ncol=2,
        frameon=False,
        fontsize=8,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.82), pad=1.0)
    OUT.parent.mkdir(exist_ok=True)
    fig.savefig(OUT, metadata={"Software": None})
    plt.close(fig)
    record(
        "fig.release-comparison",
        str(OUT.relative_to(PAPER)),
        kind="figure",
        inputs=[str(SRC.relative_to(PAPER))],
        desc="Matched workload comparison between upstream Sage and Sage Plus",
    )
    print(f"wrote {OUT.relative_to(PAPER)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
