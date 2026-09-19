#!/usr/bin/env python3
"""Plot mass offset search cost, index size, memory, and error behavior."""
from __future__ import annotations

from pathlib import Path

from matplotlib.patches import Patch

from _assets import record
from _figure_style import MUTED, panel, plt, save_figure
from _mass_offset import MASS_OFFSET_INPUTS, load, scale_rows

HERE = Path(__file__).resolve().parent
PAPER = HERE.parent.parent
OUT = PAPER / "figures" / "mass_offset.png"

# One ordinal blue ramp for the number of configured offsets, with the indexed
# baseline in neutral ink. Measures in panel D are two categorical slots.
BASELINE = "#8d99a6"
RAMP = ("#86b6ef", "#3987e5", "#1c5cab")
MEASURE = ("#2a78d6", "#eb6834")
ORDER = ("indexed", "offsets-1", "offsets-2", "offsets-3")
TICKS = ("None\n(indexed)", "One", "Two", "Three")


def bars(ax, values, fmt):
    drawn = ax.bar(TICKS, values, color=(BASELINE,) + RAMP, width=0.62)
    ax.bar_label(drawn, labels=[fmt.format(value) for value in values],
                 padding=3, fontsize=9, color=MUTED)
    ax.margins(y=0.26)
    ax.set_xlabel("Configured offsets")


def main() -> int:
    summary = load()
    rows = scale_rows(summary)

    fig, axes = plt.subplots(2, 2, figsize=(7.1, 5.4), layout="constrained")

    panel(axes[0][0], "A", "Search stage")
    bars(axes[0][0], [rows[key]["search_stage_seconds"] for key in ORDER], "{:.1f}")
    axes[0][0].set_ylabel("Seconds")

    panel(axes[0][1], "B", "Indexed peptides")
    bars(axes[0][1], [rows[key]["database_peptides"] / 1e6 for key in ORDER], "{:.2f}")
    axes[0][1].set_ylabel("Millions of peptides")

    panel(axes[1][0], "C", "Peak resident memory")
    bars(axes[1][0], [rows[key]["peak_rss_mib"] / 1024 for key in ORDER], "{:.2f}")
    axes[1][0].set_ylabel("Peak RSS (GiB)")

    ax = axes[1][1]
    panel(ax, "D", "Error behavior")
    modes = ("indexed", "offset")
    site_error, entrapment = [], []
    for mode in modes:
        correct = sum(row[mode]["correct_site_events"] for row in summary["localization"])
        incorrect = sum(row[mode]["incorrect_site_events"] for row in summary["localization"])
        site_error.append(100 * incorrect / (correct + incorrect))
        threshold = next(r for r in summary["entrapment"][mode] if r["nominal_q"] == 0.01)
        entrapment.append(100 * threshold["combined_fdp"])
    width = 0.34
    positions = range(len(modes))
    for offset, values, color, label in (
        (-width / 2, site_error, MEASURE[0], "Inconsistent sites"),
        (width / 2, entrapment, MEASURE[1], "Entrapment FDP"),
    ):
        drawn = ax.bar([p + offset for p in positions], values, width, color=color, label=label)
        ax.bar_label(drawn, labels=[f"{value:.1f}" for value in values],
                     padding=3, fontsize=9, color=MUTED)
    ax.set_xticks(list(positions))
    ax.set_xticklabels(["Indexed", "Offset"])
    ax.set_xlabel("Search mode")
    ax.set_ylabel("Percent")
    ax.margins(y=0.30)
    ax.legend(handles=[Patch(facecolor=color, label=label) for color, label in
                       zip(MEASURE, ("Inconsistent sites", "Entrapment FDP"), strict=True)],
              loc="upper left", fontsize=9)

    save_figure(fig, OUT)
    record(
        "fig.mass-offset",
        str(OUT.relative_to(PAPER)),
        kind="figure",
        inputs=MASS_OFFSET_INPUTS,
        desc="Mass offset search cost, index size, memory, and error behavior.",
    )
    vector = OUT.parent / "vector" / OUT.with_suffix(".svg").name
    record(
        "fig.mass-offset-vector",
        str(vector.relative_to(PAPER)),
        kind="figure",
        inputs=MASS_OFFSET_INPUTS,
        desc="Scalable companion to the mass offset figure.",
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
