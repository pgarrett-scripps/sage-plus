#!/usr/bin/env python3
"""Plot localization error against known phosphosites and the margin trade-off."""
from __future__ import annotations

from pathlib import Path

from _assets import record
from _figure_style import MUTED, panel, plt, save_figure
from _localization import LOCALIZATION_INPUTS, load, modification, point

HERE = Path(__file__).resolve().parent
PAPER = HERE.parent.parent
OUT = PAPER / "figures" / "localization.png"

# One ordinal ramp for the separating-ion margin; the default margin is the
# darkest so it reads as the reference series.
RAMP = {0: "#b9d3f2", 1: "#1c5cab", 2: "#5e9be6", 3: "#8db8ee"}
WIDTH = {0: 1.4, 1: 2.4, 2: 1.4, 3: 1.4}
MEASURE = ("#2a78d6", "#eb6834", "#7a8a99")


def main() -> int:
    summary = load()
    curves = summary["known_sites"]
    margins = summary["margins"]
    fig, axes = plt.subplots(1, 3, figsize=(10.4, 3.7), layout="constrained")

    ax = axes[0]
    panel(ax, "A", "Estimated vs true FLR", grid="both")
    top = 0
    for margin in margins:
        rows = [row for row in curves[str(margin)]["points"] if row["passing"]]
        x = [100 * row["estimated_flr"] for row in rows]
        y = [100 * row["true_flr"] for row in rows]
        top = max(top, *x, *y)
        ax.plot(x, y, marker="o", markersize=3.5, color=RAMP[margin], linewidth=WIDTH[margin],
                label=f"Margin {margin}" + (" (default)" if margin == 1 else ""))
    top *= 1.08
    ax.plot([0, top], [0, top], color=MUTED, linestyle=":", linewidth=1)
    ax.set_xlim(0, top)
    ax.set_ylim(0, top)
    ax.set_xlabel("Estimated FLR (%)")
    ax.set_ylabel("True FLR (%)")
    ax.legend(loc="upper left", fontsize=8.5)

    ax = axes[1]
    panel(ax, "B", "Margin trade-off at 1%")
    passing = [point(curves[str(m)])["passing"] for m in margins]
    drawn = ax.bar([str(m) for m in margins], passing, color=[RAMP[m] for m in margins], width=0.62)
    ax.bar_label(drawn, labels=[f"{point(curves[str(m)])['true_flr'] * 100:.2f}%" for m in margins],
                 padding=3, fontsize=9, color=MUTED)
    ax.set_xlabel("Separating-ion margin")
    ax.set_ylabel("Passing localizations")
    ax.margins(y=0.18)

    ax = axes[2]
    panel(ax, "C", "Biological phosphoproteome")
    per_type = summary["pxd007058"]
    rows = (("Phospho", modification(per_type, "Phospho", "79.96")),
            ("Oxidation", modification(per_type, "Oxidation", "15.99")))
    measures = (("localizations_q01", "Localizations"), ("certain", "Certain"),
                ("protein_sites_q01", "Protein sites"))
    width = 0.26
    for index, (key, label) in enumerate(measures):
        values = [row[key] for _, row in rows]
        drawn = ax.bar([p + (index - 1) * width for p in range(len(rows))], values, width,
                       color=MEASURE[index], label=label)
        ax.bar_label(drawn, labels=[f"{value:,}" for value in values], padding=2,
                     fontsize=7.5, color=MUTED, rotation=90)
    ax.set_xticks(range(len(rows)))
    ax.set_xticklabels([name for name, _ in rows])
    ax.set_ylabel("Count at 1%")
    ax.margins(y=0.32)
    ax.legend(loc="upper right", fontsize=8.5)
    for axis in axes:
        axis.tick_params(colors=MUTED)
    fig.get_layout_engine().set(w_pad=0.1)

    save_figure(fig, OUT)
    record("fig.localization", str(OUT.relative_to(PAPER)), kind="figure",
           inputs=LOCALIZATION_INPUTS,
           desc="Known-site localization error, margin trade-off, and phosphoproteome counts.")
    vector = OUT.parent / "vector" / OUT.with_suffix(".svg").name
    record("fig.localization-vector", str(vector.relative_to(PAPER)), kind="figure",
           inputs=LOCALIZATION_INPUTS, desc="Scalable companion to the localization figure.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
