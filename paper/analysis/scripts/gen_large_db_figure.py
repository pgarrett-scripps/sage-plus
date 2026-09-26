#!/usr/bin/env python3
"""Plot large-database search cost, prefilter retention, and identifications."""
from __future__ import annotations

from pathlib import Path

from matplotlib.lines import Line2D

from _assets import record
from _figure_style import GRID, MUTED, panel, plt, save_figure
from _large_db import LARGE_DB_INPUTS, MULTIPLE_ORDER, by_name, completed, load, ordered

HERE = Path(__file__).resolve().parent
PAPER = HERE.parent.parent
OUT = PAPER / "figures" / "large_db.png"

# Search modes are categorical: prefilter, unfiltered, and prefilter with
# monoisotopic precursors only. Stopped or refused runs are drawn as crosses at
# the memory limit so the axis still shows where each mode ended.
MODES = (
    ("scaling", True, "Prefilter", "#2a78d6", "-", "o"),
    ("scaling", False, "Unfiltered", "#eb6834", (0, (4, 2)), "s"),
    ("narrow", True, "Prefilter, monoisotopic", "#1c8a5a", (0, (1.5, 1.5)), "D"),
)
TICKS = ("0", "1", "3", "10", "30", "100", "all")
LIMIT_GIB = 16
HUMAN, ENTRAPMENT = "#3987e5", "#eb6834"
ANNOTATED, UNANNOTATED = "#3987e5", "#1c5cab"


def position(multiple):
    return MULTIPLE_ORDER.index(multiple)


def mode_lines(ax, summary, value, limit=None):
    for offset, (series, prefilter, _, color, style, marker) in zip((-0.15, 0, 0.15), MODES):
        rows = ordered(summary.get(series, []), prefilter)
        done = [row for row in rows if completed(row)]
        ax.plot([position(row["multiple"]) for row in done], [value(row) for row in done],
                color=color, linestyle=style, marker=marker, markersize=5,
                markerfacecolor="white", linewidth=1.6)
        stopped = [row for row in rows if not completed(row)]
        if limit is not None and stopped:
            ax.scatter([position(row["multiple"]) + offset for row in stopped], [limit] * len(stopped),
                       marker="x", color=color, s=34, linewidths=1.5, zorder=3)
    ax.set_xticks(range(len(TICKS)))
    ax.set_xticklabels(TICKS, fontsize=7.5)
    ax.set_xlim(-0.4, len(TICKS) - 0.6)
    ax.set_xlabel("Catalog added (× human)")


def main() -> int:
    summary = load()
    fig, axes = plt.subplots(2, 3, figsize=(7.4, 5.3), layout="constrained")

    ax = axes[0][0]
    panel(ax, "A", "Peak resident memory")
    mode_lines(ax, summary, lambda row: row["peak_rss_gib"], limit=LIMIT_GIB)
    ax.axhline(LIMIT_GIB, color=MUTED, linewidth=0.8, linestyle=(0, (2, 2)))
    ax.set_ylabel("GiB")

    ax = axes[0][1]
    panel(ax, "B", "Wall time")
    mode_lines(ax, summary, lambda row: row["wall_seconds"] / 60)
    ax.set_ylabel("Minutes")

    ax = axes[0][2]
    panel(ax, "C", "Peptides searched")
    for series, prefilter, _, color, style, marker in MODES:
        if not prefilter:
            continue
        rows = [row for row in ordered(summary.get(series, []), True)
                if completed(row) and row.get("prefilter_counts")]
        ax.plot([position(row["multiple"]) for row in rows],
                [100 * row["prefilter_counts"]["retained"] / row["prefilter_counts"]["streamed"]
                 for row in rows],
                color=color, linestyle=style, marker=marker, markersize=5,
                markerfacecolor="white", linewidth=1.6)
    ax.set_xticks(range(len(TICKS)))
    ax.set_xticklabels(TICKS, fontsize=7.5)
    ax.set_xlim(-0.4, len(TICKS) - 0.6)
    ax.set_ylim(0, 100)
    ax.set_xlabel("Catalog added (× human)")
    ax.set_ylabel("Retained (%)")

    ax = axes[1][0]
    panel(ax, "D", "Accepted peptides")
    rows = [row for row in ordered(summary["scaling"], True) if completed(row)]
    x = [position(row["multiple"]) for row in rows]
    human = [row["human_peptides"] / 1e3 for row in rows]
    entrap = [row["entrapment_peptides"] / 1e3 for row in rows]
    ax.bar(x, human, color=HUMAN, width=0.62, label="Human")
    ax.bar(x, entrap, bottom=human, color=ENTRAPMENT, width=0.62, label="Catalog only")
    for xi, row, top in zip(x, rows, (h + e for h, e in zip(human, entrap)), strict=True):
        if "combined_fdp" in row:
            ax.annotate(f'{row["combined_fdp"]:.1f}%', (xi, top), textcoords="offset points",
                        xytext=(0, 3), ha="center", fontsize=8, color=MUTED)
    ax.set_xticks(x)
    ax.set_xticklabels([TICKS[i] for i in x])
    ax.set_xlabel("Catalog added (× human)")
    ax.set_ylabel("Peptides (thousands)")
    ax.set_ylim(0, 13)
    ax.legend(loc="upper center", ncol=2, fontsize=7.5, frameon=False)

    ax = axes[1][1]
    panel(ax, "E", "Six-frame genomes")
    rows = by_name(summary["six_frame"])
    pairs = (("Annotated", rows["lfq-annotated-prefilter"]),
             ("Six-frame", rows["lfq-six-frame-prefilter"]))
    labels = [label for label, _ in pairs]
    annotated = [row["microbial_annotated"] / 1e3 for _, row in pairs]
    novel = [row["microbial_unannotated"] / 1e3 for _, row in pairs]
    ax.bar(labels, annotated, color=ANNOTATED, width=0.55, label="In annotated proteome")
    ax.bar(labels, novel, bottom=annotated, color=UNANNOTATED, width=0.55,
           label="Not annotated", hatch="///", edgecolor="white", linewidth=0)
    ax.set_xlabel("Microbial database")
    ax.set_ylabel("Microbial (thousands)")
    ax.set_ylim(0, 14.5)
    ax.legend(loc="upper center", fontsize=7.5, frameon=False)

    ax = axes[1][2]
    panel(ax, "F", "Fecal metaproteomes")
    rows = [row for row in summary["metaproteome"] if row["prefilter"] and completed(row)]
    labels = [row["sample"] for row in rows]
    microbial = [row["microbial_peptides"] / 1e3 for row in rows]
    human = [row["human_peptides"] / 1e3 for row in rows]
    ax.bar(labels, microbial, color="#1c8a5a", width=0.55, label="Microbial")
    ax.bar(labels, human, bottom=microbial, color=HUMAN, width=0.55, label="Human")
    ax.set_xlabel("CAMPI sample")
    ax.set_ylabel("Peptides (thousands)")
    ax.set_ylim(0, 1.3 * max(m + h for m, h in zip(microbial, human, strict=True)))
    ax.legend(loc="upper center", ncol=2, fontsize=7.5, frameon=False)

    handles = [Line2D([], [], color=color, linestyle=style, marker=marker, markersize=5,
                      markerfacecolor="white", label=label)
               for _, _, label, color, style, marker in MODES]
    handles.append(Line2D([], [], color=MUTED, marker="x", linestyle="none",
                          label="Refused or stopped"))
    fig.legend(handles=handles, loc="outside upper center", ncol=len(handles), fontsize=8.5)
    for row_axes in axes:
        for ax in row_axes:
            ax.grid(axis="y", color=GRID, linewidth=0.65)

    save_figure(fig, OUT)
    record("fig.large-db", str(OUT.relative_to(PAPER)), kind="figure",
           inputs=LARGE_DB_INPUTS,
           desc="Large-database search cost, prefilter retention, and identifications.")
    vector = OUT.parent / "vector" / OUT.with_suffix(".svg").name
    record("fig.large-db-vector", str(vector.relative_to(PAPER)), kind="figure",
           inputs=LARGE_DB_INPUTS, desc="Scalable companion to the large-database figure.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
