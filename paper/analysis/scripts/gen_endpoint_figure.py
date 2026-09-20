"""Plot comparative endpoints while retaining exact values in companion tables."""
import numpy as np
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from matplotlib.ticker import PercentFormatter, StrMethodFormatter
from _assets import record
from _figure_style import (
    plt, COLORS, INK, MUTED, engine_style, panel, save_figure,
)
from _scientific import PAPER, INPUTS, ENGINE, load, timing
from _report import report


def finish(fig, name, inputs, desc):
    target = PAPER / "figures" / f"endpoint-{name}.png"
    save_figure(fig, target)
    record(f"fig.endpoint-{name}", str(target.relative_to(PAPER)), kind="figure",
           inputs=inputs, desc=desc)
    vector = target.parent / "vector" / target.with_suffix(".svg").name
    record(f"fig.endpoint-{name}-vector", str(vector.relative_to(PAPER)), kind="figure",
           inputs=inputs, desc=f"Scalable companion: {desc}")


def public_timing():
    pilot, _ = load()
    summaries = timing(pilot)
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.7), sharey=True, layout="constrained")
    positions = [0, .65, 1.85, 2.5]
    for ax, metric, raw, title, letter in zip(
        axes, ("seconds", "rss"), ("wall_seconds", "peak_rss_mib"),
        ("Wall time", "Peak resident memory"), "AB",
    ):
        for y, row in zip(positions, summaries):
            engine = row["engine"]
            values = [r[raw] for r in pilot["jobs"] if r["suite"] == "public-timing"
                      and r["study"] == row["study"] and r["engine"] == engine
                      and not r["warmup"]]
            low, high = min(values), max(values)
            ax.hlines(y, low, high, color=COLORS[engine], linewidth=1.4)
            ax.vlines([low, high], y - .055, y + .055, color=COLORS[engine], linewidth=1)
            ax.scatter(values, np.array([-.10, 0, .10]) + y, s=13,
                       color=COLORS[engine], alpha=.5, zorder=3)
            ax.plot(row[metric], y, **{**engine_style(engine), "linestyle": "none"}, zorder=4)
            label = f"{row[metric]:.2f}" if metric == "seconds" else f"{row[metric]:,.1f}"
            ax.annotate(label, (high, y), xytext=(9, 0), textcoords="offset points",
                        ha="left", va="center", fontsize=10, color=INK)
        ax.set_yticks(positions, ["HEK · Sage", "HEK · Sage Plus",
                                 "Mixture · Sage", "Mixture · Sage Plus"])
        ax.set_ylim(2.9, -.4)
        ax.set_xlim(0, 35 if metric == "seconds" else 4050)
        ax.set_xticks([0, 10, 20, 30] if metric == "seconds" else [0, 1000, 2000, 3000, 4000])
        ax.xaxis.set_major_formatter(StrMethodFormatter("{x:,.0f}"))
        ax.set_xlabel("Seconds" if metric == "seconds" else "MiB")
        ax.spines["left"].set_visible(False)
        ax.tick_params(axis="y", length=0, pad=7)
        panel(ax, letter, title, grid="x")
    finish(fig, "timing", INPUTS, "Public timing and memory medians with measured trials and ranges")


def agreement():
    pilot, _ = load()
    rows = [r for r in pilot["public_identifications"] if r["status"] == "complete"]
    names = ["HEK 1", "HEK 2", "A Alpha", "B Alpha", "A Beta", "B Beta"]
    positions = np.array([0, 1, 2.5, 3.5, 4.5, 5.5])
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 4.25), sharey=True, layout="constrained")
    shared_color = "#CBD5DC"
    totals = np.array([r["shared_target_psms"] + r["baseline_only"] + r["candidate_only"] for r in rows])
    left = np.zeros(len(rows))
    for key, color in (("baseline_only", COLORS["upstream"]),
                       ("shared_target_psms", shared_color), ("candidate_only", COLORS["plus"])):
        values = np.array([r[key] for r in rows]) / totals
        axes[0].barh(positions, values, left=left, height=.66, color=color)
        left += values
    for y, row in zip(positions, rows):
        axes[0].text(.5, y, f"{100 * row['jaccard']:.2f}% shared", ha="center", va="center",
                     fontsize=9.5, color=INK)
    axes[0].set_xlim(0, 1)
    axes[0].set_xticks([0, .25, .5, .75, 1])
    axes[0].xaxis.set_major_formatter(PercentFormatter(1))
    axes[0].set_xlabel("Share of accepted PSM union")
    panel(axes[0], "A", "Shared and engine-only PSMs", grid=None)
    for key, engine, offset in (("baseline_only", "upstream", -.18),
                                ("candidate_only", "plus", .18)):
        values = [r[key] for r in rows]
        bars = axes[1].barh(positions + offset, values, height=.29, color=COLORS[engine])
        axes[1].bar_label(bars, labels=[f"{x:,}" for x in values], padding=4, fontsize=9)
    axes[1].set_xlim(0, 1320)
    axes[1].set_xticks([0, 500, 1000])
    axes[1].xaxis.set_major_formatter(StrMethodFormatter("{x:,.0f}"))
    axes[1].set_xlabel("Engine-only accepted PSMs")
    panel(axes[1], "B", "Engine-only counts", grid="x")
    for ax in axes:
        ax.set_yticks(positions, names)
        ax.set_ylim(6.1, -.6)
        ax.spines["left"].set_visible(False)
        ax.tick_params(axis="y", length=0, pad=7)
    fig.legend(handles=[Patch(facecolor=COLORS["upstream"], label="Sage only"),
                        Patch(facecolor=shared_color, label="Shared"),
                        Patch(facecolor=COLORS["plus"], label="Sage Plus only")],
               loc="outside upper center", ncol=3)
    finish(fig, "agreement", INPUTS, "Shared PSM fraction and engine-only counts for each public file")


def lfq_endpoints():
    rows = {r["engine"]: r["species"] for r in report("lfq")["engines"]}
    species = ("human", "yeast", "ecoli")
    fig, axes = plt.subplots(2, 2, figsize=(7.2, 5.35), sharey=True, layout="constrained")
    settings = [
        ("median_log2_bias", 1, "Ratio bias", "Median log₂(B/A) − expected", (-.25, .25), [-.2, 0, .2]),
        ("median_absolute_log2_error", 1, "Absolute ratio error", "Median absolute log₂ error", (0, .55), [0, .25, .5]),
        ("median_preparation_cv", 100, "Preparation variability", "Median preparation CV (%)", (0, 35), [0, 10, 20, 30]),
        ("missing_fraction_observed_union", 100, "Missingness", "Missing from observed union (%)", (0, 1), [0, .5, 1]),
    ]
    for ax, letter, (metric, scale, title, xlabel, limits, ticks) in zip(axes.flat, "ABCD", settings):
        for engine, offset in (("upstream", -.11), ("plus", .11)):
            values = [scale * rows[engine][sp][metric] for sp in species]
            ax.plot(values, np.arange(3) + offset,
                    **{**engine_style(engine), "linestyle": "none", "markersize": 6},
                    clip_on=False)
        ax.set_yticks(range(3), ["Human", "Yeast", r"$\it{E.\ coli}$"])
        ax.set_ylim(2.5, -.5)
        ax.set_xlim(limits)
        ax.set_xticks(ticks)
        ax.set_xlabel(xlabel)
        ax.spines["left"].set_visible(False)
        ax.tick_params(axis="y", length=0, pad=6)
        panel(ax, letter, title, grid="x")
        if metric == "median_log2_bias":
            ax.axvline(0, color=MUTED, linewidth=.9, linestyle=":")
    handles = [Line2D([], [], label=ENGINE[e],
                     **{**engine_style(e), "linestyle": "none"}) for e in ("upstream", "plus")]
    fig.legend(handles=handles, loc="outside upper center", ncol=2)
    finish(fig, "lfq", ["analysis/data/report-extension/lfq.json"],
           "Species-level LFQ bias, absolute error, preparation variability, and missingness")


if __name__ == "__main__":
    public_timing()
    agreement()
    lfq_endpoints()
