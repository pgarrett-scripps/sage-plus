"""Figures for the head-to-head. Run analyze.py first.

uv run --no-project --with pyarrow --with pandas --with matplotlib python \
    plot.py
"""
import json
import os
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

HERE = Path(os.environ.get("HEADTOHEAD_DIR", Path(__file__).resolve().parent))
OUT = Path(__file__).resolve().parents[2] / "figures" / "headtohead"
OUT.mkdir(parents=True, exist_ok=True)

ENGINES = ["upstream", "plus"]
LABEL = {"upstream": "Sage v0.15.0-beta.2", "plus": "Sage Plus Beta 16 RC"}
COLOR = {"upstream": "#2a78d6", "plus": "#eb6834"}  # blue / orange, CVD-separable
INK, INK2, GRID = "#0b0b0b", "#52514e", "#e4e3df"
DATASET = "PXD028735 HYE, Orbitrap QE-HFX DDA, 4 files (A/B x 2)"

plt.rcParams.update({
    "font.size": 10, "axes.edgecolor": INK2, "axes.labelcolor": INK, "xtick.color": INK2,
    "ytick.color": INK2, "axes.spines.top": False, "axes.spines.right": False,
    "figure.facecolor": "white", "axes.facecolor": "white", "svg.fonttype": "none",
})

m = json.loads((HERE / "metrics.json").read_text())


def save(fig, name):
    fig.savefig(OUT / f"{name}.png", dpi=150, bbox_inches="tight")
    fig.savefig(OUT / f"{name}.svg", bbox_inches="tight")
    plt.close(fig)


# Figure 1: identifications + resources, small multiples (one scale per panel).
panels = [
    ("PSMs\n1% spectrum q", lambda e: m[e]["ids"]["total"]["psms"], None, "{:,.0f}"),
    ("Peptides\n1% peptide q", lambda e: m[e]["ids"]["total"]["peptides"], None, "{:,.0f}"),
    ("Protein groups\n1% protein-group q", lambda e: m[e]["ids"]["total"]["protein_groups"], None, "{:,.0f}"),
    ("Wall time (s)\nmedian of 3, lower is better", lambda e: m[e]["wall_s_median"],
     lambda e: [t["wall_s"] for t in m[e]["timings"].values()], "{:.0f} s"),
    ("Peak RSS (GiB)\nmedian of 3, lower is better", lambda e: m[e]["peak_rss_gib_median"],
     lambda e: [t["peak_rss_gib"] for t in m[e]["timings"].values()], "{:.1f}"),
]
fig, axes = plt.subplots(1, 5, figsize=(13, 3.6), gridspec_kw={"wspace": 0.45})
for ax, (title, val, runs, fmt) in zip(axes, panels):
    xs = np.arange(len(ENGINES))
    vals = [val(e) for e in ENGINES]
    ax.bar(xs, vals, width=0.62, color=[COLOR[e] for e in ENGINES], edgecolor="white", linewidth=2)
    if runs:
        for x, e in zip(xs, ENGINES):
            r = runs(e)
            ax.scatter([x] * len(r), r, s=14, color=INK, zorder=3, linewidths=0)
    top = max(vals)
    for x, v in zip(xs, vals):
        ax.text(x, v + top * 0.02, fmt.format(v), ha="center", va="bottom", color=INK, fontsize=9)
    ax.set_ylim(0, top * 1.15)
    ax.set_xticks([])
    ax.set_title(title, fontsize=9.5, color=INK, loc="left")
    ax.yaxis.grid(True, color=GRID, linewidth=0.8)
    ax.set_axisbelow(True)
    ax.tick_params(axis="y", labelsize=8)
    ax.yaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(
        lambda v, _: f"{v/1000:.0f}k" if v >= 10000 else f"{v:g}"))
handles = [plt.Rectangle((0, 0), 1, 1, color=COLOR[e]) for e in ENGINES]
fig.legend(handles, [LABEL[e] for e in ENGINES], loc="upper center", ncol=2, frameon=False,
           bbox_to_anchor=(0.5, 1.08))
fig.text(0.5, -0.04, f"{DATASET}; one search per engine over all 4 files; target counts from each "
         "engine's own q-values. Dots = individual runs.", ha="center", color=INK2, fontsize=8.5)
save(fig, "fig1_ids_resources")

# Figure 2: LFQ log2(A/B) per species.
r = pd.read_parquet(HERE / "lfq_ratios.parquet")
species = [("human", "Human", 0.0), ("yeast", "Yeast", 1.0), ("ecoli", "E. coli", -2.0)]
fig, ax = plt.subplots(figsize=(7.5, 4.6))
w = 0.34
for i, (s, name, exp) in enumerate(species):
    ax.hlines(exp, i - 0.48, i + 0.48, colors=INK2, linestyles="--", linewidth=1.2, zorder=1)
    for j, e in enumerate(ENGINES):
        v = r.loc[(r.engine == e) & (r.species == s), "log2_a_over_b"].to_numpy()
        pos = i + (j - 0.5) * (w + 0.04)
        vp = ax.violinplot([v], positions=[pos], widths=w, showextrema=False, points=200)
        for b in vp["bodies"]:
            b.set_facecolor(COLOR[e]); b.set_alpha(0.35); b.set_edgecolor("none")
        ax.boxplot([v], positions=[pos], widths=w * 0.35, showfliers=False, patch_artist=True,
                   medianprops={"color": INK, "linewidth": 1.6},
                   boxprops={"facecolor": COLOR[e], "edgecolor": INK, "linewidth": 0.8},
                   whiskerprops={"color": INK, "linewidth": 0.8}, capprops={"color": INK, "linewidth": 0.8})
        sm = m[e]["lfq"]["species"][s]
        ax.text(pos, -0.02, f"n={sm['n']:,}\nmed {sm['median']:+.2f}", ha="center", va="top",
                fontsize=7.5, color=INK2,
                transform=matplotlib.transforms.blended_transform_factory(ax.transData, ax.transAxes))
ax.set_xticks(range(len(species)))
ax.set_xticklabels([f"{n}\n(expected {e:+.0f})" for _, n, e in species])
ax.set_ylabel("log2(A / B), mean of 2 replicates per condition")
ax.set_ylim(-4.5, 3.0)
ax.tick_params(axis="x", pad=30, length=0)
ax.yaxis.grid(True, color=GRID, linewidth=0.8)
ax.set_axisbelow(True)
handles = [plt.Rectangle((0, 0), 1, 1, color=COLOR[e]) for e in ENGINES]
handles.append(plt.Line2D([0], [0], color=INK2, linestyle="--"))
ax.legend(handles, [LABEL[e] for e in ENGINES] + ["expected ratio"], loc="upper right",
          frameon=False, fontsize=8.5)
ax.set_title(f"LFQ ratio accuracy, {DATASET}", fontsize=9.5, loc="left", color=INK)
fig.text(0.5, -0.08, "Peptides (charge states combined) at LFQ q <= 0.01 with intensity in all 4 files, "
         "single-species only; no normalization. Box = IQR, whiskers 1.5 IQR.",
         ha="center", color=INK2, fontsize=8)
save(fig, "fig2_lfq_ratios")
print("wrote", sorted(p.name for p in OUT.iterdir()))
