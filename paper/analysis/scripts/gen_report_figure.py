"""Generate comparative release figures from the audited report snapshots."""
import numpy as np
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from matplotlib.ticker import PercentFormatter, ScalarFormatter, StrMethodFormatter
from _figure_style import (
    plt, COLORS, INK, MUTED, GRID, TEAL, PURPLE, engine_style,
    engine_legend, panel, save_figure,
)
from _assets import record
from _scientific import PAPER, INPUTS, ENGINE
from _report import report, workloads


def finish(fig, name, desc):
    target = PAPER / "figures" / f"report-{name}.png"
    save_figure(fig, target)
    sections = {"identification": "public", "disagreement": "public",
                "scaling": "scaling", "ptm": "ptm", "lfq": "lfq", "control": "lfq"}
    inputs = INPUTS if name == "workloads" else [
        f"analysis/data/report-extension/{sections[name]}.json"]
    record(f"fig.report-{name}", str(target.relative_to(PAPER)), kind="figure",
           inputs=inputs, desc=desc)
    vector = target.parent / "vector" / target.with_suffix(".svg").name
    record(f"fig.report-{name}-vector", str(vector.relative_to(PAPER)), kind="figure",
           inputs=inputs, desc=f"Scalable companion: {desc}")


def identification():
    rows = report("public")
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.4), sharey=True, layout="constrained")
    for ax, study, letter in zip(axes, ("PXD001468", "PXD028735"), "AB"):
        selected = [r for r in rows if study in r["pair"]]
        qs = (.001, .005, .01, .02, .05)
        for metric, marker, text, color in (
            ("target_psms", "o", "PSMs", TEAL),
            ("target_peptidoforms", "s", "Peptidoforms", PURPLE),
        ):
            curves = np.array([[100 * (r["engines"]["plus"][str(q)][metric] /
                r["engines"]["upstream"][str(q)][metric] - 1) for q in qs] for r in selected])
            ax.fill_between(np.array(qs) * 100, curves.min(axis=0), curves.max(axis=0),
                            color=color, alpha=.12, linewidth=0)
            ax.plot(np.array(qs) * 100, curves.mean(axis=0), marker=marker,
                    markersize=5, color=color, linewidth=1.6, label=text)
        ax.axhline(0, color=MUTED, linewidth=.9)
        ax.set_xscale("log")
        ax.set_xticks([.1, .5, 1, 2, 5], ["0.1", "0.5", "1", "2", "5"])
        ax.minorticks_off()
        ax.set_xlabel("Nominal q-value (%)")
        panel(ax, letter, "HEK" if study == "PXD001468" else "Mixture")
    axes[0].set_ylabel("Sage Plus yield change (%)")
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="outside upper center", ncol=2)
    finish(fig, "identification", "Public PSM and peptidoform nominal-threshold response")


def disagreement():
    rows = report("public")
    names = ["HEK 1", "HEK 2", "A Alpha", "B Alpha", "A Beta", "B Beta"]
    positions = np.array([0, 1, 2.5, 3.5, 4.5, 5.5])
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.8), sharey=True, layout="constrained")
    categories = [
        ("same_assignment_above_threshold", "Same assignment,\nabove threshold", "#75A7B2"),
        ("different_assignment", "Different rank-one\nassignment", "#CFA34C"),
        ("no_matching_rank_one_spectrum_charge", "No matched spectrum\nand charge", "#9B86AD"),
    ]
    for ax, engine, letter in zip(axes, ("upstream", "plus"), "AB"):
        totals = np.array([sum(r["disagreement"][engine].values()) for r in rows])
        left = np.zeros(len(rows))
        for key, text, color in categories:
            values = np.array([r["disagreement"][engine].get(key, 0) for r in rows]) / totals
            ax.barh(positions, values, left=left, height=.67, color=color,
                    edgecolor="white", linewidth=.3, label=text)
            left += values
        ax.set_yticks(positions, names)
        ax.set_ylim(6.05, -.65)
        ax.set_xlim(0, 1.26)
        ax.set_xticks([0, .25, .5, .75, 1])
        ax.xaxis.set_major_formatter(PercentFormatter(1))
        ax.spines["bottom"].set_bounds(0, 1)
        ax.spines["left"].set_visible(False)
        ax.tick_params(axis="y", length=0, pad=6)
        for position, total in zip(positions, totals):
            ax.text(1.035, position, f"n={total:,}", va="center", fontsize=9, color=MUTED)
        panel(ax, letter, f"{ENGINE[engine]} only", grid=None)
        ax.set_xlabel("Share of engine-only PSMs")
    handles = [Patch(facecolor=color, label=text) for _, text, color in categories]
    fig.legend(handles=handles, loc="outside lower center", ncol=3, fontsize=9,
               handlelength=1.2, columnspacing=1.5)
    finish(fig, "disagreement", "Decomposition of accepted PSM disagreement using raw rank-one output")


def tradeoff():
    rows = [r for r in workloads() if r["kind"] != "Contextual local timing"]
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.45), sharey=True, layout="constrained")
    for ax, key, title, letter in zip(axes, ("seconds", "rss"),
                                     ("Wall time", "Peak resident memory"), "AB"):
        values = [r["engines"]["plus"][key] / r["engines"]["upstream"][key] for r in rows]
        for i, value in enumerate(values):
            color, marker = (TEAL, "o") if i < 2 else (PURPLE, "D")
            ax.plot([1, value], [i, i], color=color, linewidth=1.8, zorder=2)
            ax.scatter(value, i, color=color, marker=marker, s=45, zorder=3)
            ax.annotate(f"{value:.2f}", (value, i), xytext=(0, 9),
                        textcoords="offset points", ha="center", fontsize=10, color=INK)
        ax.set_yticks(range(len(rows)), [r["label"] for r in rows])
        ax.set_ylim(3.5, -.7)
        ax.axvline(1, color=MUTED, linewidth=1, linestyle=(0, (3, 3)))
        ax.set_xlim(.6, 1.3)
        ax.set_xticks([.6, .8, 1, 1.2])
        ax.set_xlabel("Sage Plus / Sage")
        ax.spines["left"].set_visible(False)
        ax.tick_params(axis="y", length=0, pad=7)
        panel(ax, letter, title, grid="x")
    fig.legend(handles=[
        Line2D([], [], marker="o", linestyle="none", color=TEAL, label="Repeated public searches"),
        Line2D([], [], marker="D", linestyle="none", color=PURPLE, label="Across files and seeds"),
    ], loc="outside lower center", ncol=2, fontsize=9.5)
    finish(fig, "workloads", "Within-workload performance ratios with distinct replication designs")


def scaling():
    rows = report("scaling")
    workers = (1, 2, 4, 8)
    fig, axes = plt.subplots(1, 3, figsize=(7.5, 3.2), layout="constrained")
    for engine in ("upstream", "plus"):
        baseline = np.median([r["seconds"] for r in rows if r["engine"] == engine
                              and r["threads"] == 1 and not r["warmup"]])
        for ax, key in zip(axes, ("seconds", "rss", "speedup")):
            groups = [[r for r in rows if r["engine"] == engine and r["threads"] == t
                       and not r["warmup"]] for t in workers]
            values = [np.median([r["seconds" if key == "speedup" else key] for r in group])
                      for group in groups]
            if key == "speedup":
                values = baseline / np.array(values)
            elif key == "rss":
                values = np.array(values) / 1024
            if key != "speedup":
                for t, group in zip(workers, groups):
                    ax.scatter([t] * len(group), [r[key] / (1024 if key == "rss" else 1)
                               for r in group], s=12, color=COLORS[engine], alpha=.35, zorder=2)
            ax.plot(workers, values, **engine_style(engine), zorder=3)
    for ax, letter, title, ylabel in zip(axes, "ABC",
        ("Wall time", "Peak memory", "Speedup"),
        ("Seconds", "Resident memory (GiB)", "Relative to one worker")):
        panel(ax, letter, title)
        ax.set_xscale("log", base=2)
        ax.set_xticks(workers)
        ax.xaxis.set_major_formatter(ScalarFormatter())
        ax.minorticks_off()
        ax.set_xlim(.85, 9.4)
        ax.set_xlabel("Workers")
        ax.set_ylabel(ylabel)
    axes[0].set_ylim(0, 44)
    axes[0].set_yticks([0, 10, 20, 30, 40])
    axes[1].set_ylim(1.3, 2.2)
    axes[1].set_yticks([1.4, 1.6, 1.8, 2, 2.2])
    axes[2].set_ylim(.8, 4.4)
    axes[2].set_yticks([1, 2, 3, 4], ["1×", "2×", "3×", "4×"])
    engine_legend(fig)
    finish(fig, "scaling", "Matched release thread scaling with warmups excluded")


def ptm():
    rows = report("ptm")
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.25), sharey=True, layout="constrained")
    for ax, library, letter in zip(axes, (1, 2), "AB"):
        for engine, offset in (("upstream", -.19), ("plus", .19)):
            row = next(r for r in rows if r["engine"] == engine and r["library"] == library)
            values = [row["spectrum_accepted"], row["joint_accepted"]]
            bars = ax.bar(np.arange(2) + offset, values, .32, color=COLORS[engine])
            ax.bar_label(bars, labels=[f"{value:,}" for value in values], padding=5, fontsize=10)
            ax.scatter([1 + offset], [0], marker="_", color=COLORS[engine], s=100, clip_on=False)
        ax.set_xticks([0, 1], ["Spectrum only", "Spectrum + peptide"])
        ax.set_ylim(0, max(r["spectrum_accepted"] for r in rows) * 1.18)
        ax.set_yticks([0, 1000, 2000, 3000])
        ax.yaxis.set_major_formatter(StrMethodFormatter("{x:,.0f}"))
        ax.set_xlabel("Acceptance filters at q ≤ 1%")
        panel(ax, letter, f"Synthetic HCD {library}")
    axes[0].set_ylabel("Accepted target PSMs")
    fig.legend(handles=[Patch(facecolor=COLORS[e], label=ENGINE[e])
                        for e in ("upstream", "plus")], loc="outside upper center", ncol=2)
    finish(fig, "ptm", "Synthetic PTM acceptance under spectrum-only and joint peptide filters")


def quantification():
    quant = report("lfq")
    fig, axes = plt.subplots(1, 3, figsize=(7.5, 3.15), sharey=True, layout="constrained")
    for ax, species, expected, letter, title in zip(axes, ("human", "yeast", "ecoli"),
        (0, -1, 2), "ABC", ("Human", "Yeast", r"$\it{E.\ coli}$")):
        for row in quant["engines"]:
            values = list(row["ratios"][species].values())
            ax.hist(values, bins=np.linspace(-5, 6, 89), density=False,
                    weights=np.ones(len(values)) / len(values), histtype="step",
                    **engine_style(row["engine"], markers=False))
        ax.axvline(expected, color=MUTED, linestyle=":", linewidth=1.2, zorder=1)
        ax.set_xlim(expected - 2, expected + 2)
        ax.set_xticks(np.arange(expected - 2, expected + 3))
        ax.set_ylim(0, .20)
        ax.set_yticks([0, .05, .10, .15, .20])
        ax.yaxis.set_major_formatter(PercentFormatter(1, decimals=0))
        panel(ax, letter, title)
    axes[0].set_ylabel("Fraction of ratio pairs per bin")
    fig.supxlabel("Observed log₂(B/A)", fontsize=10.5)
    engine_legend(fig, markers=False, extra=[Line2D([], [], color=MUTED,
                  linestyle=":", linewidth=1.2, label="Expected ratio")])
    finish(fig, "lfq", "Matched LFQ species-ratio distributions")


def control():
    rows = report("lfq")["engines"]
    fig, ax = plt.subplots(figsize=(5.8, 3.35), layout="constrained")
    for engine, offset in (("upstream", -.16), ("plus", .16)):
        counts = next(r["control"] for r in rows if r["engine"] == engine)
        fractions = [counts["foreign"] / counts["quantified"],
                     counts["foreign_without_strict_ms2"] / counts["without_strict_ms2"]]
        bars = ax.bar(np.arange(2) + offset, fractions, .27, color=COLORS[engine])
        ax.bar_label(bars, labels=[f"{value:.1%}" for value in fractions], padding=5,
                     fontsize=11)
    ax.set_xticks([0, 1], ["All accepted\ncontrol rows", "Without direct\nMS2 support"])
    ax.set_ylabel("Foreign-species fraction")
    ax.set_xlim(-.6, 1.6)
    ax.set_ylim(0, .8)
    ax.set_yticks([0, .2, .4, .6, .8])
    ax.yaxis.set_major_formatter(PercentFormatter(1))
    ax.grid(axis="y", color=GRID, linewidth=.65)
    ax.set_axisbelow(True)
    fig.legend(handles=[Patch(facecolor=COLORS[e], label=ENGINE[e])
                        for e in ("upstream", "plus")], loc="outside upper center", ncol=2)
    finish(fig, "control", "Foreign-species fractions in all accepted and MS2-unsupported human-only control rows")


if __name__ == "__main__":
    for function in (identification, disagreement, tradeoff, scaling, ptm, quantification, control):
        function()
