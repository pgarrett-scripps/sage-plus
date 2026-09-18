"""Generate public timing and expanded entrapment release comparisons."""
import shutil
import numpy as np
from matplotlib.lines import Line2D
from _figure_style import (
    plt, COLORS, MUTED, GRID, TEAL, PURPLE, panel, engine_style,
    engine_legend, save_figure,
)
from _assets import record
from _scientific import PAPER, INPUTS, load


def finish(fig, name, desc):
    target = PAPER / "figures" / f"scientific-entrapment-{name}.png"
    save_figure(fig, target)
    record(f"fig.scientific-entrapment-{name}", str(target.relative_to(PAPER)),
           kind="figure", inputs=INPUTS, desc=desc)
    vector = target.parent / "vector" / target.with_suffix(".svg").name
    record(f"fig.scientific-entrapment-{name}-vector", str(vector.relative_to(PAPER)),
           kind="figure", inputs=INPUTS, desc=f"Scalable companion: {desc}")


def main():
    source = PAPER.parent / "benchmarks/scientific-results/20260914/figures/public-timing-pilot.png"
    target = PAPER / "figures/scientific-public-timing-pilot.png"
    shutil.copyfile(source, target)
    record("fig.scientific-public-timing-pilot", str(target.relative_to(PAPER)), kind="figure",
           inputs=["../benchmarks/scientific-results/20260914/figures/public-timing-pilot.png",
                   "../benchmarks/scientific-results/20260914/figures/figures.json"],
           desc="Frozen public release timing comparison")

    pilot, _ = load()
    fig, axes = plt.subplots(2, 2, figsize=(7.2, 6.7), layout="constrained")
    for col, study in enumerate(("human", "hye")):
        ax = axes[0, col]
        for engine in ("upstream", "plus"):
            cells = [r for r in pilot["calibration"]
                     if r["suite"].startswith(f"entrapment-{study}-") and r["job"].endswith("-" + engine)]
            x = np.array([p["nominal_q"] for p in cells[0]["thresholds"]]) * 100
            values = np.array([[p["paired_fdp_tie_max"] * 100 for p in r["thresholds"]] for r in cells])
            for row in values:
                ax.plot(x, row, color=COLORS[engine], alpha=.16, linewidth=.7)
            ax.plot(x, values.mean(axis=0), **engine_style(engine))
        ax.plot([.035, 7], [.035, 7], linestyle=":", color=MUTED, linewidth=1.1, zorder=1)
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xticks([.1, .5, 1, 2, 5], ["0.1", "0.5", "1", "2", "5"])
        ax.set_yticks([.1, .5, 1, 2, 5], ["0.1", "0.5", "1", "2", "5"])
        ax.minorticks_off()
        ax.set_xlim(.08, 6)
        ax.set_ylim(.035, 7)
        ax.set_xlabel("Nominal peptide q-value (%)")
        ax.set_ylabel("Paired FDP (%)")
        panel(ax, "AB"[col], "HEK calibration" if study == "human" else "Mixture calibration", grid="both")

        ax = axes[1, col]
        deltas = []
        for file in (0, 1):
            for seed in (20260914, 20260915, 20260916):
                values = {engine: next(p["paired_fdp_tie_max"] for r in pilot["calibration"]
                    if r["suite"] == f"entrapment-{study}-{seed}" and r["job"] == f"file-{file}-{engine}"
                    for p in r["thresholds"] if p["nominal_q"] == .01) for engine in ("upstream", "plus")}
                deltas.append(100 * (values["plus"] - values["upstream"]))
        positions = [0, 1, 2, 3.5, 4.5, 5.5]
        ax.scatter(deltas, positions, color=TEAL, s=27, zorder=3)
        summary = next(r for r in pilot["calibration_uncertainty"] if r["study"] == study)
        mean = 100 * (summary["means"]["plus"] - summary["means"]["upstream"])
        lo, hi = np.array(summary["percentile_intervals"]["paired_difference"]) * 100
        ax.errorbar(mean, 7, xerr=[[mean - lo], [hi - mean]], fmt="D", color=PURPLE,
                    capsize=3, markersize=5, elinewidth=1.5, zorder=3)
        ax.set_yticks(positions + [7], [f"File {file + 1}, seed {i + 1}"
                     for file in (0, 1) for i in range(3)] + ["Mean + interval"])
        ax.set_ylim(7.8, -.7)
        ax.axvline(0, color=MUTED, linestyle=(0, (3, 3)), linewidth=1)
        ax.set_xlim(-.14, .09)
        ax.set_xticks([-.10, -.05, 0, .05], ["−0.10", "−0.05", "0.00", "+0.05"])
        ax.set_xlabel("Sage Plus − Sage FDP (pp)")
        ax.spines["left"].set_visible(False)
        ax.tick_params(axis="y", length=0)
        panel(ax, "CD"[col], "HEK differences" if study == "human" else "Mixture differences", grid="x")
    engine_legend(fig, extra=[Line2D([], [], color=MUTED, linestyle=":", label="FDP = nominal q")])
    finish(fig, "detail", "Threshold calibration and file-by-seed paired differences for the frozen release pair")

    fig, ax = plt.subplots(figsize=(5.9, 2.5), layout="constrained")
    for i, study in enumerate(("human", "hye")):
        summary = next(r for r in pilot["calibration_uncertainty"] if r["study"] == study)
        mean = 100 * (summary["means"]["plus"] - summary["means"]["upstream"])
        lo, hi = np.array(summary["percentile_intervals"]["paired_difference"]) * 100
        ax.errorbar(mean, i, xerr=[[mean - lo], [hi - mean]], fmt="D", color=TEAL,
                    capsize=4, markersize=6.5, elinewidth=1.7, zorder=3)
        ax.annotate(f"{mean:+.3f}".replace("-", "−"), (mean, i), xytext=(0, 11),
                    textcoords="offset points", ha="center", fontsize=10, color=TEAL)
    ax.set_yticks([0, 1], ["HEK", "Mixture"])
    ax.set_ylim(1.5, -.65)
    ax.set_xlim(-.10, .10)
    ax.set_xticks([-.10, -.05, 0, .05, .10], ["−0.10", "−0.05", "0.00", "+0.05", "+0.10"])
    ax.axvline(0, color=MUTED, linestyle=(0, (3, 3)), linewidth=1)
    ax.set_xlabel("Sage Plus − Sage FDP (percentage points)")
    ax.grid(axis="x", color=GRID, linewidth=.65)
    ax.set_axisbelow(True)
    ax.spines["left"].set_visible(False)
    ax.tick_params(axis="y", length=0, pad=8)
    finish(fig, "pilot", "Study mean paired FDP differences and conditional bootstrap intervals at one percent")


if __name__ == "__main__":
    main()
