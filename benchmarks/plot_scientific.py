#!/usr/bin/env python3
"""Generate pilot figures directly from the machine-readable scientific summary."""

import argparse
import json
import statistics
from collections import defaultdict
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from provenance import atomic_json, sha256


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("summary", type=Path)
    p.add_argument("--output", type=Path, required=True)
    args = p.parse_args()
    data = json.loads(args.summary.read_text())
    args.output.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10,
                         "axes.spines.top": False, "axes.spines.right": False,
                         "svg.hashsalt": "sage-plus-scientific-pilot", "savefig.dpi": 250})
    outputs = {}

    def save(fig, name, caption):
        fig.tight_layout()
        for suffix in ("png", "svg"):
            path = args.output / f"{name}.{suffix}"
            fig.savefig(path, metadata={"Date": None} if suffix == "svg" else {"Software": None})
            outputs[path.name] = {"sha256": sha256(path), "caption": caption}
        plt.close(fig)

    if data["calibration"]:
        fig, axes = plt.subplots(1, 2, figsize=(8, 3.4))
        for axis, study in zip(axes, ("human", "hye")):
            for engine, color, label in (("upstream", "#276B9B", "Upstream Sage"), ("plus", "#C85B29", "Sage Plus beta.3")):
                grouped = defaultdict(list)
                for result in data["calibration"]:
                    if result["suite"].startswith(f"entrapment-{study}-") and result["job"].endswith(engine):
                        for point in result["thresholds"]:
                            if point["paired_fdp_tie_max"] is not None:
                                grouped[100 * point["nominal_q"]].append(100 * point["paired_fdp_tie_max"])
                x = sorted(grouped)
                if x:
                    axis.plot(x, [statistics.mean(grouped[q]) for q in x], "o-", color=color, label=label, markersize=4)
                    for q in x:
                        axis.scatter([q] * len(grouped[q]), grouped[q], color=color, alpha=0.25, s=10)
            axis.plot([0, 5], [0, 5], color="#888888", linestyle=":", linewidth=1)
            axis.set(title="HEK DDA" if study == "human" else "Mixed-species DDA",
                     xlabel="Reported peptide threshold (%)", ylabel="Paired entrapment FDP (%)")
            axis.set_xlim(0, 5.1)
            axis.set_ylim(bottom=0)
        axes[0].legend(frameon=False, fontsize=8)
        save(fig, "entrapment-pilot", "Lines show means and points show individual selected file and seed results. The conservative end of the exact-score tie interval is used. These are pilot FDP estimates, not production calibration certificates.")

    scaling = [r for r in data["jobs"] if r["suite"] == "scaling" and r["status"] == "complete"
               and not r["warmup"] and r["workload"] == "thread-scaling-and-exact-prefilter"]
    if scaling:
        fig, axes = plt.subplots(1, 2, figsize=(8, 3.4))
        for prefilter, color in ((0, "#276B9B"), (1, "#C85B29")):
            grouped = defaultdict(list)
            for row in scaling:
                if f"prefilter-{prefilter}-" in row["id"]:
                    grouped[row["threads"]].append(row)
            x = sorted(grouped)
            for axis, metric, scale in ((axes[0], "wall_seconds", 1), (axes[1], "peak_rss_mib", 1 / 1024)):
                medians = [statistics.median(r[metric] * scale for r in grouped[n]) for n in x]
                lows = [medians[i] - min(r[metric] * scale for r in grouped[n]) for i, n in enumerate(x)]
                highs = [max(r[metric] * scale for r in grouped[n]) - medians[i] for i, n in enumerate(x)]
                axis.errorbar(x, medians, yerr=[lows, highs], fmt="o-", capsize=3, color=color,
                              label="Prefilter on" if prefilter else "Prefilter off")
                axis.set_xticks([1, 2, 4, 8])
                axis.set_xlabel("Rayon threads")
        axes[0].set_ylabel("Whole-search wall time (s)")
        axes[1].set_ylabel("Peak resident memory (GiB)")
        axes[0].legend(frameon=False)
        save(fig, "scaling-pilot", "Public HEK input and complete reviewed human reference. Points are medians and bars span measured timing repetitions. Warmups and failed jobs are excluded from these curves and retained in the result matrix.")

    timing = [r for r in data["jobs"] if r["suite"] == "public-timing"
              and r["status"] == "complete" and not r["warmup"]]
    if timing:
        fig, axes = plt.subplots(1, 2, figsize=(8, 3.8))
        for offset, engine, color, label in ((-0.18, "upstream", "#276B9B", "Upstream Sage"),
                                            (0.18, "plus", "#C85B29", "Sage Plus beta.3")):
            for axis, metric, scale in ((axes[0], "wall_seconds", 1), (axes[1], "peak_rss_mib", 1 / 1024)):
                x, medians, lows, highs = [], [], [], []
                for index, study in enumerate(("PXD001468", "PXD028735")):
                    values = [r[metric] * scale for r in timing if r["engine"] == engine and r["study"] == study]
                    if values:
                        median = statistics.median(values)
                        x.append(index + offset)
                        medians.append(median)
                        lows.append(median - min(values))
                        highs.append(max(values) - median)
                axis.bar(x, medians, width=0.35, yerr=[lows, highs], capsize=3, color=color, label=label)
                axis.set_xticks([0, 1], ["HEK DDA\nPXD001468", "Mixture DDA\nPXD028735"])
        axes[0].set_ylabel("Whole-search wall time (s)")
        axes[1].set_ylabel("Peak resident memory (GiB)")
        axes[0].legend(frameon=False, fontsize=8)
        save(fig, "public-timing-pilot", "One selected file from each public study. Bars are medians and error bars span three measured technical repetitions, after one warmup per engine. Engine order alternates. Acquisition, conversion and primary search matrices finish before this timing experiment. The mixed reference excludes seven proteins with undefined residues identically for both engines.")

    if data["ptm"]:
        fig, axes = plt.subplots(1, 2, figsize=(9, 4))
        fig.suptitle("Restricted synthetic pilot: all retained site rows have peptide q = 1", fontsize=10)
        order = ["oracle-0", "oracle-25", "oracle-50", "oracle-100", "all"]
        for library, color in ((1, "#276B9B"), (2, "#C85B29")):
            rows = {next(k for k in order if r["job"].endswith(k)): r for r in data["ptm"] if r["job"].startswith(f"library-{library}-")}
            keys = [k for k in order if k in rows]
            x = [order.index(k) for k in keys]
            axes[0].plot(x, [rows[k]["correct_site_events"] for k in keys], "o-", color=color, label=f"HCD {library}")
            valid = [k for k in keys if rows[k]["empirical_site_error_fraction"] is not None]
            axes[1].plot([order.index(k) for k in valid], [100 * rows[k]["empirical_site_error_fraction"] for k in valid], "o-", color=color)
        for axis in axes:
            axis.set_xticks(range(5), ["0%", "25%", "50%", "100%", "All\nsites"], fontsize=9)
            axis.set_xlabel("Oracle coverage or unrestricted sites")
        axes[0].set_ylabel("Synthesis-consistent site events")
        axes[1].set_ylabel("Synthesis-inconsistent site events (%)")
        axes[0].legend(frameon=False)
        save(fig, "ptm-pilot", "Synthesis consistency at 1% reported PSM and localization thresholds. The restricted synthetic searches trigger heuristic rescoring and yield peptide q-values of one. Oracle coverage uses the known truth and is not a biological annotation benchmark. Empty denominators are omitted.")

    if data["quantification"]:
        fig, axes = plt.subplots(1, 2, figsize=(8, 4))
        species = ["human", "yeast", "ecoli"]
        empty = [r["job"] for r in data["quantification"]
                 if all(r["species"][s]["observable_features"] == 0 for s in species)]
        fig.suptitle("Primary 1% LFQ threshold" + (". No accepted features: " + ", ".join(empty) if empty else ""), fontsize=10)
        for offset, result, color in zip((-0.18, 0.18), data["quantification"], ("#276B9B", "#C85B29")):
            for axis, metric in ((axes[0], "median_log2_bias"), (axes[1], "missing_fraction_observed_union")):
                x, y = [], []
                for i, s in enumerate(species):
                    value = result["species"][s][metric]
                    if value is not None:
                        x.append(i + offset)
                        y.append(value)
                axis.bar(x, y, width=0.35, color=color, label=result["job"].replace("mbr-0", "MBR off").replace("mbr-1", "MBR on"))
                axis.set_xticks(range(3), ["Human", "Yeast", "E. coli"])
        axes[0].axhline(0, color="#888888", linewidth=1)
        axes[0].set_ylabel("Median log2(B/A) bias")
        axes[1].set_ylabel("Missing fraction in observed feature union")
        axes[0].legend(frameon=False)
        save(fig, "quantification-pilot", "Ratios use matched Alpha and Beta preparations, unambiguous species assignments and no imputation. Missingness uses the observed feature union. Pure-human absence-control counts are reported separately in the summary.")
        diagnostics = [r for r in data["quantification"] if "direct_ms2_ratio_diagnostic" in r]
        if diagnostics:
            fig, axes = plt.subplots(1, 2, figsize=(8, 4))
            fig.suptitle("Exploratory MS2-supported ratios before LFQ filtering", fontsize=10)
            for offset, result, color in zip((-0.18, 0.18), diagnostics, ("#276B9B", "#C85B29")):
                for axis, metric in ((axes[0], "median_log2_bias"), (axes[1], "median_absolute_log2_error")):
                    x, y = [], []
                    for index, name in enumerate(species):
                        value = result["direct_ms2_ratio_diagnostic"]["species"][name][metric]
                        if value is not None:
                            x.append(index + offset)
                            y.append(value)
                    axis.bar(x, y, width=0.35, color=color, label=result["job"].replace("mbr-0", "MBR off").replace("mbr-1", "MBR on"))
                    axis.set_xticks(range(3), ["Human", "Yeast", "E. coli"])
            axes[0].axhline(0, color="#888888", linewidth=1)
            axes[0].set_ylabel("Median log2(B/A) bias")
            axes[1].set_ylabel("Median absolute log2 ratio error")
            axes[0].legend(frameon=False)
            save(fig, "quantification-diagnostic", "Post hoc diagnostic added after the empty MBR-off discovery set at the primary LFQ threshold. Positive intensities require direct MS2 evidence passing both 1% PSM and peptide thresholds. LFQ q-values are deliberately not filtered. Cross-species I/L ambiguity is excluded. These are not accepted 1% LFQ discoveries or evidence of calibrated LFQ confidence.")
    atomic_json(args.output / "figures.json", {"summary_sha256": sha256(args.summary),
                "generator_sha256": sha256(Path(__file__).resolve()), "outputs": outputs})
    print(f"Generated {len(outputs)} figure files")


if __name__ == "__main__":
    main()
