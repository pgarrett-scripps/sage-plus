"""Metrics for the upstream Sage vs Sage Plus head-to-head on PXD028735 HYE.

Usage (from anywhere):
  uv run --no-project --with pyarrow --with pandas --with matplotlib python \
    analyze.py

Reads <engine>/run*/ outputs, writes metrics.json next to this script.
"""
import json
import os
import math
import re
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(os.environ.get("HEADTOHEAD_DIR", Path(__file__).resolve().parent))
Q = 0.01
FILES = [
    "LFQ_Orbitrap_DDA_Condition_A_Sample_Alpha_01.mzML",
    "LFQ_Orbitrap_DDA_Condition_B_Sample_Alpha_01.mzML",
    "LFQ_Orbitrap_DDA_Condition_A_Sample_Beta_01.mzML",
    "LFQ_Orbitrap_DDA_Condition_B_Sample_Beta_01.mzML",
]
A_FILES = [f for f in FILES if "Condition_A" in f]
B_FILES = [f for f in FILES if "Condition_B" in f]
# quantification-plan.json "primary_ratios": "B over A: human 1, yeast 0.5, E. coli 4"
EXPECTED_LOG2_A_OVER_B = {"human": 0.0, "yeast": 1.0, "ecoli": -2.0}


def parse_time(path):
    txt = path.read_text()
    wall = re.search(r"Elapsed \(wall clock\) time \(h:mm:ss or m:ss\): (\S+)", txt).group(1)
    parts = [float(p) for p in wall.split(":")]
    secs = 0.0
    for p in parts:
        secs = secs * 60 + p
    rss_kb = int(re.search(r"Maximum resident set size \(kbytes\): (\d+)", txt).group(1))
    status = int(re.search(r"Exit status: (\d+)", txt).group(1))
    return {"wall_s": secs, "peak_rss_gib": rss_kb / 1024 / 1024, "exit": status}


def species_of(proteins):
    """Species set from accessions like 'human:sp|P05023|AT1A1_HUMAN' (';'-joined)."""
    out = set()
    for acc in str(proteins).split(";"):
        acc = acc.strip()
        if acc.startswith("rev_"):
            acc = acc[4:]
        out.add(acc.split(":", 1)[0])
    return out


def load_psms(engine, run_dir):
    if engine == "upstream":
        df = pd.read_csv(run_dir / "results.sage.tsv", sep="\t",
                         usecols=["peptide", "proteins", "filename", "label", "rank",
                                  "spectrum_q", "peptide_q", "protein_q",
                                  "protein_groups", "protein_group_q"])
        df["target"] = df["label"] == 1
    else:
        df = pd.read_parquet(run_dir / "results.sage.parquet",
                             columns=["peptide", "proteins", "filename", "is_decoy", "rank",
                                      "spectrum_q", "peptide_q", "protein_q",
                                      "protein_groups", "protein_group_q"])
        df["target"] = ~df["is_decoy"]
    df["filename"] = df["filename"].map(lambda s: Path(s).name)
    return df


def id_metrics(df):
    t = df[df["target"]]
    psm = t[t["spectrum_q"] <= Q]
    pep = t[t["peptide_q"] <= Q]
    prot = psm[psm["protein_q"] <= Q]
    grp = psm[psm["protein_group_q"] <= Q]
    m = {"total": {
        "psms": int(len(psm)),
        "peptides": int(pep["peptide"].nunique()),
        "proteins": int(prot["proteins"].nunique()),
        "protein_groups": int(grp["protein_groups"].nunique()),
    }, "per_file": {}}
    for f in FILES:
        m["per_file"][f] = {
            "psms": int((psm["filename"] == f).sum()),
            "peptides": int(pep.loc[pep["filename"] == f, "peptide"].nunique()),
            "proteins": int(prot.loc[prot["filename"] == f, "proteins"].nunique()),
            "protein_groups": int(grp.loc[grp["filename"] == f, "protein_groups"].nunique()),
        }
    return m


def load_lfq(engine, run_dir):
    """Wide table: peptide, proteins, q_value, one intensity column per file (NaN = missing)."""
    if engine == "upstream":
        df = pd.read_csv(run_dir / "lfq.tsv", sep="\t")
        df.columns = [Path(c).name if c.endswith(".mzML") else c for c in df.columns]
        df = df[~df["proteins"].astype(str).str.split(";").map(
            lambda a: all(x.startswith("rev_") for x in a))]
        for f in FILES:
            df[f] = df[f].where(df[f] > 0)
        return df[["peptide", "charge", "proteins", "q_value"] + FILES]
    long = pd.read_parquet(run_dir / "lfq.parquet")
    long = long[~long["is_decoy"]]
    long["filename"] = long["filename"].map(lambda s: Path(s).name)
    long["charge"] = long["charge"].fillna(-1)
    wide = (long.set_index(["peptide", "charge", "proteins", "q_value", "filename"])["intensity"]
            .unstack("filename").reset_index())
    for f in FILES:
        wide[f] = wide[f].where(wide[f] > 0) if f in wide else np.nan
    return wide[["peptide", "charge", "proteins", "q_value"] + FILES]


def lfq_metrics(df):
    out = {"precursors_q01": int((df["q_value"] <= Q).sum())}
    d = df[df["q_value"] <= Q].copy()
    complete = d[FILES].notna().all(axis=1)
    d = d[complete]
    out["precursors_q01_all4"] = int(len(d))
    sp = d["proteins"].map(species_of)
    d = d[sp.map(len) == 1]
    d["species"] = sp[sp.map(len) == 1].map(lambda s: next(iter(s)))
    d = d[d["species"].isin(EXPECTED_LOG2_A_OVER_B)]
    d["log2_a_over_b"] = np.log2(d[A_FILES].mean(axis=1) / d[B_FILES].mean(axis=1))
    out["species"] = {}
    for s, exp in EXPECTED_LOG2_A_OVER_B.items():
        v = d.loc[d["species"] == s, "log2_a_over_b"]
        q1, q3 = np.percentile(v, [25, 75])
        med = float(np.median(v))
        out["species"][s] = {
            "n": int(len(v)), "expected": exp, "median": med,
            "median_error": med - exp, "iqr": float(q3 - q1),
            "mad": float(np.median(np.abs(v - med))),
        }
    hum = out["species"]["human"]["median"]
    for s in out["species"]:
        out["species"][s]["median_error_human_centered"] = (
            out["species"][s]["median"] - hum - EXPECTED_LOG2_A_OVER_B[s])
    return out, d[["peptide", "species", "log2_a_over_b"]]


def main():
    results = {}
    ratio_frames = []
    for engine in ["upstream", "plus"]:
        runs = sorted((HERE / engine).glob("run[0-9]*"))
        timings = {}
        for r in runs:
            if (r / "time.txt").exists() and (r / "finished_at.txt").exists():
                timings[r.name] = parse_time(r / "time.txt")
        first = runs[0]
        psms = load_psms(engine, first)
        idm = id_metrics(psms)
        lfq, ratios = lfq_metrics(load_lfq(engine, first))
        ratios["engine"] = engine
        ratio_frames.append(ratios)
        walls = [t["wall_s"] for t in timings.values()]
        rss = [t["peak_rss_gib"] for t in timings.values()]
        # Identification counts must not depend on the run; check run-to-run identity.
        repeat = {}
        for r in runs[1:]:
            if (r / "finished_at.txt").exists():
                repeat[r.name] = id_metrics(load_psms(engine, r))["total"]
        results[engine] = {
            "timings": timings,
            "wall_s_median": float(np.median(walls)),
            "peak_rss_gib_median": float(np.median(rss)),
            "ids": idm, "ids_from": first.name, "repeat_run_totals": repeat,
            "lfq": lfq,
        }
    (HERE / "metrics.json").write_text(json.dumps(results, indent=2))
    pd.concat(ratio_frames).to_parquet(HERE / "lfq_ratios.parquet")
    print(json.dumps({e: {k: v for k, v in r.items() if k != "ids"} | {"ids_total": r["ids"]["total"]}
                      for e, r in results.items()}, indent=1))


if __name__ == "__main__":
    main()
