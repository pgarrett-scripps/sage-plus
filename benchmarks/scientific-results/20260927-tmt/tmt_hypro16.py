#!/usr/bin/env python3
"""HYpro16 (PXD020815) TMTpro-16 known-ratio benchmark: Sage Plus vs upstream Sage.

One command runs every stage and skips work that is already done:

    benchmarks/scientific-results/20260927-tmt/run.sh <sage-plus-binary> <label>

Stages: fetch RAW (PRIDE) -> convert mzML (ThermoRawFileParser, same build as the
20260914 LFQ pilot) -> write configs -> search (memgate-wrapped) -> analyze.

Search jobs (runs/<engine>-<method>-<input>):
  upstream-hrms2-mzml, upstream-sps-mzml                   (frozen upstream df92199)
  plus-<label>-hrms2-mzml, plus-<label>-sps-mzml           (candidate binary)
  plus-<label>-hrms2-raw,  plus-<label>-sps-raw            (candidate on Thermo RAW)

Missing reporter channels: before Beta 13 Sage Plus wrote 0.0, from Beta 13 null.
Upstream writes 0.0. The analysis treats 0, null and NaN identically as missing.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import re
import shutil
import subprocess
import sys
import time
from pathlib import Path

import numpy as np
import pyarrow.parquet as pq

HERE = Path(__file__).resolve().parent
DEFAULT_ROOT = Path("/mnt/data1/sage-plus-scientific/tmt-hypro16")
PILOT = Path("/mnt/data1/sage-plus-scientific/20260914")
TRFP = PILOT / "tools/thermo/ThermoRawFileParser"
MEMGATE = Path("/mnt/data1/explore-data/memgate.sh")
UPSTREAM_BIN = Path("/home/ty/Repos/sage-plus/benchmarks/.work/targets/baseline-df9219951cc9a54c/release/sage")
UPSTREAM_COMMIT = "df9219951cc9a54c (df92199), sage 0.15.0-beta.2; same binary as the 20260914 LFQ pilot"
PRIDE = "https://ftp.pride.ebi.ac.uk/pride/data/archive/2021/08/PXD020815"

# Deposit file sizes (PRIDE API, 2026-09-27) for download integrity.
FILES = {
    "hrms2": {"a18203": 277859801, "a18209": 262285211, "a18212": 255589225},
    "sps": {"a18207": 251612736, "a18210": 248201395, "a18213": 245163665},
}

CHANNELS = ["126", "127N", "127C", "128N", "128C", "129N", "129C", "130N",
            "130C", "131N", "131C", "132N", "132C", "133N", "133C", "134N"]
# Yeast amount relative to the 1x channels (133N, 133C, 134N). Frozen after checking
# Navarrete-Perea et al. JASMS 2021 32:247 (doi:10.1021/jasms.0c00299, PMC8210950):
# "a single empty channel, followed by triplicates channels representing a 20:1, 10:1,
# 4:1, and 2:1 decrease in the amount of yeast peptides from the final three channels";
# Methods list yeast as 0:(1/20)x3:(1/10)x3:(1/4)x3:(1/2)..:1:1:1 (the printed list drops
# one 1/2 entry; 16 channels = 1 empty + 5 triplicates). Human is 1:1 in all 16.
# Matches PXD020815_community_annotated.sdrf.tsv. See metadata/channel-map-source.md.
YEAST_LEVEL = [0.0] + [1 / 20] * 3 + [1 / 10] * 3 + [1 / 4] * 3 + [1 / 2] * 3 + [1.0] * 3
LEVELS = {"1/20": [1, 2, 3], "1/10": [4, 5, 6], "1/4": [7, 8, 9], "1/2": [10, 11, 12], "1": [13, 14, 15]}
REF = [13, 14, 15]
ZERO = 0
# Paper median observed yeast ratios (1x : level), Lumos + FAIMS, Sequest, S/N >= 10/channel.
# SPS-MS3 values are only stated for 2:1 and 20:1; 4:1 and 10:1 are given for RTS-MS3 only.
PAPER = {
    "hrms2": {"1/2": 1.81, "1/4": 3.16, "1/10": 5.67, "1/20": 7.69},
    "sps": {"1/2": 1.9, "1/4": None, "1/10": None, "1/20": 13.55},
    "rts_reference": {"1/2": 1.9, "1/4": 3.64, "1/10": 8.01, "1/20": 13.55},
}
Q = 0.01


def base_config(method: str, fasta: Path) -> dict:
    hr = method == "hrms2"
    # The research plan's config (tmt-known-ratio-candidates.md section 3), identical JSON for
    # both engines; everything else is left at each engine's own default.
    return {
        "database": {
            "enzyme": {"missed_cleavages": 2, "cleave_at": "KR", "restrict": "P", "min_len": 7, "max_len": 50},
            "static_mods": {"C": 57.021464, "K": 304.207146, "^": 304.207146},
            "variable_mods": {"M": [15.994915]},
            "max_variable_mods": 2,
            "decoy_tag": "rev_",
            "generate_decoys": True,
            "fasta": str(fasta),
        },
        "precursor_tol": {"ppm": [-10.0, 10.0]},
        # SPS-MS3 identification scans are ion-trap CID: wide fragment window, no deisotoping.
        "fragment_tol": {"ppm": [-20.0, 20.0]} if hr else {"da": [-0.5, 0.5]},
        "isotope_errors": [-1, 3],
        "deisotope": hr,
        "quant": {"tmt": "Tmt16", "tmt_settings": {"level": 2 if hr else 3, "sn": False}},
        "mzml_paths": [],
    }


# ---------------------------------------------------------------- utilities

def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 22), b""):
            h.update(chunk)
    return h.hexdigest()


def log(msg: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def version(binary: Path) -> str:
    out = subprocess.run([str(binary), "--version"], capture_output=True, text=True)
    return (out.stdout or out.stderr).strip()


# ---------------------------------------------------------------- stages

def fetch(root: Path) -> None:
    (root / "raw").mkdir(parents=True, exist_ok=True)
    (root / "metadata").mkdir(parents=True, exist_ok=True)
    for name in ["2020_Hyp.xlsx", "PXD020815_community_annotated.sdrf.tsv"]:
        if not (root / "metadata" / name).exists():
            subprocess.run(["curl", "-sSfL", "-o", str(root / "metadata" / name), f"{PRIDE}/{name}"], check=True)
    for files in FILES.values():
        for stem, size in files.items():
            path = root / "raw" / f"{stem}.raw"
            if path.exists() and path.stat().st_size == size:
                continue
            log(f"download {stem}.raw")
            subprocess.run(["curl", "-sSfL", "-C", "-", "--retry", "5", "-o", str(path), f"{PRIDE}/{stem}.raw"], check=True)
            if path.stat().st_size != size:
                sys.exit(f"{path}: size {path.stat().st_size} != deposit {size}")


def convert(root: Path) -> None:
    (root / "mzml").mkdir(exist_ok=True)
    for files in FILES.values():
        for stem in files:
            out = root / "mzml" / f"{stem}.mzML"
            if out.exists() and out.stat().st_size > 0:
                continue
            log(f"convert {stem}.raw")
            subprocess.run([str(TRFP), f"-i={root / 'raw' / (stem + '.raw')}", f"-b={out}", "-f=1", "-m=0"],
                           check=True, stdout=subprocess.DEVNULL)


def jobs(root: Path, plus_label: str, plus_bin: Path, upstream_bin: Path) -> list[dict]:
    out = []
    for method in FILES:
        out.append(dict(name=f"upstream-{method}-mzml", engine="upstream", method=method, input="mzml", binary=upstream_bin))
        for inp in ("mzml", "raw"):
            out.append(dict(name=f"plus-{plus_label}-{method}-{inp}", engine="plus", method=method, input=inp, binary=plus_bin))
    for j in out:
        ext = ".mzML" if j["input"] == "mzml" else ".raw"
        j["spectra"] = [str(root / j["input"] / f"{s}{ext}") for s in FILES[j["method"]]]
        j["dir"] = root / "runs" / j["name"]
    return out


def search(job: dict, fasta: Path, gb: int, force: bool) -> None:
    d: Path = job["dir"]
    binary: Path = job["binary"]
    stamp = d / "job.json"
    bin_sha = sha256(binary)
    if not force and stamp.exists() and (d / "results.sage.parquet").exists():
        prev = json.loads(stamp.read_text())
        if prev.get("binary_sha256") == bin_sha and prev.get("returncode") == 0:
            log(f"skip {job['name']} (done with same binary)")
            return
    if d.exists():
        shutil.rmtree(d)
    d.mkdir(parents=True)
    cfg = base_config(job["method"], fasta)
    cfg["mzml_paths"] = job["spectra"]
    cfg["output_directory"] = str(d)
    (d / "config.json").write_text(json.dumps(cfg, indent=2))
    cmd = [str(MEMGATE), str(gb), "/usr/bin/time", "-v", str(binary), str(d / "config.json"),
           "--output_directory", str(d), "--batch-size", "3", "--disable-telemetry-i-dont-want-to-improve-sage"]
    if job["engine"] == "upstream":
        cmd.append("--parquet")  # Sage Plus always writes parquet and has no --parquet flag
    log(f"search {job['name']}")
    t0 = time.time()
    with open(d / "search.log", "w") as fh:
        rc = subprocess.run(cmd, stdout=fh, stderr=subprocess.STDOUT).returncode
    text = (d / "search.log").read_text(errors="replace")
    m = re.search(r"Maximum resident set size \(kbytes\): (\d+)", text)
    stamp.write_text(json.dumps({
        "name": job["name"], "engine": job["engine"], "method": job["method"], "input": job["input"],
        "binary": str(binary), "binary_sha256": bin_sha, "binary_version": version(binary),
        "command": cmd, "returncode": rc, "wall_seconds": round(time.time() - t0, 1),
        "max_rss_gib": round(int(m.group(1)) / 2**20, 2) if m else None,
        "config_sha256": sha256(d / "config.json"),
    }, indent=2))
    if rc != 0:
        sys.exit(f"{job['name']} failed (rc={rc}); see {d / 'search.log'}")


# ---------------------------------------------------------------- analysis

def load_psms(path: Path) -> dict:
    t = pq.read_table(path)
    cols = set(t.column_names)
    keep = [c for c in ["filename", "scannr", "peptide", "stripped_peptide", "proteins", "is_decoy", "rank",
                        "spectrum_q", "peptide_q", "protein_q", "reporter_ion_intensity", "charge"] if c in cols]
    d = t.select(keep).to_pydict()
    n = len(d["peptide"])
    ok = [(not d["is_decoy"][i]) and d["spectrum_q"][i] <= Q and d["peptide_q"][i] <= Q
          and d.get("rank", [1] * n)[i] == 1 for i in range(n)]
    rows = {k: [v[i] for i in range(n) if ok[i]] for k, v in d.items()}
    rep = np.full((len(rows["peptide"]), 16), np.nan)
    for i, r in enumerate(rows["reporter_ion_intensity"]):
        if r is None:
            continue
        v = np.array([np.nan if x is None else x for x in r], dtype=float)
        v[~np.isfinite(v) | (v <= 0)] = np.nan  # 0, null, NaN, inf => missing
        rep[i, : len(v)] = v[:16]
    rows["rep"] = rep
    rows["has_reporter"] = np.array([r is not None for r in rows["reporter_ion_intensity"]])
    del rows["reporter_ion_intensity"]
    if "stripped_peptide" not in rows:
        rows["stripped_peptide"] = [re.sub(r"[^A-Z]", "", re.sub(r"\[[^\]]*\]", "", p)) for p in rows["peptide"]]
    rows["file"] = [Path(f).name.split(".")[0] for f in rows["filename"]]
    rows["scan"] = [int(re.findall(r"(\d+)", str(s))[-1]) for s in rows["scannr"]]
    rows["protein_list"] = [p if isinstance(p, list) else str(p).split(";") for p in rows["proteins"]]
    return rows


class Species:
    """Species from accession suffix, then an I/L-normalized substring check against the
    other proteome so that peptides shared across species (as in the LFQ pilot) are excluded."""

    def __init__(self, fasta_dir: Path):
        self.seq = {}
        for sp, name in [("human", "human.fasta"), ("yeast", "yeast.fasta")]:
            parts, cur = [], []
            for line in open(PILOT / "references" / name):
                if line.startswith(">"):
                    if cur:
                        parts.append("".join(cur))
                    cur = []
                else:
                    cur.append(line.strip())
            parts.append("".join(cur))
            self.seq[sp] = "|".join(parts).replace("I", "L")
        self.cache = {}

    def of(self, proteins: list[str], stripped: str) -> str | None:
        tags = set()
        for p in proteins:
            p = p.replace("rev_", "")
            tags.add("human" if p.endswith("_HUMAN") else "yeast" if p.endswith("_YEAST") else "other")
        if len(tags) != 1 or "other" in tags:
            return None
        sp = tags.pop()
        key = stripped.replace("I", "L")
        if key not in self.cache:
            self.cache[key] = {s for s, seq in self.seq.items() if key in seq}
        return sp if self.cache[key] == {sp} else None


def nanmedian(x):
    x = np.asarray(x, dtype=float)
    x = x[np.isfinite(x)]
    return float(np.median(x)) if len(x) else None


def r(x, nd=4):
    return None if x is None or (isinstance(x, float) and not math.isfinite(x)) else round(float(x), nd)


def ratio_levels(mask, l2, stripped, method):
    """Per-peptide median log2(channel / 1x mean) over PSMs in mask, pooled per ratio level."""
    peptides = {}
    for i in np.where(mask)[0]:
        peptides.setdefault(stripped[i], []).append(i)
    pep_l2 = np.array([np.nanmedian(l2[idx], axis=0) for idx in peptides.values()]) if peptides else np.empty((0, 16))
    levels, xs, ys = {}, [], []
    for name, chans in LEVELS.items():
        expected = math.log2(YEAST_LEVEL[chans[0]])
        vals = pep_l2[:, chans].ravel()
        vals = vals[np.isfinite(vals)]
        err = vals - expected
        med = nanmedian(vals)
        levels[name] = {
            "expected_log2": r(expected), "n_values": int(len(vals)),
            "median_log2": r(med), "median_bias": r(nanmedian(err)),
            "mae": r(float(np.mean(np.abs(err)))) if len(err) else None,
            "median_abs_error": r(nanmedian(np.abs(err))),
            "observed_fold_1x_over_level": r(2 ** -med, 3) if med is not None and name != "1" else None,
            "paper_fold": PAPER[method].get(name) if name != "1" else None,
        }
        xs += [expected] * len(vals)
        ys += list(vals)
    slope = intercept = None
    if len(xs) > 2:
        slope, intercept = np.polyfit(np.array(xs), np.array(ys), 1)
    # Slope on the four spiked levels only (1x channels are part of the reference).
    xs4 = [x for x in xs if x != 0.0]
    ys4 = [y for x, y in zip(xs, ys) if x != 0.0]
    slope4 = float(np.polyfit(xs4, ys4, 1)[0]) if len(xs4) > 2 else None
    return levels, slope, intercept, slope4, peptides


def analyze_run(path: Path, species: Species, method: str) -> dict:
    rows = load_psms(path)
    sp = np.array([species.of(p, s) for p, s in zip(rows["protein_list"], rows["stripped_peptide"])], dtype=object)
    rep = rows["rep"]
    files = np.array(rows["file"])
    observed_any = np.isfinite(rep).any(axis=1)
    quant = rows["has_reporter"] & observed_any
    human = quant & (sp == "human")
    yeast = quant & (sp == "yeast")

    # Per-file channel normalization to equal human medians (human is 1:1 by design).
    norm = rep.copy()
    for f in np.unique(files):
        m = human & (files == f)
        med = np.array([nanmedian(rep[m, c]) or np.nan for c in range(16)])
        norm[files == f] = rep[files == f] / (med / np.nanmean(med))

    # Missing channels among quantified PSMs, channel 126 excluded from the headline rate.
    miss = ~np.isfinite(rep)
    missing = {
        "all_species_excl_126": r(miss[quant][:, 1:].mean()) if quant.any() else None,
        "human_excl_126": r(miss[human][:, 1:].mean()) if human.any() else None,
        "yeast_excl_126": r(miss[yeast][:, 1:].mean()) if yeast.any() else None,
        "yeast_126_missing": r(miss[yeast][:, 0].mean()) if yeast.any() else None,
        "per_channel_all": {CHANNELS[c]: r(miss[quant][:, c].mean()) for c in range(16)} if quant.any() else {},
        "psms_with_any_missing_excl_126": r(miss[quant][:, 1:].any(axis=1).mean()) if quant.any() else None,
    }

    # Yeast log2(channel / mean of 1x channels); per peptide median over PSMs.
    ref = np.nanmean(norm[:, REF], axis=1)
    ref_ok = np.isfinite(norm[:, REF]).all(axis=1)
    ym = yeast & ref_ok
    with np.errstate(divide="ignore", invalid="ignore"):
        l2 = np.log2(norm / ref[:, None])
    levels, slope, intercept, slope4, peptides = ratio_levels(ym, l2, rows["stripped_peptide"], method)
    # Sensitivity view closer to the paper's S/N >= 10 per channel / summed S/N >= 200 filter
    # (not reproducible from intensities): all 15 non-126 channels observed and summed reporter
    # intensity in the top half of those yeast PSMs.
    total = np.nansum(rep, axis=1)
    full = ym & np.isfinite(rep[:, 1:]).all(axis=1)
    cut = float(np.median(total[full])) if full.any() else np.inf
    hi = full & (total >= cut)
    lv_hi, _, _, slope4_hi, pep_hi = ratio_levels(hi, l2, rows["stripped_peptide"], method)
    sensitivity = {"filter": "yeast PSMs with 127N-134N all observed and summed intensity >= median",
                   "psms": int(hi.sum()), "peptides": len(pep_hi), "compression_slope_spiked_levels": r(slope4_hi),
                   "levels": {k: {x: v[x] for x in ["median_bias", "mae", "observed_fold_1x_over_level", "paper_fold"]}
                              for k, v in lv_hi.items()}}

    # Leakage: yeast signal in the zero channel 126 relative to the 1x mean.
    leak = norm[ym, 0] / ref[ym]
    leak_obs = leak[np.isfinite(leak)]
    leak_all = np.where(np.isfinite(leak), leak, 0.0)
    leakage = {
        "yeast_psms": int(ym.sum()),
        "fraction_126_observed": r(np.isfinite(leak).mean()) if ym.any() else None,
        "median_126_over_1x_missing_as_0": r(nanmedian(leak_all)),
        "mean_126_over_1x_missing_as_0": r(float(leak_all.mean())) if len(leak_all) else None,
        "median_126_over_1x_observed_only": r(nanmedian(leak_obs)),
        "median_126_over_1_20_level_observed_only": r(nanmedian(
            norm[ym, 0][np.isfinite(leak)] / np.nanmean(norm[ym][:, LEVELS["1/20"]], axis=1)[np.isfinite(leak)])),
    }

    # Human precision: CV across 16 normalized channels, PSMs with all channels observed.
    hc = norm[human & np.isfinite(norm).all(axis=1)]
    cv = np.std(hc, axis=1, ddof=1) / np.mean(hc, axis=1) if len(hc) else np.array([])
    hpep = {}
    for i in np.where(human & np.isfinite(norm).all(axis=1))[0]:
        hpep.setdefault(rows["stripped_peptide"][i], []).append(i)
    pcv = [np.std(v, ddof=1) / np.mean(v) for v in (np.nanmedian(norm[idx], axis=0) for idx in hpep.values())]

    def counts(mask):
        prot = {";".join(sorted(rows["protein_list"][i])) for i in np.where(mask)[0]
                if rows.get("protein_q", [0] * len(mask))[i] <= Q}
        return {"psms": int(mask.sum()), "peptides": len({rows["stripped_peptide"][i] for i in np.where(mask)[0]}),
                "proteins_q01": len(prot)}

    per_file = {f: {"quantified_psms": int((quant & (files == f)).sum()),
                    "yeast_psms": int((yeast & (files == f)).sum()),
                    "human_psms": int((human & (files == f)).sum())} for f in np.unique(files)}
    return {
        "accepted_psms_q01": int(len(rows["peptide"])),
        "quantified": {"all": counts(quant), "human": counts(human), "yeast": counts(yeast)},
        "cross_or_unassigned_species_excluded": int((quant & (sp == None)).sum()),  # noqa: E711
        "levels": levels,
        "compression_slope_all_levels": r(slope), "intercept": r(intercept),
        "compression_slope_spiked_levels": r(slope4),
        "yeast_peptides_in_ratio_analysis": int(len(peptides)),
        "sensitivity_high_signal": sensitivity,
        "leakage_126": leakage,
        "missing_channels": missing,
        "human_cv": {"psms_all16": int(len(hc)), "median_psm_cv": r(nanmedian(cv)),
                     "median_peptide_cv": r(nanmedian(pcv)), "peptides": len(hpep)},
        "per_file": per_file,
        "_keys": {(rows["file"][i], rows["scan"][i]): (rows["stripped_peptide"][i], rep[i]) for i in np.where(quant)[0]},
    }


def parity(a: dict, b: dict) -> dict:
    """Sage Plus RAW vs mzML: shared quantified spectra, same-peptide rate, reporter agreement."""
    ka, kb = a["_keys"], b["_keys"]
    shared = set(ka) & set(kb)
    same = [k for k in shared if ka[k][0] == kb[k][0]]
    diffs, both_missing, one_missing, n = [], 0, 0, 0
    for k in same:
        x, y = ka[k][1], kb[k][1]
        fx, fy = np.isfinite(x), np.isfinite(y)
        both_missing += int((~fx & ~fy).sum())
        one_missing += int((fx ^ fy).sum())
        n += 16
        m = fx & fy
        diffs.extend(np.abs(np.log2(x[m] / y[m])))
    diffs = np.array(diffs)
    return {"raw_quantified_psms": len(ka), "mzml_quantified_psms": len(kb), "shared_spectra": len(shared),
            "same_peptide": len(same), "channel_values_compared": n,
            "median_abs_log2_diff": r(nanmedian(diffs), 6),
            "p99_abs_log2_diff": r(float(np.quantile(diffs, 0.99)), 6) if len(diffs) else None,
            "fraction_exact_equal": r(float((diffs < 1e-6).mean()), 4) if len(diffs) else None,
            "channel_missing_in_one_only": one_missing, "channel_missing_in_both": both_missing}


def analyze(root: Path, all_jobs: list[dict], label: str) -> dict:
    species = Species(root / "fasta")
    out = {"label": label, "runs": {}, "parity_raw_vs_mzml": {}, "paired_upstream_vs_plus_per_file": {}}
    res = {}
    for j in all_jobs:
        p = j["dir"] / "results.sage.parquet"
        if not p.exists():
            continue
        log(f"analyze {j['name']}")
        res[j["name"]] = analyze_run(p, species, j["method"])
        stamp = json.loads((j["dir"] / "job.json").read_text())
        res[j["name"]]["job"] = {k: stamp[k] for k in ["engine", "method", "input", "binary", "binary_sha256",
                                                       "binary_version", "wall_seconds", "max_rss_gib", "config_sha256"]}
    for method in FILES:
        a, b = res.get(f"plus-{label}-{method}-raw"), res.get(f"plus-{label}-{method}-mzml")
        if a and b:
            out["parity_raw_vs_mzml"][method] = parity(a, b)
        u, pl = res.get(f"upstream-{method}-mzml"), res.get(f"plus-{label}-{method}-mzml")
        if u and pl:
            out["paired_upstream_vs_plus_per_file"][method] = {
                f: {"upstream": u["per_file"].get(f), "plus": pl["per_file"].get(f)} for f in FILES[method]}
    for k, v in res.items():
        v.pop("_keys", None)
        out["runs"][k] = v
    return out


def markdown(summary: dict) -> str:
    lines = ["| Run | Quant PSMs (H/Y) | Yeast pep | Bias 1/2 | 1/4 | 1/10 | 1/20 | MAE 1/20 | Fold 20:1 (paper) | Slope | 126 leak (obs med) | Missing excl 126 | Human CV |",
             "|---|---|---|---|---|---|---|---|---|---|---|---|---|"]
    for name, v in summary["runs"].items():
        L = v["levels"]
        paper = L["1/20"]["paper_fold"]
        lines.append(
            f"| {name} | {v['quantified']['human']['psms']}/{v['quantified']['yeast']['psms']} | "
            f"{v['quantified']['yeast']['peptides']} | "
            + " | ".join(str(L[k]["median_bias"]) for k in ["1/2", "1/4", "1/10", "1/20"])
            + f" | {L['1/20']['mae']} | {L['1/20']['observed_fold_1x_over_level']} ({paper}) | "
            f"{v['compression_slope_spiked_levels']} | {v['leakage_126']['median_126_over_1x_observed_only']} | "
            f"{v['missing_channels']['all_species_excl_126']} | {v['human_cv']['median_psm_cv']} |")
    lines += ["", "High-signal sensitivity (all 127N-134N observed, summed intensity >= median): observed fold 1x:level (paper)", "",
              "| Run | PSMs | 2:1 | 4:1 | 10:1 | 20:1 | Slope |", "|---|---|---|---|---|---|---|"]
    for name, v in summary["runs"].items():
        h = v["sensitivity_high_signal"]
        lines.append(f"| {name} | {h['psms']} | " + " | ".join(
            f"{h['levels'][k]['observed_fold_1x_over_level']} ({h['levels'][k]['paper_fold']})" for k in ["1/2", "1/4", "1/10", "1/20"])
            + f" | {h['compression_slope_spiked_levels']} |")
    return "\n".join(lines) + "\n"


def provenance(root: Path, all_jobs: list[dict]) -> dict:
    prov_path = root / "metadata/input-sha256.json"
    prov = json.loads(prov_path.read_text()) if prov_path.exists() else {}
    for sub, ext in [("raw", ".raw"), ("mzml", ".mzML")]:
        for files in FILES.values():
            for stem in files:
                p = root / sub / f"{stem}{ext}"
                key = f"{sub}/{stem}{ext}"
                if p.exists() and (key not in prov or prov[key]["bytes"] != p.stat().st_size):
                    prov[key] = {"bytes": p.stat().st_size, "sha256": sha256(p)}
    for name in ["human.fasta", "yeast.fasta"]:
        p = PILOT / "references" / name
        prov[f"fasta/{name}"] = {"bytes": p.stat().st_size, "sha256": sha256(p),
                                 "receipt": str(PILOT / "references" / name.replace(".fasta", ".receipt.json"))}
    p = root / "fasta/human-yeast.fasta"
    prov["fasta/human-yeast.fasta"] = {"bytes": p.stat().st_size, "sha256": sha256(p), "made_by": "cat human.fasta yeast.fasta"}
    prov_path.write_text(json.dumps(prov, indent=2, sort_keys=True))
    return {"inputs": prov, "thermorawfileparser": str(TRFP), "trfp_zip_sha256": sha256(PILOT / "tools/ThermoRawFileParser-v.2.0.0-dev-linux.zip"),
            "upstream": UPSTREAM_COMMIT, "channel_map_source": "metadata/channel-map-source.md"}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--plus-bin", type=Path, required=True)
    ap.add_argument("--plus-label", required=True, help="e.g. b12, b13; used in run directory names")
    ap.add_argument("--upstream-bin", type=Path, default=UPSTREAM_BIN)
    ap.add_argument("--root", type=Path, default=DEFAULT_ROOT)
    ap.add_argument("--memgate-gb", type=int, default=12)
    ap.add_argument("--force", action="store_true", help="rerun searches even if outputs exist")
    ap.add_argument("--analyze-only", action="store_true")
    ap.add_argument("--no-raw", action="store_true", help="skip the Sage Plus RAW-input runs")
    a = ap.parse_args()
    root = a.root
    fasta = root / "fasta/human-yeast.fasta"
    if not fasta.exists():
        fasta.parent.mkdir(parents=True, exist_ok=True)
        with open(fasta, "wb") as out:
            for name in ["human.fasta", "yeast.fasta"]:
                out.write((PILOT / "references" / name).read_bytes())
    all_jobs = jobs(root, a.plus_label, a.plus_bin.resolve(), a.upstream_bin.resolve())
    if a.no_raw:
        all_jobs = [j for j in all_jobs if j["input"] != "raw"]
    if not a.analyze_only:
        fetch(root)
        convert(root)
        for j in all_jobs:
            search(j, fasta, a.memgate_gb, a.force)
    summary = analyze(root, all_jobs, a.plus_label)
    summary["provenance"] = provenance(root, all_jobs)
    summary["settings"] = {"q_value": Q, "channels": CHANNELS, "yeast_level": YEAST_LEVEL, "paper": PAPER,
                           "config_hrms2": base_config("hrms2", fasta), "config_sps": base_config("sps", fasta)}
    outdir = root / "analysis" / a.plus_label
    outdir.mkdir(parents=True, exist_ok=True)
    (outdir / "summary.json").write_text(json.dumps(summary, indent=2))
    (outdir / "table.md").write_text(markdown(summary))
    shutil.copy(outdir / "summary.json", HERE / f"summary-{a.plus_label}.json")
    shutil.copy(outdir / "table.md", HERE / f"table-{a.plus_label}.md")
    print(markdown(summary))


if __name__ == "__main__":
    main()
