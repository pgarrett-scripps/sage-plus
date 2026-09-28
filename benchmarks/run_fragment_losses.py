#!/usr/bin/env python3
"""Generic fragment-loss (water, ammonia) benchmark.

Searches each dataset with a baseline executable (no fragment losses) and a
candidate executable in several integration variants, then reports PSMs,
peptides and proteins at 1% q-value, decoys at 1%, species-entrapment FDP
(for human samples searched against human + yeast + E. coli), wall time and
peak RSS.

Variants (see FRAGMENT_LOSSES.md):
  base  baseline executable, no `fragment_losses` key
  A     candidate, no `fragment_losses` key (must equal base byte for byte)
  B     candidate, losses scored inside the hyperscore
  C     candidate, losses as separate LDA features (the shipped behaviour)
  D     candidate, as C but chimeric peak removal keeps loss-only peaks
The `-chimera` suffix runs a variant with `chimera: true`.

B and D are selected with the experiment hook SAGE_PLUS_FRAGMENT_LOSS_SCORING,
which exists only in the experiment commit named in FRAGMENT_LOSSES.md.

usage: run_fragment_losses.py --base BIN --candidate BIN --output DIR
       [--datasets NAME ...] [--variants NAME ...]
Requires the `duckdb` Python package.
"""

from __future__ import annotations

import argparse
import copy
import json
import os
import re
import subprocess
import sys
from pathlib import Path

import duckdb

sys.path.insert(0, str(Path(__file__).parent))
from picked_fdr_entrapment import is_target_list, ratios  # noqa: E402
from provenance import sha256  # noqa: E402

TIME = "/usr/bin/time"
MEMGATE = "/mnt/data1/explore-data/memgate.sh"
SCI = "/mnt/data1/sage-plus-scientific/20260914"
HYE = f"{SCI}/references/hye-irt-defined.fasta"
PROFILE = "/mnt/data1/scan-profile-bench"
REPO_DATA = "/home/ty/Repos/sage-plus/data/silac-k6r6"

# Configuration approved for Beta 14.
LOSSES = {
    "fragment_losses": {
        "Water": {"mass": 18.010565, "sites": ["S", "T", "E", "D"],
                  "ion_kinds": ["b", "y"], "allow_modified": False},
        "Ammonia": {"mass": 17.026549, "sites": ["R", "K", "N", "Q"], "ion_kinds": ["y"]},
    },
    "max_fragment_losses": 1,
}

TEMPLATE = {
    "database": {
        "bucket_size": 16384,
        "enzyme": {"missed_cleavages": 1, "min_len": 7, "max_len": 50,
                   "cleave_at": "KR", "restrict": "P"},
        "peptide_min_mass": 500.0,
        "peptide_max_mass": 5000.0,
        "static_mods": {"C": 57.021464},
        "variable_mods": {"M": [15.994915]},
        "max_variable_mods": 2,
        "decoy_tag": "rev_",
        "generate_decoys": True,
        "fasta": HYE,
    },
    "precursor_tol": {"ppm": [-10.0, 10.0]},
    "fragment_tol": {"ppm": [-20.0, 20.0]},
    "isotope_errors": [-1, 3],
    "min_matched_peaks": 6,
    "report_psms": 1,
    "output_filter": {"psm_q_value": 1.0},
    "predict_rt": True,
    "protein_grouping": True,
    "max_memory_gb": 10.0,
    "batch_size": 1,
}

# name -> (description, entrapment, overrides)
DATASETS = {
    "hek-silac": ("HEK SILAC K6R6, Orbitrap HCD, human FASTA", False, {
        "database": {"fasta": f"{REPO_DATA}/human-reviewed.fasta"},
        "mzml_paths": [f"{REPO_DATA}/HEK_SILAC-K6R6.mzML"],
    }),
    "pxd028735-human01": ("PXD028735 LFQ human-only, Orbitrap HCD", True, {
        "mzml_paths": [f"{SCI}/converted/PXD028735/LFQ_Orbitrap_DDA_Human_01.mzML"],
    }),
    "pxd001468-hek293t": ("PXD001468 HEK293T, 2 files, Orbitrap HCD", True, {
        "mzml_paths": [
            f"{SCI}/converted/PXD001468/b1906_293T_proteinID_01A_QE3_122212.mgf",
            f"{SCI}/converted/PXD001468/b1937_293T_proteinID_01B_QE3_122212.mgf",
        ],
    }),
    "pxd011070-itcid": ("PXD011070 HeLa, ion-trap CID (0.5 Da)", True, {
        "fragment_tol": {"da": [-0.5, 0.5]},
        "deisotope": False,
        "mzml_paths": [f"{PROFILE}/PXD011070/mzml/ch_23Aug2018_HeLa_Std_1.mzML"],
    }),
    "pxd004447-etcid": ("PXD004447 Jurkat, Orbitrap ETciD, b/c/y/z-dot", True, {
        "database": {"ion_kinds": ["b", "c", "y", "z_dot"]},
        "mzml_paths": [
            f"{PROFILE}/PXD004447/mzml/F1_Agilent5_20150701_LT_JurTryp_ETciD30_2_2.mzML"],
    }),
}

VARIANTS = {
    # name -> (engine, losses, scoring hook, chimera)
    "base": ("base", False, None, False),
    "A": ("candidate", False, None, False),
    "B": ("candidate", True, "hyperscore", False),
    "C": ("candidate", True, "features", False),
    "base-chimera": ("base", False, None, True),
    "C-chimera": ("candidate", True, "features", True),
    "D-chimera": ("candidate", True, "features_keep_peaks", True),
}


def merge(dst: dict, src: dict) -> dict:
    for key, value in src.items():
        if isinstance(value, dict) and isinstance(dst.get(key), dict):
            merge(dst[key], value)
        else:
            dst[key] = copy.deepcopy(value)
    return dst


def config(dataset: str, variant: str) -> dict:
    _, losses, _, chimera = VARIANTS[variant]
    cfg = merge(copy.deepcopy(TEMPLATE), DATASETS[dataset][2])
    if losses:
        merge(cfg["database"], LOSSES)
    if chimera:
        cfg["chimera"] = True
    return cfg


def run(binary: str, cfg_path: Path, out: Path, hook: str | None) -> dict:
    env = dict(os.environ)
    env.pop("SAGE_PLUS_FRAGMENT_LOSS_SCORING", None)
    if hook:
        env["SAGE_PLUS_FRAGMENT_LOSS_SCORING"] = hook
    # Time inside the memory gate so queueing for memory is not counted.
    cmd = [MEMGATE, "12", TIME, "-v", binary, str(cfg_path), "-o", str(out), "--overwrite",
           "--disable-telemetry-i-dont-want-to-improve-sage"]
    proc = subprocess.run(cmd, env=env, capture_output=True, text=True)
    (out.parent / f"{out.name}.log").write_text(proc.stderr)
    if proc.returncode:
        raise SystemExit(f"{cfg_path}: exit {proc.returncode}\n{proc.stderr[-2000:]}")
    wall = re.search(r"Elapsed \(wall clock\) time.*: (\S+)", proc.stderr).group(1)
    seconds = sum(float(x) * 60**i for i, x in enumerate(reversed(wall.split(":"))))
    rss = int(re.search(r"Maximum resident set size \(kbytes\): (\d+)", proc.stderr).group(1))
    return {"wall_s": round(seconds, 1), "peak_rss_gb": round(rss / 1024**2, 2)}


def metrics(parquet: Path, fasta: str, entrapment: bool, r: dict | None) -> dict:
    con = duckdb.connect()
    con.execute(f"create view t as select * from read_parquet('{parquet}')")
    q = 0.01
    out = {
        "psms": con.execute(
            f"select count(*) from t where not is_decoy and spectrum_q <= {q}").fetchone()[0],
        "decoy_psms": con.execute(
            f"select count(*) from t where is_decoy and spectrum_q <= {q}").fetchone()[0],
        "peptides": con.execute(
            f"select count(distinct peptide) from t where not is_decoy and peptide_q <= {q}"
        ).fetchone()[0],
        "proteins": con.execute(
            f"select count(distinct proteins) from t where not is_decoy and num_proteins = 1 "
            f"and protein_q <= {q}").fetchone()[0],
    }
    if entrapment:
        def fdp(rows, ratio):
            n = len(rows)
            n_e = n - sum(is_target_list(p) for p in rows)
            return round(n_e * (1 + 1 / ratio) / n, 4) if n else None
        psm = [p for (p,) in con.execute(
            f"select proteins from t where not is_decoy and spectrum_q <= {q}").fetchall()]
        pep = [p for _, p in con.execute(
            f"select distinct peptide, proteins from t where not is_decoy and peptide_q <= {q}"
        ).fetchall()]
        prot = [p for (p,) in con.execute(
            f"select distinct proteins from t where not is_decoy and num_proteins = 1 "
            f"and protein_q <= {q}").fetchall()]
        out["fdp_psm"] = fdp(psm, r["peptide_r"])
        out["fdp_peptide"] = fdp(pep, r["peptide_r"])
        out["fdp_protein"] = fdp(prot, r["protein_r"])
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--base", required=True)
    ap.add_argument("--candidate", required=True)
    ap.add_argument("--output", required=True, type=Path)
    ap.add_argument("--datasets", nargs="+", default=list(DATASETS))
    ap.add_argument("--variants", nargs="+", default=["base", "A", "B", "C"])
    args = ap.parse_args()
    engines = {"base": args.base, "candidate": args.candidate}
    args.output.mkdir(parents=True, exist_ok=True)
    summary_path = args.output / "summary.json"
    summary = json.loads(summary_path.read_text()) if summary_path.exists() else {}
    summary["engines"] = {k: {"path": v, "sha256": sha256(Path(v))} for k, v in engines.items()}
    ratio_cache: dict[str, dict] = {}
    for dataset in args.datasets:
        _, entrapment, _ = DATASETS[dataset]
        for variant in args.variants:
            engine, _, hook, _ = VARIANTS[variant]
            cfg = config(dataset, variant)
            name = f"{dataset}.{variant}"
            cfg_path = args.output / f"{name}.json"
            cfg_path.write_text(json.dumps(cfg, indent=1))
            out = args.output / name
            print(f"running {name}", flush=True)
            row = run(engines[engine], cfg_path, out, hook)
            parquet = out / "results.sage.parquet"
            fasta = cfg["database"]["fasta"]
            if entrapment and fasta not in ratio_cache:
                ratio_cache[fasta] = ratios(fasta)
            row.update(metrics(parquet, fasta, entrapment, ratio_cache.get(fasta)))
            row["results_sha256"] = sha256(parquet)
            summary.setdefault("runs", {})[name] = row
            summary_path.write_text(json.dumps(summary, indent=1))
            print(f"  {row}", flush=True)


if __name__ == "__main__":
    main()
