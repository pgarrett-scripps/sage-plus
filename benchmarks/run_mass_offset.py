#!/usr/bin/env python3
"""Mass offset evaluation matrix with recorded identities.

Compares searching a modification as a search-time mass offset against the
equivalent database expansion: cost and index size as offsets are added,
agreement of accepted identifications, localization against synthesis-defined
sites, and entrapment error at the peptide level.

Every job records the executable, configuration, and input identities, the
command, process limits, timings, and output identities.
"""

from __future__ import annotations

import argparse
import copy
import json
import platform
import re
import shutil
import subprocess
import time
from pathlib import Path

from provenance import atomic_json, file_identities, sha256

REPO = Path(__file__).resolve().parents[1]


def offset(mass: float, name: str) -> dict:
    return {"mass": mass, "name": name, "max_count": 1, "search_mode": "mass_offset"}


def indexed(mass: float, name: str) -> dict:
    return {"mass": mass, "name": name, "max_count": 1}


OXIDATION, PHOSPHO, DEAMIDATION = 15.994915, 79.966331, 0.984016


def scale_jobs(spectra: Path, fasta: Path, repeats: int) -> list[dict]:
    """Cost against the number of offsets, which are never combined."""
    acetyl = [indexed(42.010565, "Acetyl")]
    modes = {
        "indexed": {"^": acetyl, "M": [indexed(OXIDATION, "Oxidation")]},
        "offsets-1": {"^": acetyl, "M": [offset(OXIDATION, "Oxidation")]},
        "offsets-2": {"^": acetyl, "M": [offset(OXIDATION, "Oxidation")],
                      **{site: [offset(PHOSPHO, "Phospho")] for site in "STY"}},
        "offsets-3": {"^": acetyl, "M": [offset(OXIDATION, "Oxidation")],
                      **{site: [offset(PHOSPHO, "Phospho")] for site in "STY"},
                      **{site: [offset(DEAMIDATION, "Deamidation")] for site in "NQ"}},
    }
    jobs = []
    for name, mods in modes.items():
        for repeat in range(1, repeats + 1):
            config = {
                "database": {
                    "bucket_size": 16384,
                    "enzyme": {"missed_cleavages": 1, "cleave_at": "KR", "restrict": "P",
                               "min_len": 7, "max_len": 50},
                    "static_mods": {"C": {"mass": 57.021464, "name": "Carbamidomethyl"}},
                    "variable_mods": mods,
                    "max_variable_mods": 1,
                    "max_total_variable_mods": 1,
                    "decoy_tag": "rev_", "generate_decoys": True, "fasta": str(fasta),
                },
                "precursor_tol": {"ppm": [-10.0, 10.0]},
                "fragment_tol": {"ppm": [-20.0, 20.0]},
                "isotope_errors": [-1, 3], "deisotope": True, "chimera": False,
                "max_fragment_charge": 2, "min_matched_peaks": 6, "report_psms": 1,
                "output_filter": {"psm_q_value": 0.05},
                "max_memory_gb": 22.0, "min_free_memory_gb": 2.0, "batch_size": 1,
                "mzml_paths": [str(spectra)], "score_type": "SageHyperScore",
            }
            jobs.append({"suite": "scale", "id": f"{name}-r{repeat}", "config": config})
    return jobs


def localization_jobs(fasta: Path, files: list[Path]) -> list[dict]:
    """Synthesis-defined phosphosites, expansion against one offset."""
    jobs = []
    for spectra in files:
        for name, mods in (("indexed", indexed(PHOSPHO, "Phospho")),
                           ("offset", offset(PHOSPHO, "Phospho"))):
            config = {
                "database": {
                    "bucket_size": 16384,
                    "enzyme": {"missed_cleavages": 0, "cleave_at": "$", "restrict": "",
                               "min_len": 7, "max_len": 50},
                    "static_mods": {"C": 57.021464},
                    "variable_mods": {site: [dict(mods)] for site in "STY"},
                    "max_variable_mods": 1, "max_total_variable_mods": 1,
                    "max_combinations": 128,
                    "decoy_tag": "rev_", "generate_decoys": True, "fasta": str(fasta),
                },
                "precursor_tol": {"ppm": [-10.0, 10.0]},
                "fragment_tol": {"ppm": [-20.0, 20.0]},
                "isotope_errors": [-1, 3], "deisotope": True, "chimera": False,
                "max_fragment_charge": 2, "min_matched_peaks": 6, "report_psms": 1,
                "output_filter": {"psm_q_value": 1.0},
                "ptm_localization": {"enabled": True, "psm_q_value": 1.0,
                                     "localization_q_value": 1.0},
                "max_memory_gb": 22.0, "min_free_memory_gb": 2.0, "batch_size": 1,
                "mzml_paths": [str(spectra)], "score_type": "SageHyperScore",
            }
            jobs.append({"suite": "localization",
                         "id": f"{spectra.name.split('.')[0]}-{name}", "config": config})
    return jobs


def entrapment_jobs(paired_fasta: Path, spectra: Path) -> list[dict]:
    """Paired one-to-one target and entrapment peptides."""
    jobs = []
    for name, mods in (("indexed", indexed(OXIDATION, "Oxidation")),
                       ("offset", offset(OXIDATION, "Oxidation"))):
        config = {
            "database": {
                "bucket_size": 16384,
                "enzyme": {"missed_cleavages": 0, "cleave_at": "$", "restrict": "",
                           "min_len": 7, "max_len": 50},
                "static_mods": {"C": 57.021464},
                "variable_mods": {"M": [dict(mods)], "^": [42.010565]},
                "max_variable_mods": 2, "max_total_variable_mods": 2,
                "decoy_tag": "rev_", "generate_decoys": True, "fasta": str(paired_fasta),
            },
            "precursor_tol": {"ppm": [-10.0, 10.0]},
            "fragment_tol": {"ppm": [-20.0, 20.0]},
            "isotope_errors": [-1, 3], "deisotope": True, "chimera": False,
            "max_fragment_charge": 2, "min_matched_peaks": 6, "report_psms": 1,
            "output_filter": {"psm_q_value": 1.0},
            "max_memory_gb": 22.0, "min_free_memory_gb": 2.0, "batch_size": 1,
            "mzml_paths": [str(spectra)], "score_type": "SageHyperScore",
        }
        jobs.append({"suite": "entrapment", "id": name, "config": config})
    return jobs


def execute(job: dict, sage: Path, root: Path, threads: int) -> dict:
    # Each job owns its directory, and earlier attempts remain intact.
    directory = root / job["suite"] / job["id"]
    if directory.exists():
        raise FileExistsError(f'Refusing to replace existing evidence: {directory}')
    directory.mkdir(parents=True)
    config = copy.deepcopy(job["config"])
    config["output_directory"] = str(directory)
    config_path = directory / "config.json"
    config_path.write_text(json.dumps(config, indent=2, sort_keys=True) + "\n")
    inputs = [Path(p) for p in config["mzml_paths"]] + [Path(config["database"]["fasta"])]
    # `/usr/bin/time` reports this job's own peak, unlike a cumulative
    # getrusage high-water mark shared by every child of this process.
    command = ["/usr/bin/time", "-v", str(sage), str(config_path), "--disable-telemetry-i-dont-want-to-improve-sage"]
    started = time.time()
    with (directory / "search.log").open("w") as log:
        completed = subprocess.run(command, stdout=subprocess.DEVNULL, stderr=log,
                                   cwd=REPO, env={"PATH": "/usr/bin:/bin",
                                                  "RAYON_NUM_THREADS": str(threads)})
    wall = time.time() - started
    text = (directory / "search.log").read_text()
    stage = lambda pattern: (int(match.group(1)) if (match := re.search(pattern, text)) else None)
    summary_path = directory / "run-summary.json"
    summary = json.loads(summary_path.read_text()) if summary_path.exists() else None
    outputs = sorted(p for p in directory.iterdir() if p.is_file() and p.name != "job.json")
    record = {
        "suite": job["suite"], "id": job["id"], "exit_status": completed.returncode,
        "command": command, "threads": threads,
        "wall_seconds": round(wall, 3),
        "max_rss_bytes": (stage(r"Maximum resident set size \(kbytes\): (\d+)") or 0) * 1024,
        "search_stage_ms": stage(r"search:\s+(\d+) ms"),
        "index_build_ms": stage(r"peptides in (\d+)\.\d+ms") ,
        "config": file_identities([config_path]),
        "inputs": file_identities(inputs),
        "outputs": file_identities(outputs),
        "summary": summary,
    }
    atomic_json(directory / "job.json", record)
    status = "ok" if completed.returncode == 0 else f"exit {completed.returncode}"
    print(f"{job['suite']}/{job['id']}: {status} wall {wall:.0f}s", flush=True)
    return record


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--sage", type=Path, required=True)
    parser.add_argument("--corpus", type=Path,
                        default=Path("/data/sage-plus-scientific/20260914"))
    parser.add_argument("--repeats", type=int, default=2)
    parser.add_argument("--threads", type=int, default=8)
    parser.add_argument("--paired-fasta", type=Path, default=None,
                        help="FDRBench paired target/entrapment reference; "
                             "defaults to inputs/paired.fasta under --root")
    args = parser.parse_args()

    spectra = sorted((args.corpus / "converted/PXD001468").glob("*.mgf"))[0]
    human = args.corpus / "references/human.fasta"
    synthetic = args.corpus / "ptm-truth/synthetic.fasta"
    phospho_files = sorted((args.corpus / "inputs/PXD000138").glob("HCD_*.mgf"))
    paired = args.paired_fasta or args.root / "inputs/paired.fasta"

    jobs = (scale_jobs(spectra, human, args.repeats)
            + localization_jobs(synthetic, phospho_files)
            + entrapment_jobs(paired, spectra))
    args.root.mkdir(parents=True, exist_ok=True)
    records = [execute(job, args.sage, args.root, args.threads) for job in jobs]
    atomic_json(args.root / "matrix.json", {
        "schema_version": 1,
        "executable": {str(args.sage.resolve()): sha256(args.sage)},
        "environment": {
            "platform": platform.platform(),
            "processor_count": len(__import__("os").sched_getaffinity(0)),
            "threads": args.threads,
        },
        "jobs": records,
    })
    failed = [r["id"] for r in records if r["exit_status"] != 0]
    print(f"completed {len(records) - len(failed)}/{len(records)} jobs", flush=True)


if __name__ == "__main__":
    main()
