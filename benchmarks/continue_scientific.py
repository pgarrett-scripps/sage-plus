#!/usr/bin/env python3
"""Continue the selected pilot sequentially as acquisition completes."""

import argparse
import copy
import json
import os
import signal
import subprocess
import sys
import time
from pathlib import Path

from provenance import atomic_json
from run_scientific import REPO, base_config


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--root", type=Path, required=True)
    p.add_argument("--resume-after-reference-repair", action="store_true")
    args = p.parse_args()
    root = args.root.resolve()
    steps = []
    log_root = root / ("continuation-reference-repair" if args.resume_after_reference_repair else "continuation")
    log_root.mkdir(exist_ok=True)

    def stage(name, command, timeout=3600):
        start = time.monotonic()
        with (log_root / f"{name}.log").open("w") as log:
            process = subprocess.Popen(command, cwd=REPO, stdout=log, stderr=subprocess.STDOUT)
            try:
                code = process.wait(timeout=timeout)
            except subprocess.TimeoutExpired:
                def stop_tree(pid):
                    children = Path(f"/proc/{pid}/task/{pid}/children")
                    try:
                        descendants = children.read_text().split()
                    except FileNotFoundError:
                        return
                    for child in descendants:
                        stop_tree(int(child))
                    try:
                        os.kill(pid, signal.SIGKILL)
                    except ProcessLookupError:
                        pass
                stop_tree(process.pid)
                process.wait()
                code = None
        steps.append({"name": name, "command": command, "exit_code": code,
                      "elapsed_seconds": time.monotonic() - start})
        atomic_json(log_root / "stages.json", steps)
        print(name, "complete" if code == 0 else f"failed ({code})", flush=True)
        return code == 0

    # The earlier scaling matrix owns the search resource until every cell finishes.
    deadline = time.monotonic() + 4 * 3600
    while True:
        state = root / "runs/scaling/matrix-status.json"
        if state.exists() and not json.loads(state.read_text())["pending"]:
            break
        if time.monotonic() > deadline:
            raise TimeoutError("Scaling did not finish within the pilot budget")
        time.sleep(5)

    hek = sorted((root / "converted/PXD001468").glob("*.mgf"))
    if len(hek) != 2:
        raise ValueError("Exactly two prespecified converted HEK files are required")
    for seed in (20260914, 20260915, 20260916):
        if args.resume_after_reference_repair:
            continue
        stage(f"entrapment-human-{seed}", [sys.executable, "benchmarks/scientific_entrapment.py",
              "--reference", str(root / "references/human.fasta"), "--spectra", *map(str, hek),
              "--output", str(root / f"runs/entrapment-human-{seed}"), "--seed", str(seed)])

    names = [f"LFQ_Orbitrap_DDA_Condition_{condition}_Sample_{sample}_01.raw"
             for sample in ("Alpha", "Beta") for condition in ("A", "B")]
    names.append("LFQ_Orbitrap_DDA_Human_01.raw")
    converted = []
    for index, name in enumerate(names):
        source = root / "inputs/PXD028735" / name
        receipt = source.with_name(source.name + ".receipt.json")
        while not receipt.exists():
            if time.monotonic() > deadline:
                raise TimeoutError(f"Acquisition did not finish: {name}")
            time.sleep(5)
        if not args.resume_after_reference_repair and not stage(f"convert-lfq-{index}", [sys.executable, "benchmarks/prepare_scientific.py",
                     str(source), "--root", str(root)], timeout=1900):
            raise RuntimeError(f"Conversion failed: {name}")
        converted.append(root / "converted/PXD028735" / (source.stem + ".mzML"))
        if not converted[-1].exists():
            raise ValueError(f"Missing converted input: {converted[-1]}")

    reference = root / ("references/hye-irt-defined.fasta" if args.resume_after_reference_repair else "references/hye-irt.fasta")
    if args.resume_after_reference_repair:
        original = root / "runs/public-comparison/matrix-status.json"
        while json.loads(original.read_text())["pending"]:
            if time.monotonic() > deadline:
                raise TimeoutError("The original public matrix did not finish")
            time.sleep(5)
    # Paired raw identification comparisons, with identical complete input files.
    jobs = []
    for study, spectra, fasta in (("PXD001468", hek, root / "references/human.fasta"),
                                   ("PXD028735", converted[:4], reference)):
        if args.resume_after_reference_repair and study == "PXD001468":
            continue
        for index, spectrum in enumerate(spectra):
            for engine in ("upstream", "plus"):
                jobs.append({"id": f"{study}-{index}-{engine}", "study": study,
                             "workload": "public-identification", "engine": engine,
                             "threads": 8, "config": base_config(fasta, [spectrum])})
    suite = "public-comparison-v2" if args.resume_after_reference_repair else "public-comparison"
    plan = root / f"{suite}-plan.json"
    atomic_json(plan, {"status": "pilot", "jobs": jobs})
    stage("public-comparison", [sys.executable, "benchmarks/run_scientific.py", str(plan),
          "--output", str(root / "runs" / suite)])

    jobs = []
    for mbr in (False, True):
        config = base_config(reference, converted)
        config["quant"] = {"lfq": True, "lfq_settings": {"mbr": mbr}}
        jobs.append({"id": f"mbr-{int(mbr)}", "study": "PXD028735", "workload": "known-ratio-lfq",
                     "engine": "plus", "threads": 8, "config": config})
    plan = root / "quantification-plan.json"
    atomic_json(plan, {"status": "pilot", "jobs": jobs,
                "primary_ratios": "B over A: human 1, yeast 0.5, E. coli 4",
                "control": "Pure human. Species attribution must exclude peptides shared with human."})
    stage("quantification", [sys.executable, "benchmarks/run_scientific.py", str(plan),
          "--output", str(root / "runs/quantification"), "--timeout", "1800"], timeout=4000)

    for seed in (20260914, 20260915, 20260916):
        stage(f"entrapment-hye-{seed}", [sys.executable, "benchmarks/scientific_entrapment.py",
              "--reference", str(reference), "--spectra", *map(str, converted[:2]),
              "--output", str(root / f"runs/entrapment-hye-{seed}"), "--seed", str(seed)])
    if args.resume_after_reference_repair:
        stage("public-timing", [sys.executable, "benchmarks/run_scientific.py",
              str(root / "public-repeat-v2-plan.json"), "--output", str(root / "runs/public-timing")])
    stage("analysis", [sys.executable, "benchmarks/analyze_scientific.py", "--root", str(root),
                       "--output", str(root / "pilot-summary.json")])
    print("Sequential pilot stages finished. Inspect stages.json and result matrices for failures.", flush=True)


if __name__ == "__main__":
    main()
