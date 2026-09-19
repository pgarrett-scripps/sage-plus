#!/usr/bin/env python3
"""Run a frozen sequential search matrix and retain every success or failure."""

import argparse
import copy
import json
import os
import platform
import signal
import subprocess
import time
from pathlib import Path

from provenance import atomic_json, file_identities, sha256
from scientific_metrics import identification_metrics, normalize_psms, read_table


REPO = Path(__file__).resolve().parents[1]
BASELINE = REPO / "benchmarks/.work/targets/baseline-df9219951cc9a54c/release/sage"
CANDIDATE = REPO / "benchmarks/results/beta3-publication-20260910/hosted-gnu/unpacked/sage-plus-v0.1.0-beta.3-x86_64-unknown-linux-gnu/sage"


def base_config(fasta, spectra):
    return {"database": {"bucket_size": 16384, "enzyme": {
        "missed_cleavages": 1, "cleave_at": "KR", "restrict": "P", "min_len": 7, "max_len": 50},
        "static_mods": {"C": 57.021464}, "decoy_tag": "rev_", "generate_decoys": True, "fasta": str(fasta)},
        "precursor_tol": {"ppm": [-10, 10]}, "fragment_tol": {"ppm": [-20, 20]},
        "isotope_errors": [-1, 3], "deisotope": True, "chimera": False,
        "max_fragment_charge": 2, "min_matched_peaks": 6, "report_psms": 1,
        "output_filter": {"psm_q_value": 1.0}, "batch_size": 1,
        "mzml_paths": [str(p) for p in spectra], "score_type": "SageHyperScore"}


def config_inputs(config):
    paths = [Path(config["database"]["fasta"])]
    paths.extend(Path(p) for p in config["mzml_paths"])
    if library := config["database"].get("ptm_library"):
        paths.append(Path(library["path"]))
    return paths


def run_job(job, destination, timeout):
    output = destination / job["id"]
    output.mkdir(parents=True, exist_ok=True)
    settings = copy.deepcopy(job["config"])
    if job["engine"] == "plus":
        settings.update(max_memory_gb=10, min_free_memory_gb=2)
    binary = Path(job.get("binary") or (CANDIDATE if job["engine"] == "plus" else BASELINE))
    config = output / "config.json"
    atomic_json(config, settings)
    result = output / "result.json"
    inputs = file_identities([binary, config, *config_inputs(settings), Path(__file__).resolve(),
                              REPO / "benchmarks/scientific_metrics.py"])
    signature = {"job": job, "inputs": inputs}
    if result.exists():
        previous = json.loads(result.read_text())
        if previous["signature"] != signature:
            raise RuntimeError(f"Refusing to reuse changed job {job['id']}")
        if previous["status"] == "complete" and previous["outputs"] == file_identities([Path(p) for p in previous["outputs"]]):
            return previous
        raise RuntimeError(f"Previous job failed or outputs changed: {job['id']}. Use a new experiment ID.")
    available = next(int(line.split()[1]) * 1024 for line in Path("/proc/meminfo").read_text().splitlines() if line.startswith("MemAvailable:"))
    record = {"signature": signature, "status": "running", "available_memory_before": available}
    atomic_json(result, record)
    if available < 8 * 1024**3:
        record.update(status="resource_unavailable", reason="Less than 8 GiB available before search")
        atomic_json(result, record)
        return record
    environment = os.environ.copy()
    environment["RAYON_NUM_THREADS"] = str(job["threads"])
    timing = output / "timing.tsv"
    limit = int(job.get("address_space_gib", 16) * 1024**3)
    command = ["/usr/bin/time", "-f", "%e\t%M\t%x", "-o", str(timing),
               "/usr/bin/prlimit", f"--as={limit}:{limit}", str(binary), str(config),
               "--output_directory", str(output), "--batch-size", str(settings["batch_size"]),
               "--disable-telemetry-i-dont-want-to-improve-sage"]
    record["command"] = command
    start = time.monotonic()
    with (output / "search.log").open("w") as log:
        process = subprocess.Popen(command, cwd=REPO, env=environment, stdout=log,
                                   stderr=subprocess.STDOUT, start_new_session=True)
        try:
            code = process.wait(timeout=timeout)
        except subprocess.TimeoutExpired:
            os.killpg(process.pid, signal.SIGKILL)
            process.wait()
            code = None
    record.update(elapsed_seconds=time.monotonic() - start, exit_code=code)
    if timing.exists():
        lines = timing.read_text().splitlines()
        fields = lines[-1].split("\t") if lines else []
        if len(fields) == 3:
            record.update(wall_seconds=float(fields[0]), peak_rss_mib=int(fields[1]) / 1024)
    if code != 0:
        record["status"] = "timeout" if code is None else "failed"
    else:
        source = output / ("results.sage.parquet" if job["engine"] == "plus" else "results.sage.tsv")
        try:
            rows = normalize_psms(read_table(source))
            atomic_json(output / "identification-metrics.json", identification_metrics(rows))
            record.update(status="complete", rank_one_psms=len(rows),
                          decoy_rows_retained=sum(r["is_decoy"] for r in rows))
        except (ValueError, OSError, KeyError, subprocess.CalledProcessError) as error:
            record.update(status="analysis_failed", reason=str(error))
    record["outputs"] = file_identities([p for p in output.iterdir() if p.is_file() and p != result])
    if file_identities([Path(p) for p in inputs]) != inputs:
        record.update(status="invalid", reason="Inputs changed during execution")
    atomic_json(result, record)
    print(job["id"], record["status"], round(record["elapsed_seconds"], 2), flush=True)
    return record


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("plan", type=Path)
    p.add_argument("--output", type=Path, required=True)
    p.add_argument("--timeout", type=int, default=900)
    args = p.parse_args()
    plan = json.loads(args.plan.read_text())
    jobs = plan["jobs"]
    ids = [j["id"] for j in jobs]
    if len(ids) != len(set(ids)) or any(Path(i).name != i or i in (".", "..") for i in ids):
        raise ValueError("Job IDs must be unique single directory names")
    args.output.mkdir(parents=True, exist_ok=True)
    frozen = args.output / "plan.json"
    if frozen.exists() and json.loads(frozen.read_text()) != plan:
        raise RuntimeError("A different plan is already frozen here")
    atomic_json(frozen, plan)
    atomic_json(args.output / "environment.json", {"platform": platform.platform(),
                "python": platform.python_version(), "plan_sha256": sha256(frozen),
                "baseline_sha256": sha256(BASELINE), "candidate_sha256": sha256(CANDIDATE),
                "cpu": Path("/proc/cpuinfo").read_text(), "timeout_seconds_per_job": args.timeout})
    results = []
    for job in jobs:
        results.append(run_job(job, args.output, args.timeout))
        atomic_json(args.output / "matrix-status.json", {
            "expected_jobs": ids, "completed": [r["signature"]["job"]["id"] for r in results if r["status"] == "complete"],
            "failures": [{"id": r["signature"]["job"]["id"], "status": r["status"]} for r in results if r["status"] != "complete"],
            "pending": ids[len(results):]})
    if any(r["status"] != "complete" for r in results):
        raise SystemExit(1)


if __name__ == "__main__":
    main()
