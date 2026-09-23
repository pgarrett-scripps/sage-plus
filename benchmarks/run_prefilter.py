#!/usr/bin/env python3
"""Paired prefilter benchmark with recorded identities.

Runs each configuration with a baseline and a candidate executable, both with
prefiltering enabled, and optionally the candidate with prefiltering disabled.
Engine order alternates between repeats. Every run records executable and
configuration hashes, wall time, peak RSS, the prefilter phase from log
timestamps, identification counts from `run-summary.json`, and the
`results.sage.parquet` hash.

Byte-identical result hashes establish identical output. For different hashes,
use `compare_runs.py` to separate search features from rescoring differences.
"""

from __future__ import annotations

import argparse
import copy
import json
import os
import re
import statistics
import subprocess
from datetime import datetime
from pathlib import Path

from provenance import atomic_json, sha256

TIME = "/usr/bin/time"


def prefilter_seconds(log: str) -> float | None:
    """Seconds from database chunking to the end of the prefilter stage."""
    lines = log.splitlines()
    start = next((line for line in lines if "db chunks of size" in line), None)
    ends = [
        line
        for line in lines
        if "pre-filtered peptides for fasta chunk" in line or "peptides streamed" in line
    ]
    if not start or not ends:
        return None
    stamp = lambda line: datetime.strptime(line[1:20], "%Y-%m-%dT%H:%M:%S")
    return (stamp(ends[-1]) - stamp(start)).total_seconds()


def run(executable: Path, config: Path, output: Path, env: dict[str, str]) -> dict:
    output.mkdir(parents=True, exist_ok=True)
    command = [TIME, "-v", str(executable), str(config), "-o", str(output),
               "--overwrite", "--disable-telemetry-i-dont-want-to-improve-sage"]
    process = subprocess.run(command, capture_output=True, text=True,
                             env={**os.environ, **env})
    (output / "stderr.log").write_text(process.stderr)
    wall = re.search(r"Elapsed \(wall clock\) time.*: (\S+)", process.stderr)
    rss = re.search(r"Maximum resident set size \(kbytes\): (\d+)", process.stderr)
    seconds = 0.0
    for part in (wall.group(1).split(":") if wall else []):
        seconds = seconds * 60 + float(part)
    record = {
        "command": command,
        "exit_code": process.returncode,
        "wall_seconds": seconds,
        "peak_rss_gib": int(rss.group(1)) / 2**20 if rss else None,
        "prefilter_seconds": prefilter_seconds(process.stderr),
    }
    summary = output / "run-summary.json"
    results = output / "results.sage.parquet"
    if process.returncode == 0 and summary.exists() and results.exists():
        data = json.loads(summary.read_text())
        record |= {
            "results_sha256": sha256(results),
            "database_peptides": data["peptides_in_database"],
            "psms_at_one_percent_fdr": data["psms_at_one_percent_fdr"],
            "peptides_at_one_percent_fdr": data["peptides_at_one_percent_fdr"],
        }
    return record


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--baseline", type=Path, required=True)
    parser.add_argument("--candidate", type=Path, required=True)
    parser.add_argument("--config", type=Path, action="append", required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--repeats", type=int, default=3)
    parser.add_argument("--prefilter-off", action="store_true",
                        help="also run the candidate with prefiltering disabled")
    parser.add_argument("--env", action="append", default=[],
                        help="KEY=VALUE applied to candidate runs")
    args = parser.parse_args()

    env = dict(item.split("=", 1) for item in args.env)
    engines = {"baseline": args.baseline.resolve(), "candidate": args.candidate.resolve()}
    report = {
        "executables": {name: {"path": str(path), "sha256": sha256(path)}
                        for name, path in engines.items()},
        "candidate_env": env,
        "configs": {},
    }
    for config in args.config:
        name = config.stem
        base = json.loads(config.read_text())
        base["database"]["prefilter"] = True
        variants = {"baseline": (engines["baseline"], base, {}),
                    "candidate": (engines["candidate"], base, env)}
        if args.prefilter_off:
            off = copy.deepcopy(base)
            off["database"]["prefilter"] = False
            variants["candidate-off"] = (engines["candidate"], off, env)
        runs: dict[str, list[dict]] = {variant: [] for variant in variants}
        for repeat in range(1, args.repeats + 1):
            order = list(variants) if repeat % 2 else list(reversed(variants))
            for variant in order:
                executable, content, variant_env = variants[variant]
                directory = args.output / name / f"{variant}-r{repeat}"
                path = directory / "config.json"
                directory.mkdir(parents=True, exist_ok=True)
                path.write_text(json.dumps(content, indent=1))
                record = run(executable, path, directory, variant_env)
                record["config_sha256"] = sha256(path)
                runs[variant].append(record)
                print(f"{name} {variant} r{repeat}: exit {record['exit_code']}, "
                      f"{record['wall_seconds']:.1f} s", flush=True)

        def median(variant: str, key: str) -> float | None:
            values = [r[key] for r in runs[variant] if r["exit_code"] == 0 and r.get(key) is not None]
            return statistics.median(values) if values else None

        hashes = {variant: sorted({r.get("results_sha256") for r in records})
                  for variant, records in runs.items()}
        report["configs"][name] = {
            "runs": runs,
            "median": {variant: {key: median(variant, key)
                                 for key in ("wall_seconds", "peak_rss_gib", "prefilter_seconds")}
                       for variant in variants},
            "results_sha256": hashes,
            "byte_identical": len({h for values in hashes.values() for h in values}) == 1,
        }
        atomic_json(args.output / "prefilter-benchmark.json", report)


if __name__ == "__main__":
    main()
