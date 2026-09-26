#!/usr/bin/env python3
"""Run the large-database searches and record their cost.

Reads a plan of named runs, each a complete Sage Plus configuration, and runs
them one at a time under GNU Time. Every run records the executable and
configuration hashes, exit status, wall time, peak resident memory, the
database preflight estimates, the prefilter counts from the log, and the
identification counts from `run-summary.json`. A run that already has a
record is skipped, so an interrupted plan resumes where it stopped.
"""

from __future__ import annotations

import argparse
import json
import os
import re
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from provenance import atomic_json, sha256  # noqa: E402

TIME = "/usr/bin/time"
GATE = "/mnt/data1/explore-data/memgate.sh"

PATTERNS = {
    "preflight": re.compile(
        r"database preflight: (\d+) unmodified peptides \(([\d.]+) GiB peak\), up to (\d+) "
        r"modified peptides \(([\d.]+) GiB peak\), up to (\d+) fragments \(([\d.]+) GiB index peak\)"
    ),
    "final": re.compile(
        r"final database preflight: (\d+) peptides, (\d+) fragments, estimated ([\d.]+) GiB peak"
    ),
    "prefilter": re.compile(
        r"prefilter search:\s+(\d+) ms \((\d+) peptides streamed, (\d+) retained\)"
    ),
    "spectra": re.compile(r"indexed (\d+) spectra: .* ([\d.]+) MiB, window depth"),
}


def parse_log(log: str) -> dict:
    record: dict = {}
    if match := PATTERNS["preflight"].search(log):
        record["preflight"] = {
            "unmodified_peptides": int(match[1]),
            "unmodified_gib": float(match[2]),
            "modified_peptides": int(match[3]),
            "modified_gib": float(match[4]),
            "fragments": int(match[5]),
            "fragment_index_gib": float(match[6]),
        }
    if match := PATTERNS["final"].search(log):
        record["final_preflight"] = {
            "peptides": int(match[1]),
            "fragments": int(match[2]),
            "gib": float(match[3]),
        }
    if match := PATTERNS["prefilter"].search(log):
        record["prefilter"] = {
            "seconds": int(match[1]) / 1000,
            "streamed": int(match[2]),
            "retained": int(match[3]),
        }
    if match := PATTERNS["spectra"].search(log):
        record["spectrum_index"] = {"spectra": int(match[1]), "mib": float(match[2])}
    errors = [line for line in log.splitlines()
              if any(word in line for word in ("Error", "exceeds", "memory allocation", "panicked"))]
    if errors:
        record["error_lines"] = errors[-5:]
    return record


def run(executable: Path, config: Path, output: Path, gate_gb: int, env: dict) -> dict:
    command = [GATE, str(gate_gb), TIME, "-v", str(executable), str(config), "-o",
               str(output), "--overwrite", "--disable-telemetry-i-dont-want-to-improve-sage"]
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
        "env": env,
        "exit_code": process.returncode,
        "wall_seconds": seconds,
        "peak_rss_gib": int(rss.group(1)) / 2**20 if rss else None,
        **parse_log(process.stderr),
    }
    summary = output / "run-summary.json"
    results = output / "results.sage.parquet"
    if process.returncode == 0 and summary.exists() and results.exists():
        data = json.loads(summary.read_text())
        record |= {
            "results_sha256": sha256(results),
            "database_peptides": data.get("peptides_in_database"),
            "psms_at_one_percent_fdr": data.get("psms_at_one_percent_fdr"),
            "peptides_at_one_percent_fdr": data.get("peptides_at_one_percent_fdr"),
        }
    return record


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--executable", type=Path, required=True)
    parser.add_argument("--plan", type=Path, required=True,
                        help="JSON list of {name, config, gate_gb, env}")
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--only", nargs="*", help="run only these names")
    args = parser.parse_args()

    executable = args.executable.resolve()
    report_path = args.output / "large-db-runs.json"
    report = json.loads(report_path.read_text()) if report_path.exists() else {"runs": {}}
    report["executable"] = {"path": str(executable), "sha256": sha256(executable)}
    for item in json.loads(args.plan.read_text()):
        name = item["name"]
        if args.only and name not in args.only:
            continue
        if name in report["runs"]:
            print(f"{name}: recorded, skipped", flush=True)
            continue
        directory = args.output / name
        directory.mkdir(parents=True, exist_ok=True)
        config = directory / "config.json"
        config.write_text(json.dumps(item["config"], indent=1))
        record = run(executable, config, directory, item.get("gate_gb", 18),
                     item.get("env", {}))
        record["config_sha256"] = sha256(config)
        record["meta"] = item.get("meta", {})
        report["runs"][name] = record
        atomic_json(report_path, report)
        print(f"{name}: exit {record['exit_code']}, {record['wall_seconds']:.1f} s, "
              f"{record['peak_rss_gib'] or 0:.2f} GiB, "
              f"{record.get('psms_at_one_percent_fdr')} PSMs", flush=True)


if __name__ == "__main__":
    main()
