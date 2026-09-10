#!/usr/bin/env python3
"""Repeat flagged paired timing cases with process CPU and host load diagnostics."""

import argparse
import json
import os
from pathlib import Path
import statistics
import subprocess

from provenance import atomic_json, file_identities


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--baseline", type=Path, required=True)
    parser.add_argument("--candidate", type=Path, required=True)
    parser.add_argument("--config", type=Path, action="append", required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--trials", type=int, default=5)
    parser.add_argument("--threads", type=int, default=8)
    args = parser.parse_args()
    if args.trials < 1 or args.threads < 1:
        parser.error("trials and threads must be positive")
    root = args.output.resolve()
    if root.exists() and any(root.iterdir()):
        parser.error("output must be a fresh directory")
    root.mkdir(parents=True, exist_ok=True)
    binaries = {key: getattr(args, key).resolve(strict=True) for key in ("baseline", "candidate")}
    configs = [path.resolve(strict=True) for path in args.config]
    if len({path.stem for path in configs}) != len(configs):
        parser.error("configuration stems must be distinct")
    inputs = list(binaries.values()) + configs + [Path("/usr/bin/time"), Path("/usr/bin/prlimit"), Path(__file__)]
    for config in configs:
        value = json.loads(config.read_text())
        inputs.extend(Path(path) for path in value["mzml_paths"])
        inputs.append(Path(value["database"]["fasta"]))
    manifest = {"status": "running", "inputs": file_identities(inputs), "threads": args.threads,
                "trials": args.trials, "warmups": 1, "host": list(os.uname()), "runs": []}
    atomic_json(root / "manifest.json", manifest)
    for config in configs:
        for trial in range(args.trials + 1):
            for engine in (list(binaries) if trial % 2 == 0 else list(reversed(binaries))):
                directory = root / config.stem / engine / str(trial)
                directory.mkdir(parents=True)
                timing = directory / "resources.tsv"
                command = ["/usr/bin/time", "-f", "%e\t%M\t%x\t%U\t%S", "-o", str(timing),
                           "/usr/bin/prlimit", f"--as={20 * 1024**3}", str(binaries[engine]),
                           str(config), "--output_directory", str(directory / "output"),
                           "--disable-telemetry-i-dont-want-to-improve-sage"]
                row = {"case": config.stem, "engine": engine, "trial": trial,
                       "command": command, "load_before": os.getloadavg()}
                print(config.stem, engine, trial, flush=True)
                with (directory / "run.log").open("wb") as log:
                    subprocess.run(command, env={**os.environ, "RAYON_NUM_THREADS": str(args.threads)},
                                   stdout=log, stderr=subprocess.STDOUT, check=True)
                wall, rss, status, user, system = timing.read_text().strip().split("\t")
                if status != "0":
                    raise RuntimeError("search timing indicates failure")
                row.update(wall_seconds=float(wall), peak_rss_mib=int(rss) / 1024,
                           cpu_seconds=float(user) + float(system), load_after=os.getloadavg())
                manifest["runs"].append(row)
                atomic_json(root / "manifest.json", manifest)
    manifest["medians"] = {}
    for config in configs:
        manifest["medians"][config.stem] = {}
        for engine in binaries:
            rows = [row for row in manifest["runs"] if row["case"] == config.stem and row["engine"] == engine and row["trial"]]
            manifest["medians"][config.stem][engine] = {
                key: statistics.median(row[key] for row in rows)
                for key in ("wall_seconds", "peak_rss_mib", "cpu_seconds")}
    if file_identities(inputs) != manifest["inputs"]:
        raise RuntimeError("benchmark inputs changed during execution")
    manifest["status"] = "completed"
    atomic_json(root / "manifest.json", manifest)
    print(json.dumps(manifest["medians"], indent=2))


if __name__ == "__main__":
    main()
