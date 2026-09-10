#!/usr/bin/env python3
"""Run paired hardening workloads using already frozen Sage Plus binaries."""

import argparse
import json
import os
from pathlib import Path
import statistics
import subprocess

from compare_runs import compare, load_rows, markdown
from provenance import atomic_json, file_identities, sha256
from run_fdrbench_validation import timing_values


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--baseline", type=Path, required=True)
    parser.add_argument("--candidate", type=Path, required=True)
    parser.add_argument("--standard", type=Path, required=True)
    parser.add_argument("--modified", type=Path, required=True)
    parser.add_argument("--feature", type=Path, required=True)
    parser.add_argument("--duckdb", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--threads", type=int, default=8)
    parser.add_argument("--trials", type=int, default=3)
    args = parser.parse_args()
    if args.threads < 1 or args.trials < 1:
        parser.error("threads and trials must be positive")
    root = args.output.resolve()
    if root.exists() and any(root.iterdir()):
        parser.error("output must be a fresh directory")
    root.mkdir(parents=True, exist_ok=True)
    binaries = {name: getattr(args, name).resolve(strict=True) for name in ("baseline", "candidate")}
    configurations = {name: json.loads(getattr(args, name).read_text()) for name in ("standard", "modified", "feature")}
    configurations["prefilter"] = json.loads(json.dumps(configurations["standard"]))
    configurations["prefilter"]["database"]["prefilter"] = True
    configurations["prefilter"]["database"]["prefilter_chunk_size"] = 1000
    inputs = list(binaries.values()) + [args.duckdb.resolve(strict=True), Path("/usr/bin/time"), Path("/usr/bin/prlimit")]
    for name, config in configurations.items():
        config.pop("output_directory", None)
        config["batch_size"] = 1
        for key in ("fasta", "peptides", "custom_cleavage_sites"):
            if config["database"].get(key):
                path = Path(config["database"][key]).resolve(strict=True)
                config["database"][key] = str(path)
                inputs.append(path)
        config["mzml_paths"] = [str(Path(value).resolve(strict=True)) for value in config["mzml_paths"]]
        inputs.extend(map(Path, config["mzml_paths"]))
        path = root / f"{name}.json"
        atomic_json(path, config)
        inputs.append(path)
    manifest = {"schema_version": 1, "status": "running", "inputs": file_identities(inputs),
                "threads": args.threads, "warmups": 1, "trials": args.trials, "runs": []}
    atomic_json(root / "manifest.json", manifest)
    environment = {**os.environ, "RAYON_NUM_THREADS": str(args.threads)}
    for case in configurations:
        # Alternate engine order between trials to reduce systematic order effects.
        for trial in range(args.trials + 1):
            for engine in (list(binaries) if trial % 2 == 0 else list(reversed(binaries))):
                directory = root / case / engine / ("warmup" if trial == 0 else f"trial-{trial}")
                directory.mkdir(parents=True)
                output, timing = directory / "output", directory / "timing.tsv"
                command = ["/usr/bin/time", "-f", "%e\t%M\t%x", "-o", str(timing),
                           "/usr/bin/prlimit", f"--as={20 * 1024**3}", str(binaries[engine]),
                           str(root / f"{case}.json"), "--output_directory", str(output),
                           "--disable-telemetry-i-dont-want-to-improve-sage"]
                print(f"{case} {engine} trial {trial}", flush=True)
                with (directory / "run.log").open("wb") as log:
                    subprocess.run(command, env=environment, stdout=log, stderr=subprocess.STDOUT, check=True)
                summary = json.loads((output / "run-summary.json").read_text())
                row = {"case": case, "engine": engine, "trial": trial, "command": command,
                       **timing_values(timing), "summary": summary,
                       "outputs": file_identities(sorted(path for path in output.iterdir() if path.is_file()))}
                manifest["runs"].append(row)
                atomic_json(root / "manifest.json", manifest)
    comparisons = {}
    for case in configurations:
        paths = [root / case / engine / "trial-1/output/results.sage.parquet" for engine in binaries]
        report = compare(*(load_rows(path, args.duckdb) for path in paths), 0.01)
        report["inputs"] = file_identities(paths)
        report["performance"] = {}
        for engine in binaries:
            rows = [row for row in manifest["runs"] if row["case"] == case and row["engine"] == engine and row["trial"]]
            report["performance"][engine] = {key: statistics.median(row[key] for row in rows) for key in ("wall_seconds", "peak_rss_mib")}
        report["review_flags"] = [key for key in ("wall_seconds", "peak_rss_mib")
                                  if report["performance"]["candidate"][key] > 1.1 * report["performance"]["baseline"][key]]
        if report["accepted_target_psms"]["candidate"] < 0.99 * report["accepted_target_psms"]["baseline"]:
            report["review_flags"].append("accepted_target_psm_loss")
        atomic_json(root / f"{case}-comparison.json", report)
        (root / f"{case}-comparison.md").write_text(markdown(report))
        comparisons[case] = report
    paths = [root / case / "candidate/trial-1/output/results.sage.parquet" for case in ("standard", "prefilter")]
    exact = compare(*(load_rows(path, args.duckdb) for path in paths), 0.01)
    exact["byte_identical"] = sha256(paths[0]) == sha256(paths[1])
    atomic_json(root / "prefilter-equivalence.json", exact)
    if file_identities(inputs) != manifest["inputs"]:
        raise RuntimeError("benchmark inputs changed during execution")
    manifest["status"] = "completed"
    manifest["comparisons"] = file_identities(list(root.glob("*-comparison.json")) + [root / "prefilter-equivalence.json"])
    atomic_json(root / "manifest.json", manifest)
    print(root)


if __name__ == "__main__":
    main()
