#!/usr/bin/env python3
"""Run the deterministic charge-aware preprocessing and search benchmark."""

from __future__ import annotations

import argparse
import json
import os
import statistics
import subprocess
from pathlib import Path
from typing import Any

import benchmark as common


def parse_value(value: str) -> Any:
    if value == "true":
        return True
    if value == "false":
        return False
    try:
        return int(value)
    except ValueError:
        try:
            return float(value)
        except ValueError:
            return value


def parse_output(output: str) -> dict[str, Any]:
    values = {}
    for line in output.splitlines():
        if "=" not in line:
            continue
        key, value = line.split("=", 1)
        values[key] = parse_value(value)
    return values


def build_example() -> Path:
    target = common.WORK_ROOT / "targets" / "candidate"
    environment = os.environ.copy()
    environment["CARGO_TARGET_DIR"] = str(target)
    subprocess.run(
        [
            "cargo",
            "build",
            "--release",
            "--locked",
            "-p",
            "sage-core",
            "--example",
            "charge_matching_benchmark",
        ],
        cwd=common.REPO,
        env=environment,
        check=True,
    )
    binary = target / "release" / "examples" / "charge_matching_benchmark"
    if not binary.is_file():
        raise RuntimeError(f"benchmark example was not created: {binary}")
    return binary


def median(records: list[dict[str, Any]], field: str) -> float:
    return statistics.median(float(record[field]) for record in records)


def percent_delta(candidate: float, baseline: float) -> float:
    return (candidate - baseline) / baseline * 100


def render_report(metadata: dict[str, Any], records: list[dict[str, Any]]) -> str:
    modes = {}
    for record in records:
        modes.setdefault(record["mode"], []).append(record)
    rows = {}
    for mode, trials in modes.items():
        rows[mode] = {
            "trials": len(trials),
            "wall_seconds": median(trials, "wall_seconds"),
            "peak_rss_kib": median(trials, "peak_rss_kib"),
            "preprocess_ns_per_spectrum": median(trials, "preprocess_ns_per_spectrum"),
            "ns_per_search": median(trials, "ns_per_search"),
            "known_charges": median(trials, "known_charges"),
            "psms": median(trials, "psms"),
            "matched_peaks": median(trials, "matched_peaks"),
            "deterministic": len(
                {
                    (
                        trial["peptide_checksum"],
                        trial["score_checksum"],
                        trial["matched_peaks"],
                    )
                    for trial in trials
                }
            )
            == 1,
        }
    lines = [
        "# Charge-aware preprocessing benchmark",
        "",
        f"Created: {metadata['created_at_utc']}",
        "",
        f"- CPU: {metadata['cpu']}",
        f"- Candidate: `{metadata['candidate_commit']}`",
        f"- Rust: {metadata['rustc']}",
        "",
        "| Fragment charge source | Trials | Wall time | Peak RSS | Preprocessing per spectrum | Search per spectrum | Known charges | PSMs | Matched peaks | Deterministic |",
        "|---|---:|---:|---:|---:|---:|---:|---:|---:|---|",
    ]
    labels = {"inferred": "Inferred", "supplied": "Supplied array"}
    for mode in ("inferred", "supplied"):
        row = rows[mode]
        lines.append(
            f"| {labels[mode]} | {row['trials']} | {row['wall_seconds']:.2f} s | "
            f"{row['peak_rss_kib'] / 1024:.1f} MiB | "
            f"{row['preprocess_ns_per_spectrum'] / 1000:.2f} us | "
            f"{row['ns_per_search'] / 1000:.2f} us | {row['known_charges']:,.0f} | "
            f"{row['psms']:,.0f} | {row['matched_peaks']:,.0f} | "
            f"{'yes' if row['deterministic'] else 'no'} |"
        )
    inferred = rows["inferred"]
    supplied = rows["supplied"]
    preprocess_change = percent_delta(
        supplied["preprocess_ns_per_spectrum"], inferred["preprocess_ns_per_spectrum"]
    )
    search_change = percent_delta(supplied["ns_per_search"], inferred["ns_per_search"])
    lines.extend(
        [
            "",
            "## Summary",
            "",
            f"Supplied fragment charges changed preprocessing time by {preprocess_change:.1f} percent and search time by {search_change:.1f} percent. Both modes returned {supplied['psms']:,.0f} PSMs across the fixed synthetic workload.",
            "",
        ]
    )
    return "\n".join(lines)


def run(args: argparse.Namespace) -> Path:
    if args.repeats < 1:
        raise RuntimeError("--repeats must be at least one")
    if args.warmups < 0:
        raise RuntimeError("--warmups must not be negative")
    common.ensure_tools()
    binary = build_example()
    result_directory = common.timestamped_result_directory("charge")
    metadata = {
        "created_at_utc": common.dt.datetime.now(common.dt.timezone.utc).isoformat(),
        "cpu": common.cpu_model(),
        "candidate_commit": common.git_commit("HEAD"),
        "rustc": common.run_text(["rustc", "--version"]),
        "repeats": args.repeats,
        "warmups": args.warmups,
    }
    common.write_json(result_directory / "metadata.json", metadata)
    records = []
    for mode in ("inferred", "supplied"):
        for index in range(args.warmups + args.repeats):
            warmup = index < args.warmups
            measured = index - args.warmups + 1
            name = f"warmup-{index + 1}" if warmup else f"trial-{measured}"
            trial = result_directory / "runs" / mode / name
            trial.mkdir(parents=True)
            timing = trial / "timing.tsv"
            environment = os.environ.copy()
            if mode == "supplied":
                environment["SAGE_USE_FRAGMENT_CHARGES"] = "1"
            else:
                environment.pop("SAGE_USE_FRAGMENT_CHARGES", None)
            completed = subprocess.run(
                [
                    str(common.GNU_TIME),
                    "-f",
                    "%e\t%M\t%x",
                    "-o",
                    str(timing),
                    str(binary),
                ],
                cwd=common.REPO,
                env=environment,
                text=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
            )
            (trial / "stdout.log").write_text(completed.stdout)
            (trial / "stderr.log").write_text(completed.stderr)
            if completed.returncode != 0:
                raise RuntimeError(f"charge benchmark failed in {mode}/{name}")
            if not warmup:
                record = {
                    "mode": mode,
                    "trial": measured,
                    **common.parse_timing(timing),
                    **parse_output(completed.stdout),
                }
                records.append(record)
                common.write_json(result_directory / "records.json", records)
                (result_directory / "report.md").write_text(
                    render_report(metadata, records)
                    if {record["mode"] for record in records} == {"inferred", "supplied"}
                    else "# Charge-aware preprocessing benchmark\n\nBenchmark is still running.\n"
                )
    (result_directory / "report.md").write_text(render_report(metadata, records))
    print(result_directory)
    return result_directory


def parser() -> argparse.ArgumentParser:
    value = argparse.ArgumentParser(description=__doc__)
    value.add_argument("--repeats", type=int, default=3)
    value.add_argument("--warmups", type=int, default=1)
    return value


def main() -> int:
    try:
        run(parser().parse_args())
    except (OSError, RuntimeError, subprocess.CalledProcessError) as error:
        print(f"charge benchmark failed: {error}", file=common.sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
