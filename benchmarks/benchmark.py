#!/usr/bin/env python3
"""Small repeatable benchmark harness for Sage Plus."""

from __future__ import annotations

import argparse
import datetime as dt
import hashlib
import json
import os
import platform
import shutil
import statistics
import subprocess
import sys
import tarfile
from pathlib import Path, PurePosixPath
from typing import Any


REPO = Path(__file__).resolve().parents[1]
BENCHMARK_ROOT = REPO / "benchmarks"
WORK_ROOT = BENCHMARK_ROOT / ".work"
RESULTS_ROOT = BENCHMARK_ROOT / "results"
GNU_TIME = Path("/usr/bin/time")
DEFAULT_BASELINE = "v0.1.0-beta.1"
SUMMARY_FIELDS = (
    "runtime_secs",
    "files",
    "peptides_in_database",
    "fragments_in_database",
    "psms_at_one_percent_fdr",
    "peptides_at_one_percent_fdr",
    "proteins_at_one_percent_fdr",
    "protein_groups_at_one_percent_fdr",
)


def run_text(command: list[str], cwd: Path = REPO) -> str:
    result = subprocess.run(
        command,
        cwd=cwd,
        check=True,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
    )
    return result.stdout.strip()


def git_commit(reference: str) -> str:
    return run_text(["git", "rev-parse", f"{reference}^{{commit}}"])


def candidate_label() -> str:
    commit = git_commit("HEAD")[:12]
    dirty = subprocess.run(["git", "diff", "--quiet"], cwd=REPO).returncode != 0
    return f"{commit}-dirty" if dirty else commit


def safe_reference_name(commit: str) -> str:
    return commit[:16]


def ensure_tools() -> None:
    if not GNU_TIME.is_file():
        raise RuntimeError("GNU /usr/bin/time is required")
    for tool in ("cargo", "git", "rustc"):
        if shutil.which(tool) is None:
            raise RuntimeError(f"required tool is not available: {tool}")


def extract_baseline(commit: str) -> Path:
    destination = WORK_ROOT / "sources" / safe_reference_name(commit)
    marker = destination / ".benchmark-source"
    if marker.is_file() and marker.read_text().strip() == commit:
        return destination
    if destination.exists():
        raise RuntimeError(f"incomplete baseline source directory exists: {destination}")
    destination.mkdir(parents=True)
    process = subprocess.Popen(
        ["git", "archive", "--format=tar", commit],
        cwd=REPO,
        stdout=subprocess.PIPE,
    )
    if process.stdout is None:
        raise RuntimeError("could not read git archive output")
    with tarfile.open(fileobj=process.stdout, mode="r|") as archive:
        for member in archive:
            path = PurePosixPath(member.name)
            if path.is_absolute() or ".." in path.parts:
                process.kill()
                raise RuntimeError(f"unsafe path in git archive: {member.name}")
            archive.extract(member, destination)
    return_code = process.wait()
    if return_code != 0:
        raise RuntimeError(f"git archive failed with exit code {return_code}")
    marker.write_text(f"{commit}\n")
    return destination


def cargo_build(source: Path, target: Path) -> Path:
    environment = os.environ.copy()
    environment["CARGO_TARGET_DIR"] = str(target)
    command = ["cargo", "build", "--release", "--locked", "-p", "sage-cli"]
    subprocess.run(command, cwd=source, env=environment, check=True)
    binary = target / "release" / "sage"
    if not binary.is_file():
        raise RuntimeError(f"Sage binary was not created: {binary}")
    return binary


def build_binaries(baseline_ref: str, include_baseline: bool = True) -> dict[str, Any]:
    ensure_tools()
    baseline_commit = git_commit(baseline_ref)
    targets = WORK_ROOT / "targets"
    candidate_binary = cargo_build(REPO, targets / "candidate")
    baseline_binary = None
    if include_baseline:
        baseline_source = extract_baseline(baseline_commit)
        baseline_binary = cargo_build(
            baseline_source,
            targets / f"baseline-{safe_reference_name(baseline_commit)}",
        )
    return {
        "baseline_commit": baseline_commit,
        "baseline_binary": baseline_binary,
        "candidate_commit": git_commit("HEAD"),
        "candidate_binary": candidate_binary,
        "candidate_label": candidate_label(),
    }


def total_memory_bytes() -> int | None:
    meminfo = Path("/proc/meminfo")
    if not meminfo.is_file():
        return None
    for line in meminfo.read_text().splitlines():
        if line.startswith("MemTotal:"):
            return int(line.split()[1]) * 1024
    return None


def cpu_model() -> str:
    cpuinfo = Path("/proc/cpuinfo")
    if cpuinfo.is_file():
        for line in cpuinfo.read_text(errors="replace").splitlines():
            if line.startswith("model name"):
                return line.split(":", 1)[1].strip()
    return platform.processor() or "unknown"


def machine_metadata(build: dict[str, Any], baseline_ref: str) -> dict[str, Any]:
    status = run_text(["git", "status", "--short"])
    return {
        "created_at_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        "platform": platform.platform(),
        "cpu": cpu_model(),
        "logical_cpus": os.cpu_count(),
        "total_memory_bytes": total_memory_bytes(),
        "python": platform.python_version(),
        "rustc": run_text(["rustc", "--version"]),
        "cargo": run_text(["cargo", "--version"]),
        "baseline_ref": baseline_ref,
        "baseline_commit": build["baseline_commit"],
        "candidate_commit": build["candidate_commit"],
        "candidate_dirty": bool(status),
        "candidate_status": status.splitlines(),
    }


def timestamped_result_directory(suite: str) -> Path:
    stamp = dt.datetime.now().strftime("%Y%m%d-%H%M%S")
    result = RESULTS_ROOT / f"{stamp}-{suite}"
    suffix = 1
    while result.exists():
        result = RESULTS_ROOT / f"{stamp}-{suite}-{suffix}"
        suffix += 1
    result.mkdir(parents=True)
    return result


def read_json(path: Path) -> dict[str, Any]:
    with path.open() as handle:
        value = json.load(handle)
    if not isinstance(value, dict):
        raise RuntimeError(f"expected a JSON object in {path}")
    return value


def write_json(path: Path, value: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as handle:
        json.dump(value, handle, indent=2, sort_keys=True)
        handle.write("\n")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def parse_timing(path: Path) -> dict[str, Any]:
    lines = [line for line in path.read_text().splitlines() if line.strip()]
    if not lines:
        raise RuntimeError(f"GNU Time did not write timing data: {path}")
    fields = lines[-1].split("\t")
    if len(fields) != 3:
        raise RuntimeError(f"unexpected GNU Time output: {lines[-1]}")
    return {
        "wall_seconds": float(fields[0]),
        "peak_rss_kib": int(fields[1]),
        "exit_code": int(fields[2]),
    }


def selected_summary(summary: dict[str, Any]) -> dict[str, Any]:
    selected = {field: summary.get(field) for field in SUMMARY_FIELDS}
    for field in (
        "models",
        "quantification",
        "modifications",
        "spectral_library",
        "library_search",
    ):
        if field in summary:
            selected[field] = summary[field]
    return selected


class BenchmarkSession:
    def __init__(
        self,
        suite: str,
        build: dict[str, Any],
        baseline_ref: str,
        repeats: int,
        warmups: int,
        threads: int,
    ) -> None:
        self.result_directory = timestamped_result_directory(suite)
        self.build = build
        self.repeats = repeats
        self.warmups = warmups
        self.threads = threads
        self.records: list[dict[str, Any]] = []
        self.metadata = machine_metadata(build, baseline_ref)
        self.metadata.update(
            {"suite": suite, "repeats": repeats, "warmups": warmups, "threads": threads}
        )
        write_json(self.result_directory / "metadata.json", self.metadata)

    def copy_config(self, name: str, config: Path) -> Path:
        destination = self.result_directory / "configs" / f"{name}.json"
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(config, destination)
        return destination

    def run_search_case(
        self,
        case: str,
        version: str,
        binary: Path,
        config: Path,
    ) -> None:
        keep_outputs = os.environ.get("SAGE_BENCH_KEEP_OUTPUTS") == "1"
        total = self.warmups + self.repeats
        for index in range(total):
            warmup = index < self.warmups
            measured_index = index - self.warmups + 1
            trial_name = f"warmup-{index + 1}" if warmup else f"trial-{measured_index}"
            trial = self.result_directory / "runs" / case / version / trial_name
            output = trial / "output"
            trial.mkdir(parents=True, exist_ok=True)
            timing_path = trial / "timing.tsv"
            stdout_path = trial / "stdout.log"
            stderr_path = trial / "stderr.log"
            command = [
                str(GNU_TIME),
                "-f",
                "%e\t%M\t%x",
                "-o",
                str(timing_path),
                str(binary),
                str(config),
                "--output_directory",
                str(output),
                "--batch-size",
                "1",
                "--disable-telemetry-i-dont-want-to-improve-sage",
            ]
            environment = os.environ.copy()
            environment["RAYON_NUM_THREADS"] = str(self.threads)
            with stdout_path.open("wb") as stdout, stderr_path.open("wb") as stderr:
                completed = subprocess.run(
                    command,
                    cwd=REPO,
                    env=environment,
                    stdout=stdout,
                    stderr=stderr,
                )
            timing = parse_timing(timing_path)
            summary_path = output / "run-summary.json"
            summary = read_json(summary_path) if summary_path.is_file() else {}
            artifact_hashes = {}
            if output.is_dir():
                for artifact in sorted(output.iterdir()):
                    if artifact.is_file() and artifact.name not in {"results.json", "run-summary.json"}:
                        artifact_hashes[artifact.name] = sha256(artifact)
            write_json(trial / "artifact-hashes.json", artifact_hashes)
            if summary:
                write_json(trial / "run-summary.json", summary)
            if not keep_outputs and output.exists():
                shutil.rmtree(output)
            if completed.returncode != 0:
                failure = {
                    "case": case,
                    "version": version,
                    "trial": trial_name,
                    "status": "failed",
                    "process_return_code": completed.returncode,
                    **timing,
                }
                self.records.append(failure)
                self.save()
                raise RuntimeError(
                    f"benchmark case {case}/{version}/{trial_name} failed with "
                    f"exit code {completed.returncode}. See {stderr_path}"
                )
            if not warmup:
                record = {
                    "case": case,
                    "version": version,
                    "trial": measured_index,
                    "status": "ok",
                    "process_return_code": completed.returncode,
                    **timing,
                    **selected_summary(summary),
                    "artifact_hashes": artifact_hashes,
                }
                self.records.append(record)
                self.save()

    def save(self) -> None:
        write_json(self.result_directory / "records.json", self.records)
        (self.result_directory / "report.md").write_text(
            render_report(self.metadata, self.records)
        )


def median(values: list[float | int | None]) -> float | None:
    finite = [float(value) for value in values if isinstance(value, (int, float))]
    return statistics.median(finite) if finite else None


def aggregate(records: list[dict[str, Any]]) -> list[dict[str, Any]]:
    groups: dict[tuple[str, str], list[dict[str, Any]]] = {}
    for record in records:
        if record.get("status", "ok") != "ok":
            continue
        groups.setdefault((record["case"], record["version"]), []).append(record)
    rows = []
    for (case, version), trials in sorted(groups.items()):
        row: dict[str, Any] = {
            "case": case,
            "version": version,
            "trials": len(trials),
            "wall_seconds": median([trial.get("wall_seconds") for trial in trials]),
            "peak_rss_kib": median([trial.get("peak_rss_kib") for trial in trials]),
        }
        for field in SUMMARY_FIELDS[2:]:
            row[field] = median([trial.get(field) for trial in trials])
        hashes = [
            trial.get("artifact_hashes", {}).get("results.sage.parquet") for trial in trials
        ]
        row["result_hashes"] = sorted({value for value in hashes if value})
        rows.append(row)
    return rows


def format_number(value: Any, digits: int = 1) -> str:
    if value is None:
        return "n/a"
    if isinstance(value, float):
        return f"{value:,.{digits}f}"
    return f"{value:,}"


def format_rss(value: float | None) -> str:
    if value is None:
        return "n/a"
    return f"{value / 1024:.1f} MiB"


def percent_delta(candidate: float | None, baseline: float | None) -> float | None:
    if candidate is None or baseline in (None, 0):
        return None
    return (candidate - baseline) / baseline * 100


def result_lookup(rows: list[dict[str, Any]], case: str, version: str) -> dict[str, Any] | None:
    return next(
        (row for row in rows if row["case"] == case and row["version"] == version),
        None,
    )


def comparison_findings(rows: list[dict[str, Any]]) -> list[str]:
    findings = []
    comparisons = (
        ("standard", "baseline", "candidate", "Baseline to candidate"),
        ("memory", "baseline", "candidate", "Synthetic memory"),
        ("prefilter", "off", "on", "Prefilter off to on"),
    )
    for case, left_name, right_name, label in comparisons:
        left = result_lookup(rows, case, left_name)
        right = result_lookup(rows, case, right_name)
        if left is None or right is None:
            continue
        time_delta = percent_delta(right["wall_seconds"], left["wall_seconds"])
        rss_delta = percent_delta(right["peak_rss_kib"], left["peak_rss_kib"])
        if time_delta is None:
            findings.append(f"INFO: {label} wall-time change was not measurable")
        else:
            time_status = "REVIEW" if time_delta > 10 else "PASS"
            findings.append(
                f"{time_status}: {label} wall-time change is "
                f"{format_number(time_delta)} percent"
            )
        if rss_delta is None:
            findings.append(f"INFO: {label} peak-RSS change was not measurable")
        else:
            rss_status = "REVIEW" if rss_delta > 10 else "PASS"
            findings.append(
                f"{rss_status}: {label} peak-RSS change is "
                f"{format_number(rss_delta)} percent"
            )
        if case == "standard":
            for field, name in (
                ("psms_at_one_percent_fdr", "PSMs"),
                ("peptides_at_one_percent_fdr", "peptides"),
            ):
                change = percent_delta(right.get(field), left.get(field))
                if change is None:
                    findings.append(
                        f"INFO: Baseline to candidate one-percent FDR {name} change "
                        "was not measurable"
                    )
                else:
                    status = "REVIEW" if change < -1 else "PASS"
                    findings.append(
                        f"{status}: Baseline to candidate one-percent FDR {name} change is "
                        f"{format_number(change)} percent"
                    )
        if case == "prefilter":
            counts_equal = all(
                left.get(field) == right.get(field)
                for field in (
                    "psms_at_one_percent_fdr",
                    "peptides_at_one_percent_fdr",
                    "proteins_at_one_percent_fdr",
                    "protein_groups_at_one_percent_fdr",
                )
            )
            hashes_equal = bool(left["result_hashes"]) and (
                left["result_hashes"] == right["result_hashes"]
            )
            findings.append(
                f"{'PASS' if counts_equal else 'REVIEW'}: Prefilter one-percent FDR counts "
                f"{'match' if counts_equal else 'differ'}"
            )
            findings.append(
                f"{'PASS' if hashes_equal else 'REVIEW'}: Prefilter result hashes "
                f"{'match' if hashes_equal else 'differ'}"
            )
    return findings


def render_report(metadata: dict[str, Any], records: list[dict[str, Any]]) -> str:
    rows = aggregate(records)
    failures = [record for record in records if record.get("status") == "failed"]
    lines = [
        "# Sage Plus benchmark report",
        "",
        f"Created: {metadata['created_at_utc']}",
        "",
        "## Environment",
        "",
        f"- CPU: {metadata['cpu']}",
        f"- Logical CPUs: {metadata['logical_cpus']}",
        f"- Threads used: {metadata['threads']}",
        f"- Baseline: `{metadata['baseline_commit']}`",
        f"- Candidate: `{metadata['candidate_commit']}`",
        f"- Candidate working tree dirty: `{str(metadata['candidate_dirty']).lower()}`",
        f"- Rust: {metadata['rustc']}",
    ]
    if failures:
        lines.extend(
            [
                "",
                "## Failures",
                "",
                "| Case | Version | Trial | Process return code | Wall time | Peak RSS |",
                "|---|---|---|---:|---:|---:|",
            ]
        )
        for failure in failures:
            lines.append(
                f"| {failure['case']} | {failure['version']} | {failure['trial']} | "
                f"{failure['process_return_code']} | {format_number(failure['wall_seconds'], 2)} s | "
                f"{format_rss(failure['peak_rss_kib'])} |"
            )
    lines.extend(
        [
            "",
            "## Median results",
            "",
            "| Case | Version | Trials | Wall time | Peak RSS | PSMs at 1% FDR | Peptides at 1% FDR | Database peptides |",
            "|---|---|---:|---:|---:|---:|---:|---:|",
        ]
    )
    for row in rows:
        lines.append(
            f"| {row['case']} | {row['version']} | {row['trials']} | "
            f"{format_number(row['wall_seconds'], 2)} s | {format_rss(row['peak_rss_kib'])} | "
            f"{format_number(row.get('psms_at_one_percent_fdr'), 0)} | "
            f"{format_number(row.get('peptides_at_one_percent_fdr'), 0)} | "
            f"{format_number(row.get('peptides_in_database'), 0)} |"
        )
    lines.extend(["", "## Threshold checks", ""])
    findings = comparison_findings(rows)
    if findings:
        lines.extend(f"- {finding}" for finding in findings)
    else:
        lines.append("- No complete comparison is available yet")
    feature_records = [record for record in records if record["case"] == "feature"]
    if feature_records:
        latest = feature_records[-1]
        lines.extend(
            [
                "",
                "## Feature summary",
                "",
                "```json",
                json.dumps(
                    {
                        key: latest.get(key)
                        for key in ("models", "quantification", "spectral_library", "library_search")
                        if key in latest
                    },
                    indent=2,
                    sort_keys=True,
                ),
                "```",
            ]
        )
    lines.extend(
        [
            "",
            "## Interpretation",
            "",
            "`REVIEW` is an investigation prompt, not an automatic failure. Timing and RSS changes under five percent are usually noise on a developer workstation.",
            "",
        ]
    )
    return "\n".join(lines)


def make_prefilter_configs(session: BenchmarkSession, source: Path, chunk_size: int) -> tuple[Path, Path]:
    base = read_json(source)
    database = base.get("database")
    if not isinstance(database, dict):
        raise RuntimeError("prefilter benchmark requires a database search configuration")
    paths = []
    for enabled in (False, True):
        config = json.loads(json.dumps(base))
        config["database"]["prefilter"] = enabled
        config["database"]["prefilter_chunk_size"] = chunk_size
        path = session.result_directory / "configs" / f"prefilter-{'on' if enabled else 'off'}.json"
        write_json(path, config)
        paths.append(path)
    return paths[0], paths[1]


def synthetic_sequence(value: int, length: int = 15) -> str:
    residues = "ACDEFGHIKLMNPQRSTVWY"
    sequence = ["A"] * length
    for index in range(length - 1):
        sequence[index] = residues[value % len(residues)]
        value //= len(residues)
    sequence[-1] = "K"
    return "".join(sequence)


def make_memory_config(session: BenchmarkSession, count: int) -> Path:
    generated = WORK_ROOT / "generated"
    generated.mkdir(parents=True, exist_ok=True)
    peptides = generated / f"peptides-{count}.tsv"
    if not peptides.is_file():
        temporary = peptides.with_suffix(".tmp")
        with temporary.open("w") as handle:
            handle.write("sequence\tprotein\n")
            for index in range(count):
                handle.write(f"{synthetic_sequence(index)}\tbenchmark_{index}\n")
        temporary.replace(peptides)
    config = {
        "database": {
            "peptides": str(peptides.resolve()),
            "generate_decoys": True,
            "decoy_tag": "rev_",
            "bucket_size": 8192,
        },
        "deisotope": False,
        "chimera": False,
        "report_psms": 1,
        "output_filter": {"psm_q_value": 1.0},
        "precursor_tol": {"ppm": [-50, 50]},
        "fragment_tol": {"ppm": [-10, 10]},
        "isotope_errors": [-1, 3],
        "mzml_paths": [str((REPO / "tests/LQSRPAAPPAPGPGQLTLR.mzML").resolve())],
        "score_type": "SageHyperScore",
    }
    path = session.result_directory / "configs" / f"memory-{count}.json"
    write_json(path, config)
    return path


def add_common_arguments(parser: argparse.ArgumentParser, include_baseline: bool = True) -> None:
    if include_baseline:
        parser.add_argument("--baseline-ref", default=DEFAULT_BASELINE)
    parser.add_argument("--repeats", type=int, default=3)
    parser.add_argument("--warmups", type=int, default=1)
    parser.add_argument("--threads", type=int, default=8)


def validate_arguments(args: argparse.Namespace) -> None:
    for name in ("repeats", "threads"):
        if hasattr(args, name) and getattr(args, name) < 1:
            raise RuntimeError(f"--{name.replace('_', '-')} must be at least one")
    if hasattr(args, "warmups") and args.warmups < 0:
        raise RuntimeError("--warmups must not be negative")
    if hasattr(args, "memory_peptides") and args.memory_peptides < 1:
        raise RuntimeError("--memory-peptides must be at least one")
    if hasattr(args, "prefilter_chunk_size") and args.prefilter_chunk_size < 1:
        raise RuntimeError("--prefilter-chunk-size must be at least one")
    if hasattr(args, "config"):
        config = Path(args.config).resolve()
        if not config.is_file():
            raise RuntimeError(f"configuration does not exist: {config}")
        read_json(config)


def execute_suite(args: argparse.Namespace) -> Path:
    validate_arguments(args)
    baseline_ref = getattr(args, "baseline_ref", DEFAULT_BASELINE)
    include_baseline = args.command in {"build", "search", "memory", "all"}
    build = build_binaries(baseline_ref, include_baseline=include_baseline)
    if args.command == "build":
        print(f"candidate: {build['candidate_binary']}")
        print(f"baseline:  {build['baseline_binary']}")
        return WORK_ROOT
    session = BenchmarkSession(
        args.command,
        build,
        baseline_ref,
        args.repeats,
        args.warmups,
        args.threads,
    )
    try:
        if args.command in {"search", "all"}:
            config = Path(args.config).resolve()
            session.copy_config("standard", config)
            session.run_search_case(
                "standard", "baseline", build["baseline_binary"], config
            )
            session.run_search_case(
                "standard", "candidate", build["candidate_binary"], config
            )
        if args.command in {"prefilter", "all"}:
            config = Path(args.config).resolve()
            off, on = make_prefilter_configs(session, config, args.prefilter_chunk_size)
            session.run_search_case("prefilter", "off", build["candidate_binary"], off)
            session.run_search_case("prefilter", "on", build["candidate_binary"], on)
        if args.command == "memory":
            config = make_memory_config(session, args.memory_peptides)
            session.run_search_case(
                "memory", "baseline", build["baseline_binary"], config
            )
            session.run_search_case(
                "memory", "candidate", build["candidate_binary"], config
            )
        if args.command == "feature":
            config = Path(args.config).resolve()
            session.copy_config("feature", config)
            session.run_search_case(
                "feature", "candidate", build["candidate_binary"], config
            )
    finally:
        session.save()
    print(session.result_directory)
    return session.result_directory


def clean(confirm: str) -> None:
    if confirm != "yes":
        raise RuntimeError("cleanup requires --confirm yes")
    for path in (WORK_ROOT, RESULTS_ROOT):
        if path.exists():
            shutil.rmtree(path)
            print(f"removed {path}")


def parser() -> argparse.ArgumentParser:
    root = argparse.ArgumentParser(description=__doc__)
    commands = root.add_subparsers(dest="command", required=True)

    build = commands.add_parser("build", help="build candidate and baseline binaries")
    build.add_argument("--baseline-ref", default=DEFAULT_BASELINE)

    search = commands.add_parser("search", help="run the baseline and candidate search")
    search.add_argument("--config", required=True)
    add_common_arguments(search)

    prefilter = commands.add_parser("prefilter", help="compare candidate prefilter modes")
    prefilter.add_argument("--config", required=True)
    prefilter.add_argument("--prefilter-chunk-size", type=int, default=1000)
    add_common_arguments(prefilter, include_baseline=False)

    memory = commands.add_parser("memory", help="run the synthetic peak-memory comparison")
    memory.add_argument("--memory-peptides", type=int, default=1_000_000)
    add_common_arguments(memory)

    feature = commands.add_parser("feature", help="run a candidate-only feature configuration")
    feature.add_argument("--config", required=True)
    add_common_arguments(feature, include_baseline=False)

    all_command = commands.add_parser("all", help="run the real-search benchmark suite")
    all_command.add_argument("--config", required=True)
    all_command.add_argument("--prefilter-chunk-size", type=int, default=1000)
    add_common_arguments(all_command)

    clean_command = commands.add_parser("clean", help="remove generated benchmark data")
    clean_command.add_argument("--confirm", default="no")
    return root


def main() -> int:
    args = parser().parse_args()
    try:
        if args.command == "clean":
            clean(args.confirm)
        else:
            execute_suite(args)
    except (OSError, RuntimeError, subprocess.CalledProcessError) as error:
        print(f"benchmark failed: {error}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
