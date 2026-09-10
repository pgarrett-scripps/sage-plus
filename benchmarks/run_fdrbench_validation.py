#!/usr/bin/env python3
"""Run a sequential, memory-bounded paired entrapment validation suite."""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
import shutil
import subprocess
from pathlib import Path

from provenance import atomic_json, cached_stage, file_identities, sha256


REPO = Path(__file__).resolve().parents[1]
DEFAULT_SUBSET = REPO / "benchmarks/.work/entrapment-small/human-subset.fasta"
DEFAULT_JAR = Path("/tmp/fdrbench-runtime/fdrbench-1.1.1/fdrbench-1.1.1.jar")
DEFAULT_BASELINE = (
    REPO / "benchmarks/.work/targets/baseline-df9219951cc9a54c/release/sage"
)
DEFAULT_CANDIDATE = REPO / "benchmarks/.work/targets/candidate/release/sage"
DEFAULT_MZML = REPO / "data/silac-k6r6/HEK_SILAC-K6R6.mzML"
DEFAULT_OUTPUT = REPO / "benchmarks/results/fdrbench-hardening"
DUCKDB = Path(shutil.which("duckdb") or "duckdb")
JAVA = Path(shutil.which("java") or "java")
GNU_TIME = Path("/usr/bin/time")
PRLIMIT = Path("/usr/bin/prlimit")
GIB = 1024 ** 3
MEMORY_LIMIT_BYTES = 8 * GIB
MIN_AVAILABLE_BYTES = 8 * GIB
MAX_FASTA_ENTRIES = 350_000
Q_LIMITS = (0.01, 0.02, 0.05, 0.10, 0.20)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--subset", type=Path, default=DEFAULT_SUBSET)
    parser.add_argument("--jar", type=Path, default=DEFAULT_JAR)
    parser.add_argument("--baseline", type=Path, default=DEFAULT_BASELINE)
    parser.add_argument("--candidate", type=Path, default=DEFAULT_CANDIDATE)
    parser.add_argument("--mzml", type=Path, default=DEFAULT_MZML)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--first-seed", type=int, default=20260902)
    parser.add_argument("--seed-count", type=int, default=20)
    parser.add_argument("--threads", type=int, default=8)
    parser.add_argument("--baseline-format", choices=["tsv", "parquet"], default="tsv")
    parser.add_argument("--candidate-format", choices=["tsv", "parquet"], default="parquet")
    parser.add_argument("--baseline-label", default="Sage")
    parser.add_argument("--candidate-label", default="Sage Plus")
    parser.add_argument("--duckdb", type=Path, default=DUCKDB)
    parser.add_argument("--java", type=Path, default=JAVA)
    parser.add_argument("--time", type=Path, default=GNU_TIME)
    parser.add_argument("--prlimit", type=Path, default=PRLIMIT)
    parser.add_argument("--preflight-only", action="store_true")
    return parser.parse_args()


def available_memory_bytes() -> int:
    for line in Path("/proc/meminfo").read_text().splitlines():
        if line.startswith("MemAvailable:"):
            return int(line.split()[1]) * 1024
    raise RuntimeError("MemAvailable is missing from /proc/meminfo")


def require_available_memory() -> None:
    available = available_memory_bytes()
    if available < MIN_AVAILABLE_BYTES:
        raise RuntimeError(
            f"only {available / GIB:.2f} GiB is available, "
            f"at least {MIN_AVAILABLE_BYTES / GIB:.2f} GiB is required"
        )


def require_file(path: Path) -> None:
    if not path.is_file():
        raise RuntimeError(f"required file is missing: {path}")


def run_logged(command: list[str], log: Path, env: dict[str, str] | None = None) -> None:
    log.parent.mkdir(parents=True, exist_ok=True)
    with log.open("wb") as handle:
        completed = subprocess.run(
            command,
            cwd=REPO,
            env=env,
            stdout=handle,
            stderr=subprocess.STDOUT,
        )
    if completed.returncode != 0:
        raise RuntimeError(
            f"command failed with exit code {completed.returncode}, see {log}"
        )


def count_fasta_entries(path: Path) -> int:
    count = 0
    with path.open("rb") as handle:
        for line in handle:
            if line.startswith(b">"):
                count += 1
    return count


def generate_database(seed_dir: Path, subset: Path, jar: Path, seed: int) -> tuple[Path, Path]:
    pair = seed_dir / "paired.txt"
    fasta = seed_dir / "paired.fasta"
    command = [
        str(JAVA),
        "-Xmx2G",
        "-jar",
        str(jar),
        "-I2L",
        "-level",
        "peptide",
        "-db",
        str(subset),
        "-o",
        str(pair),
        "-fix_nc",
        "c",
        "-enzyme",
        "1",
        "-miss_c",
        "1",
        "-minLength",
        "7",
        "-maxLength",
        "50",
        "-fold",
        "1",
        "-seed",
        str(seed),
        "-check",
        "-ns",
    ]
    with cached_stage(seed_dir / "database-stage.json", "database-v1", [subset, jar, JAVA], [pair, fasta], command) as hit:
        if not hit:
            require_available_memory()
            run_logged(command, seed_dir / "database-generation.log")
            require_file(pair)
            require_file(fasta)
            entries = count_fasta_entries(fasta)
            if entries > MAX_FASTA_ENTRIES:
                raise RuntimeError(
                    f"database has {entries:,} entries, limit is {MAX_FASTA_ENTRIES:,}"
                )
    return pair, fasta


def write_config(path: Path, fasta: Path, mzml: Path) -> None:
    config = {
        "database": {
            "bucket_size": 16384,
            "enzyme": {
                "missed_cleavages": 0,
                "cleave_at": "$",
                "restrict": "",
                "min_len": 7,
                "max_len": 50,
            },
            "static_mods": {"C": 57.021464},
            "decoy_tag": "rev_",
            "generate_decoys": True,
            "fasta": str(fasta),
        },
        "precursor_tol": {"ppm": [-10.0, 10.0]},
        "fragment_tol": {"ppm": [-20.0, 20.0]},
        "isotope_errors": [-1, 3],
        "deisotope": True,
        "chimera": False,
        "max_fragment_charge": 2,
        "min_matched_peaks": 6,
        "report_psms": 1,
        "output_filter": {"psm_q_value": 0.2},
        "max_memory_gb": 8.0,
        "min_free_memory_gb": 8.0,
        "batch_size": 1,
        "mzml_paths": [str(mzml)],
        "score_type": "SageHyperScore",
    }
    path.write_text(json.dumps(config, indent=2) + "\n")


def run_search(
    binary: Path,
    config: Path,
    output: Path,
    timing: Path,
    log: Path,
    threads: int,
    output_format: str,
) -> None:
    if output_format not in {"tsv", "parquet"}:
        raise ValueError("output_format must be tsv or parquet")
    result_name = f"results.sage.{output_format}"
    output.mkdir(parents=True, exist_ok=True)
    environment = os.environ.copy()
    environment["RAYON_NUM_THREADS"] = str(threads)
    command = [
        str(GNU_TIME),
        "-f",
        "%e\t%M\t%x",
        "-o",
        str(timing),
        str(PRLIMIT),
        f"--as={MEMORY_LIMIT_BYTES}:{MEMORY_LIMIT_BYTES}",
        str(binary),
        str(config),
        "--output_directory",
        str(output),
        "--batch-size",
        "1",
        "--disable-telemetry-i-dont-want-to-improve-sage",
    ]
    settings = json.loads(config.read_text())
    database = settings["database"]
    inputs = [binary, config, GNU_TIME, PRLIMIT]
    inputs.extend(Path(value) for key in ("fasta", "peptides", "custom_cleavage_sites") if (value := database.get(key)))
    if database.get("ptm_library"):
        inputs.append(Path(database["ptm_library"]["path"]))
    inputs.extend(Path(path) for path in settings["mzml_paths"])
    outputs = [output / result_name, timing, output / "results.json"]
    if output_format == "parquet":
        outputs.append(output / "run-summary.json")
    with cached_stage(output.parent / "search-stage.json", "search-v1", inputs,
                      outputs, {"command": command, "threads": threads, "format": output_format}) as hit:
        if not hit:
            require_available_memory()
            run_logged(command, log, environment)
            require_file(output / result_name)
        timing_values(timing)


def sql_quote(path: Path) -> str:
    return str(path).replace("'", "''")


def extract_peptides(source: Path, target: Path, candidate: bool) -> None:
    source_sql = sql_quote(source)
    target_sql = sql_quote(target)
    if candidate:
        relation = f"read_parquet('{source_sql}')"
        target_filter = "NOT is_decoy"
    else:
        relation = f"read_csv('{source_sql}', delim='\\t', header=true)"
        target_filter = "label = 1"
    query = (
        "COPY ("
        "SELECT proteins AS peptide, min(peptide_q) AS q_value, "
        "max(hyperscore) AS score, any_value(proteins) AS protein "
        f"FROM {relation} "
        f"WHERE rank = 1 AND {target_filter} AND peptide_q <= 0.2 "
        "GROUP BY proteins"
        f") TO '{target_sql}' (HEADER, DELIMITER '\\t')"
    )
    with cached_stage(target.with_suffix(".stage.json"), "peptide-extraction-v1", [source, DUCKDB], [target],
                      {"query": query, "parquet": candidate}) as hit:
        if not hit:
            run_logged([str(DUCKDB), "-c", query], target.with_suffix(".duckdb.log"))


def calculate_fdp(source: Path, pair: Path, output: Path, jar: Path) -> None:
    command = [
        str(JAVA),
        "-Xmx1G",
        "-jar",
        str(jar),
        "-i",
        str(source),
        "-level",
        "peptide",
        "-pep",
        str(pair),
        "-score",
        "score:1",
        "-o",
        str(output),
    ]
    with cached_stage(output.with_suffix(".stage.json"), "fdp-v1", [source, pair, jar, JAVA], [output], command) as hit:
        if not hit:
            run_logged(command, output.with_suffix(".log"))


def read_fdp(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def point_at_q(rows: list[dict[str, str]], limit: float) -> dict[str, float] | None:
    eligible = [row for row in rows if float(row["q_value"]) <= limit]
    if not eligible:
        return None
    row = max(
        eligible,
        key=lambda value: (float(value["q_value"]), -float(value["score"])),
    )
    return {
        "observed_q": float(row["q_value"]),
        "targets": int(row["n_t"]),
        "entrapments": int(row["n_p"]),
        "paired_fdp": float(row["paired_fdp"]),
        "combined_fdp": float(row["combined_fdp"]),
        "lower_bound_fdp": float(row["lower_bound_fdp"]),
    }


def matched_yield(rows: list[dict[str, str]], limit: float) -> dict[str, float] | None:
    eligible = [row for row in rows if float(row["paired_fdp"]) <= limit]
    if not eligible:
        return None
    row = max(eligible, key=lambda value: (int(value["n_t"]), float(value["q_value"])))
    return {
        "reported_q": float(row["q_value"]),
        "targets": int(row["n_t"]),
        "entrapments": int(row["n_p"]),
        "paired_fdp": float(row["paired_fdp"]),
    }


def timing_values(path: Path) -> dict[str, float]:
    fields = path.read_text().strip().split("\t")
    if len(fields) != 3 or fields[2] != "0":
        raise RuntimeError(f"unexpected timing output in {path}")
    wall, rss = float(fields[0]), int(fields[1])
    if not math.isfinite(wall) or wall < 0 or rss <= 0:
        raise RuntimeError(f"invalid resource measurement in {path}")
    return {"wall_seconds": wall, "peak_rss_mib": rss / 1024}


def process_seed(args: argparse.Namespace, seed: int) -> dict[str, object]:
    seed_dir = args.output / "seeds" / str(seed)
    summary_path = seed_dir / "summary.json"
    seed_dir.mkdir(parents=True, exist_ok=True)
    pair, fasta = generate_database(seed_dir, args.subset, args.jar, seed)
    config = seed_dir / "search-config.json"
    write_config(config, fasta, args.mzml)
    engines = {
        args.baseline_label: (args.baseline, args.baseline_format),
        args.candidate_label: (args.candidate, args.candidate_format),
    }
    summary: dict[str, object] = {
        "seed": seed,
        "database_entries": count_fasta_entries(fasta),
        "database_sha256": sha256(fasta),
        "pair_sha256": sha256(pair),
        "engines": {},
    }
    for engine, (binary, output_format) in engines.items():
        candidate = output_format == "parquet"
        slug = "sage-plus" if engine == args.candidate_label else "sage"
        engine_dir = seed_dir / slug
        search_output = engine_dir / "search-output"
        timing = engine_dir / "timing.tsv"
        run_search(
            binary,
            config,
            search_output,
            timing,
            engine_dir / "search.log",
            args.threads,
            output_format,
        )
        source_name = "results.sage.parquet" if candidate else "results.sage.tsv"
        peptides = engine_dir / "peptides.tsv"
        fdp = engine_dir / "fdp.csv"
        extract_peptides(search_output / source_name, peptides, candidate)
        calculate_fdp(peptides, pair, fdp, args.jar)
        rows = read_fdp(fdp)
        summary["engines"][engine] = {
            "timing": timing_values(timing),
            "q_points": {
                str(limit): point_at_q(rows, limit) for limit in Q_LIMITS
            },
            "matched_yield": {
                str(limit): matched_yield(rows, limit) for limit in Q_LIMITS
            },
            "fdp_sha256": sha256(fdp),
        }
    summary["status"] = "completed"
    atomic_json(summary_path, summary)
    return summary


def write_aggregate(
    output: Path,
    summaries: list[dict[str, object]],
    args: argparse.Namespace,
) -> None:
    threshold_path = output / "threshold-points.csv"
    with threshold_path.open("w", newline="") as handle:
        fields = [
            "seed",
            "engine",
            "q_limit",
            "observed_q",
            "targets",
            "entrapments",
            "paired_fdp",
            "combined_fdp",
            "lower_bound_fdp",
        ]
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for summary in summaries:
            for engine, values in summary["engines"].items():
                for limit, point in values["q_points"].items():
                    if point is not None:
                        writer.writerow(
                            {
                                "seed": summary["seed"],
                                "engine": engine,
                                "q_limit": limit,
                                **point,
                            }
                        )
    yield_path = output / "matched-yield.csv"
    with yield_path.open("w", newline="") as handle:
        fields = [
            "seed",
            "engine",
            "fdp_limit",
            "reported_q",
            "targets",
            "entrapments",
            "paired_fdp",
        ]
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for summary in summaries:
            for engine, values in summary["engines"].items():
                for limit, point in values["matched_yield"].items():
                    if point is not None:
                        writer.writerow(
                            {
                                "seed": summary["seed"],
                                "engine": engine,
                                "fdp_limit": limit,
                                **point,
                            }
                        )
    curves_path = output / "calibration-curves.csv"
    with curves_path.open("w", newline="") as handle:
        fields = [
            "seed",
            "engine",
            "q_value",
            "targets",
            "entrapments",
            "paired_fdp",
            "combined_fdp",
            "lower_bound_fdp",
        ]
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for summary in summaries:
            seed = summary["seed"]
            for engine in summary["engines"]:
                slug = "sage-plus" if engine == args.candidate_label else "sage"
                rows = read_fdp(output / "seeds" / str(seed) / slug / "fdp.csv")
                for row in rows:
                    writer.writerow(
                        {
                            "seed": seed,
                            "engine": engine,
                            "q_value": row["q_value"],
                            "targets": row["n_t"],
                            "entrapments": row["n_p"],
                            "paired_fdp": row["paired_fdp"],
                            "combined_fdp": row["combined_fdp"],
                            "lower_bound_fdp": row["lower_bound_fdp"],
                        }
                    )
    manifest = {
        "schema_version": 2,
        "engine_labels": [args.baseline_label, args.candidate_label],
        "summary_hashes": file_identities([output / "seeds" / str(s["seed"]) / "summary.json" for s in summaries]),
        "binaries": file_identities([args.baseline, args.candidate]),
        "input_mzml_sha256": sha256(args.mzml),
        "seed_count": len(summaries),
        "seeds": [summary["seed"] for summary in summaries],
        "subset": str(args.subset.resolve()),
        "subset_sha256": sha256(args.subset),
        "fdrbench_jar_sha256": sha256(args.jar),
        "memory_limit_gib": MEMORY_LIMIT_BYTES / GIB,
        "minimum_available_memory_gib": MIN_AVAILABLE_BYTES / GIB,
        "max_fasta_entries": MAX_FASTA_ENTRIES,
        "sequential": True,
    }
    atomic_json(output / "manifest.json", manifest)


def main() -> int:
    global DUCKDB, JAVA, GNU_TIME, PRLIMIT
    args = parse_args()
    for name in ("subset", "jar", "baseline", "candidate", "mzml", "output", "duckdb", "java", "time", "prlimit"):
        setattr(args, name, getattr(args, name).resolve())
    DUCKDB, JAVA, GNU_TIME, PRLIMIT = args.duckdb, args.java, args.time, args.prlimit
    if args.threads < 1:
        raise RuntimeError("threads must be positive")
    if not args.baseline_label or not args.candidate_label or args.baseline_label == args.candidate_label:
        raise RuntimeError("engine labels must be nonempty and distinct")
    if args.seed_count < 1:
        raise RuntimeError("seed count must be positive")
    for path in (
        args.subset,
        args.jar,
        args.baseline,
        args.candidate,
        args.mzml,
        DUCKDB,
        JAVA,
        GNU_TIME,
        PRLIMIT,
    ):
        require_file(path)
    for path in (args.baseline, args.candidate, DUCKDB, JAVA, GNU_TIME, PRLIMIT):
        if not os.access(path, os.X_OK):
            raise RuntimeError(f"required tool is not executable: {path}")
    if args.preflight_only:
        print("All required inputs and tools are accessible")
        return 0
    args.output.mkdir(parents=True, exist_ok=True)
    # A new invocation must complete its requested matrix before it can be analyzed.
    (args.output / "manifest.json").unlink(missing_ok=True)
    seeds = range(args.first_seed, args.first_seed + args.seed_count)
    summaries = []
    for index, seed in enumerate(seeds, start=1):
        print(f"seed {index}/{args.seed_count}: {seed}", flush=True)
        summaries.append(process_seed(args, seed))
    write_aggregate(args.output, summaries, args)
    print(f"wrote {args.output}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
