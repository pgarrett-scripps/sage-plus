#!/usr/bin/env python3
"""Reduce the mass offset matrix to the retained summary the chapter reads.

Reads the job records written by `run_mass_offset.py` and derives cost against
the number of configured offsets, agreement between offset and database-expanded
searches, localization against synthesis-defined phosphosites, and peptide-level
entrapment error. Parquet is read through the same helper as the other
scientific analyses.
"""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
from statistics import median

from provenance import atomic_json, file_identities, sha256
from scientific_metrics import localization_metrics, read_table

REPO = Path(__file__).resolve().parents[1]
PHOSPHO = 79.966331


def rank_one(path: Path) -> dict[tuple[str, str], dict]:
    return {(row["filename"], row["scannr"]): row
            for row in read_table(path) if int(row["rank"]) == 1}


def agreement(indexed: Path, offset: Path) -> dict:
    left, right = rank_one(indexed), rank_one(offset)
    shared = set(left) & set(right)
    same = sum(left[key]["peptide"] == right[key]["peptide"] for key in shared)
    same_score = sum(
        left[key]["peptide"] == right[key]["peptide"]
        and abs(float(left[key]["hyperscore"]) - float(right[key]["hyperscore"])) < 1e-4
        for key in shared)
    accepted = lambda rows: {k for k, r in rows.items()
                             if float(r["spectrum_q"]) <= 0.01 and r["is_decoy"] in ("false", "0", False)}
    return {"shared_spectra": len(shared), "same_peptidoform": same,
            "same_peptidoform_and_score": same_score,
            "accepted_indexed": len(accepted(left)), "accepted_offset": len(accepted(right))}


def localization(job: Path, truth: dict[str, set[int]]) -> dict:
    rows = [row for row in read_table(job / "results.sage.ptm-sites.parquet")
            if float(row["peptide_q"]) <= 0.01]
    names = {row["filename"] for row in rows}
    mapping = {(name, sequence): positions for name in names for sequence, positions in truth.items()}
    metrics = localization_metrics(rows, mapping)
    return {"correct_site_events": metrics["correct_site_events"],
            "incorrect_site_events": metrics["incorrect_site_events"],
            "site_error_fraction": metrics["empirical_site_error_fraction"],
            "scope": "site events at 1% spectrum, peptide, and localization q, against "
                     "synthesis-defined sites. Includes identification error"}


def entrapment(job: Path, labels: dict[str, str]) -> list[dict]:
    best: dict[str, float] = {}
    for row in read_table(job / "results.sage.parquet"):
        if row["is_decoy"] in ("true", "1", True):
            continue
        sequence, q = row["stripped_peptide"], float(row["peptide_q"])
        best[sequence] = min(best.get(sequence, 1.0), q)
    out = []
    for threshold in (0.01, 0.05):
        selected = [s for s, q in best.items() if q <= threshold]
        unknown = [s for s in selected if s not in labels]
        entraps = sum(labels[s] == "p_target" for s in selected if s in labels)
        total = len(selected) - len(unknown)
        out.append({"nominal_q": threshold, "targets": total - entraps, "entrapments": entraps,
                    "combined_fdp": 2 * entraps / total if total else None,
                    "lower_bound_fdp": entraps / total if total else None,
                    "peptides_outside_pairing": len(unknown)})
    return out


def scale(records: list[dict]) -> list[dict]:
    rows = []
    for name in ("indexed", "offsets-1", "offsets-2", "offsets-3"):
        group = [r for r in records if r["suite"] == "scale" and r["id"].startswith(name + "-r")]
        assert group and all(r["exit_status"] == 0 for r in group), name
        summary = group[0]["summary"]
        rows.append({
            "configuration": name, "trials": len(group),
            "search_stage_seconds": median(r["search_stage_ms"] for r in group) / 1000,
            "wall_seconds": median(r["wall_seconds"] for r in group),
            "peak_rss_mib": median(r["max_rss_bytes"] for r in group) / 1024**2,
            "database_peptides": summary["peptides_in_database"],
            "database_fragments": summary["fragments_in_database"],
            "psms": summary["psms_at_one_percent_fdr"],
            "peptides": summary["peptides_at_one_percent_fdr"],
            "offset_psms": summary["modifications"]["mass_offset_psms"],
        })
    return rows


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--truth", type=Path,
                        default=Path("/data/sage-plus-scientific/20260914/ptm-truth/truth.tsv"))
    parser.add_argument("--pairs", type=Path, default=None)
    parser.add_argument("--output", type=Path,
                        default=REPO / "benchmarks/scientific-results/mass-offset-20260919/summary.json")
    args = parser.parse_args()
    pairs = args.pairs or args.root / "inputs/paired.txt"

    matrix = json.loads((args.root / "matrix.json").read_text())
    records = matrix["jobs"]
    assert all(r["exit_status"] == 0 for r in records), "matrix contains failed jobs"
    verified = {}
    for record in records:
        for name, expected in (record["inputs"] | record["config"] | record["outputs"]).items():
            if name not in verified:
                verified[name] = sha256(Path(name))
            assert verified[name] == expected, f"Changed mass-offset evidence: {name}"

    with args.truth.open() as handle:
        truth = {row["sequence"]: {int(p) for p in row["positions"].split(",")}
                 for row in csv.DictReader(handle, delimiter="\t")
                 if row["unambiguous_site"] == "True"}
    with pairs.open() as handle:
        labels = {row["sequence"]: row["peptide_type"] for row in csv.DictReader(handle, delimiter="\t")
                  if row["peptide_type"] in ("target", "p_target")}

    localization_rows = []
    for spectra in ("HCD_1", "HCD_2"):
        jobs = {mode: args.root / "localization" / f"{spectra}-{mode}" for mode in ("indexed", "offset")}
        row = {"spectra": spectra,
               "agreement": agreement(jobs["indexed"] / "results.sage.parquet",
                                      jobs["offset"] / "results.sage.parquet")}
        for mode, job in jobs.items():
            row[mode] = localization(job, truth)
            row[mode]["database_peptides"] = json.loads((job / "job.json").read_text())["summary"]["peptides_in_database"]
        localization_rows.append(row)

    summary = {
        "schema_version": 1,
        "scope": "Mass offset search against the equivalent database expansion. "
                 "Engineering evidence on public spectra. Not an independent biological study.",
        "executable": matrix["executable"],
        "environment": matrix["environment"],
        "evidence_root": str(args.root),
        "verified_evidence_sha256": verified,
        "inputs_sha256": file_identities([
            args.root / "matrix.json", args.truth, pairs,
            Path(__file__).resolve(), Path(__file__).with_name("scientific_metrics.py")]),
        "scale": scale(records),
        "localization": localization_rows,
        "entrapment": {mode: entrapment(args.root / "entrapment" / mode, labels)
                       for mode in ("indexed", "offset")},
    }
    atomic_json(args.output, summary)
    print(f"wrote {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
