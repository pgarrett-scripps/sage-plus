#!/usr/bin/env python3
"""Evaluate development fixes without modifying the frozen pilot evidence."""

import argparse
import csv
import json
import math
from collections import defaultdict
from pathlib import Path

from generate_ptm_library_benchmark import fasta_entries
from provenance import atomic_json, file_identities
from scientific_metrics import (
    boolean, compare_identifications, localization_metrics, normalize_psms, read_table, stripped,
)


def verified_rows(path):
    record = json.loads((path.parent / "result.json").read_text())
    if record["status"] != "complete":
        raise ValueError(f"Incomplete search: {path.parent}")
    actual = file_identities([path])
    if any(record["outputs"].get(name) != digest for name, digest in actual.items()):
        raise ValueError(f"Search output changed: {path}")
    return read_table(path)


def transfer_q_values(rows):
    """Experimental count estimates within recipient files and transfer candidates."""
    groups = defaultdict(list)
    result = {}
    for index, row in enumerate(rows):
        if not row.get("transfer_candidate") or not boolean(row["transfer_candidate"]):
            continue
        score = float(row["file_score"])
        intensity = float(row["intensity"])
        if not math.isfinite(score) or not math.isfinite(intensity) or intensity <= 0:
            raise ValueError("Transfer candidate has invalid signal evidence")
        groups[row["filename"]].append((score, boolean(row["is_decoy"]), index))
    for group in groups.values():
        group.sort(key=lambda row: -row[0])
        decoys, targets, start = 1, 0, 0
        thresholds = []
        while start < len(group):
            end = start + 1
            while end < len(group) and group[end][0] == group[start][0]:
                end += 1
            decoys += sum(row[1] for row in group[start:end])
            targets += sum(not row[1] for row in group[start:end])
            thresholds.append((start, end, min(1.0, decoys / targets) if targets else 1.0))
            start = end
        minimum = 1.0
        for start, end, q in reversed(thresholds):
            minimum = min(minimum, q)
            for _, _, index in group[start:end]:
                result[index] = minimum
    return result


def control_metrics(rows, estimates, references):
    memberships = {}
    result = []
    for endpoint in ("legacy_precursor_1pct", "experimental_transfer_1pct", "experimental_transfer_5pct"):
        quantified = foreign = ambiguous = 0
        for index, row in enumerate(rows):
            if "DDA_Human_" not in row["filename"] or boolean(row["is_decoy"]):
                continue
            if not row.get("intensity") or float(row["intensity"]) <= 0:
                continue
            if endpoint == "legacy_precursor_1pct":
                accepted = float(row["q_value"]) <= 0.01 and not boolean(row["ms2_confirmed_strict"])
            else:
                cutoff = 0.01 if endpoint.endswith("1pct") else 0.05
                accepted = index in estimates and estimates[index] <= cutoff
            if not accepted:
                continue
            sequence = stripped(row["peptide"]).replace("I", "L")
            if sequence not in memberships:
                memberships[sequence] = {species for species, reference in references.items() if sequence in reference}
            membership = memberships[sequence]
            if len(membership) != 1:
                ambiguous += 1
                continue
            quantified += 1
            foreign += bool(membership & {"yeast", "ecoli"})
        fraction = foreign / quantified if quantified else None
        result.append({"endpoint": endpoint, "unconfirmed_control_rows": quantified,
                       "foreign_rows": foreign, "foreign_fraction": fraction,
                       "ambiguous_or_unmapped_excluded": ambiguous,
                       "interpretation": "No calibration claim. Foreign assignments are a diagnostic under purity and reference assumptions."})
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("root", type=Path)
    parser.add_argument("--pilot", type=Path, required=True)
    args = parser.parse_args()
    ptm = []
    inputs = [Path(__file__).resolve(), Path(__file__).with_name("SCIENTIFIC_HARDENING.md")]
    truth_path = args.pilot / "ptm-truth/truth.tsv"
    inputs.append(truth_path)
    with truth_path.open() as handle:
        truth = {row["sequence"]: {int(p) for p in row["positions"].split(",")}
                 for row in csv.DictReader(handle, delimiter="\t") if row["unambiguous_site"] == "True"}
    ptm_paths = list((args.root / "runs").glob("ptm-*/results.sage.ptm-sites.parquet"))
    ptm_paths += list((args.root / "runs-ptm-ties").glob("ptm-*/results.sage.ptm-sites.parquet"))
    for path in sorted(ptm_paths):
        rows = verified_rows(path)
        inputs.append(path)
        mapping = {(Path(row["filename"]).name, sequence): positions
                   for row in rows[:1] for sequence, positions in truth.items()}
        accepted = [row for row in rows if float(row["peptide_q"]) <= 0.01]
        ptm.append({"job": path.parent.name, "site_rows": len(rows),
                    "joint_1pct_synthesis_consistency": localization_metrics(accepted, mapping),
                    "peptide_q_one_rows": sum(float(row["peptide_q"]) == 1 for row in rows),
                    "confidence_fallback_logged": "using target-decoy counts with a +1 correction" in (path.parent / "search.log").read_text()})
    references = {}
    for species in ("human", "yeast", "ecoli"):
        path = args.pilot / f"references/{species}.fasta"
        inputs.append(path)
        references[species] = "\x00".join(sequence.replace("I", "L") for _, sequence in fasta_entries(path))
    lfq = []
    for path in sorted((args.root / "runs").glob("mbr-*/lfq.parquet")):
        inputs.append(path)
        rows = verified_rows(path)
        estimates = transfer_q_values(rows)
        output = path.parent / "experimental-transfer-q.tsv"
        with output.open("w") as handle:
            writer = csv.writer(handle, delimiter="\t")
            writer.writerow(["peptide", "charge", "filename", "is_decoy", "file_score", "experimental_transfer_q"])
            for index, q in estimates.items():
                row = rows[index]
                writer.writerow([row["peptide"], row["charge"], row["filename"], row["is_decoy"], row["file_score"], q])
        by_file = {}
        for name in sorted({row["filename"] for row in rows}):
            indices = [i for i in estimates if rows[i]["filename"] == name]
            by_file[name] = {"target_candidates": sum(not boolean(rows[i]["is_decoy"]) for i in indices),
                             "decoy_candidates": sum(boolean(rows[i]["is_decoy"]) for i in indices),
                             "target_transfers_1pct": sum(not boolean(rows[i]["is_decoy"]) and estimates[i] <= 0.01 for i in indices),
                             "target_transfers_5pct": sum(not boolean(rows[i]["is_decoy"]) and estimates[i] <= 0.05 for i in indices)}
        lfq.append({"job": path.parent.name, "by_file": by_file,
                    "pure_human": control_metrics(rows, estimates, references),
                    "estimate_table": file_identities([output])})
    regression = []
    for path in sorted((args.root / "runs-ptm-ties").glob("hek-regression-*/results.sage.parquet")):
        index = path.parent.name.rsplit("-", 1)[1]
        baseline = args.pilot / f"runs/public-comparison/PXD001468-{index}-plus/results.sage.parquet"
        inputs.extend([path, baseline])
        comparison = compare_identifications(normalize_psms(verified_rows(baseline)),
                                             normalize_psms(verified_rows(path)))
        regression.append({"job": path.parent.name, **comparison})
    result = {"status": "development_evaluation", "inputs": file_identities(inputs),
              "ptm": ptm, "lfq": lfq,
              "hek_regression": regression,
              "scope": "Previously inspected pilot data. File scores and transfer estimates remain experimental. PTM truth mapping remains unaudited."}
    atomic_json(args.root / "evaluation.json", result)


if __name__ == "__main__":
    main()
