#!/usr/bin/env python3
"""Full-reference entrapment pilot with an independent threshold count audit."""

import argparse
import csv
import json
import math
import subprocess
from pathlib import Path

from provenance import atomic_json, cached_stage, file_identities
from run_scientific import base_config, run_job
from scientific_metrics import THRESHOLDS, normalize_psms, peptide_representatives, read_table


JAR = Path("/tmp/fdrbench-runtime/fdrbench-1.1.1/fdrbench-1.1.1.jar")
JAVA = Path("/usr/bin/java").resolve()


def load_pairs(path):
    labels, groups = {}, {}
    with path.open() as f:
        for row in csv.DictReader(f, delimiter="\t"):
            if row["peptide_type"] not in ("target", "p_target"):
                continue
            sequence, kind, pair = row["sequence"], row["peptide_type"], row["peptide_pair_index"]
            if sequence in labels:
                raise ValueError("Duplicate target or entrapment sequence")
            labels[sequence] = kind
            group = groups.setdefault(pair, {})
            if kind in group:
                raise ValueError("Nonunique pair member")
            group[kind] = sequence
    if any(set(g) != {"target", "p_target"} for g in groups.values()):
        raise ValueError("Incomplete target-entrapment pairing")
    partners = {g["p_target"]: g["target"] for g in groups.values()}
    return labels, partners


def threshold_fdp(rows, labels, partners, threshold):
    selected = {r["sequence"]: r for r in rows if r["peptide_q"] <= threshold}
    if len(selected) != sum(r["peptide_q"] <= threshold for r in rows):
        raise ValueError("Peptide sequence is not unique")
    if any(s not in labels for s in selected):
        raise ValueError("An identified peptide is missing from the pairing file")
    nt = sum(labels[s] == "target" for s in selected)
    ne = len(selected) - nt
    absent = earlier = tied = 0
    for seq, row in selected.items():
        if labels[seq] != "p_target":
            continue
        partner = selected.get(partners[seq])
        if partner is None:
            absent += 1
            continue
        a, b = (row["peptide_q"], -row["score"]), (partner["peptide_q"], -partner["score"])
        earlier += a < b
        tied += a == b
    n = nt + ne
    return {"nominal_q": threshold, "targets": nt, "entrapments": ne,
            "combined_fdp": 2 * ne / n if n else None,
            "lower_bound_fdp": ne / n if n else None,
            "paired_fdp_tie_min": (ne + 2 * earlier + absent) / n if n else None,
            "paired_fdp_tie_max": (ne + 2 * (earlier + tied) + absent) / n if n else None,
            "paired_score_ties": tied, "denominator": n}


def generate(reference, output, seed):
    output.mkdir(parents=True, exist_ok=True)
    pair, fasta = output / "paired.txt", output / "paired.fasta"
    command = [str(JAVA), "-XX:ActiveProcessorCount=8", "-Xmx3G", "-jar", str(JAR),
               "-I2L", "-level", "peptide", "-db", str(reference), "-o", str(pair),
               "-fix_nc", "c", "-enzyme", "1", "-miss_c", "1", "-minLength", "7",
               "-maxLength", "50", "-fold", "1", "-seed", str(seed), "-check", "-ns"]
    with cached_stage(output / "generation.json", "full-reference-entrapment-v1",
                      [reference, JAR, JAVA], [pair, fasta], command) as hit:
        if not hit:
            with (output / "generation.log").open("w") as log:
                subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, check=True, timeout=1800)
    return pair, fasta


def evaluate(source, pair, output):
    rows = [r for r in peptide_representatives(normalize_psms(read_table(source))) if not r["is_decoy"]]
    labels, partners = load_pairs(pair)
    if len({r["sequence"] for r in rows}) != len(rows):
        raise ValueError("This pilot supports one static peptidoform per sequence only")
    extracted = output / "fdrbench-input.tsv"
    with extracted.open("w") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["peptide", "q_value", "score", "protein"])
        for r in sorted(rows, key=lambda r: (r["peptide_q"], -r["score"], r["sequence"])):
            w.writerow([r["sequence"], r["peptide_q"], r["score"], r["proteins"]])
    independent = [threshold_fdp(rows, labels, partners, q) for q in THRESHOLDS]
    official = output / "fdrbench.csv"
    command = [str(JAVA), "-XX:ActiveProcessorCount=2", "-Xmx2G", "-jar", str(JAR),
               "-i", str(extracted), "-level", "peptide", "-pep", str(pair),
               "-score", "score:1", "-o", str(official)]
    with (output / "fdrbench.log").open("w") as log:
        subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, check=True, timeout=1200)
    with official.open() as f:
        reference_rows = list(csv.DictReader(f))
    for point in independent:
        eligible = [r for r in reference_rows if float(r["q_value"]) <= point["nominal_q"]]
        if not eligible:
            if point["denominator"]:
                raise ValueError("Missing official threshold output")
            continue
        row = max(eligible, key=lambda r: (float(r["q_value"]), -float(r["score"])))
        if int(row["n_t"]) != point["targets"] or int(row["n_p"]) != point["entrapments"]:
            raise ValueError("Independent counts disagree with FDRBench")
        if not math.isclose(float(row["combined_fdp"]), point["combined_fdp"], abs_tol=1e-10):
            raise ValueError("Combined estimator disagrees with FDRBench")
        paired = float(row["paired_fdp"])
        if not point["paired_fdp_tie_min"] - 1e-10 <= paired <= point["paired_fdp_tie_max"] + 1e-10:
            raise ValueError("Paired estimator disagrees with audited tie interval")
        point["fdrbench_paired_fdp"] = paired
    atomic_json(output / "calibration.json", {
        "status": "independently_checked", "thresholds": independent,
        "scope": "Conditional pilot FDP, not a production FDR certification",
        "target_sequences": sum(v == "target" for v in labels.values()),
        "inputs": file_identities([source, pair, JAR, JAVA, Path(__file__).resolve()]),
        "outputs": file_identities([extracted, official])})


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--reference", type=Path, required=True)
    p.add_argument("--spectra", nargs="+", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    p.add_argument("--seed", type=int, default=20260914)
    p.add_argument("--generate-only", action="store_true")
    args = p.parse_args()
    pair, fasta = generate(args.reference.resolve(), args.output, args.seed)
    if args.generate_only:
        return
    for index, spectra in enumerate(args.spectra):
        config = base_config(fasta, [spectra.resolve()])
        config["database"]["enzyme"].update(cleave_at="$", restrict="", missed_cleavages=0)
        for engine in ("upstream", "plus"):
            job = {"id": f"file-{index}-{engine}", "engine": engine, "threads": 8,
                   "seed": args.seed, "config": config}
            result = run_job(job, args.output, 900)
            if result["status"] != "complete":
                continue
            out = args.output / job["id"]
            source = out / ("results.sage.parquet" if engine == "plus" else "results.sage.tsv")
            evaluate(source, pair, out)


if __name__ == "__main__":
    main()
