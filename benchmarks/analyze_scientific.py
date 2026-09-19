#!/usr/bin/env python3
"""Summarize the pilot without hiding missing cells or pooling pseudo-replicates."""

import argparse
import csv
import json
import random
import statistics
from collections import defaultdict
from functools import partial
from pathlib import Path

from provenance import atomic_json, file_identities
from scientific_diagnostics import direct_ratio_diagnostic, lfq_threshold_yields
from scientific_metrics import (compare_identifications, localization_metrics,
                                normalize_psms, quantification_metrics, read_table)


def completed_output(path):
    manifest = path.parent / "result.json"
    return manifest.exists() and json.loads(manifest.read_text())["status"] == "complete"


def summarize_jobs(root):
    jobs = []
    for path in sorted((root / "runs").glob("*/*/result.json")):
        record = json.loads(path.read_text())
        job = record["signature"]["job"]
        row = {"suite": path.parent.parent.name, "id": job["id"], "engine": job["engine"],
               "workload": job.get("workload"), "study": job.get("study"),
               "warmup": job.get("warmup", False), "trial": job.get("trial"),
               "threads": job["threads"], "status": record["status"],
               "wall_seconds": record.get("wall_seconds"), "peak_rss_mib": record.get("peak_rss_mib"),
               "result_path": str(path.relative_to(root))}
        if record["status"] == "complete":
            metrics = json.loads((path.parent / "identification-metrics.json").read_text())
            row.update(metrics["0.01"])
        jobs.append(row)
    return jobs


def compare_local(root):
    comparisons = []
    directory = root / "runs/local-paired"
    for workload in ("standard", "common-mods", "broad-ptm"):
        for trial in range(1, 4):
            paths = [directory / f"{workload}-{e}-{trial}" for e in ("upstream", "plus")]
            if not all((p / "result.json").exists() and json.loads((p / "result.json").read_text())["status"] == "complete" for p in paths):
                continue
            a = normalize_psms(read_table(paths[0] / "results.sage.tsv"))
            b = normalize_psms(read_table(paths[1] / "results.sage.parquet"))
            comparisons.append({"workload": workload, "trial": trial, **compare_identifications(a, b)})
    return comparisons


def check_prefilter(root):
    directory = root / "runs/scaling"
    comparisons = []
    for threads in (1, 2, 4, 8):
        for trial in range(1, 4):
            paths = [directory / f"threads-{threads}-prefilter-{mode}-{trial}" for mode in (0, 1)]
            if not all((p / "result.json").exists() and json.loads((p / "result.json").read_text())["status"] == "complete" for p in paths):
                continue
            a, b = [normalize_psms(read_table(p / "results.sage.parquet")) for p in paths]
            comparison = compare_identifications(a, b)
            comparison["all_scores_and_qvalues_equal"] = {
                r["key"]: (r["score"], r["hyperscore"], r["spectrum_q"], r["peptide_q"]) for r in a
            } == {r["key"]: (r["score"], r["hyperscore"], r["spectrum_q"], r["peptide_q"]) for r in b}
            comparisons.append({"threads": threads, "trial": trial, **comparison})
    return comparisons


def compare_public(root, suite=None):
    """Compare matched public inputs and expose every missing engine pair."""
    if suite is None:
        suites = sorted(p.name for p in (root / "runs").glob("public-comparison*")
                        if (p / "plan.json").exists())
        if suites:
            return [row for name in suites for row in compare_public(root, name)]
        suite = "public-comparison"
    plan = root / f"{suite}-plan.json"
    if not plan.exists():
        return []
    groups = defaultdict(dict)
    for job in json.loads(plan.read_text())["jobs"]:
        pair = job["id"].rsplit("-", 1)[0]
        if job["engine"] in groups[pair]:
            raise ValueError("Duplicate public comparison engine")
        groups[pair][job["engine"]] = job
    comparisons = []
    for pair, engines in sorted(groups.items()):
        paths = {engine: root / "runs" / suite / job["id"]
                 for engine, job in engines.items()}
        states = {}
        for engine in ("upstream", "plus"):
            result = paths.get(engine, root / "missing") / "result.json"
            states[engine] = json.loads(result.read_text())["status"] if result.exists() else "missing"
        row = {"pair": f"{suite}/{pair}", "engine_status": states}
        if all(state == "complete" for state in states.values()):
            configs = [engines[e]["config"] for e in ("upstream", "plus")]
            if configs[0] != configs[1]:
                raise ValueError(f"Public comparison configs differ: {pair}")
            a = normalize_psms(read_table(paths["upstream"] / "results.sage.tsv"))
            b = normalize_psms(read_table(paths["plus"] / "results.sage.parquet"))
            row.update(status="complete", **compare_identifications(a, b))
        else:
            row["status"] = "incomplete"
        comparisons.append(row)
    return comparisons


def ptm_results(root):
    truth_path = root / "ptm-truth/truth.tsv"
    if not truth_path.exists():
        return []
    with truth_path.open() as f:
        truth = {r["sequence"]: {int(p) for p in r["positions"].split(",")}
                 for r in csv.DictReader(f, delimiter="\t") if r["unambiguous_site"] == "True"}
    results = []
    for suite in ("ptm", "ptm-oracle-v2"):
        for path in sorted((root / "runs" / suite).glob("*/results.sage.ptm-sites.parquet")):
            if not completed_output(path):
                continue
            rows = read_table(path)
            names = {r["filename"] for r in rows}
            mapping = {(name, seq): positions for name in names for seq, positions in truth.items()}
            result = localization_metrics(rows, mapping)
            result["joint_psm_peptide_localization_1pct"] = localization_metrics(
                [r for r in rows if float(r["peptide_q"]) <= 0.01], mapping)
            result["job"] = path.parent.name
            result["suite"] = suite
            result["scope"] = "Consistency with synthesis-defined sites across all libraries. File-to-library mapping remains unaudited."
            result["retained_site_rows"] = len(rows)
            result["peptide_q_one_fraction"] = sum(float(r["peptide_q"]) == 1 for r in rows) / len(rows) if rows else None
            log = path.parent / "search.log"
            result["heuristic_rescoring_fallback"] = "falling back to heuristic discriminant score" in log.read_text()
            results.append(result)
    return results


def quant_results(root):
    from scientific_metrics import boolean, stripped, canonical_peptide
    from generate_ptm_library_benchmark import fasta_entries
    results = []
    paths = [p for p in sorted((root / "runs/quantification").glob("*/lfq.parquet"))
             if completed_output(p)]
    if not paths:
        return results
    references = {name: "\x00".join(s.replace("I", "L") for _, s in
                  fasta_entries(root / f"references/{name}.fasta"))
                  for name in ("human", "yeast", "ecoli")}
    sequence_species = {}
    for path in paths:
        rows = read_table(path)
        strict_ms2 = accepted_ms2_keys(normalize_psms(read_table(path.parent / "results.sage.parquet")))
        species = json.loads((root / "references/hye-species.json").read_text())
        filenames = {r["filename"] for r in rows}
        design = {}
        for name in filenames:
            if "Condition_" in name:
                condition = name.split("Condition_")[1][0]
                preparation = name.split("Sample_")[1].split("_")[0]
                design[name] = {"condition": condition, "preparation": preparation}
        mixture_rows = [r for r in rows if r["filename"] in design]
        unambiguous = []
        cross_species = 0
        for row in mixture_rows:
            if boolean(row["is_decoy"]):
                continue
            sequence = stripped(row["peptide"]).replace("I", "L")
            if sequence not in sequence_species:
                sequence_species[sequence] = reference_membership(sequence, references)
            membership = sequence_species[sequence]
            if len(membership) > 1:
                cross_species += 1
                continue
            unambiguous.append(row)
        result = quantification_metrics(unambiguous, design, species, {"human": 0, "yeast": -1, "ecoli": 2})
        result["direct_ms2_ratio_diagnostic"] = direct_ratio_diagnostic(
            unambiguous, design, species, strict_ms2, {"human": 0, "yeast": -1, "ecoli": 2})
        result["lfq_threshold_yields"] = lfq_threshold_yields(rows)
        result["confidence_unit"] = "One precursor-peak q-value is repeated across all files. It is not an individual transfer q-value. The engine's discovery counter uses 5%, while the pilot's primary LFQ endpoint uses 1%."
        result["cross_species_il_rows_excluded"] = cross_species
        controls = [r for r in rows if "DDA_Human_" in r["filename"]]
        # Search engines can assign I/L-indistinguishable shared peptides differently.
        # Exclude any foreign assignment also occurring in the human reference.
        human = references["human"]
        candidates = {stripped(r["peptide"]).replace("I", "L") for r in controls
                      if {species.get(p) for p in r["proteins"].split(chr(59))} & {"yeast", "ecoli"}}
        shared = {s for s in candidates if s in human}
        absent = transferred = quantified = without_ms2 = strict_without_ms2 = strict_foreign = 0
        for row in controls:
            if boolean(row["is_decoy"]) or float(row["q_value"]) > 0.01 or not row.get("intensity") or float(row["intensity"]) <= 0:
                continue
            source_species = {species.get(p) for p in row["proteins"].split(chr(59))}
            if len(source_species) != 1 or None in source_species:
                continue
            if source_species & {"yeast", "ecoli"} and stripped(row["peptide"]).replace("I", "L") in shared:
                continue
            quantified += 1
            without_ms2 += not boolean(row["ms2_confirmed"])
            key = (Path(row["filename"]).name, canonical_peptide(row["peptide"]), str(row.get("charge") or ""))
            lacks_strict_ms2 = key not in strict_ms2
            strict_without_ms2 += lacks_strict_ms2
            if source_species & {"yeast", "ecoli"}:
                absent += 1
                transferred += not boolean(row["ms2_confirmed"])
                strict_foreign += lacks_strict_ms2
        result.update(job=path.parent.name, pure_human_control_quantified=quantified,
                      absent_species_quantified=absent, absent_species_without_ms2=transferred,
                      pure_human_control_without_ms2=without_ms2,
                      absent_species_fraction_among_unconfirmed=transferred / without_ms2 if without_ms2 else None,
                      control_without_jointly_accepted_ms2=strict_without_ms2,
                      foreign_without_jointly_accepted_ms2=strict_foreign,
                      foreign_fraction_without_jointly_accepted_ms2=strict_foreign / strict_without_ms2 if strict_without_ms2 else None,
                      absent_species_fraction=absent / quantified if quantified else None,
                      il_indistinguishable_foreign_sequences_excluded=len(shared))
        result["false_transfer_note"] = "The pure-human control measures foreign-species assignments among unconfirmed quantified features. It does not identify every incorrect human-to-human transfer or certify a total transfer FDR."
        results.append(result)
    return results


def accepted_ms2_keys(rows):
    """Direct evidence requires both PSM and peptide q-values at most 1%."""
    keys = set()
    for row in rows:
        if not row["is_decoy"] and row["spectrum_q"] <= 0.01 and row["peptide_q"] <= 0.01:
            filename, _, charge, peptide, _ = row["key"]
            keys.add((filename, peptide, str(charge)))
            keys.add((filename, peptide, ""))
    return keys


def reference_membership(sequence, references):
    """Conservatively identify any I/L-equivalent occurrence in each proteome."""
    normalized = sequence.replace("I", "L")
    return {species for species, proteins in references.items() if normalized in proteins}


def calibration_results(root):
    results = []
    for path in sorted((root / "runs").glob("entrapment-*/*/calibration.json")):
        data = json.loads(path.read_text())
        results.append({"suite": path.parent.parent.name, "job": path.parent.name,
                        "status": data["status"], "target_sequences": data["target_sequences"],
                        "thresholds": data["thresholds"]})
    return results


def calibration_uncertainty(results, replicates=10000, seed=20260914):
    """Resample files and shared entrapment seeds jointly for both engines."""
    summaries = []
    for study in ("human", "hye"):
        cells = {}
        for result in results:
            if not result["suite"].startswith(f"entrapment-{study}-"):
                continue
            construction = int(result["suite"].rsplit("-", 1)[1])
            _, file_id, engine = result["job"].split("-", 2)
            point = next(p for p in result["thresholds"] if p["nominal_q"] == 0.01)
            key = (int(file_id), construction, engine)
            if key in cells:
                raise ValueError("Duplicate calibration cell")
            cells[key] = point["paired_fdp_tie_max"]
        expected = {(file_id, construction, engine) for file_id in (0, 1)
                    for construction in (20260914, 20260915, 20260916)
                    for engine in ("upstream", "plus")}
        missing = sorted(expected - cells.keys())
        undefined = sum(v is None for v in cells.values())
        if missing or undefined:
            summaries.append({"study": study, "status": "insufficient_complete_evidence",
                              "missing_cells": missing, "undefined_fdp_cells": undefined})
            continue
        rng = random.Random(seed)
        samples = defaultdict(list)
        for _ in range(replicates):
            files = rng.choices((0, 1), k=2)
            seeds = rng.choices((20260914, 20260915, 20260916), k=3)
            values = {}
            for engine in ("upstream", "plus"):
                values[engine] = statistics.mean(cells[f, s, engine] for f in files for s in seeds)
                samples[engine].append(values[engine])
            samples["paired_difference"].append(values["plus"] - values["upstream"])
        intervals = {}
        for key, values in samples.items():
            ordered = sorted(values)
            intervals[key] = [ordered[int(0.025 * replicates)], ordered[int(0.975 * replicates)]]
        means = {engine: statistics.mean(cells[f, s, engine] for f in (0, 1)
                 for s in (20260914, 20260915, 20260916)) for engine in ("upstream", "plus")}
        summaries.append({"study": study, "status": "descriptive_pilot_interval",
                          "means": means, "percentile_intervals": intervals,
                          "bootstrap_replicates": replicates, "bootstrap_seed": seed,
                          "scope": "Conditional on two selected files and three shared construction seeds. Engines remain paired. This is not a between-study confidence interval or a power calculation."})
    return summaries


def main():
    global read_table
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--root", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    p.add_argument("--duckdb", type=Path, help="Use a retained DuckDB executable when replaying an archive")
    args = p.parse_args()
    if args.duckdb:
        read_table = partial(read_table, duckdb=args.duckdb.resolve())
    jobs = summarize_jobs(args.root)
    summary = {"status": "pilot_results", "jobs": jobs,
               "paired_identifications": compare_local(args.root),
               "public_identifications": compare_public(args.root),
               "exact_prefilter": check_prefilter(args.root),
               "ptm": ptm_results(args.root), "quantification": quant_results(args.root),
               "calibration": calibration_results(args.root)}
    summary["calibration_uncertainty"] = calibration_uncertainty(summary["calibration"])
    environment = args.root / "runs/public-timing/environment.json"
    if not environment.exists():
        environment = args.root / "runs/scaling/environment.json"
    if environment.exists():
        recorded = json.loads(environment.read_text())
        cpu_lines = recorded["cpu"].splitlines()
        summary["environment"] = {k: recorded[k] for k in
                                  ("platform", "python", "baseline_sha256", "candidate_sha256")}
        summary["environment"]["cpu_model"] = next((line.split(":", 1)[1].strip() for line in cpu_lines
                                                      if line.startswith("model name")), "unavailable")
        summary["environment"]["logical_processors_reported"] = sum(line.startswith("processor") for line in cpu_lines)
    summary["completed_jobs"] = sum(j["status"] == "complete" for j in jobs)
    summary["running_jobs"] = sum(j["status"] == "running" for j in jobs)
    summary["failed_or_invalid_jobs"] = sum(j["status"] not in ("complete", "running") for j in jobs)
    summary["limitations"] = [
        "Pilot selection and finite precision do not establish production calibration.",
        "Local HEK origin remains unverified and local timing may overlap acquisition or conversion work.",
        "Timing comes from one Linux workstation and does not establish performance across platforms or hardware.",
        "Seven E. coli proteins with undefined X residues are excluded from the repaired mixed-reference experiments.",
        "PTM site-library controls use synthesis truth and are oracle sensitivity experiments.",
        "PTM synthesis consistency includes identification and localization error and is not automatically an arrangement-level FLR estimate.",
        "Preparation variability, technical injections and entrapment seeds are not biological replication.",
        "Contaminant coverage and pure-species sample purity need validation before a final FDR or false-transfer claim."]
    inputs = list((args.root / "runs").glob("*/*/result.json"))
    inputs += list((args.root / "runs").glob("entrapment-*/*/calibration.json"))
    summary["input_manifests"] = file_identities(inputs)
    summary["analysis_sources"] = file_identities([Path(__file__).resolve(),
        Path(__file__).with_name("scientific_metrics.py"), Path(__file__).with_name("scientific_diagnostics.py")])
    atomic_json(args.output, summary)
    print(f"{summary['completed_jobs']} complete jobs, {summary['failed_or_invalid_jobs']} failures retained")


if __name__ == "__main__":
    main()
