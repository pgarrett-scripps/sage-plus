"""Common scientific endpoints independent of engine summary counters."""

import csv
import io
import math
import re
import statistics
import subprocess
from collections import defaultdict
from decimal import Decimal
from pathlib import Path


THRESHOLDS = (0.001, 0.005, 0.01, 0.02, 0.05)
DUCKDB = Path("/snap/duckdb/9/duckdb")


def canonical_peptide(value):
    """Preserve modification positions and masses, normalizing numeric spelling."""
    return re.sub(r"\[([+-]?[0-9.]+)\]",
                  lambda m: "[" + str(Decimal(m[1]).normalize()) + "]", value)


def stripped(value):
    return re.sub(r"\[[^]]*\]", "", value).replace("-", "")


def read_table(path, duckdb=DUCKDB):
    path = Path(path)
    if path.suffix != ".parquet":
        with path.open(newline="") as f:
            return list(csv.DictReader(f, delimiter="\t"))
    escaped = str(path.resolve()).replace("'", "''")
    command = [str(duckdb), "-csv", "-c", f"SELECT * FROM read_parquet('{escaped}')"]
    result = subprocess.run(command, check=True, capture_output=True, text=True)
    return list(csv.DictReader(io.StringIO(result.stdout)))


def boolean(value):
    if value in (True, "true", "True", "1", 1):
        return True
    if value in (False, "false", "False", "0", 0):
        return False
    raise ValueError(f"Invalid boolean: {value!r}")


def probability(value):
    value = float(value)
    if not math.isfinite(value) or not 0 <= value <= 1:
        raise ValueError(f"Invalid probability: {value}")
    return value


def normalize_psms(rows):
    normalized = []
    seen = set()
    for row in rows:
        if int(row["rank"]) != 1:
            continue
        peptide = canonical_peptide(row["peptide"])
        decoy = boolean(row["is_decoy"]) if "is_decoy" in row else int(row["label"]) == -1
        if "label" in row and int(row["label"]) not in (-1, 1):
            raise ValueError("Unknown target/decoy label")
        key = (Path(row["filename"]).name, row["scannr"], int(row["charge"]), peptide, decoy)
        if key in seen:
            raise ValueError(f"Duplicate rank-one PSM identity: {key}")
        seen.add(key)
        score = float(row["sage_discriminant_score"])
        if not math.isfinite(score):
            raise ValueError("Nonfinite ranking score")
        normalized.append({"key": key, "peptide": peptide,
                           "sequence": stripped(peptide), "is_decoy": decoy,
                           "spectrum_q": probability(row["spectrum_q"]),
                           "peptide_q": probability(row["peptide_q"]),
                           "score": score, "hyperscore": float(row["hyperscore"]),
                           "proteins": row["proteins"]})
    return normalized


def peptide_representatives(rows):
    """Select one actual PSM by the engine's discriminant ranking, never mix rows."""
    groups = defaultdict(list)
    for row in rows:
        groups[(row["peptide"], row["is_decoy"])].append(row)
    result = []
    for group in groups.values():
        qs = {r["peptide_q"] for r in group}
        if max(qs) - min(qs) > 1e-6:
            raise ValueError("Inconsistent q-values within a peptidoform")
        result.append(min(group, key=lambda r: (-r["score"], r["key"])))
    return result


def identification_metrics(rows):
    peptides = peptide_representatives(rows)
    return {str(q): {
        "target_psms": sum(not r["is_decoy"] and r["spectrum_q"] <= q for r in rows),
        "decoy_psms": sum(r["is_decoy"] and r["spectrum_q"] <= q for r in rows),
        "target_peptidoforms": sum(not r["is_decoy"] and r["peptide_q"] <= q for r in peptides),
        "target_sequences": len({r["sequence"] for r in peptides if not r["is_decoy"] and r["peptide_q"] <= q}),
    } for q in THRESHOLDS}


def compare_identifications(left, right, q=0.01):
    def identities(rows):
        return {r["key"] for r in rows if not r["is_decoy"] and r["spectrum_q"] <= q}
    a, b = identities(left), identities(right)
    maps = [{r["key"]: r for r in rows} for rows in (left, right)]
    common = maps[0].keys() & maps[1].keys()
    return {"q": q, "shared_target_psms": len(a & b), "baseline_only": len(a - b),
            "candidate_only": len(b - a),
            "jaccard": len(a & b) / len(a | b) if a | b else None,
            "all_psm_identities_equal": maps[0].keys() == maps[1].keys(),
            "maximum_score_difference": max((abs(maps[0][k]["score"] - maps[1][k]["score"]) for k in common), default=None),
            "maximum_peptide_q_difference": max((abs(maps[0][k]["peptide_q"] - maps[1][k]["peptide_q"]) for k in common), default=None)}


def quantification_metrics(rows, design, species_by_protein, expected_log2_b_over_a):
    """No imputation. Restrict ratios to unambiguous species and matched preparations."""
    values = defaultdict(dict)
    observed = defaultdict(set)
    ambiguous = 0
    transfer_count = 0
    for row in rows:
        if boolean(row["is_decoy"]) or probability(row["q_value"]) > 0.01:
            continue
        filename = Path(row["filename"]).name
        if filename not in design:
            raise ValueError(f"Missing sample design for {filename}")
        species = {species_by_protein.get(p) for p in row["proteins"].split(";")}
        if len(species) != 1 or None in species:
            ambiguous += 1
            continue
        species = species.pop()
        if species not in expected_log2_b_over_a:
            continue
        intensity = float(row["intensity"]) if row.get("intensity") else None
        if intensity is None or not math.isfinite(intensity) or intensity <= 0:
            continue
        key = (species, canonical_peptide(row["peptide"]), row.get("charge", ""))
        if filename in values[key]:
            raise ValueError("Duplicate LFQ feature and file")
        values[key][filename] = intensity
        observed[species].add(key)
        transfer_count += not boolean(row["ms2_confirmed"])
    result = {}
    preparations = sorted({v["preparation"] for v in design.values()})
    for species, expected in expected_log2_b_over_a.items():
        errors, ratios, cvs = [], [], []
        slots = len(observed[species]) * len(design)
        present = 0
        for key in observed[species]:
            v = values[key]
            present += len(v)
            for preparation in preparations:
                pair = {}
                for filename, info in design.items():
                    if info["preparation"] == preparation and filename in v:
                        if info["condition"] in pair:
                            raise ValueError("Multiple injections per preparation require explicit aggregation")
                        pair[info["condition"]] = v[filename]
                if set(pair) == {"A", "B"}:
                    ratio = math.log2(pair["B"] / pair["A"])
                    ratios.append(ratio)
                    errors.append(ratio - expected)
            for condition in ("A", "B"):
                replicates = [v[f] for f, info in design.items() if info["condition"] == condition and f in v]
                if len(replicates) >= 2:
                    cvs.append(statistics.stdev(replicates) / statistics.mean(replicates))
        result[species] = {"expected_log2_b_over_a": expected, "ratio_pairs": len(ratios),
                           "median_log2_ratio": statistics.median(ratios) if ratios else None,
                           "median_log2_bias": statistics.median(errors) if errors else None,
                           "median_absolute_log2_error": statistics.median(map(abs, errors)) if errors else None,
                           "median_preparation_cv": statistics.median(cvs) if cvs else None,
                           "missing_fraction_observed_union": 1 - present / slots if slots else None,
                           "observable_features": len(observed[species])}
    return {"species": result, "ambiguous_or_unmapped_rows": ambiguous,
            "quantified_rows_without_ms2_confirmation": transfer_count,
            "false_transfer_rate": None,
            "false_transfer_note": "Requires species-absent controls. Missing MS2 confirmation alone is not a false transfer."}


def localization_metrics(rows, truth, q=0.01):
    """Truth maps (file, unmodified peptide) to allowed phosphorylated positions."""
    correct = incorrect = unassessable = 0
    seen = set()
    for row in rows:
        if probability(row["spectrum_q"]) > q or probability(row["localization_q_value"]) > q:
            continue
        if abs(float(row["modification_mass"]) - 79.966331) > 0.001:
            continue
        key = (Path(row["filename"]).name, stripped(row["peptide"]))
        identity = (key, row["scannr"], int(row["position"]))
        if identity in seen:
            raise ValueError("Duplicate localization event")
        seen.add(identity)
        if key not in truth:
            unassessable += 1
        elif int(row["position"]) in truth[key]:
            correct += 1
        else:
            incorrect += 1
    n = correct + incorrect
    return {"correct_site_events": correct, "incorrect_site_events": incorrect,
            "unassessable_site_events": unassessable,
            "empirical_site_error_fraction": incorrect / n if n else None,
            "unit": "site events, distinct from arrangement-level FLR"}
