"""Exploratory endpoints that remain separate from the primary confidence filter."""

import math
import statistics
from collections import defaultdict
from pathlib import Path

from scientific_metrics import boolean, canonical_peptide


def lfq_threshold_yields(rows, thresholds=(0.01, 0.05)):
    result = {}
    for threshold in thresholds:
        precursors = set()
        positive = 0
        for row in rows:
            if boolean(row["is_decoy"]) or float(row["q_value"]) > threshold:
                continue
            precursors.add((canonical_peptide(row["peptide"]), str(row.get("charge") or "")))
            intensity = float(row["intensity"]) if row.get("intensity") else 0
            positive += math.isfinite(intensity) and intensity > 0
        result[str(threshold)] = {"target_precursors": len(precursors), "positive_target_file_rows": positive}
    return result


def direct_ratio_diagnostic(rows, design, species_by_protein, accepted_ms2, expected):
    values = defaultdict(dict)
    for row in rows:
        if boolean(row["is_decoy"]):
            continue
        filename = Path(row["filename"]).name
        peptide = canonical_peptide(row["peptide"])
        charge = str(row.get("charge") or "")
        if filename not in design or (filename, peptide, charge) not in accepted_ms2:
            continue
        species = {species_by_protein.get(p) for p in row["proteins"].split(chr(59))}
        if len(species) != 1 or not species <= expected.keys():
            continue
        intensity = float(row["intensity"]) if row.get("intensity") else 0
        if not math.isfinite(intensity) or intensity <= 0:
            continue
        key = (species.pop(), peptide, charge)
        if filename in values[key]:
            raise ValueError("Duplicate diagnostic feature and file")
        values[key][filename] = intensity
    errors = defaultdict(list)
    for (species, _, _), observed in values.items():
        for preparation in sorted({info["preparation"] for info in design.values()}):
            pair = {}
            for filename, info in design.items():
                if info["preparation"] == preparation and filename in observed:
                    if info["condition"] in pair:
                        raise ValueError("Diagnostic requires one injection per preparation and condition")
                    pair[info["condition"]] = observed[filename]
            if set(pair) == {"A", "B"}:
                errors[species].append(math.log2(pair["B"]) - math.log2(pair["A"]) - expected[species])
    return {"status": "posthoc_diagnostic",
            "selection": "Positive intensities with direct target MS2 evidence passing both 1% PSM and peptide thresholds. The LFQ q-value filter is deliberately not applied. Cross-species I/L ambiguity is excluded by the caller.",
            "scope": "Explains extraction behavior before LFQ confidence filtering. These are not accepted 1% LFQ discoveries or a calibration claim.",
            "species": {species: {"ratio_pairs": len(errors[species]),
                        "median_log2_bias": statistics.median(errors[species]) if errors[species] else None,
                        "median_absolute_log2_error": statistics.median(map(abs, errors[species])) if errors[species] else None}
                        for species in expected}}
