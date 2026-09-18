"""Evaluate complete peptide q-value steps at a common entrapment FDP ceiling."""
import csv
import gzip
import hashlib
import itertools
import json
import math
import sys
from pathlib import Path
from statistics import mean

PAPER = Path(__file__).resolve().parents[2]
REPO = PAPER.parent
sys.path.insert(0, str(REPO / "benchmarks"))
from scientific_entrapment import threshold_fdp

ROOT = Path("/data/sage-plus-scientific/20260914/runs")
OUT = PAPER / "analysis/data/matched-fdp.json"
CEILING = 0.01


def digest(path):
    with path.open("rb") as handle:
        return hashlib.file_digest(handle, "sha256").hexdigest()


def curve(rows, labels, partners):
    """Include complete q-value ties and conservatively count exact score ties."""
    if len({r["sequence"] for r in rows}) != len(rows):
        raise ValueError("Repeated peptide sequence")
    if any(r["sequence"] not in labels for r in rows):
        raise ValueError("Unlabeled peptide")
    if any(not math.isfinite(r["score"]) or not 0 <= r["peptide_q"] <= 1 for r in rows):
        raise ValueError("Invalid score or q-value")
    reverse = {target: entrapment for entrapment, target in partners.items()}
    selected = {}
    numerator = targets = entrapments = 0

    def contribution(sequence):
        row = selected.get(sequence)
        if row is None:
            return 0
        partner = selected.get(partners[sequence])
        if partner is None:
            return 2
        return 1 + 2 * ((row["peptide_q"], -row["score"]) <=
                        (partner["peptide_q"], -partner["score"]))

    result = []
    for q, group in itertools.groupby(sorted(rows, key=lambda r: r["peptide_q"]),
                                      key=lambda r: r["peptide_q"]):
        group = list(group)
        affected = set()
        for row in group:
            seq = row["sequence"]
            if labels[seq] == "p_target":
                affected.add(seq)
            elif seq in reverse:
                affected.add(reverse[seq])
        numerator -= sum(contribution(seq) for seq in affected)
        for row in group:
            selected[row["sequence"]] = row
            targets += labels[row["sequence"]] == "target"
            entrapments += labels[row["sequence"]] == "p_target"
        numerator += sum(contribution(seq) for seq in affected)
        result.append(dict(nominal_q=q, targets=targets, entrapments=entrapments,
                           denominator=len(selected), numerator=numerator,
                           paired_fdp=numerator / len(selected)))
    return result


def select(points):
    eligible = [p for p in points if p["targets"] > 0 and p["paired_fdp"] <= CEILING]
    return max(eligible, key=lambda p: (p["targets"], p["nominal_q"])) if eligible else None


def main():
    inputs = {}
    runs = []
    for study in ("human", "hye"):
        for seed in (20260914, 20260915, 20260916):
            directory = ROOT / f"entrapment-{study}-{seed}"
            extracted = {}
            audits = {}
            for file in (0, 1):
                for engine in ("upstream", "plus"):
                    run = directory / f"file-{file}-{engine}"
                    source = run / "fdrbench-input.tsv"
                    audit_path = run / "calibration.json"
                    audit = json.loads(audit_path.read_text())
                    inputs[str(source)] = digest(source)
                    inputs[str(audit_path)] = digest(audit_path)
                    assert audit["status"] == "independently_checked"
                    assert inputs[str(source)] == audit["outputs"][str(source)]
                    with source.open() as handle:
                        extracted[file, engine] = [dict(sequence=r["peptide"],
                            peptide_q=float(r["q_value"]), score=float(r["score"]))
                            for r in csv.DictReader(handle, delimiter="\t")]
                    audits[file, engine] = audit
            observed = {r["sequence"] for rows in extracted.values() for r in rows}
            pairing = directory / "paired.txt"
            inputs[str(pairing)] = digest(pairing)
            assert all(a["inputs"][str(pairing)] == inputs[str(pairing)] for a in audits.values())
            labels, pair_ids = {}, {}
            with pairing.open() as handle:
                for row in csv.DictReader(handle, delimiter="\t"):
                    seq, kind = row["sequence"], row["peptide_type"]
                    if seq in observed and kind in ("target", "p_target"):
                        assert seq not in labels
                        labels[seq] = kind
                        pair_ids.setdefault(row["peptide_pair_index"], {})[kind] = seq
            partners = {members["p_target"]: members.get("target")
                        for members in pair_ids.values() if "p_target" in members}
            for (file, engine), rows in extracted.items():
                points = curve(rows, labels, partners)
                chosen = select(points)
                assert chosen is not None, "No nonempty set below the FDP ceiling"
                for original in audits[file, engine]["thresholds"]:
                    q = original["nominal_q"]
                    independent = threshold_fdp(rows, labels, partners, q)
                    assert independent["targets"] == original["targets"]
                    assert independent["entrapments"] == original["entrapments"]
                    assert independent["paired_fdp_tie_max"] == original["paired_fdp_tie_max"]
                    eligible = [p for p in points if p["nominal_q"] <= q]
                    if eligible:
                        assert eligible[-1]["paired_fdp"] == independent["paired_fdp_tie_max"]
                independent = threshold_fdp(rows, labels, partners, chosen["nominal_q"])
                assert chosen["paired_fdp"] == independent["paired_fdp_tie_max"]
                runs.append(dict(study=study, seed=seed, file=file, engine=engine,
                                 selected=chosen, curve=points))
            print(f"Audited {study} seed {seed}", flush=True)
    summaries = []
    for study in ("human", "hye"):
        group = [r for r in runs if r["study"] == study]
        engines = {}
        for engine in ("upstream", "plus"):
            points = [r["selected"] for r in group if r["engine"] == engine]
            engines[engine] = dict(mean_targets=mean(p["targets"] for p in points),
                mean_fdp=mean(p["paired_fdp"] for p in points),
                fdp_min=min(p["paired_fdp"] for p in points),
                fdp_max=max(p["paired_fdp"] for p in points),
                q_min=min(p["nominal_q"] for p in points),
                q_max=max(p["nominal_q"] for p in points))
        paired = []
        for seed in (20260914, 20260915, 20260916):
            for file in (0, 1):
                counts = {r["engine"]: r["selected"]["targets"] for r in group
                          if r["seed"] == seed and r["file"] == file}
                paired.append(100 * (counts["plus"] / counts["upstream"] - 1))
        summaries.append(dict(study=study, engines=engines,
            mean_paired_percent_change=mean(paired),
            min_paired_percent_change=min(paired), max_paired_percent_change=max(paired)))
    sources = [Path(__file__).resolve(), REPO / "benchmarks/scientific_entrapment.py"]
    curves_path = PAPER / "analysis/results/matched-fdp-curves.json.gz"
    curves_path.parent.mkdir(parents=True, exist_ok=True)
    curves = [dict(study=r["study"], seed=r["seed"], file=r["file"],
                   engine=r["engine"], curve=r.pop("curve")) for r in runs]
    curves_path.write_bytes(gzip.compress(json.dumps(curves, separators=(",", ":")).encode(), mtime=0))
    output = dict(ceiling=CEILING, selection="Largest target count among complete observed q-value steps at or below the conservative paired FDP ceiling. Ties choose the largest q-value.",
        scope="Post hoc threshold selection on the same entrapment observations used to estimate FDP. Descriptive yield at a common estimated ceiling, not held-out error control or a sensitivity test.",
        inputs_sha256=inputs, analysis_sources_sha256={str(p): digest(p) for p in sources},
        curve_file=str(curves_path.relative_to(PAPER)), curve_sha256=digest(curves_path),
        summaries=summaries, runs=runs)
    OUT.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n")
    print(json.dumps(summaries, indent=2))


if __name__ == "__main__":
    main()
