#!/usr/bin/env python3
"""Summarize the large-database runs into one committed JSON file.

Reads `large-db-runs.json` from each run series and the `results.sage.parquet`
of each completed run. Accepted identifications are target PSMs with spectrum
q-value at most 0.01, and target peptides with peptide q-value at most 0.01.

Scaling series: peptides are split by the proteins that contain them. A
peptide in any human reference protein counts as human; a peptide found only in
catalog proteins is an entrapment identification, because the spectra are from
a human cell line. The combined entrapment estimate is
N_e (1 + 1/r) / (N_t + N_e), with r the catalog-to-human ratio of residues.

Six-frame series: microbial peptides from the translated genomes are checked
against the annotated E. coli and yeast proteomes, with I and L equated.

Narrow series: the same subsets searched with monoisotopic precursors only.
Metaproteome series: fecal runs against the sample-specific
metagenome database with the human reference.

Requires the `duckdb` Python package.
"""

from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

import duckdb

HUMAN_PREFIXES = ("sp|", "human:", "irt:")


def load_runs(path: Path) -> dict:
    return json.loads((path / "large-db-runs.json").read_text())["runs"]


def results(directory: Path) -> str:
    """A relation with stripped_peptide, proteins, is_decoy and both q-values.

    Upstream Sage writes only the tab-separated report, whose peptide column
    carries modification masses in brackets."""
    parquet = directory / "results.sage.parquet"
    if parquet.exists():
        path = str(parquet).replace("'", "''")
        return f"read_parquet('{path}')"
    path = str(directory / "results.sage.tsv").replace("'", "''")
    return ("(select regexp_replace(peptide, '\\[[^]]*\\]', '', 'g') as stripped_peptide, "
            "proteins, label = -1 as is_decoy, spectrum_q, peptide_q "
            f"from read_csv('{path}', delim='\t', header=true))")


def accepted(directory: Path) -> tuple[int, list[tuple[str, str]]]:
    """Accepted PSM count and (stripped peptide, proteins) for accepted peptides."""
    relation = results(directory)
    psms = duckdb.sql(
        f"select count(*) from {relation} "
        "where not is_decoy and spectrum_q <= 0.01"
    ).fetchone()[0]
    peptides = duckdb.sql(
        f"select stripped_peptide, any_value(proteins) from {relation} "
        "where not is_decoy and peptide_q <= 0.01 group by stripped_peptide"
    ).fetchall()
    return psms, peptides


def is_human(proteins: str) -> bool:
    return any(protein.startswith(HUMAN_PREFIXES) for protein in proteins.split(";"))


def outcome(record: dict, directory: Path) -> str:
    """completed, refused (a memory preflight declined to start), guard (the
    runtime memory guard stopped the search), or ceiling (an allocation failed
    at the address-space ceiling)."""
    if record["exit_code"] == 0:
        return "completed"
    log = (directory / "stderr.log").read_text(errors="replace")
    if "would exceed `max_memory_gb`" in log or "currently available while preserving" in log:
        return "refused"
    if "reached its configured memory limit" in log:
        return "guard"
    if "memory allocation of" in log or "capacity overflow" in log:
        return "ceiling"
    return "failed"


def refused_estimate(record: dict) -> float | None:
    """Total GiB a preflight refusal says the search needed: the estimated
    additional peak plus Sage's resident memory at the time."""
    for line in record.get("error_lines", []):
        if match := re.search(r"peak is ([\d.]+) GiB in addition to Sage's current ([\d.]+) GiB",
                              line):
            return float(match[1]) + float(match[2])
    return None


def resources(record: dict) -> dict:
    keep = ("exit_code", "wall_seconds", "peak_rss_gib", "database_peptides",
            "psms_at_one_percent_fdr", "peptides_at_one_percent_fdr", "results_sha256",
            "preflight", "final_preflight", "spectrum_index", "error_lines")
    row = {key: record[key] for key in keep if key in record}
    # The log's prefilter counts; the row's `prefilter` is the search mode.
    if "prefilter" in record:
        row["prefilter_counts"] = record["prefilter"]
    return row


def scaling(root: Path, receipt: dict,
            reference: set[str] | None = None) -> tuple[list[dict], set[str] | None]:
    """Rows for one series and the human-only reference peptide set, which is
    taken from this series when it contains the human-only prefilter run."""
    if not (root / "large-db-runs.json").exists():
        return [], reference
    runs = load_runs(root)
    rows = []
    for name, record in runs.items():
        meta = record["meta"]
        row = {"name": name, "multiple": meta["multiple"], "prefilter": meta["prefilter"],
               "outcome": outcome(record, root / name),
               "refused_estimate_gib": refused_estimate(record), **resources(record)}
        if meta["multiple"]:
            subset = receipt["subsets"]["1000" if meta["multiple"] == "full" else meta["multiple"].rstrip("x")]
            row["catalog_proteins"] = subset["proteins"]
            row["catalog_residues"] = subset["residues"]
            row["residue_ratio"] = subset["residues"] / receipt["human_residues"]
        else:
            row["catalog_proteins"] = row["catalog_residues"] = 0
            row["residue_ratio"] = 0.0
        if record["exit_code"] == 0:
            psms, peptides = accepted(root / name)
            human = {peptide for peptide, proteins in peptides if is_human(proteins)}
            entrapment = len(peptides) - len(human)
            row |= {"accepted_psms": psms, "accepted_peptides": len(peptides),
                    "human_peptides": len(human), "entrapment_peptides": entrapment}
            if row["residue_ratio"]:
                r = row["residue_ratio"]
                row["combined_fdp"] = 100 * entrapment * (1 + 1 / r) / len(peptides)
            row["_human"] = human
            if meta["multiple"] == 0 and meta["prefilter"]:
                reference = human
        rows.append(row)
    for row in rows:
        human = row.pop("_human", None)
        if human is not None and reference:
            row["reference_peptides_retained"] = len(human & reference)
            row["reference_peptides"] = len(reference)
    return rows, reference


def normalize(sequence: str) -> str:
    return sequence.replace("I", "L")


def proteome_text(paths: list[Path]) -> str:
    parts = []
    for path in paths:
        sequence = []
        for line in path.read_text().splitlines():
            if line.startswith(">"):
                parts.append("".join(sequence))
                sequence = []
            else:
                sequence.append(line.strip())
        parts.append("".join(sequence))
    return normalize("#".join(parts))


def six_frame(root: Path, annotated_microbes: list[Path]) -> list[dict]:
    runs = load_runs(root)
    microbes = proteome_text(annotated_microbes)
    rows = []
    microbial_sets = {}
    for name, record in runs.items():
        meta = record["meta"]
        row = {"name": name, "database": meta["database"], "prefilter": meta["prefilter"],
               "outcome": outcome(record, root / name), **resources(record)}
        if record["exit_code"] == 0:
            psms, peptides = accepted(root / name)
            microbial = {peptide for peptide, proteins in peptides if not is_human(proteins)}
            annotated = {peptide for peptide in microbial if normalize(peptide) in microbes}
            row |= {"accepted_psms": psms, "accepted_peptides": len(peptides),
                    "microbial_peptides": len(microbial),
                    "microbial_annotated": len(annotated),
                    "microbial_unannotated": len(microbial) - len(annotated)}
            if meta["prefilter"]:
                microbial_sets[meta["database"]] = {normalize(p) for p in microbial}
        rows.append(row)
    if {"annotated", "six-frame"} <= microbial_sets.keys():
        shared = len(microbial_sets["annotated"] & microbial_sets["six-frame"])
        for row in rows:
            row["microbial_shared_with_other_database"] = shared
    return rows


def metaproteome(root: Path) -> list[dict]:
    if not (root / "large-db-runs.json").exists():
        return []
    rows = []
    for name, record in load_runs(root).items():
        meta = record["meta"]
        row = {"name": name, "sample": meta["sample"], "database": meta["database"],
               "prefilter": meta["prefilter"], "outcome": outcome(record, root / name),
               "refused_estimate_gib": refused_estimate(record), **resources(record)}
        if record["exit_code"] == 0:
            psms, peptides = accepted(root / name)
            human = sum(is_human(proteins) for _, proteins in peptides)
            row |= {"accepted_psms": psms, "accepted_peptides": len(peptides),
                    "human_peptides": human, "microbial_peptides": len(peptides) - human}
        rows.append(row)
    return rows


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    runs = args.root / "runs"
    receipt = json.loads((args.root / "databases/scaling-databases.json").read_text())
    references = Path("/data/sage-plus-scientific/20260914/references")
    scaling_rows, reference = scaling(runs / "scaling", receipt)
    summary = {
        "executable": json.loads((runs / "scaling/large-db-runs.json").read_text())["executable"],
        "databases": receipt,
        "campi_databases": json.loads((args.root / "databases/campi-databases.json").read_text()),
        "scaling": scaling_rows,
        "narrow": scaling(runs / "narrow", receipt, reference)[0],
        "six_frame": six_frame(runs / "six-frame",
                               [references / "ecoli.fasta", references / "yeast.fasta"]),
        "metaproteome": metaproteome(runs / "metaproteome"),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(summary, indent=1, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()
