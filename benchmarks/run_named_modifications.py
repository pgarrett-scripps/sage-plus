#!/usr/bin/env python3
"""Exercise named modifications through discovery, library reuse, and iteration."""
import argparse
import csv
import json
from pathlib import Path
import re
import subprocess
import tempfile


def run(binary, work):
    root = Path(__file__).resolve().parents[1]
    source = (root / "crates/sage/src/mass.rs").read_text()
    table = re.search(r"MONOISOTOPIC_MASSES: \[f32; 26\] = \[(.*?)\];", source, re.S).group(1)
    masses = [float(value) for value in table.split(",") if value.strip()]
    sequence = "KSTAGVQK"
    work.mkdir(parents=True, exist_ok=True)
    fasta = work / "proteins.fasta"
    fasta.write_text(f">P1\n{sequence}\n")
    cases = [
        ("first", ["first_residue:K"], [0], 1, "residue", 1),
        ("nterm", ["peptide_n_term:K"], ["N"], 1, "peptide_n_term", 1),
        ("last", ["last_residue:K"], [7], 1, "residue", 8),
        ("cterm", ["peptide_c_term:K"], ["C"], 1, "peptide_c_term", 8),
        ("both", ["peptide_n_term", "peptide_c_term"], ["N", "C"], 2, None, None),
        ("ambiguous", ["first_residue:K", "peptide_n_term:K"], ["N"], 0, None, None),
    ]
    evidence = []
    for mode in ("database", "mass_offset"):
        for name, sites, attachments, expected_count, expected_attachment, expected_position in cases:
            if name == "both" and mode == "mass_offset":
                continue
            case = work / f"{mode}-{name}"
            case.mkdir(exist_ok=True)
            residue_masses = [masses[ord(aa) - 65] + (42.010565 if i in attachments else 0) for i, aa in enumerate(sequence)]
            nterm = 42.010565 if "N" in attachments else 0
            cterm = 42.010565 if "C" in attachments else 0
            total = sum(residue_masses) + 18.010565 + nterm + cterm
            peaks = []
            for cut in range(1, len(sequence)):
                peaks.extend([sum(residue_masses[:cut]) + nterm + 1.0072764,
                              sum(residue_masses[cut:]) + cterm + 18.010565 + 1.0072764])
            mgf = case / "spectra.mgf"
            mgf.write_text("\n".join(["BEGIN IONS", f"TITLE={name}", "RTINSECONDS=60", f"PEPMASS={total / 2 + 1.0072764}", "CHARGE=2+"] + [f"{mz:.6f} 1000" for mz in sorted(peaks)] + ["END IONS", ""]))
            config = {
                "database": {
                    "fasta": str(fasta), "generate_decoys": True,
                    "enzyme": {"cleave_at": "$", "min_len": 5, "max_len": 50},
                    "variable_mods": {"TestMod": {"mass": 42.010565, "sites": sites,
                        "max_count": len(attachments), "search_mode": mode}},
                    "max_variable_mods": len(attachments), "max_total_variable_mods": len(attachments),
                },
                "mzml_paths": [str(mgf)], "precursor_tol": {"ppm": [-10, 10]},
                "fragment_tol": {"ppm": [-10, 10]}, "deisotope": False,
                "min_peaks": 5, "min_matched_peaks": 4, "predict_rt": False,
                "output_filter": {"psm_q_value": 1.0},
                "ptm_localization": {"enabled": True, "psm_q_value": 1.0, "localization_q_value": 1.0},
            }
            previous = None
            phases = ["discovery"] if name == "ambiguous" else ["discovery", "guided", "iteration"]
            for phase in phases:
                if previous:
                    config["database"]["ptm_library"] = {"path": str(previous), "strict": True}
                    config["database"]["variable_mods"]["TestMod"]["site_mode"] = "library" if phase == "guided" else "both"
                path = case / f"{phase}.json"
                path.write_text(json.dumps(config, indent=2) + "\n")
                output = case / phase
                with (case / f"{phase}.log").open("w") as log:
                    subprocess.run([str(binary), str(path), "-o", str(output), "--write-pin", "--overwrite", "--disable-telemetry-i-dont-want-to-improve-sage"], stdout=log, stderr=log, check=True)
                library = output / "results.sage.ptm-library.tsv"
                with library.open() as handle:
                    rows = list(csv.DictReader(handle, delimiter="\t"))
                assert len(rows) == expected_count, (mode, name, phase, rows)
                if expected_attachment:
                    assert rows[0]["attachment"] == expected_attachment, rows
                    assert int(rows[0]["position"]) == expected_position, rows
                if name == "both":
                    assert {row["attachment"] for row in rows} == {"peptide_n_term", "peptide_c_term"}, rows
                assert all(row["modification"] == "TestMod" for row in rows)
                if previous:
                    assert library.read_text() == previous.read_text(), (mode, name, phase)
                previous = library
                evidence.append({"mode": mode, "case": name, "phase": phase, "library_sites": len(rows), "passed": True})
    return evidence


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sage", type=Path, required=True)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--work", type=Path)
    args = parser.parse_args()
    binary = args.sage.resolve()
    if args.work:
        results = run(binary, args.work.resolve())
    else:
        with tempfile.TemporaryDirectory(prefix="sage-named-mods-") as directory:
            results = run(binary, Path(directory))
    text = json.dumps({"checks": len(results), "results": results}, indent=2) + "\n"
    if args.output:
        args.output.write_text(text)
    print(text)
