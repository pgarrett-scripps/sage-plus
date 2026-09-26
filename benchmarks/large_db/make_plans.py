#!/usr/bin/env python3
"""Write the run plans for the large-database section.

Every plan starts from the paper's public search settings (tryptic, one missed
cleavage, 7 to 50 residues, fixed carbamidomethylation, 10 ppm precursor and
20 ppm fragment tolerances) with an explicit 16 GiB memory limit and eight
workers. Each database is searched with the prefilter and, where the memory
preflight allows it, without. `--min-matched` and `--max-peaks` set the
prefilter's match threshold and peak cap; the threshold sweep sets its own.
"""

from __future__ import annotations

import argparse
import copy
import json
from pathlib import Path

ROOT = Path("/data/sage-plus-scientific")
LARGE = ROOT / "large-db-20260925"
HEK = ROOT / "20260914/converted/PXD001468/b1906_293T_proteinID_01A_QE3_122212.mgf"
LFQ = ROOT / "20260914/converted/PXD028735/LFQ_Orbitrap_DDA_Condition_A_Sample_Alpha_01.mzML"
REFERENCES = ROOT / "20260914/references"
MULTIPLES = ["1x", "3x", "10x", "30x", "100x", "full"]
ENV = {"RAYON_NUM_THREADS": "8"}

BASE = {
    "batch_size": 1,
    "chimera": False,
    "database": {
        "bucket_size": 16384,
        "decoy_tag": "rev_",
        "enzyme": {"cleave_at": "KR", "max_len": 50, "min_len": 7,
                   "missed_cleavages": 1, "restrict": "P"},
        "generate_decoys": True,
        "static_mods": {"C": 57.021464},
        "prefilter_chunk_size": 0,
    },
    "deisotope": True,
    "fragment_tol": {"ppm": [-20, 20]},
    "isotope_errors": [-1, 3],
    "max_fragment_charge": 2,
    "max_memory_gb": 16,
    "min_free_memory_gb": 2,
    "min_matched_peaks": 6,
    "output_filter": {"psm_q_value": 1.0},
    "precursor_tol": {"ppm": [-10, 10]},
    "report_psms": 1,
    "score_type": "SageHyperScore",
}


PREFILTER: dict = {}


def config(fasta: Path, spectra: list[Path], prefilter: bool, oxidation: bool = False,
           isotope_errors: list[int] | None = None) -> dict:
    content = copy.deepcopy(BASE)
    if prefilter:
        content["database"] |= PREFILTER
    if isotope_errors is not None:
        content["isotope_errors"] = isotope_errors
    content["database"]["fasta"] = str(fasta)
    content["database"]["prefilter"] = prefilter
    if oxidation:
        content["database"]["variable_mods"] = {
            "M": [{"mass": 15.994915, "max_count": 1, "name": "Oxidation"}]
        }
        content["database"]["max_variable_mods"] = 1
    content["mzml_paths"] = [str(path) for path in spectra]
    return content


def pair(name: str, fasta: Path, spectra: list[Path], meta: dict,
         modes: tuple[str, ...] = ("prefilter", "full"), **kwargs) -> list[dict]:
    return [
        {"name": f"{name}-{mode}", "config": config(fasta, spectra, mode == "prefilter", **kwargs),
         "env": ENV, "gate_gb": 18, "meta": {**meta, "prefilter": mode == "prefilter"}}
        for mode in modes
    ]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--min-matched", type=int, help="prefilter_min_matched_peaks")
    parser.add_argument("--max-peaks", type=int, help="prefilter_max_peaks")
    args = parser.parse_args()
    if args.min_matched is not None:
        PREFILTER["prefilter_min_matched_peaks"] = args.min_matched
    if args.max_peaks is not None:
        PREFILTER["prefilter_max_peaks"] = args.max_peaks
    args.out.mkdir(parents=True, exist_ok=True)

    scaling = pair("hek-human", REFERENCES / "human.fasta", [HEK],
                   {"series": "scaling", "multiple": 0})
    for multiple in MULTIPLES:
        scaling += pair(f"hek-igc-{multiple}", LARGE / f"databases/human-igc-{multiple}.fasta",
                        [HEK], {"series": "scaling", "multiple": multiple})

    genomes = LARGE / "inputs/genomes"
    six_frame = pair("lfq-annotated", REFERENCES / "hye-irt-defined.fasta", [LFQ],
                     {"series": "six-frame", "database": "annotated"})
    six_frame += pair("lfq-six-frame", genomes / "human-irt-6frame.fasta", [LFQ],
                      {"series": "six-frame", "database": "six-frame"})

    # The precursor search space sets how many peptides the exact prefilter
    # keeps, so the same subsets are repeated with monoisotopic precursors only.
    narrow = []
    for multiple in ("3x", "10x", "30x"):
        narrow += pair(f"hek-igc-{multiple}-mono", LARGE / f"databases/human-igc-{multiple}.fasta",
                       [HEK], {"series": "narrow", "multiple": multiple},
                       modes=("prefilter",), isotope_errors=[0, 0])

    # Fecal metaproteomes against the sample-specific metagenome database.
    campi = LARGE / "inputs/campi"
    metaproteome = []
    for sample in ("F06", "F05"):
        metaproteome += pair(f"campi-{sample}", LARGE / "databases/gut-db2mg-human.fasta",
                             [campi / f"{sample}.mgf"],
                             {"series": "metaproteome", "sample": sample,
                              "database": "sample-specific"})

    # Prefilter threshold and peak cap on one subset, prefilter arm only.
    sweep = []
    for minimum in (1, 2, 3, 4, 6):
        for peaks in (None, 75, 50):
            PREFILTER.clear()
            PREFILTER["prefilter_min_matched_peaks"] = minimum
            if peaks is not None:
                PREFILTER["prefilter_max_peaks"] = peaks
            sweep += pair(f"hek-igc-10x-n{minimum}-p{peaks or 'all'}",
                          LARGE / "databases/human-igc-10x.fasta", [HEK],
                          {"series": "threshold-sweep", "multiple": "10x",
                           "min_matched": minimum, "max_peaks": peaks},
                          modes=("prefilter",))

    for name, plan in (("scaling", scaling), ("narrow", narrow), ("six-frame", six_frame),
                       ("metaproteome", metaproteome), ("threshold-sweep", sweep)):
        (args.out / f"{name}.json").write_text(json.dumps(plan, indent=1))


if __name__ == "__main__":
    main()
