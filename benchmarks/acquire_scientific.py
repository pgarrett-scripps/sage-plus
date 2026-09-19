#!/usr/bin/env python3
"""Acquire a prespecified public pilot, with byte limits and content receipts."""

import argparse
import hashlib
import json
import shutil
import urllib.request
from pathlib import Path

from provenance import atomic_json, sha256


REPO = Path(__file__).resolve().parents[1]
METADATA = REPO / "benchmarks/results/dissertation-20260914/metadata"
SELECTION = {
    "PXD001468": [
        "b1906_293T_proteinID_01A_QE3_122212.mzXML",
        "b1937_293T_proteinID_01B_QE3_122212.mzXML",
    ],
    "PXD028735": [
        f"LFQ_Orbitrap_DDA_Condition_{condition}_Sample_{sample}_01.raw"
        for sample in ("Alpha", "Beta") for condition in ("A", "B")
    ],
    "PXD000138": [
        "HCD_1.raw.-1.mgf", "HCD_2.raw.-1.mgf",
        "009606_IPI_v3_72_with_P_Lib_062111.fasta",
    ],
}


def selected_files(metadata):
    selected = []
    for accession, names in SELECTION.items():
        inventory = metadata / f"{accession}-files-all.json"
        if not inventory.exists():
            inventory = metadata / f"{accession}-files.json"
        entries = {row["fileName"]: row for row in json.loads(inventory.read_text())}
        project = json.loads((metadata / f"{accession}-project.json").read_text())
        for name in names:
            row = entries[name]
            ftp = next(p["value"] for p in row["publicFileLocations"]
                       if p["name"] == "FTP Protocol")
            selected.append({
                "study": accession, "name": name,
                "url": ftp.replace("ftp://", "https://", 1),
                "bytes": int(row["fileSizeBytes"]), "license": project["license"],
                "provider_checksum": row.get("checksum") or None,
                "role": "pilot",
            })
    return selected


def ensure_metadata(destination):
    """Keep public API responses with the experiment, including all inventory pages."""
    destination.mkdir(parents=True, exist_ok=True)
    for accession in SELECTION:
        project = destination / f"{accession}-project.json"
        files = destination / f"{accession}-files-all.json"
        if not project.exists():
            cached = METADATA / project.name
            if cached.exists():
                shutil.copyfile(cached, project)
            else:
                url = f"https://www.ebi.ac.uk/pride/ws/archive/v2/projects/{accession}"
                with urllib.request.urlopen(url, timeout=60) as response:
                    atomic_json(project, json.load(response))
        if not files.exists():
            cached = METADATA / files.name
            if accession == "PXD001468" and not cached.exists():
                cached = METADATA / f"{accession}-files.json"
            if cached.exists():
                shutil.copyfile(cached, files)
            else:
                rows = []
                for page in range(100):
                    url = f"https://www.ebi.ac.uk/pride/ws/archive/v2/projects/{accession}/files?page={page}&pageSize=100"
                    with urllib.request.urlopen(url, timeout=60) as response:
                        current = json.load(response)
                    rows.extend(current)
                    if len(current) < 100:
                        break
                else:
                    raise RuntimeError("Dataset inventory exceeded pagination bound")
                if len({r["fileName"] for r in rows}) != len(rows):
                    raise ValueError("Duplicate file names in paginated inventory")
                atomic_json(files, rows)
    return destination


def acquire(item, root, reserve_bytes):
    path = root / "inputs" / item["study"] / item["name"]
    receipt = path.with_name(path.name + ".receipt.json")
    if receipt.exists() and path.exists():
        previous = json.loads(receipt.read_text())
        if previous["input"] == item and previous["sha256"] == sha256(path):
            print("Verified cached input", item["name"], flush=True)
            return previous
        raise RuntimeError(f"Existing input or receipt changed: {path}")
    if path.exists():
        raise RuntimeError(f"Input without receipt already exists: {path}")
    if shutil.disk_usage(root).free < item["bytes"] + reserve_bytes:
        raise RuntimeError("Insufficient free space including the configured reserve")
    path.parent.mkdir(parents=True, exist_ok=True)
    partial = path.with_name(path.name + ".partial")
    h = hashlib.sha256()
    size = 0
    request = urllib.request.Request(item["url"], headers={"User-Agent": "SagePlus-scientific-benchmark/1"})
    with urllib.request.urlopen(request, timeout=120) as response, partial.open("wb") as out:
        headers = dict(response.headers)
        final_url = response.url
        while block := response.read(1024 * 1024):
            size += len(block)
            if size > item["bytes"]:
                raise RuntimeError(f"Download exceeds declared size: {path}")
            out.write(block)
            h.update(block)
    if size != item["bytes"]:
        raise RuntimeError(f"Incomplete input: {path}, {size} bytes")
    partial.rename(path)
    record = {"input": item, "sha256": h.hexdigest(), "bytes": size,
              "final_url": final_url, "headers": headers,
              "checksum_scope": "locally computed, provider checksum may be unavailable"}
    atomic_json(receipt, record)
    print("Acquired", item["name"], size, h.hexdigest(), flush=True)
    return record


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--metadata", type=Path)
    parser.add_argument("--download-gib", type=float, default=20)
    parser.add_argument("--reserve-gib", type=float, default=20)
    parser.add_argument("--study", choices=list(SELECTION))
    parser.add_argument("--plan-only", action="store_true")
    args = parser.parse_args()
    root = args.root.resolve()
    root.mkdir(parents=True, exist_ok=True)
    metadata = args.metadata or ensure_metadata(root / "metadata")
    items = selected_files(metadata)
    total = sum(item["bytes"] for item in items)
    if total > args.download_gib * 1024 ** 3:
        raise RuntimeError("Prespecified selection exceeds total download budget")
    plan = {"status": "pilot_selection_before_outcomes", "files": items,
            "total_bytes": total, "selection_basis":
            "First named HEK fraction A/B, first injections of Alpha/Beta mixture preparations, first two HCD synthetic libraries. Remaining files reserved for later evaluation.",
            "split_note": "A/B in HEK names are fractions or acquisition series, not independent studies. Alpha/Beta are mixture preparations, not biological replicates.",
            "download_budget_gib": args.download_gib}
    plan_path = root / "acquisition-plan.json"
    if plan_path.exists() and json.loads(plan_path.read_text()) != plan:
        raise RuntimeError("Refusing to replace a different acquisition plan")
    atomic_json(plan_path, plan)
    print(f"Prespecified {len(items)} inputs, {total / 1024**3:.2f} GiB", flush=True)
    if args.plan_only:
        return
    for item in items:
        if args.study is None or item["study"] == args.study:
            acquire(item, root, int(args.reserve_gib * 1024 ** 3))


if __name__ == "__main__":
    main()
