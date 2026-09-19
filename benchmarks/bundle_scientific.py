#!/usr/bin/env python3
"""Archive scientific evidence and verify every archived member by SHA-256."""

import argparse
import hashlib
import json
import platform
import subprocess
import tarfile
import tempfile
from pathlib import Path

from provenance import atomic_json, sha256
from run_scientific import BASELINE, CANDIDATE
from scientific_entrapment import JAR
from scientific_metrics import DUCKDB


REPO = Path(__file__).resolve().parents[1]
LARGE_INPUT_SUFFIXES = {".raw", ".mzxml", ".mzml", ".mgf"}


def select_files(root):
    included, external = {}, {}
    for path in sorted(root.rglob("*")):
        if not path.is_file() or path.is_symlink() or "artifacts" in path.relative_to(root).parts:
            continue
        relative = path.relative_to(root)
        if path.name.endswith(".partial"):
            raise ValueError(f"Acquisition is incomplete: {relative}")
        if "__pycache__" in relative.parts:
            continue
        if path.suffix.lower() in LARGE_INPUT_SUFFIXES:
            external[str(relative)] = {"sha256": sha256(path), "bytes": path.stat().st_size}
        else:
            included["evidence/" + str(relative)] = path
    for path in sorted((REPO / "benchmarks").glob("*.py")):
        included["code/benchmarks/" + path.name] = path
    for path in sorted((REPO / "benchmarks/tests").glob("test_*.py")):
        included["code/benchmarks/tests/" + path.name] = path
    for name in ("SCIENTIFIC_PROTOCOL.md", "SCIENTIFIC_PILOT.md", "datasets.json"):
        path = REPO / "benchmarks" / name
        if path.exists():
            included["code/benchmarks/" + name] = path
    for path in sorted((REPO / "benchmarks/results/dissertation-20260914").rglob("*")):
        if path.is_file():
            # This preliminary page fetch is an unrelated drug-pipeline article.
            if path.relative_to(REPO / "benchmarks/results/dissertation-20260914").as_posix() == "metadata/ptm-paper.txt":
                continue
            included["evidence/repository-metadata/" + str(path.relative_to(REPO / "benchmarks/results/dissertation-20260914"))] = path
    included["binaries/sage-upstream"] = BASELINE
    included["binaries/sage-plus-beta3"] = CANDIDATE
    included["tools/fdrbench-1.1.1.jar"] = JAR
    included["tools/duckdb"] = DUCKDB
    for result in sorted((root / "runs").glob("*/*/result.json")):
        identities = json.loads(result.read_text())["signature"]["inputs"]
        for name, digest in identities.items():
            path = Path(name)
            if path.is_relative_to(root):
                continue
            if path.suffix.lower() in LARGE_INPUT_SUFFIXES:
                external[str(path)] = {"sha256": digest, "bytes": path.stat().st_size,
                                       "scope": "Existing local engineering input. Acquisition identity and redistribution terms remain unverified."}
            elif path.suffix.lower() in {".fasta", ".fa", ".tsv"}:
                included[f"evidence/engineering-inputs/{digest[:16]}-{path.name}"] = path
    return included, external


def verify(archive, expected):
    seen = set()
    with tarfile.open(archive, "r:gz") as bundle:
        for member in bundle:
            if not member.isfile() or member.name not in expected or member.name in seen:
                raise ValueError(f"Unexpected archive member: {member.name}")
            seen.add(member.name)
            stream = bundle.extractfile(member)
            h = hashlib.sha256()
            for block in iter(lambda: stream.read(1024 * 1024), b""):
                h.update(block)
            if h.hexdigest() != expected[member.name]:
                raise ValueError(f"Archive checksum mismatch: {member.name}")
    if seen != set(expected):
        raise ValueError("Archive is incomplete")


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--root", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    args = p.parse_args()
    root, output = args.root.resolve(), args.output.resolve()
    if output.exists():
        raise ValueError("Use a new archive name to preserve earlier evidence snapshots")
    output.parent.mkdir(parents=True, exist_ok=True)
    included, external = select_files(root)
    expected = {name: sha256(path) for name, path in included.items()}
    with tempfile.TemporaryDirectory(dir=output.parent) as temporary:
        temporary = Path(temporary)
        for label, revision in (("sage-plus-beta3", "v0.1.0-beta.3"), ("sage-upstream", "df9219951cc9a54c")):
            source = temporary / f"{label}-source.tar.gz"
            with source.open("wb") as stream:
                subprocess.run(["git", "archive", "--format=tar.gz", revision], cwd=REPO, stdout=stream, check=True)
            name = "sources/" + source.name
            included[name], expected[name] = source, sha256(source)
        manifest = temporary / "manifest.json"
        atomic_json(manifest, {"schema_version": 1, "status": "local_evidence_archive",
                    "archive_members": expected, "external_spectra": external,
                    "original_data_root": str(root), "python": platform.python_version(),
                    "scope": "Public raw and converted spectra remain on /data with acquisition receipts and hashes. Existing local engineering spectra remain at their original paths with unverified acquisition identity. The archive contains results, references, executable binaries, source, configs and analysis scripts. Publication and external deposition have not occurred."})
        included["manifest.json"], expected["manifest.json"] = manifest, sha256(manifest)
        with tarfile.open(output, "w:gz", compresslevel=3) as archive:
            for name, path in sorted(included.items()):
                info = archive.gettarinfo(str(path), arcname=name)
                info.uid, info.gid, info.uname, info.gname = 0, 0, "", ""
                with path.open("rb") as stream:
                    archive.addfile(info, stream)
        verify(output, expected)
    atomic_json(output.with_suffix(".verification.json"), {
        "status": "verified", "archive_sha256": sha256(output),
        "members_verified": len(expected), "external_spectra": len(external),
        "archive_bytes": output.stat().st_size})
    print(output, flush=True)


if __name__ == "__main__":
    main()
