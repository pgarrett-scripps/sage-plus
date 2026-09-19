#!/usr/bin/env python3
"""Recheck retained input and output hashes before freezing pilot evidence."""

import argparse
import json
from pathlib import Path

from provenance import atomic_json, sha256


def audit(root):
    manifests = sorted((root / "runs").glob("*/*/result.json"))
    manifests += sorted((root / "runs").glob("*/*/calibration.json"))
    manifests += sorted((root / "runs").glob("*/generation.json"))
    manifests += sorted((root / "converted").glob("*/*.conversion.json"))
    manifests += [root / "ptm-truth/manifest.json"]
    manifests += sorted((root / "references").glob("*reference.json"))
    actual = {}
    errors = []
    checks = 0

    def check(path, expected, manifest):
        nonlocal checks
        checks += 1
        path = Path(path)
        if path not in actual:
            actual[path] = sha256(path) if path.is_file() else None
        if actual[path] != expected:
            errors.append({"manifest": str(manifest), "file": str(path),
                           "expected": expected, "actual": actual[path]})

    for manifest in manifests:
        if not manifest.exists():
            errors.append({"manifest": str(manifest), "reason": "missing"})
            continue
        data = json.loads(manifest.read_text())
        if data.get("status") == "running":
            errors.append({"manifest": str(manifest), "reason": "still running"})
        for identities in (data.get("signature", {}).get("inputs", {}),
                           data.get("inputs", {}), data.get("outputs", {})):
            for path, expected in identities.items():
                check(path, expected, manifest)
    for receipt in sorted((root / "inputs").glob("*/*.receipt.json")):
        data = json.loads(receipt.read_text())
        path = receipt.with_name(receipt.name.removesuffix(".receipt.json"))
        check(path, data["sha256"], receipt)
        if not path.exists() or path.stat().st_size != data["bytes"]:
            errors.append({"manifest": str(receipt), "reason": "input size mismatch"})
    for matrix in sorted((root / "runs").glob("*/matrix-status.json")):
        if json.loads(matrix.read_text())["pending"]:
            errors.append({"manifest": str(matrix), "reason": "pending jobs"})
    partials = list(root.rglob("*.partial"))
    errors.extend({"file": str(p), "reason": "incomplete acquisition"} for p in partials)
    return {"status": "verified" if not errors else "failed", "hash_checks": checks,
            "unique_files_hashed": len(actual), "search_and_conversion_manifests": len(manifests),
            "errors": errors,
            "scope": "Evidence integrity and recorded matrix completion. Preserved failed experiments remain failures. This does not certify scientific conclusions."}


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--root", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    args = p.parse_args()
    result = audit(args.root.resolve())
    atomic_json(args.output, result)
    print(result["status"], result["hash_checks"], "checks", flush=True)
    if result["errors"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
