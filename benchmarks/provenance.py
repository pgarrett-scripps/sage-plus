"""Content-verified stage reuse for scientific benchmark jobs."""

from __future__ import annotations

import hashlib
import json
import os
import tempfile
from contextlib import contextmanager
from pathlib import Path


def sha256(path: Path) -> str:
    with path.open("rb") as handle:
        digest = hashlib.sha256()
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
        return digest.hexdigest()


def atomic_json(path: Path, value: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = None
    try:
        with tempfile.NamedTemporaryFile(mode="w", dir=path.parent, delete=False) as handle:
            temporary = Path(handle.name)
            json.dump(value, handle, indent=2, sort_keys=True, allow_nan=False)
            handle.write("\n")
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, path)
    finally:
        if temporary is not None:
            temporary.unlink(missing_ok=True)


def file_identities(paths: list[Path]) -> dict[str, str]:
    return {str(path.resolve(strict=True)): sha256(path) for path in paths}


@contextmanager
def cached_stage(manifest: Path, stage: str, inputs: list[Path], outputs: list[Path], parameters: object):
    """Reuse only completed matching stages with intact outputs.

    Callers run stages sequentially within a private experiment directory.
    A miss removes only the explicitly declared generated outputs.
    """
    signature = {
        "stage": stage,
        "inputs": file_identities(inputs),
        "parameters": parameters,
    }
    expected = {str(path.resolve()) for path in outputs}
    try:
        previous = json.loads(manifest.read_text())
        hit = (
            previous.get("schema_version") == 1
            and previous.get("signature") == signature
            and set(previous.get("outputs", {})) == expected
            and previous["outputs"] == file_identities(outputs)
        )
    except (OSError, ValueError, TypeError, AttributeError):
        hit = False
    if hit:
        yield True
        return
    manifest.unlink(missing_ok=True)
    for path in outputs:
        path.unlink(missing_ok=True)
    yield False
    if file_identities(inputs) != signature["inputs"]:
        raise RuntimeError(f"stage {stage} inputs changed during execution")
    for path in outputs:
        if path.stat().st_size == 0:
            raise RuntimeError(f"stage {stage} produced an empty file: {path}")
    atomic_json(manifest, {
        "schema_version": 1,
        "signature": signature,
        "outputs": file_identities(outputs),
    })
