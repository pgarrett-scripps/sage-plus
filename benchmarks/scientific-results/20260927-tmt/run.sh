#!/usr/bin/env bash
# One-command rerun: run.sh <sage-plus-binary> <label> [extra tmt_hypro16.py args]
# e.g. run.sh /mnt/data1/build-cache/b13-tmt-target/release/sage b13
set -euo pipefail
here="$(cd "$(dirname "$0")" && pwd)"
export UV_CACHE_DIR="${UV_CACHE_DIR:-/mnt/data1/build-cache/uv-cache}"
exec uv run --no-project --with pyarrow==25.0.1 --with numpy==2.5.3 \
    python "$here/tmt_hypro16.py" --plus-bin "$1" --plus-label "$2" "${@:3}"
