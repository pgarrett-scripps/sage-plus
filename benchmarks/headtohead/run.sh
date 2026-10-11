#!/usr/bin/env bash
# Head-to-head: upstream Sage v0.15.0-beta.2 vs Sage Plus Beta 16 RC on PXD028735 HYE (4 files).
# Usage: ./run.sh <upstream|plus> <run number> [memgate GB]
# Engines must be run sequentially, never in parallel.
set -euo pipefail
# Runs go to $HEADTOHEAD_DIR (default: this directory). config.json must point at your local
# PXD028735 mzMLs and FASTA.
script_dir="$(cd "$(dirname "$0")" && pwd)"
here="${HEADTOHEAD_DIR:-$script_dir}"
engine="$1"; n="$2"; gb="${3:-12}"
UPSTREAM="${SAGE_UPSTREAM_BINARY:?set SAGE_UPSTREAM_BINARY to the upstream v0.15.0-beta.2 sage}"
PLUS="${SAGE_PLUS_BINARY:?set SAGE_PLUS_BINARY to the Sage Plus sage}"
out="$here/$engine/run$n"
mkdir -p "$out"
case "$engine" in
  upstream) cmd=("$UPSTREAM" --disable-telemetry-i-dont-want-to-improve-sage -o "$out" "$script_dir/config.json") ;;
  plus)     cmd=("$PLUS" --overwrite -o "$out" "$script_dir/config.json") ;;
  *) echo "engine must be upstream or plus" >&2; exit 2 ;;
esac
cp "$script_dir/config.json" "$out/config.used.json"
printf '%q ' "${cmd[@]}" > "$out/command.txt"; echo >> "$out/command.txt"
date -Is > "$out/started_at.txt"
${SAGE_MEMGATE:+$SAGE_MEMGATE "$gb"} /usr/bin/time -v -o "$out/time.txt" "${cmd[@]}" > "$out/stdout.log" 2> "$out/stderr.log"
date -Is > "$out/finished_at.txt"
