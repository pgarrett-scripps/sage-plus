#!/bin/bash
# Run configs sequentially, record wall time.
O=/mnt/data1/scan-profile-bench/sim
BIN=/mnt/data1/build-cache/beta9/release/sage
export RAYON_NUM_THREADS=12
for name in "$@"; do
  s=$(date +%s)
  "$BIN" --overwrite "$O/cfg/$name.json" > "$O/out/$name.log" 2>&1
  echo "$name status=$? wall=$(( $(date +%s) - s ))s" >> "$O/timings.txt"
done
echo done >> "$O/timings.txt"
