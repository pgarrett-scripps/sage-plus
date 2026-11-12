#!/usr/bin/env bash
# Public head-to-head: Sage Plus reading vendor .raw against upstream Sage
# reading mzML, with the conversion the mzML arm needs counted in its time.
#
#   ./run_public_benchmark.sh --arm plus-raw     --out DIR file1.raw [file2.raw ...]
#   ./run_public_benchmark.sh --arm upstream-mzml --out DIR file1.raw [file2.raw ...]
#
# Both arms take the same vendor files and the same configuration; the
# upstream arm converts them first. Every path, binary and version used is
# written to DIR/receipt.json. See README.md for how to obtain the inputs.
set -uo pipefail

arm=""
out=""
config="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/pxd028735.json"
fasta="${SAGE_BENCH_FASTA:-}"
plus_bin="${SAGE_BENCH_PLUS:-sage}"
upstream_bin="${SAGE_BENCH_UPSTREAM:-sage-upstream}"
converter="${SAGE_BENCH_CONVERTER:-docker}"
threads="${SAGE_BENCH_THREADS:-${AJS_CPU:-}}"
mzml_cache="${SAGE_BENCH_MZML_CACHE:-}"

while [ $# -gt 0 ]; do
  case "$1" in
    --arm) arm=$2; shift 2 ;;
    --out) out=$2; shift 2 ;;
    --config) config=$2; shift 2 ;;
    --fasta) fasta=$2; shift 2 ;;
    --mzml-cache) mzml_cache=$2; shift 2 ;;
    --) shift; break ;;
    -*) echo "unknown option: $1" >&2; exit 2 ;;
    *) break ;;
  esac
done
files=("$@")

[ -n "$arm" ] && [ -n "$out" ] && [ -n "$fasta" ] && [ ${#files[@]} -gt 0 ] || {
  echo "usage: $0 --arm {plus-raw|upstream-mzml} --out DIR --fasta FASTA file.raw ..." >&2
  exit 2
}
case "$arm" in plus-raw|upstream-mzml) ;; *) echo "unknown arm: $arm" >&2; exit 2 ;; esac

mkdir -p "$out"
out=$(cd "$out" && pwd)
search_dir=$out/search
rm -rf "$search_dir"; mkdir -p "$search_dir"

# ---------------------------------------------------------------- conversion
inputs=()
convert_seconds=0
convert_source=none
if [ "$arm" = upstream-mzml ]; then
  # The conversion is part of this arm's cost, so it is timed once and the
  # time is recorded beside the mzML. A rerun pointed at the same cache
  # reuses both the files and that recorded time instead of paying it again.
  mzml_dir=${mzml_cache:-$out/mzml}
  mkdir -p "$mzml_dir"
  record=$mzml_dir/conversion.json
  missing=0
  for raw in "${files[@]}"; do
    [ -s "$mzml_dir/$(basename "${raw%.*}").mzML" ] || missing=1
  done
  if [ "$missing" = 0 ] && [ -s "$record" ]; then
    convert_seconds=$(python3 -c 'import json,sys;print(json.load(open(sys.argv[1]))["convert_seconds"])' "$record")
    convert_source=recorded
    for raw in "${files[@]}"; do
      inputs+=("$mzml_dir/$(basename "$(readlink -f "${raw%.*}")").mzML")
    done
    echo "reusing $mzml_dir (conversion took ${convert_seconds}s when it was measured)"
  else
  convert_start=$(date +%s.%N)
  for raw in "${files[@]}"; do
    raw=$(readlink -f "$raw")
    case "$converter" in
      docker)
        # Digest-pinned so the converter is the same build for everyone.
        docker run --rm --network none --user "$(id -u):$(id -g)" \
          --volume "$(dirname "$raw")":"$(dirname "$raw")":ro \
          --volume "$mzml_dir":"$mzml_dir":rw \
          quay.io/biocontainers/thermorawfileparser@sha256:465b623c81a8f300d2802196246307e4280e9135ebfcc4c4526a2beac7e7ea8f \
          ThermoRawFileParser --input="$raw" --output_directory="$mzml_dir" --format=1 \
          > "$out/convert.$(basename "$raw").log" 2>&1 || exit 1
        ;;
      *)
        # Any local ThermoRawFileParser with the same flags.
        "$converter" --input="$raw" --output_directory="$mzml_dir" --format=1 \
          > "$out/convert.$(basename "$raw").log" 2>&1 || exit 1
        ;;
    esac
    inputs+=("$mzml_dir/$(basename "${raw%.*}").mzML")
  done
  convert_seconds=$(echo "$(date +%s.%N) - $convert_start" | bc)
  convert_source=measured
  printf '{"convert_seconds": %.2f, "converter": "%s", "files": %s}\n' \
    "$convert_seconds" "$converter" "$(printf '%s\n' "${inputs[@]}" | python3 -c 'import json,sys;print(json.dumps([l.strip() for l in sys.stdin]))')" \
    > "$record"
  fi
else
  for raw in "${files[@]}"; do inputs+=("$(readlink -f "$raw")"); done
fi

# ------------------------------------------------------------------ the search
run_config=$out/config.used.json
python3 - "$config" "$fasta" "$run_config" "${inputs[@]}" <<'PY'
import json, sys
config, fasta, dest, *inputs = sys.argv[1:]
spec = json.load(open(config))
spec["database"]["fasta"] = fasta
spec["mzml_paths"] = inputs
json.dump(spec, open(dest, "w"), indent=2)
PY

if [ "$arm" = plus-raw ]; then
  binary=$plus_bin
  extra=(--overwrite)
else
  binary=$upstream_bin
  # Upstream Sage asks to be told not to phone home, has no --overwrite (the
  # script always writes into a fresh directory), and writes TSV unless asked
  # for parquet, which is what both arms are counted from.
  extra=(--disable-telemetry-i-dont-want-to-improve-sage --parquet)
fi
# Both engines take their thread count from rayon, so this is the one knob
# that gives them the same number of cores.
[ -n "$threads" ] && export RAYON_NUM_THREADS="$threads"

/usr/bin/time -v "$binary" "$run_config" -o "$search_dir" "${extra[@]}" \
  > "$out/search.stdout" 2> "$out/search.time"
search_status=$?
if [ $search_status -ne 0 ]; then
  # Nothing is timed or counted for a search that did not run: say what
  # happened and stop, so a broken command cannot look like a slow engine.
  echo "$arm: the search failed with exit $search_status; its first error was:" >&2
  grep -m3 -v '^\s' "$out/search.time" >&2
  printf '{"arm": "%s", "exit_code": %d, "binary": "%s", "searched": false}\n' \
    "$arm" "$search_status" "$binary" > "$out/receipt.json"
  exit $search_status
fi
search_seconds=$(awk -F': ' '/Elapsed \(wall clock\)/ {split($2,t,":");
  print (length(t)==3 ? t[1]*3600+t[2]*60+t[3] : t[1]*60+t[2])}' "$out/search.time")
peak_rss_kb=$(awk -F': ' '/Maximum resident set size/ {print $2}' "$out/search.time")

# ----------------------------------------------------------------- the numbers
python3 - "$arm" "$out" "$search_dir" "$convert_seconds" "${search_seconds:-0}" \
         "${peak_rss_kb:-0}" "$search_status" "$binary" "$run_config" "$convert_source" <<'PY'
import json, subprocess, sys, pathlib
(arm, out, search_dir, convert_s, search_s, rss_kb, status, binary, config,
 convert_source) = sys.argv[1:]
counts = {}
path = pathlib.Path(search_dir, "results.sage.parquet")
if path.exists():
    import pyarrow.parquet as pq
    table = pq.read_table(path, columns=None)
    names = table.column_names
    # Sage Plus writes `is_decoy`; upstream Sage writes `label` (1 / -1).
    if "is_decoy" in names:
        target = [not d for d in table["is_decoy"].to_pylist()]
    else:
        target = [v == 1 for v in table["label"].to_pylist()]
    def at_one_percent(column, key=None):
        if column not in names:
            return None
        q = table[column].to_pylist()
        rows = [i for i in range(len(q)) if target[i] and q[i] is not None and q[i] <= 0.01]
        if key is None:
            return len(rows)
        values = table[key].to_pylist()
        return len({values[i] for i in rows})
    counts = {
        "psms_at_1pct_fdr": at_one_percent("spectrum_q"),
        "peptides_at_1pct_fdr": at_one_percent("peptide_q", "peptide"),
        "proteins_at_1pct_fdr": at_one_percent("protein_q", "proteins"),
    }
receipt = {
    "arm": arm,
    "exit_code": int(status),
    "binary": binary,
    "binary_version": subprocess.run([binary, "--version"], capture_output=True,
                                     text=True).stdout.strip(),
    "config": config,
    "searched": True,
    "convert_seconds": round(float(convert_s), 2),
    "convert_time_source": convert_source,
    "search_seconds": round(float(search_s), 2),
    "total_seconds": round(float(convert_s) + float(search_s), 2),
    "peak_rss_gib": round(int(rss_kb) / 1048576, 3),
    **counts,
}
json.dump(receipt, open(pathlib.Path(out, "receipt.json"), "w"), indent=2)
print(json.dumps(receipt, indent=2))
PY
exit $search_status
