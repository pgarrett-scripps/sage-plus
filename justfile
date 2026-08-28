set dotenv-load := true

baseline_ref := env_var_or_default("BASELINE_REF", "a38387b")
repeats := env_var_or_default("REPEATS", "3")
warmups := env_var_or_default("WARMUPS", "1")
threads := env_var_or_default("THREADS", "8")
memory_peptides := env_var_or_default("MEMORY_PEPTIDES", "1000000")
prefilter_chunk_size := env_var_or_default("PREFILTER_CHUNK_SIZE", "1000")

default:
    @just --list

# Build the candidate working tree and the pinned baseline.
bench-build:
    python3 benchmarks/benchmark.py build --baseline-ref "{{baseline_ref}}"

# Run the real search regression and prefilter benchmarks.
bench config:
    python3 benchmarks/benchmark.py all --config "{{config}}" --baseline-ref "{{baseline_ref}}" --repeats "{{repeats}}" --warmups "{{warmups}}" --threads "{{threads}}"

# Run the core suite with the local 219 MB HEK SILAC dataset.
bench-local:
    python3 benchmarks/benchmark.py all --config benchmarks/configs/local-standard.json --baseline-ref "{{baseline_ref}}" --repeats "{{repeats}}" --warmups "{{warmups}}" --threads "{{threads}}"

# Exercise SILAC, LFQ, and spectral-library export on the local real dataset.
bench-local-feature:
    python3 benchmarks/benchmark.py feature --config data/silac-k6r6/config.json --repeats "{{repeats}}" --warmups "{{warmups}}" --threads "{{threads}}"

# Compare a bounded variable-modification search on the local real dataset.
bench-local-mods:
    python3 benchmarks/benchmark.py search --config benchmarks/configs/local-modifications.json --baseline-ref "{{baseline_ref}}" --repeats "{{repeats}}" --warmups "{{warmups}}" --threads "{{threads}}"

# Measure scored deisotoping with inferred and supplied fragment charges.
bench-charge:
    python3 benchmarks/charge_benchmark.py --repeats "{{repeats}}" --warmups "{{warmups}}"

# Compare the baseline and candidate with a normal search configuration.
bench-search config:
    python3 benchmarks/benchmark.py search --config "{{config}}" --baseline-ref "{{baseline_ref}}" --repeats "{{repeats}}" --warmups "{{warmups}}" --threads "{{threads}}"

# Compare exact prefiltering off and on in the candidate working tree.
bench-prefilter config:
    python3 benchmarks/benchmark.py prefilter --config "{{config}}" --repeats "{{repeats}}" --warmups "{{warmups}}" --threads "{{threads}}" --prefilter-chunk-size "{{prefilter_chunk_size}}"

# Compare peak memory using a generated pre-digested peptide database.
bench-memory:
    python3 benchmarks/benchmark.py memory --baseline-ref "{{baseline_ref}}" --repeats "{{repeats}}" --warmups "{{warmups}}" --threads "{{threads}}" --memory-peptides "{{memory_peptides}}"

# Run a candidate-only configuration that enables the features being checked.
bench-feature config:
    python3 benchmarks/benchmark.py feature --config "{{config}}" --repeats "{{repeats}}" --warmups "{{warmups}}" --threads "{{threads}}"

# Check the benchmark harness without compiling Sage or running a search.
bench-check:
    python3 -m unittest discover -s benchmarks/tests -v
    python3 benchmarks/benchmark.py --help > /dev/null
    python3 benchmarks/charge_benchmark.py --help > /dev/null

# Remove generated builds and reports after an explicit confirmation flag.
bench-clean confirm="no":
    python3 benchmarks/benchmark.py clean --confirm "{{confirm}}"
