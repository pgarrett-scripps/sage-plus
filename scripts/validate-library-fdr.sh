#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 3 || $# -gt 4 ]]; then
  echo "usage: $0 RESULTS_PARQUET LIBRARY_PARQUET ENTRAPMENT_PROTEIN_REGEX [NON_ENTRAPMENT_TO_ENTRAPMENT_RATIO]" >&2
  exit 2
fi

results_path="$1"
library_path="$2"
entrapment_regex="$3"
ratio_override="${4:-}"

if [[ ! -f "$results_path" ]]; then
  echo "results file does not exist: $results_path" >&2
  exit 2
fi
if [[ ! -f "$library_path" ]]; then
  echo "library file does not exist: $library_path" >&2
  exit 2
fi
if ! command -v duckdb >/dev/null 2>&1; then
  echo "duckdb is required to inspect Sage Parquet output" >&2
  exit 2
fi
if [[ -n "$ratio_override" && ! "$ratio_override" =~ ^[0-9]+([.][0-9]+)?$ ]]; then
  echo "NON_ENTRAPMENT_TO_ENTRAPMENT_RATIO must be a non-negative number" >&2
  exit 2
fi

results_sql=${results_path//\'/\'\'}
library_sql=${library_path//\'/\'\'}
regex_sql=${entrapment_regex//\'/\'\'}
ratio_sql="${ratio_override:-NULL}"

duckdb -csv -c "
WITH thresholds(reported_q) AS (
    VALUES (0.001), (0.005), (0.010), (0.020), (0.050)
), psms AS (
    SELECT
        spectrum_q,
        is_decoy,
        NOT is_decoy AND regexp_matches(proteins, '${regex_sql}') AS is_entrapment
    FROM read_parquet('${results_sql}')
    WHERE rank = 1
), library_entries AS (
    SELECT
        library_entry_id,
        bool_or(regexp_matches(proteins, '${regex_sql}')) AS is_entrapment
    FROM read_parquet('${library_sql}')
    GROUP BY library_entry_id
), candidate_counts AS (
    SELECT
        count(*) FILTER (WHERE NOT is_entrapment) AS genuine_library_entries,
        count(*) FILTER (WHERE is_entrapment) AS entrapment_library_entries
    FROM library_entries
), candidate_ratio AS (
    SELECT
        genuine_library_entries,
        entrapment_library_entries,
        coalesce(
            ${ratio_sql},
            genuine_library_entries::DOUBLE / nullif(entrapment_library_entries, 0)
        ) AS non_entrapment_to_entrapment_ratio
    FROM candidate_counts
), counts AS (
    SELECT
        thresholds.reported_q,
        count(*) FILTER (WHERE NOT is_decoy AND spectrum_q <= thresholds.reported_q) AS reported_targets,
        count(*) FILTER (WHERE is_decoy AND spectrum_q <= thresholds.reported_q) AS reported_decoys,
        count(*) FILTER (WHERE is_entrapment AND spectrum_q <= thresholds.reported_q) AS entrapment_psms,
        count(*) FILTER (
            WHERE NOT is_decoy AND NOT is_entrapment AND spectrum_q <= thresholds.reported_q
        ) AS genuine_target_psms
    FROM thresholds
    CROSS JOIN psms
    GROUP BY thresholds.reported_q
)
SELECT
    reported_q,
    reported_targets,
    reported_decoys,
    entrapment_psms,
    genuine_target_psms,
    genuine_library_entries,
    entrapment_library_entries,
    non_entrapment_to_entrapment_ratio,
    CASE WHEN genuine_target_psms = 0 THEN NULL ELSE
        entrapment_psms * non_entrapment_to_entrapment_ratio / genuine_target_psms
    END AS scaled_entrapment_fdr
FROM counts
CROSS JOIN candidate_ratio
ORDER BY reported_q;
"
