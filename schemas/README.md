# Sage analytical schemas

These files are the versioned, machine-readable Parquet message schemas for Sage's canonical analytical outputs.

- `results.sage.v1.parquet.schema` describes `results.sage.parquet`.
- `lfq.v1.parquet.schema` describes the separate long-form `lfq.parquet` table.
- `scores.v1.md` defines the score and evidence fields used by those schemas.

Within a schema major version, fields may be added only when existing readers can safely ignore them. Removing a field, changing its physical type or nullability, changing row granularity, or changing a score's meaning requires a new schema major version. Files embed `sage.schema.name` and `sage.schema.version` in their Parquet key-value metadata.

`results.sage.parquet` and `matched_fragments.sage.parquet` also embed the inclusive output cutoff as `sage.output_filter.spectrum_q_max`. Every fragment row belongs to a PSM retained in `results.sage.parquet` under that cutoff.
