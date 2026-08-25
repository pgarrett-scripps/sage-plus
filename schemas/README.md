# Sage analytical schemas

These files are the versioned, machine-readable Parquet message schemas for Sage's canonical analytical outputs.

- `results.sage.v1.parquet.schema` describes `results.sage.parquet`.
- `results.sage.v2.parquet.schema` adds precursor label channel and group identity.
- `lfq.v1.parquet.schema` describes the separate long-form `lfq.parquet` table.
- `lfq.v2.parquet.schema` adds label identity and reference-channel ratios.
- `spectral_library.sage.v1.parquet.schema` describes the empirical, long-form
  `spectral_library.sage.parquet` transition table.
- `spectral_library.sage.v2.parquet.schema` preserves label channel, group, and reference metadata.
- `scores.v1.md` defines the score and evidence fields used by those schemas.

Within a schema major version, fields may be added only when existing readers can safely ignore them. Removing a field, changing its physical type or nullability, changing row granularity, or changing a score's meaning requires a new schema major version. Files embed `sage.schema.name` and `sage.schema.version` in their Parquet key-value metadata.

`results.sage.parquet` and `matched_fragments.sage.parquet` also embed the inclusive output cutoff as `sage.output_filter.spectrum_q_max`. Every fragment row belongs to a PSM retained in `results.sage.parquet` under that cutoff.

The spectral-library table contains one row per selected fragment transition. Its Parquet metadata
records the schema version, selection strategy, PSM and peptide q-value cutoffs, minimum consensus
support, and minimum fragment frequency. `library_entry_id` groups transitions belonging to the
same exact peptidoform and precursor charge.

Unlabeled searches continue to write version 1 schemas. A configured precursor-label search writes
version 2 results and LFQ schemas. Spectral libraries use version 2 only when labeled entries are
present.
