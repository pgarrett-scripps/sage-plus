# Grounding rules: tests and benchmarks

This maps each applicable hard constraint (HC) of the proteomics grounding spec to the tests or
benchmarks that check it, and lists the gaps. It reflects the state after Beta 12 and the Beta 13
provenance work. Differential statistics (HC-STAT-01 to 03) and SRM/PRM (HC-QUANT-04) do not
apply to a search engine.

Test paths are relative to the repository. Unit tests in `crates/sage/tests/unit/` are compiled
into the `sage-core` crate. `integration.rs` means `crates/sage-cli/tests/integration.rs`.

| Rule | Covered by | Gaps |
| --- | --- | --- |
| HC-FDR-01 separate FDR per level | `fdr.rs`: `picked_peptide_assigns_one_to_orphaned_competition_twins`, `picked_protein_counts_unique_proteins_and_skips_shared_peptides`, `decoy_protein_groups_compete_with_their_target_group`; `ml/qvalue.rs`: `groups_are_estimated_separately` | Protein-site reports carry PSM-level q only; no site-level FDR test (no site FDR exists). |
| HC-FDR-02 decoy strategy | `database.rs`: `filtered_decoy_generation_drops_reversals_that_collide_with_targets`, `generated_decoys_that_are_twins_of_a_target_are_dropped`, `reversed_decoys_of_modified_targets_mirror_sites_and_keep_mass`, `generated_decoys_of_n_terminal_peptides_stay_balanced_and_n_terminal` | FASTA-supplied decoys are not balance-checked (documented, by design). |
| HC-FDR-03 target-decoy q, not p-value corrections | `fdr.rs`: `count_confidence_is_identical_for_ties_and_keeps_plus_one`, `sparse_decoys_use_finite_count_based_confidence`, `empty_and_decoy_only_count_confidence_are_conservative`; `ml/qvalue.rs`: `equal_scores_receive_the_same_q_value`, `tied_score_order_does_not_change_q_values` | None known. |
| HC-FDR-04 filter on q, then report | `ptm.rs`: `target_decoy_q_values_are_monotonic`, `tied_target_arrangements_cannot_inherit_confident_competition_q_values`; `sage-cloudpath` `parquet.rs` writes `sage.output_filter.spectrum_q_max` | No test that a protein-site row's q is re-estimated after filtering (it is not). |
| HC-FDR-05 DIA FDR | `sage-dia` unit tests cover pseudo-spectrum construction only. | Gap: no DIA decoy model, peak-group scoring or DIA entrapment benchmark. |
| HC-FDR-06 MBR transfer FDR | `lfq.rs`: `strict_ms2_evidence_requires_both_psm_and_peptide_acceptance`, `disabling_mbr_quantifies_each_file_against_its_own_anchor`; `fdr.rs`: `picked_precursor_assigns_count_q_values_per_charge_state`; pilot in `benchmarks/scientific-results/20260914/` | Gap: no per-transfer q and no test that would fail on foreign-species transfers. |
| HC-QUANT-01 quantify after identification | `lfq.rs`: `strict_ms2_evidence_requires_both_psm_and_peptide_acceptance`, `one_identified_label_channel_seeds_all_channel_precursors` | None known. |
| HC-QUANT-02 missing values are not zero | `sage-cloudpath` `parquet.rs`: `lfq_preserves_missingness_and_ms2_evidence`; `lfq.rs`: `file_evidence_does_not_invent_a_score_for_a_missing_trace`; `tmt.rs`: `reporter_search_selects_the_most_intense_peak_and_marks_missing_channels` | Gap on this branch: TMT missing channels are written as 0.0; no Parquet test asserts null reporter values. |
| HC-QUANT-03 intensity semantics labelled | Docs only: `schemas/scores.v1.md`, DOCS.md "What quantification does not do". | No test (documentation rule). |
| HC-RPT-01 reproducible outputs | `integration.rs`: `parquet_footers_and_run_summary_record_provenance`; `sage-cli` unit `provenance.rs` (FASTA hash stable, spectrum hashing off by default); `sage-cloudpath` `parquet.rs`: `provenance_follows_the_schema_keys_in_the_footer`; `config_schema.rs`: `committed_config_schema_matches_rust_types` | The git commit omits uncommitted changes. Spectrum hashes are opt-in. |
| HC-RPT-02 score definitions and directions | Docs: DOCS.md PSM column list and `schemas/scores.v1.md`. `sage-cloudpath` `parquet.rs` round-trip tests pin column names. | No test that DOCS lists every schema column. |
| HC-EFF-01 complexity stated | Docs: DOCS.md "Performance and complexity". Measured runtime and RSS in `benchmarks/RESULTS.md`, `benchmarks/large_db`. | The complexity statements are analytical; no scaling benchmark over N or S. |
| HC-EFF-02 bounded mod combinatorics | `database.rs`: `max_total_count_caps_library_and_max_count_caps_new_sites`, `memory_estimate_counts_per_modification_caps_exactly`, `estimates_variable_modification_expansion_before_allocation`, `preflight_variant_counts_match_generated_peptides` | None known. |
| HC-EFF-03 indexed search | `database.rs`: `fragment_index_preserves_exact_ids_masses_and_ranges`, `interleaved_bucket_ranges_match_bucket_search_on_random_queries`; `integration.rs` prefilter equivalence tests (`ptm_library_sites_match_with_and_without_prefilter` and others) | None known. |
| HC-EFF-04 index built once | Structural (one `IndexedDatabase` per `Runner`). | No test asserts a single build. |
| HC-EFF-05 memory bounds | `sage-cli` unit `memory.rs`; `database.rs`: `memory_estimate_totals_are_additive_over_protein_chunks`; `benchmarks/large_db` | Peak RSS is recorded in benchmarks, not per run. |
| HC-INTER-01 open, versioned formats | `config_schema.rs`: `committed_config_schema_matches_rust_types`, `generated_config_schema_has_expected_contract`; Parquet schemas in `schemas/` are loaded by the writers and round-tripped in `sage-cloudpath` `parquet.rs` | Gap: `run-summary.json`, `results.json`, `digestion.tsv` and `diagnostic_ions.tsv` have no published schema. |
| HC-INTER-02 controlled vocabularies | `peff.rs` and `unimod.rs` unit tests (input side). | Gap: outputs carry no UNIMOD accessions or CV terms. |
| HC-INTER-03 database identity | `sage-cli` unit `provenance.rs`: `fasta_hash_is_stable_and_matches_sha256sum`, `headers_count_decoys_sources_and_organisms`; `integration.rs`: `parquet_footers_and_run_summary_record_provenance` | No release is recorded; UniProt FASTA headers carry none. |
| HC-INTER-04 documented API | `sage-cli` unit `api.rs`. | Error conditions are not systematically documented; no published rustdoc. |
| HC-TEST-01 rules mapped to tests | This file. | Failure-mode tests are missing for HC-QUANT-02 (TMT) and HC-FDR-06. |
| HC-TEST-02 FDR validated by entrapment | `benchmarks/archive/HARDENING_RESULTS.md` (FDRBench, 20 seeds, peptide level) | One file and a 2,000-protein subset. Protein groups, PTM sites, MBR transfers and DIA are not certified. |
| HC-TEST-03 quantification validated on known ratios | `benchmarks/scientific-results/20260914/SCIENTIFIC_REPORT.md` (three-species LFQ) | Run on Beta 3, not rerun since. No TMT known-ratio set. |
| HC-TEST-04 end-to-end on open data | `integration.rs` (open mzML and FASTA in `tests/`); CI synthetic smoke run in `.github/workflows/rust.yml` | None known. |
