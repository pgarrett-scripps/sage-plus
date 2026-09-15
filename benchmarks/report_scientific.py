#!/usr/bin/env python3
"""Write a reviewable scientific pilot report from the verified result summary."""

import argparse
import json
import statistics
from collections import defaultdict
from pathlib import Path

from provenance import sha256


def value(number, digits=2):
    return "undefined" if number is None else f"{number:.{digits}f}"


def percent(number):
    return "undefined" if number is None else f"{100 * number:.2f}%"


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("summary", type=Path)
    p.add_argument("--output", type=Path, required=True)
    args = p.parse_args()
    d = json.loads(args.summary.read_text())
    environment = d.get("environment", {})
    lines = ["# Sage Plus scientific pilot", "",
             "The released beta.3 binary and upstream Sage commit df9219951cc9a54cf4cd55d76541af24b687bd3d are frozen. These measurements assess the selected public pilot and do not certify production calibration.", "",
             f"Recorded jobs: {d['completed_jobs']} complete, {d['failed_or_invalid_jobs']} failed or invalid, {d.get('running_jobs', 0)} running.", "",
             f"Host CPU: {environment.get('cpu_model', 'unavailable')}. Logical processors reported: {environment.get('logical_processors_reported', 'unavailable')}. Search worker limits are recorded per job.", "",
             "## Paired engineering comparisons", "",
             "Counts come from the underlying PSM files. These local HEK timing observations can overlap acquisition work and must not be presented as controlled final timing experiments.", "",
             "| Workload | Engine | Measured trials | Median seconds | Median peak MiB | Median target PSMs at 1% | Median peptidoforms at 1% |",
             "| --- | --- | ---: | ---: | ---: | ---: | ---: |"]
    grouped = defaultdict(list)
    for row in d["jobs"]:
        if row["suite"] == "local-paired" and not row["warmup"] and row["status"] == "complete":
            grouped[row["workload"], row["engine"]].append(row)
    for (workload, engine), rows in sorted(grouped.items()):
        med = lambda k: statistics.median(r[k] for r in rows)
        lines.append(f"| {workload} | {engine} | {len(rows)} | {med('wall_seconds'):.2f} | {med('peak_rss_mib'):.1f} | {med('target_psms'):.0f} | {med('target_peptidoforms'):.0f} |")
    lines += ["", "## Public identification agreement", "",
              "Matched files use identical search configurations. Overlap is measured among target PSMs passing each engine's reported 1% spectrum q-value threshold.", "",
              "The original mixed reference contains seven E. coli proteins with unknown X residues and is rejected by beta.3. The v2 comparison, mixed-reference entrapment, quantification and repeated timing use a documented variant excluding those same seven proteins for both engines. Their known subsequences are also excluded. The original reference and failed attempts are retained.", "",
              "| Input pair | Shared PSMs | Upstream only | Plus only | Jaccard | Status |",
              "| --- | ---: | ---: | ---: | ---: | --- |"]
    for result in d.get("public_identifications", []):
        if result["status"] == "complete":
            lines.append(f"| {result['pair']} | {result['shared_target_psms']} | {result['baseline_only']} | {result['candidate_only']} | {value(result['jaccard'], 4)} | complete |")
        else:
            lines.append(f"| {result['pair']} | unavailable | unavailable | unavailable | unavailable | incomplete |")
    lines += ["", "These are single paired public-file observations. They do not estimate timing variability.", "",
              "| Public search | Seconds | Peak MiB | Target PSMs at 1% | Peptidoforms at 1% |",
              "| --- | ---: | ---: | ---: | ---: |"]
    for row in d["jobs"]:
        if row["suite"].startswith("public-comparison") and row["status"] == "complete":
            lines.append(f"| {row['suite']}/{row['id']} | {value(row['wall_seconds'])} | {value(row['peak_rss_mib'], 1)} | {row['target_psms']} | {row['target_peptidoforms']} |")
    lines += ["", "## Repeated public timing", "",
              "One selected file per study, one warmup and three measured trials per engine. These runs follow acquisition, conversion and the primary search matrices. Parentheses show the observed minimum and maximum, not a confidence interval.", "",
              "| Study | Engine | Trials | Median seconds (range) | Median peak MiB | Median target PSMs at 1% |",
              "| --- | --- | ---: | ---: | ---: | ---: |"]
    repeated = defaultdict(list)
    for row in d["jobs"]:
        if row["suite"] == "public-timing" and row["status"] == "complete" and not row["warmup"]:
            repeated[row["study"], row["engine"]].append(row)
    for (study, engine), rows in sorted(repeated.items()):
        times = [r["wall_seconds"] for r in rows]
        lines.append(f"| {study} | {engine} | {len(rows)} | {statistics.median(times):.2f} ({min(times):.2f}, {max(times):.2f}) | {statistics.median(r['peak_rss_mib'] for r in rows):.1f} | {statistics.median(r['target_psms'] for r in rows):.0f} |")
    lines += ["", "## Independent entrapment pilot", "",
              "The table reports each file and construction seed. FDP estimators describe realized discovery sets. A small estimate or similarity between engines alone is not proof of valid FDR control.", "",
              "| Experiment | File and engine | Target peptides | Entrapments | Conservative paired FDP at nominal 1% | Seconds | Peak MiB |",
              "| --- | --- | ---: | ---: | ---: | ---: | ---: |"]
    measurements = {(r["suite"], r["id"]): r for r in d["jobs"]}
    for result in d["calibration"]:
        point = next(p for p in result["thresholds"] if p["nominal_q"] == 0.01)
        estimate = point["paired_fdp_tie_max"]
        measurement = measurements.get((result["suite"], result["job"]), {})
        lines.append(f"| {result['suite']} | {result['job']} | {point['targets']} | {point['entrapments']} | {percent(estimate)} | {value(measurement.get('wall_seconds'))} | {value(measurement.get('peak_rss_mib'), 1)} |")
    lines += ["", "Conditional resampling summaries:", ""]
    for result in d["calibration_uncertainty"]:
        if result["status"] != "descriptive_pilot_interval":
            lines.append(f"- {result['study']}: incomplete or undefined evidence, no interval reported.")
        else:
            lo, hi = result["percentile_intervals"]["paired_difference"]
            mean = result["means"]["plus"] - result["means"]["upstream"]
            lines.append(f"- {result['study']}: mean candidate-minus-upstream difference {mean * 100:.3f} percentage points, descriptive interval [{lo * 100:.3f}, {hi * 100:.3f}]. Two files and three shared seeds do not support broad generalization.")
    lines += ["", "## Exact prefilter", ""]
    checks = d["exact_prefilter"]
    equal = sum(c["all_scores_and_qvalues_equal"] for c in checks)
    lines.append(f"Exact PSM, score and q-value equality passed in {equal} of {len(checks)} completed measured pairs. Missing or failed pairs are not counted as passes.")
    lines += ["", "## PTM synthesis consistency", "",
              "These are site events passing reported 1% PSM and localization thresholds. The synthesis database is restricted, file-to-library mapping is not independently established, and the searches trigger heuristic rescoring. The table explicitly retains the peptide-level acceptance limitation.", "",
              "| Search | Consistent sites | Inconsistent sites | Inconsistent fraction | Site rows with peptide q = 1 |",
              "| --- | ---: | ---: | ---: | ---: |"]
    for r in d["ptm"]:
        error = r["empirical_site_error_fraction"]
        q = r["peptide_q_one_fraction"]
        lines.append(f"| {r['job']} | {r['correct_site_events']} | {r['incorrect_site_events']} | {percent(error)} | {percent(q)} |")
    joint = [r["joint_psm_peptide_localization_1pct"] for r in d["ptm"]]
    if joint and all(sum(r[k] for k in ("correct_site_events", "incorrect_site_events", "unassessable_site_events")) == 0 for r in joint):
        lines += ["", "No site events survive the joint 1% PSM, peptide and localization thresholds in this restricted synthetic pilot. Error among that empty accepted set is undefined."]
    lines += ["", "Oracle coverage is derived from synthesis truth. It measures sensitivity to missing known sites and cannot demonstrate the benefit of an independently annotated biological site library. The site library constrains peptide generation, but the released localizer reconsiders all residue-compatible positions. Full oracle coverage therefore does not force synthesis-consistent localization.", "",
              "## Known-ratio quantification", "",
              "B/A ground truth is human 1, yeast 0.5 and E. coli 4. Ratios pair Alpha and Beta preparations without imputation. CVs measure preparation variability, not technical-repeat precision.", "",
              "| Search | Species | Ratio pairs | Median log2 bias | Median absolute log2 error | Median preparation CV | Missing fraction |",
              "| --- | --- | ---: | ---: | ---: | ---: | ---: |"]
    for result in d["quantification"]:
        for species, r in result["species"].items():
            lines.append(f"| {result['job']} | {species} | {r['ratio_pairs']} | {value(r['median_log2_bias'], 3)} | {value(r['median_absolute_log2_error'], 3)} | {value(r['median_preparation_cv'], 3)} | {value(r['missing_fraction_observed_union'], 6)} |")
    for result in d["quantification"]:
        lines.append("")
        lines.append(f"Pure-human control, {result['job']}: {result['absent_species_quantified']} quantified foreign-species rows. Among {result['pure_human_control_without_ms2']} quantified rows without direct MS2 confirmation, {result['absent_species_without_ms2']} are assigned exclusively to foreign species. I/L-indistinguishable matches to human are excluded. This measures the detectable foreign component of transfer error. Sample purity and reference completeness remain assumptions.")
        lines.append(f"An independent check requires both PSM and peptide q-values at most 1% for direct MS2 evidence: {result['foreign_without_jointly_accepted_ms2']} foreign-species rows among {result['control_without_jointly_accepted_ms2']} control rows lacking that evidence. The engine's exported MS2 flag uses its LFQ peptide threshold without a separate PSM q-value requirement, so the two counts need not agree.")
        lines.append("")
    if not d["quantification"]:
        lines.append("No completed quantification results are available.")
    if d["quantification"]:
        lines += ["", "The serializer repeats one precursor-peak q-value across all files. It does not report an individual transfer q-value. Foreign assignments in the control therefore cannot be interpreted as a direct calibration test at the global precursor unit. The engine's discovery counter uses 5%, while this pilot retains its primary 1% LFQ endpoint. The following counts expose that threshold distinction before species filtering.", "",
                  "| Search | LFQ threshold | Target precursors | Positive target file rows |",
                  "| --- | ---: | ---: | ---: |"]
        for result in d["quantification"]:
            for threshold, row in result["lfq_threshold_yields"].items():
                lines.append(f"| {result['job']} | {100 * float(threshold):.0f}% | {row['target_precursors']} | {row['positive_target_file_rows']} |")
    diagnostics = [r for r in d["quantification"] if "direct_ms2_ratio_diagnostic" in r]
    if diagnostics:
        lines += ["", "### Exploratory ratio diagnostic before LFQ filtering", "",
                  "This diagnostic was added after observing that MBR off accepts no features at the primary 1% LFQ threshold. It uses positive intensities supported by direct MS2 evidence passing both 1% PSM and peptide thresholds, with cross-species I/L ambiguity excluded. It deliberately omits the LFQ q-value filter and does not change the primary result.", "",
                  "| Search | Species | Ratio pairs | Median log2 bias | Median absolute log2 error |",
                  "| --- | --- | ---: | ---: | ---: |"]
        for result in diagnostics:
            for species, row in result["direct_ms2_ratio_diagnostic"]["species"].items():
                lines.append(f"| {result['job']} | {species} | {row['ratio_pairs']} | {value(row['median_log2_bias'], 3)} | {value(row['median_absolute_log2_error'], 3)} |")
    lines += ["", "## Failed and invalid jobs", "",
              "Original rejected configurations and resource failures are retained. Configuration repairs use separate output directories.", "",
              "The broad-PTM stress case exceeds the selected resource budgets. Plus estimates a 23.8 GiB additional modified-database peak against its 10 GiB guard, while upstream fails allocation under the 16 GiB address-space ceiling. Initial oracle configurations omit the required max_count field and are corrected in ptm-oracle-v2. The 2 GiB address-space stress case fails with the prefilter off. Original mixed-reference failures are the undefined-residue case described above.", ""]
    for row in d["jobs"]:
        if row["status"] not in ("complete", "running"):
            lines.append(f"- {row['suite']}/{row['id']}: {row['status']}.")
    lines += ["", "## Follow-up priorities", "",
              "- Establish and validate confidence for individual MBR transfers. Keep precursor-level confidence and transfer-level confidence distinct in outputs and manuscript claims.",
              "- Resolve the restricted PTM pilot's peptide acceptance failure, independently establish file-to-library truth, and test localization separately from identification and oracle candidate generation.",
              "- Freeze a separate evaluation after the pilot, including held-out learned-model comparisons and a sample-size justification. Preserve the enlarged-reference runtime costs and explicit reference exclusions in any broad performance or coverage claims.",
              "", "## Reproducibility and remaining scope", ""]
    lines.extend("- " + item for item in d["limitations"])
    lines += ["", "The base RT and final LDA fits reuse scored observations. Grouped additive PTM-offset folds do not make the entire pipeline cross-fitted. This requires a separate model-validation experiment before strong learned-model claims.", "",
              f"Machine-readable source SHA-256: `{sha256(args.summary)}`.",
              "The evidence bundle includes configs, result files, reference sequences, binary hashes, source snapshots, extraction scripts and source receipts. Raw and converted spectra remain on /data. No external deposition has been performed.", ""]
    figures = sorted((args.summary.parent / "figures").glob("*.svg"))
    if figures:
        lines += ["## Figures", ""]
        lines.extend(f"- [{path.stem.replace('-', ' ')}](figures/{path.name})" for path in figures)
        lines += ["", "Captions and image checksums are recorded in [figures.json](figures/figures.json).", ""]
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text("\n".join(lines))
    print(args.output)


if __name__ == "__main__":
    main()
