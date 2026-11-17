//! `database.variant_fasta` through whole runs on the repository's
//! one-spectrum fixture, whose peptide LQSRPAAPPAPGPGQLTLR is put on a variant
//! protein only, so its PSM is a variant identification.

use super::*;

const PEPTIDE: &str = "LQSRPAAPPAPGPGQLTLR";

fn workspace() -> String {
    concat!(env!("CARGO_MANIFEST_DIR"), "/../..").to_string()
}

fn directory(name: &str) -> std::path::PathBuf {
    let unique = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .unwrap()
        .as_nanos();
    let directory = std::env::temp_dir().join(format!(
        "sage-variant-{name}-{}-{unique}",
        std::process::id()
    ));
    std::fs::create_dir_all(&directory).unwrap();
    directory
}

/// The configuration of a run over the fixture with `fasta` as the reference
/// proteome and, when given, `variant_fasta`.
fn input(
    directory: &std::path::Path,
    fasta: &str,
    variant_fasta: Option<&str>,
) -> crate::input::Input {
    input_with(
        directory,
        fasta,
        variant_fasta,
        serde_json::json!({}),
        serde_json::json!({}),
    )
}

/// [`input`] with `database` keys and top-level keys added.
fn input_with(
    directory: &std::path::Path,
    fasta: &str,
    variant_fasta: Option<&str>,
    database_extra: serde_json::Value,
    extra: serde_json::Value,
) -> crate::input::Input {
    let reference = directory.join("reference.fasta");
    std::fs::write(&reference, fasta).unwrap();
    let mut database = serde_json::json!({
        "fasta": reference.to_string_lossy(),
        "enzyme": {"cleave_at": "KR", "restrict": "P", "missed_cleavages": 0, "min_len": 5, "max_len": 50, "semi_enzymatic": false},
        "static_mods": {},
        "ion_kinds": ["b", "y"]
    });
    if let Some(text) = variant_fasta {
        let path = directory.join("variants.fasta");
        std::fs::write(&path, text).unwrap();
        database["variant_fasta"] = serde_json::json!(path.to_string_lossy());
    }
    for (key, value) in database_extra.as_object().unwrap() {
        database[key] = value.clone();
    }
    let mut config = serde_json::json!({
        "database": database,
        "precursor_tol": { "ppm": [-50, 50] },
        "fragment_tol": { "ppm": [-10, 10] },
        "isotope_errors": [-1, 3],
        "min_peaks": 15,
        "min_matched_peaks": 4,
        "mzml_paths": [format!("{}/tests/LQSRPAAPPAPGPGQLTLR.mzML", workspace())],
        "predict_rt": false,
        "output_filter": {"psm_q_value": 1.0},
        "output_directory": directory.join("out").to_string_lossy(),
    });
    for (key, value) in extra.as_object().unwrap() {
        config[key] = value.clone();
    }
    serde_json::from_value(config).unwrap()
}

const REFERENCE: &str = ">sp|P1|REFERENCE\nMSAGGWGGKAAAEELLGGK\n";

fn schema_version(path: &std::path::Path) -> String {
    use parquet::file::reader::{FileReader, SerializedFileReader};
    let reader = SerializedFileReader::new(std::fs::File::open(path).unwrap()).unwrap();
    reader
        .metadata()
        .file_metadata()
        .key_value_metadata()
        .unwrap()
        .iter()
        .find(|entry| entry.key == "sage.schema.version")
        .and_then(|entry| entry.value.clone())
        .unwrap()
}

/// Each row as a map from column name to its value, as text.
fn rows(path: &std::path::Path) -> Vec<HashMap<String, String>> {
    use parquet::file::reader::{FileReader, SerializedFileReader};
    let reader = SerializedFileReader::new(std::fs::File::open(path).unwrap()).unwrap();
    reader
        .get_row_iter(None)
        .unwrap()
        .map(|row| {
            row.unwrap()
                .get_column_iter()
                .map(|(name, value)| (name.clone(), value.to_string()))
                .collect()
        })
        .collect()
}

#[test]
fn a_variant_protein_that_reuses_an_accession_is_refused_before_any_spectrum_is_read() {
    let directory = directory("shared-accession");
    let input = input(
        &directory,
        REFERENCE,
        Some(&format!(">sp|P1|REFERENCE\nGGWGGK{PEPTIDE}\n")),
    );
    let error = match Runner::new(input.build().unwrap(), 1) {
        Ok(_) => panic!("a variant protein with a reference accession must be refused"),
        Err(error) => format!("{error:#}"),
    };
    assert!(
        error.contains("protein `sp|P1|REFERENCE` is already in another FASTA"),
        "{error}"
    );
    assert!(error.contains("variants.fasta"), "{error}");
    assert!(!directory.join("out").join("results.sage.parquet").exists());
    std::fs::remove_dir_all(directory).ok();
}

#[test]
fn a_whole_run_marks_the_variant_peptide_and_keeps_it_out_of_protein_inference() {
    // Off: the version 7 file, without variant columns or a summary block.
    let off_directory = directory("off");
    let off = Runner::new(
        input(
            &off_directory,
            &format!("{REFERENCE}>sp|P2|CARRIER\nGGWGGK{PEPTIDE}\n"),
            None,
        )
        .build()
        .unwrap(),
        1,
    )
    .unwrap()
    .run(1)
    .unwrap();
    assert!(off.variants.is_none());
    let off_results = off_directory.join("out").join("results.sage.parquet");
    assert_eq!(schema_version(&off_results), "7");
    let off_rows = rows(&off_results);
    assert!(!off_rows.is_empty());
    assert!(off_rows.iter().all(|row| !row.contains_key("is_variant")));

    // On: the fixture peptide is only on a variant protein.
    let on_directory = directory("on");
    let on = Runner::new(
        input(
            &on_directory,
            REFERENCE,
            Some(&format!(">sp|P2_V|CARRIER_VARIANT\nGGWGGK{PEPTIDE}\n")),
        )
        .build()
        .unwrap(),
        1,
    )
    .unwrap()
    .run(1)
    .unwrap();
    let on_results = on_directory.join("out").join("results.sage.parquet");
    assert_eq!(schema_version(&on_results), "11");
    let on_rows = rows(&on_results);
    let variant = on_rows
        .iter()
        .find(|row| row["peptide"].trim_matches('"') == PEPTIDE)
        .expect("the fixture spectrum is matched to its peptide");
    assert_eq!(variant["is_variant"], "true");
    assert_ne!(variant["variant_spectrum_q"], "null");
    assert_ne!(variant["variant_peptide_q"], "null");
    assert_eq!(variant["in_source_rescued"], "null");
    // Reported as kept out of protein inference. On one spectrum these values
    // hold either way; the exclusion itself is tested in sage-core's
    // `a_variant_decoy_cannot_raise_the_reference_proteins_q_values`.
    assert_eq!(variant["protein_q"], "1.0");
    assert_eq!(variant["protein_group_q"], "1.0");
    assert_eq!(
        variant["protein_groups"].trim_matches('"'),
        "sp|P2_V|CARRIER_VARIANT"
    );

    let stats = on.variants.expect("a variant FASTA was searched");
    assert_eq!(stats.variant_proteins, 1);
    // GGWGGK and the fixture peptide, both only on the variant protein.
    assert_eq!(stats.variant_peptides, 2);
    assert_eq!(stats.shared_peptides, 0);
    assert_eq!(stats.unpaired_peptides, 0);
    // One spectrum, one target PSM of the class, no decoy: q is the floor 1/1.
    let one = sage_core::variant::VariantClassCounts {
        targets: 1,
        decoys: 0,
        passing_variant_q: 0,
        passing_global_q: 0,
        q_floor: Some(1.0),
        lowest_q: Some(1.0),
    };
    assert_eq!(stats.psms, one);
    assert_eq!(stats.peptides, one);
    // One spectrum either way, so as many rows; byte identity of a run with the
    // option off is checked against main on real data with two binaries.
    assert_eq!(on_rows.len(), off_rows.len());
    std::fs::remove_dir_all(off_directory).ok();
    std::fs::remove_dir_all(on_directory).ok();
}

#[test]
fn a_variant_fasta_without_a_target_protein_is_refused() {
    // Every entry carries the decoy tag, so with generated decoys none is
    // read: the run must not quietly search as if the option were unset.
    let directory = directory("empty-variants");
    let input = input(
        &directory,
        REFERENCE,
        Some(&format!(">rev_sp|P2_V|CARRIER_VARIANT\nGGWGGK{PEPTIDE}\n")),
    );
    let error = match Runner::new(input.build().unwrap(), 1) {
        Ok(_) => panic!("a variant FASTA with no target protein must be refused"),
        Err(error) => format!("{error:#}"),
    };
    assert!(error.contains("Failed to read variant FASTA"), "{error}");
    std::fs::remove_dir_all(directory).ok();
}

#[test]
fn a_peptide_the_reference_explains_is_reported_on_the_reference_protein_only() {
    // The fixture peptide is in the reference protein and in its variant copy,
    // which differs elsewhere: a reference peptide, on the reference protein.
    let directory = directory("reference-mapped");
    let summary = Runner::new(
        input(
            &directory,
            &format!(">sp|P2|CARRIER\nGGWGGK{PEPTIDE}\n"),
            Some(&format!(">sp|P2_V|CARRIER_G1A\nAGWGGK{PEPTIDE}\n")),
        )
        .build()
        .unwrap(),
        1,
    )
    .unwrap()
    .run(1)
    .unwrap();
    let rows = rows(&directory.join("out").join("results.sage.parquet"));
    let row = rows
        .iter()
        .find(|row| row["peptide"].trim_matches('"') == PEPTIDE)
        .expect("the fixture spectrum is matched to its peptide");
    assert_eq!(row["is_variant"], "false");
    assert_eq!(row["variant_spectrum_q"], "null");
    assert_eq!(row["proteins"].trim_matches('"'), "sp|P2|CARRIER");
    assert_eq!(row["num_proteins"], "1");
    let stats = summary.variants.expect("a variant FASTA was searched");
    assert!(stats.shared_peptides >= 1);
    std::fs::remove_dir_all(directory).ok();
}

/// A whole run with the fixture peptide on a variant protein only, and the
/// given extra settings; returns the summary and the results rows.
fn variant_run(
    name: &str,
    database_extra: serde_json::Value,
    extra: serde_json::Value,
) -> (
    crate::runner::RunSummary,
    String,
    Vec<HashMap<String, String>>,
) {
    let directory = directory(name);
    let summary = Runner::new(
        input_with(
            &directory,
            REFERENCE,
            Some(&format!(">sp|P2_V|CARRIER_VARIANT\nGGWGGK{PEPTIDE}\n")),
            database_extra,
            extra,
        )
        .build()
        .unwrap(),
        1,
    )
    .unwrap()
    .run(1)
    .unwrap();
    let results = directory.join("out").join("results.sage.parquet");
    let version = schema_version(&results);
    let rows = rows(&results);
    std::fs::remove_dir_all(directory).ok();
    (summary, version, rows)
}

fn variant_row(rows: &[HashMap<String, String>]) -> &HashMap<String, String> {
    rows.iter()
        .find(|row| row["peptide"].trim_matches('"') == PEPTIDE)
        .expect("the fixture spectrum is matched to its peptide")
}

#[test]
fn the_variant_class_is_the_same_with_the_prefilter() {
    let (plain, _, plain_rows) = variant_run("plain", serde_json::json!({}), serde_json::json!({}));
    let (prefiltered, version, rows) = variant_run(
        "prefilter",
        serde_json::json!({"prefilter": true, "prefilter_min_matched_peaks": 1}),
        serde_json::json!({}),
    );
    assert_eq!(version, "11");
    assert_eq!(variant_row(&rows)["is_variant"], "true");
    assert_eq!(
        variant_row(&rows)["variant_spectrum_q"],
        variant_row(&plain_rows)["variant_spectrum_q"]
    );
    assert_eq!(
        prefiltered.variants.unwrap().psms,
        plain.variants.unwrap().psms
    );
}

#[test]
fn precursor_labels_write_schema_12_with_the_variant_columns() {
    let (summary, version, rows) = variant_run(
        "labels",
        serde_json::json!({"static_mods": {"Arg10": {
            "mass": 0.0,
            "channel_offsets": {"light": 0.0, "heavy": 10.008269},
            "sites": ["R"]
        }}}),
        serde_json::json!({}),
    );
    assert_eq!(version, "12");
    let row = variant_row(&rows);
    assert_eq!(row["is_variant"], "true");
    assert_ne!(row["variant_spectrum_q"], "null");
    assert_eq!(summary.variants.unwrap().psms.targets, 1);
}

#[test]
fn with_the_rescue_pass_on_the_rescue_column_is_written_beside_the_variant_columns() {
    let (_, version, rows) = variant_run(
        "rescue",
        serde_json::json!({}),
        serde_json::json!({"in_source": {"rescue": true}}),
    );
    assert_eq!(version, "11");
    let row = variant_row(&rows);
    // A search row: not rescued, classed.
    assert_eq!(row["in_source_rescued"], "false");
    assert_eq!(row["is_variant"], "true");
}

#[test]
fn entrapment_and_contaminant_fastas_keep_their_own_columns_beside_the_variants() {
    let directory = directory("classes");
    let entrapment = directory.join("entrapment.fasta");
    std::fs::write(&entrapment, ">AT1G01010\nMNEEKAWVVRSTTPPHGGQK\n").unwrap();
    let contaminant = directory.join("contaminant.fasta");
    std::fs::write(&contaminant, ">CON_P02768\nMKWVTFISLLFLFSSAYSR\n").unwrap();
    let (summary, version, rows) = variant_run(
        "classes",
        serde_json::json!({
            "entrapment_fasta": entrapment.to_string_lossy(),
            "contaminant_fasta": contaminant.to_string_lossy()
        }),
        serde_json::json!({}),
    );
    assert_eq!(version, "11");
    let row = variant_row(&rows);
    assert_eq!(row["is_variant"], "true");
    assert_eq!(row["is_entrapment"], "false");
    assert_eq!(row["is_contaminant"], "false");
    assert!(summary.entrapment.is_some());
    assert!(summary.contaminants.is_some());
    let variants = summary.variants.unwrap();
    assert_eq!(variants.variant_proteins, 1);
    std::fs::remove_dir_all(directory).ok();
}
