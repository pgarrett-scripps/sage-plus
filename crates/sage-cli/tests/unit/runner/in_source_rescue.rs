//! `Runner::rescue_file` and `Runner::in_source_rescue` on the repository's
//! one-spectrum fixture, whose peptide LQSRPAAPPAPGPGQLTLR is made the suffix
//! product of a longer parent that the database holds while the product is left
//! out of it, as a fully enzymatic search would leave out a truncation.

use super::*;
use sage_core::in_source::{InSourceAnnotations, InSourceSettings, Placement};
use sage_core::in_source_rescue::InSourceRescueFileStats;

const PRODUCT: &str = "LQSRPAAPPAPGPGQLTLR";
/// The parent: one missed cleavage, the product as its C-terminal piece.
const PARENT: &str = "GGWGGKLQSRPAAPPAPGPGQLTLR";

fn workspace() -> String {
    concat!(env!("CARGO_MANIFEST_DIR"), "/../..").to_string()
}

fn fixture() -> String {
    format!("{}/tests/LQSRPAAPPAPGPGQLTLR.mzML", workspace())
}

/// A runner over `files` copies of the fixture (or `paths` when given), whose
/// database holds the parent and, with `keep_product`, the product too.
fn fixture_runner(paths: &[String], keep_product: bool) -> (Runner, std::path::PathBuf) {
    let unique = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .unwrap()
        .as_nanos();
    let directory =
        std::env::temp_dir().join(format!("sage-rescue-test-{}-{unique}", std::process::id()));
    std::fs::create_dir_all(&directory).unwrap();
    let fasta = directory.join("parent.fasta");
    std::fs::write(&fasta, format!(">sp|P1|PARENT\n{PARENT}\n")).unwrap();
    let input: crate::input::Input = serde_json::from_value(serde_json::json!({
        "database": {
            "fasta": fasta.to_string_lossy(),
            "enzyme": {"cleave_at": "KR", "restrict": "P", "missed_cleavages": 1, "min_len": 5, "max_len": 50, "semi_enzymatic": false},
            "static_mods": {},
            "ion_kinds": ["b", "y"]
        },
        "precursor_tol": { "ppm": [-50, 50] },
        "fragment_tol": { "ppm": [-10, 10] },
        "isotope_errors": [-1, 3],
        "deisotope": true,
        "min_peaks": 15,
        "min_matched_peaks": 4,
        "mzml_paths": paths,
        "in_source": {"rescue": true},
        "output_directory": directory.to_string_lossy(),
    }))
    .unwrap();
    let mut runner = Runner::new(input.build().unwrap(), 1).unwrap();
    // The main database: the full digest, less the product unless kept.
    let fasta =
        sage_core::fasta::Fasta::parse(std::fs::read_to_string(&fasta).unwrap(), "rev_", true)
            .unwrap();
    let peptides = runner
        .database_parameters
        .digest(&fasta)
        .into_iter()
        .filter(|peptide| {
            keep_product || {
                let sequence = std::str::from_utf8(peptide.sequence.as_bytes()).unwrap();
                sequence != PRODUCT && sequence != reversed_interior(PRODUCT)
            }
        })
        .collect::<Vec<_>>();
    runner.database = runner
        .database_parameters
        .clone()
        .build_from_peptides(peptides);
    (runner, directory)
}

fn reversed_interior(sequence: &str) -> String {
    let bytes = sequence.as_bytes();
    let mut out = vec![bytes[0]];
    out.extend(bytes[1..bytes.len() - 1].iter().rev());
    out.push(bytes[bytes.len() - 1]);
    String::from_utf8(out).unwrap()
}

/// The fixture spectrum's id and retention time.
fn fixture_spectrum() -> (String, f32) {
    let spectra =
        sage_cloudpath::util::read_mzml(&sage_cloudpath::to_url(&fixture()).unwrap(), 0, None)
            .unwrap();
    (spectra[0].id.clone(), spectra[0].scan_start_time)
}

/// A confident PSM of the parent in `file_id`, on another spectrum, at `rt`.
fn parent_psm(runner: &Runner, file_id: usize, psm_id: usize, rt: f32) -> Feature {
    let parent = (0..runner.database.peptides.len())
        .map(|index| sage_core::database::PeptideIx(index as u32))
        .find(|&index| {
            let peptide = &runner.database[index];
            !peptide.decoy && peptide.sequence.as_bytes() == PARENT.as_bytes()
        })
        .expect("the parent is in the database");
    Feature {
        psm_id,
        peptide_idx: parent,
        file_id,
        spec_id: format!("parent-scan-{psm_id}"),
        rt,
        label: 1,
        rank: 1,
        spectrum_q: 0.001,
        peptide_q: 0.001,
        ..Feature::default()
    }
}

fn settings() -> InSourceSettings {
    InSourceSettings {
        rt_tolerance: Some(0.5),
        apex_fwhm_fraction: None,
        rescue: Some(true),
    }
}

#[test]
fn an_unidentified_spectrum_of_a_coeluting_parents_truncation_is_rescued() {
    let (runner, directory) = fixture_runner(&[fixture()], false);
    let (_, rt) = fixture_spectrum();
    let features = [parent_psm(&runner, 0, 0, rt)];
    let (hits, stats) = runner.rescue_file(
        0,
        &features,
        &settings(),
        &InSourceAnnotations::default(),
        None,
    );
    assert_eq!(stats.parents, 1);
    assert_eq!(stats.searched_spectra, 1);
    assert!(stats.candidates > 0);
    assert_eq!(stats.error, None);
    assert_eq!(hits.len(), 1, "{stats:?}");
    let (feature, origin, peptide) = &hits[0];
    assert_eq!(peptide.sequence.as_bytes(), PRODUCT.as_bytes());
    assert!(!peptide.decoy);
    assert_eq!(origin.parent, features[0].peptide_idx);
    assert_eq!(origin.kind.placement, Placement::Suffix);
    assert_eq!(origin.parent_len, PARENT.len());
    // A lone target competes with no decoy: (0 + 1) / 1.
    assert_eq!(feature.spectrum_q, 1.0);
    assert_eq!((feature.peptide_q, feature.protein_q), (1.0, 1.0));
    assert_eq!(feature.label, 1);
    assert_eq!(stats.target_hits, 1);
    assert_eq!(stats.decoy_hits, 0);
    // Missed cleavages are the product's own: none, not the parent's one.
    assert_eq!(feature.missed_cleavages, 0);
    std::fs::remove_dir_all(directory).ok();
}

#[test]
fn a_parent_that_does_not_coelute_rescues_nothing() {
    let (runner, directory) = fixture_runner(&[fixture()], false);
    let (_, rt) = fixture_spectrum();
    for parent_rt in [rt + 5.0, rt - 5.0] {
        let features = [parent_psm(&runner, 0, 0, parent_rt)];
        let (hits, stats) = runner.rescue_file(
            0,
            &features,
            &settings(),
            &InSourceAnnotations::default(),
            None,
        );
        assert!(hits.is_empty(), "{stats:?}");
        assert_eq!(stats.searched_spectra, 1);
    }
    // The null control's shift moves the window off the spectrum as well.
    let features = [parent_psm(&runner, 0, 0, rt)];
    let (hits, stats) = runner.rescue_file(
        0,
        &features,
        &settings(),
        &InSourceAnnotations::default(),
        Some(5.0),
    );
    assert!(hits.is_empty());
    assert_eq!(stats.null_shift_widths, Some(5.0));
    std::fs::remove_dir_all(directory).ok();
}

#[test]
fn identified_spectra_and_database_peptides_are_not_searched_again() {
    // A spectrum with a confident PSM is not searched.
    let (runner, directory) = fixture_runner(&[fixture()], false);
    let (id, rt) = fixture_spectrum();
    let mut identified = parent_psm(&runner, 0, 1, rt);
    identified.spec_id = id;
    let features = [parent_psm(&runner, 0, 0, rt), identified];
    let (hits, stats) = runner.rescue_file(
        0,
        &features,
        &settings(),
        &InSourceAnnotations::default(),
        None,
    );
    assert!(hits.is_empty());
    assert_eq!(stats.searched_spectra, 0);
    std::fs::remove_dir_all(directory).ok();

    // A truncation the database already holds is not a candidate: the main
    // search scored it.
    let (with_product, directory) = fixture_runner(&[fixture()], true);
    let (without_product, other) = fixture_runner(&[fixture()], false);
    let candidates = |runner: &Runner| {
        let parent = parent_psm(runner, 0, 0, rt).peptide_idx;
        let mut parents = fnv::FnvHashMap::default();
        parents.insert(parent, rt);
        let mut stats = InSourceRescueFileStats::default();
        let (origins, peptides) = runner.rescue_candidates(&parents, &mut stats);
        (origins, peptides, stats)
    };
    let (origins, peptides, stats) = candidates(&with_product);
    let (_, all_peptides, all_stats) = candidates(&without_product);
    assert_eq!(stats.candidates + 1, all_stats.candidates);
    assert!(!origins.contains_key(&(PRODUCT.to_string(), false)));
    // Every target keeps its decoy.
    for list in [&peptides, &all_peptides] {
        assert_eq!(
            list.iter().filter(|p| p.decoy).count(),
            list.iter().filter(|p| !p.decoy).count()
        );
    }
    std::fs::remove_dir_all(directory).ok();
    std::fs::remove_dir_all(other).ok();
}

#[test]
fn a_file_that_cannot_be_read_again_is_recorded_as_an_error() {
    // The file exists when the search starts and is gone when the rescue pass
    // reads it again.
    let staging = std::env::temp_dir().join(format!(
        "sage-rescue-vanishing-{}-{}.mzML",
        std::process::id(),
        std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .unwrap()
            .as_nanos()
    ));
    std::fs::copy(fixture(), &staging).unwrap();
    let (runner, directory) = fixture_runner(&[staging.to_string_lossy().into_owned()], false);
    std::fs::remove_file(&staging).unwrap();
    let features = [parent_psm(&runner, 0, 0, 10.0)];
    let (hits, stats) = runner.rescue_file(
        0,
        &features,
        &settings(),
        &InSourceAnnotations::default(),
        None,
    );
    assert!(hits.is_empty());
    let reason = stats.error.expect("the file is recorded as failed");
    assert!(reason.contains("sage-rescue-vanishing"), "{reason}");
    std::fs::remove_dir_all(directory).ok();
}

#[test]
fn rescued_psms_of_several_files_share_one_registered_peptide_and_stay_out_of_the_search() {
    let (mut runner, directory) = fixture_runner(&[fixture(), fixture()], false);
    let (_, rt) = fixture_spectrum();
    let features = [parent_psm(&runner, 0, 0, rt), parent_psm(&runner, 1, 7, rt)];
    let before = features.clone();
    let indexed = runner.database.peptides.len();
    let mut annotations = InSourceAnnotations::default();
    let outcome = runner
        .in_source_rescue(&features, &settings(), &mut annotations, None)
        .unwrap();
    assert_eq!(outcome.error, None);
    let (rescued, stats) = (outcome.features, outcome.stats);
    assert_eq!(stats.len(), 2);
    assert_eq!(rescued.len(), 2);
    // One product, registered once, beside the indexed peptides.
    assert_eq!(rescued[0].peptide_idx, rescued[1].peptide_idx);
    assert_eq!(rescued[0].peptide_idx.0 as usize, indexed);
    assert_eq!(
        runner.database[rescued[0].peptide_idx].sequence.as_bytes(),
        PRODUCT.as_bytes()
    );
    // Numbered after every PSM of the search, each one recorded as rescued
    // against its parent, and the search's own PSMs untouched.
    assert_eq!(rescued[0].psm_id, 8);
    assert_eq!(rescued[1].psm_id, 9);
    assert_eq!(rescued[0].file_id, 0);
    assert_eq!(rescued[1].file_id, 1);
    for feature in &rescued {
        assert!(annotations.is_rescued(feature.psm_id));
        assert_eq!(
            annotations.psm(feature.psm_id).unwrap().parent,
            features[0].peptide_idx
        );
    }
    assert!(!annotations.is_rescued(0) && !annotations.is_rescued(7));
    assert_eq!(annotations.flagged_psms(), 0);
    for (a, b) in features.iter().zip(&before) {
        assert_eq!(
            (a.psm_id, a.peptide_idx, a.spectrum_q),
            (b.psm_id, b.peptide_idx, b.spectrum_q)
        );
    }
    std::fs::remove_dir_all(directory).ok();
}

#[test]
fn an_open_precursor_tolerance_is_refused() {
    // At the configuration...
    let workspace = workspace();
    let config = |tolerance: serde_json::Value, extra: serde_json::Value| {
        let mut config = serde_json::json!({
            "database": {
                "fasta": format!("{workspace}/tests/Q99536.fasta"),
                "enzyme": {"cleave_at": "KR", "restrict": "P", "missed_cleavages": 0, "semi_enzymatic": false},
                "static_mods": {},
                "ion_kinds": ["b", "y"]
            },
            "precursor_tol": tolerance,
            "fragment_tol": { "ppm": [-10, 10] },
            "mzml_paths": [fixture()],
            "in_source": {"rescue": true},
        });
        for (key, value) in extra.as_object().unwrap() {
            config[key] = value.clone();
        }
        serde_json::from_value::<crate::input::Input>(config)
            .unwrap()
            .build()
    };
    let open = config(
        serde_json::json!({"da": [-500, 100]}),
        serde_json::json!({}),
    );
    assert!(format!("{:#}", open.err().unwrap()).contains("precursor window"));
    let wide = config(
        serde_json::json!({"ppm": [-10, 10]}),
        serde_json::json!({"wide_window": true}),
    );
    assert!(format!("{:#}", wide.err().unwrap()).contains("DDA spectra only"));
    assert!(config(serde_json::json!({"ppm": [-10, 10]}), serde_json::json!({})).is_ok());

    // ...and by the pass itself, for a runner built some other way.
    let (mut runner, directory) = fixture_runner(&[fixture()], false);
    runner.parameters.precursor_tol = sage_core::mass::Tolerance::Da(-500.0, 100.0);
    let (_, rt) = fixture_spectrum();
    let features = [parent_psm(&runner, 0, 0, rt)];
    let outcome = runner
        .in_source_rescue(
            &features,
            &settings(),
            &mut InSourceAnnotations::default(),
            None,
        )
        .unwrap();
    assert!(outcome.features.is_empty());
    assert!(outcome.error.unwrap().contains("precursor window"));
    std::fs::remove_dir_all(directory).ok();
}

#[test]
fn a_failure_on_one_file_keeps_every_files_rescue_out() {
    let (mut runner, directory) = fixture_runner(&[fixture(), fixture()], false);
    let (_, rt) = fixture_spectrum();
    let features = [parent_psm(&runner, 0, 0, rt), parent_psm(&runner, 1, 7, rt)];
    let indexed = runner.database.peptides.len();
    let mut annotations = InSourceAnnotations::default();
    super::FAIL_REREAD.with(|fail| fail.set(Some(1)));
    let outcome = runner.in_source_rescue(&features, &settings(), &mut annotations, None);
    super::FAIL_REREAD.with(|fail| fail.set(None));
    let outcome = outcome.unwrap();
    // File 0 found its hit, but nothing of it is kept.
    assert!(outcome.features.is_empty());
    assert_eq!(outcome.stats.len(), 2);
    assert_eq!(outcome.stats[0].error, None);
    assert_eq!(outcome.stats[0].target_hits, 1);
    assert!(outcome.stats[1].error.is_some());
    let error = outcome.error.expect("the run is told the pass failed");
    assert!(error.contains("file 1"), "{error}");
    assert_eq!(runner.database.peptides.len(), indexed);
    assert!(!annotations.is_rescued(8) && !annotations.is_rescued(9));
    std::fs::remove_dir_all(directory).ok();
}

/// A whole run whose rescue fails on its second file writes the main results
/// as a run without the pass writes them, records the failure in
/// run-summary.json, and exits with an error naming the file.
#[test]
fn a_run_whose_rescue_fails_writes_the_main_results_and_exits_with_an_error() {
    let workspace = workspace();
    let run = |name: &str, rescue: bool, fail: Option<usize>| {
        let unique = std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .unwrap()
            .as_nanos();
        let directory = std::env::temp_dir().join(format!(
            "sage-rescue-run-{name}-{}-{unique}",
            std::process::id()
        ));
        let input: crate::input::Input = serde_json::from_value(serde_json::json!({
            "database": {
                "fasta": format!("{workspace}/tests/Q99536.fasta"),
                "enzyme": {"cleave_at": "KR", "restrict": "P", "missed_cleavages": 0, "semi_enzymatic": false},
                "static_mods": {"C": 57.021464},
                "ion_kinds": ["b", "y"]
            },
            "precursor_tol": { "ppm": [-50, 50] },
            "fragment_tol": { "ppm": [-10, 10] },
            "mzml_paths": [fixture(), fixture()],
            "in_source": {"rescue": rescue},
            "predict_rt": false,
            "output_directory": directory.to_string_lossy(),
        }))
        .unwrap();
        let runner = Runner::new(input.build().unwrap(), 1).unwrap();
        super::FAIL_REREAD.with(|cell| cell.set(fail));
        let result = runner.run(1);
        super::FAIL_REREAD.with(|cell| cell.set(None));
        (result, directory)
    };
    let rows = |directory: &std::path::Path| {
        use parquet::file::reader::{FileReader, SerializedFileReader};
        let file = std::fs::File::open(directory.join("results.sage.parquet")).unwrap();
        let reader = SerializedFileReader::new(file).unwrap();
        let metadata = reader.metadata().file_metadata();
        let fields = metadata
            .schema()
            .get_fields()
            .iter()
            .map(|field| field.name().to_string())
            .collect::<Vec<_>>();
        (metadata.num_rows(), fields)
    };
    let (without, without_directory) = run("off", false, None);
    without.unwrap();
    let (failed, failed_directory) = run("failed", true, Some(1));
    let error = format!("{:#}", failed.unwrap_err());
    assert!(error.contains("file 1"), "{error}");
    assert_eq!(rows(&failed_directory), rows(&without_directory));
    assert!(!rows(&failed_directory)
        .1
        .iter()
        .any(|f| f == "in_source_rescued"));
    let summary: serde_json::Value =
        serde_json::from_slice(&std::fs::read(failed_directory.join("run-summary.json")).unwrap())
            .unwrap();
    let recorded = summary["in_source"]["rescue_error"].as_str().unwrap();
    assert!(recorded.contains("file 1"), "{recorded}");
    assert!(summary["in_source"]["rescue"][1]["error"].is_string());
    // The same run without the failure writes the rescue column.
    let (completed, completed_directory) = run("completed", true, None);
    completed.unwrap();
    assert!(rows(&completed_directory)
        .1
        .iter()
        .any(|f| f == "in_source_rescued"));
    for directory in [without_directory, failed_directory, completed_directory] {
        std::fs::remove_dir_all(directory).ok();
    }
}

#[test]
fn palindromic_candidates_are_dropped_with_their_decoys() {
    let (mut runner, directory) = fixture_runner(&[fixture()], false);
    // A parent whose every piece reads the same with its interior reversed.
    let fasta = sage_core::fasta::Fasta::parse(">P2\nAAAAAAAAK\n".into(), "rev_", true).unwrap();
    let mut peptides = runner.database_parameters.digest(&fasta);
    peptides.retain(|p| !p.decoy);
    runner.database = runner
        .database_parameters
        .clone()
        .build_from_peptides(peptides);
    let mut parents = fnv::FnvHashMap::default();
    parents.insert(sage_core::database::PeptideIx(0), 10.0);
    let mut stats = InSourceRescueFileStats::default();
    let (origins, peptides) = runner.rescue_candidates(&parents, &mut stats);
    assert_eq!(stats.candidates, 0);
    assert!(stats.dropped_decoy_collisions > 0);
    assert!(origins.is_empty() && peptides.is_empty());
    std::fs::remove_dir_all(directory).ok();
}
