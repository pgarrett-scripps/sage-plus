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
    let excluded = if keep_product {
        Vec::new()
    } else {
        vec![PRODUCT.to_string()]
    };
    runner_with(
        paths,
        &[PARENT.to_string()],
        &excluded,
        serde_json::json!({}),
    )
}

/// A runner over `paths` whose database is the digest of one protein per
/// sequence in `proteins`, less `excluded` peptides and their decoys, with
/// `variable_mods` added to the configuration.
fn runner_with(
    paths: &[String],
    proteins: &[String],
    excluded: &[String],
    variable_mods: serde_json::Value,
) -> (Runner, std::path::PathBuf) {
    let unique = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .unwrap()
        .as_nanos();
    let directory =
        std::env::temp_dir().join(format!("sage-rescue-test-{}-{unique}", std::process::id()));
    std::fs::create_dir_all(&directory).unwrap();
    let fasta = directory.join("parent.fasta");
    let text = proteins
        .iter()
        .enumerate()
        .map(|(index, sequence)| format!(">sp|P{index}|PARENT{index}\n{sequence}\n"))
        .collect::<String>();
    std::fs::write(&fasta, text).unwrap();
    let input: crate::input::Input = serde_json::from_value(serde_json::json!({
        "database": {
            "fasta": fasta.to_string_lossy(),
            "enzyme": {"cleave_at": "KR", "restrict": "P", "missed_cleavages": 1, "min_len": 5, "max_len": 50, "semi_enzymatic": false},
            "static_mods": {},
            "variable_mods": variable_mods,
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
    // The main database: the full digest, less the excluded peptides.
    let fasta =
        sage_core::fasta::Fasta::parse(std::fs::read_to_string(&fasta).unwrap(), "rev_", true)
            .unwrap();
    let dropped = excluded
        .iter()
        .flat_map(|sequence| [sequence.clone(), reversed_interior(sequence)])
        .collect::<HashSet<_>>();
    let peptides = runner
        .database_parameters
        .digest(&fasta)
        .into_iter()
        .filter(|peptide| {
            !dropped.contains(std::str::from_utf8(peptide.sequence.as_bytes()).unwrap())
        })
        .collect::<Vec<_>>();
    let offsets = runner.database.mass_offsets.clone();
    runner.database = runner
        .database_parameters
        .clone()
        .build_from_peptides(peptides);
    assert_eq!(runner.database.mass_offsets.len(), offsets.len());
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
    psm_of(runner, PARENT, file_id, psm_id, rt)
}

/// A confident PSM of the database peptide `sequence`, on another spectrum.
fn psm_of(runner: &Runner, sequence: &str, file_id: usize, psm_id: usize, rt: f32) -> Feature {
    let parent = (0..runner.database.peptides.len())
        .map(|index| sage_core::database::PeptideIx(index as u32))
        .find(|&index| {
            let peptide = &runner.database[index];
            !peptide.decoy && peptide.sequence.as_bytes() == sequence.as_bytes()
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
        &HashMap::new(),
        None,
    );
    assert_eq!(stats.parents, 1);
    assert_eq!(stats.searched_spectra, 1);
    assert!(stats.candidates > 0);
    assert_eq!(stats.error, None);
    assert_eq!(hits.len(), 1, "{stats:?}");
    let ((feature, origin, peptide), occurrence) = &hits[0];
    assert_eq!(*occurrence, 0);
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
            &HashMap::new(),
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
        &HashMap::new(),
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
        &HashMap::new(),
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
        &HashMap::new(),
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
        .in_source_rescue(
            &features,
            &settings(),
            &mut annotations,
            None,
            &HashMap::new(),
        )
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
            &HashMap::new(),
        )
        .unwrap();
    assert!(outcome.features.is_empty());
    assert!(outcome.error.unwrap().message.contains("precursor window"));
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
    let outcome = runner.in_source_rescue(
        &features,
        &settings(),
        &mut annotations,
        None,
        &HashMap::new(),
    );
    super::FAIL_REREAD.with(|fail| fail.set(None));
    let outcome = outcome.unwrap();
    // File 0 found its hit, but nothing of it is kept, not even its counts.
    assert!(outcome.features.is_empty());
    assert!(outcome.stats.is_empty());
    assert!(outcome.occurrences.is_empty());
    let error = outcome.error.expect("the run is told the pass failed");
    assert_eq!(error.file_id, Some(1));
    assert!(error.to_string().contains("file 1"), "{error}");
    assert_eq!(runner.database.peptides.len(), indexed);
    assert!(!annotations.is_rescued(8) && !annotations.is_rescued(9));
    std::fs::remove_dir_all(directory).ok();
}

/// The parent of the whole-run tests: the fixture's peptide after a residue
/// that is not a cleavage site, so the product is no tryptic peptide and the
/// fully enzymatic search cannot find it.
const WHOLE_RUN_PARENT: &str = "GGWGGALQSRPAAPPAPGPGQLTLR";

struct WholeRun {
    result: anyhow::Result<crate::runner::RunSummary>,
    directory: std::path::PathBuf,
    events: Vec<serde_json::Value>,
}

/// A whole run over `files` copies of the fixture, with a confident PSM of
/// the parent injected in every file at the fixture spectrum's retention
/// time, the rescue pass `rescue` (the in-source annotation is on either way),
/// and the pass made to fail on file `fail`. Every row is written.
fn whole_run(name: &str, files: usize, rescue: bool, fail: Option<usize>) -> WholeRun {
    let unique = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .unwrap()
        .as_nanos();
    let directory = std::env::temp_dir().join(format!(
        "sage-rescue-run-{name}-{}-{unique}",
        std::process::id()
    ));
    std::fs::create_dir_all(&directory).unwrap();
    let fasta = directory.join("parent.fasta");
    std::fs::write(&fasta, format!(">sp|P1|PARENT\n{WHOLE_RUN_PARENT}\n")).unwrap();
    let output = directory.join("out");
    let (_, rt) = fixture_spectrum();
    // The fixture plus MS1 scans around its retention time with the parent's
    // isotope envelope at charge 3, so LFQ has something to quantify.
    let parent_mass = {
        let probe = fixture_runner(&[fixture()], false);
        std::fs::remove_dir_all(probe.1).ok();
        let fasta = sage_core::fasta::Fasta::parse(
            format!(">sp|P1|PARENT\n{WHOLE_RUN_PARENT}\n"),
            "rev_",
            true,
        )
        .unwrap();
        probe
            .0
            .database_parameters
            .digest(&fasta)
            .into_iter()
            .find(|p| !p.decoy && p.sequence.as_bytes() == WHOLE_RUN_PARENT.as_bytes())
            .unwrap()
            .monoisotopic
    };
    let spectra_path = directory.join("with-ms1.mzML");
    std::fs::write(&spectra_path, with_ms1_scans(rt, parent_mass, 3)).unwrap();
    let input: crate::input::Input = serde_json::from_value(serde_json::json!({
        "database": {
            "fasta": fasta.to_string_lossy(),
            "enzyme": {"cleave_at": "KR", "restrict": "P", "missed_cleavages": 0, "min_len": 5, "max_len": 50, "semi_enzymatic": false},
            "static_mods": {},
            "ion_kinds": ["b", "y"]
        },
        "precursor_tol": { "ppm": [-50, 50] },
        "fragment_tol": { "ppm": [-10, 10] },
        "isotope_errors": [-1, 3],
        "min_peaks": 15,
        "min_matched_peaks": 4,
        "mzml_paths": vec![spectra_path.to_string_lossy().into_owned(); files],
        "in_source": {"rescue": rescue, "rt_tolerance": 0.5},
        "predict_rt": false,
        "quant": {"lfq": true},
        "write_pin": true,
        "output_filter": {"psm_q_value": 1.0},
        "output_directory": output.to_string_lossy(),
    }))
    .unwrap();
    let events_path = directory.join("events.jsonl");
    let runner = Runner::new_with_control(
        input.build().unwrap(),
        1,
        crate::events::EventEmitter::from_writer(std::fs::File::create(&events_path).unwrap()),
        crate::events::CancellationToken::default(),
    )
    .unwrap();
    let parents = (0..files)
        .map(|file_id| {
            let mut psm = psm_of(&runner, WHOLE_RUN_PARENT, file_id, 1_000_000 + file_id, rt);
            psm.spec_id = format!("injected-parent-{file_id}");
            psm.protein_q = 0.001;
            psm.protein_group_q = 0.001;
            psm.charge = 3;
            psm.expmass = parent_mass;
            psm.calcmass = parent_mass;
            psm
        })
        .collect::<Vec<_>>();
    super::INJECT_PSMS.with(|cell| *cell.borrow_mut() = parents);
    super::FAIL_REREAD.with(|cell| cell.set(fail));
    let result = runner.run(1);
    super::FAIL_REREAD.with(|cell| cell.set(None));
    super::INJECT_PSMS.with(|cell| cell.borrow_mut().clear());
    let events = std::fs::read_to_string(&events_path)
        .unwrap()
        .lines()
        .map(|line| serde_json::from_str(line).unwrap())
        .collect();
    WholeRun {
        result,
        directory: output,
        events,
    }
}

fn base64(bytes: &[u8]) -> String {
    const TABLE: &[u8; 64] = b"ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";
    let mut out = String::new();
    for chunk in bytes.chunks(3) {
        let b = [
            chunk[0],
            *chunk.get(1).unwrap_or(&0),
            *chunk.get(2).unwrap_or(&0),
        ];
        let n = (b[0] as u32) << 16 | (b[1] as u32) << 8 | b[2] as u32;
        for i in 0..4 {
            if i <= chunk.len() {
                out.push(TABLE[(n >> (18 - 6 * i) & 63) as usize] as char);
            } else {
                out.push('=');
            }
        }
    }
    out
}

/// The fixture with nine centroided MS1 scans, 0.1 min apart around `rt`,
/// each carrying the first four isotopes of a precursor of `mass` at
/// `charge` with a Gaussian elution profile.
fn with_ms1_scans(rt: f32, mass: f32, charge: u8) -> String {
    let mzml = std::fs::read_to_string(fixture()).unwrap();
    let ms1 = |index: usize, scan_rt: f32| {
        let mz = (mass as f64 + charge as f64 * 1.007276) / charge as f64;
        let height = 1e8 * (-((scan_rt - rt) as f64 / 0.15).powi(2)).exp();
        let mzs = (0..4)
            .map(|k| mz + k as f64 * 1.003355 / charge as f64)
            .collect::<Vec<_>>();
        let intensities = [1.0f32, 0.9, 0.5, 0.2].map(|f| (height * f as f64) as f32);
        let mz_bytes = mzs.iter().flat_map(|v| v.to_le_bytes()).collect::<Vec<_>>();
        let int_bytes = intensities
            .iter()
            .flat_map(|v| v.to_le_bytes())
            .collect::<Vec<_>>();
        let (mz64, int64) = (base64(&mz_bytes), base64(&int_bytes));
        format!(
            r#"<spectrum index="{index}" id="controllerType=0 controllerNumber=1 scan={scan}" defaultArrayLength="4">
<cvParam cvRef="MS" accession="MS:1000579" name="MS1 spectrum" value=""/>
<cvParam cvRef="MS" accession="MS:1000511" name="ms level" value="1"/>
<cvParam cvRef="MS" accession="MS:1000127" name="centroid spectrum" value=""/>
<scanList count="1"><scan>
<cvParam cvRef="MS" accession="MS:1000016" name="scan start time" value="{scan_rt}" unitCvRef="UO" unitAccession="UO:0000031" unitName="minute"/>
</scan></scanList>
<binaryDataArrayList count="2">
<binaryDataArray encodedLength="{l1}">
<cvParam cvRef="MS" accession="MS:1000523" name="64-bit float" value=""/>
<cvParam cvRef="MS" accession="MS:1000576" name="no compression" value=""/>
<cvParam cvRef="MS" accession="MS:1000514" name="m/z array" value=""/>
<binary>{mz64}</binary>
</binaryDataArray>
<binaryDataArray encodedLength="{l2}">
<cvParam cvRef="MS" accession="MS:1000521" name="32-bit float" value=""/>
<cvParam cvRef="MS" accession="MS:1000576" name="no compression" value=""/>
<cvParam cvRef="MS" accession="MS:1000515" name="intensity array" value=""/>
<binary>{int64}</binary>
</binaryDataArray>
</binaryDataArrayList>
</spectrum>
"#,
            scan = 40000 + index,
            l1 = mz64.len(),
            l2 = int64.len(),
        )
    };
    let start = mzml.find("<spectrum ").unwrap();
    let end = mzml.find("</spectrum>").unwrap() + "</spectrum>".len();
    let offsets = [-0.4f32, -0.3, -0.2, -0.1, 0.05, 0.15, 0.25, 0.35, 0.45];
    let before = offsets
        .iter()
        .filter(|&&d| d < 0.0)
        .enumerate()
        .map(|(i, d)| ms1(i + 100, rt + d))
        .collect::<String>();
    let after = offsets
        .iter()
        .filter(|&&d| d > 0.0)
        .enumerate()
        .map(|(i, d)| ms1(i + 200, rt + d))
        .collect::<String>();
    format!(
        "{}{before}{}\n{after}{}",
        &mzml[..start],
        &mzml[start..end],
        &mzml[end..]
    )
    .replacen(
        r#"<spectrumList count="1""#,
        r#"<spectrumList count="10""#,
        1,
    )
}

/// Each row of a parquet file, as text, with the rescue column dropped.
fn parquet_rows(path: &std::path::Path) -> Vec<String> {
    use parquet::file::reader::{FileReader, SerializedFileReader};
    let reader = SerializedFileReader::new(std::fs::File::open(path).unwrap()).unwrap();
    reader
        .get_row_iter(None)
        .unwrap()
        .map(|row| {
            row.unwrap()
                .get_column_iter()
                .filter(|(name, _)| name.as_str() != "in_source_rescued")
                .map(|(name, value)| format!("{name}={value}"))
                .collect::<Vec<_>>()
                .join("\t")
        })
        .collect()
}

fn parquet_fields(path: &std::path::Path) -> Vec<String> {
    use parquet::file::reader::{FileReader, SerializedFileReader};
    let reader = SerializedFileReader::new(std::fs::File::open(path).unwrap()).unwrap();
    reader
        .metadata()
        .file_metadata()
        .schema()
        .get_fields()
        .iter()
        .map(|field| field.name().to_string())
        .collect()
}

fn event_names(run: &WholeRun) -> Vec<String> {
    run.events
        .iter()
        .map(|event| event["event"].as_str().unwrap().to_string())
        .collect()
}

/// A rescued PSM is written to results.sage.parquet and to nothing else: the
/// PIN file, lfq.parquet and every quality-control table are exactly those of
/// the same run without the pass.
#[test]
fn a_rescued_psm_reaches_the_results_and_no_other_output() {
    let off = whole_run("leak-off", 1, false, None);
    off.result.as_ref().unwrap();
    let on = whole_run("leak-on", 1, true, None);
    on.result.as_ref().unwrap();
    let results = on.directory.join("results.sage.parquet");
    assert!(parquet_fields(&results)
        .iter()
        .any(|f| f == "in_source_rescued"));
    let rows = parquet_rows(&off.directory.join("results.sage.parquet"));
    let rows_on = parquet_rows(&results);
    let rescued = rows_on.len() - rows.len();
    assert!(rescued >= 1, "the run wrote no rescue row");
    assert_eq!(
        rows_on[..rows.len()],
        rows[..],
        "the search's own rows changed"
    );
    assert!(
        rows_on[rows.len()..]
            .iter()
            .all(|row| (row.contains(&format!("\tpeptide=\"{PRODUCT}\""))
                && row.contains(&format!("in_source_parent_peptide=\"{WHOLE_RUN_PARENT}\"")))
                || row.contains("is_decoy=true")),
        "{:?}",
        &rows_on[rows.len()..]
    );
    // No other output changes.
    let mut compared = 0;
    for entry in std::fs::read_dir(&off.directory).unwrap() {
        let name = entry.unwrap().file_name().into_string().unwrap();
        match name.as_str() {
            "results.sage.parquet" | "results.json" | "run-summary.json" => continue,
            "lfq.parquet" => assert_eq!(
                parquet_rows(&off.directory.join(&name)),
                parquet_rows(&on.directory.join(&name)),
                "{name}"
            ),
            _ => assert_eq!(
                std::fs::read(off.directory.join(&name)).unwrap(),
                std::fs::read(on.directory.join(&name)).unwrap(),
                "{name}"
            ),
        }
        compared += 1;
    }
    for required in ["results.sage.pin", "lfq.parquet", "qc.tsv", "digestion.tsv"] {
        assert!(
            on.directory.join(required).exists(),
            "{required} was not written"
        );
    }
    assert!(compared >= 4);
    // LFQ is built before the rescue pass, and this fixture is too small
    // for it to quantify anything, so lfq.parquet is compared here as
    // written; the evidence run compares it row for row on real data.
    // Byte for byte the PIN file of the run without the pass (above). The
    // parent's sequence contains the product's, so a text search would prove
    // nothing.
    assert_eq!(
        event_names(&on).last().map(String::as_str),
        Some("job_completed")
    );
    for run in [off, on] {
        std::fs::remove_dir_all(run.directory.parent().unwrap()).ok();
    }
}

/// A whole run whose rescue fails on its second file, when the first file
/// has a rescue row, writes the main results exactly as a run without the
/// pass, records only the failure and the file in run-summary.json, never
/// reports the job completed, and exits with an error naming the file.
#[test]
fn a_run_whose_rescue_fails_writes_the_main_results_and_exits_with_an_error() {
    let off = whole_run("fail-off", 2, false, None);
    off.result.as_ref().unwrap();
    let completed = whole_run("fail-completed", 2, true, None);
    completed.result.as_ref().unwrap();
    let failed = whole_run("fail-failed", 2, true, Some(1));
    let error = format!("{:#}", failed.result.as_ref().err().unwrap());
    assert!(error.contains("file 1"), "{error}");
    // The completed run has rescue rows, so the failed run could have leaked.
    let rows = parquet_rows(&off.directory.join("results.sage.parquet"));
    assert!(parquet_rows(&completed.directory.join("results.sage.parquet")).len() > rows.len());
    let results = failed.directory.join("results.sage.parquet");
    assert_eq!(parquet_rows(&results), rows);
    assert_eq!(
        parquet_fields(&results),
        parquet_fields(&off.directory.join("results.sage.parquet"))
    );
    let summary: serde_json::Value =
        serde_json::from_slice(&std::fs::read(failed.directory.join("run-summary.json")).unwrap())
            .unwrap();
    let recorded = &summary["in_source"]["rescue_error"];
    assert_eq!(recorded["file_id"], 1);
    assert!(recorded["path"]
        .as_str()
        .unwrap()
        .ends_with("with-ms1.mzML"));
    assert!(
        summary["in_source"].get("rescue").is_none(),
        "partial counts were written"
    );
    // The caller reports the failure; the runner never reports success first.
    let names = event_names(&failed);
    assert!(
        !names.iter().any(|name| name == "job_completed"),
        "{names:?}"
    );
    for run in [off, completed, failed] {
        std::fs::remove_dir_all(run.directory.parent().unwrap()).ok();
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

/// The fixture's peptide with the residues at `i` and `i + 1` exchanged, for
/// each pair given: same composition, so the same precursor mass.
fn swapped(pairs: &[usize]) -> String {
    let mut bytes = PRODUCT.as_bytes().to_vec();
    for &i in pairs {
        bytes.swap(i, i + 1);
    }
    String::from_utf8(bytes).unwrap()
}

/// Second review G: co-elution restricts a spectrum's candidates before they
/// compete. Eleven candidates of parents eluting elsewhere outscore the one
/// co-eluting candidate on the fixture spectrum; with only the top ten hits
/// kept and filtered afterwards, the co-eluting candidate was never seen.
#[test]
fn candidates_of_parents_eluting_elsewhere_cannot_crowd_out_a_coeluting_one() {
    let elsewhere = std::iter::once(PRODUCT.to_string())
        .chain([6, 8, 9, 10, 11, 12, 13, 14, 15, 16].map(|i| swapped(&[i])))
        .collect::<Vec<_>>();
    let coeluting = swapped(&[1, 6, 10, 14]);
    let products = elsewhere
        .iter()
        .chain([&coeluting])
        .cloned()
        .collect::<Vec<_>>();
    let parents = products
        .iter()
        .map(|product| format!("GGWGGK{product}"))
        .collect::<Vec<_>>();
    let (runner, directory) = runner_with(&[fixture()], &parents, &products, serde_json::json!({}));
    let (_, rt) = fixture_spectrum();
    let mut features = parents[..elsewhere.len()]
        .iter()
        .enumerate()
        .map(|(index, parent)| psm_of(&runner, parent, 0, index, rt + 30.0))
        .collect::<Vec<_>>();
    features.push(psm_of(
        &runner,
        parents.last().unwrap(),
        0,
        elsewhere.len(),
        rt,
    ));

    // The setup: on this spectrum the co-eluting candidate ranks below ten.
    let mut stats = InSourceRescueFileStats::default();
    let medians = features
        .iter()
        .map(|feature| (feature.peptide_idx, feature.rt))
        .collect::<fnv::FnvHashMap<_, _>>();
    let (_, peptides) = runner.rescue_candidates(&medians, &mut stats);
    let database = runner.rescue_database(peptides);
    let mut scorer = runner.scorer();
    scorer.db = &database;
    scorer.report_psms = 10;
    scorer.chimera = false;
    let spectra =
        sage_cloudpath::util::read_mzml(&sage_cloudpath::to_url(&fixture()).unwrap(), 0, None)
            .unwrap();
    let (_, processor) = runner.spectrum_reading();
    let spectrum = processor.process(spectra[0].clone());
    let top = scorer
        .score(&spectrum)
        .into_iter()
        .map(|feature| {
            std::str::from_utf8(database[feature.peptide_idx].sequence.as_bytes())
                .unwrap()
                .to_string()
        })
        .collect::<Vec<_>>();
    assert_eq!(top.len(), 10);
    assert!(
        !top.contains(&coeluting),
        "the setup needs it outside the top ten: {top:?}"
    );
    // Third review N1: the restriction comes before every cut, so even a
    // depth of one finds the co-eluting candidate, across the isotope
    // windows (-1..3) and charges this search merges.
    let coeluting_idx = (0..database.peptides.len())
        .map(|i| sage_core::database::PeptideIx(i as u32))
        .find(|&i| database[i].sequence.as_bytes() == coeluting.as_bytes() && !database[i].decoy)
        .unwrap();
    scorer.report_psms = 1;
    let plain = scorer.score(&spectrum);
    assert_ne!(plain[0].peptide_idx, coeluting_idx);
    let keep = |peptide: sage_core::database::PeptideIx| {
        database[peptide].sequence.as_bytes() == coeluting.as_bytes()
    };
    let restricted = scorer.score_restricted(&spectrum, &keep);
    assert_eq!(restricted.len(), 1);
    assert_eq!(restricted[0].peptide_idx, coeluting_idx);
    assert_eq!(restricted[0].delta_best, 0.0);

    let (hits, stats) = runner.rescue_file(
        0,
        &features,
        &settings(),
        &InSourceAnnotations::default(),
        &HashMap::new(),
        None,
    );
    assert_eq!(hits.len(), 1, "{stats:?}");
    let ((_, origin, peptide), _) = &hits[0];
    assert_eq!(peptide.sequence.as_bytes(), coeluting.as_bytes());
    assert_eq!(origin.parent, features.last().unwrap().peptide_idx);
    std::fs::remove_dir_all(directory).ok();
}

/// Second review B: the rescue database has no mass offsets, so no offset
/// hypothesis takes a scorer slot from an unshifted product.
#[test]
fn the_rescue_database_carries_no_mass_offsets() {
    let (runner, directory) = runner_with(
        &[fixture()],
        &[PARENT.to_string()],
        &[PRODUCT.to_string()],
        serde_json::json!({
            "Phospho": {"mass": 79.966331, "sites": ["S", "T", "Y"], "search_mode": "mass_offset"}
        }),
    );
    assert!(
        !runner.database.mass_offsets.is_empty(),
        "the search has an offset"
    );
    let (_, rt) = fixture_spectrum();
    let mut parents = fnv::FnvHashMap::default();
    parents.insert(parent_psm(&runner, 0, 0, rt).peptide_idx, rt);
    let mut stats = InSourceRescueFileStats::default();
    let (_, peptides) = runner.rescue_candidates(&parents, &mut stats);
    let database = runner.rescue_database(peptides);
    assert!(database.mass_offsets.is_empty());
    assert!(database.offset_library.is_none());
    // And the product is still rescued with the offset configured.
    let features = [parent_psm(&runner, 0, 0, rt)];
    let (hits, stats) = runner.rescue_file(
        0,
        &features,
        &settings(),
        &InSourceAnnotations::default(),
        &HashMap::new(),
        None,
    );
    assert_eq!(hits.len(), 1, "{stats:?}");
    assert!(hits[0].0 .0.mass_offset.is_none());
    std::fs::remove_dir_all(directory).ok();
}

/// Second review E: a spectrum is its ID and its occurrence among the file's
/// spectra of that ID. With the fixture's spectrum written twice under one ID,
/// a confident PSM on the first copy leaves the second searched, and the
/// rescued PSM records the second occurrence, so its reporter ions and
/// fragments are looked up on its own spectrum.
#[test]
fn a_repeated_spectrum_id_is_told_apart_by_its_occurrence() {
    let mzml = std::fs::read_to_string(fixture()).unwrap();
    let start = mzml.find("<spectrum ").unwrap();
    let end = mzml.find("</spectrum>").unwrap() + "</spectrum>".len();
    let spectrum = &mzml[start..end];
    let doubled = format!("{}{spectrum}\n{spectrum}{}", &mzml[..start], &mzml[end..]).replacen(
        "<spectrumList count=\"1\"",
        "<spectrumList count=\"2\"",
        1,
    );
    let path = std::env::temp_dir().join(format!(
        "sage-rescue-repeated-{}-{}.mzML",
        std::process::id(),
        std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .unwrap()
            .as_nanos()
    ));
    std::fs::write(&path, doubled).unwrap();
    let (mut runner, directory) = fixture_runner(&[path.to_string_lossy().into_owned()], false);
    let (id, rt) = fixture_spectrum();
    let mut first = parent_psm(&runner, 0, 1, rt);
    first.spec_id = id;
    let features = [parent_psm(&runner, 0, 0, rt), first];
    // The first copy has the confident PSM (occurrence 0, the default).
    let (hits, stats) = runner.rescue_file(
        0,
        &features,
        &settings(),
        &InSourceAnnotations::default(),
        &HashMap::new(),
        None,
    );
    assert_eq!(stats.searched_spectra, 1, "{stats:?}");
    assert_eq!(hits.len(), 1);
    assert_eq!(
        hits[0].1, 1,
        "the rescued spectrum is the second occurrence"
    );
    // The whole pass records it against the rescued PSM's id.
    let outcome = runner
        .in_source_rescue(
            &features,
            &settings(),
            &mut InSourceAnnotations::default(),
            None,
            &HashMap::new(),
        )
        .unwrap();
    assert_eq!(outcome.features.len(), 1);
    assert_eq!(
        outcome.occurrences.get(&outcome.features[0].psm_id),
        Some(&1)
    );
    // With the confident PSM on the second copy instead, the first is searched.
    let occurrences = HashMap::from([(1usize, 1usize)]);
    let (hits, stats) = runner.rescue_file(
        0,
        &features,
        &settings(),
        &InSourceAnnotations::default(),
        &occurrences,
        None,
    );
    assert_eq!(stats.searched_spectra, 1);
    assert_eq!(hits[0].1, 0);
    std::fs::remove_file(path).ok();
    std::fs::remove_dir_all(directory).ok();
}

/// Third review N2: a scan without a precursor charge is not searched, since
/// every charge tried would read a parent's own spectrum at a fraction of its
/// mass; and the charge override is refused at configuration.
#[test]
fn a_scan_without_a_precursor_charge_is_not_searched() {
    let mzml = std::fs::read_to_string(fixture()).unwrap();
    let chargeless = mzml.replacen(
        r#"<cvParam cvRef="MS" accession="MS:1000041" name="charge state" value="3"/>"#,
        "",
        1,
    );
    assert_ne!(chargeless, mzml);
    let path = std::env::temp_dir().join(format!(
        "sage-rescue-chargeless-{}-{}.mzML",
        std::process::id(),
        std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .unwrap()
            .as_nanos()
    ));
    std::fs::write(&path, chargeless).unwrap();
    let (runner, directory) = fixture_runner(&[path.to_string_lossy().into_owned()], false);
    let (_, rt) = fixture_spectrum();
    let features = [parent_psm(&runner, 0, 0, rt)];
    let (hits, stats) = runner.rescue_file(
        0,
        &features,
        &settings(),
        &InSourceAnnotations::default(),
        &HashMap::new(),
        None,
    );
    assert!(hits.is_empty());
    assert_eq!(stats.searched_spectra, 0);
    assert_eq!(stats.skipped_without_charge, 1);
    std::fs::remove_file(path).ok();
    std::fs::remove_dir_all(directory).ok();

    let workspace = workspace();
    let input: crate::input::Input = serde_json::from_value(serde_json::json!({
        "database": {
            "fasta": format!("{workspace}/tests/Q99536.fasta"),
            "enzyme": {"cleave_at": "KR", "restrict": "P", "missed_cleavages": 0, "semi_enzymatic": false},
            "static_mods": {},
            "ion_kinds": ["b", "y"]
        },
        "precursor_tol": { "ppm": [-10, 10] },
        "fragment_tol": { "ppm": [-10, 10] },
        "mzml_paths": [fixture()],
        "override_precursor_charge": true,
        "in_source": {"rescue": true},
    }))
    .unwrap();
    let error = format!("{:#}", input.build().err().unwrap());
    assert!(error.contains("override_precursor_charge"), "{error}");
}

/// Third review N6: a product that several parents give lies in all of
/// their proteins, within a file and across files.
#[test]
fn a_product_of_several_parents_lies_in_all_of_their_proteins() {
    let (runner, directory) = runner_with(
        &[fixture()],
        &[PARENT.to_string(), format!("MMWWK{PRODUCT}")],
        &[PRODUCT.to_string()],
        serde_json::json!({}),
    );
    let (_, rt) = fixture_spectrum();
    let mut parents = fnv::FnvHashMap::default();
    parents.insert(psm_of(&runner, PARENT, 0, 0, rt).peptide_idx, rt);
    parents.insert(
        psm_of(&runner, &format!("MMWWK{PRODUCT}"), 0, 1, rt).peptide_idx,
        rt,
    );
    let mut stats = InSourceRescueFileStats::default();
    let (origins, peptides) = runner.rescue_candidates(&parents, &mut stats);
    let product = peptides
        .iter()
        .find(|p| !p.decoy && p.sequence.as_bytes() == PRODUCT.as_bytes())
        .unwrap();
    assert_eq!(product.proteins.len(), 2, "{:?}", product.proteins);
    assert_eq!(origins[&(PRODUCT.to_string(), false)].len(), 2);
    std::fs::remove_dir_all(directory).ok();
}
