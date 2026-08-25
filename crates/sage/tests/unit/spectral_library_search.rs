use super::*;
use crate::ion_series::Kind;
use crate::spectrum::Precursor;

fn fragment(mz: f32, relative_intensity: f32) -> LibraryFragment {
    LibraryFragment {
        kind: Kind::Y,
        ordinal: 1,
        charge: 1,
        neutral_loss: 0.0,
        mz,
        relative_intensity,
    }
}

fn entry(id: &str, charge: u8, fragments: Vec<LibraryFragment>) -> DdaLibraryEntry {
    DdaLibraryEntry {
        library_entry_id: id.into(),
        source_file: String::new(),
        source_spectrum: String::new(),
        proforma: id.into(),
        stripped_peptide: id.into(),
        proteins: "P1".into(),
        label_channel: None,
        label_group: None,
        label_reference: None,
        precursor_charge: charge,
        precursor_neutral_mass: 1_000.0,
        precursor_mz: 1_000.0 / charge as f32 + PROTON,
        retention_time_minutes: 10.0,
        ion_mobility: 1.1,
        source_spectrum_q: 0.001,
        is_decoy: false,
        fragments,
    }
}

fn query(peaks: &[(f32, f32)]) -> ProcessedSpectrum {
    ProcessedSpectrum {
        level: 2,
        scan_start_time: 10.5,
        precursors: vec![Precursor {
            inverse_ion_mobility: Some(1.2),
            ..Default::default()
        }],
        masses: peaks.iter().map(|(mz, _)| mz - PROTON).collect(),
        intensities: peaks.iter().map(|(_, intensity)| *intensity).collect(),
        charges: vec![1; peaks.len()],
        ..Default::default()
    }
}

#[test]
fn matching_intensity_pattern_beats_mass_matched_distractor() {
    let index = DdaLibraryIndex::new(vec![
        entry(
            "matching",
            2,
            vec![fragment(200.0, 1.0), fragment(300.0, 0.25)],
        ),
        entry(
            "distractor",
            2,
            vec![fragment(200.0, 0.1), fragment(300.0, 1.0)],
        ),
    ])
    .unwrap();
    let query = query(&[(200.0, 100.0), (300.0, 25.0)]);
    let matches = index.search(
        &query,
        1_000.0,
        2,
        DdaLibrarySearchParameters {
            min_matched_peaks: 2,
            max_hits: 2,
            annotate_matches: true,
            ..Default::default()
        },
    );
    assert_eq!(matches.len(), 2);
    assert_eq!(
        index.entries()[matches[0].entry_index].library_entry_id,
        "matching"
    );
    assert!((matches[0].spectral_angle - 1.0).abs() < 1e-5);
    assert!(matches[0].spectral_angle > matches[1].spectral_angle);
    assert_eq!(matches[0].ion_mobility_delta, Some(0.100_000_024));
    assert_eq!(matches[0].fragments.as_ref().unwrap().kinds.len(), 2);
    assert!(matches[0].average_fragment_ppm.abs() < 1e-4);
    assert!(matches[0].signed_fragment_ppm.abs() < 1e-4);
}

#[test]
fn precursor_charge_filters_candidates() {
    let index = DdaLibraryIndex::new(vec![
        entry("charge-2", 2, vec![fragment(200.0, 1.0)]),
        entry("charge-3", 3, vec![fragment(200.0, 1.0)]),
    ])
    .unwrap();
    let query = query(&[(200.0, 100.0)]);
    let matches = index.search(
        &query,
        1_000.0,
        3,
        DdaLibrarySearchParameters {
            precursor_tolerance: Tolerance::Da(-0.1, 0.1),
            min_matched_peaks: 1,
            max_hits: 10,
            ..Default::default()
        },
    );
    assert_eq!(matches.len(), 1);
    assert_eq!(
        index.entries()[matches[0].entry_index].library_entry_id,
        "charge-3"
    );
}

#[test]
fn one_query_peak_cannot_match_two_library_peaks() {
    let index = DdaLibraryIndex::new(vec![entry(
        "overlap",
        2,
        vec![fragment(200.000, 1.0), fragment(200.001, 0.5)],
    )])
    .unwrap();
    let query = query(&[(200.0005, 100.0)]);
    let matches = index.search(
        &query,
        1_000.0,
        2,
        DdaLibrarySearchParameters {
            precursor_tolerance: Tolerance::Da(-0.1, 0.1),
            fragment_tolerance: Tolerance::Da(-0.01, 0.01),
            min_matched_peaks: 1,
            max_hits: 1,
            ..Default::default()
        },
    );
    assert_eq!(matches[0].matched_peaks, 1);
}

#[test]
fn isotope_error_corrects_precursor_mass_and_is_reported() {
    let index = DdaLibraryIndex::new(vec![entry(
        "isotope-plus-one",
        2,
        vec![fragment(200.0, 1.0)],
    )])
    .unwrap();
    let query = query(&[(200.0, 100.0)]);
    let observed_mass = 1_000.0 + NEUTRON;

    let without_isotope = index.search(
        &query,
        observed_mass,
        2,
        DdaLibrarySearchParameters {
            precursor_tolerance: Tolerance::Da(-0.01, 0.01),
            min_matched_peaks: 1,
            ..Default::default()
        },
    );
    assert!(without_isotope.is_empty());

    let with_isotope = index.search(
        &query,
        observed_mass,
        2,
        DdaLibrarySearchParameters {
            precursor_tolerance: Tolerance::Da(-0.01, 0.01),
            min_matched_peaks: 1,
            min_isotope_error: 0,
            max_isotope_error: 1,
            ..Default::default()
        },
    );
    assert_eq!(with_isotope.len(), 1);
    assert_eq!(with_isotope[0].isotope_error, 1);
    assert!(with_isotope[0].precursor_ppm.abs() < 1e-4);
}

#[test]
fn overlapping_isotope_windows_do_not_duplicate_entries() {
    let index =
        DdaLibraryIndex::new(vec![entry("wide-window", 2, vec![fragment(200.0, 1.0)])]).unwrap();
    let matches = index.search(
        &query(&[(200.0, 100.0)]),
        1_000.0 + NEUTRON,
        2,
        DdaLibrarySearchParameters {
            precursor_tolerance: Tolerance::Da(-2.0, 2.0),
            min_matched_peaks: 1,
            max_hits: 10,
            min_isotope_error: 0,
            max_isotope_error: 1,
            ..Default::default()
        },
    );
    assert_eq!(matches.len(), 1);
    assert_eq!(matches[0].isotope_error, 1);
}

#[test]
fn systematic_mass_offsets_are_corrected_during_search() {
    let index =
        DdaLibraryIndex::new(vec![entry("calibrated", 2, vec![fragment(200.0, 1.0)])]).unwrap();
    let ppm = 5.0;
    let shifted_precursor = 1_000.0 * (1.0 + ppm * 1e-6);
    let shifted_fragment = 200.0 * (1.0 + ppm * 1e-6);
    let query = query(&[(shifted_fragment, 100.0)]);

    let matches = index.search(
        &query,
        shifted_precursor,
        2,
        DdaLibrarySearchParameters {
            precursor_tolerance: Tolerance::Ppm(-1.0, 1.0),
            fragment_tolerance: Tolerance::Ppm(-1.0, 1.0),
            min_matched_peaks: 1,
            precursor_offset_ppm: ppm,
            fragment_offset_ppm: ppm,
            ..Default::default()
        },
    );
    assert_eq!(matches.len(), 1);
    assert!((matches[0].raw_precursor_ppm - ppm).abs() < 0.1);
    assert!(matches[0].precursor_ppm.abs() < 0.1);
    assert!((matches[0].signed_fragment_ppm - ppm).abs() < 0.1);
    assert!(matches[0].aligned_average_fragment_ppm < 0.1);
}

#[test]
fn invalid_entries_are_rejected() {
    let invalid = entry("empty", 2, Vec::new());
    assert!(DdaLibraryIndex::new(vec![invalid])
        .unwrap_err()
        .contains("no fragment peaks"));
}

#[test]
fn parses_mass_and_unimod_proforma() {
    let parsed = parse_proforma("[+42.010565]-AC[UNIMOD:4]DM[Oxidation]-[+1.0]/2").unwrap();
    assert_eq!(parsed.sequence, b"ACDM");
    assert!((parsed.nterm.unwrap() - 42.010_565).abs() < 1e-5);
    assert!((parsed.modifications[1] - 57.021_465).abs() < 1e-4);
    assert!((parsed.modifications[3] - 15.994_915).abs() < 1e-4);
    assert_eq!(parsed.cterm, Some(1.0));
}

#[test]
fn shuffled_decoys_preserve_precursor_and_intensities() {
    let target = DdaLibraryEntry {
        library_entry_id: "PEPTIDER/2".into(),
        source_file: "sample.mzML".into(),
        source_spectrum: "scan=42".into(),
        proforma: "PEPTIDER".into(),
        stripped_peptide: "PEPTIDER".into(),
        proteins: "P1;P2".into(),
        precursor_charge: 2,
        precursor_neutral_mass: 955.461,
        precursor_mz: 478.7378,
        fragments: vec![
            LibraryFragment {
                kind: Kind::B,
                ordinal: 3,
                charge: 1,
                neutral_loss: 0.0,
                mz: 324.155,
                relative_intensity: 1.0,
            },
            LibraryFragment {
                kind: Kind::Y,
                ordinal: 4,
                charge: 2,
                neutral_loss: 0.0,
                mz: 250.0,
                relative_intensity: 0.4,
            },
        ],
        ..DdaLibraryEntry::default()
    };
    let entries = generate_decoys(
        vec![target.clone()],
        &LibrarySearchSettings {
            path: "library.parquet".into(),
            ..LibrarySearchSettings::default()
        },
    )
    .unwrap();
    let decoy = entries.iter().find(|entry| entry.is_decoy).unwrap();
    assert_ne!(decoy.proforma, target.proforma);
    assert_eq!(decoy.precursor_neutral_mass, target.precursor_neutral_mass);
    assert_eq!(decoy.proteins, target.proteins);
    assert_eq!(decoy.source_file, target.source_file);
    assert_eq!(decoy.source_spectrum, target.source_spectrum);
    let mut intensities = decoy
        .fragments
        .iter()
        .map(|fragment| fragment.relative_intensity)
        .collect::<Vec<_>>();
    intensities.sort_unstable_by(f32::total_cmp);
    assert_eq!(intensities, vec![0.4, 1.0]);
    assert!(decoy.fragments.iter().any(|decoy| {
        target.fragments.iter().any(|target| {
            decoy.kind == target.kind && decoy.ordinal == target.ordinal && decoy.mz != target.mz
        })
    }));
}

#[test]
fn parses_mzspeclib_with_protein_mappings() {
    let text = r#"<mzSpecLib>

<Spectrum=1>
MS:1003061|library spectrum name=PEPTIDER/2
MS:1003208|experimental precursor monoisotopic m/z=478.7378

<Analyte=1>
MS:1003270|proforma peptidoform ion notation=PEPTIDER/2
MS:1001117|theoretical mass=955.461
[1]MS:1000885|protein accession=P1
[2]MS:1000885|protein accession=P2

<Peaks>
300.0 10000.0 y3
400.0 5000.0 b4^2
"#;
    let entries = deserialize_mzspeclib(text).unwrap();
    assert_eq!(entries.len(), 1);
    assert_eq!(entries[0].proteins, "P1;P2");
    assert_eq!(entries[0].fragments.len(), 2);
    assert_eq!(entries[0].fragments[0].relative_intensity, 1.0);
    assert_eq!(entries[0].fragments[1].charge, 2);
}
