use super::*;
use crate::database::PeptideIx;
use crate::scoring::Fragments;
use std::sync::Arc;

fn feature(peptide_idx: u32, psm_id: usize, spectrum_q: f32, score: f32) -> Feature {
    Feature {
        peptide_idx: PeptideIx(peptide_idx),
        psm_id,
        label: 1,
        rank: 1,
        charge: 2,
        matched_peaks: 8,
        spectrum_q,
        peptide_q: 0.005,
        discriminant_score: score,
        calcmass: 1_000.0,
        rt: 12.0,
        aligned_rt: 11.5,
        fragments: Some(Fragments {
            kinds: vec![Kind::B, Kind::Y, Kind::Y],
            fragment_ordinals: vec![2, 3, 4],
            charges: vec![1, 1, 2],
            neutral_losses: vec![0.0, 18.010_565, 42.0],
            mz_calculated: vec![200.0, 300.0, 400.0],
            mz_experimental: vec![200.01, 300.01, 400.01],
            intensities: vec![25.0, 100.0, 0.5],
            peak_indices: vec![0, 1, 2],
        }),
        ..Feature::default()
    }
}

fn database() -> IndexedDatabase {
    let mut database = IndexedDatabase::default();
    database.peptides.push(Peptide {
        sequence: (&b"PEPTIDE"[..]).into(),
        proteins: vec![Arc::from("P12345")].into(),
        ..Peptide::default()
    });
    database
}

#[test]
fn best_psm_is_deterministic_and_counts_support() {
    let database = database();
    let settings = SpectralLibrarySettings {
        enabled: true,
        ..Default::default()
    };
    let features = vec![
        feature(0, 3, 0.005, 10.0),
        feature(0, 2, 0.001, 5.0),
        feature(0, 1, 0.001, 8.0),
    ];
    let selected = select_best_psms(&features, &database, &settings);
    assert_eq!(selected.len(), 1);
    assert_eq!(selected[0].feature_index, 2);
    assert_eq!(selected[0].supporting_psms, 3);
    assert_eq!(selected[0].feature_indices, vec![2]);
}

#[test]
fn builds_normalized_top_fragments_and_mzspeclib() {
    let database = database();
    let settings = SpectralLibrarySettings {
        enabled: true,
        max_fragments: 2,
        ..Default::default()
    };
    let features = vec![feature(0, 1, 0.001, 8.0)];
    let selected = select_best_psms(&features, &database, &settings);
    let mut entries = build_entries(
        &features,
        &database,
        &["sample.mzML".into()],
        &selected,
        &settings,
    )
    .unwrap();
    assert_eq!(entries.len(), 1);
    assert_eq!(entries[0].fragments.len(), 2);
    assert_eq!(entries[0].fragments[0].relative_intensity, 0.25);
    assert_eq!(entries[0].fragments[1].relative_intensity, 1.0);
    assert_eq!(entries[0].fragments[1].annotation(), "y3-H2O");
    entries[0].label_channel = Some("heavy".into());
    entries[0].label_group = Some("PEPTIDE".into());
    entries[0].label_reference = Some("light".into());

    let text = String::from_utf8(serialize_mzspeclib(
        &entries,
        "0.1.0-beta.1",
        SpectralLibraryStrategy::BestPsm,
    ))
    .unwrap();
    assert!(text.starts_with("<mzSpecLib>\n"));
    assert!(text.contains("<Spectrum=1>"));
    assert!(text.contains("MS:1003270|proforma peptidoform ion notation=PEPTIDE/2"));
    assert!(text.contains("300.000000\t10000.000000\ty3-H2O"));
    assert!(text.contains("SAGE:1000001|label channel=heavy"));
    assert!(text.contains("SAGE:1000002|label group=PEPTIDE"));
    assert!(text.contains("SAGE:1000003|label reference=light"));
}

#[test]
fn consensus_uses_median_properties_and_reproducible_fragments() {
    let database = database();
    let settings = SpectralLibrarySettings {
        enabled: true,
        strategy: SpectralLibraryStrategy::Consensus,
        min_fragment_frequency: 0.66,
        ..Default::default()
    };
    let mut features = vec![
        feature(0, 1, 0.001, 8.0),
        feature(0, 2, 0.002, 7.0),
        feature(0, 3, 0.003, 6.0),
    ];
    features[0].rt = 10.0;
    features[0].aligned_rt = 9.0;
    features[0].ims = 1.0;
    features[0].fragments.as_mut().unwrap().intensities = vec![25.0, 100.0, 0.0];
    features[1].rt = 12.0;
    features[1].aligned_rt = 11.0;
    features[1].ims = 1.2;
    features[1].fragments.as_mut().unwrap().intensities = vec![50.0, 100.0, 0.0];
    features[2].rt = 30.0;
    features[2].aligned_rt = 25.0;
    features[2].ims = 1.4;
    features[2].fragments.as_mut().unwrap().intensities = vec![75.0, 0.0, 100.0];

    let selected = select_psms(&features, &database, &settings);
    assert_eq!(selected[0].feature_indices, vec![0, 1, 2]);
    let entries = build_entries(
        &features,
        &database,
        &["one.mzML".into()],
        &selected,
        &settings,
    )
    .unwrap();
    assert_eq!(entries[0].supporting_psms, 3);
    assert_eq!(entries[0].retention_time_minutes, 12.0);
    assert_eq!(entries[0].aligned_retention_time_minutes, 11.0);
    assert_eq!(entries[0].ion_mobility, 1.2);
    assert_eq!(entries[0].fragments.len(), 2);
    assert_eq!(entries[0].fragments[0].relative_intensity, 0.5);
    assert_eq!(entries[0].fragments[1].relative_intensity, 1.0);

    let text = String::from_utf8(serialize_mzspeclib(
        &entries,
        "0.1.0-beta.1",
        SpectralLibraryStrategy::Consensus,
    ))
    .unwrap();
    assert!(text.contains("MS:1003067|consensus spectrum"));
    assert!(text.contains("strategy=consensus"));
}

#[test]
fn settings_reject_invalid_values_and_duplicate_formats() {
    let mut settings = SpectralLibrarySettings {
        enabled: true,
        ..Default::default()
    };
    settings.max_fragments = 0;
    assert!(settings.validate().is_err());
    settings.max_fragments = 10;
    settings.min_fragment_frequency = 0.0;
    assert!(settings.validate().is_err());
    settings.min_fragment_frequency = 0.5;
    settings.formats.push(SpectralLibraryFormat::MzSpecLib);
    assert!(settings.validate().is_err());
}

#[test]
fn mass_delta_proforma_does_not_depend_on_modification_names() {
    let peptide = Peptide {
        sequence: (&b"ACD"[..]).into(),
        modifications: crate::peptide::CompactModifications::from_dense([0.0, 57.021_465, 0.0]),
        nterm: Some(42.010_565),
        ..Peptide::default()
    };
    assert_eq!(
        mass_delta_proforma(&peptide),
        "[+42.010567]-AC[+57.021465]D"
    );
}

#[test]
fn mass_delta_proforma_renders_negative_and_c_terminal_deltas() {
    let peptide = Peptide {
        sequence: (&b"MK"[..]).into(),
        modifications: crate::peptide::CompactModifications::from_dense([-1.5, 0.0]),
        cterm: Some(-0.984_016),
        ..Peptide::default()
    };
    assert_eq!(mass_delta_proforma(&peptide), "M[-1.5]K-[-0.984016]");

    let plain = Peptide {
        sequence: (&b"PEPTIDE"[..]).into(),
        ..Peptide::default()
    };
    assert_eq!(mass_delta_proforma(&plain), "PEPTIDE");
}

#[test]
fn every_invalid_setting_has_its_own_message() {
    let base = SpectralLibrarySettings {
        enabled: true,
        ..Default::default()
    };
    assert_eq!(base.validate(), Ok(()));
    let cases: Vec<(SpectralLibrarySettings, &str)> = vec![
        (
            SpectralLibrarySettings {
                psm_q_value: 1.5,
                ..base.clone()
            },
            "psm_q_value",
        ),
        (
            SpectralLibrarySettings {
                psm_q_value: f32::NAN,
                ..base.clone()
            },
            "psm_q_value",
        ),
        (
            SpectralLibrarySettings {
                peptide_q_value: -0.1,
                ..base.clone()
            },
            "peptide_q_value",
        ),
        (
            SpectralLibrarySettings {
                min_matched_peaks: 0,
                ..base.clone()
            },
            "min_matched_peaks",
        ),
        (
            SpectralLibrarySettings {
                min_relative_intensity: 2.0,
                ..base.clone()
            },
            "min_relative_intensity",
        ),
        (
            SpectralLibrarySettings {
                min_consensus_psms: 0,
                ..base.clone()
            },
            "min_consensus_psms",
        ),
        (
            SpectralLibrarySettings {
                min_fragment_frequency: 1.01,
                ..base.clone()
            },
            "min_fragment_frequency",
        ),
        (
            SpectralLibrarySettings {
                formats: Vec::new(),
                ..base.clone()
            },
            "formats must not be empty",
        ),
    ];
    for (settings, field) in cases {
        let err = settings.validate().unwrap_err();
        assert!(err.contains(field), "expected {field}: {err}");
    }
    // Empty formats are fine while the library is disabled.
    let disabled = SpectralLibrarySettings {
        enabled: false,
        formats: Vec::new(),
        ..Default::default()
    };
    assert_eq!(disabled.validate(), Ok(()));
}

#[test]
fn writes_requires_enabled_and_listed_format() {
    let mut settings = SpectralLibrarySettings {
        enabled: true,
        formats: vec![SpectralLibraryFormat::MzSpecLib],
        ..Default::default()
    };
    assert!(settings.writes(SpectralLibraryFormat::MzSpecLib));
    assert!(!settings.writes(SpectralLibraryFormat::SageParquet));
    settings.enabled = false;
    assert!(!settings.writes(SpectralLibraryFormat::MzSpecLib));
}

#[test]
fn fragment_annotations_follow_mzpaf() {
    let fragment = |kind: Kind, ordinal: i32, charge: i32, neutral_loss: f32| LibraryFragment {
        kind,
        ordinal,
        charge,
        neutral_loss,
        mz: 0.0,
        relative_intensity: 1.0,
    };
    assert_eq!(fragment(Kind::B, 3, 1, 0.0).annotation(), "b3");
    assert_eq!(fragment(Kind::Y, 7, 2, 0.0).annotation(), "y7^2");
    assert_eq!(fragment(Kind::A, 2, 1, 0.0).annotation(), "a2");
    assert_eq!(fragment(Kind::C, 4, 1, 0.0).annotation(), "c4");
    assert_eq!(fragment(Kind::X, 5, 1, 0.0).annotation(), "x5");
    assert_eq!(fragment(Kind::Z, 6, 1, 0.0).annotation(), "z6");
    assert_eq!(fragment(Kind::ZDot, 6, 3, 0.0).annotation(), "z.6^3");
    assert_eq!(fragment(Kind::Y, 5, 1, 18.010_565).annotation(), "y5-H2O");
    assert_eq!(fragment(Kind::B, 5, 2, 17.026_55).annotation(), "b5-NH3^2");
    assert_eq!(fragment(Kind::Y, 9, 1, 97.976_9).annotation(), "y9-H3PO4");
    // Unknown losses are not invented.
    assert_eq!(fragment(Kind::Y, 9, 1, 44.0).annotation(), "?");
}

#[test]
fn median_ignores_non_finite_values() {
    assert_eq!(median(vec![3.0, 1.0, 2.0]), Some(2.0));
    assert_eq!(median(vec![4.0, 1.0, 3.0, 2.0]), Some(2.5));
    assert_eq!(median(vec![f32::NAN, 5.0, f32::INFINITY]), Some(5.0));
    assert_eq!(median(vec![f32::NAN]), None);
    assert_eq!(median(Vec::new()), None);
}
