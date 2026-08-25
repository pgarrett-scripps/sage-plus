use super::*;
use crate::spectrum::Precursor;

#[test]
fn predefined_tags_expose_expected_channels_and_masses() {
    for (tag, channels, modification_mass) in [
        (Isobaric::Tmt6, 6, 229.162932),
        (Isobaric::Tmt10, 10, 229.162932),
        (Isobaric::Tmt11, 11, 229.162932),
        (Isobaric::Tmt16, 16, 304.2071),
        (Isobaric::Tmt18, 18, 304.2135),
    ] {
        assert_eq!(tag.reporter_masses().len(), channels);
        assert_eq!(tag.headers().len(), channels);
        assert_eq!(tag.modification_mass(), Some(modification_mass));
    }
}

#[test]
fn user_tags_have_stable_headers_and_no_modification_mass() {
    let tag = Isobaric::User(vec![100.0, 110.0, 120.0]);

    assert_eq!(tag.reporter_masses(), &[100.0, 110.0, 120.0]);
    assert_eq!(tag.headers(), vec!["user_1", "user_2", "user_3"]);
    assert_eq!(tag.modification_mass(), None);
}

#[test]
fn reporter_search_selects_the_most_intense_peak_and_marks_missing_channels() {
    let label = Isobaric::Tmt6.reporter_masses()[0];
    let masses = vec![label - PROTON - 0.005, label - PROTON + 0.002, 200.0];
    let intensities = vec![5.0, 25.0, 1000.0];

    let found = find_reporter_ions(
        &masses,
        &intensities,
        &[label, 150.0],
        Tolerance::Da(-0.01, 0.01),
    );

    assert_eq!(found, vec![Some(25.0), None]);
}

fn spectrum(level: u8, id: &str, parent: Option<&str>, intensity: f32) -> ProcessedSpectrum {
    let label = Isobaric::Tmt6.reporter_masses()[0];
    ProcessedSpectrum {
        level,
        id: id.into(),
        file_id: level as usize,
        ion_injection_time: 12.5,
        precursors: parent
            .map(|spectrum_ref| Precursor {
                spectrum_ref: Some(spectrum_ref.into()),
                ..Default::default()
            })
            .into_iter()
            .collect(),
        masses: vec![label - PROTON],
        intensities: vec![intensity],
        ..Default::default()
    }
}

#[test]
fn ms2_quantification_uses_the_spectrum_identifier() {
    let spectra = vec![
        spectrum(1, "ms1", None, 1.0),
        spectrum(2, "ms2", Some("ms1"), 42.0),
        spectrum(3, "ms3", Some("ms2"), 84.0),
    ];

    let quant = quantify(&spectra, &Isobaric::Tmt6, Tolerance::Da(-0.01, 0.01), 2);

    assert_eq!(quant.len(), 1);
    assert_eq!(quant[0].spec_id, "ms2");
    assert_eq!(quant[0].file_id, 2);
    assert_eq!(quant[0].peaks[0], 42.0);
    assert_eq!(quant[0].peaks.len(), 6);
}

#[test]
fn msn_quantification_uses_the_parent_spectrum_identifier() {
    let spectra = vec![
        spectrum(3, "ms3-with-parent", Some("ms2-parent"), 84.0),
        spectrum(3, "ms3-without-parent", None, 21.0),
    ];

    let quant = quantify(&spectra, &Isobaric::Tmt6, Tolerance::Da(-0.01, 0.01), 3);

    assert_eq!(quant.len(), 2);
    assert_eq!(quant[0].spec_id, "ms2-parent");
    assert_eq!(quant[1].spec_id, "");
}

#[test]
fn ms1_quantification_is_explicitly_disabled() {
    let spectra = vec![spectrum(1, "ms1", None, 99.0)];

    assert!(quantify(&spectra, &Isobaric::Tmt6, Tolerance::Da(-0.01, 0.01), 1).is_empty());
}
