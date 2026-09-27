use super::*;
use crate::spectrum::Precursor;

#[test]
fn predefined_tags_expose_expected_channels_and_masses() {
    for (tag, channels, modification_mass) in [
        (Isobaric::Tmt6, 6, 229.162932),
        (Isobaric::Tmt10, 10, 229.162932),
        (Isobaric::Tmt11, 11, 229.162932),
        (Isobaric::Tmt16, 16, 304.2071),
        (Isobaric::Tmt18, 18, 304.2071),
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
    assert_eq!(quant[0].peaks[0], Some(42.0));
    assert_eq!(quant[0].peaks.len(), 6);
}

#[test]
fn absent_reporter_channels_are_missing_not_zero() {
    let spectra = vec![spectrum(2, "ms2", Some("ms1"), 42.0)];

    let quant = quantify(&spectra, &Isobaric::Tmt6, Tolerance::Da(-0.01, 0.01), 2);

    // Only the 126 channel has a peak; the other five are missing.
    assert_eq!(
        quant[0].peaks,
        vec![Some(42.0), None, None, None, None, None]
    );
}

#[test]
fn zero_and_non_finite_reporter_peaks_are_not_observations() {
    let labels = &Isobaric::Tmt6.reporter_masses()[..3];
    let at = |idx: usize, offset: f32| labels[idx] - PROTON + offset;
    // Channel 0: only a zero-intensity centroid. Channel 1: an infinite S/N
    // value (zero noise) next to a real peak. Channel 2: NaN (0 / 0) only.
    let masses = vec![at(0, 0.0), at(1, -0.002), at(1, 0.002), at(2, 0.0)];
    let intensities = vec![0.0, f32::INFINITY, 7.5, f32::NAN];

    let found = find_reporter_ions(&masses, &intensities, labels, Tolerance::Da(-0.01, 0.01));

    assert_eq!(found, vec![None, Some(7.5), None]);
}

#[test]
fn channel_summary_skips_missing_channels() {
    let quant = |peaks: Vec<Option<f32>>| TmtQuant {
        spec_id: String::new(),
        file_id: 0,
        occurrence: 0,
        ion_injection_time: 0.0,
        peaks,
    };
    let quant = vec![
        quant(vec![Some(10.0), None, None]),
        quant(vec![Some(30.0), Some(4.0), None]),
        quant(vec![None, Some(8.0), None]),
    ];
    let headers = Isobaric::User(vec![1.0, 2.0, 3.0]).headers();

    let summary = summarize_channels(&quant, &headers);

    let expected = |channel: &str, observed, missing, median| ReporterChannelSummary {
        channel: channel.into(),
        observed,
        missing,
        median_intensity: median,
    };
    assert_eq!(
        summary,
        vec![
            // A zero fill would give a median of 10 for both observed channels.
            expected("user_1", 2, 1, Some(20.0)),
            expected("user_2", 2, 1, Some(6.0)),
            expected("user_3", 0, 3, None),
        ]
    );
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
