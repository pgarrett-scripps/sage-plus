use super::*;

fn buffer(peaks: Vec<ImsPeak>) -> PeakBuffer {
    let mut order = (0..peaks.len()).collect::<Vec<_>>();
    order.sort_unstable_by(|left, right| {
        peaks[*right]
            .intensity
            .partial_cmp(&peaks[*left].intensity)
            .unwrap()
    });
    PeakBuffer {
        peaks,
        order,
        agg_buff: Vec::new(),
    }
}

#[test]
fn parses_precursor_metadata_without_losing_optional_values() {
    use timsrust::core::{Charge, FrameIndex, Im, Mz, Rt, ScanIndex};

    let precursor = timsrust::core::Precursor::new(
        Mz::from(500.25),
        Im::from(1.15),
        Rt::from(12.0),
        ScanIndex::try_from(5).unwrap(),
        Some(Charge::try_from(3).unwrap()),
        Some(1234.5),
        7,
        FrameIndex::try_from(42).unwrap(),
    );
    let parsed = TdfReader::parse_precursor(&precursor);

    assert_eq!(parsed.mz, 500.25);
    assert_eq!(parsed.charge, Some(3));
    assert_eq!(parsed.intensity, Some(1234.5));
    assert_eq!(parsed.spectrum_ref.as_deref(), Some("42"));
    assert_eq!(parsed.inverse_ion_mobility, Some(1.15));
}

#[test]
fn parses_real_bruker_directory() {
    let path = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/bruker/example_dia.d");
    let url = crate::Url::from_file_path(path).unwrap();
    let spectra =
        crate::util::read_spectra(&url, 3, None, BrukerProcessingConfig::default(), true).unwrap();

    assert!(!spectra.is_empty());
    assert!(spectra.iter().all(|spectrum| {
        spectrum.file_id == 3 && spectrum.mz.len() == spectrum.intensity.len()
    }));
    assert!(spectra.iter().any(|spectrum| spectrum.ms_level == 1));
    assert!(spectra.iter().any(|spectrum| spectrum.ms_level > 1));
    assert!(spectra
        .iter()
        .filter(|spectrum| spectrum.ms_level == 1)
        .all(|spectrum| spectrum
            .mobility
            .as_ref()
            .is_some_and(|values| { !values.is_empty() && values.len() == spectrum.mz.len() })));
    assert!(spectra
        .iter()
        .filter(|spectrum| spectrum.ms_level > 1)
        .all(|spectrum| !spectrum.precursors.is_empty()));
}

#[test]
fn mobility_offsets_expand_run_lengths_and_skip_empty_scans() {
    struct LinearImConverter;

    impl Converter<ScanIndex, Im> for LinearImConverter {
        fn convert(&self, scan_index: ScanIndex) -> Im {
            Im::from(2.0 - f64::from(scan_index) * 0.25)
        }
    }

    let converter = LinearImConverter;
    let mobility = PeakBuffer::expand_mobility_iter(&[0, 2, 2, 5], &converter).collect::<Vec<_>>();

    assert_eq!(mobility, vec![2.0, 2.0, 1.5, 1.5, 1.5]);
}

#[test]
fn centroiding_combines_nearby_mass_and_mobility_peaks() {
    let mut buffer = buffer(vec![
        ImsPeak {
            mz: 100.0,
            intensity: 10.0,
            im: 1.0,
        },
        ImsPeak {
            mz: 100.0004,
            intensity: 5.0,
            im: 1.01,
        },
        ImsPeak {
            mz: 101.0,
            intensity: 7.0,
            im: 1.2,
        },
    ]);

    let (mz, (intensity, mobility)) = buffer.fastcentroid_frame(5.0, 2.0);

    assert_eq!(mz, vec![100.0, 101.0]);
    assert_eq!(intensity, vec![15.0, 7.0]);
    assert_eq!(mobility, vec![1.0, 1.2]);
}

#[test]
fn centroiding_keeps_close_masses_separate_when_mobility_differs() {
    let mut buffer = buffer(vec![
        ImsPeak {
            mz: 100.0,
            intensity: 10.0,
            im: 1.0,
        },
        ImsPeak {
            mz: 100.0004,
            intensity: 5.0,
            im: 1.2,
        },
    ]);

    let (mz, (intensity, mobility)) = buffer.fastcentroid_frame(5.0, 2.0);

    assert_eq!(mz, vec![100.0, 100.0004]);
    assert_eq!(intensity, vec![10.0, 5.0]);
    assert_eq!(mobility, vec![1.0, 1.2]);
}

#[test]
fn clear_resets_all_reusable_storage() {
    let mut buffer = buffer(vec![ImsPeak {
        mz: 100.0,
        intensity: 10.0,
        im: 1.0,
    }]);
    buffer.agg_buff.push(buffer.peaks[0]);

    buffer.clear();

    assert!(buffer.peaks.is_empty());
    assert!(buffer.order.is_empty());
    assert!(buffer.agg_buff.is_empty());
}
