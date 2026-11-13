use super::release_unread_ms1_columns;
use sage_core::spectrum::ProcessedSpectrum;

/// The columns LFQ reads survive; the charge columns it never reads are freed.
#[test]
fn releasing_ms1_columns_keeps_what_lfq_reads() {
    let mut spectra = vec![ProcessedSpectrum {
        level: 1,
        file_id: 3,
        scan_start_time: 12.5,
        masses: vec![400.0, 500.0, 600.0],
        intensities: vec![1.0, 2.0, 3.0],
        charges: vec![1, 1, 1],
        charge_is_known: vec![false, false, false],
        mobilities: vec![0.8, 0.9, 1.0],
        ..Default::default()
    }];

    release_unread_ms1_columns(&mut spectra);

    let spectrum = &spectra[0];
    assert_eq!(spectrum.file_id, 3);
    assert_eq!(spectrum.scan_start_time, 12.5);
    assert_eq!(spectrum.masses, vec![400.0, 500.0, 600.0]);
    assert_eq!(spectrum.intensities, vec![1.0, 2.0, 3.0]);
    assert_eq!(spectrum.mobilities, vec![0.8, 0.9, 1.0]);
    assert!(spectrum.charges.is_empty());
    assert!(spectrum.charge_is_known.is_empty());
    assert_eq!(spectrum.charges.capacity(), 0);
    assert_eq!(spectrum.charge_is_known.capacity(), 0);
}

/// A spectrum whose charge columns were already empty is left alone.
#[test]
fn releasing_ms1_columns_is_idempotent() {
    let mut spectra = vec![ProcessedSpectrum {
        level: 1,
        masses: vec![400.0],
        intensities: vec![1.0],
        ..Default::default()
    }];

    release_unread_ms1_columns(&mut spectra);
    release_unread_ms1_columns(&mut spectra);

    assert_eq!(spectra[0].masses, vec![400.0]);
    assert!(spectra[0].charges.is_empty());
}
