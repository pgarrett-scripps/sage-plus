use super::{ms2_without_precursor_warnings, release_unread_ms1_columns};
use sage_core::spectrum::{Precursor, ProcessedSpectrum};

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

fn ms2(file_id: usize, peaks: usize, precursor: bool) -> ProcessedSpectrum {
    ProcessedSpectrum {
        level: 2,
        file_id,
        masses: vec![100.0; peaks],
        intensities: vec![1.0; peaks],
        precursors: match precursor {
            true => vec![Precursor {
                mz: 500.0,
                ..Default::default()
            }],
            false => Vec::new(),
        },
        ..Default::default()
    }
}

/// One warning per file, in file order, counting only the scans dropped for
/// want of a precursor.
#[test]
fn precursor_less_ms2_scans_are_counted_per_file() {
    let spectra = vec![
        // File 2 first in the list, to show the messages come out in file order.
        ms2(2, 20, false),
        ms2(0, 20, false),
        ms2(0, 20, false),
        ms2(0, 20, false),
        ms2(0, 20, true),
        // Dropped for too few peaks, which is a different reason.
        ms2(1, 3, false),
        // An MS1 scan has no precursor by definition.
        ProcessedSpectrum {
            level: 1,
            file_id: 1,
            masses: vec![100.0; 20],
            intensities: vec![1.0; 20],
            ..Default::default()
        },
    ];

    let warnings =
        ms2_without_precursor_warnings(&spectra, 15, |file_id| format!("file-{file_id}"));

    assert_eq!(
        warnings,
        vec![
            "file-0: 3 MS2 scans were not searched because no precursor was recorded for them"
                .to_string(),
            "file-2: 1 MS2 scan was not searched because no precursor was recorded for it"
                .to_string(),
        ]
    );
}

/// Nothing to say about a file whose MS2 scans all record a precursor.
#[test]
fn precursor_less_ms2_warnings_are_empty_when_every_scan_has_one() {
    let spectra = vec![ms2(0, 20, true), ms2(1, 20, true)];

    assert!(
        ms2_without_precursor_warnings(&spectra, 15, |file_id| format!("file-{file_id}"))
            .is_empty()
    );
    assert!(ms2_without_precursor_warnings(&[], 15, |_| String::new()).is_empty());
}
