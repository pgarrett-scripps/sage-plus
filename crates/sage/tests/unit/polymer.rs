use super::*;

/// Sorted peaks from `(mz, intensity)` pairs.
fn spectrum(mut peaks: Vec<(f64, f32)>) -> (Vec<f32>, Vec<f32>) {
    peaks.sort_by(|a, b| a.0.total_cmp(&b.0));
    peaks.into_iter().map(|(mz, i)| (mz as f32, i)).unzip()
}

fn ladder(polymer: Polymer, charge: f64, adduct: f64, n: std::ops::Range<usize>) -> Vec<f64> {
    n.map(|n| (polymer.end_group + n as f64 * polymer.repeat + charge * adduct) / charge)
        .collect()
}

#[test]
fn peg_ladder_is_found_and_counted_once() {
    let mut peaks = ladder(PEG, 1.0, 1.007_276, 8..14)
        .into_iter()
        .map(|mz| (mz, 100.0))
        .collect::<Vec<_>>();
    // Doubly sodiated PEG ladder.
    peaks.extend(
        ladder(PEG, 2.0, 22.989_218, 20..26)
            .into_iter()
            .map(|mz| (mz, 50.0)),
    );
    // Peptide-like peaks that are not on a ladder.
    peaks.extend([(501.27, 1000.0), (733.41, 500.0)]);
    let (mz, intensity) = spectrum(peaks);

    let scanner = PolymerScanner::default();
    let found = scanner.scan(&mz, &intensity);
    assert_eq!(found.len(), 3);
    assert!(
        (found[0] - (6.0 * 100.0 + 6.0 * 50.0)).abs() < 1e-6,
        "{found:?}"
    );
    assert_eq!(found[1], 0.0);
    assert_eq!(found[2], 0.0);
}

#[test]
fn short_ladders_do_not_count() {
    let (mz, intensity) = spectrum(
        ladder(POLYSILOXANE, 1.0, 1.007_276, 5..8)
            .into_iter()
            .map(|mz| (mz, 10.0))
            .collect(),
    );
    assert_eq!(
        PolymerScanner::default().scan(&mz, &intensity),
        vec![0.0; 3]
    );
    // One more member reaches the four-peak minimum.
    let (mz, intensity) = spectrum(
        ladder(POLYSILOXANE, 1.0, 1.007_276, 5..9)
            .into_iter()
            .map(|mz| (mz, 10.0))
            .collect(),
    );
    assert_eq!(PolymerScanner::default().scan(&mz, &intensity)[2], 40.0);
    // D6 cyclic siloxane is the familiar 445.12 lock mass.
    assert!((ladder(POLYSILOXANE, 1.0, 1.007_276, 6..7)[0] - 445.1200).abs() < 1e-3);
}

#[test]
fn file_stats_report_percent_of_tic() {
    let mut stats = PolymerFileStats::new(&DEFAULT_POLYMERS);
    stats.add(1000.0, &[100.0, 0.0, 0.0]);
    stats.add(1000.0, &[0.0, 0.0, 50.0]);
    stats.finish();
    assert_eq!(stats.ms1_spectra, 2);
    assert!((stats.polymers[0].tic_pct - 5.0).abs() < 1e-9);
    assert!((stats.polymers[2].tic_pct - 2.5).abs() < 1e-9);
    assert_eq!(PolymerScanner::default().scan(&[], &[]), vec![0.0; 3]);
}
