use super::*;
use crate::in_source::{FALLBACK_RT_FRACTION, MIN_SPREAD_OBSERVATIONS};
use crate::lfq::{FileEvidence, Peak};

fn stats(charges: &[u8]) -> ChromatographyStats {
    let mut stats = ChromatographyStats::default();
    for &charge in charges {
        stats.add_charge(charge);
    }
    stats
}

#[test]
fn charge_four_and_above_share_the_last_bucket() {
    let stats = stats(&[1, 2, 2, 3, 4, 5, 9]);
    assert_eq!(stats.confident_precursors, 7);
    assert_eq!(stats.charges, [1, 2, 1, 3]);
    assert_eq!(stats.charge_pct(2), 100.0 * 2.0 / 7.0);
    assert_eq!(stats.charge_pct(4), 100.0 * 3.0 / 7.0);
    // Asking past the last bucket is zero, not a panic.
    assert_eq!(stats.charge_pct(5), 0.0);
    assert_eq!(stats.charge_unknown, 0);
}

#[test]
fn an_unstated_charge_is_counted_apart_so_the_shares_still_sum() {
    let stats = stats(&[0, 0, 2, 2]);
    assert_eq!(stats.charge_unknown, 2);
    assert_eq!(stats.charge_unknown_pct(), 50.0);
    let sum = (1..=4).map(|c| stats.charge_pct(c)).sum::<f64>() + stats.charge_unknown_pct();
    assert!((sum - 100.0).abs() < 1e-9, "{sum}");
}

#[test]
fn an_empty_file_divides_by_nothing() {
    let stats = ChromatographyStats::default();
    assert_eq!(stats.charge_pct(2), 0.0);
    assert_eq!(stats.charge_unknown_pct(), 0.0);
    assert_eq!(stats.transferred_pct(), None);
}

#[test]
fn the_transferred_share_needs_quantification_and_a_denominator() {
    let mut stats = ChromatographyStats::default();
    assert_eq!(stats.transferred_pct(), None);
    stats.quantified_precursors = Some(0);
    stats.transferred_precursors = Some(0);
    // Nothing quantified is not a share of zero percent.
    assert_eq!(stats.transferred_pct(), None);
    stats.quantified_precursors = Some(8);
    stats.transferred_precursors = Some(2);
    assert_eq!(stats.transferred_pct(), Some(25.0));
}

#[test]
fn the_spread_is_the_median_pairwise_half_distance() {
    let mut rts = RepeatRts::default();
    // Enough pairs to derive rather than fall back: 20 peptides seen twice,
    // each pair 1.0 apart, so every half distance is 0.5.
    for index in 0..MIN_SPREAD_OBSERVATIONS {
        let rt = index as f32;
        rts.insert(PeptideIx(index as u32), vec![rt, rt + 1.0]);
    }
    let (spread, kind) = repeat_psm_spread(&rts).unwrap();
    assert_eq!(kind, PeakWidth::RepeatPsmSpread);
    assert_eq!(spread, 0.5);
}

#[test]
fn too_few_repeats_fall_back_to_a_fraction_of_the_gradient() {
    let mut rts = RepeatRts::default();
    rts.insert(PeptideIx(0), vec![0.0, 2.0]);
    rts.insert(PeptideIx(1), vec![100.0]);
    let (spread, kind) = repeat_psm_spread(&rts).unwrap();
    assert_eq!(kind, PeakWidth::RepeatPsmSpreadFallback);
    assert_eq!(spread, 100.0 * FALLBACK_RT_FRACTION);
}

#[test]
fn a_file_with_no_psm_has_no_spread() {
    assert_eq!(repeat_psm_spread(&RepeatRts::default()), None);
}

#[test]
fn the_peak_width_says_which_quantity_it_is() {
    assert_eq!(PeakWidth::LfqMedianFwhm.label(), "lfq_median_fwhm");
    assert_eq!(PeakWidth::RepeatPsmSpread.label(), "repeat_psm_spread");
    assert_eq!(
        PeakWidth::RepeatPsmSpreadFallback.label(),
        "repeat_psm_spread_fallback"
    );
}

fn evidence(transfer_candidate: bool, fwhm: Option<f32>) -> FileEvidence {
    FileEvidence {
        transfer_candidate,
        fwhm,
        ..Default::default()
    }
}

fn quantified(
    intensities: Vec<Option<f64>>,
    file_evidence: Vec<Option<FileEvidence>>,
) -> QuantifiedPeak {
    QuantifiedPeak {
        peak: Peak::default(),
        intensities,
        ms2_confirmed: Vec::new(),
        ms2_confirmed_strict: Vec::new(),
        file_evidence,
        paired_decoy_evidence: Vec::new(),
    }
}

#[test]
fn quantification_counts_transfers_and_takes_the_median_width() {
    let mut areas = std::collections::HashMap::new();
    // File 0 identified both precursors; file 1 only got them by transfer, and
    // one of them has no integrated intensity there at all.
    areas.insert(
        (PrecursorId::Combined(PeptideIx(0)), false),
        quantified(
            vec![Some(10.0), Some(5.0)],
            vec![
                Some(evidence(false, Some(0.2))),
                Some(evidence(true, Some(0.4))),
            ],
        ),
    );
    areas.insert(
        (PrecursorId::Combined(PeptideIx(1)), false),
        quantified(
            vec![Some(20.0), None],
            vec![Some(evidence(false, Some(0.6))), None],
        ),
    );
    // A decoy is not a quantified precursor.
    areas.insert(
        (PrecursorId::Combined(PeptideIx(2)), true),
        quantified(
            vec![Some(99.0), Some(99.0)],
            vec![Some(evidence(true, Some(9.0))), Some(evidence(true, None))],
        ),
    );
    let by_file = quantification_by_file(&areas, 2);
    assert_eq!(by_file[0], (2, 0, Some(0.4)));
    assert_eq!(by_file[1], (1, 1, Some(0.4)));
}

#[test]
fn a_peak_with_no_usable_width_still_counts_as_quantified() {
    let mut areas = std::collections::HashMap::new();
    areas.insert(
        (PrecursorId::Combined(PeptideIx(0)), false),
        quantified(
            vec![Some(10.0)],
            vec![Some(evidence(false, Some(f32::NAN)))],
        ),
    );
    areas.insert(
        (PrecursorId::Combined(PeptideIx(1)), false),
        quantified(vec![Some(10.0)], vec![Some(evidence(false, Some(0.0)))]),
    );
    let by_file = quantification_by_file(&areas, 1);
    assert_eq!(by_file[0], (2, 0, None));
}
