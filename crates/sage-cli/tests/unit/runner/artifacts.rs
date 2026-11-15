use super::*;

use sage_core::ptm::ModLocalization;

fn localization(decoy_winner: bool, localization_q_value: f32) -> ModLocalization {
    ModLocalization {
        mass: 79.966_33,
        label: Some("Phospho".into()),
        site_count: 1,
        candidate_sites: 2,
        site_determining_ions: 10,
        site_determining_matched: 6,
        delta_score: 12.0,
        target_decoy_score: 4.0,
        decoy_winner,
        competition_eligible: true,
        separating_margin: 2,
        localization_q_value,
        best_sites: Vec::new(),
        all_sites: Vec::new(),
    }
}

/// The two reasons a localization reaches no site report, and nothing else.
#[test]
fn a_localization_is_declined_for_a_decoy_win_or_a_high_q_value() {
    assert_eq!(site_report_decline(&localization(false, 0.004), 0.01), None);
    assert_eq!(site_report_decline(&localization(false, 0.01), 0.01), None);
    assert_eq!(
        site_report_decline(&localization(false, 0.011), 0.01),
        Some(Decline::BelowLocalizationQValue)
    );
    assert_eq!(
        site_report_decline(&localization(true, 0.0), 0.01),
        Some(Decline::DecoyArrangementWon)
    );
    // A decoy win is not a threshold: loosening the cutoff does not report it.
    assert_eq!(
        site_report_decline(&localization(true, 0.0), 1.0),
        Some(Decline::DecoyArrangementWon)
    );
    // A row that fails both is counted once, as the decoy win it is.
    assert_eq!(
        site_report_decline(&localization(true, 0.9), 0.01),
        Some(Decline::DecoyArrangementWon)
    );
}

/// `localization_declines` must count exactly the modifications
/// `collect_site_rows` leaves out, which is why both read the same predicate.
#[test]
fn the_counts_are_what_the_site_report_leaves_out() {
    let cutoff = 0.01;
    let mods = vec![
        localization(false, 0.001),
        localization(false, 0.5),
        localization(true, 0.0),
        localization(true, 0.5),
        localization(false, 0.009),
    ];
    let declined = mods
        .iter()
        .filter(|m| site_report_decline(m, cutoff).is_some())
        .count();
    let kept = mods.len() - declined;
    assert_eq!(kept, 2, "two localizations are reportable");
    let wins = mods
        .iter()
        .filter(|m| site_report_decline(m, cutoff) == Some(Decline::DecoyArrangementWon))
        .count();
    let below = mods
        .iter()
        .filter(|m| site_report_decline(m, cutoff) == Some(Decline::BelowLocalizationQValue))
        .count();
    assert_eq!((wins, below), (2, 1));
    assert_eq!(
        wins + below,
        declined,
        "every decline has exactly one reason"
    );
}
