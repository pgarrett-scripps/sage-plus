use super::{RunSummary, ToleranceRecommendation};
use sage_core::scoring::Feature;

#[test]
fn older_run_summaries_receive_compatible_defaults() {
    let summary: RunSummary = serde_json::from_value(serde_json::json!({
        "runtime_secs": 1,
        "files": 1,
        "peptides_in_database": 10,
        "fragments_in_database": 20,
        "psms_at_one_percent_fdr": 2,
        "peptides_at_one_percent_fdr": 1,
        "proteins_at_one_percent_fdr": 1,
        "protein_groups_at_one_percent_fdr": 1,
        "output_paths": []
    }))
    .unwrap();

    assert_eq!(summary.schema_version, 1);
    assert_eq!(
        summary.recommended_tolerances,
        ToleranceRecommendation::default()
    );
    assert!(!summary.ptm_localization.enabled);
    assert_eq!(summary.quantification.lfq_features, 0);
}

#[test]
fn tolerance_recommendation_uses_confident_signed_errors() {
    let confident = |i: usize| Feature {
        rank: 1,
        label: 1,
        spectrum_q: 0.001,
        // Precursor: -2 ppm bias, errors -3/-2/-1.
        delta_mass: -2.0 + (i % 3) as f32 - 1.0,
        // Fragment: +3 ppm bias; each PSM's ions spread by 2 ppm.
        signed_fragment_ppm: 3.0,
        average_ppm: 3.0,
        fragment_ppm_sd: 2.0,
        ..Feature::default()
    };
    let mut features = (0..150).map(confident).collect::<Vec<_>>();
    // Decoys, lower ranks and non-confident PSMs are ignored.
    for excluded in [
        Feature {
            label: -1,
            delta_mass: 90.0,
            ..confident(0)
        },
        Feature {
            rank: 2,
            delta_mass: 90.0,
            ..confident(0)
        },
        Feature {
            spectrum_q: 0.5,
            delta_mass: 90.0,
            ..confident(0)
        },
    ] {
        features.extend(std::iter::repeat_n(excluded, 60));
    }

    let recommendation = ToleranceRecommendation::from_features(&features);
    assert_eq!(recommendation.psms, 150);
    assert!(recommendation.skipped.is_none());
    let precursor = recommendation.precursor.unwrap();
    assert_eq!(precursor.bias_ppm, -2.0);
    // |-2| + 4 * 1.4826 = 7.93 -> 10 ppm.
    assert_eq!(precursor.recommended_ppm, Some(10.0));
    let fragment = recommendation.fragment.unwrap();
    assert_eq!(fragment.bias_ppm, 3.0);
    // |3| + 4 * 2 = 11 -> 20 ppm.
    assert_eq!(fragment.recommended_ppm, Some(20.0));
    assert_eq!(
        recommendation.log_line(),
        "recommended tolerances: precursor ±10 ppm, fragment ±20 ppm (from 150 PSMs)"
    );

    let few = ToleranceRecommendation::from_features(&features[..99]);
    assert_eq!(few.psms, 99);
    assert!(few.precursor.is_none() && few.fragment.is_none());
    assert!(few.log_line().contains("skipped (99 confident PSMs"));
}
