use super::*;

fn feature(label: i32, discriminant_score: f32) -> Feature {
    Feature {
        label,
        discriminant_score,
        ..Feature::default()
    }
}

#[test]
fn equal_scores_receive_the_same_q_value() {
    let mut features = vec![
        feature(1, 10.0),
        feature(-1, 10.0),
        feature(1, 9.0),
        feature(1, 9.0),
    ];

    spectrum_q_value(&mut features);

    assert_eq!(features[0].spectrum_q, features[1].spectrum_q);
}

#[test]
fn tied_score_order_does_not_change_q_values() {
    let mut target_first = vec![
        feature(1, 10.0),
        feature(-1, 10.0),
        feature(1, 9.0),
        feature(1, 9.0),
    ];
    let mut decoy_first = vec![
        feature(-1, 10.0),
        feature(1, 10.0),
        feature(1, 9.0),
        feature(1, 9.0),
    ];

    spectrum_q_value(&mut target_first);
    spectrum_q_value(&mut decoy_first);

    let target_first_q = target_first
        .iter()
        .map(|feature| feature.spectrum_q)
        .collect::<Vec<_>>();
    let decoy_first_q = decoy_first
        .iter()
        .map(|feature| feature.spectrum_q)
        .collect::<Vec<_>>();
    assert_eq!(target_first_q, decoy_first_q);
}

#[test]
fn alternate_score_supports_provisional_fdr() {
    let mut features = vec![feature(1, 0.0), feature(-1, 0.0)];
    features[0].poisson = -5.0;
    features[1].poisson = -4.0;

    spectrum_q_value_by(&mut features, |feature| feature.poisson);

    assert!(features
        .iter()
        .all(|feature| feature.spectrum_q.is_finite()));
}
