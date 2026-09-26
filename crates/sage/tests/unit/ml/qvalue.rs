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

fn discriminant(feature: &Feature) -> f64 {
    f64::from(feature.discriminant_score)
}

fn is_decoy(feature: &Feature) -> bool {
    feature.label == -1
}

#[test]
fn one_group_matches_ungrouped_q_values() {
    // Deterministic pseudo-random scores on a coarse grid, so ties occur
    let mut state = 0x2545_f491_u64;
    let unsorted = (0..500)
        .map(|_| {
            state = state
                .wrapping_mul(6_364_136_223_846_793_005)
                .wrapping_add(1);
            let label = if (state >> 40) % 3 == 0 { -1 } else { 1 };
            feature(label, ((state >> 20) % 60) as f32 / 4.0)
        })
        .collect::<Vec<_>>();

    // Grouped input is unsorted; the ungrouped function needs sorted input
    let grouped = grouped_q_values(&unsorted, discriminant, is_decoy, |_| ());

    let mut order = (0..unsorted.len()).collect::<Vec<_>>();
    order.sort_by(|&a, &b| discriminant(&unsorted[b]).total_cmp(&discriminant(&unsorted[a])));
    let mut sorted = order
        .iter()
        .map(|&index| unsorted[index].clone())
        .collect::<Vec<_>>();
    spectrum_q_value(&mut sorted);

    assert!(sorted.iter().any(|feature| feature.spectrum_q < 1.0));
    for (feature, &index) in sorted.iter().zip(&order) {
        assert_eq!(feature.spectrum_q, grouped[index]);
    }
}

#[test]
fn groups_are_estimated_separately() {
    // (group, label, score), interleaved across groups
    let items = [
        ('a', 1, 10.0),
        ('b', 1, 5.0),
        ('a', 1, 9.0),
        ('c', 1, 3.0),
        ('b', -1, 4.0),
        ('a', 1, 8.0),
        ('c', -1, 3.0),
        ('b', -1, 3.0),
        ('a', 1, 7.0),
        ('c', 1, 2.0),
        ('b', 1, 2.0),
        ('a', -1, 6.0),
        ('c', 1, 1.0),
    ];
    let q = grouped_q_values(&items, |item| item.2, |item| item.1 == -1, |item| item.0);

    // a: 10T 9T 8T 7T 6D -> FDR 1, 1/2, 1/3, 1/4, 2/4 -> q 1/4 for targets
    // b: 5T 4D 3D 2T -> FDR 1, 2, 3, 3/2 -> q 1 throughout
    // c: {3T 3D} 2T 1T -> FDR 2, 1, 2/3 -> q 2/3 throughout, ties shared
    let expected: [f32; 13] = [
        0.25,
        1.0,
        0.25,
        2.0 / 3.0,
        1.0,
        0.25,
        2.0 / 3.0,
        1.0,
        0.25,
        2.0 / 3.0,
        1.0,
        0.5,
        2.0 / 3.0,
    ];
    assert_eq!(q, expected);

    // Pooled, group b's targets would borrow group a's confidence
    let pooled = grouped_q_values(&items, |item| item.2, |item| item.1 == -1, |_| ());
    assert!(pooled[1] < q[1]);
}

#[test]
fn empty_input_and_singleton_groups() {
    let empty: [(u8, bool, f64); 0] = [];
    assert!(grouped_q_values(&empty, |item| item.2, |item| item.1, |item| item.0).is_empty());

    // A lone target has FDR (0 + 1) / 1
    let lone = [(0u8, false, 1.0), (1u8, false, 2.0)];
    let q = grouped_q_values(&lone, |item| item.2, |item| item.1, |item| item.0);
    assert_eq!(q, [1.0f32, 1.0]);
}

#[test]
fn all_decoy_group_gets_q_one_and_leaves_other_groups_alone() {
    let items = [
        (0u8, true, 9.0),
        (1u8, false, 8.0),
        (0u8, true, 7.0),
        (1u8, false, 6.0),
        (1u8, false, 5.0),
        (1u8, false, 4.0),
    ];
    let q = grouped_q_values(&items, |item| item.2, |item| item.1, |item| item.0);
    assert_eq!(q, [1.0f32, 0.25, 1.0, 0.25, 0.25, 0.25]);
    assert!(q.iter().all(|q| q.is_finite()));
}
