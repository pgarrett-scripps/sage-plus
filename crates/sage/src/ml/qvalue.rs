use crate::scoring::Feature;
use std::cmp::Ordering;

/// Assign q_values in place to a set of PSMs, returning the number of PSMs
/// q <= 0.01
///
/// # Invariants
/// * `scores` must be sorted in descending order (e.g. best PSM is first)
pub fn spectrum_q_value(scores: &mut [Feature]) -> usize {
    spectrum_q_value_by(scores, |feature| f64::from(feature.discriminant_score))
}

/// Assign spectrum q-values while treating adjacent equal scores as one group.
/// `scores` must already be sorted from best to worst by `score`.
pub fn spectrum_q_value_by(scores: &mut [Feature], score: impl Fn(&Feature) -> f64) -> usize {
    // FDR Calculation:
    // * Sort by score, descending
    // * Estimate FDR after each complete tied-score group
    // * Calculate q-value

    let mut decoy = 1;
    let mut target = 0;

    let mut start = 0;
    while start < scores.len() {
        let tied_score = score(&scores[start]);
        let mut end = start + 1;
        while end < scores.len() && score(&scores[end]).total_cmp(&tied_score) == Ordering::Equal {
            end += 1;
        }
        for feature in &scores[start..end] {
            if feature.label == -1 {
                decoy += 1;
            } else {
                target += 1;
            }
        }
        let fdr = decoy as f32 / target as f32;
        for feature in &mut scores[start..end] {
            feature.spectrum_q = fdr;
        }
        start = end;
    }

    // Reverse slice, and calculate the cumulative minimum
    let mut q_min = 1.0f32;
    let mut passing = 0;
    for score in scores.iter_mut().rev() {
        q_min = q_min.min(score.spectrum_q);
        score.spectrum_q = q_min;
        if q_min <= 0.01 && score.label != -1 {
            passing += 1;
        }
    }
    passing
}

#[cfg(test)]
mod tests {
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
}
