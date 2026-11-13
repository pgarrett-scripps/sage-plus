//! Small numeric helpers shared by the runner and its submodules.
//!
//! Moved out of `runner.rs` unchanged.

use super::*;

pub(super) fn median_finite(values: impl IntoIterator<Item = f32>) -> Option<f32> {
    let mut values = values
        .into_iter()
        .filter(|value| value.is_finite())
        .collect::<Vec<_>>();
    values.sort_by(f32::total_cmp);
    let middle = values.len() / 2;
    match values.len() {
        0 => None,
        length if length % 2 == 0 => Some((values[middle - 1] + values[middle]) / 2.0),
        _ => Some(values[middle]),
    }
}

pub(super) fn average_finite(values: impl IntoIterator<Item = f32>) -> Option<f32> {
    let values = values
        .into_iter()
        .filter(|value| value.is_finite())
        .collect::<Vec<_>>();
    (!values.is_empty()).then(|| values.iter().sum::<f32>() / values.len() as f32)
}

pub(super) fn normalize_finite(values: Vec<f64>) -> Vec<f64> {
    if values.is_empty() {
        return values;
    }
    let min = values.iter().copied().fold(f64::INFINITY, f64::min);
    let max = values.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    let range = max - min;
    if !range.is_finite() || range <= f64::EPSILON {
        return vec![0.0; values.len()];
    }
    values
        .into_iter()
        .map(|value| (value - min) / range)
        .collect()
}

pub(super) fn labeled_finite_values(
    features: &[Feature],
    value: impl Fn(&Feature) -> f64,
) -> (Vec<f64>, Vec<i32>) {
    features
        .iter()
        .filter_map(|feature| {
            let value = value(feature);
            (value.is_finite() && matches!(feature.label, -1 | 1)).then_some((value, feature.label))
        })
        .unzip()
}

#[cfg(test)]
#[path = "../../tests/unit/runner/util.rs"]
mod test;
