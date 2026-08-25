use super::*;

#[test]
fn fit_perfect_line() {
    // y = 2 x + 1, with intercept embedded as the last column.
    let items: Vec<(f64, f64)> = (0..50).map(|i| (i as f64, 2.0 * i as f64 + 1.0)).collect();
    let lr =
        LinearRegression::fit::<_, 2>(&items, |_| true, |&(x, _)| [x, 1.0], |&(_, y)| y).unwrap();
    assert!((lr.beta[0] - 2.0).abs() < 1e-9, "slope: {}", lr.beta[0]);
    assert!((lr.beta[1] - 1.0).abs() < 1e-9, "intercept: {}", lr.beta[1]);
    assert!((lr.r2 - 1.0).abs() < 1e-9, "r2: {}", lr.r2);
}

#[test]
fn fit_with_noise() {
    // y ~= 3 x + 2 with a deterministic perturbation; r^2 should be high.
    let items: Vec<(f64, f64)> = (0..200)
        .map(|i| {
            let x = i as f64 / 10.0;
            let noise = ((i as f64) * 0.7).sin() * 0.1;
            (x, 3.0 * x + 2.0 + noise)
        })
        .collect();
    let lr =
        LinearRegression::fit::<_, 2>(&items, |_| true, |&(x, _)| [x, 1.0], |&(_, y)| y).unwrap();
    assert!((lr.beta[0] - 3.0).abs() < 0.05, "slope: {}", lr.beta[0]);
    assert!((lr.beta[1] - 2.0).abs() < 0.1, "intercept: {}", lr.beta[1]);
    assert!(lr.r2 > 0.99, "r2: {}", lr.r2);
}

#[test]
fn empty_filter_returns_none() {
    let items: Vec<f64> = vec![1.0, 2.0, 3.0];
    let lr = LinearRegression::fit::<_, 1>(&items, |_| false, |_| [1.0], |&y| y);
    assert!(lr.is_none());
}
