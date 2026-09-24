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

#[test]
fn fit_is_bitwise_deterministic_across_parallel_runs() {
    // Ill-conditioned, large-magnitude rows make floating-point summation
    // order visible in the fitted coefficients.
    let items: Vec<[f64; 4]> = (0..200_000)
        .map(|i| {
            let x = i as f64;
            [
                (x * 0.37).sin() * 1e3,
                (x * 0.011).cos() * 1e-2 + x * 1e-5,
                ((x * 0.7).sin() * 1e4).fract(),
                x.sqrt() * 1e2 + (x * 1.3).sin(),
            ]
        })
        .collect();
    let fit = || {
        LinearRegression::fit::<_, 4>(
            &items,
            |row| row[2] > -0.9,
            |row| [row[0], row[1], row[2], 1.0],
            |row| row[3],
        )
        .unwrap()
    };
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(8)
        .build()
        .unwrap();
    let bits = |lr: &LinearRegression| {
        (
            lr.beta.iter().map(|b| b.to_bits()).collect::<Vec<_>>(),
            lr.r2.to_bits(),
        )
    };
    let expected = bits(&pool.install(fit));
    for _ in 0..20 {
        assert_eq!(bits(&pool.install(fit)), expected);
    }
    // Independent of the number of worker threads, too.
    let serial = rayon::ThreadPoolBuilder::new()
        .num_threads(1)
        .build()
        .unwrap();
    assert_eq!(bits(&serial.install(fit)), expected);
}
