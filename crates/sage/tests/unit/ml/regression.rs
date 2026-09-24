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

#[test]
fn zero_variance_and_constant_columns_do_not_break_the_fit() {
    // y = 2 x - 3 z + 1, plus an all-zero column and a second constant column.
    let items: Vec<(f64, f64)> = (0..500)
        .map(|i| {
            let x = i as f64 / 50.0;
            (x, ((i * 7) % 13) as f64)
        })
        .collect();
    let y = |&(x, z): &(f64, f64)| 2.0 * x - 3.0 * z + 1.0;
    let lr = LinearRegression::fit::<_, 5>(&items, |_| true, |&(x, z)| [0.0, x, 5.0, z, 1.0], y)
        .unwrap();
    assert_eq!(lr.beta[0], 0.0);
    assert_eq!(lr.beta[4], 0.0);
    assert!((lr.beta[1] - 2.0).abs() < 1e-9, "{:?}", lr.beta);
    assert!((lr.beta[3] + 3.0).abs() < 1e-9, "{:?}", lr.beta);
    // The first constant column carries the intercept.
    assert!((5.0 * lr.beta[2] - 1.0).abs() < 1e-9, "{:?}", lr.beta);
    assert!((lr.r2 - 1.0).abs() < 1e-12);

    // The intercept is placed on the first constant column, scaled by its value.
    let lr = LinearRegression::fit::<_, 3>(&items, |_| true, |&(x, z)| [x, 4.0, z], y).unwrap();
    assert_eq!(lr.beta.len(), 3);
    assert!((4.0 * lr.beta[1] - 1.0).abs() < 1e-9, "{:?}", lr.beta);
}

#[test]
fn fit_without_intercept_goes_through_the_origin() {
    let items: Vec<(f64, f64)> = (1..200).map(|i| (i as f64, (i % 7) as f64)).collect();
    let lr = LinearRegression::fit::<_, 2>(
        &items,
        |_| true,
        |&(x, z)| [x, z],
        |&(x, z)| 0.5 * x + 2.0 * z,
    )
    .unwrap();
    assert!((lr.beta[0] - 0.5).abs() < 1e-9, "{:?}", lr.beta);
    assert!((lr.beta[1] - 2.0).abs() < 1e-9, "{:?}", lr.beta);
}

#[test]
fn duplicated_columns_share_the_coefficient_and_keep_predictions() {
    // Exact duplicates are a pure null-space direction: the minimum-norm
    // standardized solution splits the effect evenly.
    let items: Vec<f64> = (0..300).map(|i| i as f64 / 10.0).collect();
    let single =
        LinearRegression::fit::<_, 2>(&items, |_| true, |&x| [x, 1.0], |&x| 3.0 * x + 1.0).unwrap();
    let double =
        LinearRegression::fit::<_, 3>(&items, |_| true, |&x| [x, x, 1.0], |&x| 3.0 * x + 1.0)
            .unwrap();
    eprintln!("duplicate split: {:?} vs {:?}", double.beta, single.beta);
    assert!((double.beta[0] - double.beta[1]).abs() < 1e-6);
    assert!((double.beta[0] + double.beta[1] - single.beta[0]).abs() < 1e-9);
    assert!((double.beta[2] - single.beta[1]).abs() < 1e-9);
}

#[test]
fn cholesky_solves_positive_definite_systems_and_rejects_others() {
    // [[4, 2], [2, 3]] x = [10, 8] -> x = [1.75, 1.5]
    let x = cholesky_solve(vec![4.0, 2.0, 2.0, 3.0], vec![10.0, 8.0]).unwrap();
    assert!(
        (x[0] - 1.75).abs() < 1e-12 && (x[1] - 1.5).abs() < 1e-12,
        "{x:?}"
    );
    assert!(cholesky_solve(vec![1.0, 1.0, 1.0, 1.0], vec![1.0, 1.0]).is_none());
    assert!(cholesky_solve(vec![-1.0], vec![1.0]).is_none());
    assert!(cholesky_solve(vec![f64::NAN], vec![1.0]).is_none());
    assert_eq!(cholesky_solve(vec![], vec![]), Some(vec![]));
}

#[test]
fn collinear_fit_is_bitwise_deterministic_across_thread_counts() {
    let items: Vec<[f64; 3]> = (0..100_000)
        .map(|i| {
            let x = i as f64;
            [
                (x * 0.37).sin() * 1e3,
                (x * 0.011).cos(),
                x.sqrt() + (x * 1.3).sin(),
            ]
        })
        .collect();
    // Column 2 = 3 * column 0 - column 1 + 2 * intercept: exactly collinear.
    let fit = || {
        LinearRegression::fit::<_, 4>(
            &items,
            |_| true,
            |row| [row[0], row[1], 3.0 * row[0] - row[1] + 2.0, 1.0],
            |row| row[2],
        )
        .unwrap()
    };
    let bits = |lr: &LinearRegression| {
        (
            lr.beta.iter().map(|b| b.to_bits()).collect::<Vec<_>>(),
            lr.r2.to_bits(),
        )
    };
    let expected = bits(&fit());
    for threads in [1, 2, 3, 8, 16] {
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(threads)
            .build()
            .unwrap();
        for _ in 0..3 {
            assert_eq!(bits(&pool.install(fit)), expected, "{threads} threads");
        }
    }
}
