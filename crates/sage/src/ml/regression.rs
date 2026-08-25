//! Streaming OLS linear regression.
//!
//! Fits `beta = (X^T X)^-1 X^T y` without materializing the `n x D` design
//! matrix: one parallel fold/reduce pass accumulates `X^T X`, `X^T y`,
//! `sum(y)`, `sum(y^2)`, `n` from per-row outer products. A second pass
//! evaluates `sum((X beta - y)^2)` to report r^2 that matches the
//! materialized formula numerically.
//!
//! Per-worker scratch is `O(D^2)` (the cov accumulator), independent of `n`.

use super::{gauss::Gauss, matrix::Matrix};
use rayon::prelude::*;

pub struct LinearRegression {
    pub beta: Vec<f64>,
    pub r2: f64,
}

struct Acc {
    cov: Vec<f64>, // D*D row-major
    b: Vec<f64>,   // D
    sum_y: f64,
    sum_y2: f64,
    n: usize,
}

impl Acc {
    fn zero(d: usize) -> Self {
        Self {
            cov: vec![0.0; d * d],
            b: vec![0.0; d],
            sum_y: 0.0,
            sum_y2: 0.0,
            n: 0,
        }
    }

    fn add_row(&mut self, row: &[f64], y: f64) {
        let d = row.len();
        for j in 0..d {
            let rj = row[j];
            self.b[j] += rj * y;
            let off = j * d;
            for k in 0..d {
                self.cov[off + k] += rj * row[k];
            }
        }
        self.sum_y += y;
        self.sum_y2 += y * y;
        self.n += 1;
    }

    fn merge(mut self, other: Acc) -> Acc {
        for i in 0..self.cov.len() {
            self.cov[i] += other.cov[i];
        }
        for i in 0..self.b.len() {
            self.b[i] += other.b[i];
        }
        self.sum_y += other.sum_y;
        self.sum_y2 += other.sum_y2;
        self.n += other.n;
        self
    }
}

impl LinearRegression {
    /// Fit OLS over `items` with predicate `filter`. `embed(item)` produces a
    /// design row of length `D`; `target(item)` produces the response.
    ///
    /// Returns `None` if no items pass the filter or `X^T X` is singular.
    pub fn fit<T: Sync, const D: usize>(
        items: &[T],
        filter: impl Fn(&T) -> bool + Sync,
        embed: impl Fn(&T) -> [f64; D] + Sync,
        target: impl Fn(&T) -> f64 + Sync,
    ) -> Option<Self> {
        let acc = items
            .par_iter()
            .filter(|x| filter(x))
            .fold(
                || Acc::zero(D),
                |mut acc, x| {
                    let row = embed(x);
                    acc.add_row(&row, target(x));
                    acc
                },
            )
            .reduce(|| Acc::zero(D), Acc::merge);

        if acc.n == 0 {
            return None;
        }

        let nf = acc.n as f64;
        let y_mean = acc.sum_y / nf;
        let y_var = acc.sum_y2 - nf * y_mean * y_mean;

        let cov = Matrix::new(acc.cov, D, D);
        let b_mat = Matrix::col_vector(acc.b);
        let beta = Gauss::solve(cov, b_mat)?.take();

        // Streaming pass for SSE = sum((X beta - y)^2). O(N*D), small vs fit.
        let sse: f64 = items
            .par_iter()
            .filter(|x| filter(x))
            .map(|x| {
                let row = embed(x);
                let pred: f64 = row.iter().zip(&beta).map(|(v, w)| v * w).sum();
                let act = target(x);
                (pred - act).powi(2)
            })
            .sum();

        let r2 = 1.0 - sse / y_var;
        Some(Self { beta, r2 })
    }
}

#[cfg(test)]
#[path = "../../tests/unit/ml/regression.rs"]
mod test;
