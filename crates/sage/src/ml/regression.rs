//! Streaming linear regression on standardized features.
//!
//! Fits `y ~ X beta` without materializing the `n x D` design matrix. The
//! sequence embeddings used by the retention-time and ion-mobility models are
//! exactly collinear (residue counts sum to the peptide length, terminal
//! one-hots sum to a multiple of the intercept, residue-class counts are sums
//! of residue counts, ...), so `X^T X` is singular and plain OLS has no unique
//! solution. Instead of perturbing the raw normal equations until elimination
//! "succeeds", we solve a well-posed ridge problem on standardized columns:
//!
//! 1. Column means `mu`, `mean(y)` (pass 1) and the centered cross products
//!    `Xc^T Xc`, `Xc^T yc` (pass 2).
//! 2. Columns with zero variance (the intercept, or residues that never occur)
//!    are removed from the system; the remaining columns are scaled to unit
//!    variance, `Z = Xc / sigma`, so `diag(Z^T Z) = n`.
//! 3. `(Z^T Z + lambda I) b = Z^T yc` with `lambda = RIDGE_PER_ROW * n`, solved
//!    by Cholesky factorization followed by a fixed number of iterative
//!    refinement steps against the unregularized equations. Along identifiable
//!    directions this converges to OLS; along exactly collinear directions the
//!    residual is pure rounding, so the coefficients stay at the (tiny)
//!    minimum-norm ridge solution. Fitted values therefore do not depend on
//!    how redundant columns are laid out.
//! 4. Back-transform: `beta_j = b_j / sigma_j`, and the intercept
//!    `mean(y) - sum_j beta_j mu_j` is placed on the first constant column.
//!    Without a constant column the uncentered (through-origin) problem is
//!    solved with the same scaling.
//!
//! A final pass evaluates `sum((X beta - y)^2)` for r^2. Every pass sums
//! fixed-size chunks sequentially and combines the chunk partials in index
//! order, so results are bitwise reproducible regardless of thread count or
//! scheduling. Scratch is `O(D^2)` per chunk.

use rayon::prelude::*;

/// Rows per sequentially-summed chunk. Fixed (not derived from the thread
/// count) so the floating-point summation order never changes.
const CHUNK_SIZE: usize = 4096;

/// Ridge penalty per training row on the standardized scale, i.e. relative to
/// the unit diagonal of the correlation matrix `Z^T Z / n`.
///
/// Its job is to make the exactly singular directions of the peptide
/// embeddings well posed: it must dominate their rounding noise (about
/// `D * f64::EPSILON ~ 1e-14` relative) so the Cholesky factorization is
/// stable and their coefficients stay at the small minimum-norm solution. The
/// ridge bias on identifiable directions (smallest genuine eigenvalues of the
/// embeddings are ~1e-4, e.g. monoisotopic mass versus residue counts) is then
/// removed by `REFINEMENT_STEPS` of iterative refinement, so the fit equals
/// OLS on every identifiable direction to near machine precision. On
/// synthetic basic-mobility data predictions agree to 1e-10 for penalties
/// between 1e-12 and 1e-6 and degrade only from 1e-5 upwards.
pub const RIDGE_PER_ROW: f64 = 1e-8;

/// Iterative-refinement steps applied after the ridge solve (see
/// `solve_standardized`). Each step multiplies the remaining bias along an
/// eigendirection with eigenvalue `e` by `lambda / (e + lambda)`.
const REFINEMENT_STEPS: usize = 3;

/// Relative standard deviation below which a column is treated as constant.
const CONSTANT_COLUMN_TOLERANCE: f64 = 1e-12;

pub struct LinearRegression {
    pub beta: Vec<f64>,
    pub r2: f64,
}

/// Pass 1: row count and column sums.
struct Sums {
    x: Vec<f64>,
    y: f64,
    n: usize,
}

impl Sums {
    fn zero(d: usize) -> Self {
        Self {
            x: vec![0.0; d],
            y: 0.0,
            n: 0,
        }
    }

    fn merge(mut self, other: Sums) -> Sums {
        for (a, b) in self.x.iter_mut().zip(other.x) {
            *a += b;
        }
        self.y += other.y;
        self.n += other.n;
        self
    }
}

/// Pass 2: centered cross products.
struct Acc {
    cov: Vec<f64>, // D*D row-major, sum (x - mu)(x - mu)^T
    b: Vec<f64>,   // D, sum (x - mu)(y - mean_y)
    syy: f64,      // sum (y - mean_y)^2
}

impl Acc {
    fn zero(d: usize) -> Self {
        Self {
            cov: vec![0.0; d * d],
            b: vec![0.0; d],
            syy: 0.0,
        }
    }

    fn add_row(&mut self, row: &[f64], y: f64) {
        let d = row.len();
        for j in 0..d {
            let rj = row[j];
            self.b[j] += rj * y;
            let off = j * d;
            for (k, &rk) in row.iter().enumerate() {
                self.cov[off + k] += rj * rk;
            }
        }
        self.syy += y * y;
    }

    fn merge(mut self, other: Acc) -> Acc {
        for i in 0..self.cov.len() {
            self.cov[i] += other.cov[i];
        }
        for i in 0..self.b.len() {
            self.b[i] += other.b[i];
        }
        self.syy += other.syy;
        self
    }
}

impl LinearRegression {
    /// Fit a linear model over `items` with predicate `filter`. `embed(item)`
    /// produces a design row of length `D`; `target(item)` produces the
    /// response. Collinear and constant columns are allowed (see module docs).
    ///
    /// Returns `None` if no items pass the filter or the fit is not finite.
    pub fn fit<T: Sync, const D: usize>(
        items: &[T],
        filter: impl Fn(&T) -> bool + Sync,
        embed: impl Fn(&T) -> [f64; D] + Sync,
        target: impl Fn(&T) -> f64 + Sync,
    ) -> Option<Self> {
        Self::fit_with_ridge(items, filter, embed, target, RIDGE_PER_ROW)
    }

    /// [`LinearRegression::fit`] with an explicit standardized ridge penalty per
    /// training row.
    pub fn fit_with_ridge<T: Sync, const D: usize>(
        items: &[T],
        filter: impl Fn(&T) -> bool + Sync,
        embed: impl Fn(&T) -> [f64; D] + Sync,
        target: impl Fn(&T) -> f64 + Sync,
        ridge_per_row: f64,
    ) -> Option<Self> {
        let sums = items
            .par_chunks(CHUNK_SIZE)
            .map(|chunk| {
                let mut sums = Sums::zero(D);
                for x in chunk.iter().filter(|x| filter(x)) {
                    for (s, v) in sums.x.iter_mut().zip(embed(x)) {
                        *s += v;
                    }
                    sums.y += target(x);
                    sums.n += 1;
                }
                sums
            })
            .collect::<Vec<_>>()
            .into_iter()
            .fold(Sums::zero(D), Sums::merge);

        if sums.n == 0 {
            return None;
        }
        let nf = sums.n as f64;
        let mu = sums.x.iter().map(|s| s / nf).collect::<Vec<_>>();
        let y_mean = sums.y / nf;

        let acc = items
            .par_chunks(CHUNK_SIZE)
            .map(|chunk| {
                let mut acc = Acc::zero(D);
                let mut centered = [0.0; D];
                for x in chunk.iter().filter(|x| filter(x)) {
                    for ((c, v), m) in centered.iter_mut().zip(embed(x)).zip(&mu) {
                        *c = v - m;
                    }
                    acc.add_row(&centered, target(x) - y_mean);
                }
                acc
            })
            .collect::<Vec<_>>()
            .into_iter()
            .fold(Acc::zero(D), Acc::merge);

        let beta = solve_standardized(acc.cov, acc.b, &mu, y_mean, nf, ridge_per_row)?;

        // Streaming pass for SSE = sum((X beta - y)^2). O(N*D), small vs fit.
        let sse: f64 = items
            .par_chunks(CHUNK_SIZE)
            .map(|chunk| {
                chunk
                    .iter()
                    .filter(|x| filter(x))
                    .map(|x| {
                        let row = embed(x);
                        let pred: f64 = row.iter().zip(&beta).map(|(v, w)| v * w).sum();
                        let act = target(x);
                        (pred - act).powi(2)
                    })
                    .sum::<f64>()
            })
            .collect::<Vec<_>>()
            .into_iter()
            .sum();

        let r2 = 1.0 - sse / acc.syy;
        Some(Self { beta, r2 })
    }
}

/// Solve the standardized ridge problem from centered moments and map the
/// coefficients back to the original columns.
fn solve_standardized(
    mut cov: Vec<f64>,
    mut b: Vec<f64>,
    mu: &[f64],
    y_mean: f64,
    n: f64,
    ridge_per_row: f64,
) -> Option<Vec<f64>> {
    let d = mu.len();
    let is_constant = |j: usize, second_moment: f64| {
        let sd = (second_moment.max(0.0) / n).sqrt();
        sd <= CONSTANT_COLUMN_TOLERANCE * mu[j].abs().max(1.0)
    };
    let intercept = (0..d).find(|&j| mu[j] != 0.0 && is_constant(j, cov[j * d + j]));

    if intercept.is_none() {
        // No constant column: the model goes through the origin, so undo the
        // centering and work with raw second moments.
        for j in 0..d {
            b[j] += n * mu[j] * y_mean;
            for k in 0..d {
                cov[j * d + k] += n * mu[j] * mu[k];
            }
        }
    }
    let active = (0..d)
        .filter(|&j| {
            let moment = cov[j * d + j];
            moment > 0.0
                && if intercept.is_some() {
                    !is_constant(j, moment)
                } else {
                    (moment / n).sqrt() > CONSTANT_COLUMN_TOLERANCE
                }
        })
        .collect::<Vec<_>>();
    let scale = active
        .iter()
        .map(|&j| (cov[j * d + j] / n).sqrt())
        .collect::<Vec<_>>();

    let m = active.len();
    let mut gram = vec![0.0; m * m];
    let mut rhs = vec![0.0; m];
    for (p, (&j, &sj)) in active.iter().zip(&scale).enumerate() {
        rhs[p] = b[j] / sj;
        for (q, (&k, &sk)) in active.iter().zip(&scale).enumerate() {
            gram[p * m + q] = cov[j * d + k] / (sj * sk);
        }
    }
    let mut regularized = gram.clone();
    for p in 0..m {
        regularized[p * m + p] += ridge_per_row * n;
    }
    let factor = cholesky_factor(regularized)?;
    let mut standardized = rhs.clone();
    cholesky_substitute(&factor, &mut standardized);
    // Iterated Tikhonov refinement against the unregularized normal
    // equations: each step shrinks the ridge bias along an eigendirection with
    // eigenvalue `e` by `lambda / (e + lambda)`, while exactly collinear
    // directions (where the residual is pure rounding) stay at the ridge
    // (minimum-norm) solution.
    for _ in 0..REFINEMENT_STEPS {
        let mut residual = (0..m)
            .map(|p| {
                rhs[p]
                    - gram[p * m..(p + 1) * m]
                        .iter()
                        .zip(&standardized)
                        .map(|(g, x)| g * x)
                        .sum::<f64>()
            })
            .collect::<Vec<_>>();
        cholesky_substitute(&factor, &mut residual);
        for (x, dx) in standardized.iter_mut().zip(residual) {
            *x += dx;
        }
    }

    let mut beta = vec![0.0; d];
    for ((&j, &sj), coef) in active.iter().zip(&scale).zip(standardized) {
        beta[j] = coef / sj;
    }
    if let Some(c) = intercept {
        let offset = y_mean - (0..d).map(|j| beta[j] * mu[j]).sum::<f64>();
        beta[c] = offset / mu[c];
    }
    beta.iter().all(|b| b.is_finite()).then_some(beta)
}

/// Solve `A x = b` for a symmetric positive-definite `A` (row-major, `n x n`)
/// by Cholesky factorization. Deterministic; returns `None` if `A` is not
/// numerically positive definite or the inputs are not finite.
pub fn cholesky_solve(a: Vec<f64>, mut b: Vec<f64>) -> Option<Vec<f64>> {
    assert_eq!(a.len(), b.len() * b.len(), "cholesky_solve: shape mismatch");
    let factor = cholesky_factor(a)?;
    cholesky_substitute(&factor, &mut b);
    b.iter().all(|x| x.is_finite()).then_some(b)
}

/// In-place lower-triangular factor `L` with `A = L L^T` (upper triangle is
/// left untouched and ignored).
fn cholesky_factor(mut a: Vec<f64>) -> Option<Vec<f64>> {
    let n = (a.len() as f64).sqrt() as usize;
    debug_assert_eq!(n * n, a.len());
    for j in 0..n {
        let mut diag = a[j * n + j];
        for k in 0..j {
            diag -= a[j * n + k] * a[j * n + k];
        }
        if !(diag > 0.0 && diag.is_finite()) {
            return None;
        }
        let diag = diag.sqrt();
        a[j * n + j] = diag;
        for i in j + 1..n {
            let mut value = a[i * n + j];
            for k in 0..j {
                value -= a[i * n + k] * a[j * n + k];
            }
            a[i * n + j] = value / diag;
        }
    }
    Some(a)
}

/// Overwrite `b` with `(L L^T)^-1 b`.
fn cholesky_substitute(l: &[f64], b: &mut [f64]) {
    let n = b.len();
    for i in 0..n {
        let mut value = b[i];
        for k in 0..i {
            value -= l[i * n + k] * b[k];
        }
        b[i] = value / l[i * n + i];
    }
    for i in (0..n).rev() {
        let mut value = b[i];
        for k in i + 1..n {
            value -= l[k * n + i] * b[k];
        }
        b[i] = value / l[i * n + i];
    }
}

#[cfg(test)]
#[path = "../../tests/unit/ml/regression.rs"]
mod test;
