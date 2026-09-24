//! Mass-error residuals as a mixture of a location-scale signal and a
//! uniform background over the searched window.
//!
//! Confident PSMs and the fragment ions they match contain wrong matches that
//! land anywhere in the searched window. Their residuals follow
//! `pi * f(x) + (1 - pi) / (hi - lo)` on the window `[lo, hi]`, where `f` is
//! the error distribution of correct matches. Fitting this by EM separates the
//! two, so a tolerance can be sized from `f` alone instead of from all
//! residuals, which fill whatever window was searched.
//!
//! Three signal forms are fitted: a Gaussian, a Student-t with four degrees of
//! freedom, and two Gaussians with a shared center. Held-out log-likelihood
//! chooses between them; a heavier form replaces a simpler one only if it
//! gains at least `min_gain_per_point` nats per held-out residual.

use serde::{Deserialize, Serialize};

/// Shape of the signal part of the mixture.
#[derive(Copy, Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum SignalForm {
    Gaussian,
    /// Student-t with four degrees of freedom.
    StudentT4,
    /// Two Gaussians with a shared center.
    TwoGaussian,
}

/// Settings for [`fit_signal_mixture`].
#[derive(Copy, Clone, Debug, PartialEq)]
pub struct MixtureOptions {
    /// Two-sided signal mass inside the reported half-width.
    pub signal_quantile: f64,
    pub max_iterations: usize,
    /// Relative log-likelihood change that ends EM.
    pub tolerance: f64,
    /// Held-out log-likelihood gain per point a heavier form needs.
    pub min_gain_per_point: f64,
    /// Smallest allowed scale, in ppm.
    pub min_scale_ppm: f64,
    /// Residuals used per fit (and per held-out evaluation); larger sets are
    /// thinned evenly.
    pub max_points: usize,
}

impl Default for MixtureOptions {
    fn default() -> Self {
        Self {
            signal_quantile: 0.99,
            max_iterations: 500,
            tolerance: 1e-8,
            min_gain_per_point: 0.005,
            min_scale_ppm: 0.02,
            max_points: 40_000,
        }
    }
}

/// Held-out score of one fitted form.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct FormScore {
    pub form: SignalForm,
    pub signal_fraction: f32,
    pub center_ppm: f32,
    pub scale_ppm: f32,
    pub converged: bool,
    /// Mean held-out log-likelihood per residual.
    pub validation_log_likelihood: f32,
}

/// A fitted signal + uniform mixture.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct SignalFit {
    pub form: SignalForm,
    /// Mixture weight of the signal part (`pi`).
    pub signal_fraction: f32,
    pub center_ppm: f32,
    /// Gaussian sigma, Student-t scale, or the narrow Gaussian's sigma.
    pub scale_ppm: f32,
    /// Sigma of the wide Gaussian (two-Gaussian form only).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub wide_scale_ppm: Option<f32>,
    /// Weight of the narrow Gaussian within the signal (two-Gaussian form).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub narrow_weight: Option<f32>,
    /// Half-width around `center_ppm` holding `signal_quantile` of the signal.
    pub half_width_ppm: f32,
    pub signal_quantile: f32,
    pub iterations: usize,
    pub converged: bool,
    /// A scale ended on its lower or upper bound.
    pub scale_at_bound: bool,
    /// Searched window the background is uniform over.
    pub window_ppm: (f32, f32),
    pub fit_points: usize,
    pub validation_points: usize,
    /// Mean held-out log-likelihood per residual.
    pub validation_log_likelihood: f32,
    /// Every form tried, in the order considered.
    pub candidates: Vec<FormScore>,
    #[serde(skip)]
    params: Params,
}

#[derive(Copy, Clone, Debug, Default, PartialEq)]
struct Params {
    form: Option<SignalForm>,
    pi: f64,
    mu: f64,
    s1: f64,
    s2: f64,
    w: f64,
    log_uniform: f64,
}

const LN_SQRT_2PI: f64 = 0.918_938_533_204_672_8;
/// `ln Gamma(5/2) - ln Gamma(2) - ln(sqrt(4 pi))`.
const LN_T4_NORM: f64 = -0.980_829_253_011_726_2;

fn ln_gauss(x: f64, mu: f64, s: f64) -> f64 {
    let d = (x - mu) / s;
    -0.5 * d * d - s.ln() - LN_SQRT_2PI
}

fn ln_add(a: f64, b: f64) -> f64 {
    let m = a.max(b);
    if m == f64::NEG_INFINITY {
        return m;
    }
    m + ((a - m).exp() + (b - m).exp()).ln()
}

impl Params {
    fn ln_signal(&self, x: f64) -> f64 {
        match self.form.expect("fitted form") {
            SignalForm::Gaussian => ln_gauss(x, self.mu, self.s1),
            SignalForm::StudentT4 => {
                let d = (x - self.mu) / self.s1;
                LN_T4_NORM - self.s1.ln() - 2.5 * (d * d / 4.0).ln_1p()
            }
            SignalForm::TwoGaussian => ln_add(
                self.w.ln() + ln_gauss(x, self.mu, self.s1),
                (1.0 - self.w).ln() + ln_gauss(x, self.mu, self.s2),
            ),
        }
    }

    /// Log mixture density and the signal responsibility.
    fn evaluate(&self, x: f64) -> (f64, f64) {
        let signal = self.pi.ln() + self.ln_signal(x);
        let background = (1.0 - self.pi).ln() + self.log_uniform;
        let total = ln_add(signal, background);
        (total, (signal - total).exp())
    }

    fn mean_log_likelihood(&self, xs: &[f64]) -> f64 {
        if xs.is_empty() {
            return f64::NAN;
        }
        xs.iter().map(|&x| self.evaluate(x).0).sum::<f64>() / xs.len() as f64
    }

    /// Two-sided signal tail mass beyond `h` from the center.
    fn tail(&self, h: f64) -> f64 {
        let gauss_tail = |s: f64| erfc(h / (s * std::f64::consts::SQRT_2));
        match self.form.expect("fitted form") {
            SignalForm::Gaussian => gauss_tail(self.s1),
            SignalForm::StudentT4 => {
                // For nu = 4, F(t) = 1/2 + sin(a) (1 + cos(a)^2 / 2) / 2 with
                // a = atan(t / 2).
                let t = h / self.s1;
                let sin = t / (t * t + 4.0).sqrt();
                let cos2 = 1.0 - sin * sin;
                (1.0 - sin * (1.0 + 0.5 * cos2)).max(0.0)
            }
            SignalForm::TwoGaussian => {
                self.w * gauss_tail(self.s1) + (1.0 - self.w) * gauss_tail(self.s2)
            }
        }
    }

    fn half_width(&self, quantile: f64) -> f64 {
        let target = 1.0 - quantile;
        let (mut lo, mut hi) = (0.0, 1.0);
        while self.tail(hi) > target && hi < 1e7 {
            hi *= 2.0;
        }
        for _ in 0..100 {
            let mid = 0.5 * (lo + hi);
            if self.tail(mid) > target {
                lo = mid;
            } else {
                hi = mid;
            }
        }
        hi
    }
}

/// Complementary error function (Numerical Recipes `erfcc`, fractional error
/// below 1.2e-7).
fn erfc(x: f64) -> f64 {
    let z = x.abs();
    let t = 1.0 / (1.0 + 0.5 * z);
    let poly = -z * z - 1.265_512_23
        + t * (1.000_023_68
            + t * (0.374_091_96
                + t * (0.096_784_18
                    + t * (-0.186_288_06
                        + t * (0.278_868_07
                            + t * (-1.135_203_98
                                + t * (1.488_515_87 + t * (-0.822_152_23 + t * 0.170_872_77))))))));
    let r = t * poly.exp();
    if x >= 0.0 {
        r
    } else {
        2.0 - r
    }
}

fn median(values: &mut [f64]) -> f64 {
    let middle = values.len() / 2;
    let (_, m, _) = values.select_nth_unstable_by(middle, f64::total_cmp);
    *m
}

fn thin(values: &[f32], max: usize) -> Vec<f64> {
    let stride = values.len().div_ceil(max.max(1)).max(1);
    values
        .iter()
        .step_by(stride)
        .map(|&v| v as f64)
        .filter(|v| v.is_finite())
        .collect()
}

struct Fitted {
    params: Params,
    iterations: usize,
    converged: bool,
    at_bound: bool,
}

fn fit_form(xs: &[f64], form: SignalForm, lo: f64, hi: f64, options: &MixtureOptions) -> Fitted {
    let mut sorted = xs.to_vec();
    let mu0 = median(&mut sorted);
    let mut deviations = xs.iter().map(|x| (x - mu0).abs()).collect::<Vec<_>>();
    let max_scale = 0.5 * (hi - lo);
    let min_scale = options.min_scale_ppm;
    let sigma0 = (1.4826 * median(&mut deviations)).clamp(min_scale, max_scale);
    let mut p = Params {
        form: Some(form),
        pi: 0.8,
        mu: mu0,
        s1: sigma0,
        s2: (3.0 * sigma0).min(max_scale),
        w: 0.6,
        log_uniform: -(hi - lo).ln(),
    };
    if form == SignalForm::TwoGaussian {
        p.s1 = (0.7 * sigma0).max(min_scale);
    }
    let n = xs.len() as f64;
    let mut previous = f64::NEG_INFINITY;
    let mut converged = false;
    let mut iterations = 0;
    let mut at_bound = false;
    while iterations < options.max_iterations {
        iterations += 1;
        // E step, accumulating the sufficient statistics of each form.
        let mut ll = 0.0;
        let (mut g_sum, mut gu_sum, mut gux_sum) = (0.0, 0.0, 0.0);
        let (mut r1_sum, mut r1x, mut r2_sum, mut r2x) = (0.0, 0.0, 0.0, 0.0);
        let mut responsibilities = Vec::with_capacity(xs.len());
        for &x in xs {
            let (total, g) = p.evaluate(x);
            ll += total;
            g_sum += g;
            match form {
                SignalForm::Gaussian => {
                    gux_sum += g * x;
                    responsibilities.push((g, 1.0));
                }
                SignalForm::StudentT4 => {
                    let d = (x - p.mu) / p.s1;
                    let u = 5.0 / (4.0 + d * d);
                    gu_sum += g * u;
                    gux_sum += g * u * x;
                    responsibilities.push((g, u));
                }
                SignalForm::TwoGaussian => {
                    let a = p.w.ln() + ln_gauss(x, p.mu, p.s1);
                    let b = (1.0 - p.w).ln() + ln_gauss(x, p.mu, p.s2);
                    let h = (a - ln_add(a, b)).exp();
                    let (r1, r2) = (g * h, g * (1.0 - h));
                    r1_sum += r1;
                    r1x += r1 * x;
                    r2_sum += r2;
                    r2x += r2 * x;
                    responsibilities.push((r1, r2));
                }
            }
        }
        if g_sum <= 0.0 {
            break;
        }
        // M step.
        p.pi = (g_sum / n).clamp(1e-4, 1.0 - 1e-9);
        let clamp = |s: f64, at_bound: &mut bool| {
            let c = s.clamp(min_scale, max_scale);
            *at_bound = c != s || !s.is_finite();
            if s.is_finite() {
                c
            } else {
                max_scale
            }
        };
        match form {
            SignalForm::Gaussian => {
                p.mu = gux_sum / g_sum;
                let ss = xs
                    .iter()
                    .zip(&responsibilities)
                    .map(|(x, (g, _))| g * (x - p.mu).powi(2))
                    .sum::<f64>();
                p.s1 = clamp((ss / g_sum).sqrt(), &mut at_bound);
            }
            SignalForm::StudentT4 => {
                p.mu = gux_sum / gu_sum;
                let ss = xs
                    .iter()
                    .zip(&responsibilities)
                    .map(|(x, (g, u))| g * u * (x - p.mu).powi(2))
                    .sum::<f64>();
                p.s1 = clamp((ss / g_sum).sqrt(), &mut at_bound);
            }
            SignalForm::TwoGaussian => {
                let (a1, a2) = (r1_sum / p.s1.powi(2), r2_sum / p.s2.powi(2));
                p.mu = (r1x / p.s1.powi(2) + r2x / p.s2.powi(2)) / (a1 + a2);
                let (ss1, ss2) =
                    xs.iter()
                        .zip(&responsibilities)
                        .fold((0.0, 0.0), |(s1, s2), (x, (r1, r2))| {
                            let d2 = (x - p.mu).powi(2);
                            (s1 + r1 * d2, s2 + r2 * d2)
                        });
                let mut bound1 = false;
                let mut bound2 = false;
                p.s1 = clamp((ss1 / r1_sum.max(1e-300)).sqrt(), &mut bound1);
                p.s2 = clamp((ss2 / r2_sum.max(1e-300)).sqrt(), &mut bound2);
                at_bound = bound1 || bound2;
                p.w = (r1_sum / g_sum).clamp(0.01, 0.99);
                if p.s1 > p.s2 {
                    std::mem::swap(&mut p.s1, &mut p.s2);
                    p.w = 1.0 - p.w;
                }
            }
        }
        if (ll - previous).abs() <= options.tolerance * ll.abs() {
            converged = true;
            break;
        }
        previous = ll;
    }
    Fitted {
        params: p,
        iterations,
        converged,
        at_bound,
    }
}

/// Fit the mixture on `fit` residuals over the searched window `[lo, hi]`
/// and choose a signal form by held-out log-likelihood on `validation`.
/// Returns `None` for an empty fit set or an empty window.
pub fn fit_signal_mixture(
    fit: &[f32],
    validation: &[f32],
    lo: f32,
    hi: f32,
    options: &MixtureOptions,
) -> Option<SignalFit> {
    let (lo, hi) = (lo as f64, hi as f64);
    let fit_xs = thin(fit, options.max_points);
    let validation_xs = thin(validation, options.max_points);
    if fit_xs.len() < 10 || hi.partial_cmp(&lo) != Some(std::cmp::Ordering::Greater) {
        return None;
    }
    let score_on = if validation_xs.is_empty() {
        &fit_xs
    } else {
        &validation_xs
    };
    let mut candidates = Vec::new();
    let mut best: Option<(Fitted, f64)> = None;
    for form in [
        SignalForm::Gaussian,
        SignalForm::StudentT4,
        SignalForm::TwoGaussian,
    ] {
        let fitted = fit_form(&fit_xs, form, lo, hi, options);
        let ll = fitted.params.mean_log_likelihood(score_on);
        candidates.push(FormScore {
            form,
            signal_fraction: fitted.params.pi as f32,
            center_ppm: fitted.params.mu as f32,
            scale_ppm: fitted.params.s1 as f32,
            converged: fitted.converged,
            validation_log_likelihood: ll as f32,
        });
        let better = match &best {
            None => true,
            Some((_, best_ll)) => {
                fitted.converged && ll.is_finite() && ll - best_ll >= options.min_gain_per_point
            }
        };
        if better {
            best = Some((fitted, ll));
        }
    }
    let (fitted, ll) = best?;
    let p = fitted.params;
    let two = p.form == Some(SignalForm::TwoGaussian);
    Some(SignalFit {
        form: p.form.expect("fitted form"),
        signal_fraction: p.pi as f32,
        center_ppm: p.mu as f32,
        scale_ppm: p.s1 as f32,
        wide_scale_ppm: two.then_some(p.s2 as f32),
        narrow_weight: two.then_some(p.w as f32),
        half_width_ppm: p.half_width(options.signal_quantile) as f32,
        signal_quantile: options.signal_quantile as f32,
        iterations: fitted.iterations,
        converged: fitted.converged,
        scale_at_bound: fitted.at_bound,
        window_ppm: (lo as f32, hi as f32),
        fit_points: fit_xs.len(),
        validation_points: validation_xs.len(),
        validation_log_likelihood: ll as f32,
        candidates,
        params: p,
    })
}

impl SignalFit {
    /// Posterior probability that residual `x` belongs to the signal.
    pub fn signal_probability(&self, x: f32) -> f32 {
        self.params.evaluate(x as f64).1 as f32
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn noise(i: u64) -> f64 {
        let mut x = i.wrapping_mul(0x9e3779b97f4a7c15);
        x ^= x >> 31;
        x = x.wrapping_mul(0xbf58476d1ce4e5b9);
        x ^= x >> 29;
        (x % 1_000_003) as f64 / 1_000_003.0
    }

    /// Standard normal draw by Box-Muller.
    fn normal(i: u64) -> f64 {
        let u1 = noise(2 * i + 1).max(1e-12);
        let u2 = noise(2 * i + 2);
        (-2.0 * u1.ln()).sqrt() * (2.0 * std::f64::consts::PI * u2).cos()
    }

    fn sample(n: u64, pi: f64, mu: f64, sigma: f64, lo: f64, hi: f64) -> Vec<f32> {
        (0..n)
            .map(|i| {
                if noise(i + 5_000_000) < pi {
                    (mu + sigma * normal(i)) as f32
                } else {
                    (lo + (hi - lo) * noise(i + 9_000_000)) as f32
                }
            })
            .collect()
    }

    #[test]
    fn erfc_and_t4_tails_match_tables() {
        assert!((erfc(0.0) - 1.0).abs() < 1e-6);
        assert!((erfc(2.5758 / std::f64::consts::SQRT_2) - 0.01).abs() < 1e-5);
        assert!((erfc(-1.0) - 1.842_700_79).abs() < 1e-6);
        let t4 = Params {
            form: Some(SignalForm::StudentT4),
            s1: 1.0,
            ..Params::default()
        };
        // t(0.995, 4) = 4.6041.
        assert!((t4.half_width(0.99) - 4.6041).abs() < 1e-3);
        let gauss = Params {
            form: Some(SignalForm::Gaussian),
            s1: 2.0,
            ..Params::default()
        };
        assert!((gauss.half_width(0.99) - 2.0 * 2.5758).abs() < 1e-3);
    }

    #[test]
    fn separates_gaussian_signal_from_uniform_background() {
        let data = sample(20_000, 0.6, 1.0, 2.0, -50.0, 50.0);
        let (fit, validation): (Vec<_>, Vec<_>) =
            data.iter().enumerate().partition(|(i, _)| i % 10 >= 3);
        let fit = fit.into_iter().map(|(_, &x)| x).collect::<Vec<_>>();
        let validation = validation.into_iter().map(|(_, &x)| x).collect::<Vec<_>>();
        let result = fit_signal_mixture(&fit, &validation, -50.0, 50.0, &MixtureOptions::default())
            .expect("fit");
        assert_eq!(result.form, SignalForm::Gaussian);
        assert!(result.converged);
        assert!((result.signal_fraction - 0.6).abs() < 0.02);
        assert!((result.center_ppm - 1.0).abs() < 0.1);
        assert!((result.scale_ppm - 2.0).abs() < 0.1);
        assert!((result.half_width_ppm - 2.5758 * result.scale_ppm).abs() < 1e-3);
        assert!(result.signal_probability(1.0) > 0.95);
        assert!(result.signal_probability(40.0) < 1e-6);
    }

    #[test]
    fn heavy_tailed_signal_prefers_a_heavier_form() {
        // Mixture of a 1 ppm and a 4 ppm Gaussian sharing a center, over a
        // sparse background.
        let narrow = sample(15_000, 0.97, 0.0, 1.0, -20.0, 20.0);
        let wide = sample(5_000, 0.97, 0.0, 4.0, -20.0, 20.0);
        let data = narrow.into_iter().chain(wide).collect::<Vec<_>>();
        let fit = data.iter().step_by(2).copied().collect::<Vec<_>>();
        let validation = data.iter().skip(1).step_by(2).copied().collect::<Vec<_>>();
        let result = fit_signal_mixture(&fit, &validation, -20.0, 20.0, &MixtureOptions::default())
            .expect("fit");
        assert_ne!(result.form, SignalForm::Gaussian);
        // The true 99% half-width of 0.75 N(0,1) + 0.25 N(0,4) is ~8.5 ppm.
        assert!(
            result.half_width_ppm > 7.0 && result.half_width_ppm < 11.0,
            "{result:?}"
        );
    }
}
