//! One-dimensional signal helpers: apex climbing, peak bounds, smoothing.
//!
//! Moved out of `lfq.rs` unchanged; every item keeps its `lfq::` path.

use super::*;

/// Fraction of the summed apex height at which apex-mode integration stops,
/// unless a valley comes first. Half height matches the window search; 20%
/// took in more interference and worsened the E. coli ratio on PXD028735.
pub(super) const APEX_BOUND_FRACTION: f64 = 0.5;

/// Climb from `seed` to the apex of the peak it sits on. The seed snaps to
/// the nearest signal within [`APEX_SNAP_BINS`]; the climb then moves to the
/// highest point reachable without the trace falling below half of the
/// current apex, so that small dips in a smoothed peak do not stop it early.
pub(super) fn climb_to_apex(trace: &[f64], seed: usize) -> Option<usize> {
    let n = trace.len();
    let mut apex = (0..=APEX_SNAP_BINS)
        .flat_map(|d| [seed.checked_sub(d), Some(seed + d).filter(|&i| i < n)])
        .flatten()
        .find(|&i| trace[i] > 0.0)?;
    loop {
        let height = trace[apex];
        let mut next = None;
        for dir in [-1isize, 1] {
            let mut i = apex as isize + dir;
            while i >= 0 && (i as usize) < n && trace[i as usize] >= 0.5 * height {
                if trace[i as usize] > trace[next.unwrap_or(apex)] {
                    next = Some(i as usize);
                }
                i += dir;
            }
        }
        match next {
            Some(higher) => apex = higher,
            None => return Some(apex),
        }
    }
}

/// Inclusive integration bounds `[left, right]` around `apex`: extend while
/// the trace stays above `fraction` of the apex intensity and keeps falling (a
/// valley ends the peak). Inclusive to match [`PeakWindow::right`], which every
/// consumer slices with `..=`.
pub(super) fn peak_bounds(trace: &[f64], apex: usize, fraction: f64) -> (usize, usize) {
    let floor = trace[apex] * fraction;
    let mut left = apex;
    while left > 0 && trace[left - 1] >= floor && trace[left - 1] <= trace[left] {
        left -= 1;
    }
    let mut right = apex;
    while right + 1 < trace.len() && trace[right + 1] >= floor && trace[right + 1] <= trace[right] {
        right += 1;
    }
    (left, right)
}

/// Width, in (fractional) bins, of `trace` at half its height at `apex`.
pub(super) fn half_maximum_width(trace: &[f64], apex: usize) -> Option<f64> {
    let half = trace[apex] / 2.0;
    if half <= 0.0 {
        return None;
    }
    let crossing = |mut i: usize, step: isize| -> Option<f64> {
        loop {
            let j = i as isize + step;
            if j < 0 || j as usize >= trace.len() {
                return None;
            }
            let j = j as usize;
            if trace[j] < half {
                // Interpolate between i (above) and j (below).
                let frac = (trace[i] - half) / (trace[i] - trace[j]);
                return Some(i as f64 + step as f64 * frac);
            }
            i = j;
        }
    };
    Some(crossing(apex, 1)? - crossing(apex, -1)?)
}

/// Create a symmetrical gaussian kernel of given standard deviation and length
pub(super) fn gaussian_kernel(sigma: f64, len: usize) -> Vec<f64> {
    let step = 2.0 / (len - 1) as f64;
    let constant = 1.0 / (sigma * (2.0 * std::f64::consts::PI).sqrt());

    let mut kernel = (0..len)
        .map(|i| {
            let x = i as f64 * step - 1.0;
            constant * (-0.5 * (x / sigma).powi(2)).exp()
        })
        .collect::<Vec<_>>();

    let sum = kernel.iter().sum::<f64>();
    kernel.iter_mut().for_each(|x| *x /= sum);
    kernel
}

/// Convolve a signal with a symmetrical kernel
/// - This should behave the same as `np.convolve(..., mode='same')`
pub(super) fn convolve(slice: &[f64], kernel: &[f64]) -> Vec<f64> {
    // Middle index of kernel
    let n = kernel.len() - (kernel.len() / 2);

    (0..slice.len())
        .map(|idx| {
            // If idx < kernel.len(), take only some subset of the kernel (at least half)
            let k = &kernel[kernel.len().saturating_sub(n + idx)..];
            // If idx < kernel.len(), then start from 0
            let w = &slice[idx.saturating_sub(n - 1)..];
            // Dot product
            w.iter().zip(k).fold(0.0, |acc, (x, y)| acc + x * y)
        })
        .collect()
}
