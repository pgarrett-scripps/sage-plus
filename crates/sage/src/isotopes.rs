use std::sync::OnceLock;

// Manually unroll the convolution loop
fn convolve(a: &[f32; 4], b: &[f32; 4]) -> [f32; 4] {
    [
        a[0] * b[0],
        a[0] * b[1] + a[1] * b[0],
        a[0] * b[2] + a[1] * b[1] + a[2] * b[0],
        a[0] * b[3] + a[1] * b[2] + a[2] * b[1] + a[3] * b[0],
        // a[0] * b[4] + a[1] * b[3] + a[2] * b[2] + a[3] * b[1] + a[4] * b[0],
    ]
}

fn carbon_isotopes(count: u16) -> [f32; 4] {
    let lambda = count as f32 * 0.011;
    let mut c13 = [0.0; 4];

    let fact = [1, 1, 2, 6];
    for k in 0..4 {
        c13[k] = lambda.powi(k as i32) * f32::exp(-lambda) / fact[k] as f32;
    }
    c13
}

fn sulfur_isotopes(count: u16) -> [f32; 4] {
    let lambda33 = count as f32 * 0.0076;
    let lambda35 = count as f32 * 0.044;
    let mut s33 = [0.0; 4];
    let s35 = [
        lambda35.powi(0) * f32::exp(-lambda35),
        0.0,
        lambda35.powi(1) * f32::exp(-lambda35),
        0.0,
        // lambda35.powi(2) * f32::exp(-lambda35) / 2.0,
    ];

    let fact = [1, 1, 2, 6];
    for k in 0..4 {
        s33[k] = lambda33.powi(k as i32) * f32::exp(-lambda33) / fact[k] as f32;
    }

    convolve(&s33, &s35)
}

pub fn peptide_isotopes(carbons: u16, sulfurs: u16) -> [f32; 3] {
    let c = carbon_isotopes(carbons);
    let s = sulfur_isotopes(sulfurs);
    let mut c = convolve(&c, &s);
    let max = c[0].max(c[1]).max(c[2]);
    c.iter_mut().for_each(|val| *val /= max);
    [c[0], c[1], c[2]]
}

/// Approximate the first four isotope abundances for a peptide-like fragment.
///
/// Element counts are estimated from neutral mass with the averagine model.
/// Natural heavy-isotope counts are approximated as independent Poisson
/// processes at nominal mass shifts of one and two daltons. The truncated
/// pattern is normalized so it can be used directly for envelope scoring.
pub fn averagine_isotopes(neutral_mass: f32) -> [f32; 4] {
    const AVERAGINE_MASS: f32 = 111.1254;
    const CARBON: f32 = 4.9384;
    const HYDROGEN: f32 = 7.7583;
    const NITROGEN: f32 = 1.3577;
    const OXYGEN: f32 = 1.4773;
    const SULFUR: f32 = 0.0417;

    let scale = neutral_mass.max(0.0) / AVERAGINE_MASS;
    let lambda_one = scale
        * (CARBON * 0.0107
            + HYDROGEN * 0.000115
            + NITROGEN * 0.00368
            + OXYGEN * 0.00038
            + SULFUR * 0.0075);
    let lambda_two = scale * (OXYGEN * 0.00205 + SULFUR * 0.0425);

    let one_zero = (-lambda_one).exp();
    let one = [
        one_zero,
        one_zero * lambda_one,
        one_zero * lambda_one.powi(2) / 2.0,
        one_zero * lambda_one.powi(3) / 6.0,
    ];
    let two_zero = (-lambda_two).exp();
    let two_one = two_zero * lambda_two;
    let mut pattern = [
        one[0] * two_zero,
        one[1] * two_zero,
        one[2] * two_zero + one[0] * two_one,
        one[3] * two_zero + one[1] * two_one,
    ];
    let sum = pattern.iter().sum::<f32>();
    if sum > 0.0 {
        pattern.iter_mut().for_each(|value| *value /= sum);
    }
    pattern
}

/// Look up a four-peak averagine pattern in 25 Da neutral-mass bins.
pub fn averagine_isotopes_cached(neutral_mass: f32) -> [f32; 4] {
    const BIN_WIDTH: f32 = 25.0;
    const MAX_MASS: usize = 10_000;
    const BINS: usize = MAX_MASS / BIN_WIDTH as usize + 1;

    static CACHE: OnceLock<Vec<[f32; 4]>> = OnceLock::new();
    if neutral_mass > MAX_MASS as f32 {
        return averagine_isotopes(neutral_mass);
    }
    let cache = CACHE.get_or_init(|| {
        (0..BINS)
            .map(|bin| averagine_isotopes(bin as f32 * BIN_WIDTH))
            .collect()
    });
    let bin = (neutral_mass.max(0.0) / BIN_WIDTH).round() as usize;
    cache[bin.min(BINS - 1)]
}

#[cfg(test)]
#[path = "../tests/unit/isotopes.rs"]
mod tests;
