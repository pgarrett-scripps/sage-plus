use super::*;
use sage_core::mass::PROTON;

const ALPHA: f32 = 1500.7342;
const BETA: f32 = 1010.5127;
const TOL: Tolerance = Tolerance::Ppm(-10.0, 10.0);

fn dsso() -> CleavableLinker {
    CleavableLinker::dsso()
}

fn dsbu() -> CleavableLinker {
    CleavableLinker::dsbu()
}

fn settings() -> PairingSettings {
    PairingSettings {
        precursor_tol: TOL,
        isotope_errors: (0, 0),
        min_chain_mass: 500.0,
    }
}

fn precursor(mass: f32) -> [PrecursorMass; 1] {
    [PrecursorMass {
        mass,
        charge: 4,
        half_width: 0.0,
    }]
}

/// Build a spectrum from (neutral mass, charge, known charge, intensity).
pub(crate) fn spectrum(peaks: &[(f32, u8, bool, f32)]) -> ProcessedSpectrum {
    let mut stored: Vec<(f32, u8, bool, f32)> = peaks
        .iter()
        .map(|&(neutral, charge, known, intensity)| {
            if known {
                (neutral, charge, true, intensity)
            } else {
                // Unknown-charge peaks are stored as m/z - proton, charge 1.
                let mz = neutral / charge as f32 + PROTON;
                (mz - PROTON, 1, false, intensity)
            }
        })
        .collect();
    stored.sort_by(|a, b| a.0.total_cmp(&b.0));
    ProcessedSpectrum {
        level: 2,
        masses: stored.iter().map(|p| p.0).collect(),
        charges: stored.iter().map(|p| p.1).collect(),
        charge_is_known: stored.iter().map(|p| p.2).collect(),
        intensities: stored.iter().map(|p| p.3).collect(),
        total_ion_current: stored.iter().map(|p| p.3).sum(),
        ..Default::default()
    }
}

pub(crate) fn signature(
    linker: &CleavableLinker,
    chain: f32,
    charge: u8,
    known: bool,
) -> [(f32, u8, bool, f32); 2] {
    [
        (chain + linker.light_stub, charge, known, 100.0),
        (chain + linker.heavy_stub, charge, known, 80.0),
    ]
}

#[test]
fn four_peak_dsso_signature_yields_one_observed_pair() {
    let linker = dsso();
    let mut peaks = Vec::new();
    peaks.extend(signature(&linker, ALPHA, 2, true));
    peaks.extend(signature(&linker, BETA, 3, false));
    // Ordinary backbone fragments that must not form doublets.
    peaks.extend([
        (500.25, 1, true, 50.0),
        (620.31, 1, true, 40.0),
        (731.4, 2, true, 30.0),
    ]);
    let spectrum = spectrum(&peaks);
    let mass = ALPHA + BETA + linker.crosslink_mass;

    let doublets = find_doublets(&spectrum, &linker, 5, TOL);
    assert_eq!(doublets.len(), 2, "{doublets:?}");

    let pairs = pair_hypotheses(&doublets, &precursor(mass), &linker, settings());
    assert_eq!(pairs.len(), 1, "{pairs:?}");
    let pair = pairs[0];
    assert!(pair.both_observed());
    assert!(TOL.contains(ALPHA, pair.alpha_mass));
    assert!(TOL.contains(BETA, pair.beta_mass));
}

#[test]
fn single_doublet_infers_partner_at_each_isotope() {
    let linker = dsso();
    let spectrum = spectrum(&signature(&linker, ALPHA, 2, true));
    // Reported monoisotopic peak is one isotope too high.
    let mass = ALPHA + BETA + linker.crosslink_mass + NEUTRON;
    let doublets = find_doublets(&spectrum, &linker, 4, TOL);
    let pairs = pair_hypotheses(
        &doublets,
        &precursor(mass),
        &linker,
        PairingSettings {
            isotope_errors: (-1, 3),
            ..settings()
        },
    );
    assert_eq!(pairs.len(), 5);
    assert!(pairs.iter().all(|p| !p.both_observed() && p.alpha_observed));
    assert!(pairs.iter().any(|p| TOL.contains(BETA, p.beta_mass)));
}

#[test]
fn observed_chains_pair_despite_wrong_monoisotopic_mass() {
    let linker = dsso();
    let mut peaks = Vec::new();
    peaks.extend(signature(&linker, ALPHA, 2, true));
    peaks.extend(signature(&linker, BETA, 2, true));
    let spectrum = spectrum(&peaks);
    let doublets = find_doublets(&spectrum, &linker, 4, TOL);
    // Reported mass is 0.4 Da off: not an isotope error.
    let mass = ALPHA + BETA + linker.crosslink_mass + 0.4;
    let narrow = pair_hypotheses(&doublets, &precursor(mass), &linker, settings());
    assert!(narrow.iter().all(|p| !p.both_observed()));
    let window = [PrecursorMass {
        half_width: 0.7,
        ..precursor(mass)[0]
    }];
    let wide = pair_hypotheses(&doublets, &window, &linker, settings());
    assert!(wide[0].both_observed());
    assert!(TOL.contains(BETA, wide[0].beta_mass));

    // Without a window, an isotope-shifted precursor still pairs the
    // observed chains; the 0.4 Da offset does not, even inside the
    // isotope span.
    let shifted = precursor(ALPHA + BETA + linker.crosslink_mass + NEUTRON);
    let isotopes = PairingSettings {
        isotope_errors: (-1, 3),
        ..settings()
    };
    let exact = pair_hypotheses(&doublets, &shifted, &linker, isotopes);
    assert!(exact[0].both_observed());
    let off = pair_hypotheses(&doublets, &precursor(mass), &linker, isotopes);
    assert!(off.iter().all(|p| !p.both_observed()));
}

#[test]
fn dsbu_spacing_is_distinct_from_dsso() {
    let spectrum = spectrum(&signature(&dsbu(), ALPHA, 2, true));
    assert_eq!(find_doublets(&spectrum, &dsbu(), 4, TOL).len(), 1);
    assert!(find_doublets(&spectrum, &dsso(), 4, TOL).is_empty());
}

#[test]
fn doublet_members_must_share_charge() {
    let [light, heavy] = signature(&dsso(), ALPHA, 2, true);
    // Same neutral spacing, but the heavy member is deisotoped at charge 3.
    let spectrum = spectrum(&[light, (heavy.0, 3, true, heavy.3)]);
    assert!(find_doublets(&spectrum, &dsso(), 5, TOL).is_empty());
}
