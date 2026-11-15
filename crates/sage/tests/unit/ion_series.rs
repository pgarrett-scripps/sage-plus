#![allow(clippy::excessive_precision)]

use super::*;
use crate::modification::{ModificationDefinition, ModificationSpecificity, NeutralLossMode};
use crate::peptide::Peptide;
use crate::{enzyme::Digest, mass::PROTON};
use std::{collections::HashMap, sync::Arc};

fn peptide(s: &str) -> Peptide {
    Peptide::try_from(Digest {
        sequence: s.into(),
        ..Default::default()
    })
    .unwrap()
}

fn peptide_with_loss(mode: NeutralLossMode) -> Peptide {
    let modification = Arc::new(ModificationDefinition {
        mass: 20.0,
        name: Some(Arc::from("TestMod")),
        neutral_losses: Arc::from([10.0]),
        site_losses: None,
        neutral_loss_mode: mode,
        channel_offsets: Arc::default(),
    });
    peptide("AMK")
        .apply(
            &[(
                ModificationSpecificity::Residue(b'M'),
                modification,
                Some(1),
            )],
            &HashMap::default(),
            1,
            None,
        )
        .into_iter()
        .find(|peptide| peptide.to_string().contains("TestMod"))
        .unwrap()
}

fn check_within<I: Iterator<Item = Ion>>(iter: I, expected_mz: &[f32]) {
    let observed = iter.map(|ion| ion.monoisotopic_mass).collect::<Vec<f32>>();
    assert_eq!(expected_mz.len(), observed.len());
    assert!(
        expected_mz
            .iter()
            .zip(observed.iter())
            .all(|(a, b)| (a - b).abs() < 0.005),
        "{:?}",
        expected_mz
            .iter()
            .zip(observed.iter())
            .map(|(a, b)| a - b)
            .collect::<Vec<_>>()
    );
}

macro_rules! ions {
    ($peptide:expr, $kind:expr, $charge:expr) => {{
        IonSeries::new($peptide, $kind).map(|mut ion| {
            ion.monoisotopic_mass = (ion.monoisotopic_mass + $charge * PROTON) / $charge;
            ion
        })
    }};
}

#[test]
fn abc_xyz() {
    let peptide = peptide("PEPTIDE");
    let expected_a = vec![70.065, 199.108, 296.160, 397.208, 510.292, 625.32];
    let expected_b = vec![98.0600, 227.1026, 324.155, 425.2030, 538.287, 653.314];
    let expected_c = vec![115.086, 244.129, 341.182, 442.229, 555.314, 670.341];
    let expected_x = vec![729.294, 600.251, 503.198, 402.151, 289.066, 174.039];
    let expected_y = vec![703.314, 574.2719, 477.219, 376.171, 263.0874, 148.0604];
    let expected_z = vec![686.288, 557.245, 460.193, 359.145, 246.061, 131.034];

    check_within(ions!(&peptide, Kind::A, 1.0), &expected_a);
    check_within(ions!(&peptide, Kind::B, 1.0), &expected_b);
    check_within(ions!(&peptide, Kind::C, 1.0), &expected_c);
    check_within(ions!(&peptide, Kind::X, 1.0), &expected_x);
    check_within(ions!(&peptide, Kind::Y, 1.0), &expected_y);
    check_within(ions!(&peptide, Kind::Z, 1.0), &expected_z);
}

#[test]
fn iterate_b_ions() {
    let peptide = peptide("PEPTIDE");

    // Charge state 2
    let expected_mz = vec![
        98.06004, 227.10263, 324.155_4, 425.203_06, 538.287_2, 653.314_1,
    ];

    check_within(ions!(&peptide, Kind::B, 1.0), &expected_mz);
}

#[test]
fn iterate_y_ions() {
    let peptide = peptide("PEPTIDE");

    // Charge state 1
    let expected_mz = vec![
        703.31447, 574.27188, 477.21912, 376.17144, 263.08737, 148.06043,
    ];

    check_within(ions!(&peptide, Kind::Y, 1.0), &expected_mz);
}

#[test]
fn y_index() {
    let peptide = peptide("PEPTIDE");
    // Charge state 1
    let expected_ion: Vec<(usize, f32)> = vec![
        (6, 703.31447),
        (5, 574.27188),
        (4, 477.21912),
        (3, 376.17144),
        (2, 263.08737),
        (1, 148.06043),
    ];
    assert!(IonSeries::new(&peptide, Kind::Y)
        .enumerate()
        .map(|(idx, ion)| (peptide.sequence.len().saturating_sub(1) - idx, ion))
        .zip(expected_ion.into_iter())
        .all(|((idx, ion), (idx_, mz))| {
            idx == idx_ && (ion.monoisotopic_mass + PROTON - mz).abs() <= 0.01
        }),)
}

#[test]
fn index_filtering() {
    let peptide = &peptide("PEPTIDE");
    let ions = IonSeries::new(peptide, Kind::B)
        .enumerate()
        .chain(IonSeries::new(peptide, Kind::Y).enumerate())
        .filter(|(ion_idx, ion)| {
            // Don't store b1, b2, y1, y2 ions for preliminary scoring
            match ion.kind {
                Kind::A | Kind::B | Kind::C => (ion_idx + 1) > 2,
                Kind::X | Kind::Y | Kind::Z | Kind::ZDot => {
                    peptide.sequence.len().saturating_sub(1) - ion_idx > 2
                }
            }
        })
        .map(|(_, mut ion)| {
            ion.monoisotopic_mass += PROTON;
            ion
        })
        .collect::<Vec<_>>();

    #[rustfmt::skip]
        let expected = vec![
            Ion { kind: Kind::B, monoisotopic_mass: 324.155397 },
            Ion { kind: Kind::B, monoisotopic_mass: 425.203076 },
            Ion { kind: Kind::B, monoisotopic_mass: 538.287140 },
            Ion { kind: Kind::B, monoisotopic_mass: 653.314083 },
            Ion { kind: Kind::Y, monoisotopic_mass: 703.314477 },
            Ion { kind: Kind::Y, monoisotopic_mass: 574.271884 },
            Ion { kind: Kind::Y, monoisotopic_mass: 477.219120 },
            Ion { kind: Kind::Y, monoisotopic_mass: 376.171441 },
        ];

    assert_eq!(expected.len(), ions.len(), "{:?}\n{:?}", ions, expected);
    assert!(
        ions.iter().zip(expected.iter()).all(|(left, right)| {
            left.kind == right.kind && (left.monoisotopic_mass - right.monoisotopic_mass) <= 0.1
        }),
        "{:?}",
        ions
    );
}

#[test]
fn decoy() {
    let peptide_ = peptide("PEPTIDE");

    // Charge state 2
    let expected_mz = vec![
        352.16087, 287.639_6, 239.11319, 188.58935, 132.04732, 74.53385,
    ];

    check_within(ions!(&peptide_, Kind::Y, 2.0), &expected_mz);

    let peptide = peptide("EDITPEP");

    // Charge state 2
    let expected_mz = vec![
        336.16596, 278.652_5, 222.110_46, 171.586_62, 123.060237, 58.538_94,
    ];

    check_within(ions!(&peptide, Kind::Y, 2.0), &expected_mz);
}

#[test]
fn nterm_mod() {
    let static_mods = [(ModificationSpecificity::PeptideN(None), 229.01)].into();
    let peptide = peptide("PEPTIDE")
        .apply(&[], &static_mods, 1, None)
        .remove(0);

    // Charge state 1, b-ions should be TMT tagged
    let expected_b = [
        98.06004, 227.10263, 324.155_4, 425.203_06, 538.287_2, 653.314_1,
    ]
    .into_iter()
    .map(|x| x + 229.01)
    .collect::<Vec<_>>();

    // y-ions shouldn't have TMT tag
    let expected_y = vec![
        703.31447, 574.27188, 477.21912, 376.17144, 263.08737, 148.06043,
    ];

    check_within(ions!(&peptide, Kind::B, 1.0), &expected_b);
    check_within(ions!(&peptide, Kind::Y, 1.0), &expected_y);
}

#[test]
fn cterm_mod() {
    let static_mods = [(ModificationSpecificity::PeptideC(None), 229.01)].into();
    let peptide = peptide("PEPTIDE")
        .apply(&[], &static_mods, 1, None)
        .remove(0);
    assert!((peptide.monoisotopic - 1028.37).abs() < 0.001);

    // b-ions should not be tagged
    let expected_b = [
        98.06004, 227.10263, 324.155_4, 425.203_06, 538.287_2, 653.314_1,
    ];

    // y-ions should be tagged
    let expected_y = vec![
        703.31447, 574.27188, 477.21912, 376.17144, 263.08737, 148.06043,
    ]
    .into_iter()
    .map(|x| x + 229.01)
    .collect::<Vec<_>>();

    check_within(ions!(&peptide, Kind::B, 1.0), &expected_b);
    check_within(ions!(&peptide, Kind::Y, 1.0), &expected_y);
}

#[test]
fn internal_mod() {
    let peptide = peptide("PEPTIDE");
    let static_mods = [(ModificationSpecificity::Residue(b'I'), 29.0)].into();
    let peptide = peptide.apply(&[], &static_mods, 1, None).remove(0);

    let expected_b = [
        98.06004,
        227.10263,
        324.155_4,
        425.203_06,
        538.287_2 + 29.0,
        653.314_1 + 29.0,
    ];

    let expected_y = vec![
        703.31447 + 29.0,
        574.27188 + 29.0,
        477.21912 + 29.0,
        376.17144 + 29.0,
        263.08737,
        148.06043,
    ];

    check_within(ions!(&peptide, Kind::B, 1.0), &expected_b);
    check_within(ions!(&peptide, Kind::Y, 1.0), &expected_y);
}

#[test]
fn optional_neutral_loss_keeps_retained_fragment() {
    let peptide = peptide_with_loss(NeutralLossMode::Optional);
    let groups = IonGroupSeries::new(&peptide, Kind::B).collect::<Vec<_>>();

    assert_eq!(groups[0].variants.len(), 1); // b1 does not contain M
    assert_eq!(groups[1].variants.len(), 2); // b2 contains M
    let retained = IonSeries::new(&peptide, Kind::B).nth(1).unwrap();
    assert!(groups[1].variants.iter().any(|variant| {
        variant.neutral_loss.is_none()
            && (variant.monoisotopic_mass - retained.monoisotopic_mass).abs() < 1e-5
    }));
    assert!(groups[1].variants.iter().any(|variant| {
        variant.neutral_loss == Some(10.0)
            && (variant.monoisotopic_mass - (retained.monoisotopic_mass - 10.0)).abs() < 1e-5
    }));
}

#[test]
fn common_ion_groups_use_inline_variant_storage() {
    let peptide = peptide("PEPTIDE");
    assert!(IonGroupSeries::new(&peptide, Kind::B).all(|group| !group.variants.spilled()));

    let peptide = peptide_with_loss(NeutralLossMode::Optional);
    assert!(IonGroupSeries::new(&peptide, Kind::B).all(|group| !group.variants.spilled()));
}

#[test]
fn required_neutral_loss_suppresses_retained_fragment() {
    let peptide = peptide_with_loss(NeutralLossMode::Required);
    let b = IonGroupSeries::new(&peptide, Kind::B).collect::<Vec<_>>();
    let y = IonGroupSeries::new(&peptide, Kind::Y).collect::<Vec<_>>();

    assert_eq!(b[0].variants.len(), 1); // b1 does not contain M
    assert_eq!(b[1].variants.len(), 1);
    assert_eq!(b[1].variants[0].neutral_loss, Some(10.0));
    assert_eq!(y[0].variants.len(), 1);
    assert_eq!(y[0].variants[0].neutral_loss, Some(10.0));
    assert_eq!(y[1].variants.len(), 1); // y1 does not contain M
    assert_eq!(y[1].variants[0].neutral_loss, None);
}

#[test]
fn neutral_losses_combine_across_multiple_modified_sites() {
    let modification = Arc::new(ModificationDefinition {
        mass: 20.0,
        name: Some(Arc::from("TestMod")),
        neutral_losses: Arc::from([10.0]),
        site_losses: None,
        neutral_loss_mode: NeutralLossMode::Optional,
        channel_offsets: Arc::default(),
    });
    let peptide = peptide("MMK")
        .apply(
            &[(
                ModificationSpecificity::Residue(b'M'),
                modification,
                Some(2),
            )],
            &HashMap::default(),
            2,
            None,
        )
        .into_iter()
        .find(|peptide| peptide.to_string().matches("TestMod").count() == 2)
        .unwrap();

    let b2 = IonGroupSeries::new(&peptide, Kind::B).nth(1).unwrap();
    assert!(b2.variants.spilled());
    let losses = b2
        .variants
        .iter()
        .map(|variant| variant.neutral_loss.unwrap_or(0.0))
        .collect::<Vec<_>>();
    assert_eq!(losses, vec![0.0, 10.0, 20.0]);
}

#[test]
fn terminal_required_losses_affect_only_containing_series() {
    let modification = Arc::new(ModificationDefinition {
        mass: 20.0,
        name: Some(Arc::from("TerminalMod")),
        neutral_losses: Arc::from([10.0]),
        site_losses: None,
        neutral_loss_mode: NeutralLossMode::Required,
        channel_offsets: Arc::default(),
    });
    let peptide = peptide("AMK")
        .apply(
            &[(
                ModificationSpecificity::PeptideN(None),
                modification,
                Some(1),
            )],
            &HashMap::default(),
            1,
            None,
        )
        .into_iter()
        .find(|peptide| peptide.to_string().contains("TerminalMod"))
        .unwrap();

    assert!(IonGroupSeries::new(&peptide, Kind::B).all(|group| group
        .variants
        .iter()
        .all(|variant| variant.neutral_loss == Some(10.0))));
    assert!(IonGroupSeries::new(&peptide, Kind::Y).all(|group| group
        .variants
        .iter()
        .all(|variant| variant.neutral_loss.is_none())));
}

#[test]
fn z_dot_ions_are_one_hydrogen_heavier_than_z() {
    // z• (radical z, from ETD/EThcD) = z + H; PEPTIDE z•1 (E) is 132.042.
    let peptide = peptide("PEPTIDE");
    let expected_z_dot = vec![687.296, 558.253, 461.201, 360.153, 247.069, 132.042];
    check_within(ions!(&peptide, Kind::ZDot, 1.0), &expected_z_dot);

    for label in ["\"z_dot\"", "\"zdot\"", "\"z.\""] {
        assert_eq!(serde_json::from_str::<Kind>(label).unwrap(), Kind::ZDot);
    }
    assert_eq!(serde_json::to_string(&Kind::ZDot).unwrap(), "\"z_dot\"");
}

/// Assert m/z to 0.2 mDa: tight enough to catch the 0.549 mDa (one electron)
/// ammonia error this file pins, loose enough for f32 rounding over a
/// six-residue running sum.
fn check_exactly<I: Iterator<Item = Ion>>(iter: I, expected_mz: &[f32]) {
    let observed = iter.map(|ion| ion.monoisotopic_mass).collect::<Vec<f32>>();
    assert_eq!(expected_mz.len(), observed.len());
    for (expected, got) in expected_mz.iter().zip(observed.iter()) {
        assert!(
            (expected - got).abs() < 0.0002,
            "expected {expected}, got {got}, off by {:.6}",
            got - expected
        );
    }
}

/// Every series of PEPTIDE against m/z computed by hand from the element
/// masses, pinned tightly enough to see the ammonia error.
///
/// Neutral ammonia is N + 3H = 17.026549 and NH2 is N + 2H = 16.018724: a c ion
/// is b + NH3, a z ion is y - NH3, a radical z-dot ion is y - NH2 and an x ion
/// is y + CO - 2H. The third hydrogen of NH3 used to be a proton (1.0072764),
/// one electron light, so c ions came out 0.549 mDa too light (4.8 ppm at c1,
/// 0.8 ppm at c6) and x, z and z-dot ions 0.549 mDa too heavy. a and b ions
/// never involved it and are pinned here to show they did not move.
#[test]
fn series_masses_pin_the_ammonia_constant() {
    let peptide = peptide("PEPTIDE");

    check_exactly(
        ions!(&peptide, Kind::A, 1.0),
        &[
            70.06512, 199.1077, 296.16046, 397.20816, 510.29221, 625.31909,
        ],
    );
    check_exactly(
        ions!(&peptide, Kind::B, 1.0),
        &[
            98.06004, 227.10262, 324.15536, 425.20306, 538.28711, 653.31403,
        ],
    );
    check_exactly(
        ions!(&peptide, Kind::C, 1.0),
        &[
            115.08659, 244.12917, 341.18192, 442.22961, 555.31366, 670.34058,
        ],
    );
    check_exactly(
        ions!(&peptide, Kind::X, 1.0),
        &[
            729.2937, 600.2511, 503.19836, 402.15067, 289.06662, 174.03967,
        ],
    );
    check_exactly(
        ions!(&peptide, Kind::Y, 1.0),
        &[
            703.31445, 574.27185, 477.21909, 376.17139, 263.08734, 148.06041,
        ],
    );
    check_exactly(
        ions!(&peptide, Kind::Z, 1.0),
        &[
            686.2879, 557.2453, 460.19254, 359.14484, 246.06081, 131.03386,
        ],
    );
    check_exactly(
        ions!(&peptide, Kind::ZDot, 1.0),
        &[
            687.29572, 558.25311, 461.20035, 360.15265, 247.06863, 132.04169,
        ],
    );
}

/// A z-dot ion is exactly one hydrogen heavier than the even-electron z ion,
/// and exactly NH2 lighter than y.
#[test]
fn z_dot_minus_z_is_exactly_one_hydrogen() {
    let peptide = peptide("PEPTIDE");
    let z = ions!(&peptide, Kind::Z, 1.0).map(|ion| ion.monoisotopic_mass);
    let z_dot = ions!(&peptide, Kind::ZDot, 1.0)
        .map(|ion| ion.monoisotopic_mass)
        .collect::<Vec<_>>();
    for (z, z_dot) in z.zip(z_dot) {
        assert!(
            (z_dot - z - 1.007825).abs() < 0.0002,
            "z {z}, z-dot {z_dot}, difference {:.6}",
            z_dot - z
        );
    }
}

/// `z_radical` is upstream's name for the same ion (`cb9934f`); we accept it in
/// a configuration and keep writing `z_dot`.
#[test]
fn z_radical_is_accepted_as_z_dot() {
    let kinds: Vec<Kind> = serde_json::from_str(r#"["c","z_radical","z"]"#).unwrap();
    assert_eq!(kinds, vec![Kind::C, Kind::ZDot, Kind::Z]);
    assert_eq!(
        serde_json::to_string(&kinds).unwrap(),
        r#"["c","z_dot","z"]"#
    );
    for alias in [r#""z_dot""#, r#""zdot""#, r#""z.""#, r#""z_radical""#] {
        assert_eq!(serde_json::from_str::<Kind>(alias).unwrap(), Kind::ZDot);
    }
}
