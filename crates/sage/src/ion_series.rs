use serde::{Deserialize, Serialize};

use crate::mass::monoisotopic;
use crate::modification::NeutralLossMode;
use crate::peptide::{Peptide, Site};

#[derive(Copy, Clone, Debug, PartialEq, Eq, Hash, PartialOrd, Deserialize, Serialize)]
#[serde(rename_all = "lowercase")]
pub enum Kind {
    A,
    B,
    C,
    X,
    Y,
    Z,
}

/// Theoretical B/Y ion
#[derive(Copy, Clone, Debug)]
pub struct Ion {
    /// B or Y ion
    pub kind: Kind,
    /// Neutral fragment mass (no charge)
    pub monoisotopic_mass: f32,
}

#[derive(Clone, Debug)]
pub struct IonVariant {
    pub kind: Kind,
    pub monoisotopic_mass: f32,
    /// Total neutral loss represented by this fragment variant. `None` is the
    /// retained (no-loss) form.
    pub neutral_loss: Option<f32>,
}

#[derive(Clone, Debug)]
pub struct IonGroup {
    pub kind: Kind,
    /// Zero-based series index used by scoring and minimum-ion filtering.
    pub series_index: usize,
    pub variants: Vec<IonVariant>,
}

/// Generate groups of mutually alternative fragment forms for each peptide
/// cleavage. A group contains the retained ion when every applicable
/// modification has optional loss behavior, plus the configured neutral-loss
/// combinations. Required losses remove the retained option for fragments
/// containing that modification.
pub struct IonGroupSeries<'p> {
    peptide: &'p Peptide,
    base: IonSeries<'p>,
    series_index: usize,
}

impl<'p> IonGroupSeries<'p> {
    pub fn new(peptide: &'p Peptide, kind: Kind) -> Self {
        Self {
            peptide,
            base: IonSeries::new(peptide, kind),
            series_index: 0,
        }
    }

    fn contains_site(&self, site: Site, series_index: usize) -> bool {
        match (self.base.kind, site) {
            (Kind::A | Kind::B | Kind::C, Site::Nterm) => true,
            (Kind::A | Kind::B | Kind::C, Site::Cterm) => false,
            (Kind::A | Kind::B | Kind::C, Site::Sequence(index)) => index as usize <= series_index,
            (Kind::X | Kind::Y | Kind::Z, Site::Nterm) => false,
            (Kind::X | Kind::Y | Kind::Z, Site::Cterm) => true,
            (Kind::X | Kind::Y | Kind::Z, Site::Sequence(index)) => index as usize > series_index,
        }
    }

    fn losses(&self, series_index: usize) -> Vec<f32> {
        let mut totals = vec![0.0f32];
        for applied in self.peptide.applied_modifications.iter().filter(|applied| {
            self.contains_site(applied.site, series_index)
                && !applied.modification.neutral_losses.is_empty()
        }) {
            let mut options = Vec::with_capacity(applied.modification.neutral_losses.len() + 1);
            if applied.modification.neutral_loss_mode == NeutralLossMode::Optional {
                options.push(0.0);
            }
            options.extend(applied.modification.neutral_losses.iter().copied());

            let mut next = Vec::with_capacity(totals.len().saturating_mul(options.len()));
            for total in &totals {
                for option in &options {
                    next.push(total + option);
                }
            }
            next.sort_unstable_by(f32::total_cmp);
            next.dedup_by(|a, b| (*a - *b).abs() < 1e-5);
            totals = next;
        }
        totals
    }
}

impl Iterator for IonGroupSeries<'_> {
    type Item = IonGroup;

    fn next(&mut self) -> Option<Self::Item> {
        let ion = self.base.next()?;
        let series_index = self.series_index;
        self.series_index += 1;

        let variants = self
            .losses(series_index)
            .into_iter()
            .filter_map(|loss| {
                let mass = ion.monoisotopic_mass - loss;
                (mass > 0.0).then_some(IonVariant {
                    kind: ion.kind,
                    monoisotopic_mass: mass,
                    neutral_loss: (loss > 0.0).then_some(loss),
                })
            })
            .collect();

        Some(IonGroup {
            kind: ion.kind,
            series_index,
            variants,
        })
    }
}

/// Generate B/Y ions for a candidate peptide under a given charge state
pub struct IonSeries<'p> {
    pub kind: Kind,
    cumulative_mass: f32,
    peptide: &'p Peptide,
    idx: usize,
}

impl<'p> IonSeries<'p> {
    /// Create a new [`IonSeries`] iterator for a specified peptide
    pub fn new(peptide: &'p Peptide, kind: Kind) -> Self {
        const C: f32 = 12.0;
        const O: f32 = 15.994914;
        const H: f32 = 1.007825;
        const PRO: f32 = 1.0072764;
        const N: f32 = 14.003074;
        const NH3: f32 = N + H * 2.0 + PRO;

        let cumulative_mass = match kind {
            Kind::A => peptide.nterm.unwrap_or_default() - (C + O),
            Kind::B => peptide.nterm.unwrap_or_default(),
            Kind::C => peptide.nterm.unwrap_or_default() + NH3,
            Kind::X => {
                peptide.monoisotopic - peptide.nterm.unwrap_or_default() + (C + O - NH3 + N + H)
            }
            Kind::Y => peptide.monoisotopic - peptide.nterm.unwrap_or_default(),
            Kind::Z => peptide.monoisotopic - peptide.nterm.unwrap_or_default() - NH3,
        };
        Self {
            kind,
            cumulative_mass,
            peptide,
            idx: 0,
        }
    }
}

impl<'p> Iterator for IonSeries<'p> {
    type Item = Ion;

    // Dynamic programming solution - memoize cumulative mass of
    // peptide fragment for fast fragment ion generation
    fn next(&mut self) -> Option<Self::Item> {
        if self.idx >= self.peptide.sequence.len() - 1 {
            return None;
        }
        let r = self.peptide.sequence[self.idx];
        let m = self.peptide.modifications[self.idx];

        self.cumulative_mass += match self.kind {
            Kind::A | Kind::B | Kind::C => monoisotopic(r) + m,
            Kind::X | Kind::Y | Kind::Z => -(monoisotopic(r) + m),
        };
        self.idx += 1;

        Some(Ion {
            kind: self.kind,
            monoisotopic_mass: self.cumulative_mass,
        })
    }
}

#[cfg(test)]
mod test {
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
            neutral_loss_mode: mode,
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
                let ion_idx_filter = match ion.kind {
                    Kind::A | Kind::B | Kind::C => (ion_idx + 1) > 2,
                    Kind::X | Kind::Y | Kind::Z => {
                        peptide.sequence.len().saturating_sub(1) - ion_idx > 2
                    }
                };
                ion_idx_filter
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
            neutral_loss_mode: NeutralLossMode::Optional,
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
            neutral_loss_mode: NeutralLossMode::Required,
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
}
