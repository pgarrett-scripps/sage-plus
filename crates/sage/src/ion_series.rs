use serde::{Deserialize, Serialize};
use smallvec::{smallvec, SmallVec};

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
    pub variants: SmallVec<[IonVariant; 2]>,
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

    fn losses(&self, series_index: usize) -> SmallVec<[f32; 4]> {
        let mut totals = smallvec![0.0f32];
        for applied in self.peptide.applied_modifications().filter(|applied| {
            self.contains_site(applied.site, series_index)
                && !applied.modification.neutral_losses.is_empty()
        }) {
            let option_count = applied.modification.neutral_losses.len()
                + usize::from(applied.modification.neutral_loss_mode == NeutralLossMode::Optional);
            let mut next: SmallVec<[f32; 4]> =
                SmallVec::with_capacity(totals.len().saturating_mul(option_count));
            for total in &totals {
                if applied.modification.neutral_loss_mode == NeutralLossMode::Optional {
                    next.push(*total);
                }
                for loss in applied.modification.neutral_losses.iter() {
                    next.push(total + loss);
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
            .collect::<SmallVec<[_; 2]>>();

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
    modification_idx: usize,
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
            modification_idx: 0,
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
        let m = self
            .peptide
            .modifications
            .mass_at_with_cursor(self.idx, &mut self.modification_idx);

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
#[path = "../tests/unit/ion_series.rs"]
mod test;
