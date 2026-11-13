//! Open search mass offsets.
//!
//! Moved out of `database.rs` unchanged; every item keeps its `database::` path.

use crate::modification::{ModificationDefinition, ModificationSpecificity, SiteMode};
use crate::peptide::Site;
use std::hash::Hash;
use std::sync::Arc;

/// Upper bound on distinct search-time offsets; candidate hypotheses encode
/// the offset in one byte, with zero reserved for "no offset".
pub const MAX_MASS_OFFSETS: usize = u8::MAX as usize - 1;

/// A modification tested at search time instead of being expanded into the
/// fragment index. At most one offset copy is placed on a peptide.
#[derive(Clone, Debug)]
pub struct MassOffset {
    pub definition: Arc<ModificationDefinition>,
    /// Every configured site rule for this definition.
    pub specificities: Vec<ModificationSpecificity>,
    pub site_mode: SiteMode,
}

impl MassOffset {
    pub fn mass(&self) -> f32 {
        self.definition.mass
    }

    /// Mass difference of the first generated fragment form containing the
    /// modification. This mirrors the preliminary fragment index, which keeps
    /// the first variant of every ion group; ion series sort losses ascending,
    /// so a required loss contributes its smallest configured mass. With
    /// per-site losses the preliminary search uses one shift for every site:
    /// the smallest required loss only when every site has a loss, since a
    /// site without losses keeps the retained form.
    pub fn fragment_shift(&self) -> f32 {
        match self.definition.neutral_loss_mode {
            crate::modification::NeutralLossMode::Required => {
                let smallest = match &self.definition.site_losses {
                    None => self
                        .definition
                        .neutral_losses
                        .iter()
                        .copied()
                        .min_by(f32::total_cmp)
                        .unwrap_or_default(),
                    Some(rules) => {
                        let losses = |specificity: &ModificationSpecificity| {
                            rules
                                .iter()
                                .find(|(rule, _)| rule == specificity)
                                .map(|(_, losses)| losses)
                        };
                        if self.specificities.iter().all(|s| losses(s).is_some()) {
                            self.specificities
                                .iter()
                                .filter_map(losses)
                                .flat_map(|losses| losses.iter().copied())
                                .min_by(f32::total_cmp)
                                .unwrap_or_default()
                        } else {
                            0.0
                        }
                    }
                };
                self.definition.mass - smallest
            }
            crate::modification::NeutralLossMode::Optional => self.definition.mass,
        }
    }
}

/// The offset placement selected for a candidate. `offset` indexes
/// [`crate::database::IndexedDatabase::mass_offsets`].
#[derive(Copy, Clone, Debug, PartialEq, Eq, Hash)]
pub struct MassOffsetAssignment {
    pub offset: u16,
    pub site: Site,
}
