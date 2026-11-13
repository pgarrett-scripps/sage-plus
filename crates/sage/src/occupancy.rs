//! Occupancy of the searched variable modifications: how often the sites a
//! modification could occupy actually carry it.
//!
//! One row per searched definition and site rule, counted over distinct
//! confident peptidoforms rather than PSMs, so an abundant peptide does not
//! weigh more than a scarce one. Nothing here knows the name of any
//! modification: every rule comes from the configuration, and sites come from
//! the same resolver the search and the localizer use
//! ([`crate::database::Parameters::localization_rules`]), so a site counts as
//! eligible exactly when the search could have modified it.

use crate::database::Parameters;
use crate::modification::{ModificationDefinition, ModificationSpecificity, SearchMode};
use crate::peptide::Peptide;
use std::sync::Arc;

/// One searched variable modification definition at one of its site rules.
#[derive(Debug, Clone)]
pub struct OccupancyRule {
    pub site_rule: ModificationSpecificity,
    pub definition: Arc<ModificationDefinition>,
    pub search_mode: SearchMode,
}

impl OccupancyRule {
    /// The site rule as it is written in a configuration key.
    pub fn site_rule_name(&self) -> String {
        self.site_rule.explicit_name()
    }

    pub fn mass(&self) -> f32 {
        self.definition.mass
    }

    /// The configured name, absent for a definition given as a bare mass.
    pub fn name(&self) -> Option<&str> {
        self.definition.name.as_deref()
    }

    /// `database` for a definition expanded into the fragment index,
    /// `mass_offset` for one tested at search time.
    pub fn search_mode_name(&self) -> &'static str {
        match self.search_mode {
            SearchMode::Database => "database",
            SearchMode::MassOffset => "mass_offset",
        }
    }

    fn key(&self) -> (ModificationSpecificity, &ModificationDefinition) {
        (self.site_rule, self.definition.as_ref())
    }
}

/// Sites one rule could occupy and the sites that carry it.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct OccupancyCounts {
    /// Peptidoforms with at least one eligible site.
    pub peptidoforms: usize,
    /// Peptidoforms carrying the definition on at least one eligible site.
    pub modified_peptidoforms: usize,
    /// Sites the rule selects, summed over the peptidoforms.
    pub eligible_sites: usize,
    /// Eligible sites carrying the definition.
    pub modified_sites: usize,
}

impl OccupancyCounts {
    /// Modified sites as a percent of eligible sites; zero without an
    /// eligible site.
    pub fn occupancy_pct(&self) -> f64 {
        if self.eligible_sites == 0 {
            0.0
        } else {
            100.0 * self.modified_sites as f64 / self.eligible_sites as f64
        }
    }
}

/// Every searched variable modification definition and site rule, in the
/// order [`Parameters::localization_rules`] resolves them, with a definition
/// configured twice at one rule counted once.
pub fn rules(parameters: &Parameters) -> Vec<OccupancyRule> {
    let mut rules = parameters
        .variable_mods
        .iter()
        .flat_map(|(site_rule, entries)| {
            entries.iter().map(move |entry| OccupancyRule {
                site_rule: *site_rule,
                definition: Arc::new(entry.definition()),
                search_mode: entry.search_mode(),
            })
        })
        .collect::<Vec<_>>();
    rules.sort_by(|left, right| {
        left.site_rule
            .cmp(&right.site_rule)
            .then_with(|| left.definition.cmp(&right.definition))
    });
    rules.dedup_by(|left, right| left.key() == right.key());
    rules
}

/// Count each rule's eligible and modified sites over `peptidoforms`, which
/// the caller has already filtered and deduplicated. The returned counts are
/// parallel to `rules`.
///
/// A site counts as modified when the modification applied there is this
/// definition: the same mass, name, losses and channel offsets. A peptidoform
/// without an eligible site for a rule does not count towards that rule.
pub fn count<'a, I>(
    parameters: &Parameters,
    rules: &[OccupancyRule],
    peptidoforms: I,
) -> Vec<OccupancyCounts>
where
    I: IntoIterator<Item = &'a Peptide>,
{
    let mut counts = vec![OccupancyCounts::default(); rules.len()];
    for peptide in peptidoforms {
        for resolved in parameters.localization_rules(peptide) {
            if resolved.sites.is_empty() {
                continue;
            }
            let key = (resolved.specificity, resolved.definition.as_ref());
            let Ok(index) = rules.binary_search_by(|rule| {
                rule.site_rule
                    .cmp(&key.0)
                    .then_with(|| rule.definition.as_ref().cmp(key.1))
            }) else {
                // `rules` and `localization_rules` read the same
                // configuration, so every resolved rule is one of ours.
                continue;
            };
            let modified = peptide
                .applied_modifications()
                .filter(|applied| {
                    applied.modification == resolved.definition.as_ref()
                        && resolved.sites.contains(&applied.site)
                })
                .count();
            let counts = &mut counts[index];
            counts.peptidoforms += 1;
            counts.eligible_sites += resolved.sites.len();
            counts.modified_sites += modified;
            if modified > 0 {
                counts.modified_peptidoforms += 1;
            }
        }
    }
    counts
}

#[cfg(test)]
#[path = "../tests/unit/occupancy.rs"]
mod test;
