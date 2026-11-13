//! Enumeration of modified peptidoforms from a peptide and its rules.
//!
//! Moved out of `peptide.rs` unchanged; every item keeps its `peptide::` path.

use std::{collections::HashMap, fmt::Debug, sync::Arc};

use super::*;
use crate::modification::{ModificationDefinition, ModificationSpecificity, SiteMode};
use fnv::FnvHashSet;

pub trait ModificationSource {
    fn definition(&self) -> Arc<ModificationDefinition>;
}

#[derive(Clone)]
pub(crate) struct VariableRule {
    pub specificity: ModificationSpecificity,
    pub modification: Arc<ModificationDefinition>,
    /// Per-peptide limit on new (non-library) placements of this modification.
    pub max_count: Option<usize>,
    /// Per-peptide limit on all placements, library-supported included.
    /// Defaults to `max_count`.
    pub max_total_count: Option<usize>,
    pub site_mode: SiteMode,
    /// Rules with the same named modification share one occurrence counter.
    pub count_group: usize,
}

impl VariableRule {
    pub(crate) fn total_limit(&self) -> Option<usize> {
        self.max_total_count.or(self.max_count)
    }
}

#[derive(Clone, Debug)]
pub(crate) struct LibrarySite {
    pub attachment: crate::ptm_library::Attachment,
    pub position: u32,
    pub modification: Arc<str>,
}

#[derive(Copy, Clone)]
pub(super) struct ModificationCandidate {
    pub(super) site: Site,
    pub(super) rule: usize,
    pub(super) library_supported: bool,
}

pub(crate) struct ModificationPlan<'a> {
    pub(super) variable_mods: &'a [VariableRule],
    pub(super) static_mods: &'a HashMap<ModificationSpecificity, Arc<ModificationDefinition>>,
    pub(super) lookup: Arc<ModificationLookup>,
    max_exhaustive_mods: usize,
    pub(super) max_total_mods: usize,
    max_combinations: Option<usize>,
}

impl<'a> ModificationPlan<'a> {
    pub(crate) fn new(
        variable_mods: &'a [VariableRule],
        static_mods: &'a HashMap<ModificationSpecificity, Arc<ModificationDefinition>>,
        lookup: Arc<ModificationLookup>,
        max_exhaustive_mods: usize,
        max_total_mods: usize,
        max_combinations: Option<usize>,
    ) -> Self {
        Self {
            variable_mods,
            static_mods,
            lookup,
            max_exhaustive_mods,
            max_total_mods,
            max_combinations,
        }
    }
}

impl ModificationSource for f32 {
    fn definition(&self) -> Arc<ModificationDefinition> {
        Arc::new(ModificationDefinition::bare(*self))
    }
}

impl ModificationSource for Arc<ModificationDefinition> {
    fn definition(&self) -> Arc<ModificationDefinition> {
        self.clone()
    }
}

pub(super) fn nonzero_mass(mass: f32) -> Option<f32> {
    (mass.abs() >= 1e-5).then_some(mass)
}

pub(super) struct ModificationEnumeration<'a, 'p> {
    peptide: &'p Peptide,
    modifications: &'p [ModificationCandidate],
    plan: &'p ModificationPlan<'a>,
    occupied: FnvHashSet<Site>,
    pub(super) selected: Vec<ModificationCandidate>,
    new_counts: Vec<usize>,
    total_counts: Vec<usize>,
    pub(super) output: Vec<Peptide>,
}

impl<'a, 'p> ModificationEnumeration<'a, 'p> {
    pub(super) fn new(
        peptide: &'p Peptide,
        modifications: &'p [ModificationCandidate],
        plan: &'p ModificationPlan<'a>,
        group_count: usize,
    ) -> Self {
        Self {
            peptide,
            modifications,
            plan,
            occupied: FnvHashSet::default(),
            selected: Vec::new(),
            new_counts: vec![0; group_count],
            total_counts: vec![0; group_count],
            output: vec![peptide.clone()],
        }
    }

    pub(super) fn enumerate(
        &mut self,
        start: usize,
        remaining: usize,
        exhaustive_count: usize,
    ) -> bool {
        for idx in start..self.modifications.len() {
            let modification = self.modifications[idx];
            if !self.occupied.insert(modification.site) {
                continue;
            }

            let rule = &self.plan.variable_mods[modification.rule];
            let group = rule.count_group;
            let is_new = !modification.library_supported;
            let next_exhaustive = exhaustive_count + usize::from(is_new);
            let exceeds_limit = next_exhaustive > self.plan.max_exhaustive_mods
                || (is_new
                    && rule
                        .max_count
                        .is_some_and(|limit| self.new_counts[group] >= limit))
                || rule
                    .total_limit()
                    .is_some_and(|limit| self.total_counts[group] >= limit);
            if exceeds_limit {
                self.occupied.remove(&modification.site);
                continue;
            }

            self.new_counts[group] += usize::from(is_new);
            self.total_counts[group] += 1;
            self.selected.push(modification);
            let should_continue = if remaining == 1 {
                self.push_variant()
            } else {
                self.enumerate(idx + 1, remaining - 1, next_exhaustive)
            };
            self.selected.pop();
            self.new_counts[group] -= usize::from(is_new);
            self.total_counts[group] -= 1;
            self.occupied.remove(&modification.site);
            if !should_continue {
                return false;
            }
        }
        true
    }

    fn push_variant(&mut self) -> bool {
        if self
            .plan
            .max_combinations
            .is_some_and(|cap| self.output.len() >= cap)
        {
            return false;
        }
        let mut modified = self.peptide.clone();
        for selected in self.selected.iter().copied() {
            modified.apply_site(
                selected.site,
                self.plan.variable_mods[selected.rule].modification.clone(),
            );
        }
        self.output.push(modified);
        true
    }
}
