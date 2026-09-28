//! Generic fragment neutral losses, such as water and ammonia.
//!
//! These losses depend on the residues a fragment contains, not on a
//! modification. They are configured under `database.fragment_losses` and are
//! separate from a modification's `neutral_losses`. Loss ions are never part
//! of the preliminary fragment index: they are matched only when a candidate
//! is fully scored, and they feed the rescoring model as their own features
//! rather than the hyperscore (see `benchmarks/FRAGMENT_LOSSES.md`).

use crate::ion_series::Kind;
use crate::modification::ModificationSpecificity;
use crate::peptide::{Peptide, Site};
use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;
use std::sync::Arc;

/// One configured fragment loss. The map key is a label only.
#[derive(Clone, Debug, Serialize, Deserialize, PartialEq, schemars::JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct FragmentLossEntry {
    /// Neutral mass removed from the fragment, in daltons. Positive and finite.
    pub mass: f32,
    /// Where the loss can occur, in the modification site vocabulary (for
    /// example `"S"`, `"first_residue:E"` or `"peptide_c_term"`). A fragment
    /// is eligible when it contains at least one listed site.
    #[schemars(length(min = 1))]
    pub sites: Vec<String>,
    /// Ion kinds that can carry the loss. A subset of `database.ion_kinds`.
    #[schemars(length(min = 1))]
    pub ion_kinds: Vec<Kind>,
    /// Whether a residue or terminus carrying any modification still counts
    /// as a loss site (default false).
    #[serde(default, skip_serializing_if = "std::ops::Not::not")]
    pub allow_modified: bool,
}

/// A validated fragment loss.
#[derive(Clone, Debug, PartialEq)]
pub struct FragmentLoss {
    pub name: Arc<str>,
    pub mass: f32,
    pub sites: Vec<ModificationSpecificity>,
    pub ion_kinds: Vec<Kind>,
    pub allow_modified: bool,
}

/// Validated fragment-loss settings used during full candidate scoring.
#[derive(Clone, Debug, PartialEq)]
pub struct FragmentLosses {
    pub losses: Vec<FragmentLoss>,
    /// Most generic losses stacked on one fragment.
    pub max_losses: usize,
    /// Experiment hook: how matched loss ions enter scoring.
    pub scoring: LossScoring,
}

/// Where matched loss ions enter scoring (experiment hook).
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub enum LossMode {
    /// Hyperscore ignores loss ions.
    Features,
    /// B: loss forms are alternatives of their cleavage (most intense wins).
    Alternatives,
    /// B-int: as B, but a cleavage matched only by a loss form adds its
    /// intensity and no matched-ion count.
    AlternativesIntensity,
    /// B-intact: as B, but the intact form wins whenever it matched.
    AlternativesIntactFirst,
    /// Parent-supported: a loss peak counts only when its intact form
    /// matched at the same cleavage and charge; it adds `weight` to the
    /// hyperscore's matched-ion count and its intensity to the summed
    /// intensity. `matched_peaks` and the other features are unchanged.
    Parent,
    /// Comet-like: every matched loss peak, with or without its intact
    /// form, adds `weight` to the hyperscore's matched-ion count and
    /// `weight` times its intensity to the summed intensity.
    Weighted,
}

/// Experiment hook, selected with `SAGE_PLUS_FRAGMENT_LOSS_*` variables.
#[derive(Copy, Clone, Debug, PartialEq)]
pub struct LossScoring {
    pub mode: LossMode,
    /// Parent mode: matched-ion count weight of a supported loss peak.
    pub weight: f64,
    /// Parent mode: whether a supported loss peak adds its intensity.
    pub intensity: bool,
    /// Emit `matched_loss_peaks` and `loss_intensity_pct` LDA features.
    pub features: bool,
    /// Features count only loss peaks whose intact form matched.
    pub feature_parent: bool,
}

impl LossScoring {
    fn from_env() -> Self {
        let var = |name: &str| std::env::var(name).ok();
        let flag = |name: &str, default: bool| {
            var(name).map_or(default, |value| value == "1" || value == "true")
        };
        let mode = match var("SAGE_PLUS_FRAGMENT_LOSS_SCORING").as_deref() {
            Some("hyperscore") => LossMode::Alternatives,
            Some("hyperscore_intensity") => LossMode::AlternativesIntensity,
            Some("intact_first") => LossMode::AlternativesIntactFirst,
            Some("parent") => LossMode::Parent,
            Some("weighted") => LossMode::Weighted,
            _ => LossMode::Features,
        };
        LossScoring {
            mode,
            weight: var("SAGE_PLUS_FRAGMENT_LOSS_WEIGHT")
                .and_then(|value| value.parse().ok())
                .unwrap_or(1.0),
            intensity: flag("SAGE_PLUS_FRAGMENT_LOSS_PARENT_INTENSITY", true),
            features: flag(
                "SAGE_PLUS_FRAGMENT_LOSS_FEATURES",
                mode == LossMode::Features,
            ),
            feature_parent: flag("SAGE_PLUS_FRAGMENT_LOSS_FEATURE_PARENT", false),
        }
    }
}

/// Default for `database.max_fragment_losses`.
pub const DEFAULT_MAX_FRAGMENT_LOSSES: usize = 1;

/// Validate the `database.fragment_losses` and `database.max_fragment_losses`
/// settings against the searched `ion_kinds`. Returns `None` when fragment
/// losses are not configured.
pub fn resolve(
    entries: Option<&BTreeMap<String, FragmentLossEntry>>,
    max_losses: Option<usize>,
    ion_kinds: &[Kind],
) -> Result<Option<FragmentLosses>, String> {
    let Some(entries) = entries else {
        if max_losses.is_some() {
            return Err(
                "`database.max_fragment_losses` requires `database.fragment_losses`".into(),
            );
        }
        return Ok(None);
    };
    if entries.is_empty() {
        return Err(
            "`database.fragment_losses` must define at least one loss; remove the key to disable fragment losses"
                .into(),
        );
    }
    let max_losses = max_losses.unwrap_or(DEFAULT_MAX_FRAGMENT_LOSSES);
    if max_losses == 0 {
        return Err("`database.max_fragment_losses` must be at least 1".into());
    }
    let mut losses = Vec::with_capacity(entries.len());
    for (name, entry) in entries {
        if name.trim().is_empty() || name.trim() != name || name.chars().any(char::is_control) {
            return Err(
                "fragment loss names must be nonempty and have no surrounding whitespace or control characters"
                    .into(),
            );
        }
        if !entry.mass.is_finite() || entry.mass <= 0.0 {
            return Err(format!(
                "fragment loss `{name}` needs a positive, finite `mass` (got {})",
                entry.mass
            ));
        }
        if entry.sites.is_empty() {
            return Err(format!("fragment loss `{name}` has no sites"));
        }
        let mut sites = Vec::with_capacity(entry.sites.len());
        for site in &entry.sites {
            let specificity = site
                .parse::<ModificationSpecificity>()
                .map_err(|_| format!("invalid site `{site}` for fragment loss `{name}`"))?;
            if specificity.is_motif() {
                return Err(format!(
                    "fragment loss `{name}` cannot use motif site `{site}`; use residue or terminal sites"
                ));
            }
            if specificity.explicit_name() != *site {
                return Err(format!(
                    "use explicit site `{}` instead of `{site}` for fragment loss `{name}`",
                    specificity.explicit_name()
                ));
            }
            if !sites.contains(&specificity) {
                sites.push(specificity);
            }
        }
        if entry.ion_kinds.is_empty() {
            return Err(format!("fragment loss `{name}` has no `ion_kinds`"));
        }
        let mut kinds = Vec::with_capacity(entry.ion_kinds.len());
        for kind in &entry.ion_kinds {
            if !ion_kinds.contains(kind) {
                return Err(format!(
                    "fragment loss `{name}` lists ion kind `{}`, which is not in `database.ion_kinds`",
                    kind_name(*kind)
                ));
            }
            if !kinds.contains(kind) {
                kinds.push(*kind);
            }
        }
        losses.push(FragmentLoss {
            name: Arc::from(name.as_str()),
            mass: entry.mass,
            sites,
            ion_kinds: kinds,
            allow_modified: entry.allow_modified,
        });
    }
    Ok(Some(FragmentLosses {
        losses,
        max_losses,
        scoring: LossScoring::from_env(),
    }))
}

fn kind_name(kind: Kind) -> String {
    serde_json::to_value(kind)
        .ok()
        .and_then(|value| value.as_str().map(str::to_owned))
        .unwrap_or_else(|| format!("{kind:?}"))
}

/// Loss sites of one configured loss on one peptide.
#[derive(Clone, Debug, Default)]
pub(crate) struct LossSites {
    pub mass: f32,
    /// Whether the loss applies to the ion kind being generated.
    pub applies: bool,
    pub nterm: bool,
    pub cterm: bool,
    /// Sorted zero-based residue indices.
    pub residues: Vec<u32>,
}

impl LossSites {
    pub(crate) fn new(loss: &FragmentLoss, peptide: &Peptide, kind: Kind) -> Self {
        let mut sites = LossSites {
            mass: loss.mass,
            applies: loss.ion_kinds.contains(&kind),
            ..Default::default()
        };
        if !sites.applies {
            return sites;
        }
        for specificity in &loss.sites {
            for site in specificity.sites(&peptide.sequence, peptide.position) {
                if !loss.allow_modified && is_modified(peptide, site) {
                    continue;
                }
                match site {
                    Site::Nterm => sites.nterm = true,
                    Site::Cterm => sites.cterm = true,
                    Site::Sequence(index) => sites.residues.push(index),
                }
            }
        }
        sites.residues.sort_unstable();
        sites.residues.dedup();
        sites
    }

    /// Number of loss sites in the fragment of `kind` at `series_index`.
    pub(crate) fn count(&self, kind: Kind, series_index: usize) -> usize {
        if !self.applies {
            return 0;
        }
        // Residues 0..=series_index belong to the N-terminal fragment.
        let split = self
            .residues
            .partition_point(|&index| index as usize <= series_index);
        match kind {
            Kind::A | Kind::B | Kind::C => split + usize::from(self.nterm),
            Kind::X | Kind::Y | Kind::Z | Kind::ZDot => {
                self.residues.len() - split + usize::from(self.cterm)
            }
        }
    }
}

fn is_modified(peptide: &Peptide, site: Site) -> bool {
    match site {
        Site::Nterm if peptide.nterm.is_some() => return true,
        Site::Cterm if peptide.cterm.is_some() => return true,
        _ => {}
    }
    peptide
        .applied_modifications()
        .any(|applied| applied.site == site)
}

/// Totals of every combination of at most `max_losses` generic losses the
/// fragment can carry, each loss used no more often than it has sites.
/// Excludes the empty combination.
pub(crate) fn loss_totals(
    counts: &[(f32, usize)],
    max_losses: usize,
) -> smallvec::SmallVec<[f32; 4]> {
    fn walk(
        counts: &[(f32, usize)],
        remaining: usize,
        total: f32,
        used: usize,
        out: &mut smallvec::SmallVec<[f32; 4]>,
    ) {
        let Some(((mass, available), rest)) = counts.split_first() else {
            if used > 0 {
                out.push(total);
            }
            return;
        };
        for n in 0..=(*available).min(remaining) {
            walk(rest, remaining - n, total + mass * n as f32, used + n, out);
        }
    }
    let mut out = smallvec::SmallVec::new();
    walk(counts, max_losses, 0.0, 0, &mut out);
    out
}

#[cfg(test)]
#[path = "../tests/unit/fragment_loss.rs"]
mod tests;
