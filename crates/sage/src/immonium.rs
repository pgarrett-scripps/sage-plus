//! Per-PSM immonium-ion evidence.
//!
//! An immonium ion is the single-residue internal fragment `H2N=CH-R+`, with
//! m/z equal to the residue mass minus CO plus a proton. Two kinds are used:
//!
//! - **Residue ions** of unmodified F, Y, W, P, H, V and L/I, whose intensity
//!   tracks the presence of the residue in the peptide (Hohmann et al. 2008,
//!   Anal. Chem. 80:5596). `residue_ions` switches all of them on or off.
//! - **Modified-residue ions** declared on a modification's `immonium_ions`,
//!   such as phosphotyrosine at m/z 216.0420 (Steen et al. 2001, Anal. Chem.
//!   73:1440) or acetyl-lysine at m/z 126.0913 (Trelle & Jensen 2008, Anal.
//!   Chem. 80:3422). There are no built-in modified ions.
//!
//! For each PSM the ions are looked up in the processed spectrum and split
//! into ions explained by the peptide, ions the peptide should give but that
//! are absent, and ions observed whose residue is not in the peptide. The
//! counts are reported and, unless `rescore` is false, added to the linear
//! discriminant. They never enter the fragment index or the hyperscore.

use std::collections::{BTreeMap, BTreeSet, HashMap};

use crate::mass::{monoisotopic, Tolerance, PROTON};
use crate::modification::{ModificationSpecificity, StaticModEntry, VarModEntry};
use crate::peptide::{Peptide, Site};
use crate::spectrum::ProcessedSpectrum;
use serde::{Deserialize, Serialize};

/// Monoisotopic mass of CO, lost from a residue to form its immonium ion.
const CO: f32 = 27.994915;

/// Largest number of modified-residue ions, so observed ions fit a bit mask.
pub const MAX_MODIFIED_IONS: usize = 32;

/// An unnamed modification (for example one parsed from a ProForma mass
/// delta) is the declared modification when its mass is within this many Da.
const MODIFICATION_MASS_TOLERANCE: f32 = 0.01;

/// Residue ions: label and the residues that produce the ion. Leucine and
/// isoleucine give the same ion.
pub const RESIDUE_IONS: [(&str, &[u8]); 7] = [
    ("P", b"P"),
    ("V", b"V"),
    ("L/I", b"LIJ"),
    ("H", b"H"),
    ("F", b"F"),
    ("Y", b"Y"),
    ("W", b"W"),
];

/// Singly charged m/z of the immonium ion of an unmodified residue.
pub fn residue_immonium_mz(residue: u8) -> f32 {
    monoisotopic(residue) - CO + PROTON
}

/// A modified-residue immonium ion, declared on a modification's
/// `immonium_ions`.
#[derive(Debug, Clone, PartialEq, Serialize)]
pub struct ModifiedImmoniumIon {
    /// Label reported in the output: the modification name and its residues,
    /// such as `Phospho@Y` or `Phospho@S/T`.
    pub label: String,
    /// Name of the modification that declares the ion.
    pub modification: String,
    /// Mass of the modification, used only for unnamed modifications.
    #[serde(skip)]
    pub mass: f32,
    /// Residues that carry the modification when it gives this ion.
    #[serde(skip)]
    pub residues: Vec<u8>,
    /// The ion also comes from a site without a fixed residue (a terminus).
    #[serde(skip)]
    pub any_residue: bool,
    /// Singly charged m/z of the ion.
    pub mz: f32,
}

/// Residue a modification rule attaches to, or `None` for a terminus
/// without a residue.
fn rule_residues(specificity: ModificationSpecificity) -> Option<Vec<u8>> {
    use ModificationSpecificity::*;
    match specificity {
        Residue(r) | Internal(r) | PeptideNTerm(r) | PeptideCTerm(r) | ProteinNTerm(r)
        | ProteinCTerm(r) => Some(vec![r]),
        PeptideN(r) | PeptideC(r) | ProteinN(r) | ProteinC(r) => r.map(|r| vec![r]),
        Motif(motif) => Some(motif.site_residues()),
    }
}

/// Modified-residue ions declared by the search's modifications. One ion
/// per modification and m/z; a modification that declares the same m/z on
/// several sites gives one ion whose label lists their residues.
pub fn modified_ions(
    static_mods: &HashMap<ModificationSpecificity, StaticModEntry>,
    variable_mods: &HashMap<ModificationSpecificity, Vec<VarModEntry>>,
) -> Vec<ModifiedImmoniumIon> {
    #[derive(Default)]
    struct Group {
        mass: f32,
        residues: BTreeSet<u8>,
        sites: BTreeSet<String>,
        any_residue: bool,
    }
    let declared = static_mods
        .iter()
        .map(|(site, entry)| (*site, entry.definition(), entry.immonium_ions()))
        .chain(variable_mods.iter().flat_map(|(site, entries)| {
            entries
                .iter()
                .map(move |entry| (*site, entry.definition(), entry.immonium_ions()))
        }));
    let mut groups = BTreeMap::<(String, u32), Group>::new();
    for (specificity, definition, ions) in declared {
        let site = specificity.explicit_name();
        let name = definition
            .name
            .as_deref()
            .map_or_else(|| format!("{:+.4}", definition.mass), str::to_owned);
        for &mz in ions.for_site(&site) {
            let group = groups.entry((name.clone(), mz.to_bits())).or_default();
            group.mass = definition.mass;
            match rule_residues(specificity) {
                Some(residues) => {
                    for residue in residues {
                        group.residues.insert(residue);
                        group.sites.insert((residue as char).to_string());
                    }
                }
                None => {
                    group.any_residue = true;
                    group.sites.insert(site.clone());
                }
            }
        }
    }
    let mut ions = groups
        .into_iter()
        .map(|((modification, mz), group)| ModifiedImmoniumIon {
            label: format!(
                "{modification}@{}",
                group.sites.into_iter().collect::<Vec<_>>().join("/")
            ),
            modification,
            mass: group.mass,
            residues: group.residues.into_iter().collect(),
            any_residue: group.any_residue,
            mz: f32::from_bits(mz),
        })
        .collect::<Vec<_>>();
    ions.sort_by(|a, b| a.mz.total_cmp(&b.mz).then_with(|| a.label.cmp(&b.label)));
    ions
}

/// `immonium` configuration: `true` for the defaults, or an object.
#[derive(Debug, Clone, PartialEq, Deserialize, schemars::JsonSchema)]
#[serde(untagged)]
pub enum ImmoniumConfig {
    Enabled(bool),
    Settings(ImmoniumOptions),
}

#[derive(Debug, Clone, Default, PartialEq, Deserialize, schemars::JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct ImmoniumOptions {
    /// Add the immonium counts to the linear discriminant (default true);
    /// false only reports them.
    pub rescore: Option<bool>,
    /// Global switch for the unmodified-residue ions of P, V, L/I, H, F, Y
    /// and W (default true). Modified-residue ions come from each
    /// modification's `immonium_ions`.
    pub residue_ions: Option<bool>,
    /// Matching tolerance; the search's `fragment_tol` when omitted.
    pub tolerance: Option<Tolerance>,
}

impl ImmoniumConfig {
    /// Settings to use, with the modified-residue ions the modifications
    /// declare; `None` when disabled.
    pub fn resolve(
        self,
        fragment_tol: Tolerance,
        modified: Vec<ModifiedImmoniumIon>,
    ) -> Option<ImmoniumSettings> {
        let options = match self {
            Self::Enabled(false) => return None,
            Self::Enabled(true) => ImmoniumOptions::default(),
            Self::Settings(options) => options,
        };
        Some(ImmoniumSettings {
            rescore: options.rescore.unwrap_or(true),
            residue_ions: options.residue_ions.unwrap_or(true),
            modified,
            tolerance: options.tolerance.unwrap_or(fragment_tol),
        })
    }
}

/// Resolved immonium settings.
#[derive(Debug, Clone, PartialEq, Serialize)]
pub struct ImmoniumSettings {
    pub rescore: bool,
    pub residue_ions: bool,
    /// Modified-residue ions declared on the search's modifications.
    pub modified: Vec<ModifiedImmoniumIon>,
    pub tolerance: Tolerance,
}

impl ImmoniumSettings {
    /// Check the settings; the error names the offending ion.
    pub fn validate(&self) -> Result<(), String> {
        if !self.residue_ions && self.modified.is_empty() {
            return Err(
                "`immonium` is on but has no ions: `residue_ions` is false and no modification declares `immonium_ions`"
                    .into(),
            );
        }
        if self.modified.len() > MAX_MODIFIED_IONS {
            return Err(format!(
                "modifications declare {} immonium ions; at most {MAX_MODIFIED_IONS} are allowed",
                self.modified.len()
            ));
        }
        for ion in &self.modified {
            if ion.label.contains(['\t', '\n', '\r', ',']) {
                return Err(format!(
                    "modification `{}` declares immonium ions, so its name must be free of tabs, newlines and commas",
                    ion.modification
                ));
            }
            if !ion.mz.is_finite() || ion.mz <= 0.0 {
                return Err(format!(
                    "immonium ion `{}` must have a positive m/z",
                    ion.label
                ));
            }
        }
        let (low, high) = self.tolerance.bounds(100.0);
        if !(low.is_finite() && high.is_finite() && low <= 100.0 && 100.0 <= high) {
            return Err("`immonium.tolerance` must contain zero error".into());
        }
        Ok(())
    }

    /// Immonium evidence for `peptide` in `spectrum`.
    pub fn evaluate(&self, spectrum: &ProcessedSpectrum, peptide: &Peptide) -> ImmoniumEvidence {
        let mut evidence = ImmoniumEvidence::default();
        let lowest = spectrum.masses.first().copied().unwrap_or(f32::INFINITY);
        let unmodified = |position: usize| {
            peptide.modification_at(position) == 0.0
                && (position > 0 || peptide.nterm.unwrap_or(0.0) == 0.0)
        };

        if self.residue_ions {
            for (index, (_, residues)) in RESIDUE_IONS.iter().enumerate() {
                let in_peptide = peptide
                    .sequence
                    .iter()
                    .enumerate()
                    .any(|(position, aa)| residues.contains(aa) && unmodified(position));
                let mass = residue_immonium_mz(residues[0]) - PROTON;
                let observed = self.observed(spectrum, mass);
                if observed {
                    evidence.residue_observed |= 1 << index;
                }
                match (in_peptide, observed) {
                    (true, true) => evidence.explained += 1,
                    // An ion below the lowest retained peak may be outside
                    // the scan range; its absence is not evidence.
                    (true, false) if mass >= lowest => evidence.missing += 1,
                    (false, true) => evidence.unexplained += 1,
                    _ => {}
                }
            }
        }

        for (index, ion) in self.modified.iter().enumerate() {
            if self.observed(spectrum, ion.mz - PROTON) {
                evidence.modified_observed |= 1 << index;
                if carries(peptide, ion) {
                    evidence.modified_explained += 1;
                } else {
                    evidence.modified_unexplained += 1;
                }
            }
        }
        evidence
    }

    /// Whether a singly charged peak lies within tolerance of the neutral
    /// `mass` in the processed spectrum.
    fn observed(&self, spectrum: &ProcessedSpectrum, mass: f32) -> bool {
        let (low, high) = self.tolerance.bounds(mass);
        let start = spectrum.masses.partition_point(|&value| value < low);
        (start..spectrum.masses.len())
            .take_while(|&peak| spectrum.masses[peak] <= high)
            .any(|peak| spectrum.charges.get(peak).is_none_or(|&charge| charge <= 1))
    }

    /// Names of the modified ions set in `mask`, comma-separated.
    pub fn modified_names(&self, mask: u32) -> String {
        self.modified
            .iter()
            .enumerate()
            .filter(|(index, _)| mask & (1 << index) != 0)
            .map(|(_, ion)| ion.label.as_str())
            .collect::<Vec<_>>()
            .join(",")
    }
}

/// Does `peptide` carry the modification that gives `ion`, on one of the
/// ion's residues?
fn carries(peptide: &Peptide, ion: &ModifiedImmoniumIon) -> bool {
    let residue_at = |site: Site| match site {
        Site::Sequence(index) => peptide.sequence.get(index as usize).copied(),
        Site::Nterm => peptide.sequence.first().copied(),
        Site::Cterm => peptide.sequence.last().copied(),
    };
    peptide.applied_modifications().any(|applied| {
        let definition = applied.modification;
        let same = match definition.name.as_deref() {
            Some(name) => name == ion.modification,
            None => (definition.mass - ion.mass).abs() <= MODIFICATION_MASS_TOLERANCE,
        };
        same && (ion.any_residue
            || residue_at(applied.site).is_some_and(|residue| ion.residues.contains(&residue)))
    })
}

/// Immonium ions of one PSM.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct ImmoniumEvidence {
    /// Residue ions observed whose unmodified residue is in the peptide.
    pub explained: u8,
    /// Residue ions of the peptide that are not observed, counted only above
    /// the lowest retained peak.
    pub missing: u8,
    /// Residue ions observed whose unmodified residue is not in the peptide.
    pub unexplained: u8,
    /// Modified-residue ions observed whose modified residue is in the peptide.
    pub modified_explained: u8,
    /// Modified-residue ions observed whose modified residue is not in the
    /// peptide.
    pub modified_unexplained: u8,
    /// Bit `i` is set when [`RESIDUE_IONS`]`[i]` is observed.
    pub residue_observed: u8,
    /// Bit `i` is set when modified ion `i` is observed.
    pub modified_observed: u32,
}

/// Number of linear-discriminant columns added by [`ImmoniumEvidence::lda_row`].
pub const LDA_FEATURES: usize = 5;

impl ImmoniumEvidence {
    /// Feature columns for the linear discriminant.
    pub fn lda_row(&self) -> [f64; LDA_FEATURES] {
        [
            self.explained as f64,
            self.missing as f64,
            self.unexplained as f64,
            self.modified_explained as f64,
            self.modified_unexplained as f64,
        ]
    }

    /// Labels of the observed residue ions, comma-separated.
    pub fn residue_names(&self) -> String {
        RESIDUE_IONS
            .iter()
            .enumerate()
            .filter(|(index, _)| self.residue_observed & (1 << index) != 0)
            .map(|(_, (name, _))| *name)
            .collect::<Vec<_>>()
            .join(",")
    }
}

#[cfg(test)]
#[path = "../tests/unit/immonium.rs"]
mod test;
