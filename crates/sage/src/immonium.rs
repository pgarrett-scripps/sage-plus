//! Per-PSM immonium-ion evidence.
//!
//! An immonium ion is the single-residue internal fragment `H2N=CH-R+`, with
//! m/z equal to the residue mass minus CO plus a proton. Two kinds are used:
//!
//! - **Residue ions** of unmodified F, Y, W, P, H, V and L/I, whose intensity
//!   tracks the presence of the residue in the peptide (Hohmann et al. 2008,
//!   Anal. Chem. 80:5596).
//! - **Modified-residue ions**, such as phosphotyrosine at m/z 216.0420
//!   (Steen et al. 2001, Anal. Chem. 73:1440) and acetyl-lysine at m/z
//!   126.0913 (Trelle & Jensen 2008, Anal. Chem. 80:3422).
//!
//! For each PSM the ions are looked up in the processed spectrum and split
//! into ions explained by the peptide, ions the peptide should give but that
//! are absent, and ions observed whose residue is not in the peptide. The
//! counts are reported and, when `rescore` is on, added to the linear
//! discriminant. They never enter the fragment index or the hyperscore.

use crate::mass::{monoisotopic, Tolerance, PROTON};
use crate::peptide::Peptide;
use crate::spectrum::ProcessedSpectrum;
use serde::{Deserialize, Serialize};

/// Monoisotopic mass of CO, lost from a residue to form its immonium ion.
const CO: f32 = 27.994915;

/// Largest number of modified-residue ions, so observed ions fit a bit mask.
pub const MAX_MODIFIED_IONS: usize = 32;

/// Modification masses within this many Da of a configured ion's
/// `modification` are the same modification.
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

/// A modified-residue immonium ion.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize, schemars::JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct ModifiedImmoniumIon {
    /// Label reported in the output.
    pub name: String,
    /// Modified residue, as a one-letter code.
    pub residue: char,
    /// Mass of the modification on that residue, matched within 0.01 Da.
    pub modification: f32,
    /// Singly charged m/z of the ion.
    pub mz: f32,
}

impl ModifiedImmoniumIon {
    fn new(name: &str, residue: char, modification: f32, mz: f32) -> Self {
        Self {
            name: name.into(),
            residue,
            modification,
            mz,
        }
    }
}

/// Built-in modified-residue ions: phosphotyrosine (Steen et al. 2001) and
/// acetyl-lysine after ammonia loss (Trelle & Jensen 2008, 98% specific;
/// the intact 143.118 ion is not specific and is left out).
pub fn default_modified_ions() -> Vec<ModifiedImmoniumIon> {
    vec![
        ModifiedImmoniumIon::new("pY", 'Y', 79.96633, 216.042),
        ModifiedImmoniumIon::new("acK", 'K', 42.010565, 126.0913),
    ]
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
    /// Add the immonium counts to the linear discriminant (default false).
    pub rescore: Option<bool>,
    /// Look for the unmodified F, Y, W, P, H, V and L/I ions (default true).
    pub residues: Option<bool>,
    /// Modified-residue ions; the built-in pY and acK ions when omitted, none
    /// when empty.
    pub modified: Option<Vec<ModifiedImmoniumIon>>,
    /// Matching tolerance; the search's `fragment_tol` when omitted.
    pub tolerance: Option<Tolerance>,
}

impl ImmoniumConfig {
    /// Settings to use; `None` when disabled.
    pub fn resolve(self, fragment_tol: Tolerance) -> Option<ImmoniumSettings> {
        let options = match self {
            Self::Enabled(false) => return None,
            Self::Enabled(true) => ImmoniumOptions::default(),
            Self::Settings(options) => options,
        };
        Some(ImmoniumSettings {
            rescore: options.rescore.unwrap_or(false),
            residues: options.residues.unwrap_or(true),
            modified: options.modified.unwrap_or_else(default_modified_ions),
            tolerance: options.tolerance.unwrap_or(fragment_tol),
        })
    }
}

/// Resolved immonium settings.
#[derive(Debug, Clone, PartialEq, Serialize)]
pub struct ImmoniumSettings {
    pub rescore: bool,
    pub residues: bool,
    pub modified: Vec<ModifiedImmoniumIon>,
    pub tolerance: Tolerance,
}

impl ImmoniumSettings {
    /// Check the modified-ion list; the error names the offending ion.
    pub fn validate(&self) -> Result<(), String> {
        if self.modified.len() > MAX_MODIFIED_IONS {
            return Err(format!(
                "`immonium.modified` allows at most {MAX_MODIFIED_IONS} ions"
            ));
        }
        for ion in &self.modified {
            if ion.name.is_empty() || ion.name.contains(['\t', '\n', '\r', ',']) {
                return Err(
                    "immonium ion names must be non-empty and free of tabs, newlines and commas"
                        .into(),
                );
            }
            if !ion.residue.is_ascii_uppercase() || !monoisotopic(ion.residue as u8).is_normal() {
                return Err(format!(
                    "immonium ion `{}` needs a one-letter amino-acid residue",
                    ion.name
                ));
            }
            if !ion.modification.is_finite() || ion.modification == 0.0 {
                return Err(format!(
                    "immonium ion `{}` needs a finite, non-zero modification mass",
                    ion.name
                ));
            }
            if !ion.mz.is_finite() || ion.mz <= 0.0 {
                return Err(format!(
                    "immonium ion `{}` must have a positive m/z",
                    ion.name
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

        if self.residues {
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
            let residue = ion.residue as u8;
            let in_peptide = peptide.sequence.iter().enumerate().any(|(position, &aa)| {
                let mut mass = peptide.modification_at(position);
                if position == 0 {
                    mass += peptide.nterm.unwrap_or(0.0);
                }
                aa == residue && (mass - ion.modification).abs() <= MODIFICATION_MASS_TOLERANCE
            });
            if self.observed(spectrum, ion.mz - PROTON) {
                evidence.modified_observed |= 1 << index;
                if in_peptide {
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
            .map(|(_, ion)| ion.name.as_str())
            .collect::<Vec<_>>()
            .join(",")
    }
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
