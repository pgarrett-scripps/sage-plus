//! MS-cleavable crosslinkers and the `crosslink` configuration block.

use sage_core::mass::{Tolerance, H2O};
use serde::{Deserialize, Serialize};

/// An MS-cleavable crosslinker. Masses are monoisotopic and neutral.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct CleavableLinker {
    pub name: String,
    /// Mass added when two residues are linked.
    pub crosslink_mass: f32,
    /// Stub left on a chain by the lighter cleavage product.
    pub light_stub: f32,
    /// Stub left on a chain by the heavier cleavage product. The light/heavy
    /// spacing defines the signature doublet.
    pub heavy_stub: f32,
    /// Further stubs matched on fragments but not used to find doublets, for
    /// example the DSSO sulfenic acid (thiol + water).
    #[serde(default)]
    pub extra_stubs: Vec<f32>,
}

impl CleavableLinker {
    /// DSSO. Stubs are alkene (C3H2O) and unsaturated thiol (C3H2OS); the
    /// sulfenic acid stub (C3H4O2S) is the thiol plus water.
    pub fn dsso() -> Self {
        Self {
            name: "DSSO".into(),
            crosslink_mass: 158.003_77,
            light_stub: 54.010_56,
            heavy_stub: 85.982_63,
            extra_stubs: vec![85.982_63 + H2O],
        }
    }

    /// DSBU (BuUrBu). Stubs are Bu (C4H7NO) and BuUr (C5H5NO2).
    pub fn dsbu() -> Self {
        Self {
            name: "DSBU".into(),
            crosslink_mass: 196.084_8,
            light_stub: 85.052_76,
            heavy_stub: 111.032_03,
            extra_stubs: Vec::new(),
        }
    }

    /// Spacing between the two members of a signature doublet.
    pub fn doublet_spacing(&self) -> f32 {
        self.heavy_stub - self.light_stub
    }

    /// Every stub a released chain can carry on its fragments.
    pub fn fragment_stubs(&self) -> Vec<f32> {
        let mut stubs = vec![self.light_stub, self.heavy_stub];
        stubs.extend(&self.extra_stubs);
        stubs
    }
}

/// A linker named by preset (`"DSSO"`, `"DSBU"`) or given explicitly.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(untagged)]
pub enum LinkerConfig {
    Preset(String),
    Custom(CleavableLinker),
}

impl Default for LinkerConfig {
    fn default() -> Self {
        Self::Preset("DSSO".into())
    }
}

impl LinkerConfig {
    pub fn resolve(&self) -> Result<CleavableLinker, String> {
        match self {
            Self::Preset(name) => match name.to_ascii_uppercase().as_str() {
                "DSSO" => Ok(CleavableLinker::dsso()),
                "DSBU" | "BUURBU" => Ok(CleavableLinker::dsbu()),
                _ => Err(format!(
                    "unknown crosslinker `{name}`: use DSSO, DSBU, or give crosslink_mass, light_stub and heavy_stub"
                )),
            },
            Self::Custom(linker) => {
                let finite = [linker.crosslink_mass, linker.light_stub, linker.heavy_stub]
                    .iter()
                    .chain(&linker.extra_stubs)
                    .all(|mass| mass.is_finite() && *mass > 0.0);
                if !finite || linker.heavy_stub <= linker.light_stub {
                    return Err(
                        "crosslinker masses must be positive, with heavy_stub > light_stub".into(),
                    );
                }
                Ok(linker.clone())
            }
        }
    }
}

/// The `"crosslink"` block of a search configuration.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(default, deny_unknown_fields)]
pub struct CrosslinkSettings {
    pub linker: LinkerConfig,
    /// Residues the linker reacts with. A residue at the peptide C-terminus is
    /// only linkable at the protein C-terminus, since a linked lysine is not
    /// cleaved by trypsin.
    pub residues: String,
    /// Whether protein N-termini are linkable.
    pub protein_n_term: bool,
    /// Precursor tolerance for chain lookups; defaults to the search
    /// `precursor_tol`.
    pub precursor_tol: Option<Tolerance>,
    /// Precursor isotope errors tried when one chain mass is inferred from the
    /// precursor.
    pub isotope_errors: (i8, i8),
    /// Precursor charges tried when a spectrum has none.
    pub missing_charges: (u8, u8),
    /// Half width, in m/z, of the precursor window a pair of observed chains
    /// must fall in when the precursor mass itself does not match (for example
    /// a misassigned monoisotopic peak). Used when the spectrum does not
    /// record its isolation window.
    pub isolation_half_width: f32,
    /// Lightest chain considered.
    pub min_chain_mass: f32,
    /// Chain-mass pair hypotheses scored per spectrum.
    pub max_pairs: usize,
    /// Chain candidates fully scored per chain mass, from preliminary matches.
    pub preliminary_candidates: usize,
    /// Scored chain candidates per chain combined into crosslinks.
    pub chain_candidates: usize,
    /// Minimum matched fragments for each chain.
    pub min_chain_matched_peaks: u16,
    /// Crosslink-spectrum matches with a CSM q-value above this are not
    /// written.
    pub output_q_value: f32,
}

impl Default for CrosslinkSettings {
    fn default() -> Self {
        Self {
            linker: LinkerConfig::default(),
            residues: "K".into(),
            protein_n_term: true,
            precursor_tol: None,
            isotope_errors: (-1, 3),
            missing_charges: (3, 6),
            isolation_half_width: 1.0,
            min_chain_mass: 400.0,
            max_pairs: 12,
            preliminary_candidates: 10,
            chain_candidates: 5,
            min_chain_matched_peaks: 2,
            output_q_value: 1.0,
        }
    }
}

impl CrosslinkSettings {
    pub fn validate(&self) -> Result<(), String> {
        self.linker.resolve()?;
        if self.residues.is_empty() && !self.protein_n_term {
            return Err("crosslink.residues is empty and protein_n_term is false".into());
        }
        if let Some(residue) = self.residues.bytes().find(|r| !r.is_ascii_uppercase()) {
            return Err(format!(
                "crosslink.residues must be upper-case residues, found `{}`",
                residue as char
            ));
        }
        if self.isotope_errors.0 > self.isotope_errors.1 {
            return Err("crosslink.isotope_errors must be ordered [low, high]".into());
        }
        if self.missing_charges.0 == 0 || self.missing_charges.0 > self.missing_charges.1 {
            return Err("crosslink.missing_charges must be ordered positive charges".into());
        }
        if self.max_pairs == 0 || self.preliminary_candidates == 0 || self.chain_candidates == 0 {
            return Err(
                "crosslink max_pairs, preliminary_candidates and chain_candidates must be positive"
                    .into(),
            );
        }
        if !(self.isolation_half_width >= 0.0 && self.min_chain_mass >= 0.0) {
            return Err("crosslink isolation_half_width and min_chain_mass must be >= 0".into());
        }
        if !(0.0..=1.0).contains(&self.output_q_value) {
            return Err("crosslink.output_q_value must be between 0 and 1".into());
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn stub_masses_are_consistent_with_linker_mass() {
        let dsso = CleavableLinker::dsso();
        // alkene + sulfenic acid = crosslink.
        assert!((dsso.light_stub + dsso.extra_stubs[0] - dsso.crosslink_mass).abs() < 1e-3);
        let dsbu = CleavableLinker::dsbu();
        assert!((dsbu.light_stub + dsbu.heavy_stub - dsbu.crosslink_mass).abs() < 1e-3);
    }

    #[test]
    fn settings_parse_presets_and_custom_linkers() {
        let settings: CrosslinkSettings = serde_json::from_str(r#"{"linker": "dsbu"}"#).unwrap();
        assert_eq!(settings.linker.resolve().unwrap().name, "DSBU");
        assert_eq!(settings.residues, "K");

        let settings: CrosslinkSettings = serde_json::from_str(
            r#"{"linker": {"name": "X", "crosslink_mass": 100.0, "light_stub": 40.0, "heavy_stub": 60.0}}"#,
        )
        .unwrap();
        assert_eq!(settings.linker.resolve().unwrap().doublet_spacing(), 20.0);

        let unknown: CrosslinkSettings = serde_json::from_str(r#"{"linker": "BS3"}"#).unwrap();
        assert!(unknown.validate().is_err());
        assert!(serde_json::from_str::<CrosslinkSettings>(r#"{"bogus": 1}"#).is_err());
    }
}
