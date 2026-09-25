//! The `"glyco": {...}` configuration block.

use serde::{Deserialize, Serialize};

use crate::composition::{n_glycan_composition_space, GlycanComposition, GlycanLibrary};

/// Residue mass of HexNAc, the innermost glycan residue kept on fragments.
pub const HEXNAC: f64 = 203.079_373;

/// Accepted `glyco.site_fragment_forms`.
pub const SITE_FRAGMENT_FORMS: [&str; 4] = ["hexnac", "bare", "hexnac_fuc", "hexnac2"];

/// Accepted `glyco.twin_feature`.
pub const TWIN_FEATURES: [&str; 3] = ["off", "any", "opposite"];

/// Intact N-glycopeptide search settings.
///
/// Every field is optional. With no glycan list at all, a built-in
/// mammalian N-glycan composition space is searched.
#[derive(Clone, Debug, Serialize, Deserialize, PartialEq)]
#[serde(default, deny_unknown_fields)]
pub struct GlycoConfig {
    /// Glycan composition list files, one composition per line
    /// (`HexNAc(4)Hex(5)Fuc(1)NeuAc(2)` or `N(4)H(5)F(1)A(2)`).
    pub glycan_files: Vec<String>,
    /// Extra compositions given inline.
    pub glycans: Vec<String>,
    /// Add HexNAc(2)Hex(3..=20), the high-mannose and paucimannose series,
    /// including the long Hex chains of yeast glycans.
    pub high_mannose: bool,
    /// Sequon motif that carries the glycan.
    pub sequon: String,
    /// Maximum number of noncovalent ammonium (NH3) adducts on the precursor.
    pub ammonium_adducts: u8,
    /// Precursor isotope errors considered when explaining the glycan mass.
    pub isotope_errors: (i8, i8),
    /// Minimum number of oxonium ions, besides HexNAc 204.087, needed for a
    /// spectrum to be searched as a glycopeptide.
    pub min_oxonium_ions: usize,
    /// Glycan FDR threshold for `passes_fdr` in the output.
    pub glycan_fdr: f32,
    /// Peptide FDR threshold for `passes_fdr` in the output.
    pub peptide_fdr: f32,
    /// Peptide candidates kept per spectrum before glycan explanation. Most
    /// top-ranked sequon peptides have no glycan that explains the precursor
    /// delta, so a deep list is what lets a spectrum reach an explained one.
    pub report_candidates: usize,
    /// Explained peptide candidates kept per spectrum (at most
    /// `report_candidates`). The glycan evidence of each is scored, and the
    /// best-supported one represents the spectrum in FDR. 1 keeps only the
    /// highest-ranked explained candidate.
    pub explain_candidates: usize,
    /// Fragment index bucket size of the sequon peptide index.
    pub bucket_size: usize,
    /// Keep `database.variable_mods` in the sequon peptide index. Off by
    /// default: the glyco precursor window spans thousands of Daltons, so
    /// every variable-mod variant is a candidate for every spectrum, which
    /// costs search time and adds peptide decoys to compete with.
    pub variable_mods: bool,
    /// Add the peptide's core Y ions (Y0, Y1, Y1+Fuc and HexNAc(2)Hex(1..=3))
    /// to the fragment index, so that preliminary candidate retrieval counts
    /// them next to b and y ions. Y ions fix the peptide mass, which the open
    /// glycan window leaves free.
    pub index_y_ions: bool,
    /// Forms in which a b or y ion spanning the glycosylation site can
    /// match, counted once per ion: `hexnac` (the innermost HexNAc kept,
    /// required), `bare` (the glycan fully lost), `hexnac_fuc` (the core
    /// fucose kept too) and `hexnac2` (the chitobiose kept).
    pub site_fragment_forms: Vec<String>,
    /// Peaks that match a Y ion (peptide plus a core-like glycan fragment,
    /// at any charge up to the precursor's) or an oxonium ion are not
    /// counted as b/y matches of the candidate.
    pub exclude_glycan_peaks: bool,
    /// Largest b/y fragment charge, below the precursor charge, when the
    /// main `max_fragment_charge` is not set. `null` allows every charge
    /// below the precursor's.
    pub max_fragment_charge: Option<u8>,
    /// Add site-spanning fragment counts to the glyco peptide model.
    pub site_features: bool,
    /// Same-mass competitor ("twin") feature for the glyco peptide model:
    /// the hyperscore gap to the best other peptide candidate with the same
    /// bare mass. `any` compares with every such peptide, `opposite` only
    /// with those of the other target/decoy label (a target's reversed
    /// twin), `off` leaves the feature at zero.
    pub twin_feature: String,
    /// Re-pick each spectrum's candidate among its explained peptides with
    /// the peptide discriminant (b/y, Y-ion, oxonium and glycan evidence
    /// together) instead of the glycan score alone.
    pub rescore_candidates: bool,
    /// Add trunk Hex-ladder features (HexNAc(2)Hex(k) Y ions, their longest
    /// consecutive run, and a control ladder at a fixed offset) to the glyco
    /// peptide model.
    pub ladder_feature: bool,
    /// Add to the glyco peptide model how many other spectra picked the same
    /// peptide (glycoform and charge-state siblings).
    pub sibling_feature: bool,
    /// Fit the glyco peptide model and its q-values separately for
    /// candidates whose glycan has only core Y ions (such as HexNAc(1)).
    pub core_only_subgroup: bool,
}

impl Default for GlycoConfig {
    fn default() -> Self {
        GlycoConfig {
            glycan_files: Vec::new(),
            glycans: Vec::new(),
            high_mannose: false,
            sequon: "motif:N*-{P}-[ST]".into(),
            ammonium_adducts: 1,
            isotope_errors: (0, 1),
            min_oxonium_ions: 1,
            glycan_fdr: 0.01,
            peptide_fdr: 0.01,
            report_candidates: 50,
            explain_candidates: 5,
            bucket_size: 8192,
            variable_mods: false,
            index_y_ions: true,
            site_fragment_forms: vec!["hexnac".into(), "bare".into()],
            exclude_glycan_peaks: true,
            max_fragment_charge: Some(3),
            site_features: true,
            twin_feature: "off".into(),
            rescore_candidates: true,
            ladder_feature: false,
            sibling_feature: false,
            core_only_subgroup: false,
        }
    }
}

impl GlycoConfig {
    /// Parse and check the configuration block.
    pub fn from_value(value: serde_json::Value) -> Result<Self, String> {
        let config: GlycoConfig =
            serde_json::from_value(value).map_err(|error| format!("invalid `glyco`: {error}"))?;
        config.validate()?;
        Ok(config)
    }

    pub fn validate(&self) -> Result<(), String> {
        let (lo, hi) = self.isotope_errors;
        if lo > hi || lo < -1 || hi > 3 {
            return Err("`glyco.isotope_errors` must be ordered and within -1..=3".into());
        }
        if self.ammonium_adducts > 3 {
            return Err("`glyco.ammonium_adducts` must be at most 3".into());
        }
        for (name, value) in [
            ("glyco.glycan_fdr", self.glycan_fdr),
            ("glyco.peptide_fdr", self.peptide_fdr),
        ] {
            if !(value.is_finite() && (0.0..=1.0).contains(&value)) {
                return Err(format!("`{name}` must be between 0 and 1"));
            }
        }
        if self.report_candidates == 0 {
            return Err("`glyco.report_candidates` must be at least 1".into());
        }
        if self.explain_candidates == 0 {
            return Err("`glyco.explain_candidates` must be at least 1".into());
        }
        if self.bucket_size == 0 {
            return Err("`glyco.bucket_size` must be at least 1".into());
        }
        if !self.sequon.starts_with("motif:") {
            return Err("`glyco.sequon` must be a motif, e.g. `motif:N*-{P}-[ST]`".into());
        }
        for glycan in &self.glycans {
            GlycanComposition::parse(glycan)?;
        }
        for form in &self.site_fragment_forms {
            if !SITE_FRAGMENT_FORMS.contains(&form.as_str()) {
                return Err(format!(
                    "`glyco.site_fragment_forms`: unknown form `{form}`, expected one of {SITE_FRAGMENT_FORMS:?}"
                ));
            }
        }
        if !self.site_fragment_forms.iter().any(|form| form == "hexnac") {
            return Err("`glyco.site_fragment_forms` must include `hexnac`".into());
        }
        if !TWIN_FEATURES.contains(&self.twin_feature.as_str()) {
            return Err(format!(
                "`glyco.twin_feature` must be one of {TWIN_FEATURES:?}, not `{}`",
                self.twin_feature
            ));
        }
        if self.max_fragment_charge == Some(0) {
            return Err("`glyco.max_fragment_charge` must be at least 1".into());
        }
        Ok(())
    }

    /// Neutral losses of the HexNAc offset that produce the configured
    /// site-spanning fragment forms. A negative loss is a gain: the form
    /// keeps more than the innermost HexNAc.
    pub fn site_losses(&self) -> Vec<f32> {
        let mut losses: Vec<f32> = self
            .site_fragment_forms
            .iter()
            .filter_map(|form| match form.as_str() {
                "bare" => Some(HEXNAC),
                "hexnac_fuc" => Some(-crate::composition::Monosaccharide::Fuc.mass()),
                "hexnac2" => Some(-HEXNAC),
                _ => None,
            })
            .map(|loss| loss as f32)
            .collect();
        losses.sort_unstable_by(f32::total_cmp);
        losses.dedup();
        losses
    }

    /// Build the composition library from the inline list and the contents
    /// of `glycan_files` (in the same order). Unparseable file lines are
    /// skipped with a warning; an empty result falls back to the built-in
    /// composition space.
    pub fn library(&self, file_contents: &[String]) -> Result<GlycanLibrary, String> {
        let mut compositions = Vec::new();
        for glycan in &self.glycans {
            compositions.push(GlycanComposition::parse(glycan)?);
        }
        let mut skipped = 0;
        for text in file_contents {
            for line in text
                .lines()
                .map(str::trim)
                .filter(|line| !line.is_empty() && !line.starts_with('#'))
            {
                match GlycanComposition::parse(line) {
                    Ok(composition) => compositions.push(composition),
                    Err(error) => {
                        skipped += 1;
                        log::debug!("glyco: skipping glycan `{line}`: {error}");
                    }
                }
            }
        }
        if skipped > 0 {
            log::warn!("glyco: skipped {skipped} unparseable glycan compositions");
        }
        if compositions.is_empty() {
            compositions = n_glycan_composition_space();
        }
        if self.high_mannose {
            for hex in 3..=20 {
                compositions.push(GlycanComposition([2, hex, 0, 0, 0]));
            }
        }
        // A composition must at least carry the HexNAc attached to the sequon.
        compositions.retain(|composition| composition.0[0] >= 1);
        Ok(GlycanLibrary::new(compositions))
    }

    /// JSON for the HexNAc mass offset added to the user's database
    /// settings: fragments carry either nothing or the innermost HexNAc.
    pub fn hexnac_offset(&self) -> serde_json::Value {
        serde_json::json!({
            "mass": HEXNAC,
            "sites": [self.sequon],
            "max_count": 1,
            "search_mode": "mass_offset",
            "neutral_losses": [HEXNAC]
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn defaults_and_validation() {
        let config = GlycoConfig::from_value(serde_json::json!({})).unwrap();
        assert_eq!(config, GlycoConfig::default());
        assert!(GlycoConfig::from_value(serde_json::json!({"isotope_errors": [2, 1]})).is_err());
        assert!(GlycoConfig::from_value(serde_json::json!({"glycans": ["Kdn(1)"]})).is_err());
        assert!(GlycoConfig::from_value(serde_json::json!({"unknown": 1})).is_err());
        assert!(GlycoConfig::from_value(serde_json::json!({"sequon": "N"})).is_err());
        assert!(
            GlycoConfig::from_value(serde_json::json!({"site_fragment_forms": ["bare"]})).is_err()
        );
        assert!(GlycoConfig::from_value(
            serde_json::json!({"site_fragment_forms": ["hexnac", "x"]})
        )
        .is_err());
        let forms = GlycoConfig::from_value(serde_json::json!({
            "site_fragment_forms": ["hexnac", "bare", "hexnac2"]
        }))
        .unwrap();
        assert_eq!(forms.site_losses(), vec![-HEXNAC as f32, HEXNAC as f32]);
        assert_eq!(GlycoConfig::default().site_losses(), vec![HEXNAC as f32]);
        assert!(GlycoConfig::from_value(serde_json::json!({"twin_feature": "decoy"})).is_err());
        for mode in TWIN_FEATURES {
            assert!(GlycoConfig::from_value(serde_json::json!({"twin_feature": mode})).is_ok());
        }
    }

    #[test]
    fn library_merges_sources_and_requires_hexnac() {
        let config = GlycoConfig {
            glycans: vec!["HexNAc(4)Hex(5)".into()],
            high_mannose: true,
            ..Default::default()
        };
        let files = vec!["# comment\nN(2)H(9)\nHex(3)\nnot a glycan\n".to_string()];
        let library = config.library(&files).unwrap();
        // HexNAc(4)Hex(5), HexNAc(2)Hex(3..=20) (which includes N(2)H(9));
        // Hex(3) has no HexNAc and is dropped.
        assert_eq!(library.len(), 1 + 18);
        let default = GlycoConfig::default().library(&[]).unwrap();
        assert!(default.len() > 100);
    }

    #[test]
    fn hexnac_offset_is_a_valid_modification() {
        let builder: sage_core::database::Builder = serde_json::from_value(serde_json::json!({
            "variable_mods": {"HexNAc": GlycoConfig::default().hexnac_offset()}
        }))
        .unwrap();
        builder.validate_modification_keys().unwrap();
        let parameters = builder.make_parameters();
        let offsets = parameters.mass_offset_modifications();
        assert_eq!(offsets.len(), 1);
        assert!((offsets[0].mass() as f64 - HEXNAC).abs() < 1e-3);
    }
}
