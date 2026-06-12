use std::{
    collections::HashMap,
    fmt::{Display, Write},
    str::FromStr,
};

use serde::{Deserialize, Serialize};

use crate::mass::VALID_AA;
use crate::unimod;

#[derive(Copy, Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum ModificationSpecificity {
    PeptideN(Option<u8>),
    PeptideC(Option<u8>),
    ProteinN(Option<u8>),
    ProteinC(Option<u8>),
    Residue(u8),
}

impl Display for ModificationSpecificity {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let r = match self {
            ModificationSpecificity::PeptideN(r) => {
                f.write_char('^')?;
                *r
            }
            ModificationSpecificity::PeptideC(r) => {
                f.write_char('$')?;
                *r
            }
            ModificationSpecificity::ProteinN(r) => {
                f.write_char('[')?;
                *r
            }
            ModificationSpecificity::ProteinC(r) => {
                f.write_char(']')?;
                *r
            }
            ModificationSpecificity::Residue(r) => Some(*r),
        };

        if let Some(r) = r {
            f.write_char(r as char)?;
        }

        Ok(())
    }
}

impl Serialize for ModificationSpecificity {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: serde::Serializer,
    {
        serializer.serialize_str(&self.to_string())
    }
}

#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum InvalidModification {
    Empty,
    InvalidResidue(char),
    TooLong(String),
}

impl FromStr for ModificationSpecificity {
    type Err = InvalidModification;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        if s.len() > 2 {
            return Err(InvalidModification::TooLong(s.into()));
        }
        if let Some(rest) = s.strip_prefix('^') {
            return Ok(ModificationSpecificity::PeptideN(
                rest.chars().next().map(|ch| ch as u8),
            ));
        }
        if let Some(rest) = s.strip_prefix('$') {
            return Ok(ModificationSpecificity::PeptideC(
                rest.chars().next().map(|ch| ch as u8),
            ));
        }
        if let Some(rest) = s.strip_prefix('[') {
            return Ok(ModificationSpecificity::ProteinN(
                rest.chars().next().map(|ch| ch as u8),
            ));
        }
        if let Some(rest) = s.strip_prefix(']') {
            return Ok(ModificationSpecificity::ProteinC(
                rest.chars().next().map(|ch| ch as u8),
            ));
        }
        match s.chars().next() {
            Some(c) => {
                if VALID_AA.contains(&(c as u8)) {
                    Ok(ModificationSpecificity::Residue(c as u8))
                } else {
                    Err(InvalidModification::InvalidResidue(c))
                }
            }
            None => Err(InvalidModification::Empty),
        }
    }
}

/// User-facing modification value: either a numeric monoisotopic delta mass
/// or a Unimod name (e.g. `"Oxidation"`). Name strings are resolved against
/// the embedded Unimod table at validation time.
#[derive(Clone, Debug, Deserialize)]
#[serde(untagged)]
pub enum MassOrName {
    Mass(f32),
    Name(String),
}

impl MassOrName {
    /// Resolve to a numeric mass. When the input is a name, register the
    /// (canonical) Unimod label so that downstream output renders the name
    /// instead of the raw mass. Returns `None` for unknown names.
    fn resolve(&self) -> Option<f32> {
        match self {
            MassOrName::Mass(m) => Some(*m),
            MassOrName::Name(name) => {
                let mass = unimod::mass_by_name(name)?;
                let canonical = unimod::canonical_name(name).unwrap_or(name);
                unimod::register_label(mass, canonical);
                Some(mass)
            }
        }
    }
}

fn log_spec_err(action: &str, err: InvalidModification) {
    match err {
        InvalidModification::Empty => log::error!("{action}: empty modification string"),
        InvalidModification::InvalidResidue(c) => {
            log::error!("{action}: unrecognized residue ({c})")
        }
        InvalidModification::TooLong(s) => log::error!("{action}: {s} is too long"),
    }
}

pub fn validate_mods(
    input: Option<HashMap<String, MassOrName>>,
) -> HashMap<ModificationSpecificity, f32> {
    let mut output = HashMap::new();
    if let Some(input) = input {
        for (s, value) in input {
            match ModificationSpecificity::from_str(&s) {
                Ok(m) => match value.resolve() {
                    Some(mass) => {
                        output.insert(m, mass);
                    }
                    None => log::error!(
                        "Unknown Unimod modification name supplied for static mod on {}",
                        s
                    ),
                },
                Err(e) => log_spec_err("Invalid modification string", e),
            }
        }
    }
    output
}

pub fn validate_var_mods(
    input: Option<HashMap<String, Vec<MassOrName>>>,
) -> HashMap<ModificationSpecificity, Vec<f32>> {
    let mut output = HashMap::new();
    if let Some(input) = input {
        for (s, values) in input {
            match ModificationSpecificity::from_str(&s) {
                Ok(m) => {
                    let resolved: Vec<f32> = values
                        .iter()
                        .filter_map(|v| {
                            let r = v.resolve();
                            if r.is_none() {
                                if let MassOrName::Name(n) = v {
                                    log::error!(
                                        "Unknown Unimod modification name '{}' for variable mod on {}",
                                        n, s
                                    );
                                }
                            }
                            r
                        })
                        .collect();
                    if !resolved.is_empty() {
                        output.insert(m, resolved);
                    }
                }
                Err(e) => log_spec_err("Skipping invalid modification string", e),
            }
        }
    }
    output
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn massorname_resolves_names_and_passes_masses_through() {
        // Numeric form is preserved as-is.
        let m = MassOrName::Mass(15.99491).resolve().unwrap();
        assert!((m - 15.99491).abs() < 1e-5);

        // Name resolves via Unimod and registers a label for output.
        let oxi = MassOrName::Name("Oxidation".into()).resolve().unwrap();
        assert!((oxi - 15.994915).abs() < 1e-4);
        assert_eq!(
            crate::unimod::label_for(oxi).as_deref(),
            Some("Oxidation")
        );

        // Unknown name → None (caller logs and skips).
        assert!(MassOrName::Name("not-a-real-mod".into()).resolve().is_none());
    }

    #[test]
    fn validate_mods_accepts_mixed_input() {
        let mut input: HashMap<String, MassOrName> = HashMap::new();
        input.insert("M".into(), MassOrName::Name("Oxidation".into()));
        input.insert("C".into(), MassOrName::Mass(57.0216));
        let out = validate_mods(Some(input));
        assert_eq!(out.len(), 2);
        assert!((out[&ModificationSpecificity::Residue(b'C')] - 57.0216).abs() < 1e-5);
        assert!(
            (out[&ModificationSpecificity::Residue(b'M')] - 15.994915).abs() < 1e-4
        );
    }

    #[test]
    fn parse_modifications() {
        use InvalidModification::*;
        use ModificationSpecificity::*;
        assert_eq!("[".parse::<ModificationSpecificity>(), Ok(ProteinN(None)));
        assert_eq!(
            "[M".parse::<ModificationSpecificity>(),
            Ok(ProteinN(Some(b'M')))
        );
        assert_eq!(
            "]M".parse::<ModificationSpecificity>(),
            Ok(ProteinC(Some(b'M')))
        );
        assert_eq!("M".parse::<ModificationSpecificity>(), Ok(Residue(b'M')));
        assert_eq!(
            "Z".parse::<ModificationSpecificity>(),
            Err(InvalidResidue('Z'))
        );
    }
}
