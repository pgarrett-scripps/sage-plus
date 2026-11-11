//! The configuration wire format: hand-written Deserialize, Serialize, Display and FromStr.
//!
//! Moved out of `modification.rs` unchanged; every item keeps its `modification::` path.

use std::{
    collections::{BTreeMap, HashMap},
    fmt::{Display, Write},
    str::FromStr,
};

use serde::{
    de::{self, value::MapAccessDeserializer, MapAccess, Visitor},
    Deserialize, Deserializer, Serialize,
};

use super::*;
use crate::site_map::SiteMap;

impl<'de> Deserialize<'de> for StaticModification {
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: Deserializer<'de>,
    {
        #[derive(Deserialize)]
        #[serde(deny_unknown_fields)]
        struct Raw {
            mass: f32,
            #[serde(default)]
            name: Option<String>,
            #[serde(default)]
            neutral_losses: SiteMap<f32>,
            #[serde(default)]
            neutral_loss_mode: NeutralLossMode,
            #[serde(default)]
            channel_offsets: BTreeMap<String, f32>,
            #[serde(default)]
            immonium_ions: SiteMap<f32>,
        }

        let raw = Raw::deserialize(deserializer)?;
        validate_details::<D::Error>(
            raw.mass,
            &raw.name,
            &raw.neutral_losses,
            raw.neutral_loss_mode,
            &raw.channel_offsets,
            &raw.immonium_ions,
        )?;
        Ok(Self {
            mass: raw.mass,
            name: raw.name,
            neutral_losses: raw.neutral_losses,
            neutral_loss_mode: raw.neutral_loss_mode,
            channel_offsets: raw.channel_offsets,
            immonium_ions: raw.immonium_ions,
        })
    }
}

impl<'de> Deserialize<'de> for VariableModification {
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: Deserializer<'de>,
    {
        #[derive(Deserialize)]
        #[serde(deny_unknown_fields)]
        struct Raw {
            mass: f32,
            #[serde(default)]
            max_count: Option<usize>,
            #[serde(default)]
            max_total_count: Option<usize>,
            #[serde(default)]
            name: Option<String>,
            #[serde(default)]
            neutral_losses: SiteMap<f32>,
            #[serde(default)]
            neutral_loss_mode: NeutralLossMode,
            #[serde(default)]
            site_mode: SiteMode,
            #[serde(default)]
            search_mode: SearchMode,
            #[serde(default)]
            channel_offsets: BTreeMap<String, f32>,
            #[serde(default)]
            immonium_ions: SiteMap<f32>,
        }

        let raw = Raw::deserialize(deserializer)?;
        validate_details::<D::Error>(
            raw.mass,
            &raw.name,
            &raw.neutral_losses,
            raw.neutral_loss_mode,
            &raw.channel_offsets,
            &raw.immonium_ions,
        )?;
        if raw.max_total_count == Some(0) {
            return Err(de::Error::custom("max_total_count must be positive"));
        }
        if let (Some(max_count), Some(max_total_count)) = (raw.max_count, raw.max_total_count) {
            if max_total_count < max_count {
                return Err(de::Error::custom(
                    "max_total_count must be at least max_count",
                ));
            }
        }
        if raw.search_mode == SearchMode::MassOffset {
            if raw.mass.abs() < 1e-5 {
                return Err(de::Error::custom(
                    "search_mode `mass_offset` requires a non-zero modification mass",
                ));
            }
            if !raw.channel_offsets.is_empty() {
                return Err(de::Error::custom(
                    "search_mode `mass_offset` does not support channel_offsets",
                ));
            }
            if !raw.immonium_ions.is_empty() {
                return Err(de::Error::custom(
                    "search_mode `mass_offset` does not support immonium_ions: mass-offset PSMs are checked against the unmodified peptide",
                ));
            }
        }
        Ok(Self {
            mass: raw.mass,
            max_count: raw.max_count,
            max_total_count: raw.max_total_count,
            name: raw.name,
            neutral_losses: raw.neutral_losses,
            neutral_loss_mode: raw.neutral_loss_mode,
            site_mode: raw.site_mode,
            search_mode: raw.search_mode,
            channel_offsets: raw.channel_offsets,
            immonium_ions: raw.immonium_ions,
        })
    }
}

impl<'de> Deserialize<'de> for StaticModEntry {
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: Deserializer<'de>,
    {
        struct StaticModEntryVisitor;

        impl<'de> Visitor<'de> for StaticModEntryVisitor {
            type Value = StaticModEntry;

            fn expecting(&self, formatter: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
                formatter.write_str("a modification mass or a structured modification object")
            }

            fn visit_f64<E: de::Error>(self, value: f64) -> Result<Self::Value, E> {
                if !value.is_finite() || value.abs() > f32::MAX as f64 {
                    return Err(E::custom("static modification mass is out of range"));
                }
                Ok(StaticModEntry::Mass(value as f32))
            }

            fn visit_f32<E: de::Error>(self, value: f32) -> Result<Self::Value, E> {
                self.visit_f64(value as f64)
            }

            fn visit_i64<E: de::Error>(self, value: i64) -> Result<Self::Value, E> {
                self.visit_f64(value as f64)
            }

            fn visit_u64<E: de::Error>(self, value: u64) -> Result<Self::Value, E> {
                self.visit_f64(value as f64)
            }

            fn visit_map<A: MapAccess<'de>>(self, map: A) -> Result<Self::Value, A::Error> {
                StaticModification::deserialize(MapAccessDeserializer::new(map))
                    .map(StaticModEntry::Detailed)
            }
        }

        deserializer.deserialize_any(StaticModEntryVisitor)
    }
}

impl<'de> Deserialize<'de> for VarModEntry {
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: Deserializer<'de>,
    {
        struct VarModEntryVisitor;

        impl VarModEntryVisitor {
            fn mass<E>(value: f64) -> Result<VarModEntry, E>
            where
                E: de::Error,
            {
                if !value.is_finite() || value.abs() > f32::MAX as f64 {
                    return Err(E::custom("variable modification mass is out of range"));
                }
                Ok(VarModEntry::Mass(value as f32))
            }
        }

        impl<'de> Visitor<'de> for VarModEntryVisitor {
            type Value = VarModEntry;

            fn expecting(&self, formatter: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
                formatter.write_str("a modification mass or a structured modification object")
            }

            fn visit_f64<E>(self, value: f64) -> Result<Self::Value, E>
            where
                E: de::Error,
            {
                Self::mass(value)
            }

            fn visit_f32<E>(self, value: f32) -> Result<Self::Value, E>
            where
                E: de::Error,
            {
                Self::mass(value as f64)
            }

            fn visit_i64<E>(self, value: i64) -> Result<Self::Value, E>
            where
                E: de::Error,
            {
                Self::mass(value as f64)
            }

            fn visit_u64<E>(self, value: u64) -> Result<Self::Value, E>
            where
                E: de::Error,
            {
                Self::mass(value as f64)
            }

            fn visit_map<A>(self, map: A) -> Result<Self::Value, A::Error>
            where
                A: MapAccess<'de>,
            {
                VariableModification::deserialize(MapAccessDeserializer::new(map))
                    .map(VarModEntry::Detailed)
            }
        }

        deserializer.deserialize_any(VarModEntryVisitor)
    }
}

impl Display for ModificationSpecificity {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        if matches!(
            self,
            Self::PeptideNTerm(_)
                | Self::PeptideCTerm(_)
                | Self::ProteinNTerm(_)
                | Self::ProteinCTerm(_)
                | Self::Motif(_)
        ) {
            return f.write_str(&self.explicit_name());
        }
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
            ModificationSpecificity::Internal(r) => {
                f.write_char('~')?;
                Some(*r)
            }
            _ => unreachable!(),
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
    Motif(String),
}

impl FromStr for ModificationSpecificity {
    type Err = InvalidModification;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        if let Some(pattern) = s.strip_prefix("motif:") {
            return crate::motif::SiteMotif::intern(pattern)
                .map(Self::Motif)
                .map_err(InvalidModification::Motif);
        }
        let (kind, residue) = s.split_once(':').map_or((s, None), |(k, r)| (k, Some(r)));
        let aa = match residue {
            Some(r) if r.len() == 1 && VALID_AA.contains(&r.as_bytes()[0]) => Some(r.as_bytes()[0]),
            Some(r) if r.chars().count() == 1 => {
                return Err(InvalidModification::InvalidResidue(
                    r.chars().next().unwrap(),
                ))
            }
            Some(_) => return Err(InvalidModification::TooLong(s.into())),
            None => None,
        };
        let explicit = match (kind, aa) {
            ("peptide_n_term", None) => Some(Self::PeptideN(None)),
            ("peptide_c_term", None) => Some(Self::PeptideC(None)),
            ("protein_n_term", None) => Some(Self::ProteinN(None)),
            ("protein_c_term", None) => Some(Self::ProteinC(None)),
            ("peptide_n_term", Some(r)) => Some(Self::PeptideNTerm(r)),
            ("peptide_c_term", Some(r)) => Some(Self::PeptideCTerm(r)),
            ("protein_n_term", Some(r)) => Some(Self::ProteinNTerm(r)),
            ("protein_c_term", Some(r)) => Some(Self::ProteinCTerm(r)),
            ("first_residue", Some(r)) => Some(Self::PeptideN(Some(r))),
            ("last_residue", Some(r)) => Some(Self::PeptideC(Some(r))),
            ("internal_residue", Some(r)) => Some(Self::Internal(r)),
            ("protein_first", Some(r)) => Some(Self::ProteinN(Some(r))),
            ("protein_last", Some(r)) => Some(Self::ProteinC(Some(r))),
            _ => None,
        };
        if let Some(explicit) = explicit {
            return Ok(explicit);
        }
        let bytes = s.as_bytes();
        if bytes.is_empty() {
            return Err(InvalidModification::Empty);
        }
        if bytes.len() > 2 {
            return Err(InvalidModification::TooLong(s.into()));
        }
        let (prefix, residue) = match bytes {
            [b'^' | b'$' | b'[' | b']'] => (bytes[0], None),
            [residue] if VALID_AA.contains(residue) => {
                return Ok(Self::Residue(*residue));
            }
            [prefix @ (b'^' | b'$' | b'[' | b']' | b'~'), residue]
                if VALID_AA.contains(residue) =>
            {
                (*prefix, Some(*residue))
            }
            _ => {
                return Err(InvalidModification::InvalidResidue(
                    s.chars().last().unwrap(),
                ))
            }
        };
        Ok(match (prefix, residue) {
            (b'^', residue) => Self::PeptideN(residue),
            (b'$', residue) => Self::PeptideC(residue),
            (b'[', residue) => Self::ProteinN(residue),
            (b']', residue) => Self::ProteinC(residue),
            (b'~', Some(residue)) => Self::Internal(residue),
            // `~` is only accepted above with a residue, so a residue-less `~`
            // is a malformed key, not an internal invariant: report it instead
            // of unwrapping `None`, which would panic on user config if another
            // single-byte prefix is ever added to the match above.
            (b'~', None) => return Err(InvalidModification::InvalidResidue('~')),
            _ => unreachable!(),
        })
    }
}

pub trait ModMapValue: serde::de::DeserializeOwned {
    fn named(value: serde_json::Value) -> Result<Self, String>;
    fn merge(&mut self, other: Self) -> Result<(), String>;
}

impl ModMapValue for StaticModEntry {
    fn named(value: serde_json::Value) -> Result<Self, String> {
        serde_json::from_value(value).map_err(|e| e.to_string())
    }
    fn merge(&mut self, other: Self) -> Result<(), String> {
        if self.definition() != other.definition() {
            return Err("different static modifications target the same site rule".into());
        }
        Ok(())
    }
}

impl ModMapValue for Vec<VarModEntry> {
    fn named(value: serde_json::Value) -> Result<Self, String> {
        serde_json::from_value(value)
            .map(|v| vec![v])
            .map_err(|e| e.to_string())
    }
    fn merge(&mut self, other: Self) -> Result<(), String> {
        self.extend(other);
        Ok(())
    }
}

/// Symbol-keyed maps accept only upstream Sage syntax: a residue or a terminal
/// symbol (`^ $ [ ]`, optionally followed by a residue) mapped to masses. The
/// Sage Plus-only extensions that named definitions replaced in Beta 6 (`~K`
/// and explicit-site keys, and object values) are rejected.
fn reject_legacy_extension(
    id: &str,
    value: &serde_json::Value,
    specificity: ModificationSpecificity,
) -> Result<(), String> {
    const SEE: &str = "See \"Modifications\" in DOCS.md";
    let site = specificity.explicit_name();
    let upstream_key = matches!(id.as_bytes(), [_] | [b'^' | b'$' | b'[' | b']', _]);
    if !upstream_key {
        return Err(format!(
            "modification key `{id}` is no longer supported in symbol-keyed maps. \
             Use a named definition, e.g. {{\"MyMod\": {{\"mass\": 42.0106, \"sites\": [\"{site}\"]}}}}. {SEE}"
        ));
    }
    let has_object = match value {
        serde_json::Value::Object(_) => true,
        serde_json::Value::Array(items) => items.iter().any(serde_json::Value::is_object),
        _ => false,
    };
    if has_object {
        return Err(format!(
            "symbol key `{id}` must map to a mass, not an object. Move `name`, `max_count`, \
             neutral losses and other fields into a named definition, e.g. \
             {{\"MyMod\": {{\"mass\": 42.0106, \"sites\": [\"{site}\"]}}}}. {SEE}"
        ));
    }
    Ok(())
}

/// Normalize named definitions and upstream symbol-keyed configurations into one engine model.
pub fn deserialize_mod_map<'de, D, T>(
    deserializer: D,
) -> Result<Option<HashMap<String, T>>, D::Error>
where
    D: Deserializer<'de>,
    T: ModMapValue,
{
    let Some(input) = Option::<BTreeMap<String, serde_json::Value>>::deserialize(deserializer)?
    else {
        return Ok(None);
    };
    let named = input.values().any(|value| value.get("sites").is_some());
    let mut result: HashMap<String, T> = HashMap::new();
    for (id, mut value) in input {
        if !named {
            if id.starts_with("motif:") {
                return Err(de::Error::custom(format!(
                    "motif site `{id}` requires a named definition, e.g. {{\"HexNAc\": {{\"mass\": 203.079373, \"sites\": [\"{id}\"]}}}}"
                )));
            }
            let specificity = id.parse::<ModificationSpecificity>().map_err(|_| de::Error::custom(format!("invalid modification key `{id}`. Named definitions require a nonempty `sites` array")))?;
            reject_legacy_extension(&id, &value, specificity).map_err(de::Error::custom)?;
            let entry: T = serde_json::from_value(value).map_err(de::Error::custom)?;
            let key = specificity.to_string();
            if let Some(existing) = result.get_mut(&key) {
                existing.merge(entry).map_err(de::Error::custom)?
            } else {
                result.insert(key, entry);
            }
            continue;
        }
        if id.trim().is_empty() || id.trim() != id || id.chars().any(char::is_control) {
            return Err(de::Error::custom("modification IDs must be nonempty and have no surrounding whitespace or control characters"));
        }
        let object = value.as_object_mut().ok_or_else(|| {
            de::Error::custom(
                "named and legacy modification declarations cannot be mixed in one section",
            )
        })?;
        let sites: Vec<String> =
            serde_json::from_value(object.remove("sites").ok_or_else(|| {
                de::Error::custom(format!("modification `{id}` requires `sites`"))
            })?)
            .map_err(de::Error::custom)?;
        if sites.is_empty() {
            return Err(de::Error::custom(format!(
                "modification `{id}` has no sites"
            )));
        }
        if let Some(immonium_ions) = object.get("immonium_ions") {
            // Shape and value errors are reported when the entry is parsed.
            if let Ok(immonium_ions) = SiteMap::<f32>::deserialize(immonium_ions) {
                immonium_ions
                    .validate_sites("immonium_ions", &id, &sites)
                    .map_err(de::Error::custom)?;
            }
        }
        if let Some(name) = object.get("name") {
            if name.as_str() != Some(id.as_str()) {
                return Err(de::Error::custom(format!(
                    "modification `{id}` must not override its identity with a different name"
                )));
            }
        }
        object.insert("name".into(), serde_json::Value::String(id.clone()));
        validate_site_losses(&id, &sites, object).map_err(de::Error::custom)?;
        if object.get("max_count").and_then(|v| v.as_u64()) == Some(0) {
            return Err(de::Error::custom("max_count must be positive"));
        }
        let mut seen = std::collections::HashSet::new();
        for site in sites {
            let specificity = site.parse::<ModificationSpecificity>().map_err(|error| {
                let detail = match error {
                    InvalidModification::Motif(detail) => format!(": {detail}"),
                    _ => String::new(),
                };
                de::Error::custom(format!(
                    "invalid site `{site}` for modification `{id}`{detail}"
                ))
            })?;
            if specificity.explicit_name() != site {
                return Err(de::Error::custom(format!(
                    "use explicit site `{}` instead of `{site}`",
                    specificity.explicit_name()
                )));
            }
            if !seen.insert(specificity) {
                continue;
            }
            let entry = T::named(value.clone()).map_err(de::Error::custom)?;
            let key = specificity.to_string();
            if let Some(existing) = result.get_mut(&key) {
                existing.merge(entry).map_err(de::Error::custom)?
            } else {
                result.insert(key, entry);
            }
        }
    }
    Ok(Some(result))
}

pub(super) trait NamedEntry {
    fn json(&self) -> serde_json::Value;
}
impl NamedEntry for StaticModEntry {
    fn json(&self) -> serde_json::Value {
        serde_json::to_value(self).unwrap()
    }
}
impl NamedEntry for VarModEntry {
    fn json(&self) -> serde_json::Value {
        serde_json::to_value(self).unwrap()
    }
}

pub fn serialize_static_mods<S: serde::Serializer>(
    mods: &HashMap<ModificationSpecificity, StaticModEntry>,
    serializer: S,
) -> Result<S::Ok, S::Error> {
    if mods.values().all(|m| m.definition().name.is_some()) {
        named_modifications(mods.iter().map(|(s, m)| (*s, m as &dyn NamedEntry)))
            .serialize(serializer)
    } else {
        mods.serialize(serializer)
    }
}

pub fn serialize_variable_mods<S: serde::Serializer>(
    mods: &HashMap<ModificationSpecificity, Vec<VarModEntry>>,
    serializer: S,
) -> Result<S::Ok, S::Error> {
    if mods
        .values()
        .flatten()
        .all(|m| m.definition().name.is_some())
    {
        named_modifications(
            mods.iter()
                .flat_map(|(s, entries)| entries.iter().map(move |m| (*s, m as &dyn NamedEntry))),
        )
        .serialize(serializer)
    } else {
        mods.serialize(serializer)
    }
}
