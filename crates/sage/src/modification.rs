use std::{
    cmp::Ordering,
    collections::{BTreeMap, HashMap},
    fmt::{Display, Write},
    str::FromStr,
    sync::Arc,
};

use serde::{
    de::{self, value::MapAccessDeserializer, MapAccess, Visitor},
    Deserialize, Deserializer, Serialize,
};

#[derive(
    Copy,
    Clone,
    Debug,
    Default,
    Deserialize,
    Serialize,
    schemars::JsonSchema,
    PartialEq,
    Eq,
    PartialOrd,
    Ord,
)]
#[serde(rename_all = "snake_case")]
pub enum NeutralLossMode {
    #[default]
    Optional,
    Required,
}

/// Controls where candidates for a variable modification are generated.
#[derive(
    Copy, Clone, Debug, Default, Deserialize, Serialize, schemars::JsonSchema, PartialEq, Eq,
)]
#[serde(rename_all = "snake_case")]
pub enum SiteMode {
    /// Generate the modification at every compatible site from the config.
    #[default]
    Exhaustive,
    /// Generate the modification only at sites listed in the PTM library.
    Library,
    /// Generate exhaustive candidates and allow library sites as targeted additions.
    Both,
}

/// Controls how a variable modification enters the search.
#[derive(
    Copy, Clone, Debug, Default, Deserialize, Serialize, schemars::JsonSchema, PartialEq, Eq,
)]
#[serde(rename_all = "snake_case")]
pub enum SearchMode {
    /// Expand modified peptidoforms while building the fragment index.
    #[default]
    Database,
    /// Keep the index unmodified and test at most one copy of this
    /// modification per peptide as a search-time precursor and fragment offset.
    MassOffset,
}

fn is_optional(mode: &NeutralLossMode) -> bool {
    *mode == NeutralLossMode::Optional
}

fn validate_details<E: de::Error>(
    mass: f32,
    name: &Option<String>,
    neutral_losses: &[f32],
    neutral_loss_mode: NeutralLossMode,
    channel_offsets: &BTreeMap<String, f32>,
) -> Result<(), E> {
    if !mass.is_finite() {
        return Err(E::custom("modification mass must be finite"));
    }
    if matches!(name.as_ref(), Some(name) if name.trim().is_empty()) {
        return Err(E::custom("modification name must not be empty"));
    }
    if neutral_losses
        .iter()
        .any(|loss| !loss.is_finite() || *loss <= 0.0)
    {
        return Err(E::custom(
            "neutral loss masses must be finite and greater than zero",
        ));
    }
    if neutral_loss_mode == NeutralLossMode::Required && neutral_losses.is_empty() {
        return Err(E::custom(
            "neutral_loss_mode `required` requires at least one neutral loss",
        ));
    }
    for (channel, offset) in channel_offsets {
        if channel.is_empty() || channel.trim() != channel {
            return Err(E::custom(
                "channel offset names must be non-empty and contain no surrounding whitespace",
            ));
        }
        if channel
            .chars()
            .any(|character| matches!(character, '\r' | '\n' | '=' | ';'))
        {
            return Err(E::custom(format!(
                "channel offset name `{channel}` contains an unsupported character"
            )));
        }
        if !offset.is_finite() {
            return Err(E::custom(format!(
                "channel offset `{channel}` must be finite"
            )));
        }
    }
    Ok(())
}

/// A structured static modification. Numeric static modifications remain
/// supported through [`StaticModEntry::Mass`].
#[derive(Clone, Debug, Serialize, schemars::JsonSchema, PartialEq)]
#[serde(deny_unknown_fields)]
pub struct StaticModification {
    pub mass: f32,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub name: Option<String>,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub neutral_losses: Vec<f32>,
    #[serde(default, skip_serializing_if = "is_optional")]
    pub neutral_loss_mode: NeutralLossMode,
    #[serde(default, skip_serializing_if = "BTreeMap::is_empty")]
    pub channel_offsets: BTreeMap<String, f32>,
}

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
            neutral_losses: Vec<f32>,
            #[serde(default)]
            neutral_loss_mode: NeutralLossMode,
            #[serde(default)]
            channel_offsets: BTreeMap<String, f32>,
        }

        let raw = Raw::deserialize(deserializer)?;
        validate_details::<D::Error>(
            raw.mass,
            &raw.name,
            &raw.neutral_losses,
            raw.neutral_loss_mode,
            &raw.channel_offsets,
        )?;
        Ok(Self {
            mass: raw.mass,
            name: raw.name,
            neutral_losses: raw.neutral_losses,
            neutral_loss_mode: raw.neutral_loss_mode,
            channel_offsets: raw.channel_offsets,
        })
    }
}

/// A variable modification with optional per-peptide occurrence limit and
/// optional fragment neutral-loss behavior.
#[derive(Clone, Debug, Serialize, schemars::JsonSchema, PartialEq)]
#[serde(deny_unknown_fields)]
pub struct VariableModification {
    pub mass: f32,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub max_count: Option<usize>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub name: Option<String>,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub neutral_losses: Vec<f32>,
    #[serde(default, skip_serializing_if = "is_optional")]
    pub neutral_loss_mode: NeutralLossMode,
    #[serde(default, skip_serializing_if = "is_exhaustive")]
    pub site_mode: SiteMode,
    #[serde(default, skip_serializing_if = "is_database")]
    pub search_mode: SearchMode,
    #[serde(default, skip_serializing_if = "BTreeMap::is_empty")]
    pub channel_offsets: BTreeMap<String, f32>,
}

fn is_exhaustive(mode: &SiteMode) -> bool {
    *mode == SiteMode::Exhaustive
}

fn is_database(mode: &SearchMode) -> bool {
    *mode == SearchMode::Database
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
            name: Option<String>,
            #[serde(default)]
            neutral_losses: Vec<f32>,
            #[serde(default)]
            neutral_loss_mode: NeutralLossMode,
            #[serde(default)]
            site_mode: SiteMode,
            #[serde(default)]
            search_mode: SearchMode,
            #[serde(default)]
            channel_offsets: BTreeMap<String, f32>,
        }

        let raw = Raw::deserialize(deserializer)?;
        validate_details::<D::Error>(
            raw.mass,
            &raw.name,
            &raw.neutral_losses,
            raw.neutral_loss_mode,
            &raw.channel_offsets,
        )?;
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
        }
        Ok(Self {
            mass: raw.mass,
            max_count: raw.max_count,
            name: raw.name,
            neutral_losses: raw.neutral_losses,
            neutral_loss_mode: raw.neutral_loss_mode,
            site_mode: raw.site_mode,
            search_mode: raw.search_mode,
            channel_offsets: raw.channel_offsets,
        })
    }
}

#[derive(Clone, Debug)]
pub struct ModificationDefinition {
    pub mass: f32,
    pub name: Option<Arc<str>>,
    pub neutral_losses: Arc<[f32]>,
    pub neutral_loss_mode: NeutralLossMode,
    pub channel_offsets: Arc<BTreeMap<Arc<str>, f32>>,
}

impl ModificationDefinition {
    pub fn bare(mass: f32) -> Self {
        Self {
            mass,
            name: None,
            neutral_losses: Arc::from([]),
            neutral_loss_mode: NeutralLossMode::Optional,
            channel_offsets: Arc::default(),
        }
    }

    fn detailed(
        mass: f32,
        name: &Option<String>,
        neutral_losses: &[f32],
        neutral_loss_mode: NeutralLossMode,
        channel_offsets: &BTreeMap<String, f32>,
    ) -> Self {
        Self {
            mass,
            name: name.as_deref().map(Arc::from),
            neutral_losses: Arc::from(neutral_losses),
            neutral_loss_mode,
            channel_offsets: Arc::new(
                channel_offsets
                    .iter()
                    .map(|(channel, offset)| (Arc::from(channel.as_str()), *offset))
                    .collect(),
            ),
        }
    }

    pub fn with_mass(&self, mass: f32) -> Self {
        Self {
            mass,
            name: self.name.clone(),
            neutral_losses: self.neutral_losses.clone(),
            neutral_loss_mode: self.neutral_loss_mode,
            channel_offsets: self.channel_offsets.clone(),
        }
    }
}

impl PartialEq for ModificationDefinition {
    fn eq(&self, other: &Self) -> bool {
        self.cmp(other) == Ordering::Equal
    }
}

impl Eq for ModificationDefinition {}

impl PartialOrd for ModificationDefinition {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for ModificationDefinition {
    fn cmp(&self, other: &Self) -> Ordering {
        self.mass
            .total_cmp(&other.mass)
            .then_with(|| self.name.cmp(&other.name))
            .then_with(|| {
                self.neutral_losses
                    .iter()
                    .map(|loss| loss.to_bits())
                    .cmp(other.neutral_losses.iter().map(|loss| loss.to_bits()))
            })
            .then_with(|| self.neutral_loss_mode.cmp(&other.neutral_loss_mode))
            .then_with(|| {
                self.channel_offsets
                    .iter()
                    .map(|(channel, offset)| (channel, offset.to_bits()))
                    .cmp(
                        other
                            .channel_offsets
                            .iter()
                            .map(|(channel, offset)| (channel, offset.to_bits())),
                    )
            })
    }
}

#[derive(Clone, Debug, Serialize, schemars::JsonSchema, PartialEq)]
#[serde(untagged)]
pub enum StaticModEntry {
    Mass(f32),
    Detailed(StaticModification),
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

impl StaticModEntry {
    pub fn definition(&self) -> ModificationDefinition {
        match self {
            Self::Mass(mass) => ModificationDefinition::bare(*mass),
            Self::Detailed(modification) => ModificationDefinition::detailed(
                modification.mass,
                &modification.name,
                &modification.neutral_losses,
                modification.neutral_loss_mode,
                &modification.channel_offsets,
            ),
        }
    }
}

/// A variable modification entry may use the existing bare-mass syntax or the
/// extensible object syntax.
#[derive(Clone, Debug, Serialize, schemars::JsonSchema, PartialEq)]
#[serde(untagged)]
pub enum VarModEntry {
    Mass(f32),
    Detailed(VariableModification),
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

impl VarModEntry {
    pub fn mass(&self) -> f32 {
        match self {
            VarModEntry::Mass(m) => *m,
            VarModEntry::Detailed(modification) => modification.mass,
        }
    }

    pub fn max_count(&self) -> Option<usize> {
        match self {
            VarModEntry::Mass(_) => None,
            VarModEntry::Detailed(modification) => modification.max_count,
        }
    }

    pub fn definition(&self) -> ModificationDefinition {
        match self {
            Self::Mass(mass) => ModificationDefinition::bare(*mass),
            Self::Detailed(modification) => ModificationDefinition::detailed(
                modification.mass,
                &modification.name,
                &modification.neutral_losses,
                modification.neutral_loss_mode,
                &modification.channel_offsets,
            ),
        }
    }

    pub fn site_mode(&self) -> SiteMode {
        match self {
            VarModEntry::Mass(_) => SiteMode::Exhaustive,
            VarModEntry::Detailed(modification) => modification.site_mode,
        }
    }

    pub fn search_mode(&self) -> SearchMode {
        match self {
            VarModEntry::Mass(_) => SearchMode::Database,
            VarModEntry::Detailed(modification) => modification.search_mode,
        }
    }

    pub fn channel_offsets(&self) -> &BTreeMap<String, f32> {
        match self {
            VarModEntry::Mass(_) => {
                static EMPTY: std::sync::OnceLock<BTreeMap<String, f32>> =
                    std::sync::OnceLock::new();
                EMPTY.get_or_init(BTreeMap::new)
            }
            VarModEntry::Detailed(modification) => &modification.channel_offsets,
        }
    }
}

impl StaticModEntry {
    pub fn channel_offsets(&self) -> &BTreeMap<String, f32> {
        match self {
            StaticModEntry::Mass(_) => {
                static EMPTY: std::sync::OnceLock<BTreeMap<String, f32>> =
                    std::sync::OnceLock::new();
                EMPTY.get_or_init(BTreeMap::new)
            }
            StaticModEntry::Detailed(modification) => &modification.channel_offsets,
        }
    }
}

use crate::mass::VALID_AA;

#[derive(Copy, Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum ModificationSpecificity {
    PeptideN(Option<u8>),
    PeptideC(Option<u8>),
    ProteinN(Option<u8>),
    ProteinC(Option<u8>),
    Residue(u8),
    /// Residue that is neither the first nor the last residue of the peptide
    Internal(u8),
}

impl ModificationSpecificity {
    /// Is `index` an internal position of a peptide with `len` residues?
    pub fn is_internal(index: usize, len: usize) -> bool {
        index > 0 && index < len.saturating_sub(1)
    }
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
            ModificationSpecificity::Internal(r) => {
                f.write_char('~')?;
                Some(*r)
            }
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
        Ok(match prefix {
            b'^' => Self::PeptideN(residue),
            b'$' => Self::PeptideC(residue),
            b'[' => Self::ProteinN(residue),
            b']' => Self::ProteinC(residue),
            b'~' => Self::Internal(residue.unwrap()),
            _ => unreachable!(),
        })
    }
}

/// Deserialize modification maps without dropping malformed site rules.
pub fn deserialize_mod_map<'de, D, T>(
    deserializer: D,
) -> Result<Option<HashMap<String, T>>, D::Error>
where
    D: Deserializer<'de>,
    T: Deserialize<'de>,
{
    let input = Option::<HashMap<String, T>>::deserialize(deserializer)?;
    if let Some(entries) = &input {
        for key in entries.keys() {
            key.parse::<ModificationSpecificity>().map_err(|_| {
                de::Error::custom(format!("invalid modification key `{key}`. Use a residue, ^, $, [, ], or a residue prefixed by ^, $, [, ], ~"))
            })?;
        }
    }
    Ok(input)
}

pub fn validate_mods(
    input: Option<HashMap<String, StaticModEntry>>,
) -> HashMap<ModificationSpecificity, StaticModEntry> {
    parse_mod_map(input)
}

pub fn validate_var_mods(
    input: Option<HashMap<String, Vec<VarModEntry>>>,
) -> HashMap<ModificationSpecificity, Vec<VarModEntry>> {
    parse_mod_map(input)
}

fn parse_mod_map<T>(input: Option<HashMap<String, T>>) -> HashMap<ModificationSpecificity, T> {
    input
        .unwrap_or_default()
        .into_iter()
        .map(|(key, value)| {
            let specificity = key
                .parse()
                .unwrap_or_else(|_| panic!("invalid modification key `{key}`"));
            (specificity, value)
        })
        .collect()
}

#[cfg(test)]
#[path = "../tests/unit/modification.rs"]
mod test;
