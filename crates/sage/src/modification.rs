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

#[derive(Copy, Clone, Debug, Default, Deserialize, Serialize, PartialEq, Eq, PartialOrd, Ord)]
#[serde(rename_all = "snake_case")]
pub enum NeutralLossMode {
    #[default]
    Optional,
    Required,
}

/// Controls where candidates for a variable modification are generated.
#[derive(Copy, Clone, Debug, Default, Deserialize, Serialize, PartialEq, Eq)]
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
#[derive(Clone, Debug, Serialize, PartialEq)]
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
#[derive(Clone, Debug, Serialize, PartialEq)]
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
    #[serde(default, skip_serializing_if = "BTreeMap::is_empty")]
    pub channel_offsets: BTreeMap<String, f32>,
}

fn is_exhaustive(mode: &SiteMode) -> bool {
    *mode == SiteMode::Exhaustive
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
            max_count: raw.max_count,
            name: raw.name,
            neutral_losses: raw.neutral_losses,
            neutral_loss_mode: raw.neutral_loss_mode,
            site_mode: raw.site_mode,
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

#[derive(Clone, Debug, Serialize, PartialEq)]
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
#[derive(Clone, Debug, Serialize, PartialEq)]
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

pub fn validate_mods(
    input: Option<HashMap<String, StaticModEntry>>,
) -> HashMap<ModificationSpecificity, StaticModEntry> {
    let mut output = HashMap::new();
    if let Some(input) = input {
        for (s, mass) in input {
            match ModificationSpecificity::from_str(&s) {
                Ok(m) => {
                    output.insert(m, mass);
                }
                Err(InvalidModification::Empty) => {
                    log::error!("Invalid modification string: empty")
                }
                Err(InvalidModification::InvalidResidue(c)) => {
                    log::error!("Invalid modification string: unrecognized residue ({})", c)
                }
                Err(InvalidModification::TooLong(s)) => {
                    log::error!("Invalid modification string: {} is too long", s)
                }
            }
        }
    }
    output
}

pub fn validate_var_mods(
    input: Option<HashMap<String, Vec<VarModEntry>>>,
) -> HashMap<ModificationSpecificity, Vec<VarModEntry>> {
    let mut output = HashMap::new();
    if let Some(input) = input {
        for (s, entries) in input {
            match ModificationSpecificity::from_str(&s) {
                Ok(m) => {
                    output.insert(m, entries);
                }
                Err(InvalidModification::Empty) => {
                    log::error!("Skipping invalid modification string: empty")
                }
                Err(InvalidModification::InvalidResidue(c)) => {
                    log::error!(
                        "Skipping invalid modification string: unrecognized residue ({})",
                        c
                    )
                }
                Err(InvalidModification::TooLong(s)) => {
                    log::error!("Skipping invalid modification string: {} is too long", s)
                }
            }
        }
    }
    output
}

#[cfg(test)]
mod test {
    use super::*;

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

    #[test]
    fn var_mod_entry_bare_mass() {
        let entry = VarModEntry::Mass(15.9949);
        assert_eq!(entry.mass(), 15.9949);
        assert_eq!(entry.max_count(), None);
    }

    #[test]
    fn var_mod_entry_detailed_with_limit() {
        let entry = VarModEntry::Detailed(VariableModification {
            mass: 15.9949,
            max_count: Some(1),
            name: None,
            neutral_losses: vec![],
            neutral_loss_mode: NeutralLossMode::Optional,
            site_mode: SiteMode::Exhaustive,
            channel_offsets: Default::default(),
        });
        assert_eq!(entry.mass(), 15.9949);
        assert_eq!(entry.max_count(), Some(1));
    }

    #[test]
    fn deserialize_var_mod_entries() {
        let entries: Vec<VarModEntry> = serde_json::from_str(
            r#"[15.9949, {"mass": 42.0106, "max_count": 1}, {"mass": 14.0157}]"#,
        )
        .unwrap();

        assert_eq!(entries.len(), 3);
        assert!((entries[0].mass() - 15.9949).abs() < 1e-4);
        assert_eq!(entries[0].max_count(), None);
        assert!((entries[1].mass() - 42.0106).abs() < 1e-4);
        assert_eq!(entries[1].max_count(), Some(1));
        assert!((entries[2].mass() - 14.0157).abs() < 1e-4);
        assert_eq!(entries[2].max_count(), None);

        let serialized = serde_json::to_value(&entries).unwrap();
        assert!(serialized[0].is_number());
        assert!(serialized[1].is_object());
        assert_eq!(serialized[1]["max_count"], 1);
        assert!(serialized[2].is_object());
        assert!(serialized[2].get("max_count").is_none());

        let round_tripped: Vec<VarModEntry> = serde_json::from_value(serialized).unwrap();
        assert_eq!(entries, round_tripped);
    }

    #[test]
    fn deserialize_named_neutral_loss_modifications() {
        let entry: VarModEntry = serde_json::from_str(
            r#"{
                "mass": 79.9663,
                "max_count": 2,
                "name": "Phospho",
                "neutral_losses": [97.9769],
                "neutral_loss_mode": "required",
                "site_mode": "both"
            }"#,
        )
        .unwrap();

        let VarModEntry::Detailed(entry) = &entry else {
            panic!("expected structured modification")
        };
        assert_eq!(entry.name.as_deref(), Some("Phospho"));
        assert_eq!(entry.neutral_losses, vec![97.9769]);
        assert_eq!(entry.neutral_loss_mode, NeutralLossMode::Required);
        assert_eq!(entry.site_mode, SiteMode::Both);

        let round_trip: VarModEntry =
            serde_json::from_value(serde_json::to_value(entry).unwrap()).unwrap();
        assert_eq!(round_trip, VarModEntry::Detailed(entry.clone()));
    }

    #[test]
    fn static_mods_accept_numeric_and_structured_entries() {
        let numeric: StaticModEntry = serde_json::from_str("57.0215").unwrap();
        assert!((numeric.definition().mass - 57.0215).abs() < 1e-4);

        let detailed: StaticModEntry = serde_json::from_str(
            r#"{
                "mass": 57.0215,
                "name": "Carbamidomethyl",
                "neutral_losses": [18.0106]
            }"#,
        )
        .unwrap();
        let definition = detailed.definition();
        assert_eq!(definition.name.as_deref(), Some("Carbamidomethyl"));
        assert_eq!(&*definition.neutral_losses, &[18.0106]);
        assert_eq!(definition.neutral_loss_mode, NeutralLossMode::Optional);
    }

    #[test]
    fn channel_offsets_round_trip_on_static_and_variable_modifications() {
        let static_entry: StaticModEntry = serde_json::from_value(serde_json::json!({
            "mass": 0.0,
            "name": "SILAC-K",
            "channel_offsets": {"light": 0.0, "heavy": 8.014199}
        }))
        .unwrap();
        assert_eq!(static_entry.channel_offsets()["heavy"], 8.014199);

        let variable_entry: VarModEntry = serde_json::from_value(serde_json::json!({
            "mass": 0.0,
            "name": "Optional-Lys8",
            "max_count": 2,
            "channel_offsets": {"light": 0.0, "heavy": 8.014199}
        }))
        .unwrap();
        assert_eq!(variable_entry.channel_offsets()["light"], 0.0);
        let serialized = serde_json::to_value(&variable_entry).unwrap();
        let heavy = serialized["channel_offsets"]["heavy"].as_f64().unwrap();
        assert!((heavy - 8.014199).abs() < 1e-6);
    }

    #[test]
    fn channel_offsets_reject_invalid_names_and_masses() {
        for json in [
            r#"{"mass": 0.0, "channel_offsets": {" heavy": 8.0, "light": 0.0}}"#,
            r#"{"mass": 0.0, "channel_offsets": {"heavy": 1e999, "light": 0.0}}"#,
        ] {
            assert!(serde_json::from_str::<VarModEntry>(json).is_err());
        }
    }

    #[test]
    fn reject_invalid_neutral_loss_configuration() {
        for json in [
            r#"{"mass": 79.9663, "name": ""}"#,
            r#"{"mass": 79.9663, "neutral_losses": [0.0]}"#,
            r#"{"mass": 79.9663, "neutral_losses": [-18.0]}"#,
            r#"{"mass": 79.9663, "neutral_loss_mode": "required"}"#,
            r#"{"mass": 79.9663, "neutral_loss_mode": "sometimes"}"#,
        ] {
            assert!(
                serde_json::from_str::<VarModEntry>(json).is_err(),
                "accepted invalid configuration: {json}"
            );
        }
    }

    #[test]
    fn reject_unsupported_var_mod_shapes_and_fields() {
        assert!(serde_json::from_str::<Vec<VarModEntry>>(r#"[[42.0106, 1]]"#).is_err());
        assert!(serde_json::from_str::<Vec<VarModEntry>>(
            r#"[{"mass": 42.0106, "max_counts": 1}]"#
        )
        .is_err());
    }

    #[test]
    fn validate_var_mods_mixed() {
        use ModificationSpecificity::*;
        // Mix bare masses and detailed entries
        let mut raw = HashMap::new();
        raw.insert(
            "M".to_string(),
            vec![
                VarModEntry::Mass(15.9949),
                VarModEntry::Detailed(VariableModification {
                    mass: 15.9949,
                    max_count: Some(1),
                    name: None,
                    neutral_losses: vec![],
                    neutral_loss_mode: NeutralLossMode::Optional,
                    site_mode: SiteMode::Exhaustive,
                    channel_offsets: Default::default(),
                }),
            ],
        );
        raw.insert(
            "C".to_string(),
            vec![VarModEntry::Detailed(VariableModification {
                mass: 57.0215,
                max_count: Some(2),
                name: None,
                neutral_losses: vec![],
                neutral_loss_mode: NeutralLossMode::Optional,
                site_mode: SiteMode::Exhaustive,
                channel_offsets: Default::default(),
            })],
        );
        let result = validate_var_mods(Some(raw));

        let m_entries = result.get(&Residue(b'M')).unwrap();
        assert_eq!(m_entries.len(), 2);
        assert!((m_entries[0].mass() - 15.9949).abs() < 1e-4);
        assert_eq!(m_entries[0].max_count(), None);
        assert!((m_entries[1].mass() - 15.9949).abs() < 1e-4);
        assert_eq!(m_entries[1].max_count(), Some(1));

        let c_entries = result.get(&Residue(b'C')).unwrap();
        assert_eq!(c_entries.len(), 1);
        assert!((c_entries[0].mass() - 57.0215).abs() < 1e-4);
        assert_eq!(c_entries[0].max_count(), Some(2));
    }

    #[test]
    fn validate_var_mods_invalid_residue_skipped() {
        let mut raw = HashMap::new();
        raw.insert("Z".to_string(), vec![VarModEntry::Mass(15.9949)]);
        raw.insert("M".to_string(), vec![VarModEntry::Mass(15.9949)]);
        let result = validate_var_mods(Some(raw));
        // Z is invalid — only M should survive
        assert_eq!(result.len(), 1);
        assert!(result.contains_key(&ModificationSpecificity::Residue(b'M')));
    }
}
