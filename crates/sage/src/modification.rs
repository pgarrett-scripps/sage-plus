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
    PeptideNTerm(u8),
    PeptideCTerm(u8),
    ProteinNTerm(u8),
    ProteinCTerm(u8),
    /// Residue selected by a sequence motif evaluated with protein context.
    Motif(&'static crate::motif::SiteMotif),
}

impl ModificationSpecificity {
    /// Explicit public spelling of one complete attachment rule.
    pub fn explicit_name(self) -> String {
        match self {
            Self::Residue(r) => (r as char).to_string(),
            Self::Internal(r) => format!("internal_residue:{}", r as char),
            Self::PeptideN(Some(r)) => format!("first_residue:{}", r as char),
            Self::PeptideC(Some(r)) => format!("last_residue:{}", r as char),
            Self::ProteinN(Some(r)) => format!("protein_first:{}", r as char),
            Self::ProteinC(Some(r)) => format!("protein_last:{}", r as char),
            Self::PeptideN(None) => "peptide_n_term".into(),
            Self::PeptideC(None) => "peptide_c_term".into(),
            Self::ProteinN(None) => "protein_n_term".into(),
            Self::ProteinC(None) => "protein_c_term".into(),
            Self::PeptideNTerm(r) => format!("peptide_n_term:{}", r as char),
            Self::PeptideCTerm(r) => format!("peptide_c_term:{}", r as char),
            Self::ProteinNTerm(r) => format!("protein_n_term:{}", r as char),
            Self::ProteinCTerm(r) => format!("protein_c_term:{}", r as char),
            Self::Motif(motif) => motif.canonical().to_string(),
        }
    }

    pub fn is_motif(self) -> bool {
        matches!(self, Self::Motif(_))
    }

    /// Select sites with known protein flanking residues. Only motifs use the
    /// flanks; every other rule is decided by the peptide and its position.
    pub fn sites_with_flanks(
        self,
        sequence: &[u8],
        position: crate::enzyme::Position,
        left: &[u8],
        right: &[u8],
    ) -> Vec<crate::peptide::Site> {
        use crate::enzyme::Position;
        match self {
            Self::Motif(motif) => motif
                .sites(
                    sequence,
                    crate::motif::MotifContext {
                        left,
                        right,
                        left_boundary: left.is_empty()
                            && matches!(position, Position::Nterm | Position::Full),
                        right_boundary: right.is_empty()
                            && matches!(position, Position::Cterm | Position::Full),
                    },
                )
                .into_iter()
                .map(crate::peptide::Site::Sequence)
                .collect(),
            _ => self.sites(sequence, position),
        }
    }

    /// Select physical attachment sites using one matcher for every consumer.
    pub fn sites(
        self,
        sequence: &[u8],
        position: crate::enzyme::Position,
    ) -> Vec<crate::peptide::Site> {
        use crate::enzyme::Position;
        use crate::peptide::Site;
        let Some(last) = sequence.len().checked_sub(1) else {
            return Vec::new();
        };
        let protein_n = matches!(position, Position::Nterm | Position::Full);
        let protein_c = matches!(position, Position::Cterm | Position::Full);
        match self {
            Self::Residue(r) | Self::Internal(r) => sequence
                .iter()
                .enumerate()
                .filter(|(i, observed)| {
                    **observed == r
                        && (!matches!(self, Self::Internal(_))
                            || Self::is_internal(*i, sequence.len()))
                })
                .map(|(i, _)| Site::Sequence(i as u32))
                .collect(),
            Self::PeptideN(None) => vec![Site::Nterm],
            Self::PeptideC(None) => vec![Site::Cterm],
            Self::ProteinN(None) if protein_n => vec![Site::Nterm],
            Self::ProteinC(None) if protein_c => vec![Site::Cterm],
            Self::PeptideN(Some(r)) if sequence[0] == r => vec![Site::Sequence(0)],
            Self::PeptideC(Some(r)) if sequence[last] == r => vec![Site::Sequence(last as u32)],
            Self::ProteinN(Some(r)) if protein_n && sequence[0] == r => vec![Site::Sequence(0)],
            Self::ProteinC(Some(r)) if protein_c && sequence[last] == r => {
                vec![Site::Sequence(last as u32)]
            }
            Self::PeptideNTerm(r) if sequence[0] == r => vec![Site::Nterm],
            Self::PeptideCTerm(r) if sequence[last] == r => vec![Site::Cterm],
            Self::ProteinNTerm(r) if protein_n && sequence[0] == r => vec![Site::Nterm],
            Self::ProteinCTerm(r) if protein_c && sequence[last] == r => vec![Site::Cterm],
            Self::Motif(_) => self.sites_with_flanks(sequence, position, &[], &[]),
            _ => Vec::new(),
        }
    }

    pub fn overlaps(self, other: Self) -> bool {
        use crate::enzyme::Position;
        if let Some(motif) = [self, other].into_iter().find_map(|rule| match rule {
            Self::Motif(motif) => Some(motif),
            _ => None,
        }) {
            // Conservative: motif rules overlap any residue rule sharing a residue.
            let residues = |rule| match rule {
                Self::Motif(m) => m.site_residues(),
                Self::Residue(r)
                | Self::Internal(r)
                | Self::PeptideN(Some(r))
                | Self::PeptideC(Some(r))
                | Self::ProteinN(Some(r))
                | Self::ProteinC(Some(r)) => vec![r],
                _ => Vec::new(),
            };
            let other = if self.is_motif() { other } else { self };
            let mine = motif.site_residues();
            return residues(other).iter().any(|r| mine.contains(r));
        }
        let residue = |rule| match rule {
            Self::Residue(r)
            | Self::Internal(r)
            | Self::PeptideNTerm(r)
            | Self::PeptideCTerm(r)
            | Self::ProteinNTerm(r)
            | Self::ProteinCTerm(r) => Some(r),
            Self::PeptideN(r) | Self::PeptideC(r) | Self::ProteinN(r) | Self::ProteinC(r) => r,
            Self::Motif(_) => unreachable!(),
        };
        let mut alphabet = vec![b'A'];
        alphabet.extend(residue(self));
        alphabet.extend(residue(other));
        alphabet.sort_unstable();
        alphabet.dedup();
        for &first in &alphabet {
            for &middle in &alphabet {
                for &last in &alphabet {
                    let sequence = [first, middle, last];
                    for length in 1..=3 {
                        let left = self.sites(&sequence[..length], Position::Full);
                        if other
                            .sites(&sequence[..length], Position::Full)
                            .iter()
                            .any(|site| left.contains(site))
                        {
                            return true;
                        }
                    }
                }
            }
        }
        false
    }

    /// Is `index` an internal position of a peptide with `len` residues?
    pub fn is_internal(index: usize, len: usize) -> bool {
        index > 0 && index < len.saturating_sub(1)
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

/// Public definition keyed by modification identity, with explicit attachment rules.
#[derive(schemars::JsonSchema)]
pub struct NamedStaticModification {
    pub sites: Vec<String>,
    #[serde(flatten)]
    pub details: StaticModification,
}

#[derive(schemars::JsonSchema)]
pub struct NamedVariableModification {
    pub sites: Vec<String>,
    #[serde(flatten)]
    pub details: VariableModification,
}

#[derive(schemars::JsonSchema)]
#[serde(untagged)]
pub enum StaticModConfig {
    Named(BTreeMap<String, NamedStaticModification>),
    Legacy(HashMap<String, StaticModEntry>),
}

#[derive(schemars::JsonSchema)]
#[serde(untagged)]
pub enum VariableModConfig {
    Named(BTreeMap<String, NamedVariableModification>),
    Legacy(HashMap<String, Vec<VarModEntry>>),
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

/// Normalize named definitions and legacy configurations into one engine model.
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
            let specificity = id.parse::<ModificationSpecificity>().map_err(|_| de::Error::custom(format!("invalid modification key `{id}`. Named definitions require a nonempty `sites` array")))?;
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
        if let Some(name) = object.get("name") {
            if name.as_str() != Some(id.as_str()) {
                return Err(de::Error::custom(format!(
                    "modification `{id}` must not override its identity with a different name"
                )));
            }
        }
        object.insert("name".into(), serde_json::Value::String(id.clone()));
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

fn named_modifications<'a>(
    entries: impl Iterator<Item = (ModificationSpecificity, &'a dyn NamedEntry)>,
) -> serde_json::Value {
    let mut named = BTreeMap::<String, serde_json::Value>::new();
    for (site, entry) in entries {
        let mut value = entry.json();
        let name = value
            .get("name")
            .and_then(|v| v.as_str())
            .expect("named serialization requires names")
            .to_string();
        value.as_object_mut().unwrap().remove("name");
        let object = named.entry(name).or_insert_with(|| {
            value["sites"] = serde_json::json!([]);
            value
        });
        object["sites"]
            .as_array_mut()
            .unwrap()
            .push(serde_json::json!(site.explicit_name()));
    }
    for value in named.values_mut() {
        let sites = value["sites"].as_array_mut().unwrap();
        sites.sort_by(|a, b| a.as_str().cmp(&b.as_str()));
        sites.dedup();
    }
    serde_json::to_value(named).unwrap()
}

trait NamedEntry {
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
