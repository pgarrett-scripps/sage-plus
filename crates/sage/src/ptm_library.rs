use std::collections::{HashMap, HashSet};
use std::sync::Arc;

/// Attachment identity is independent of the adjacent residue coordinate.
#[derive(
    Clone,
    Copy,
    Debug,
    Default,
    PartialEq,
    Eq,
    Hash,
    PartialOrd,
    Ord,
    serde::Serialize,
    serde::Deserialize,
)]
#[serde(rename_all = "snake_case")]
pub enum Attachment {
    #[default]
    Residue,
    PeptideNTerm,
    PeptideCTerm,
    ProteinNTerm,
    ProteinCTerm,
}

impl Attachment {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::Residue => "residue",
            Self::PeptideNTerm => "peptide_n_term",
            Self::PeptideCTerm => "peptide_c_term",
            Self::ProteinNTerm => "protein_n_term",
            Self::ProteinCTerm => "protein_c_term",
        }
    }
    pub fn parse(value: &str) -> Result<Self, String> {
        match value {
            "residue" => Ok(Self::Residue),
            "peptide_n_term" => Ok(Self::PeptideNTerm),
            "peptide_c_term" => Ok(Self::PeptideCTerm),
            "protein_n_term" => Ok(Self::ProteinNTerm),
            "protein_c_term" => Ok(Self::ProteinCTerm),
            _ => Err(format!("invalid PTM attachment `{value}`")),
        }
    }
    pub fn site(
        self,
        index: u32,
        length: usize,
        position: crate::enzyme::Position,
    ) -> Option<crate::peptide::Site> {
        use crate::enzyme::Position;
        use crate::peptide::Site;
        if index as usize >= length {
            return None;
        }
        match self {
            Self::Residue => Some(Site::Sequence(index)),
            Self::PeptideNTerm if index == 0 => Some(Site::Nterm),
            Self::PeptideCTerm if index as usize + 1 == length => Some(Site::Cterm),
            Self::ProteinNTerm
                if index == 0 && matches!(position, Position::Nterm | Position::Full) =>
            {
                Some(Site::Nterm)
            }
            Self::ProteinCTerm
                if index as usize + 1 == length
                    && matches!(position, Position::Cterm | Position::Full) =>
            {
                Some(Site::Cterm)
            }
            _ => None,
        }
    }
    pub fn from_site(site: crate::peptide::Site) -> Self {
        match site {
            crate::peptide::Site::Nterm => Self::PeptideNTerm,
            crate::peptide::Site::Cterm => Self::PeptideCTerm,
            crate::peptide::Site::Sequence(_) => Self::Residue,
        }
    }
}

/// One observed PTM location loaded from a site-library file.
#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub struct PtmLibrarySite {
    pub attachment: Attachment,
    pub protein: Arc<str>,
    /// Zero-based protein position. File representations are one-based.
    pub position: u32,
    pub residue: u8,
    /// Name of a variable modification defined in the search configuration.
    pub modification: Arc<str>,
}

/// Sites indexed by protein accession for database expansion.
#[derive(Clone, Debug, Default)]
pub struct PtmLibrary {
    sites: HashMap<Arc<str>, Vec<PtmLibrarySite>>,
    len: usize,
}

impl PtmLibrary {
    pub fn new(input: Vec<PtmLibrarySite>) -> Self {
        let mut seen = HashSet::new();
        let mut sites: HashMap<Arc<str>, Vec<PtmLibrarySite>> = HashMap::new();
        for site in input {
            if seen.insert(site.clone()) {
                sites.entry(site.protein.clone()).or_default().push(site);
            }
        }
        for protein_sites in sites.values_mut() {
            protein_sites.sort_unstable_by(|a, b| {
                a.position
                    .cmp(&b.position)
                    .then_with(|| a.modification.cmp(&b.modification))
                    .then_with(|| a.attachment.cmp(&b.attachment))
            });
        }
        let len = seen.len();
        Self { sites, len }
    }

    pub fn len(&self) -> usize {
        self.len
    }

    pub fn is_empty(&self) -> bool {
        self.len == 0
    }

    pub fn sites_for(&self, protein: &str) -> &[PtmLibrarySite] {
        self.sites
            .get(protein)
            .map(Vec::as_slice)
            .unwrap_or_default()
    }

    pub fn iter(&self) -> impl Iterator<Item = &PtmLibrarySite> {
        self.sites.values().flatten()
    }

    /// Parse a tab-separated PTM library. Required columns may appear in any
    /// order and additional evidence columns are ignored.
    pub fn from_tsv(contents: &str) -> Result<Self, String> {
        let mut reader = csv::ReaderBuilder::new()
            .delimiter(b'\t')
            .trim(csv::Trim::All)
            .from_reader(contents.as_bytes());
        let headers = reader.headers().map_err(|error| error.to_string())?;
        let column = |name: &str| {
            headers
                .iter()
                .position(|header| header.trim_start_matches('\u{feff}') == name)
                .ok_or_else(|| format!("PTM library is missing required column `{name}`"))
        };
        let protein = column("protein")?;
        let position = column("position")?;
        let residue = column("residue")?;
        let modification = column("modification")?;
        let attachment = headers.iter().position(|header| header == "attachment");

        let mut sites = Vec::new();
        for (index, record) in reader.records().enumerate() {
            let row = index + 2;
            let record =
                record.map_err(|error| format!("invalid PTM library row {row}: {error}"))?;
            let field = |column: usize, name: &str| {
                record
                    .get(column)
                    .filter(|value| !value.is_empty())
                    .ok_or_else(|| format!("PTM library row {row} has an empty `{name}`"))
            };
            let protein = field(protein, "protein")?;
            let position = field(position, "position")?
                .parse::<u32>()
                .map_err(|_| format!("PTM library row {row} has an invalid `position`"))?
                .checked_sub(1)
                .ok_or_else(|| {
                    format!("PTM library row {row} has position 0; positions are one-based")
                })?;
            let residue = field(residue, "residue")?.as_bytes();
            if residue.len() != 1 || !residue[0].is_ascii_alphabetic() {
                return Err(format!(
                    "PTM library row {row} has an invalid one-letter `residue`"
                ));
            }
            let modification = field(modification, "modification")?;
            sites.push(PtmLibrarySite {
                attachment: attachment
                    .map(|column| field(column, "attachment").and_then(Attachment::parse))
                    .transpose()?
                    .unwrap_or_default(),
                protein: Arc::from(protein),
                position,
                residue: residue[0].to_ascii_uppercase(),
                modification: Arc::from(modification),
            });
        }
        Ok(Self::new(sites))
    }
}

/// TSV inputs are selected by extension, including transparently compressed
/// `.tsv.gz` files.
pub fn is_tsv_path(path: &str) -> bool {
    let path = path.to_ascii_lowercase();
    path.strip_suffix(".gz").unwrap_or(&path).ends_with(".tsv")
}

#[cfg(test)]
#[path = "../tests/unit/ptm_library.rs"]
mod tests;
