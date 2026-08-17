use std::collections::{HashMap, HashSet};
use std::sync::Arc;

/// One observed PTM location loaded from a site-library file.
#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub struct PtmLibrarySite {
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
mod tests {
    use super::*;

    #[test]
    fn parses_tsv_with_extra_columns_and_deduplicates() {
        let contents = concat!(
            "score\tmodification\tresidue\tposition\tprotein\n",
            "0.99\tPhospho\tS\t3\tP12345\n",
            "0.95\tPhospho\tS\t3\tP12345\n",
            "0.90\tOxidation\tm\t7\tP12345\n",
        );
        let library = PtmLibrary::from_tsv(contents).unwrap();
        assert_eq!(library.len(), 2);
        let sites = library.sites_for("P12345");
        assert_eq!(sites[0].position, 2);
        assert_eq!(sites[0].residue, b'S');
        assert_eq!(sites[1].position, 6);
        assert_eq!(sites[1].residue, b'M');
    }

    #[test]
    fn rejects_zero_based_tsv_position() {
        let error = PtmLibrary::from_tsv(
            "protein\tposition\tresidue\tmodification\nP12345\t0\tS\tPhospho\n",
        )
        .unwrap_err();
        assert!(error.contains("positions are one-based"));
    }

    #[test]
    fn detects_plain_and_compressed_tsv_paths() {
        assert!(is_tsv_path("sites.TSV"));
        assert!(is_tsv_path("s3://bucket/sites.tsv.gz"));
        assert!(!is_tsv_path("sites.parquet"));
    }
}
