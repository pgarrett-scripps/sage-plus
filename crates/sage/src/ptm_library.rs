use std::collections::{HashMap, HashSet};
use std::sync::Arc;

/// One observed PTM location loaded from the site-library Parquet file.
#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub struct PtmLibrarySite {
    pub protein: Arc<str>,
    /// Zero-based protein position. The Parquet representation is one-based.
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
}
