use std::collections::HashMap;
use std::sync::OnceLock;

use serde::Deserialize;

const UNIMOD_JSON: &str = include_str!("../data/unimodifications.json");

#[derive(Deserialize)]
struct Entry {
    id: String,
    monoisotopic_mass: f64,
}

#[derive(Deserialize)]
struct File {
    unimodifications: Vec<Entry>,
}

static TABLE: OnceLock<HashMap<u32, f32>> = OnceLock::new();

fn table() -> &'static HashMap<u32, f32> {
    TABLE.get_or_init(|| {
        let parsed: File = serde_json::from_str(UNIMOD_JSON)
            .expect("embedded unimodifications.json is malformed");
        parsed
            .unimodifications
            .into_iter()
            .filter_map(|e| e.id.parse::<u32>().ok().map(|id| (id, e.monoisotopic_mass as f32)))
            .collect()
    })
}

/// Look up the monoisotopic delta mass (Da) for a Unimod accession.
pub fn delta_mass(accession: u32) -> Option<f32> {
    table().get(&accession).copied()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn known_accessions() {
        let oxidation = delta_mass(35).expect("UNIMOD:35 missing");
        assert!((oxidation - 15.994915).abs() < 1e-4, "got {oxidation}");
        let carbamidomethyl = delta_mass(4).expect("UNIMOD:4 missing");
        assert!((carbamidomethyl - 57.021464).abs() < 1e-4);
        let phospho = delta_mass(21).expect("UNIMOD:21 missing");
        assert!((phospho - 79.966331).abs() < 1e-4);
    }

    #[test]
    fn unknown_accession() {
        assert!(delta_mass(9_999_999).is_none());
    }
}
