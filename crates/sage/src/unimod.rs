use std::collections::HashMap;
use std::sync::{Mutex, OnceLock};

use serde::Deserialize;

const UNIMOD_JSON: &str = include_str!("../data/unimodifications.json");

#[derive(Deserialize)]
struct Entry {
    id: String,
    name: String,
    monoisotopic_mass: f64,
}

#[derive(Deserialize)]
struct File {
    unimodifications: Vec<Entry>,
}

struct Tables {
    by_accession: HashMap<u32, (f32, String)>,
    /// Lower-cased name -> mass.
    by_name: HashMap<String, f32>,
    /// Lower-cased name -> canonical capitalization.
    canonical_name: HashMap<String, String>,
}

static TABLES: OnceLock<Tables> = OnceLock::new();

fn tables() -> &'static Tables {
    TABLES.get_or_init(|| {
        let parsed: File = serde_json::from_str(UNIMOD_JSON)
            .expect("embedded unimodifications.json is malformed");
        let mut by_accession = HashMap::new();
        let mut by_name = HashMap::new();
        let mut canonical_name = HashMap::new();
        for e in parsed.unimodifications {
            let mass = e.monoisotopic_mass as f32;
            if let Ok(id) = e.id.parse::<u32>() {
                by_accession.insert(id, (mass, e.name.clone()));
            }
            let lower = e.name.to_ascii_lowercase();
            by_name.entry(lower.clone()).or_insert(mass);
            canonical_name.entry(lower).or_insert(e.name);
        }
        Tables {
            by_accession,
            by_name,
            canonical_name,
        }
    })
}

/// Look up the monoisotopic delta mass (Da) for a Unimod accession.
pub fn delta_mass(accession: u32) -> Option<f32> {
    tables().by_accession.get(&accession).map(|(m, _)| *m)
}

/// Look up the monoisotopic delta mass for a Unimod modification name
/// (case-insensitive). Returns None for unknown names.
pub fn mass_by_name(name: &str) -> Option<f32> {
    tables().by_name.get(&name.to_ascii_lowercase()).copied()
}

/// Return the canonical (case-preserving) Unimod name for the supplied
/// case-insensitive `name`, if known.
pub fn canonical_name(name: &str) -> Option<&'static str> {
    tables()
        .canonical_name
        .get(&name.to_ascii_lowercase())
        .map(|s| s.as_str())
}

fn key(mass: f32) -> i64 {
    (mass * 1e5).round() as i64
}

static LABELS: OnceLock<Mutex<HashMap<i64, String>>> = OnceLock::new();

fn labels() -> &'static Mutex<HashMap<i64, String>> {
    LABELS.get_or_init(|| Mutex::new(HashMap::new()))
}

/// Associate a delta mass with a human-readable Unimod name, so that any
/// peptide carrying that exact mass (within 5 decimal places) renders as
/// `[Name]` instead of `[+mass]` in output. First registration wins for a
/// given mass; later calls with a different name are ignored.
pub fn register_label(mass: f32, name: &str) {
    let mut guard = labels().lock().expect("unimod label registry poisoned");
    guard.entry(key(mass)).or_insert_with(|| name.to_string());
}

/// Look up the registered display label for a delta mass, if any.
pub fn label_for(mass: f32) -> Option<String> {
    let guard = labels().lock().expect("unimod label registry poisoned");
    guard.get(&key(mass)).cloned()
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

    #[test]
    fn names_resolve_case_insensitively() {
        let m = mass_by_name("Oxidation").expect("Oxidation missing");
        assert!((m - 15.994915).abs() < 1e-4);
        let m2 = mass_by_name("oxidation").expect("lowercase Oxidation missing");
        assert_eq!(m, m2);
        assert!(mass_by_name("totally-not-a-mod").is_none());
        assert_eq!(canonical_name("phospho"), Some("Phospho"));
    }

    #[test]
    fn labels_first_write_wins() {
        // Use a unique mass so this test doesn't collide with other suites.
        let m = 1234.56789_f32;
        register_label(m, "MyMod");
        register_label(m, "OtherMod");
        assert_eq!(label_for(m).as_deref(), Some("MyMod"));
        assert!(label_for(0.000_001).is_none());
    }
}
