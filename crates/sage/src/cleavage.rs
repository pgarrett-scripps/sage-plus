use crate::fasta::Fasta;
use crate::mass::VALID_AA;
use std::collections::{BTreeMap, HashMap};
use std::error::Error;
use std::fmt::{Display, Formatter};

#[derive(Clone, Debug, PartialEq, Eq)]
struct CustomCleavageSite {
    position: usize,
    context: Option<String>,
}

/// Protein-specific cleavage sites parsed from a TSV file.
///
/// `position` is the zero-based index of the residue immediately before the
/// cleavage boundary. An optional `context` uses `|` to mark that boundary,
/// for example `KLGF|APQT`.
#[derive(Clone, Debug, Default)]
pub struct CustomCleavageLibrary {
    sites: HashMap<String, Vec<CustomCleavageSite>>,
}

/// A cleavage library whose accessions, coordinates, and contexts have been
/// checked against the FASTA used for the search.
#[derive(Clone, Debug, Default)]
pub struct ValidatedCustomCleavageLibrary {
    boundaries: HashMap<String, Vec<usize>>,
    pub total_sites: usize,
    pub matched_sites: usize,
    pub unmatched_sites: usize,
    pub sites_without_context: usize,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct CustomCleavageError(String);

impl CustomCleavageError {
    fn new(message: impl Into<String>) -> Self {
        Self(message.into())
    }
}

impl Display for CustomCleavageError {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        f.write_str(&self.0)
    }
}

impl Error for CustomCleavageError {}

impl CustomCleavageLibrary {
    pub fn from_tsv(content: &str) -> Result<Self, CustomCleavageError> {
        let mut lines = content.lines().enumerate().filter_map(|(index, line)| {
            let line = line.trim_end_matches('\r');
            (!line.trim().is_empty()).then_some((index + 1, line))
        });
        let (header_line, header) = lines
            .next()
            .ok_or_else(|| CustomCleavageError::new("custom cleavage-site TSV is empty"))?;
        let columns = header.split('\t').map(str::trim).collect::<Vec<_>>();
        let protein_col = columns
            .iter()
            .position(|column| *column == "protein")
            .ok_or_else(|| {
                CustomCleavageError::new(format!(
                    "custom cleavage-site TSV header on line {header_line} is missing required `protein` column"
                ))
            })?;
        let position_col = columns
            .iter()
            .position(|column| *column == "position")
            .ok_or_else(|| {
                CustomCleavageError::new(format!(
                    "custom cleavage-site TSV header on line {header_line} is missing required `position` column"
                ))
            })?;
        let context_col = columns.iter().position(|column| *column == "context");

        let mut sites: HashMap<String, BTreeMap<usize, Option<String>>> = HashMap::new();
        for (line_number, line) in lines {
            let fields = line.split('\t').map(str::trim).collect::<Vec<_>>();
            let protein = fields.get(protein_col).copied().unwrap_or_default();
            if protein.is_empty() {
                return Err(CustomCleavageError::new(format!(
                    "custom cleavage-site TSV line {line_number} has an empty `protein`"
                )));
            }
            let raw_position = fields.get(position_col).copied().unwrap_or_default();
            let position = raw_position.parse::<usize>().map_err(|_| {
                CustomCleavageError::new(format!(
                    "custom cleavage-site TSV line {line_number} has invalid `position` `{raw_position}`; expected a zero-based residue index"
                ))
            })?;
            let context = context_col
                .and_then(|column| fields.get(column))
                .copied()
                .filter(|context| !context.is_empty())
                .map(|context| {
                    validate_context(
                        context,
                        &format!("custom cleavage-site TSV line {line_number}"),
                    )
                })
                .transpose()?;
            insert_site(&mut sites, protein.to_string(), position, context)?;
        }

        Self::from_site_map(sites, "custom cleavage-site TSV contains no data rows")
    }

    /// Build a library from typed records, as used by columnar input formats.
    pub fn from_records<I>(records: I) -> Result<Self, CustomCleavageError>
    where
        I: IntoIterator<Item = (String, usize, Option<String>)>,
    {
        let mut sites: HashMap<String, BTreeMap<usize, Option<String>>> = HashMap::new();
        for (index, (protein, position, context)) in records.into_iter().enumerate() {
            let row = index + 1;
            if protein.trim().is_empty() {
                return Err(CustomCleavageError::new(format!(
                    "custom cleavage-site record {row} has an empty `protein`"
                )));
            }
            let context = context
                .filter(|context| !context.trim().is_empty())
                .map(|context| {
                    validate_context(
                        context.trim(),
                        &format!("custom cleavage-site record {row}"),
                    )
                })
                .transpose()?;
            insert_site(&mut sites, protein.trim().to_string(), position, context)?;
        }

        Self::from_site_map(sites, "custom cleavage-site input contains no records")
    }

    fn from_site_map(
        sites: HashMap<String, BTreeMap<usize, Option<String>>>,
        empty_error: &str,
    ) -> Result<Self, CustomCleavageError> {
        if sites.is_empty() {
            return Err(CustomCleavageError::new(empty_error));
        }
        Ok(Self {
            sites: sites
                .into_iter()
                .map(|(protein, sites)| {
                    (
                        protein,
                        sites
                            .into_iter()
                            .map(|(position, context)| CustomCleavageSite { position, context })
                            .collect(),
                    )
                })
                .collect(),
        })
    }

    pub fn validate(
        &self,
        fasta: &Fasta,
    ) -> Result<ValidatedCustomCleavageLibrary, CustomCleavageError> {
        let mut sequences: HashMap<&str, Vec<&str>> = HashMap::new();
        for (protein, sequence) in &fasta.targets {
            sequences.entry(protein).or_default().push(sequence);
        }

        let total_sites = self.sites.values().map(Vec::len).sum();
        let mut matched_sites = 0;
        let mut unmatched_sites = 0;
        let mut sites_without_context = 0;
        let mut boundaries = HashMap::new();

        for (protein, sites) in &self.sites {
            let Some(protein_sequences) = sequences.get(protein.as_str()) else {
                unmatched_sites += sites.len();
                continue;
            };
            let mut protein_boundaries = Vec::with_capacity(sites.len());
            for site in sites {
                let boundary = site.position.checked_add(1).ok_or_else(|| {
                    CustomCleavageError::new(format!(
                        "custom cleavage position overflows for protein `{protein}`"
                    ))
                })?;
                for sequence in protein_sequences {
                    if boundary >= sequence.len() {
                        return Err(CustomCleavageError::new(format!(
                            "custom cleavage `{protein}` position {} is not internal to the {}-residue FASTA sequence",
                            site.position,
                            sequence.len()
                        )));
                    }
                    if let Some(context) = &site.context {
                        validate_sequence_context(
                            protein,
                            site.position,
                            sequence,
                            boundary,
                            context,
                        )?;
                    }
                }
                if site.context.is_none() {
                    sites_without_context += 1;
                }
                matched_sites += 1;
                protein_boundaries.push(boundary);
            }
            boundaries.insert(protein.clone(), protein_boundaries);
        }

        if matched_sites == 0 {
            return Err(CustomCleavageError::new(
                "none of the custom cleavage-site proteins matched a FASTA accession",
            ));
        }

        Ok(ValidatedCustomCleavageLibrary {
            boundaries,
            total_sites,
            matched_sites,
            unmatched_sites,
            sites_without_context,
        })
    }
}

impl ValidatedCustomCleavageLibrary {
    pub fn boundaries_for(&self, protein: &str) -> &[usize] {
        self.boundaries
            .get(protein)
            .map(Vec::as_slice)
            .unwrap_or_default()
    }
}

fn insert_site(
    sites: &mut HashMap<String, BTreeMap<usize, Option<String>>>,
    protein: String,
    position: usize,
    context: Option<String>,
) -> Result<(), CustomCleavageError> {
    let protein_sites = sites.entry(protein.clone()).or_default();
    if let Some(existing) = protein_sites.get(&position) {
        match (existing, &context) {
            (Some(existing), Some(context)) if existing != context => {
                return Err(CustomCleavageError::new(format!(
                    "custom cleavage-site input has conflicting contexts for `{protein}` position {position}"
                )));
            }
            (None, Some(_)) => {
                protein_sites.insert(position, context);
            }
            _ => {}
        }
    } else {
        protein_sites.insert(position, context);
    }
    Ok(())
}

fn validate_context(context: &str, location: &str) -> Result<String, CustomCleavageError> {
    let mut pieces = context.split('|');
    let left = pieces.next().unwrap_or_default();
    let right = pieces.next().unwrap_or_default();
    if left.is_empty() || right.is_empty() || pieces.next().is_some() {
        return Err(CustomCleavageError::new(format!(
            "{location} has invalid context `{context}`; expected amino acids on both sides of one `|`"
        )));
    }
    if left
        .bytes()
        .chain(right.bytes())
        .any(|residue| !VALID_AA.contains(&residue))
    {
        return Err(CustomCleavageError::new(format!(
            "{location} has invalid amino acids in context `{context}`"
        )));
    }
    Ok(context.to_string())
}

fn validate_sequence_context(
    protein: &str,
    position: usize,
    sequence: &str,
    boundary: usize,
    context: &str,
) -> Result<(), CustomCleavageError> {
    let (left, right) = context
        .split_once('|')
        .expect("contexts are validated while parsing");
    let matches = sequence[..boundary].ends_with(left) && sequence[boundary..].starts_with(right);
    if !matches {
        return Err(CustomCleavageError::new(format!(
            "custom cleavage context `{context}` does not match FASTA protein `{protein}` at position {position}"
        )));
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parses_deduplicates_and_validates_sites() {
        let library = CustomCleavageLibrary::from_tsv(
            "protein\tposition\tcontext\nP1\t4\tPEPK|TIDE\nP1\t4\tPEPK|TIDE\nP2\t0\t\n",
        )
        .unwrap();
        let fasta = Fasta::parse(
            ">P1 description\nMPEPKTIDER\n>P2\nACDE\n".into(),
            "rev_",
            true,
        );
        let validated = library.validate(&fasta).unwrap();

        assert_eq!(validated.total_sites, 2);
        assert_eq!(validated.matched_sites, 2);
        assert_eq!(validated.unmatched_sites, 0);
        assert_eq!(validated.sites_without_context, 1);
        assert_eq!(validated.boundaries_for("P1"), &[5]);
        assert_eq!(validated.boundaries_for("P2"), &[1]);
    }

    #[test]
    fn reports_context_mismatch() {
        let library =
            CustomCleavageLibrary::from_tsv("protein\tposition\tcontext\nP1\t4\tPEPR|TIDE\n")
                .unwrap();
        let fasta = Fasta::parse(">P1\nMPEPKTIDER\n".into(), "rev_", true);

        let error = library.validate(&fasta).unwrap_err().to_string();
        assert!(error.contains("does not match"));
    }

    #[test]
    fn rejects_terminal_and_negative_positions() {
        let negative = CustomCleavageLibrary::from_tsv("protein\tposition\nP1\t-1\n")
            .unwrap_err()
            .to_string();
        assert!(negative.contains("zero-based residue index"));

        let library = CustomCleavageLibrary::from_tsv("protein\tposition\nP1\t3\n").unwrap();
        let fasta = Fasta::parse(">P1\nACDE\n".into(), "rev_", true);
        let terminal = library.validate(&fasta).unwrap_err().to_string();
        assert!(terminal.contains("not internal"));
    }

    #[test]
    fn allows_unmatched_library_subset_but_not_zero_matches() {
        let library =
            CustomCleavageLibrary::from_tsv("protein\tposition\nP1\t0\nMISSING\t1\n").unwrap();
        let fasta = Fasta::parse(">P1\nACDE\n".into(), "rev_", true);
        let validated = library.validate(&fasta).unwrap();
        assert_eq!(validated.matched_sites, 1);
        assert_eq!(validated.unmatched_sites, 1);

        let fasta = Fasta::parse(">OTHER\nACDE\n".into(), "rev_", true);
        assert!(library.validate(&fasta).is_err());
    }
}
