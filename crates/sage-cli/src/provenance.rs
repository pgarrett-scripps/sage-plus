//! Run provenance: the software version, effective configuration, inputs and
//! database that produced a set of outputs. The same record is written to the
//! key-value footer of every Parquet output (see [`Provenance::parquet_metadata`])
//! and to `run-summary.json` under `provenance.metadata`.

use crate::input::Search;
use rayon::prelude::*;
use sage_core::database::{IndexedDatabase, Parameters};
use serde::{Deserialize, Serialize};
use std::collections::HashMap;
use url::Url;

/// Version of the provenance record layout, written as `sage.provenance.version`.
/// Fields may be added within a version; removing or changing one bumps it.
pub const PROVENANCE_VERSION: u32 = 1;

/// Organisms listed per FASTA, most frequent first; the rest are only counted.
const MAX_ORGANISMS: usize = 20;

/// Commit Sage was built from, when the build could read it from git (or the
/// `SAGE_GIT_COMMIT` environment variable). Uncommitted changes are not recorded.
pub fn git_commit() -> Option<&'static str> {
    option_env!("SAGE_GIT_COMMIT").filter(|commit| !commit.is_empty())
}

#[derive(Debug, Clone, Default, PartialEq, Serialize, Deserialize)]
pub struct Provenance {
    /// Sage Plus crate version.
    pub version: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub git_commit: Option<String>,
    /// Effective, fully resolved configuration, as in `results.json` minus
    /// `output_paths`.
    pub config: serde_json::Value,
    /// Spectrum files, in `mzml_paths` order.
    pub inputs: Vec<InputFile>,
    /// The searched FASTA; `None` for peptide-list-only searches.
    pub fasta: Option<FastaProvenance>,
    /// Peptide lists, custom cleavage sites and PTM libraries.
    pub database_inputs: Vec<InputFile>,
    pub protein_inference: ProteinInference,
}

/// One input file. `sha256` is the hash of the bytes as stored (gzip files are
/// not decompressed), so it matches `sha256sum`.
#[derive(Debug, Clone, Default, PartialEq, Serialize, Deserialize)]
pub struct InputFile {
    /// Spectrum files: the name used in the `filename` output columns.
    /// Database files: `peptides`, `custom_cleavage_sites` or `ptm_library`.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub name: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub role: Option<String>,
    pub path: String,
    pub size_bytes: Option<u64>,
    pub sha256: Option<String>,
    /// Why `sha256` is null.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub sha256_skipped: Option<String>,
}

#[derive(Debug, Clone, Default, PartialEq, Serialize, Deserialize)]
pub struct FastaProvenance {
    pub path: String,
    pub size_bytes: Option<u64>,
    pub sha256: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub sha256_skipped: Option<String>,
    /// Entries in the file.
    pub proteins: usize,
    /// Entries whose accession does not contain the decoy tag.
    pub target_proteins: usize,
    /// Entries whose accession contains the decoy tag.
    pub decoy_proteins: usize,
    pub decoys: DecoyProvenance,
    /// UniProt header fields, when the headers carry them. UniProt FASTA
    /// headers carry no release, so none is recorded.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub uniprot: Option<UniprotHeaders>,
    /// Why the protein counts are missing, if the FASTA could not be reread.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub error: Option<String>,
}

#[derive(Debug, Clone, Default, PartialEq, Serialize, Deserialize)]
pub struct DecoyProvenance {
    /// `generated`: every target peptide is reversed with both terminal
    /// residues kept, and FASTA entries matching the tag are ignored.
    /// `supplied`: FASTA entries matching the tag are searched as decoys.
    pub strategy: String,
    /// Reversal method for generated decoys, `reversed_peptide_keep_termini`.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub method: Option<String>,
    pub decoy_tag: String,
    /// Target and decoy peptides in the searched database, from every source,
    /// after collisions are removed.
    pub target_peptides: usize,
    pub decoy_peptides: usize,
}

#[derive(Debug, Clone, Default, PartialEq, Serialize, Deserialize)]
pub struct UniprotHeaders {
    /// Target entries from Swiss-Prot (`>sp|`) and TrEMBL (`>tr|`).
    pub reviewed: usize,
    pub unreviewed: usize,
    /// Distinct `OS=`/`OX=` pairs over target entries.
    pub organism_count: usize,
    /// The most frequent organisms, at most 20.
    pub organisms: Vec<Organism>,
}

#[derive(Debug, Clone, Default, PartialEq, Eq, Serialize, Deserialize)]
pub struct Organism {
    /// `OS=` name.
    pub name: Option<String>,
    /// `OX=` NCBI taxonomy ID.
    pub taxonomy_id: Option<u64>,
    pub proteins: usize,
}

#[derive(Debug, Clone, Default, PartialEq, Serialize, Deserialize)]
pub struct ProteinInference {
    /// Short strategy name: `idpicker_parsimony` with grouping on,
    /// `protein_lists` with it off (each peptide's protein list is its group).
    pub strategy: String,
    /// Peptide q-value of the first, confident-peptide grouping pass.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub grouping_peptide_q_value: Option<f32>,
    /// How `protein_group_q` is estimated.
    pub protein_group_q: String,
    /// How `protein_q` is estimated.
    pub protein_q: String,
}

impl ProteinInference {
    pub fn from_search(search: &Search) -> Self {
        let grouping = search.protein_grouping;
        Self {
            strategy: if grouping {
                "idpicker_parsimony"
            } else {
                "protein_lists"
            }
            .into(),
            grouping_peptide_q_value: grouping.then_some(search.protein_grouping_peptide_fdr),
            protein_group_q: "picked_group_fdr_unique_peptides".into(),
            protein_q: "picked_protein_fdr_unique_peptides".into(),
        }
    }
}

impl Provenance {
    /// Collect the provenance of a finished search. Reads and hashes the FASTA
    /// and database inputs; hashes spectrum files only when
    /// `record_input_hashes` is set. Never fails: a file that cannot be read
    /// gets a null hash and a reason.
    pub fn collect(search: &Search, database: &Parameters, indexed: &IndexedDatabase) -> Self {
        let mut config = serde_json::to_value(search).unwrap_or(serde_json::Value::Null);
        if let Some(object) = config.as_object_mut() {
            object.remove("output_paths");
        }
        let inputs = search
            .mzml_paths
            .par_iter()
            .map(|url| {
                let mut file = describe_spectrum_file(url, search.record_input_hashes);
                file.name = Some(
                    sage_cloudpath::filename(url)
                        .unwrap_or_else(|| url.as_str())
                        .to_string(),
                );
                file
            })
            .collect();
        let database_inputs = [
            ("peptides", database.peptides.as_deref()),
            (
                "custom_cleavage_sites",
                database.custom_cleavage_sites.as_deref(),
            ),
            (
                "ptm_library",
                database
                    .ptm_library
                    .as_ref()
                    .map(|settings| settings.path.as_str()),
            ),
        ]
        .into_iter()
        .filter_map(|(role, path)| {
            let mut file = describe_path(path?, true);
            file.role = Some(role.into());
            Some(file)
        })
        .collect();
        let fasta = (!database.fasta.is_empty()).then(|| describe_fasta(database, indexed));
        Self {
            version: search.version.clone(),
            git_commit: git_commit().map(str::to_owned),
            config,
            inputs,
            fasta,
            database_inputs,
            protein_inference: ProteinInference::from_search(search),
        }
    }

    /// Key-value pairs for the Parquet footer. `sage.version` and
    /// `sage.git_commit` are plain strings; the others are JSON. The commit key
    /// is left out when the build did not know its commit.
    pub fn parquet_metadata(&self) -> Vec<(String, String)> {
        let json = |value: serde_json::Result<String>| value.unwrap_or_else(|_| "null".into());
        let mut pairs = vec![
            (
                "sage.provenance.version".to_string(),
                PROVENANCE_VERSION.to_string(),
            ),
            ("sage.version".into(), self.version.clone()),
        ];
        if let Some(commit) = &self.git_commit {
            pairs.push(("sage.git_commit".into(), commit.clone()));
        }
        pairs.extend([
            (
                "sage.config".into(),
                json(serde_json::to_string(&self.config)),
            ),
            (
                "sage.inputs".into(),
                json(serde_json::to_string(&self.inputs)),
            ),
            (
                "sage.fasta".into(),
                json(serde_json::to_string(&self.fasta)),
            ),
            (
                "sage.database_inputs".into(),
                json(serde_json::to_string(&self.database_inputs)),
            ),
            (
                "sage.protein_inference".into(),
                json(serde_json::to_string(&self.protein_inference)),
            ),
        ]);
        pairs
    }
}

/// Resolve `path` like the search does, then describe it.
fn describe_path(path: &str, hash: bool) -> InputFile {
    match sage_cloudpath::to_url(path) {
        Ok(url) => describe_file(&url, hash),
        Err(error) => InputFile {
            path: path.into(),
            sha256_skipped: Some(format!("path could not be resolved: {error}")),
            ..InputFile::default()
        },
    }
}

/// Size and, when `hash` is set, the SHA-256 of one file. Local files are
/// hashed by streaming them from disk. Remote spectrum files are not hashed,
/// because that would download them again; remote database files are small
/// and are streamed once more.
fn describe_file(url: &Url, hash: bool) -> InputFile {
    let mut file = InputFile {
        path: url.to_string(),
        ..InputFile::default()
    };
    let local = url.scheme() == "file";
    if local {
        let metadata = url
            .to_file_path()
            .ok()
            .and_then(|path| std::fs::metadata(path).ok());
        if let Some(metadata) = metadata.as_ref().filter(|metadata| metadata.is_file()) {
            file.size_bytes = Some(metadata.len());
        } else if metadata.is_some_and(|metadata| metadata.is_dir()) {
            file.sha256_skipped = Some("directory input; not hashed".into());
            return file;
        }
    }
    if !hash {
        file.sha256_skipped = Some("record_input_hashes is off".into());
        return file;
    }
    match sage_cloudpath::hash::sha256_url(url) {
        Ok(hashed) => {
            file.size_bytes = Some(hashed.size_bytes);
            file.sha256 = Some(hashed.sha256);
        }
        Err(error) => file.sha256_skipped = Some(format!("hashing failed: {error}")),
    }
    file
}

fn describe_spectrum_file(url: &Url, record_input_hashes: bool) -> InputFile {
    if url.scheme() != "file" && record_input_hashes {
        return InputFile {
            path: url.to_string(),
            sha256_skipped: Some("remote spectrum file; not hashed".into()),
            ..InputFile::default()
        };
    }
    describe_file(url, record_input_hashes)
}

fn describe_fasta(database: &Parameters, indexed: &IndexedDatabase) -> FastaProvenance {
    let file = describe_path(&database.fasta, true);
    let (target_peptides, decoy_peptides) =
        indexed
            .peptides
            .iter()
            .fold((0, 0), |(targets, decoys), peptide| match peptide.decoy {
                true => (targets, decoys + 1),
                false => (targets + 1, decoys),
            });
    let mut fasta = FastaProvenance {
        path: file.path,
        size_bytes: file.size_bytes,
        sha256: file.sha256,
        sha256_skipped: file.sha256_skipped,
        decoys: DecoyProvenance {
            strategy: if database.generate_decoys {
                "generated"
            } else {
                "supplied"
            }
            .into(),
            method: database
                .generate_decoys
                .then(|| "reversed_peptide_keep_termini".into()),
            decoy_tag: database.decoy_tag.clone(),
            target_peptides,
            decoy_peptides,
        },
        ..FastaProvenance::default()
    };
    let mut headers = FastaHeaders::new(&database.decoy_tag);
    let read = sage_cloudpath::to_url(&database.fasta)
        .map_err(|error| error.to_string())
        .and_then(|url| {
            sage_cloudpath::hash::for_each_line(url.as_str(), |line| headers.visit(line))
                .map_err(|error| error.to_string())
        });
    match read {
        Ok(()) => headers.finish(&mut fasta),
        Err(error) => fasta.error = Some(format!("FASTA could not be reread: {error}")),
    }
    fasta
}

/// Counts gathered from FASTA header lines.
pub(crate) struct FastaHeaders<'a> {
    decoy_tag: &'a str,
    proteins: usize,
    decoys: usize,
    reviewed: usize,
    unreviewed: usize,
    organisms: HashMap<(Option<String>, Option<u64>), usize>,
}

impl<'a> FastaHeaders<'a> {
    pub(crate) fn new(decoy_tag: &'a str) -> Self {
        Self {
            decoy_tag,
            proteins: 0,
            decoys: 0,
            reviewed: 0,
            unreviewed: 0,
            organisms: HashMap::new(),
        }
    }

    /// Count one FASTA line; only header lines matter. Decoys are recognized
    /// as the search does: the accession (first word) contains the tag.
    pub(crate) fn visit(&mut self, line: &str) {
        let Some(header) = line.trim().strip_prefix('>') else {
            return;
        };
        let accession = header.split_ascii_whitespace().next().unwrap_or_default();
        self.proteins += 1;
        if accession.contains(self.decoy_tag) {
            self.decoys += 1;
            return;
        }
        if accession.starts_with("sp|") {
            self.reviewed += 1;
        } else if accession.starts_with("tr|") {
            self.unreviewed += 1;
        }
        let name = header_field(header, "OS").map(str::to_owned);
        let taxonomy_id = header_field(header, "OX").and_then(|value| value.parse().ok());
        if name.is_some() || taxonomy_id.is_some() {
            *self.organisms.entry((name, taxonomy_id)).or_default() += 1;
        }
    }

    pub(crate) fn finish(self, fasta: &mut FastaProvenance) {
        fasta.proteins = self.proteins;
        fasta.decoy_proteins = self.decoys;
        fasta.target_proteins = self.proteins - self.decoys;
        if self.reviewed + self.unreviewed == 0 && self.organisms.is_empty() {
            return;
        }
        let organism_count = self.organisms.len();
        let mut organisms = self
            .organisms
            .into_iter()
            .map(|((name, taxonomy_id), proteins)| Organism {
                name,
                taxonomy_id,
                proteins,
            })
            .collect::<Vec<_>>();
        organisms.sort_by(|a, b| {
            b.proteins
                .cmp(&a.proteins)
                .then_with(|| a.taxonomy_id.cmp(&b.taxonomy_id))
                .then_with(|| a.name.cmp(&b.name))
        });
        organisms.truncate(MAX_ORGANISMS);
        fasta.uniprot = Some(UniprotHeaders {
            reviewed: self.reviewed,
            unreviewed: self.unreviewed,
            organism_count,
            organisms,
        });
    }
}

/// Value of a UniProt `KEY=value` header field. The value runs to the next
/// ` XX=` field (two uppercase letters) or the end of the header.
fn header_field<'h>(header: &'h str, key: &str) -> Option<&'h str> {
    let marker = format!(" {key}=");
    let start = header.find(&marker)? + marker.len();
    let rest = &header[start..];
    let bytes = rest.as_bytes();
    let end = (0..bytes.len())
        .find(|&index| {
            bytes[index] == b' '
                && bytes.len() > index + 3
                && bytes[index + 1].is_ascii_uppercase()
                && bytes[index + 2].is_ascii_uppercase()
                && bytes[index + 3] == b'='
        })
        .unwrap_or(bytes.len());
    Some(rest[..end].trim()).filter(|value| !value.is_empty())
}

#[cfg(test)]
#[path = "../tests/unit/provenance.rs"]
mod tests;
