//! Protein grouping and inference using IDPicker-style bipartite graph analysis
//!
//! Proteins with identical peptide evidence are collapsed into groups, then a
//! greedy set cover finds an (almost) minimal set of protein groups that
//! explains all observed peptides.
//!
//! Reference: Zhang, B., Chambers, M. C., & Tabb, D. L. (2007). Proteomic
//! parsimony through bipartite graph analysis improves accuracy and transparency.
//! J. Proteome Res., 6(9), 3549-3557. https://doi.org/10.1021/pr070230d

use crate::database::{IndexedDatabase, PeptideIx};
use crate::scoring::Feature;
use fnv::{FnvHashMap, FnvHashSet};
use itertools::Itertools;
use log::info;
use rayon::prelude::*;
use std::sync::Arc;
use std::time::Instant;

/// Compact protein identifier used during grouping
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
struct ProteinIx(u32);

impl ProteinIx {
    fn format(
        &self,
        proteins: &[(Arc<str>, bool)],
        decoy_tag: &str,
        generate_decoys: bool,
    ) -> String {
        let (name, decoy) = &proteins[self.0 as usize];
        if *decoy && generate_decoys {
            format!("{}{}", decoy_tag, name)
        } else {
            name.to_string()
        }
    }
}

/// A group of proteins with identical peptide evidence
#[derive(Debug, Default, PartialEq, Eq, Hash, Clone, PartialOrd, Ord)]
struct ProteinGroup(Vec<ProteinIx>);

impl ProteinGroup {
    fn format(
        &self,
        proteins: &[(Arc<str>, bool)],
        decoy_tag: &str,
        generate_decoys: bool,
    ) -> String {
        self.0
            .iter()
            .map(|ix| ix.format(proteins, decoy_tag, generate_decoys))
            .sorted()
            .join("/")
    }
}

/// Bipartite graph for greedy set cover of proteins (left) to peptides (right)
struct BipartiteGraph {
    edges: Vec<(u32, u32)>,
    /// Degree of each left node at construction time (tiebreaker)
    original_degree: Vec<u32>,
    /// Current degree of each left node
    left_degree: Vec<u32>,
    /// Current degree of each right node
    right_degree: Vec<u32>,
    /// Whether each left node is in the cover
    left_cover: Vec<bool>,
    /// Whether each right node is covered
    right_cover: Vec<bool>,
}

impl BipartiteGraph {
    fn new(edges: Vec<(u32, u32)>, left_count: usize, right_count: usize) -> Self {
        let mut left_degree = vec![0u32; left_count];
        let mut right_degree = vec![0u32; right_count];
        for &(l, r) in &edges {
            left_degree[l as usize] += 1;
            right_degree[r as usize] += 1;
        }
        Self {
            edges,
            original_degree: left_degree.clone(),
            left_degree,
            right_degree,
            left_cover: vec![false; left_count],
            right_cover: vec![false; right_count],
        }
    }

    /// Compute a greedy minimal cover of left nodes that explains all right nodes
    fn into_cover(mut self) -> Vec<bool> {
        while !self.edges.is_empty() {
            self.trim();
            if !self.edges.is_empty() {
                self.add_largest_to_cover();
            }
        }
        self.left_cover
    }

    /// Iteratively add left nodes connected to degree-1 right nodes (unique
    /// peptides force their protein into the cover), then remove all edges
    /// incident to covered nodes.
    fn trim(&mut self) {
        let mut prev_len = 0;
        while prev_len != self.edges.len() {
            prev_len = self.edges.len();

            // Any right node with degree 1 forces its left neighbor into the cover
            for &(l, r) in &self.edges {
                if self.right_degree[r as usize] == 1 {
                    self.left_cover[l as usize] = true;
                }
            }

            // Remove edges where the left node is now covered
            self.edges.retain(|&(l, r)| {
                if self.left_cover[l as usize] {
                    self.right_cover[r as usize] = true;
                    self.left_degree[l as usize] -= 1;
                    self.right_degree[r as usize] -= 1;
                    false
                } else {
                    true
                }
            });

            // Remove edges where the right node is already covered
            self.edges.retain(|&(l, r)| {
                if self.right_cover[r as usize] {
                    self.left_degree[l as usize] -= 1;
                    self.right_degree[r as usize] -= 1;
                    false
                } else {
                    true
                }
            });
        }
    }

    /// Add the left node with the most remaining connections to the cover.
    /// Ties broken by original degree (prefer proteins with more total evidence).
    fn add_largest_to_cover(&mut self) {
        if let Some((idx, _)) = self
            .left_degree
            .iter()
            .zip(&self.original_degree)
            .enumerate()
            .max_by_key(|(_, (remaining, original))| (*remaining, *original))
        {
            self.left_cover[idx] = true;
        }
    }
}

/// Assigns proteins to integer indices and groups them by shared peptide evidence
struct ProteinGrouper {
    /// Map from (protein_name, is_decoy) -> ProteinIx
    protein_index: FnvHashMap<(Arc<str>, bool), ProteinIx>,
    /// Protein groups discovered by collapsing identical evidence
    groups: Vec<ProteinGroup>,
    /// Edges from group index -> meta-peptide index
    edges: Vec<(u32, u32)>,
    /// Number of distinct peptides
    peptide_count: usize,
}

impl ProteinGrouper {
    fn build(db: &IndexedDatabase, peptides: FnvHashSet<PeptideIx>) -> Self {
        let mut protein_index: FnvHashMap<(Arc<str>, bool), ProteinIx> = FnvHashMap::default();

        // Map each peptide to a sorted vector of ProteinIx ("meta-peptide"),
        // deduplicating peptides that map to identical protein sets
        let meta_peptides: FnvHashSet<Vec<ProteinIx>> = peptides
            .into_iter()
            .sorted()
            .map(|pep_ix| {
                let peptide = &db[pep_ix];
                peptide
                    .proteins
                    .iter()
                    .map(|name| {
                        let key = (name.clone(), peptide.decoy);
                        let next_id = ProteinIx(protein_index.len() as u32);
                        *protein_index.entry(key).or_insert(next_id)
                    })
                    .sorted()
                    .collect()
            })
            .collect();

        info!("-  found {} meta peptides", meta_peptides.len(),);

        // Invert: group proteins that share identical meta-peptide sets
        let mut prot_to_metapeps: FnvHashMap<ProteinIx, Vec<usize>> = FnvHashMap::default();
        for (i, meta_pep) in meta_peptides.iter().sorted().enumerate() {
            for &prot_ix in meta_pep {
                prot_to_metapeps.entry(prot_ix).or_default().push(i);
            }
        }

        // Proteins with identical meta-peptide vectors form a group
        let mut evidence_to_group: FnvHashMap<Vec<usize>, ProteinGroup> = FnvHashMap::default();
        for (prot_ix, meta_peps) in prot_to_metapeps {
            evidence_to_group
                .entry(meta_peps)
                .or_default()
                .0
                .push(prot_ix);
        }

        let mut groups = Vec::new();
        let mut edges = Vec::new();
        for (group_idx, (meta_peps, group)) in evidence_to_group.into_iter().sorted().enumerate() {
            groups.push(group);
            for meta_pep_idx in meta_peps {
                edges.push((group_idx as u32, meta_pep_idx as u32));
            }
        }

        info!("-  found {} protein groups", groups.len());

        Self {
            protein_index,
            groups,
            edges,
            peptide_count: meta_peptides.len(),
        }
    }

    /// Run set cover inference and produce a lookup map for annotation
    fn into_group_map(self) -> ProteinGroupLookup {
        let group_count = self.groups.len();
        let cover = BipartiteGraph::new(self.edges, group_count, self.peptide_count).into_cover();

        // Build protein name list ordered by ProteinIx
        let proteins: Vec<(Arc<str>, bool)> = self
            .protein_index
            .into_iter()
            .sorted_by_key(|(_, ix)| ix.0)
            .map(|(key, _)| key)
            .collect();

        // Map each protein to its covered group indices
        let mut protein_to_groups: FnvHashMap<(Arc<str>, bool), Vec<u32>> = FnvHashMap::default();
        for (i, in_cover) in cover.into_iter().enumerate() {
            if !in_cover {
                continue;
            }
            for &prot_ix in &self.groups[i].0 {
                let (name, decoy) = &proteins[prot_ix.0 as usize];
                protein_to_groups
                    .entry((name.clone(), *decoy))
                    .or_default()
                    .push(i as u32);
            }
        }

        ProteinGroupLookup {
            groups: self.groups,
            proteins,
            protein_to_groups,
        }
    }
}

/// Maps individual proteins to their group strings for feature annotation
struct ProteinGroupLookup {
    groups: Vec<ProteinGroup>,
    proteins: Vec<(Arc<str>, bool)>,
    protein_to_groups: FnvHashMap<(Arc<str>, bool), Vec<u32>>,
}

impl ProteinGroupLookup {
    /// Get the sorted, semicolon-delimited protein group string for a peptide
    fn group_string(
        &self,
        peptide: &crate::peptide::Peptide,
        db: &IndexedDatabase,
    ) -> Option<String> {
        let group_set: FnvHashSet<&ProteinGroup> = peptide
            .proteins
            .iter()
            .filter_map(|name| self.protein_to_groups.get(&(name.clone(), peptide.decoy)))
            .flat_map(|indices| indices.iter().map(|&i| &self.groups[i as usize]))
            .collect();

        if group_set.is_empty() {
            return None;
        }

        Some(
            group_set
                .into_iter()
                .map(|g| g.format(&self.proteins, &db.decoy_tag, db.generate_decoys))
                .sorted()
                .join(";"),
        )
    }
}

/// Annotate features with protein group information.
///
/// When `protein_grouping` is enabled, proteins are grouped by shared peptide
/// evidence and a minimal cover is inferred. If `confident_peptide_threshold`
/// is provided, an initial grouping pass is run on only confident peptides,
/// followed by a second pass on all peptides.
///
/// Features not assigned to a group are annotated with their raw protein list.
pub fn generate_protein_groups(
    db: &IndexedDatabase,
    features: &mut [Feature],
    protein_grouping: bool,
    confident_peptide_threshold: Option<f32>,
) {
    let time = Instant::now();
    if protein_grouping {
        if confident_peptide_threshold.is_some() {
            annotate_features(features, db, confident_peptide_threshold);
        }
        annotate_features(features, db, None);
    }

    // Fallback: features not assigned by grouping get their raw protein list
    features
        .par_iter_mut()
        .filter(|f| f.protein_groups.is_none())
        .for_each(|feat| {
            let pep = &db[feat.peptide_idx];
            feat.protein_groups = Some(pep.proteins(&db.decoy_tag, db.generate_decoys));
            feat.num_protein_groups = pep.proteins.len() as u32;
        });
    info!(
        "Grouped and inferred proteins in {:?}ms",
        time.elapsed().as_millis()
    );
}

fn annotate_features(
    features: &mut [Feature],
    db: &IndexedDatabase,
    confident_peptide_threshold: Option<f32>,
) {
    let time = Instant::now();
    let threshold = confident_peptide_threshold.unwrap_or(1.0).clamp(0.0, 1.0);

    let peptides: FnvHashSet<PeptideIx> = features
        .par_iter()
        .filter(|f| f.label != -1 && f.peptide_q < threshold)
        .map(|f| f.peptide_idx)
        .collect();

    info!(
        "Protein grouping: {} unique peptides (threshold={}) in {:?}ms",
        peptides.len(),
        threshold,
        time.elapsed().as_millis()
    );

    let grouper = ProteinGrouper::build(db, peptides);
    let lookup = grouper.into_group_map();

    let annotated: u32 = features
        .par_iter_mut()
        .filter(|f| f.protein_groups.is_none())
        .map(|feat| {
            let pep = &db[feat.peptide_idx];
            match lookup.group_string(pep, db) {
                Some(groups) => {
                    feat.num_protein_groups = groups.matches(';').count() as u32 + 1;
                    feat.protein_groups = Some(groups);
                    1u32
                }
                None => 0,
            }
        })
        .sum();

    info!(
        "-  annotated {} features in {:?}ms",
        annotated,
        time.elapsed().as_millis()
    );
}

#[cfg(test)]
#[path = "../tests/unit/protein_grouping.rs"]
mod test;
