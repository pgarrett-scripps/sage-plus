//! PROSITE-style sequence motifs that restrict where a modification may attach.
//!
//! A motif site is written `motif:<pattern>`. Elements are separated by `-`:
//!
//! * `N` - one literal residue
//! * `x` - any residue
//! * `[ST]` - any listed residue
//! * `{P}` - any residue except those listed
//! * `e(n)` - element `e` repeated exactly `n` times, for example `x(2)`
//! * `<` before the first element anchors the motif at the protein N-terminus
//! * `>` after the last element anchors the motif at the protein C-terminus
//!
//! Exactly one single-width element carries a trailing `*`, marking the modified
//! residue: `N*-{P}-[ST]` is the N-glycosylation sequon, `R-x-x-[ST]*` a basophilic
//! kinase motif, and `C*-x-x-x>` a CaaX prenylation box.
//!
//! Motifs are evaluated against the full source protein of each peptide
//! occurrence. Residues outside the supplied context never match, so without
//! protein context (for example, peptide lists) only in-peptide motifs apply.

use std::collections::HashMap;
use std::fmt::Display;
use std::sync::{Mutex, OnceLock};

use crate::enzyme::Position;
use crate::mass::VALID_AA;

/// Longest accepted pattern, in residues. This bounds pattern size only; the
/// protein context a motif may span is unlimited.
pub const MAX_MOTIF_WIDTH: usize = 64;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
struct ResidueSet(u32);

impl ResidueSet {
    fn bit(residue: u8) -> Option<u32> {
        residue
            .is_ascii_uppercase()
            .then(|| 1u32 << (residue - b'A'))
    }

    fn contains(self, residue: u8) -> bool {
        Self::bit(residue).is_some_and(|bit| self.0 & bit != 0)
    }

    /// Every uppercase letter, so `x` and exclusions also match letters such
    /// as `X` that appear in FASTA sequences but are not searchable residues.
    fn any() -> Self {
        Self((1 << 26) - 1)
    }

    fn valid() -> u32 {
        VALID_AA
            .iter()
            .filter_map(|&r| Self::bit(r))
            .fold(0, |a, b| a | b)
    }

    pub fn residues(self) -> Vec<u8> {
        VALID_AA
            .iter()
            .copied()
            .filter(|&r| self.contains(r))
            .collect()
    }
}

#[derive(Debug)]
pub struct SiteMotif {
    /// Canonical `motif:` spelling, used for identity, ordering, and output.
    canonical: String,
    elements: Vec<ResidueSet>,
    site: usize,
    protein_n_anchor: bool,
    protein_c_anchor: bool,
}

impl PartialEq for SiteMotif {
    fn eq(&self, other: &Self) -> bool {
        self.canonical == other.canonical
    }
}
impl Eq for SiteMotif {}
impl PartialOrd for SiteMotif {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}
impl Ord for SiteMotif {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.canonical.cmp(&other.canonical)
    }
}
impl std::hash::Hash for SiteMotif {
    fn hash<H: std::hash::Hasher>(&self, state: &mut H) {
        self.canonical.hash(state)
    }
}

impl Display for SiteMotif {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(&self.canonical)
    }
}

/// Residues flanking a peptide in its source protein. `left` is ordered as in
/// the protein and ends immediately before the peptide; `right` starts
/// immediately after it. A `*_boundary` flag means the flank reaches the
/// protein terminus, so `<`/`>` anchors may match there.
#[derive(Clone, Copy, Debug, Default)]
pub struct MotifContext<'a> {
    pub left: &'a [u8],
    pub right: &'a [u8],
    pub left_boundary: bool,
    pub right_boundary: bool,
}

impl<'a> MotifContext<'a> {
    /// Peptide without known protein residues. Protein-terminal anchors follow
    /// the digest position.
    pub fn peptide_only(position: Position) -> Self {
        Self::neighbors(&None, &None, position)
    }

    /// One known residue on each side, as recorded by `prev_aa`/`next_aa`.
    pub fn neighbors(prev: &'a Option<u8>, next: &'a Option<u8>, position: Position) -> Self {
        Self {
            left: prev.as_ref().map(std::slice::from_ref).unwrap_or(&[]),
            right: next.as_ref().map(std::slice::from_ref).unwrap_or(&[]),
            left_boundary: prev.is_none() && matches!(position, Position::Nterm | Position::Full),
            right_boundary: next.is_none() && matches!(position, Position::Cterm | Position::Full),
        }
    }

    /// Full protein context around the span `start..start + len`.
    pub fn in_protein(protein: &'a [u8], start: usize, len: usize, position: Position) -> Self {
        let end = start + len;
        Self {
            left: &protein[..start],
            right: &protein[end..],
            left_boundary: start == 0 && matches!(position, Position::Nterm | Position::Full),
            right_boundary: end == protein.len()
                && matches!(position, Position::Cterm | Position::Full),
        }
    }
}

impl SiteMotif {
    /// Parse the pattern after the `motif:` prefix.
    pub fn parse(pattern: &str) -> Result<Self, String> {
        let mut body = pattern.trim();
        let protein_n_anchor = body.starts_with('<');
        if protein_n_anchor {
            body = body[1..].trim_start_matches('-');
        }
        let protein_c_anchor = body.ends_with('>');
        if protein_c_anchor {
            body = body[..body.len() - 1].trim_end_matches('-');
        }
        if body.is_empty() {
            return Err(format!("motif `{pattern}` has no elements"));
        }
        let mut elements = Vec::new();
        let mut site = None;
        let mut canonical = Vec::new();
        for raw in body.split('-') {
            let (token, marked) = match raw.strip_suffix('*') {
                Some(token) => (token, true),
                None => (raw, false),
            };
            let (token, repeat) = match token.split_once('(') {
                Some((token, count)) => {
                    let count = count
                        .strip_suffix(')')
                        .and_then(|c| c.parse::<usize>().ok())
                        .filter(|&c| c > 0)
                        .ok_or_else(|| format!("motif element `{raw}` has an invalid repeat"))?;
                    (token, count)
                }
                None => (token, 1),
            };
            let set = parse_element(token)
                .ok_or_else(|| format!("motif element `{raw}` is not a residue class"))?;
            if marked {
                if repeat != 1 {
                    return Err(format!("modified element `{raw}` must match one residue"));
                }
                if set.residues().is_empty() {
                    return Err(format!("modified element `{raw}` admits no residue"));
                }
                if site.replace(elements.len()).is_some() {
                    return Err(format!(
                        "motif `{pattern}` marks more than one site with `*`"
                    ));
                }
            }
            let mut spelled = canonical_element(set);
            if repeat != 1 {
                spelled.push_str(&format!("({repeat})"));
            }
            if marked {
                spelled.push('*');
            }
            canonical.push(spelled);
            elements.extend(std::iter::repeat_n(set, repeat));
            if elements.len() > MAX_MOTIF_WIDTH {
                return Err(format!(
                    "motif `{pattern}` is wider than {MAX_MOTIF_WIDTH} residues"
                ));
            }
        }
        let site = site.ok_or_else(|| {
            format!("motif `{pattern}` must mark the modified residue with `*`, e.g. N*-{{P}}-[ST]")
        })?;
        let mut canonical = canonical.join("-");
        if protein_n_anchor {
            canonical.insert(0, '<');
        }
        if protein_c_anchor {
            canonical.push('>');
        }
        Ok(Self {
            canonical: format!("motif:{canonical}"),
            elements,
            site,
            protein_n_anchor,
            protein_c_anchor,
        })
    }

    /// Parse and intern, so the result can live inside `Copy` specificities.
    pub fn intern(pattern: &str) -> Result<&'static SiteMotif, String> {
        static TABLE: OnceLock<Mutex<HashMap<String, &'static SiteMotif>>> = OnceLock::new();
        let motif = Self::parse(pattern)?;
        let mut table = TABLE.get_or_init(Default::default).lock().unwrap();
        Ok(*table
            .entry(motif.canonical.clone())
            .or_insert_with(|| Box::leak(Box::new(motif))))
    }

    pub fn canonical(&self) -> &str {
        &self.canonical
    }

    /// Residues that may carry the modification.
    pub fn site_residues(&self) -> Vec<u8> {
        self.elements[self.site].residues()
    }

    /// Residues required before and after the modified residue.
    pub fn reach(&self) -> (usize, usize) {
        (self.site, self.elements.len() - self.site - 1)
    }

    /// Zero-based peptide indices whose residue satisfies the motif.
    pub fn sites(&self, sequence: &[u8], context: MotifContext<'_>) -> Vec<u32> {
        let (before, after) = self.reach();
        let window = context.left.len() + sequence.len() + context.right.len();
        let at = |p: usize| {
            if p < context.left.len() {
                context.left[p]
            } else if p < context.left.len() + sequence.len() {
                sequence[p - context.left.len()]
            } else {
                context.right[p - context.left.len() - sequence.len()]
            }
        };
        (0..sequence.len())
            .filter(|&index| {
                if !self.elements[self.site].contains(sequence[index]) {
                    return false;
                }
                let center = context.left.len() + index;
                let Some(start) = center.checked_sub(before) else {
                    return false;
                };
                let end = center + after + 1;
                if end > window
                    || (self.protein_n_anchor && !(start == 0 && context.left_boundary))
                    || (self.protein_c_anchor && !(end == window && context.right_boundary))
                {
                    return false;
                }
                self.elements
                    .iter()
                    .enumerate()
                    .all(|(offset, set)| set.contains(at(start + offset)))
            })
            .map(|index| index as u32)
            .collect()
    }
}

fn parse_element(token: &str) -> Option<ResidueSet> {
    let residues = |inner: &str| {
        (!inner.is_empty()).then_some(())?;
        inner.bytes().try_fold(0u32, |acc, r| {
            VALID_AA
                .contains(&r)
                .then(|| acc | ResidueSet::bit(r).unwrap())
        })
    };
    let any = ResidueSet::any();
    if token == "x" || token == "X" {
        Some(any)
    } else if let Some(inner) = token.strip_prefix('[').and_then(|t| t.strip_suffix(']')) {
        residues(inner).map(ResidueSet)
    } else if let Some(inner) = token.strip_prefix('{').and_then(|t| t.strip_suffix('}')) {
        residues(inner)
            .map(|excluded| ResidueSet(any.0 & !excluded))
            .filter(|set| set.0 != 0)
    } else if token.len() == 1 {
        residues(token).map(ResidueSet)
    } else {
        None
    }
}

/// Positive classes contain only searchable residues; `x` and exclusions also
/// contain every other letter. Printing by that split reparses identically.
fn canonical_element(set: ResidueSet) -> String {
    let any = ResidueSet::any();
    let excluded = ResidueSet(any.0 & !set.0).residues();
    let negative = set.0 & !ResidueSet(ResidueSet::valid()).0 != 0;
    let listed = |residues: Vec<u8>| String::from_utf8(residues).expect("residues are ASCII");
    match (negative, excluded.is_empty()) {
        (true, true) => "x".into(),
        (true, false) => format!("{{{}}}", listed(excluded)),
        (false, _) => {
            let included = set.residues();
            if included.len() == 1 {
                listed(included)
            } else {
                format!("[{}]", listed(included))
            }
        }
    }
}

#[cfg(test)]
#[path = "../tests/unit/motif.rs"]
mod test;
