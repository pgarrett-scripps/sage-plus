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
//! Motifs are evaluated against the peptide plus whatever protein flanking
//! residues the caller supplies. Positions outside the supplied context never
//! match, so a motif that needs unavailable flanking residues is not satisfied.

use std::collections::HashMap;
use std::fmt::Display;
use std::sync::{Mutex, OnceLock};

use crate::mass::VALID_AA;

/// Longest supported motif, in residues. Keeps per-site checks trivially cheap.
pub const MAX_MOTIF_WIDTH: usize = 16;

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

    fn any() -> Self {
        Self(
            VALID_AA
                .iter()
                .filter_map(|&r| Self::bit(r))
                .fold(0, |a, b| a | b),
        )
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

/// One optional flanking residue as a context slice.
pub fn flank(aa: &Option<u8>) -> &[u8] {
    aa.as_ref().map(std::slice::from_ref).unwrap_or(&[])
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
            elements.extend(std::iter::repeat(set).take(repeat));
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

fn canonical_element(set: ResidueSet) -> String {
    let any = ResidueSet::any();
    let included = set.residues();
    if set == any {
        return "x".into();
    }
    if included.len() == 1 {
        return (included[0] as char).to_string();
    }
    let excluded = ResidueSet(any.0 & !set.0).residues();
    if excluded.len() < included.len() {
        format!("{{{}}}", String::from_utf8(excluded).unwrap())
    } else {
        format!("[{}]", String::from_utf8(included).unwrap())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn ctx<'a>(left: &'a [u8], right: &'a [u8]) -> MotifContext<'a> {
        MotifContext {
            left,
            right,
            ..Default::default()
        }
    }

    #[test]
    fn parses_and_canonicalizes() {
        let motif = SiteMotif::parse("N*-{P}-[TS]").unwrap();
        assert_eq!(motif.canonical(), "motif:N*-{P}-[ST]");
        assert_eq!(motif.reach(), (0, 2));
        assert_eq!(motif.site_residues(), b"N");
        let kinase = SiteMotif::parse("R-x(2)-[ST]*").unwrap();
        assert_eq!(kinase.canonical(), "motif:R-x(2)-[ST]*");
        assert_eq!(kinase.reach(), (3, 0));
        assert!(SiteMotif::parse("N-{P}-[ST]").is_err());
        assert!(SiteMotif::parse("N*-[ST]*").is_err());
        assert!(SiteMotif::parse("x(2)*").is_err());
        assert!(SiteMotif::parse("N*-[ZB]").is_err());
        assert!(SiteMotif::parse("N*-x(20)").is_err());
        assert_eq!(
            SiteMotif::intern("N*-{P}-[ST]").unwrap() as *const _,
            SiteMotif::intern("N*-{P}-[TS]").unwrap() as *const _
        );
    }

    #[test]
    fn sequon_matches_peptide_and_flanks() {
        let motif = SiteMotif::parse("N*-{P}-[ST]").unwrap();
        // NGT sequon, NPS is excluded by the proline rule.
        assert_eq!(motif.sites(b"ANGTANPSK", ctx(b"", b"")), vec![1]);
        // N second-to-last needs one flanking residue after the peptide.
        assert_eq!(motif.sites(b"AANK", ctx(b"", b"")), Vec::<u32>::new());
        assert_eq!(motif.sites(b"AANK", ctx(b"", b"S")), vec![2]);
        assert_eq!(motif.sites(b"AANK", ctx(b"", b"P")), Vec::<u32>::new());
    }

    #[test]
    fn kinase_motif_uses_left_flank() {
        let motif = SiteMotif::parse("R-x-x-[ST]*").unwrap();
        assert_eq!(motif.sites(b"GASPK", ctx(b"R", b"")), vec![2]);
        assert_eq!(motif.sites(b"GASPK", ctx(b"", b"")), Vec::<u32>::new());
    }

    #[test]
    fn protein_terminal_anchor() {
        let motif = SiteMotif::parse("C*-x-x-x>").unwrap();
        let terminal = MotifContext {
            right_boundary: true,
            ..Default::default()
        };
        assert_eq!(motif.sites(b"GKCVLS", terminal), vec![2]);
        assert_eq!(motif.sites(b"GKCVLS", ctx(b"", b"")), Vec::<u32>::new());
        assert_eq!(motif.sites(b"GKCVLSA", terminal), Vec::<u32>::new());
    }
}
