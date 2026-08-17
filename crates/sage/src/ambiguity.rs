//! Peptide sequence-ambiguity annotation.
//!
//! A fragmentation spectrum rarely pins down *every* residue: stretches of
//! sequence with no flanking b/y-ion cleavage evidence could be reordered (or,
//! in an open search, carry a mass shift anywhere within them) without changing
//! the set of matched peaks. This module encodes that uncertainty directly into
//! the peptide string, wrapping ambiguous runs in `(?...)` and — when the
//! precursor mass does not match the peptide's calculated mass — placing the
//! residual mass shift as `[+mass]` (localized to one residue), `(...)[+mass]`
//! (confined to a region) or a leading `{+mass}` (labile / un-localizable).
//!
//! This is a native port of the standalone `SagePeptideAmbiguityAnnotator`
//! Python tool (which itself drives `peptacular.annotate_ambiguity`). The
//! interval-construction, interval-combination and mass-shift-placement logic
//! below mirrors that reference implementation exactly so the column Sage emits
//! matches what the external tool produces; the doctests from the reference are
//! reproduced verbatim as unit tests to lock the port to that behavior.

use crate::peptide::Peptide;

/// PSMs whose `|expmass - calcmass|` is within this many ppm of zero are treated
/// as having no mass shift (matches the reference tool's `--mass_error` default).
pub const DEFAULT_MASS_SHIFT_PPM: f32 = 50.0;

/// Result of annotating a single PSM.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct Ambiguity {
    /// Peptide string with `(?...)` ambiguity intervals and any mass-shift
    /// placement applied, in Sage's `[+mass]` / `[Name]` modification notation.
    pub sequence: String,
    /// The residual mass shift that was placed (0.0 if none).
    pub mass_shift: f32,
}

/// Render a delta mass as Sage's bracketed mod tag, preferring a registered
/// Unimod name. Mirrors [`crate::peptide`]'s `fmt_mod`.
fn mod_tag(mass: f32) -> String {
    match crate::unimod::label_for(mass) {
        Some(name) => format!("[{}]", name),
        None => format!("[{:+}]", mass),
    }
}

/// Labile / un-localized mass shift, rendered with ProForma `{...}` notation.
fn labile_tag(mass: f32) -> String {
    match crate::unimod::label_for(mass) {
        Some(name) => format!("{{{}}}", name),
        None => format!("{{{:+}}}", mass),
    }
}

/// Construct intervals spanning runs of zeros in `counts`.
///
/// Scanning left-to-right (`reverse = false`), a run of zeros is recorded as an
/// inclusive `(start, end)` interval where `end` is the index of the covered
/// residue that terminates the run (or the last index, if the run reaches the
/// end). With `reverse = true` the list is processed right-to-left and the
/// resulting intervals are mapped back into forward coordinates.
fn construct_ambiguity_intervals(counts: &[u16], reverse: bool) -> Vec<(usize, usize)> {
    if reverse {
        let n = counts.len();
        let rev: Vec<u16> = counts.iter().rev().copied().collect();
        let mut intervals: Vec<(usize, usize)> = construct_ambiguity_intervals(&rev, false)
            .into_iter()
            .map(|(start, end)| (n - 1 - end, n - 1 - start))
            .collect();
        intervals.sort_by_key(|&(start, _)| start);
        return intervals;
    }

    let mut intervals: Vec<(usize, usize)> = Vec::new();
    let mut current: Option<(usize, usize)> = None;
    for (i, &cnt) in counts.iter().enumerate() {
        if cnt == 0 {
            current = Some(match current {
                Some((start, _)) => (start, i),
                None => (i, i),
            });
        } else if let Some((start, _)) = current {
            intervals.push((start, i));
            current = None;
        }
    }
    if let Some((start, _)) = current {
        intervals.push((start, counts.len() - 1));
    }
    intervals
}

/// Merge interval lists, keeping only positions ambiguous across *every* list.
///
/// Input and output tuples are `(start, end)` with `end` **exclusive**; this
/// matches the reference implementation, which expands inputs over
/// `range(start, end)` and reconstructs contiguous runs of common indices.
fn combine_ambiguity_intervals(lists: &[Vec<(usize, usize)>]) -> Vec<(usize, usize)> {
    let mut all_indices: std::collections::BTreeSet<usize> = std::collections::BTreeSet::new();
    for list in lists {
        for &(start, end) in list {
            if start != end {
                for i in start..end {
                    all_indices.insert(i);
                }
            }
        }
    }

    let mut common: Vec<usize> = all_indices
        .into_iter()
        .filter(|&idx| {
            lists
                .iter()
                .all(|list| list.iter().any(|&(s, e)| s <= idx && idx < e))
        })
        .collect();
    common.sort_unstable();

    let mut result: Vec<(usize, usize)> = Vec::new();
    if common.is_empty() {
        return result;
    }
    let mut start = common[0];
    for w in common.windows(2) {
        if w[1] > w[0] + 1 {
            result.push((start, w[0] + 1));
            start = w[1];
        }
    }
    result.push((start, *common.last().unwrap() + 1));
    result
}

/// Determine the interval (inclusive) where a mass shift belongs, based on the
/// highest forward-covered position and the lowest reverse-covered position.
/// Returns `None` when forward and reverse coverage overlap (the shift cannot be
/// localized and is therefore labile).
fn mass_shift_interval(forward: &[u16], reverse: &[u16]) -> Option<(usize, usize)> {
    let highest_forward: i64 = forward
        .iter()
        .enumerate()
        .filter(|(_, &c)| c > 0)
        .map(|(i, _)| i as i64)
        .max()
        .unwrap_or(-1);

    let lowest_reverse: i64 = reverse
        .iter()
        .enumerate()
        .filter(|(_, &c)| c > 0)
        .map(|(i, _)| i as i64)
        .min()
        .unwrap_or(reverse.len() as i64);

    if highest_forward >= lowest_reverse {
        return None;
    }
    if highest_forward == lowest_reverse - 1 {
        let p = (highest_forward + 1) as usize;
        return Some((p, p));
    }
    Some((
        (highest_forward + 1) as usize,
        (lowest_reverse - 1) as usize,
    ))
}

/// Where a mass shift is placed during rendering.
enum Shift {
    None,
    Labile(f32),
    /// `[+mass]` on a single residue (0-based index).
    Site(usize, f32),
    /// `(...)[+mass]` over an inclusive residue span.
    Span(usize, usize, f32),
}

/// Annotate `peptide` given per-residue forward (a/b/c) and reverse (x/y/z) ion
/// coverage counts, optionally placing a residual `mass_shift`.
///
/// `forward[i]` counts forward ions mapping to residue `i` (ion ordinal `i + 1`);
/// `reverse[i]` counts reverse ions mapping to residue `i` (ordinal `n - i`).
pub fn annotate(
    peptide: &Peptide,
    forward: &[u16],
    reverse: &[u16],
    mass_shift: Option<f32>,
) -> Ambiguity {
    let n = peptide.sequence.len();
    debug_assert_eq!(forward.len(), n);
    debug_assert_eq!(reverse.len(), n);

    // Ambiguity intervals: positions with neither forward nor reverse evidence.
    let forward_intervals = construct_ambiguity_intervals(forward, false);
    let reverse_intervals = construct_ambiguity_intervals(reverse, true);
    let combined = combine_ambiguity_intervals(&[forward_intervals, reverse_intervals]);
    // Reference applies `Interval(start, end + 1)`, i.e. inclusive residue range
    // `start..=end` for each combined `(start, end)` tuple.
    let amb: Vec<(usize, usize)> = combined.iter().map(|&(s, e)| (s, e)).collect();

    // Mass-shift placement.
    let mut shift = Shift::None;
    if let Some(mass) = mass_shift {
        match mass_shift_interval(forward, reverse) {
            None => shift = Shift::Labile(mass),
            Some((s, e)) if s == e => shift = Shift::Site(s, mass),
            Some((s, e)) => shift = Shift::Span(s, e, mass),
        }
    }

    let sequence = render(peptide, &amb, &shift);
    Ambiguity {
        sequence,
        mass_shift: mass_shift.unwrap_or(0.0),
    }
}

fn render(peptide: &Peptide, amb: &[(usize, usize)], shift: &Shift) -> String {
    let n = peptide.sequence.len();
    let mut out = String::new();

    // Labile shifts are written as a leading `{...}` token.
    if let Shift::Labile(mass) = shift {
        out.push_str(&labile_tag(*mass));
    }

    if let Some(m) = peptide.nterm {
        out.push_str(&mod_tag(m));
        out.push('-');
    }

    // A span mass shift that exactly coincides with an ambiguity interval is
    // attached to that interval (`(?...)[+mass]`) rather than nested separately.
    let (span_start, span_end, span_mass) = match shift {
        Shift::Span(s, e, m) => (Some(*s), Some(*e), Some(*m)),
        _ => (None, None, None),
    };
    let span_merges_amb =
        matches!((span_start, span_end), (Some(s), Some(e)) if amb.contains(&(s, e)));

    for i in 0..n {
        if amb.iter().any(|&(s, _)| s == i) {
            out.push_str("(?");
        }
        if !span_merges_amb && span_start == Some(i) {
            out.push('(');
        }

        out.push(peptide.sequence[i] as char);
        let m = peptide.modification_at(i);
        if m != 0.0 {
            out.push_str(&mod_tag(m));
        }

        if let Shift::Site(pos, mass) = shift {
            if *pos == i {
                out.push_str(&mod_tag(*mass));
            }
        }

        if !span_merges_amb && span_end == Some(i) {
            out.push(')');
            if let Some(mass) = span_mass {
                out.push_str(&mod_tag(mass));
            }
        }

        if amb.iter().any(|&(_, e)| e == i) {
            out.push(')');
            if span_merges_amb && span_end == Some(i) {
                out.push_str(&mod_tag(span_mass.unwrap()));
            }
        }
    }

    if let Some(m) = peptide.cterm {
        out.push('-');
        out.push_str(&mod_tag(m));
    }

    out
}

#[cfg(test)]
mod test {
    use super::*;
    use std::sync::Arc;

    // --- Ported reference doctests: construct_ambiguity_intervals ----------

    #[test]
    fn construct_forward() {
        assert_eq!(
            construct_ambiguity_intervals(&[0, 1, 1, 1, 0, 0, 0], false),
            vec![(0, 1), (4, 6)]
        );
        assert_eq!(
            construct_ambiguity_intervals(&[0, 1, 1, 1, 0, 0, 1], false),
            vec![(0, 1), (4, 6)]
        );
    }

    #[test]
    fn construct_reverse() {
        assert_eq!(
            construct_ambiguity_intervals(&[0, 0, 1, 1, 1, 0, 0], true),
            vec![(0, 1), (4, 6)]
        );
    }

    // --- Ported reference doctests: combine_ambiguity_intervals ------------

    #[test]
    fn combine() {
        assert_eq!(
            combine_ambiguity_intervals(&[vec![(0, 1), (4, 6)], vec![(0, 1)]]),
            vec![(0, 1)]
        );
        assert_eq!(
            combine_ambiguity_intervals(&[vec![(0, 1), (4, 6)], vec![(0, 1), (4, 5)]]),
            vec![(0, 1), (4, 5)]
        );
        assert_eq!(
            combine_ambiguity_intervals(&[vec![(0, 1), (4, 6)], vec![(0, 4), (5, 6)]]),
            vec![(0, 1), (5, 6)]
        );
        assert_eq!(
            combine_ambiguity_intervals(&[vec![(2, 5)], vec![(3, 6)]]),
            vec![(3, 5)]
        );
        assert_eq!(
            combine_ambiguity_intervals(&[vec![(0, 1)], vec![(4, 6)]]),
            Vec::<(usize, usize)>::new()
        );
    }

    // --- Ported reference doctests: mass_shift_interval --------------------

    #[test]
    fn mass_shift() {
        assert_eq!(
            mass_shift_interval(&[1, 1, 1, 0, 0, 0, 0], &[0, 0, 0, 0, 1, 1, 1]),
            Some((3, 3))
        );
        assert_eq!(
            mass_shift_interval(&[1, 1, 1, 0, 0, 0, 0], &[0, 0, 0, 1, 1, 1, 1]),
            Some((3, 3))
        );
        assert_eq!(
            mass_shift_interval(&[1, 1, 0, 0, 0, 0, 0], &[0, 0, 0, 0, 1, 1, 1]),
            Some((2, 3))
        );
        assert_eq!(
            mass_shift_interval(&[0, 0, 0, 0, 0, 0, 0], &[0, 0, 0, 0, 1, 1, 1]),
            Some((0, 3))
        );
        assert_eq!(
            mass_shift_interval(&[1, 1, 1, 0, 0, 0, 0], &[0, 0, 0, 0, 0, 0, 0]),
            Some((3, 6))
        );
        assert_eq!(
            mass_shift_interval(&[1, 1, 1, 1, 1, 0, 0], &[0, 0, 0, 0, 1, 1, 1]),
            None
        );
    }

    // --- End-to-end rendering ---------------------------------------------

    fn peptide(seq: &str) -> Peptide {
        Peptide {
            sequence: Arc::from(seq.as_bytes()),
            modifications: vec![0.0; seq.len()],
            ..Default::default()
        }
    }

    #[test]
    fn fully_covered_has_no_intervals() {
        let p = peptide("PEPTIDE");
        let cov = vec![1u16; 7];
        let a = annotate(&p, &cov, &cov, None);
        assert_eq!(a.sequence, "PEPTIDE");
        assert_eq!(a.mass_shift, 0.0);
    }

    #[test]
    fn ambiguous_nterm_is_wrapped() {
        // The first two residues have neither forward (b1/b2) nor reverse
        // (y6/y5) evidence, so their order is ambiguous and they are wrapped.
        // (reverse[0] is always 0 since no y-ion maps to the N-terminal residue.)
        let p = peptide("PEPTIDE");
        let forward = vec![0, 0, 1, 1, 1, 1, 1];
        let reverse = vec![0, 0, 1, 1, 1, 1, 1];
        let a = annotate(&p, &forward, &reverse, None);
        assert_eq!(a.sequence, "(?PE)PTIDE");
    }

    #[test]
    fn localized_mass_shift_is_bracketed() {
        // Forward stops after residue 2, reverse starts at residue 4 -> the
        // shift localizes to the single gap residue (index 3).
        let p = peptide("PEPTIDE");
        let forward = vec![1, 1, 1, 0, 0, 0, 0];
        let reverse = vec![0, 0, 0, 0, 1, 1, 1];
        let a = annotate(&p, &forward, &reverse, Some(79.96633));
        assert!(a.sequence.contains("T[+79.96633]"), "got: {}", a.sequence);
        assert_eq!(a.mass_shift, 79.96633);
    }

    #[test]
    fn labile_mass_shift_is_prefixed() {
        // Forward and reverse coverage overlap -> shift cannot be localized.
        let p = peptide("PEPTIDE");
        let forward = vec![1, 1, 1, 1, 1, 0, 0];
        let reverse = vec![0, 0, 0, 0, 1, 1, 1];
        let a = annotate(&p, &forward, &reverse, Some(100.0));
        assert!(a.sequence.starts_with("{+100}"), "got: {}", a.sequence);
    }
}
