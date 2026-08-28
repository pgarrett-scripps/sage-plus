use std::collections::HashMap;
use std::sync::Arc;

use crate::unimod;

/// A single PEFF-derived modification site on a protein.
/// `protein_pos` is **0-based** (PEFF on disk is 1-based; we subtract one at parse time).
#[derive(Clone, Debug, PartialEq)]
pub struct PeffMod {
    pub protein_pos: u32,
    pub mass: f32,
    /// Human-readable Unimod name as it appeared in the PEFF entry
    /// (e.g. `"Phospho"`); used to label this delta mass in output.
    pub name: String,
}

type PeffEntries = (Vec<(Arc<str>, String)>, HashMap<Arc<str>, Vec<PeffMod>>);

/// Returns true when the supplied file contents look like a PEFF file
/// (per spec, the first non-empty line is `# PEFF <version>`).
pub fn looks_like_peff(contents: &str) -> bool {
    for line in contents.lines() {
        let t = line.trim();
        if t.is_empty() {
            continue;
        }
        return t.starts_with("# PEFF");
    }
    false
}

/// Parse PEFF contents into the same `(accession, sequence)` shape that
/// `Fasta::parse` produces, plus a per-protein map of `\ModResUnimod`
/// modifications resolved against the embedded Unimod table.
///
/// Lines beginning with `#` (header / database metadata blocks) are skipped.
/// The description-line annotations honored here are limited to
/// `\ModResUnimod=(positions|UNIMOD:N|Name)[(...)]`. Other annotations
/// (`\ModResPsi`, `\ModRes`, `\VariantSimple`, etc.) are intentionally ignored.
pub fn parse(contents: &str, decoy_tag: &str, generate_decoys: bool) -> PeffEntries {
    let mut targets: Vec<(Arc<str>, String)> = Vec::new();
    let mut peff_mods: HashMap<Arc<str>, Vec<PeffMod>> = HashMap::new();

    let mut cur_acc: Option<Arc<str>> = None;
    let mut cur_mods: Vec<PeffMod> = Vec::new();
    let mut cur_seq = String::new();

    let flush = |acc: &mut Option<Arc<str>>,
                 mods: &mut Vec<PeffMod>,
                 seq: &mut String,
                 targets: &mut Vec<(Arc<str>, String)>,
                 peff_mods: &mut HashMap<Arc<str>, Vec<PeffMod>>| {
        if let Some(a) = acc.take() {
            let s = std::mem::take(seq);
            let m = std::mem::take(mods);
            let keep = !a.contains(decoy_tag) || !generate_decoys;
            if keep {
                if !m.is_empty() {
                    peff_mods.insert(a.clone(), m);
                }
                targets.push((a, s));
            }
        } else {
            seq.clear();
            mods.clear();
        }
    };

    for line in contents.lines() {
        if line.is_empty() {
            continue;
        }
        let line = line.trim();
        if line.is_empty() {
            continue;
        }
        // Skip PEFF header / metadata blocks (anything starting with `#`).
        if line.starts_with('#') {
            continue;
        }
        if let Some(desc) = line.strip_prefix('>') {
            // New entry: flush the previous one first.
            flush(
                &mut cur_acc,
                &mut cur_mods,
                &mut cur_seq,
                &mut targets,
                &mut peff_mods,
            );
            let acc_str = desc
                .split_ascii_whitespace()
                .next()
                .unwrap_or("")
                .to_string();
            cur_acc = Some(Arc::from(acc_str));
            cur_mods = parse_mod_res_unimod(desc);
        } else {
            cur_seq.push_str(line);
        }
    }
    flush(
        &mut cur_acc,
        &mut cur_mods,
        &mut cur_seq,
        &mut targets,
        &mut peff_mods,
    );

    (targets, peff_mods)
}

/// Locate `\ModResUnimod=(...)` in a PEFF description line and return the
/// resolved per-residue modifications. Each `(positions|UNIMOD:N|Name)` item
/// may carry multiple comma-separated positions; integer positions are
/// 1-based and become 0-based here. Non-integer positions (`N`, `C`, ...)
/// and unknown accessions are skipped.
fn parse_mod_res_unimod(description: &str) -> Vec<PeffMod> {
    let key = "\\ModResUnimod=";
    let Some(start) = description.find(key) else {
        return Vec::new();
    };
    let body = &description[start + key.len()..];

    let items = split_paren_items(body);
    let mut out = Vec::new();
    for item in items {
        let fields = split_pipe_fields(&item);
        if fields.len() < 3 {
            continue;
        }
        // The first field can carry an optional `<id>:` annotation prefix
        // before the positions; strip everything before the last `:` we
        // see at depth zero.
        let positions_str = match fields[0].rfind(':') {
            Some(idx) => &fields[0][idx + 1..],
            None => fields[0].as_str(),
        };
        let accession_str = fields[1].trim();
        let accession_num = accession_str
            .strip_prefix("UNIMOD:")
            .or_else(|| accession_str.strip_prefix("unimod:"))
            .unwrap_or(accession_str);
        let Ok(accession) = accession_num.trim().parse::<u32>() else {
            continue;
        };
        let Some(mass) = unimod::delta_mass(accession) else {
            log::warn!("PEFF: unknown UNIMOD accession {}", accession);
            continue;
        };
        let raw_name = fields[2].trim();
        // Prefer the canonical Unimod name from the embedded table; fall
        // back to the name string supplied by the PEFF entry.
        let name = unimod::canonical_name(raw_name)
            .map(|s| s.to_string())
            .unwrap_or_else(|| raw_name.to_string());
        for raw_pos in positions_str.split(',') {
            let p = raw_pos.trim();
            let Ok(pos_1based) = p.parse::<u32>() else {
                continue;
            };
            if pos_1based == 0 {
                continue;
            }
            out.push(PeffMod {
                protein_pos: pos_1based - 1,
                mass,
                name: name.clone(),
            });
        }
    }
    out
}

/// Split a string of the form `(item1)(item2)(...)` into its top-level items,
/// honoring nested parentheses (PEFF allows e.g. names that contain `(`/`)`).
fn split_paren_items(s: &str) -> Vec<String> {
    let mut out = Vec::new();
    let bytes = s.as_bytes();
    let mut i = 0;
    while i < bytes.len() {
        match bytes[i] {
            b'(' => {
                let mut depth = 1;
                let start = i + 1;
                i += 1;
                while i < bytes.len() && depth > 0 {
                    match bytes[i] {
                        b'(' => depth += 1,
                        b')' => {
                            depth -= 1;
                            if depth == 0 {
                                break;
                            }
                        }
                        _ => {}
                    }
                    i += 1;
                }
                if depth == 0 && i <= bytes.len() {
                    if let Ok(slice) = std::str::from_utf8(&bytes[start..i]) {
                        out.push(slice.to_string());
                    }
                    i += 1; // consume `)`
                } else {
                    break;
                }
            }
            // Stop at the next `\Key=` (start of a new annotation) or
            // any whitespace that isn't inside parens.
            b'\\' => break,
            b if b.is_ascii_whitespace() => {
                i += 1;
            }
            _ => {
                // Unexpected token between items; bail out.
                break;
            }
        }
    }
    out
}

/// Split an item body on `|` at paren-depth zero.
fn split_pipe_fields(s: &str) -> Vec<String> {
    let mut out = Vec::new();
    let mut depth = 0u32;
    let mut buf = String::new();
    for c in s.chars() {
        match c {
            '(' => {
                depth += 1;
                buf.push(c);
            }
            ')' => {
                depth = depth.saturating_sub(1);
                buf.push(c);
            }
            '|' if depth == 0 => {
                out.push(std::mem::take(&mut buf));
            }
            _ => buf.push(c),
        }
    }
    out.push(buf);
    out
}

#[cfg(test)]
#[path = "../tests/unit/peff.rs"]
mod tests;
