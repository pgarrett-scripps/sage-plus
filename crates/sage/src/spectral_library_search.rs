//! DDA spectrum-to-library candidate indexing, decoy generation, and scoring.

use crate::database::binary_search_slice;
use crate::mass::{Tolerance, PROTON};
use crate::spectral_library::{LibraryFragment, SpectralLibraryEntry};
use crate::spectrum::ProcessedSpectrum;
use serde::{Deserialize, Serialize};

#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(default)]
pub struct LibrarySearchSettings {
    /// Sage spectral-library Parquet or mzSpecLib text file.
    pub path: String,
    /// Prefix used when reporting internally generated decoy proteins.
    pub decoy_tag: String,
    /// Number of deterministic candidate shuffles evaluated per target.
    pub decoy_attempts: usize,
}

impl Default for LibrarySearchSettings {
    fn default() -> Self {
        Self {
            path: String::new(),
            decoy_tag: "rev_".into(),
            decoy_attempts: 32,
        }
    }
}

impl LibrarySearchSettings {
    pub fn validate(&self) -> Result<(), String> {
        if self.path.trim().is_empty() {
            return Err("library_search.path must not be empty".into());
        }
        if self.decoy_tag.is_empty() {
            return Err("library_search.decoy_tag must not be empty".into());
        }
        if self.decoy_attempts == 0 {
            return Err("library_search.decoy_attempts must be greater than zero".into());
        }
        Ok(())
    }
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct DdaLibraryEntry {
    pub library_entry_id: String,
    pub proforma: String,
    pub stripped_peptide: String,
    pub proteins: String,
    pub precursor_charge: u8,
    pub precursor_neutral_mass: f32,
    pub precursor_mz: f32,
    pub retention_time_minutes: f32,
    pub ion_mobility: f32,
    pub source_spectrum_q: f32,
    pub is_decoy: bool,
    pub fragments: Vec<LibraryFragment>,
}

impl DdaLibraryEntry {
    pub fn from_export(entry: &SpectralLibraryEntry) -> Self {
        Self {
            library_entry_id: entry.library_entry_id.clone(),
            proforma: entry.proforma.clone(),
            stripped_peptide: entry.stripped_peptide.clone(),
            proteins: entry.proteins.clone(),
            precursor_charge: entry.precursor_charge,
            precursor_neutral_mass: entry.precursor_neutral_mass,
            precursor_mz: entry.precursor_mz,
            retention_time_minutes: entry.aligned_retention_time_minutes,
            ion_mobility: entry.ion_mobility,
            source_spectrum_q: entry.spectrum_q,
            is_decoy: false,
            fragments: entry.fragments.clone(),
        }
    }
}

#[derive(Clone, Debug, PartialEq)]
struct ResidueToken {
    residue: u8,
    modification: f32,
}

#[derive(Clone, Debug, PartialEq)]
pub struct ParsedPeptidoform {
    pub sequence: Vec<u8>,
    pub modifications: Vec<f32>,
    pub nterm: Option<f32>,
    pub cterm: Option<f32>,
}

fn parse_modification(value: &str) -> Result<f32, String> {
    let value = value.trim();
    if let Ok(mass) = value.parse::<f32>() {
        return Ok(mass);
    }
    let accession = value
        .strip_prefix("UNIMOD:")
        .or_else(|| value.strip_prefix("Unimod:"))
        .and_then(|value| value.parse::<u32>().ok());
    if let Some(accession) = accession {
        return crate::unimod::delta_mass(accession)
            .ok_or_else(|| format!("unknown Unimod accession `{value}`"));
    }
    crate::unimod::mass_by_name(value)
        .ok_or_else(|| format!("unsupported ProForma modification `{value}`"))
}

fn bracket(input: &str, start: usize) -> Result<(&str, usize), String> {
    let rest = input
        .get(start + 1..)
        .ok_or_else(|| "truncated ProForma modification".to_string())?;
    let end = rest
        .find(']')
        .ok_or_else(|| "unterminated ProForma modification".to_string())?;
    Ok((&rest[..end], start + end + 2))
}

/// Parse residue, N-terminal, and C-terminal mass/Unimod modifications.
pub fn parse_proforma(input: &str) -> Result<ParsedPeptidoform, String> {
    let input = input
        .rsplit_once('/')
        .filter(|(_, charge)| charge.parse::<u8>().is_ok())
        .map_or(input, |(peptidoform, _)| peptidoform);
    let bytes = input.as_bytes();
    let mut cursor = 0usize;
    let mut nterm = None;
    let mut cterm = None;
    let mut tokens = Vec::<ResidueToken>::new();

    if bytes.first() == Some(&b'[') {
        let (value, next) = bracket(input, 0)?;
        if bytes.get(next) != Some(&b'-') {
            return Err("an initial ProForma modification must be N-terminal (`[mod]-`)".into());
        }
        nterm = Some(parse_modification(value)?);
        cursor = next + 1;
    }

    while cursor < bytes.len() {
        if bytes[cursor] == b'-' && bytes.get(cursor + 1) == Some(&b'[') {
            let (value, next) = bracket(input, cursor + 1)?;
            cterm = Some(parse_modification(value)?);
            cursor = next;
            break;
        }
        let residue = bytes[cursor];
        if crate::mass::monoisotopic(residue) == 0.0 {
            return Err(format!("unsupported ProForma token at byte {cursor}"));
        }
        cursor += 1;
        let mut modification = 0.0;
        if bytes.get(cursor) == Some(&b'[') {
            let (value, next) = bracket(input, cursor)?;
            modification = parse_modification(value)?;
            cursor = next;
        }
        tokens.push(ResidueToken {
            residue,
            modification,
        });
    }
    if cursor != bytes.len() || tokens.is_empty() {
        return Err("invalid or empty ProForma peptidoform".into());
    }
    Ok(ParsedPeptidoform {
        sequence: tokens.iter().map(|token| token.residue).collect(),
        modifications: tokens.iter().map(|token| token.modification).collect(),
        nterm,
        cterm,
    })
}

fn render_proforma(tokens: &[ResidueToken], nterm: Option<f32>, cterm: Option<f32>) -> String {
    use std::fmt::Write;
    let mut output = String::new();
    if let Some(mass) = nterm {
        write!(output, "[{mass:+}]-").expect("write to String");
    }
    for token in tokens {
        output.push(token.residue as char);
        if token.modification != 0.0 {
            write!(output, "[{:+}]", token.modification).expect("write to String");
        }
    }
    if let Some(mass) = cterm {
        write!(output, "-[{mass:+}]").expect("write to String");
    }
    output
}

fn fragment_token_mass(
    tokens: &[ResidueToken],
    kind: crate::ion_series::Kind,
    ordinal: usize,
) -> f32 {
    let take = ordinal.min(tokens.len());
    let sum = |tokens: &[ResidueToken]| {
        tokens
            .iter()
            .map(|token| crate::mass::monoisotopic(token.residue) + token.modification)
            .sum()
    };
    match kind {
        crate::ion_series::Kind::A | crate::ion_series::Kind::B | crate::ion_series::Kind::C => {
            sum(&tokens[..take])
        }
        crate::ion_series::Kind::X | crate::ion_series::Kind::Y | crate::ion_series::Kind::Z => {
            sum(&tokens[tokens.len() - take..])
        }
    }
}

fn shuffled_tokens(tokens: &[ResidueToken], attempt: usize) -> Vec<ResidueToken> {
    let mut shuffled = tokens.to_vec();
    if shuffled.len() <= 2 {
        shuffled.reverse();
        return shuffled;
    }
    let internal = &tokens[1..tokens.len() - 1];
    let shift = 1 + attempt % internal.len();
    let reverse = (attempt / internal.len()) % 2 == 1;
    for index in 0..internal.len() {
        let source = if reverse {
            internal.len() - 1 - ((index + shift) % internal.len())
        } else {
            (index + shift) % internal.len()
        };
        shuffled[index + 1] = internal[source].clone();
    }
    shuffled
}

fn make_decoy(
    entry: &DdaLibraryEntry,
    attempts: usize,
    decoy_tag: &str,
) -> Result<DdaLibraryEntry, String> {
    let parsed = parse_proforma(&entry.proforma)
        .map_err(|error| format!("library entry `{}`: {error}", entry.library_entry_id))?;
    let target = parsed
        .sequence
        .iter()
        .copied()
        .zip(parsed.modifications.iter().copied())
        .map(|(residue, modification)| ResidueToken {
            residue,
            modification,
        })
        .collect::<Vec<_>>();
    let mut best: Option<(usize, Vec<ResidueToken>, Vec<LibraryFragment>)> = None;
    for attempt in 0..attempts {
        let candidate = shuffled_tokens(&target, attempt);
        if candidate == target {
            continue;
        }
        let fragments = entry
            .fragments
            .iter()
            .map(|fragment| {
                let ordinal = usize::try_from(fragment.ordinal).unwrap_or_default();
                let delta = fragment_token_mass(&candidate, fragment.kind, ordinal)
                    - fragment_token_mass(&target, fragment.kind, ordinal);
                let mut fragment = fragment.clone();
                fragment.mz += delta / fragment.charge as f32;
                fragment
            })
            .collect::<Vec<_>>();
        let overlap = fragments
            .iter()
            .filter(|decoy| {
                entry.fragments.iter().any(|target| {
                    target.kind == decoy.kind
                        && target.charge == decoy.charge
                        && ((target.mz - decoy.mz) / target.mz * 1_000_000.0).abs() <= 20.0
                })
            })
            .count();
        if best.as_ref().is_none_or(|(score, _, _)| overlap < *score) {
            best = Some((overlap, candidate, fragments));
        }
    }
    let (_, tokens, mut fragments) = best.ok_or_else(|| {
        format!(
            "library entry `{}` cannot produce a distinct composition-preserving decoy",
            entry.library_entry_id
        )
    })?;
    fragments.sort_unstable_by(|left, right| left.mz.total_cmp(&right.mz));
    let proforma = render_proforma(&tokens, parsed.nterm, parsed.cterm);
    Ok(DdaLibraryEntry {
        library_entry_id: format!("{decoy_tag}{}", entry.library_entry_id),
        stripped_peptide: tokens.iter().map(|token| token.residue as char).collect(),
        proforma,
        is_decoy: true,
        fragments,
        ..entry.clone()
    })
}

/// Add one deterministic, mass-preserving shuffled decoy for every target.
pub fn generate_decoys(
    targets: Vec<DdaLibraryEntry>,
    settings: &LibrarySearchSettings,
) -> Result<Vec<DdaLibraryEntry>, String> {
    settings.validate()?;
    let mut entries = Vec::with_capacity(targets.len() * 2);
    for mut target in targets {
        target.is_decoy = false;
        let decoy = make_decoy(&target, settings.decoy_attempts, &settings.decoy_tag)?;
        entries.push(target);
        entries.push(decoy);
    }
    Ok(entries)
}

fn parse_peak_annotation(annotation: &str) -> Option<(crate::ion_series::Kind, i32, i32, f32)> {
    let annotation = annotation.trim_matches('"').trim();
    let kind = match annotation.as_bytes().first().copied()? {
        b'a' => crate::ion_series::Kind::A,
        b'b' => crate::ion_series::Kind::B,
        b'c' => crate::ion_series::Kind::C,
        b'x' => crate::ion_series::Kind::X,
        b'y' => crate::ion_series::Kind::Y,
        b'z' => crate::ion_series::Kind::Z,
        _ => return None,
    };
    let mut rest = &annotation[1..];
    let ordinal_end = rest
        .find(|character: char| !character.is_ascii_digit())
        .unwrap_or(rest.len());
    let ordinal = rest[..ordinal_end].parse().ok()?;
    rest = &rest[ordinal_end..];
    let (rest, charge) = match rest.rsplit_once('^') {
        Some((rest, charge)) => (rest, charge.parse().ok()?),
        None => (rest, 1),
    };
    let neutral_loss = match rest {
        "" => 0.0,
        "-H2O" => 18.010_565,
        "-NH3" => 17.026_55,
        "-H3PO4" => 97.976_9,
        _ => return None,
    };
    Some((kind, ordinal, charge, neutral_loss))
}

fn property_value<'a>(line: &'a str, accession: &str) -> Option<&'a str> {
    let property = if line.starts_with('[') {
        line.find(']').and_then(|end| line.get(end + 1..))?
    } else {
        line
    };
    property
        .starts_with(accession)
        .then(|| property.split_once('=').map(|(_, value)| value.trim()))
        .flatten()
}

fn finish_mzspeclib_entry(
    entries: &mut Vec<DdaLibraryEntry>,
    mut entry: DdaLibraryEntry,
) -> Result<(), String> {
    if entry.library_entry_id.is_empty() {
        return Ok(());
    }
    if entry.proforma.is_empty() || entry.precursor_charge == 0 {
        return Err(format!(
            "mzSpecLib entry `{}` is missing its ProForma analyte or charge",
            entry.library_entry_id
        ));
    }
    let parsed = parse_proforma(&entry.proforma)
        .map_err(|error| format!("mzSpecLib entry `{}`: {error}", entry.library_entry_id))?;
    if entry.stripped_peptide.is_empty() {
        entry.stripped_peptide = String::from_utf8(parsed.sequence)
            .map_err(|_| "mzSpecLib peptide sequence is not UTF-8".to_string())?;
    }
    if entry.precursor_neutral_mass <= 0.0 && entry.precursor_mz > 0.0 {
        entry.precursor_neutral_mass =
            (entry.precursor_mz - PROTON) * f32::from(entry.precursor_charge);
    }
    if entry.precursor_mz <= 0.0 && entry.precursor_neutral_mass > 0.0 {
        entry.precursor_mz =
            entry.precursor_neutral_mass / f32::from(entry.precursor_charge) + PROTON;
    }
    let max_intensity = entry
        .fragments
        .iter()
        .map(|fragment| fragment.relative_intensity)
        .max_by(f32::total_cmp)
        .unwrap_or_default();
    if max_intensity > 0.0 {
        for fragment in &mut entry.fragments {
            fragment.relative_intensity /= max_intensity;
        }
    }
    entry
        .fragments
        .sort_unstable_by(|left, right| left.mz.total_cmp(&right.mz));
    entries.push(entry);
    Ok(())
}

/// Parse mzSpecLib 1.0 text libraries. Unsupported peak annotations are
/// ignored; entries must still retain at least one supported fragment.
pub fn deserialize_mzspeclib(text: &str) -> Result<Vec<DdaLibraryEntry>, String> {
    let mut entries = Vec::new();
    let mut entry = DdaLibraryEntry::default();
    let mut proteins = Vec::<String>::new();
    let mut in_peaks = false;

    for raw_line in text.lines() {
        let line = raw_line.trim();
        if line.starts_with("<Spectrum=") {
            if !proteins.is_empty() {
                entry.proteins = proteins.join(";");
            }
            finish_mzspeclib_entry(&mut entries, entry)?;
            entry = DdaLibraryEntry::default();
            proteins = Vec::new();
            entry.library_entry_id = line
                .trim_start_matches("<Spectrum=")
                .trim_end_matches('>')
                .to_string();
            in_peaks = false;
            continue;
        }
        if line == "<Peaks>" {
            in_peaks = true;
            continue;
        }
        if line.starts_with('<') || line.is_empty() {
            in_peaks = false;
            continue;
        }
        if in_peaks {
            let fields = line.split_whitespace().collect::<Vec<_>>();
            if fields.len() >= 3 {
                if let (Ok(mz), Ok(intensity), Some((kind, ordinal, charge, neutral_loss))) = (
                    fields[0].parse::<f32>(),
                    fields[1].parse::<f32>(),
                    parse_peak_annotation(fields[2]),
                ) {
                    entry.fragments.push(LibraryFragment {
                        kind,
                        ordinal,
                        charge,
                        neutral_loss,
                        mz,
                        relative_intensity: intensity,
                    });
                }
            }
            continue;
        }
        if let Some(value) = property_value(line, "MS:1003061|") {
            entry.library_entry_id = value.to_string();
        } else if let Some(value) = property_value(line, "MS:1003270|") {
            let (proforma, charge) = value
                .rsplit_once('/')
                .ok_or_else(|| format!("mzSpecLib analyte `{value}` does not include a charge"))?;
            entry.proforma = proforma.to_string();
            entry.precursor_charge = charge
                .parse()
                .map_err(|_| format!("invalid mzSpecLib charge `{charge}`"))?;
        } else if let Some(value) = property_value(line, "MS:1000041|") {
            entry.precursor_charge = value
                .parse()
                .map_err(|_| format!("invalid mzSpecLib charge `{value}`"))?;
        } else if let Some(value) = property_value(line, "MS:1003208|") {
            entry.precursor_mz = value
                .parse()
                .map_err(|_| format!("invalid mzSpecLib precursor m/z `{value}`"))?;
        } else if let Some(value) = property_value(line, "MS:1003053|") {
            entry.precursor_mz = value
                .parse()
                .map_err(|_| format!("invalid mzSpecLib precursor m/z `{value}`"))?;
        } else if let Some(value) = property_value(line, "MS:1001117|") {
            entry.precursor_neutral_mass = value
                .parse()
                .map_err(|_| format!("invalid mzSpecLib precursor mass `{value}`"))?;
        } else if let Some(value) = property_value(line, "MS:1000894|") {
            entry.retention_time_minutes = value
                .parse()
                .map_err(|_| format!("invalid mzSpecLib retention time `{value}`"))?;
        } else if let Some(value) = property_value(line, "MS:1002815|") {
            entry.ion_mobility = value
                .parse()
                .map_err(|_| format!("invalid mzSpecLib ion mobility `{value}`"))?;
        } else if let Some(value) = property_value(line, "MS:1002354|") {
            entry.source_spectrum_q = value
                .parse()
                .map_err(|_| format!("invalid mzSpecLib q-value `{value}`"))?;
        } else if let Some(value) = property_value(line, "MS:1000885|") {
            proteins.push(value.to_string());
        }
    }
    if !proteins.is_empty() {
        entry.proteins = proteins.join(";");
    }
    finish_mzspeclib_entry(&mut entries, entry)?;
    DdaLibraryIndex::new(entries).map(|index| index.entries)
}

#[derive(Clone, Debug, PartialEq)]
pub struct DdaLibraryMatch {
    pub entry_index: usize,
    pub matched_peaks: usize,
    pub spectral_angle: f32,
    pub explained_library_intensity: f32,
    pub explained_query_intensity: f32,
    pub precursor_ppm: f32,
    pub retention_time_delta_minutes: f32,
    pub ion_mobility_delta: Option<f32>,
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DdaLibrarySearchParameters {
    pub precursor_tolerance: Tolerance,
    pub fragment_tolerance: Tolerance,
    pub min_matched_peaks: usize,
    pub max_hits: usize,
}

impl Default for DdaLibrarySearchParameters {
    fn default() -> Self {
        Self {
            precursor_tolerance: Tolerance::Ppm(-10.0, 10.0),
            fragment_tolerance: Tolerance::Ppm(-10.0, 10.0),
            min_matched_peaks: 6,
            max_hits: 1,
        }
    }
}

/// Precursor-mass-sorted empirical library used by the DDA scoring spike.
#[derive(Clone, Debug, Default)]
pub struct DdaLibraryIndex {
    entries: Vec<DdaLibraryEntry>,
}

impl DdaLibraryIndex {
    pub fn new(mut entries: Vec<DdaLibraryEntry>) -> Result<Self, String> {
        let mut identifiers = std::collections::HashSet::with_capacity(entries.len());
        for entry in &entries {
            if entry.library_entry_id.is_empty() {
                return Err("library entry identifiers must not be empty".into());
            }
            if !identifiers.insert(entry.library_entry_id.as_str()) {
                return Err(format!(
                    "duplicate library entry identifier `{}`",
                    entry.library_entry_id
                ));
            }
            if entry.precursor_charge == 0 {
                return Err(format!(
                    "library entry `{}` has precursor charge zero",
                    entry.library_entry_id
                ));
            }
            if !entry.precursor_neutral_mass.is_finite() || entry.precursor_neutral_mass <= 0.0 {
                return Err(format!(
                    "library entry `{}` has an invalid precursor neutral mass",
                    entry.library_entry_id
                ));
            }
            if entry.fragments.is_empty() {
                return Err(format!(
                    "library entry `{}` has no fragment peaks",
                    entry.library_entry_id
                ));
            }
            if entry.fragments.iter().any(|fragment| {
                fragment.charge <= 0
                    || !fragment.mz.is_finite()
                    || fragment.mz <= PROTON
                    || !fragment.relative_intensity.is_finite()
                    || fragment.relative_intensity <= 0.0
            }) {
                return Err(format!(
                    "library entry `{}` has an invalid fragment peak",
                    entry.library_entry_id
                ));
            }
        }
        entries.sort_unstable_by(|left, right| {
            left.precursor_neutral_mass
                .total_cmp(&right.precursor_neutral_mass)
                .then_with(|| left.precursor_charge.cmp(&right.precursor_charge))
                .then_with(|| left.library_entry_id.cmp(&right.library_entry_id))
        });
        Ok(Self { entries })
    }

    pub fn entries(&self) -> &[DdaLibraryEntry] {
        &self.entries
    }

    pub fn candidates(
        &self,
        precursor_neutral_mass: f32,
        precursor_charge: u8,
        precursor_tolerance: Tolerance,
    ) -> impl Iterator<Item = (usize, &DdaLibraryEntry)> {
        let (low, high) = precursor_tolerance.bounds(precursor_neutral_mass);
        let left = self
            .entries
            .partition_point(|entry| entry.precursor_neutral_mass < low);
        let right = self
            .entries
            .partition_point(|entry| entry.precursor_neutral_mass <= high);
        self.entries[left..right]
            .iter()
            .enumerate()
            .filter(move |(_, entry)| entry.precursor_charge == precursor_charge)
            .map(move |(offset, entry)| (left + offset, entry))
    }

    pub fn search(
        &self,
        query: &ProcessedSpectrum,
        precursor_neutral_mass: f32,
        precursor_charge: u8,
        parameters: DdaLibrarySearchParameters,
    ) -> Vec<DdaLibraryMatch> {
        if parameters.max_hits == 0 || query.masses.is_empty() || query.intensities.is_empty() {
            return Vec::new();
        }

        let query_intensity_sum = query
            .intensities
            .iter()
            .copied()
            .filter(|intensity| intensity.is_finite() && *intensity > 0.0)
            .sum::<f32>();
        if query_intensity_sum <= 0.0 {
            return Vec::new();
        }

        let query_mobility = query
            .precursors
            .first()
            .and_then(|precursor| precursor.inverse_ion_mobility)
            .filter(|value| value.is_finite() && *value > 0.0);
        let mut matches = self
            .candidates(
                precursor_neutral_mass,
                precursor_charge,
                parameters.precursor_tolerance,
            )
            .filter_map(|(entry_index, entry)| {
                score_entry(
                    query,
                    query_intensity_sum,
                    query_mobility,
                    precursor_neutral_mass,
                    entry_index,
                    entry,
                    parameters.fragment_tolerance,
                    parameters.min_matched_peaks,
                )
            })
            .collect::<Vec<_>>();
        matches.sort_unstable_by(|left, right| {
            right
                .spectral_angle
                .total_cmp(&left.spectral_angle)
                .then_with(|| right.matched_peaks.cmp(&left.matched_peaks))
                .then_with(|| {
                    left.precursor_ppm
                        .abs()
                        .total_cmp(&right.precursor_ppm.abs())
                })
                .then_with(|| left.entry_index.cmp(&right.entry_index))
        });
        matches.truncate(parameters.max_hits);
        matches
    }
}

#[allow(clippy::too_many_arguments)]
fn score_entry(
    query: &ProcessedSpectrum,
    query_intensity_sum: f32,
    query_mobility: Option<f32>,
    precursor_neutral_mass: f32,
    entry_index: usize,
    entry: &DdaLibraryEntry,
    fragment_tolerance: Tolerance,
    min_matched_peaks: usize,
) -> Option<DdaLibraryMatch> {
    // Potential assignments are sorted by their contribution and greedily
    // accepted to enforce a one-to-one mapping between library and query peaks.
    let mut assignments = Vec::new();
    for (library_index, fragment) in entry.fragments.iter().enumerate() {
        let center = (fragment.mz - PROTON) * fragment.charge as f32;
        let (low, high) = fragment_tolerance.bounds(center);
        let (left, right) = binary_search_slice(
            &query.masses,
            |mass, bound| mass.total_cmp(bound),
            low,
            high,
        );
        for query_index in left..right {
            let query_mass = query.masses[query_index];
            let query_intensity = query.intensities[query_index];
            if query_mass >= low
                && query_mass <= high
                && query_intensity.is_finite()
                && query_intensity > 0.0
            {
                assignments.push((
                    (fragment.relative_intensity * query_intensity).sqrt(),
                    library_index,
                    query_index,
                ));
            }
        }
    }
    assignments.sort_unstable_by(|left, right| right.0.total_cmp(&left.0));

    let mut used_library = vec![false; entry.fragments.len()];
    let mut used_query = vec![false; query.masses.len()];
    let mut numerator = 0.0f32;
    let mut matched_library_intensity = 0.0f32;
    let mut matched_query_intensity = 0.0f32;
    let mut matched_peaks = 0usize;
    for (contribution, library_index, query_index) in assignments {
        if used_library[library_index] || used_query[query_index] {
            continue;
        }
        used_library[library_index] = true;
        used_query[query_index] = true;
        numerator += contribution;
        matched_library_intensity += entry.fragments[library_index].relative_intensity;
        matched_query_intensity += query.intensities[query_index];
        matched_peaks += 1;
    }
    if matched_peaks < min_matched_peaks {
        return None;
    }

    let library_intensity_sum = entry
        .fragments
        .iter()
        .map(|fragment| fragment.relative_intensity)
        .sum::<f32>();
    let cosine = (numerator / (library_intensity_sum * query_intensity_sum).sqrt()).clamp(0.0, 1.0);
    let spectral_angle = 1.0 - 2.0 * cosine.acos() / std::f32::consts::PI;
    let precursor_ppm = (precursor_neutral_mass - entry.precursor_neutral_mass) * 1_000_000.0
        / entry.precursor_neutral_mass;
    let ion_mobility_delta = query_mobility.and_then(|observed| {
        (entry.ion_mobility.is_finite() && entry.ion_mobility > 0.0)
            .then_some(observed - entry.ion_mobility)
    });

    Some(DdaLibraryMatch {
        entry_index,
        matched_peaks,
        spectral_angle,
        explained_library_intensity: matched_library_intensity / library_intensity_sum,
        explained_query_intensity: matched_query_intensity / query_intensity_sum,
        precursor_ppm,
        retention_time_delta_minutes: query.scan_start_time - entry.retention_time_minutes,
        ion_mobility_delta,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ion_series::Kind;
    use crate::spectrum::Precursor;

    fn fragment(mz: f32, relative_intensity: f32) -> LibraryFragment {
        LibraryFragment {
            kind: Kind::Y,
            ordinal: 1,
            charge: 1,
            neutral_loss: 0.0,
            mz,
            relative_intensity,
        }
    }

    fn entry(id: &str, charge: u8, fragments: Vec<LibraryFragment>) -> DdaLibraryEntry {
        DdaLibraryEntry {
            library_entry_id: id.into(),
            proforma: id.into(),
            stripped_peptide: id.into(),
            proteins: "P1".into(),
            precursor_charge: charge,
            precursor_neutral_mass: 1_000.0,
            precursor_mz: 1_000.0 / charge as f32 + PROTON,
            retention_time_minutes: 10.0,
            ion_mobility: 1.1,
            source_spectrum_q: 0.001,
            is_decoy: false,
            fragments,
        }
    }

    fn query(peaks: &[(f32, f32)]) -> ProcessedSpectrum {
        ProcessedSpectrum {
            level: 2,
            scan_start_time: 10.5,
            precursors: vec![Precursor {
                inverse_ion_mobility: Some(1.2),
                ..Default::default()
            }],
            masses: peaks.iter().map(|(mz, _)| mz - PROTON).collect(),
            intensities: peaks.iter().map(|(_, intensity)| *intensity).collect(),
            charges: vec![1; peaks.len()],
            ..Default::default()
        }
    }

    #[test]
    fn matching_intensity_pattern_beats_mass_matched_distractor() {
        let index = DdaLibraryIndex::new(vec![
            entry(
                "matching",
                2,
                vec![fragment(200.0, 1.0), fragment(300.0, 0.25)],
            ),
            entry(
                "distractor",
                2,
                vec![fragment(200.0, 0.1), fragment(300.0, 1.0)],
            ),
        ])
        .unwrap();
        let query = query(&[(200.0, 100.0), (300.0, 25.0)]);
        let matches = index.search(
            &query,
            1_000.0,
            2,
            DdaLibrarySearchParameters {
                min_matched_peaks: 2,
                max_hits: 2,
                ..Default::default()
            },
        );
        assert_eq!(matches.len(), 2);
        assert_eq!(
            index.entries()[matches[0].entry_index].library_entry_id,
            "matching"
        );
        assert!((matches[0].spectral_angle - 1.0).abs() < 1e-5);
        assert!(matches[0].spectral_angle > matches[1].spectral_angle);
        assert_eq!(matches[0].ion_mobility_delta, Some(0.100_000_024));
    }

    #[test]
    fn precursor_charge_filters_candidates() {
        let index = DdaLibraryIndex::new(vec![
            entry("charge-2", 2, vec![fragment(200.0, 1.0)]),
            entry("charge-3", 3, vec![fragment(200.0, 1.0)]),
        ])
        .unwrap();
        let query = query(&[(200.0, 100.0)]);
        let matches = index.search(
            &query,
            1_000.0,
            3,
            DdaLibrarySearchParameters {
                precursor_tolerance: Tolerance::Da(-0.1, 0.1),
                min_matched_peaks: 1,
                max_hits: 10,
                ..Default::default()
            },
        );
        assert_eq!(matches.len(), 1);
        assert_eq!(
            index.entries()[matches[0].entry_index].library_entry_id,
            "charge-3"
        );
    }

    #[test]
    fn one_query_peak_cannot_match_two_library_peaks() {
        let index = DdaLibraryIndex::new(vec![entry(
            "overlap",
            2,
            vec![fragment(200.000, 1.0), fragment(200.001, 0.5)],
        )])
        .unwrap();
        let query = query(&[(200.0005, 100.0)]);
        let matches = index.search(
            &query,
            1_000.0,
            2,
            DdaLibrarySearchParameters {
                precursor_tolerance: Tolerance::Da(-0.1, 0.1),
                fragment_tolerance: Tolerance::Da(-0.01, 0.01),
                min_matched_peaks: 1,
                max_hits: 1,
            },
        );
        assert_eq!(matches[0].matched_peaks, 1);
    }

    #[test]
    fn invalid_entries_are_rejected() {
        let invalid = entry("empty", 2, Vec::new());
        assert!(DdaLibraryIndex::new(vec![invalid])
            .unwrap_err()
            .contains("no fragment peaks"));
    }

    #[test]
    fn parses_mass_and_unimod_proforma() {
        let parsed = parse_proforma("[+42.010565]-AC[UNIMOD:4]DM[Oxidation]-[+1.0]/2").unwrap();
        assert_eq!(parsed.sequence, b"ACDM");
        assert!((parsed.nterm.unwrap() - 42.010_565).abs() < 1e-5);
        assert!((parsed.modifications[1] - 57.021_465).abs() < 1e-4);
        assert!((parsed.modifications[3] - 15.994_915).abs() < 1e-4);
        assert_eq!(parsed.cterm, Some(1.0));
    }

    #[test]
    fn shuffled_decoys_preserve_precursor_and_intensities() {
        let target = DdaLibraryEntry {
            library_entry_id: "PEPTIDER/2".into(),
            proforma: "PEPTIDER".into(),
            stripped_peptide: "PEPTIDER".into(),
            proteins: "P1;P2".into(),
            precursor_charge: 2,
            precursor_neutral_mass: 955.461,
            precursor_mz: 478.7378,
            fragments: vec![
                LibraryFragment {
                    kind: Kind::B,
                    ordinal: 3,
                    charge: 1,
                    neutral_loss: 0.0,
                    mz: 324.155,
                    relative_intensity: 1.0,
                },
                LibraryFragment {
                    kind: Kind::Y,
                    ordinal: 4,
                    charge: 2,
                    neutral_loss: 0.0,
                    mz: 250.0,
                    relative_intensity: 0.4,
                },
            ],
            ..DdaLibraryEntry::default()
        };
        let entries = generate_decoys(
            vec![target.clone()],
            &LibrarySearchSettings {
                path: "library.parquet".into(),
                ..LibrarySearchSettings::default()
            },
        )
        .unwrap();
        let decoy = entries.iter().find(|entry| entry.is_decoy).unwrap();
        assert_ne!(decoy.proforma, target.proforma);
        assert_eq!(decoy.precursor_neutral_mass, target.precursor_neutral_mass);
        assert_eq!(decoy.proteins, target.proteins);
        let mut intensities = decoy
            .fragments
            .iter()
            .map(|fragment| fragment.relative_intensity)
            .collect::<Vec<_>>();
        intensities.sort_unstable_by(f32::total_cmp);
        assert_eq!(intensities, vec![0.4, 1.0]);
        assert!(decoy.fragments.iter().any(|decoy| {
            target.fragments.iter().any(|target| {
                decoy.kind == target.kind
                    && decoy.ordinal == target.ordinal
                    && decoy.mz != target.mz
            })
        }));
    }

    #[test]
    fn parses_mzspeclib_with_protein_mappings() {
        let text = r#"<mzSpecLib>

<Spectrum=1>
MS:1003061|library spectrum name=PEPTIDER/2
MS:1003208|experimental precursor monoisotopic m/z=478.7378

<Analyte=1>
MS:1003270|proforma peptidoform ion notation=PEPTIDER/2
MS:1001117|theoretical mass=955.461
[1]MS:1000885|protein accession=P1
[2]MS:1000885|protein accession=P2

<Peaks>
300.0 10000.0 y3
400.0 5000.0 b4^2
"#;
        let entries = deserialize_mzspeclib(text).unwrap();
        assert_eq!(entries.len(), 1);
        assert_eq!(entries[0].proteins, "P1;P2");
        assert_eq!(entries[0].fragments.len(), 2);
        assert_eq!(entries[0].fragments[0].relative_intensity, 1.0);
        assert_eq!(entries[0].fragments[1].charge, 2);
    }
}
