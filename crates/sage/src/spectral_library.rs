//! Build deterministic empirical spectral-library entries from confident PSMs.

use crate::database::{IndexedDatabase, PeptideIx};
use crate::ion_series::Kind;
use crate::mass::PROTON;
use crate::peptide::Peptide;
use crate::scoring::Feature;
use serde::{Deserialize, Serialize};
use std::collections::{HashMap, HashSet};
use std::fmt::Write;

#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum SpectralLibraryFormat {
    SageParquet,
    #[serde(rename = "mzspeclib")]
    MzSpecLib,
}

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum SpectralLibraryStrategy {
    #[default]
    BestPsm,
    Consensus,
}

fn default_psm_q_value() -> f32 {
    0.01
}

fn default_peptide_q_value() -> f32 {
    0.01
}

fn default_min_matched_peaks() -> u32 {
    6
}

fn default_max_fragments() -> usize {
    20
}

fn default_min_relative_intensity() -> f32 {
    0.01
}

fn default_min_consensus_psms() -> usize {
    1
}

fn default_min_fragment_frequency() -> f32 {
    0.5
}

fn default_formats() -> Vec<SpectralLibraryFormat> {
    vec![
        SpectralLibraryFormat::SageParquet,
        SpectralLibraryFormat::MzSpecLib,
    ]
}

/// Settings for an empirical library built from confidently identified PSMs.
#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(default, deny_unknown_fields)]
pub struct SpectralLibrarySettings {
    pub enabled: bool,
    pub psm_q_value: f32,
    pub peptide_q_value: f32,
    pub strategy: SpectralLibraryStrategy,
    pub min_matched_peaks: u32,
    pub max_fragments: usize,
    pub min_relative_intensity: f32,
    pub min_consensus_psms: usize,
    pub min_fragment_frequency: f32,
    pub include_chimeric: bool,
    pub formats: Vec<SpectralLibraryFormat>,
}

impl Default for SpectralLibrarySettings {
    fn default() -> Self {
        Self {
            enabled: false,
            psm_q_value: default_psm_q_value(),
            peptide_q_value: default_peptide_q_value(),
            strategy: SpectralLibraryStrategy::BestPsm,
            min_matched_peaks: default_min_matched_peaks(),
            max_fragments: default_max_fragments(),
            min_relative_intensity: default_min_relative_intensity(),
            min_consensus_psms: default_min_consensus_psms(),
            min_fragment_frequency: default_min_fragment_frequency(),
            include_chimeric: false,
            formats: default_formats(),
        }
    }
}

impl SpectralLibrarySettings {
    pub fn validate(&self) -> Result<(), String> {
        if !self.psm_q_value.is_finite() || !(0.0..=1.0).contains(&self.psm_q_value) {
            return Err("spectral_library.psm_q_value must be between 0 and 1".into());
        }
        if !self.peptide_q_value.is_finite() || !(0.0..=1.0).contains(&self.peptide_q_value) {
            return Err("spectral_library.peptide_q_value must be between 0 and 1".into());
        }
        if self.min_matched_peaks == 0 {
            return Err("spectral_library.min_matched_peaks must be greater than zero".into());
        }
        if self.max_fragments == 0 {
            return Err("spectral_library.max_fragments must be greater than zero".into());
        }
        if !self.min_relative_intensity.is_finite()
            || !(0.0..=1.0).contains(&self.min_relative_intensity)
        {
            return Err("spectral_library.min_relative_intensity must be between 0 and 1".into());
        }
        if self.min_consensus_psms == 0 {
            return Err("spectral_library.min_consensus_psms must be greater than zero".into());
        }
        if !self.min_fragment_frequency.is_finite()
            || !(0.0..=1.0).contains(&self.min_fragment_frequency)
            || self.min_fragment_frequency == 0.0
        {
            return Err(
                "spectral_library.min_fragment_frequency must be greater than zero and at most one"
                    .into(),
            );
        }
        if self.enabled && self.formats.is_empty() {
            return Err("spectral_library.formats must not be empty when enabled".into());
        }
        let unique = self.formats.iter().copied().collect::<HashSet<_>>();
        if unique.len() != self.formats.len() {
            return Err("spectral_library.formats must not contain duplicates".into());
        }
        Ok(())
    }

    pub fn writes(&self, format: SpectralLibraryFormat) -> bool {
        self.enabled && self.formats.contains(&format)
    }
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct LibrarySelection {
    pub feature_index: usize,
    pub supporting_psms: usize,
    pub feature_indices: Vec<usize>,
}

#[derive(Clone, Debug, PartialEq)]
pub struct LibraryFragment {
    pub kind: Kind,
    pub ordinal: i32,
    pub charge: i32,
    pub neutral_loss: f32,
    pub mz: f32,
    pub relative_intensity: f32,
}

impl LibraryFragment {
    pub fn annotation(&self) -> String {
        let kind = match self.kind {
            Kind::A => 'a',
            Kind::B => 'b',
            Kind::C => 'c',
            Kind::X => 'x',
            Kind::Y => 'y',
            Kind::Z => 'z',
        };
        let mut annotation = format!("{kind}{}", self.ordinal);
        if self.neutral_loss > 0.0 {
            if (self.neutral_loss - 18.010_565).abs() <= 0.01 {
                annotation.push_str("-H2O");
            } else if (self.neutral_loss - 17.026_55).abs() <= 0.01 {
                annotation.push_str("-NH3");
            } else if (self.neutral_loss - 97.976_9).abs() <= 0.01 {
                annotation.push_str("-H3PO4");
            } else {
                // Sage retains arbitrary loss masses. mzPAF requires a chemically
                // meaningful loss annotation, so do not invent one here.
                return "?".into();
            }
        }
        if self.charge > 1 {
            write!(annotation, "^{}", self.charge).expect("write to String");
        }
        annotation
    }
}

#[derive(Clone, Debug, PartialEq)]
pub struct SpectralLibraryEntry {
    pub library_entry_id: String,
    pub source_psm_id: usize,
    pub source_file: String,
    pub source_spectrum: String,
    pub modified_peptide: String,
    pub proforma: String,
    pub stripped_peptide: String,
    pub proteins: String,
    pub label_channel: Option<String>,
    pub label_group: Option<String>,
    pub label_reference: Option<String>,
    pub precursor_charge: u8,
    pub precursor_neutral_mass: f32,
    pub precursor_mz: f32,
    pub retention_time_minutes: f32,
    pub aligned_retention_time_minutes: f32,
    pub ion_mobility: f32,
    pub spectrum_q: f32,
    pub peptide_q: f32,
    pub supporting_psms: usize,
    pub fragments: Vec<LibraryFragment>,
}

fn eligible(feature: &Feature, settings: &SpectralLibrarySettings) -> bool {
    feature.label == 1
        && feature.spectrum_q <= settings.psm_q_value
        && feature.peptide_q <= settings.peptide_q_value
        && feature.matched_peaks >= settings.min_matched_peaks
        && (settings.include_chimeric || feature.rank == 1)
}

fn better_candidate(candidate: &Feature, incumbent: &Feature) -> bool {
    candidate
        .spectrum_q
        .total_cmp(&incumbent.spectrum_q)
        .then_with(|| candidate.peptide_q.total_cmp(&incumbent.peptide_q))
        .then_with(|| {
            incumbent
                .discriminant_score
                .total_cmp(&candidate.discriminant_score)
        })
        .then_with(|| incumbent.hyperscore.total_cmp(&candidate.hyperscore))
        .then_with(|| candidate.psm_id.cmp(&incumbent.psm_id))
        .is_lt()
}

/// Select deterministic supporting PSMs for each exact peptidoform and charge.
pub fn select_psms(
    features: &[Feature],
    database: &IndexedDatabase,
    settings: &SpectralLibrarySettings,
) -> Vec<LibrarySelection> {
    if !settings.enabled {
        return Vec::new();
    }

    let mut groups: HashMap<(PeptideIx, u8), Vec<usize>> = HashMap::new();
    for (feature_index, feature) in features.iter().enumerate() {
        if !eligible(feature, settings) {
            continue;
        }
        groups
            .entry((feature.peptide_idx, feature.charge))
            .or_default()
            .push(feature_index);
    }

    let mut selected = groups
        .into_values()
        .filter(|feature_indices| {
            settings.strategy == SpectralLibraryStrategy::BestPsm
                || feature_indices.len() >= settings.min_consensus_psms
        })
        .map(|feature_indices| {
            let feature_index = feature_indices
                .iter()
                .copied()
                .reduce(|best, candidate| {
                    if better_candidate(&features[candidate], &features[best]) {
                        candidate
                    } else {
                        best
                    }
                })
                .expect("library PSM group is not empty");
            let supporting_psms = feature_indices.len();
            let feature_indices = match settings.strategy {
                SpectralLibraryStrategy::BestPsm => vec![feature_index],
                SpectralLibraryStrategy::Consensus => feature_indices,
            };
            LibrarySelection {
                feature_index,
                supporting_psms,
                feature_indices,
            }
        })
        .collect::<Vec<_>>();
    selected.sort_unstable_by(|a, b| {
        let left = &features[a.feature_index];
        let right = &features[b.feature_index];
        database[left.peptide_idx]
            .to_string()
            .cmp(&database[right.peptide_idx].to_string())
            .then_with(|| left.charge.cmp(&right.charge))
            .then_with(|| left.psm_id.cmp(&right.psm_id))
    });
    selected
}

/// Compatibility wrapper for callers that previously selected best PSMs directly.
pub fn select_best_psms(
    features: &[Feature],
    database: &IndexedDatabase,
    settings: &SpectralLibrarySettings,
) -> Vec<LibrarySelection> {
    select_psms(features, database, settings)
}

/// Render an unambiguous mass-delta ProForma peptidoform.
pub fn mass_delta_proforma(peptide: &Peptide) -> String {
    let mut output = String::new();
    if let Some(mass) = peptide.nterm {
        write!(output, "[{mass:+}]-").expect("write to String");
    }
    for (index, residue) in peptide.sequence.iter().copied().enumerate() {
        output.push(residue as char);
        let mass = peptide.modification_at(index);
        if mass != 0.0 {
            write!(output, "[{mass:+}]").expect("write to String");
        }
    }
    if let Some(mass) = peptide.cterm {
        write!(output, "-[{mass:+}]").expect("write to String");
    }
    output
}

fn median(mut values: Vec<f32>) -> Option<f32> {
    values.retain(|value| value.is_finite());
    if values.is_empty() {
        return None;
    }
    values.sort_unstable_by(f32::total_cmp);
    let middle = values.len() / 2;
    if values.len().is_multiple_of(2) {
        Some((values[middle - 1] + values[middle]) / 2.0)
    } else {
        Some(values[middle])
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
struct FragmentKey {
    kind: Kind,
    ordinal: i32,
    charge: i32,
    neutral_loss_bits: u32,
}

fn build_consensus_fragments(
    features: &[&Feature],
    settings: &SpectralLibrarySettings,
) -> Result<Vec<LibraryFragment>, String> {
    let mut fragments_by_key = HashMap::<FragmentKey, (f32, Vec<f32>)>::new();
    for feature in features {
        let fragments = feature.fragments.as_ref().ok_or_else(|| {
            format!(
                "PSM {} was selected for the spectral library but was not annotated",
                feature.psm_id
            )
        })?;
        let lengths = [
            fragments.kinds.len(),
            fragments.fragment_ordinals.len(),
            fragments.charges.len(),
            fragments.neutral_losses.len(),
            fragments.mz_calculated.len(),
            fragments.intensities.len(),
        ];
        if lengths.iter().any(|length| *length != lengths[0]) {
            return Err(format!(
                "PSM {} has inconsistent fragment annotation lengths",
                feature.psm_id
            ));
        }
        let max_intensity = fragments
            .intensities
            .iter()
            .copied()
            .filter(|value| value.is_finite() && *value > 0.0)
            .max_by(f32::total_cmp)
            .ok_or_else(|| format!("PSM {} has no positive fragment intensity", feature.psm_id))?;
        let mut spectrum_fragments = HashMap::<FragmentKey, (f32, f32)>::new();
        for index in 0..lengths[0] {
            let intensity = fragments.intensities[index];
            if !intensity.is_finite() || intensity <= 0.0 {
                continue;
            }
            let key = FragmentKey {
                kind: fragments.kinds[index],
                ordinal: fragments.fragment_ordinals[index],
                charge: fragments.charges[index],
                neutral_loss_bits: fragments.neutral_losses[index].to_bits(),
            };
            let candidate = (fragments.mz_calculated[index], intensity / max_intensity);
            spectrum_fragments
                .entry(key)
                .and_modify(|current| {
                    if candidate.1 > current.1 {
                        *current = candidate;
                    }
                })
                .or_insert(candidate);
        }
        for (key, (mz, intensity)) in spectrum_fragments {
            fragments_by_key
                .entry(key)
                .or_insert_with(|| (mz, Vec::new()))
                .1
                .push(intensity);
        }
    }

    let support = features.len() as f32;
    let mut fragments = fragments_by_key
        .into_iter()
        .filter_map(|(key, (mz, intensities))| {
            let frequency = intensities.len() as f32 / support;
            (frequency >= settings.min_fragment_frequency).then(|| LibraryFragment {
                kind: key.kind,
                ordinal: key.ordinal,
                charge: key.charge,
                neutral_loss: f32::from_bits(key.neutral_loss_bits),
                mz,
                relative_intensity: median(intensities).unwrap_or_default(),
            })
        })
        .collect::<Vec<_>>();
    let max_intensity = fragments
        .iter()
        .map(|fragment| fragment.relative_intensity)
        .max_by(f32::total_cmp)
        .unwrap_or_default();
    if max_intensity <= 0.0 {
        return Ok(Vec::new());
    }
    fragments.iter_mut().for_each(|fragment| {
        fragment.relative_intensity /= max_intensity;
    });
    fragments.retain(|fragment| fragment.relative_intensity >= settings.min_relative_intensity);
    fragments.sort_unstable_by(|a, b| {
        b.relative_intensity
            .total_cmp(&a.relative_intensity)
            .then_with(|| a.mz.total_cmp(&b.mz))
            .then_with(|| a.ordinal.cmp(&b.ordinal))
            .then_with(|| a.charge.cmp(&b.charge))
    });
    fragments.truncate(settings.max_fragments);
    fragments.sort_unstable_by(|a, b| a.mz.total_cmp(&b.mz));
    Ok(fragments)
}

pub fn build_entries(
    features: &[Feature],
    database: &IndexedDatabase,
    filenames: &[String],
    selections: &[LibrarySelection],
    settings: &SpectralLibrarySettings,
) -> Result<Vec<SpectralLibraryEntry>, String> {
    let mut entries = Vec::with_capacity(selections.len());
    for selection in selections {
        let feature = features.get(selection.feature_index).ok_or_else(|| {
            format!(
                "spectral-library selection references missing feature {}",
                selection.feature_index
            )
        })?;
        let source_file = filenames
            .get(feature.file_id)
            .ok_or_else(|| format!("missing source filename for file {}", feature.file_id))?;
        let peptide = &database[feature.peptide_idx];
        let supporting_features = selection
            .feature_indices
            .iter()
            .map(|index| {
                features.get(*index).ok_or_else(|| {
                    format!("spectral-library selection references missing feature {index}")
                })
            })
            .collect::<Result<Vec<_>, _>>()?;
        let selected_fragments = build_consensus_fragments(&supporting_features, settings)?;
        if selected_fragments.is_empty() {
            continue;
        }

        let proforma = mass_delta_proforma(peptide);
        let retention_time_minutes = median(
            supporting_features
                .iter()
                .map(|feature| feature.rt)
                .filter(|rt| *rt > 0.0)
                .collect(),
        )
        .unwrap_or(feature.rt);
        let aligned_retention_time_minutes = median(
            supporting_features
                .iter()
                .map(|feature| {
                    if feature.aligned_rt.is_finite() && feature.aligned_rt > 0.0 {
                        feature.aligned_rt
                    } else {
                        feature.rt
                    }
                })
                .filter(|rt| *rt > 0.0)
                .collect(),
        )
        .unwrap_or(retention_time_minutes);
        let ion_mobility = median(
            supporting_features
                .iter()
                .map(|feature| feature.ims)
                .filter(|mobility| *mobility > 0.0)
                .collect(),
        )
        .unwrap_or_default();
        entries.push(SpectralLibraryEntry {
            library_entry_id: format!("{proforma}/{}", feature.charge),
            source_psm_id: feature.psm_id,
            source_file: source_file.clone(),
            source_spectrum: feature.spec_id.clone(),
            modified_peptide: peptide.to_string(),
            proforma,
            stripped_peptide: String::from_utf8_lossy(&peptide.sequence).into_owned(),
            proteins: peptide.proteins(&database.decoy_tag, database.generate_decoys),
            label_channel: peptide.label_channel.as_deref().map(str::to_owned),
            label_group: peptide
                .label_channel
                .as_ref()
                .map(|_| peptide.label_group()),
            label_reference: peptide
                .label_channel
                .as_ref()
                .and(database.label_reference.as_deref())
                .map(str::to_owned),
            precursor_charge: feature.charge,
            precursor_neutral_mass: feature.calcmass,
            precursor_mz: feature.calcmass / feature.charge as f32 + PROTON,
            retention_time_minutes,
            aligned_retention_time_minutes,
            ion_mobility,
            spectrum_q: feature.spectrum_q,
            peptide_q: feature.peptide_q,
            supporting_psms: selection.supporting_psms,
            fragments: selected_fragments,
        });
    }
    Ok(entries)
}

/// Serialize empirical entries using mzSpecLib 1.0 text syntax.
pub fn serialize_mzspeclib(
    entries: &[SpectralLibraryEntry],
    sage_version: &str,
    strategy: SpectralLibraryStrategy,
) -> Vec<u8> {
    let mut output = String::new();
    let (description, strategy_name) = match strategy {
        SpectralLibraryStrategy::BestPsm => ("Best-PSM", "best_psm"),
        SpectralLibraryStrategy::Consensus => ("Consensus", "consensus"),
    };
    writeln!(output, "<mzSpecLib>").unwrap();
    writeln!(output, "MS:1003186|library format version=1.0").unwrap();
    writeln!(output, "MS:1003190|library version={sage_version}").unwrap();
    writeln!(
        output,
        "MS:1003187|library identifier=sage-plus:{sage_version}"
    )
    .unwrap();
    writeln!(
        output,
        "MS:1003188|library name=Sage Plus empirical spectral library"
    )
    .unwrap();
    writeln!(
        output,
        "MS:1003189|library description={description} empirical library generated by Sage Plus"
    )
    .unwrap();
    writeln!(
        output,
        "MS:1003206|library creation log=Sage Plus {sage_version}, strategy={strategy_name}"
    )
    .unwrap();

    for (spectrum_number, entry) in entries.iter().enumerate() {
        writeln!(output).unwrap();
        writeln!(output, "<Spectrum={}>", spectrum_number + 1).unwrap();
        writeln!(
            output,
            "MS:1003061|library spectrum name={}/{}",
            entry.proforma, entry.precursor_charge
        )
        .unwrap();
        writeln!(output, "MS:1000511|ms level=2").unwrap();
        match strategy {
            SpectralLibraryStrategy::BestPsm => writeln!(
                output,
                "MS:1003065|spectrum aggregation type=MS:1003066|singleton spectrum"
            ),
            SpectralLibraryStrategy::Consensus => writeln!(
                output,
                "MS:1003065|spectrum aggregation type=MS:1003067|consensus spectrum"
            ),
        }
        .unwrap();
        writeln!(
            output,
            "MS:1003072|spectrum origin type=MS:1003073|observed spectrum"
        )
        .unwrap();
        writeln!(
            output,
            "MS:1003208|experimental precursor monoisotopic m/z={:.6}",
            entry.precursor_mz
        )
        .unwrap();
        writeln!(
            output,
            "[1]MS:1000894|retention time={:.6}",
            entry.aligned_retention_time_minutes
        )
        .unwrap();
        writeln!(output, "[1]UO:0000000|unit=UO:0000031|minute").unwrap();
        writeln!(
            output,
            "MS:1003059|number of peaks={}",
            entry.fragments.len()
        )
        .unwrap();

        writeln!(output).unwrap();
        writeln!(output, "<Analyte=1>").unwrap();
        writeln!(
            output,
            "MS:1003270|proforma peptidoform ion notation={}/{}",
            entry.proforma, entry.precursor_charge
        )
        .unwrap();
        writeln!(output, "MS:1000041|charge state={}", entry.precursor_charge).unwrap();
        writeln!(
            output,
            "MS:1003043|number of residues={}",
            entry.stripped_peptide.len()
        )
        .unwrap();
        writeln!(
            output,
            "MS:1003053|theoretical monoisotopic m/z={:.6}",
            entry.precursor_mz
        )
        .unwrap();
        writeln!(
            output,
            "MS:1001117|theoretical mass={:.6}",
            entry.precursor_neutral_mass
        )
        .unwrap();
        writeln!(
            output,
            "MS:1003243|adduct ion mass={:.6}",
            entry.precursor_mz * f32::from(entry.precursor_charge)
        )
        .unwrap();
        for (index, protein) in entry
            .proteins
            .split(';')
            .filter(|protein| !protein.is_empty())
            .enumerate()
        {
            writeln!(
                output,
                "[{}]MS:1000885|protein accession={protein}",
                index + 1
            )
            .unwrap();
        }
        if let Some(channel) = &entry.label_channel {
            writeln!(output, "SAGE:1000001|label channel={channel}").unwrap();
        }
        if let Some(group) = &entry.label_group {
            writeln!(output, "SAGE:1000002|label group={group}").unwrap();
        }
        if let Some(reference) = &entry.label_reference {
            writeln!(output, "SAGE:1000003|label reference={reference}").unwrap();
        }

        writeln!(output).unwrap();
        writeln!(output, "<Interpretation=1>").unwrap();
        writeln!(
            output,
            "MS:1002354|PSM-level q-value={:.8}",
            entry.spectrum_q
        )
        .unwrap();

        writeln!(output).unwrap();
        writeln!(output, "<Peaks>").unwrap();
        for fragment in &entry.fragments {
            writeln!(
                output,
                "{:.6}\t{:.6}\t{}",
                fragment.mz,
                fragment.relative_intensity * 10_000.0,
                fragment.annotation()
            )
            .unwrap();
        }
    }
    output.into_bytes()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::database::PeptideIx;
    use crate::scoring::Fragments;
    use std::sync::Arc;

    fn feature(peptide_idx: u32, psm_id: usize, spectrum_q: f32, score: f32) -> Feature {
        Feature {
            peptide_idx: PeptideIx(peptide_idx),
            psm_id,
            label: 1,
            rank: 1,
            charge: 2,
            matched_peaks: 8,
            spectrum_q,
            peptide_q: 0.005,
            discriminant_score: score,
            calcmass: 1_000.0,
            rt: 12.0,
            aligned_rt: 11.5,
            fragments: Some(Fragments {
                kinds: vec![Kind::B, Kind::Y, Kind::Y],
                fragment_ordinals: vec![2, 3, 4],
                charges: vec![1, 1, 2],
                neutral_losses: vec![0.0, 18.010_565, 42.0],
                mz_calculated: vec![200.0, 300.0, 400.0],
                mz_experimental: vec![200.01, 300.01, 400.01],
                intensities: vec![25.0, 100.0, 0.5],
            }),
            ..Feature::default()
        }
    }

    fn database() -> IndexedDatabase {
        let mut database = IndexedDatabase::default();
        database.peptides.push(Peptide {
            sequence: Arc::from(&b"PEPTIDE"[..]),
            proteins: vec![Arc::from("P12345")],
            ..Peptide::default()
        });
        database
    }

    #[test]
    fn best_psm_is_deterministic_and_counts_support() {
        let database = database();
        let settings = SpectralLibrarySettings {
            enabled: true,
            ..Default::default()
        };
        let features = vec![
            feature(0, 3, 0.005, 10.0),
            feature(0, 2, 0.001, 5.0),
            feature(0, 1, 0.001, 8.0),
        ];
        let selected = select_best_psms(&features, &database, &settings);
        assert_eq!(selected.len(), 1);
        assert_eq!(selected[0].feature_index, 2);
        assert_eq!(selected[0].supporting_psms, 3);
        assert_eq!(selected[0].feature_indices, vec![2]);
    }

    #[test]
    fn builds_normalized_top_fragments_and_mzspeclib() {
        let database = database();
        let settings = SpectralLibrarySettings {
            enabled: true,
            max_fragments: 2,
            ..Default::default()
        };
        let features = vec![feature(0, 1, 0.001, 8.0)];
        let selected = select_best_psms(&features, &database, &settings);
        let mut entries = build_entries(
            &features,
            &database,
            &["sample.mzML".into()],
            &selected,
            &settings,
        )
        .unwrap();
        assert_eq!(entries.len(), 1);
        assert_eq!(entries[0].fragments.len(), 2);
        assert_eq!(entries[0].fragments[0].relative_intensity, 0.25);
        assert_eq!(entries[0].fragments[1].relative_intensity, 1.0);
        assert_eq!(entries[0].fragments[1].annotation(), "y3-H2O");
        entries[0].label_channel = Some("heavy".into());
        entries[0].label_group = Some("PEPTIDE".into());
        entries[0].label_reference = Some("light".into());

        let text = String::from_utf8(serialize_mzspeclib(
            &entries,
            "0.16.0-beta.1",
            SpectralLibraryStrategy::BestPsm,
        ))
        .unwrap();
        assert!(text.starts_with("<mzSpecLib>\n"));
        assert!(text.contains("<Spectrum=1>"));
        assert!(text.contains("MS:1003270|proforma peptidoform ion notation=PEPTIDE/2"));
        assert!(text.contains("300.000000\t10000.000000\ty3-H2O"));
        assert!(text.contains("SAGE:1000001|label channel=heavy"));
        let parsed = crate::spectral_library_search::deserialize_mzspeclib(&text).unwrap();
        assert_eq!(parsed[0].label_channel.as_deref(), Some("heavy"));
        assert_eq!(parsed[0].label_group.as_deref(), Some("PEPTIDE"));
        assert_eq!(parsed[0].label_reference.as_deref(), Some("light"));
    }

    #[test]
    fn consensus_uses_median_properties_and_reproducible_fragments() {
        let database = database();
        let settings = SpectralLibrarySettings {
            enabled: true,
            strategy: SpectralLibraryStrategy::Consensus,
            min_fragment_frequency: 0.66,
            ..Default::default()
        };
        let mut features = vec![
            feature(0, 1, 0.001, 8.0),
            feature(0, 2, 0.002, 7.0),
            feature(0, 3, 0.003, 6.0),
        ];
        features[0].rt = 10.0;
        features[0].aligned_rt = 9.0;
        features[0].ims = 1.0;
        features[0].fragments.as_mut().unwrap().intensities = vec![25.0, 100.0, 0.0];
        features[1].rt = 12.0;
        features[1].aligned_rt = 11.0;
        features[1].ims = 1.2;
        features[1].fragments.as_mut().unwrap().intensities = vec![50.0, 100.0, 0.0];
        features[2].rt = 30.0;
        features[2].aligned_rt = 25.0;
        features[2].ims = 1.4;
        features[2].fragments.as_mut().unwrap().intensities = vec![75.0, 0.0, 100.0];

        let selected = select_psms(&features, &database, &settings);
        assert_eq!(selected[0].feature_indices, vec![0, 1, 2]);
        let entries = build_entries(
            &features,
            &database,
            &["one.mzML".into()],
            &selected,
            &settings,
        )
        .unwrap();
        assert_eq!(entries[0].supporting_psms, 3);
        assert_eq!(entries[0].retention_time_minutes, 12.0);
        assert_eq!(entries[0].aligned_retention_time_minutes, 11.0);
        assert_eq!(entries[0].ion_mobility, 1.2);
        assert_eq!(entries[0].fragments.len(), 2);
        assert_eq!(entries[0].fragments[0].relative_intensity, 0.5);
        assert_eq!(entries[0].fragments[1].relative_intensity, 1.0);

        let text = String::from_utf8(serialize_mzspeclib(
            &entries,
            "0.16.0-beta.1",
            SpectralLibraryStrategy::Consensus,
        ))
        .unwrap();
        assert!(text.contains("MS:1003067|consensus spectrum"));
        assert!(text.contains("strategy=consensus"));
    }

    #[test]
    fn settings_reject_invalid_values_and_duplicate_formats() {
        let mut settings = SpectralLibrarySettings {
            enabled: true,
            ..Default::default()
        };
        settings.max_fragments = 0;
        assert!(settings.validate().is_err());
        settings.max_fragments = 10;
        settings.min_fragment_frequency = 0.0;
        assert!(settings.validate().is_err());
        settings.min_fragment_frequency = 0.5;
        settings.formats.push(SpectralLibraryFormat::MzSpecLib);
        assert!(settings.validate().is_err());
    }

    #[test]
    fn mass_delta_proforma_does_not_depend_on_modification_names() {
        let mut peptide = Peptide {
            sequence: Arc::from(&b"ACD"[..]),
            modifications: vec![0.0, 57.021_465, 0.0],
            nterm: Some(42.010_565),
            ..Peptide::default()
        };
        peptide.applied_modifications = Arc::new(vec![]);
        assert_eq!(
            mass_delta_proforma(&peptide),
            "[+42.010567]-AC[+57.021465]D"
        );
    }
}
