//! Experimental DDA spectrum-to-library candidate indexing and scoring.
//!
//! This module deliberately stops short of assigning q-values. A production
//! library-search mode needs an explicit decoy-library policy before its scores
//! can be connected to Sage's target-decoy FDR pipeline.

use crate::database::binary_search_slice;
use crate::mass::{Tolerance, PROTON};
use crate::spectral_library::{LibraryFragment, SpectralLibraryEntry};
use crate::spectrum::ProcessedSpectrum;

#[derive(Clone, Debug, PartialEq)]
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
        for entry in &entries {
            if entry.library_entry_id.is_empty() {
                return Err("library entry identifiers must not be empty".into());
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
}
