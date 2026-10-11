//! TMT quantification
#![allow(clippy::excessive_precision)]
use crate::mass::{Tolerance, PROTON};
use crate::scoring::Feature;
use crate::spectrum::ProcessedSpectrum;
use rayon::prelude::*;
use serde::{Deserialize, Serialize};

#[derive(Debug, Clone, PartialEq, PartialOrd, Serialize, Deserialize, schemars::JsonSchema)]
pub enum Isobaric {
    Tmt6,
    Tmt10,
    Tmt11,
    Tmt16,
    Tmt18,
    User(Vec<f32>),
}

impl Isobaric {
    /// Return the monoisotopic mass of reporter ions
    pub fn reporter_masses(&self) -> &[f32] {
        match self {
            Isobaric::Tmt6 => &TMT6PLEX,
            Isobaric::Tmt10 => &TMT11PLEX[0..10],
            Isobaric::Tmt11 => &TMT11PLEX,
            Isobaric::Tmt16 => &TMT18PLEX[0..16],
            Isobaric::Tmt18 => &TMT18PLEX,
            Isobaric::User(labels) => labels,
        }
    }

    /// Return the monoisotopic mass of tag
    pub fn modification_mass(&self) -> Option<f32> {
        match self {
            Isobaric::Tmt6 | Isobaric::Tmt10 | Isobaric::Tmt11 => Some(229.162932),
            // TMTpro 18plex adds channels to the same TMTpro reagent
            // (Unimod 2016, 304.207146); upstream Sage listed 304.2135.
            Isobaric::Tmt16 | Isobaric::Tmt18 => Some(304.2071),
            Isobaric::User(_) => None,
        }
    }

    /// Return a column name for each tag
    pub fn headers(&self) -> Vec<String> {
        match self {
            Isobaric::User(v) => v
                .iter()
                .enumerate()
                .map(|(idx, _)| format!("user_{}", idx + 1))
                .collect(),
            _ => self
                .reporter_masses()
                .iter()
                .enumerate()
                .map(|(idx, _)| format!("tmt_{}", idx + 1))
                .collect(),
        }
    }
}

#[derive(Debug, Copy, Clone, PartialEq, PartialOrd)]
pub struct Purity {
    pub ratio: f32,
    pub correct_precursors: usize,
    pub incorrect_precursors: usize,
}

#[derive(Debug)]
pub struct Quant<'ms3> {
    /// Top hit for this MS3 spectrum
    pub hit: Feature,
    /// Top chimeric/co-fragmenting hit for this spectrum
    pub chimera: Option<Feature>,
    /// SPS precursor purity for the top hit
    pub hit_purity: Purity,
    /// SPS precursor purity for the chimeric hit
    pub chimera_purity: Option<Purity>,
    /// Quanitified TMT reporter ion intensities
    pub intensities: Vec<Option<f32>>,
    /// MS3 spectrum
    pub spectrum: &'ms3 ProcessedSpectrum,
}

/// Return, for each m/z in `labels`, the most intense peak within the given
/// tolerance window, or `None` when the channel was not observed.
///
/// Only finite, strictly positive intensities count as an observation. A
/// zero-intensity centroid (some converters keep them) is not a detected
/// peak, and a non-finite value can only come from dividing by a zero or
/// missing noise estimate in signal-to-noise mode. Such peaks are skipped, so
/// a real peak in the same window is still used and a window with none of
/// them is reported as missing instead of 0.
///
/// This function is MS-level agnostic, so it can be used for either MS2 or MS3
/// quant.
pub fn find_reporter_ions(
    masses: &[f32],
    intensities: &[f32],
    labels: &[f32],
    label_tolerance: Tolerance,
) -> Vec<Option<f32>> {
    debug_assert_eq!(masses.len(), intensities.len());
    labels
        .iter()
        .map(|&label| {
            let (lo, hi) = label_tolerance.bounds(label);
            let (lo, hi) = (lo - PROTON, hi - PROTON);
            let start = masses.partition_point(|mass| mass.total_cmp(&lo).is_lt());
            let mut best: Option<f32> = None;
            for (&mass, &intensity) in masses[start..].iter().zip(&intensities[start..]) {
                if mass.total_cmp(&hi).is_gt() {
                    break;
                }
                if mass < lo || !intensity.is_finite() || intensity <= 0.0 {
                    continue;
                }
                if best.is_none_or(|max| intensity >= max) {
                    best = Some(intensity);
                }
            }
            best
        })
        .collect()
}

/// Per-channel coverage of isobaric reporter ions across quantified spectra.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct ReporterChannelSummary {
    /// Channel header, e.g. `tmt_126`.
    pub channel: String,
    /// Spectra in which the channel was observed.
    pub observed: usize,
    /// Spectra in which the channel was not observed (written as null).
    pub missing: usize,
    /// Median intensity over the spectra where the channel was observed;
    /// `None` when it was never observed. Missing channels are skipped, never
    /// counted as 0.
    pub median_intensity: Option<f32>,
}

/// Summarise reporter-ion coverage per channel. Missing channels are counted
/// as missing and excluded from the intensity statistics.
pub fn summarize_channels(quant: &[TmtQuant], headers: &[String]) -> Vec<ReporterChannelSummary> {
    headers
        .iter()
        .enumerate()
        .map(|(channel, header)| {
            let mut observed = quant
                .iter()
                .filter_map(|q| q.peaks.get(channel).copied().flatten())
                .collect::<Vec<f32>>();
            let missing = quant.len() - observed.len();
            observed.sort_unstable_by(f32::total_cmp);
            let median_intensity = match observed.len() {
                0 => None,
                n if n % 2 == 1 => Some(observed[n / 2]),
                n => Some((observed[n / 2 - 1] + observed[n / 2]) / 2.0),
            };
            ReporterChannelSummary {
                channel: header.clone(),
                observed: observed.len(),
                missing,
                median_intensity,
            }
        })
        .collect()
}

const TMT6PLEX: [f32; 6] = [
    126.127726, 127.124761, 128.134436, 129.131471, 130.141145, 131.138180,
];

const TMT11PLEX: [f32; 11] = [
    126.127726, 127.124761, 127.131081, 128.128116, 128.134436, 129.131471, 129.137790, 130.134825,
    130.141145, 131.138180, 131.144499,
];

const TMT18PLEX: [f32; 18] = [
    126.127726, 127.124761, 127.131081, 128.128116, 128.134436, 129.131471, 129.137790, 130.134825,
    130.141145, 131.138180, 131.144500, 132.141535, 132.147855, 133.144890, 133.151210, 134.148245,
    134.154565, 135.15160,
];

/// Search MS/MS and quantify isobaric tag intensities from an SPS-MS3 spectrum
///
/// * `scorer`: used for searching/scoring precursor MS2 spectrum
/// * `spectra`: a slice (generally entire mzML) of spectra, that can be searched for precursor spectra
/// * `ms3`: The MS3 spectrum to search and quantify
/// * `isobaric_labels`: specify label m/zs to be used
/// * `isobaric_tolerance`: specify label tolerance
// pub fn quantify_sps<'a, 'b>(
//     scorer: &'a Scorer<'a>,
//     spectra: &[ProcessedSpectrum],
//     ms3: &'b ProcessedSpectrum,
//     isobaric_labels: &Isobaric,
//     isobaric_tolerance: Tolerance,
// ) -> Option<Quant<'b>> {
//     let first_precursor = ms3
//         .precursors
//         .first()
//         .expect("MS3 scan without at least one precursor!");

//     let ms2 = spectrum::find_spectrum_by_id(
//         spectra,
//         first_precursor
//             .scan
//             .expect("MS3 scan without a MS2 precursor scan ID"),
//     )
//     .expect("Couldn't locate parent MS2 scan!");

//     let ms1_charge = ms2
//         .precursors
//         .get(0)
//         .and_then(|p| p.charge)
//         .unwrap_or(2)
//         .saturating_sub(1);

//     let scores = scorer.score_chimera(ms2);
//     let hit = scores.first()?.clone();
//     let peptide = &scorer.db[hit.peptide_idx];
//     let hit_purity = purity_of_match(
//         &ms3.precursors,
//         ms2,
//         &mk_theoretical(peptide),
//         ms1_charge,
//         scorer.fragment_tol,
//     );

//     let chimera = scores.get(1).cloned();
//     let chimera_purity = chimera.as_ref().map(|score| {
//         purity_of_match(
//             &ms3.precursors,
//             ms2,
//             &mk_theoretical(&scorer.db[score.peptide_idx]),
//             ms1_charge,
//             scorer.fragment_tol,
//         )
//     });

//     Some(Quant {
//         hit,
//         hit_purity,
//         chimera,
//         chimera_purity,
//         intensities: find_reporter_ions(
//             &ms3.peaks,
//             isobaric_labels.reporter_masses(),
//             isobaric_tolerance,
//         ),
//         spectrum: ms3,
//     })
// }

#[derive(Clone)]
pub struct TmtQuant {
    pub spec_id: String,
    pub file_id: usize,
    /// Zero-based occurrence of `(file_id, spec_id)` among the MSn spectra, so
    /// repeated spectrum IDs (e.g. duplicate MGF titles) keep their own
    /// reporter ions. MS3 reporter spectra refer to their MS2 scan and use 0.
    pub occurrence: usize,
    pub ion_injection_time: f32,
    /// Reporter intensity per channel, in `Isobaric::reporter_masses` order.
    /// `None` means the channel was not observed; it is never filled with 0.
    pub peaks: Vec<Option<f32>>,
}

/// Quantify isobaric tags from an MS2 or MS3 spectrum
///
/// * `spectra`: a slice (generally entire mzML) of spectra, that can be searched
///   for precursor spectra
/// * `isobaric_labels`: specify label m/zs to be used
/// * `isobaric_tolerance`: specify label tolerance
/// * `level`: MSn level to extract isobaric peaks from
pub fn quantify(
    spectra: &[ProcessedSpectrum],
    isobaric_labels: &Isobaric,
    isobaric_tolerance: Tolerance,
    level: u8,
) -> Vec<TmtQuant> {
    let mut seen = std::collections::HashMap::with_capacity(spectra.len());
    let occurrences = spectra
        .iter()
        .map(|spectrum| {
            let count = seen
                .entry((spectrum.file_id, spectrum.id.as_str()))
                .or_insert(0usize);
            *count += 1;
            *count - 1
        })
        .collect::<Vec<_>>();
    spectra
        .par_iter()
        .zip(occurrences.par_iter())
        .filter(|(spectrum, _)| spectrum.level == level)
        .filter_map(|(spectrum, &occurrence)| {
            let (spec_id, occurrence) = match level {
                1 => return None,
                2 => (spectrum.id.clone(), occurrence),
                _ => (
                    spectrum
                        .precursors
                        .first()
                        .and_then(|precursor| precursor.spectrum_ref.clone())
                        .unwrap_or_default(),
                    0,
                ),
            };

            let peaks = find_reporter_ions(
                &spectrum.masses,
                &spectrum.intensities,
                isobaric_labels.reporter_masses(),
                isobaric_tolerance,
            );

            Some(TmtQuant {
                spec_id,
                file_id: spectrum.file_id,
                occurrence,
                ion_injection_time: spectrum.ion_injection_time,
                peaks,
            })
        })
        .collect()
}

#[cfg(test)]
#[path = "../tests/unit/tmt.rs"]
mod tests;
