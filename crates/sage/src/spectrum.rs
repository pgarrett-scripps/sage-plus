use crate::mass::{Tolerance, NEUTRON, PROTON};
use serde::{Deserialize, Serialize};
use std::borrow::Cow;

#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
#[serde(default, deny_unknown_fields)]
pub struct DeisotopeSettings {
    pub enabled: bool,
    pub ppm_tolerance: f32,
    pub max_charge: Option<u8>,
    pub min_envelope_peaks: usize,
    pub max_envelope_peaks: usize,
    pub min_score: f32,
    pub max_isotope_log2_ratio: f32,
}

impl DeisotopeSettings {
    pub fn from_enabled(enabled: bool) -> Self {
        Self {
            enabled,
            ..Self::default()
        }
    }

    pub fn validate(&self) -> Result<(), String> {
        if !self.ppm_tolerance.is_finite() || self.ppm_tolerance <= 0.0 {
            return Err("deisotope.ppm_tolerance must be finite and greater than zero".into());
        }
        if self.max_charge == Some(0) {
            return Err("deisotope.max_charge must be greater than zero".into());
        }
        if !(2..=4).contains(&self.min_envelope_peaks) {
            return Err("deisotope.min_envelope_peaks must be between 2 and 4".into());
        }
        if !(2..=4).contains(&self.max_envelope_peaks) {
            return Err("deisotope.max_envelope_peaks must be between 2 and 4".into());
        }
        if self.min_envelope_peaks > self.max_envelope_peaks {
            return Err("deisotope.min_envelope_peaks cannot exceed max_envelope_peaks".into());
        }
        if !self.min_score.is_finite() || !(0.0..=1.0).contains(&self.min_score) {
            return Err("deisotope.min_score must be between 0 and 1".into());
        }
        if !self.max_isotope_log2_ratio.is_finite() || self.max_isotope_log2_ratio < 0.0 {
            return Err("deisotope.max_isotope_log2_ratio must be finite and nonnegative".into());
        }
        Ok(())
    }
}

impl Default for DeisotopeSettings {
    fn default() -> Self {
        Self {
            enabled: true,
            ppm_tolerance: 10.0,
            max_charge: None,
            min_envelope_peaks: 2,
            max_envelope_peaks: 4,
            min_score: 0.45,
            max_isotope_log2_ratio: 1.5,
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Deserialize)]
#[serde(untagged)]
pub enum DeisotopeConfig {
    Enabled(bool),
    Settings(DeisotopeSettings),
}

impl DeisotopeConfig {
    pub fn resolve(self) -> DeisotopeSettings {
        match self {
            Self::Enabled(enabled) => DeisotopeSettings::from_enabled(enabled),
            Self::Settings(settings) => settings,
        }
    }
}

/// A de-isotoped peak, that might have some charge state information
#[derive(PartialEq, PartialOrd, Debug, Copy, Clone)]
pub struct Deisotoped {
    pub mz: f32,
    // Cumulative intensity of all isotopic peaks in the envelope higher than this one
    pub intensity: f32,
    // Assigned charge
    pub charge: Option<u8>,
    // If `Some(idx)`, idx is the index of the parent isotopic envelope
    pub envelope: Option<usize>,
}

#[derive(Debug, Clone)]
pub struct SpectrumProcessor {
    pub take_top_n: usize,
    pub min_deisotope_mz: f32,
    pub deisotope: DeisotopeSettings,
}

#[derive(Default, Debug, Clone)]
pub struct Precursor {
    pub mz: f32,
    pub intensity: Option<f32>,
    pub charge: Option<u8>,
    // pub scan: Option<usize>,
    pub spectrum_ref: Option<String>,
    pub isolation_window: Option<Tolerance>,
    pub inverse_ion_mobility: Option<f32>,
}

#[derive(Clone, Default, Debug)]
pub struct ProcessedSpectrum {
    /// MSn level
    pub level: u8,
    /// Scan ID
    pub id: String,
    /// File ID
    pub file_id: usize,
    /// Retention time in minutes
    pub scan_start_time: f32,
    /// Ion injection time
    pub ion_injection_time: f32,
    /// Selected ions for precursors, if `level > 1`
    pub precursors: Vec<Precursor>,
    /// MS peak masses, sorted in ascending order
    pub masses: Vec<f32>,
    /// MS peak intensities, parallel to `masses`
    pub intensities: Vec<f32>,
    /// MS peak charges, parallel to `masses`
    pub charges: Vec<u8>,
    /// Whether each MS peak charge was assigned from an isotope envelope.
    /// False means the charge is unknown and `charges` contains the fallback value 1.
    pub charge_is_known: Vec<bool>,
    /// Ion mobility values, parallel to `masses` for IMS spectra and empty otherwise
    pub mobilities: Vec<f32>,
    /// Total ion current
    pub total_ion_current: f32,
}

#[derive(Default, Debug, Clone)]
/// An unprocessed mass spectrum, as returned by a parser
/// *CRITICAL*: Users must set all fields manually, including `file_id`
pub struct RawSpectrum {
    pub file_id: usize,
    /// MSn level
    pub ms_level: u8,
    /// Spectrum identifier
    pub id: String,
    /// Vector of precursors associated with this spectrum
    pub precursors: Vec<Precursor>,
    /// Profile or Centroided data
    pub representation: Representation,
    /// Scan start time in minutes
    pub scan_start_time: f32,
    /// Ion injection time
    pub ion_injection_time: f32,
    /// Total ion current
    pub total_ion_current: f32,
    /// M/z array
    pub mz: Vec<f32>,
    /// Intensity array
    pub intensity: Vec<f32>,
    /// Mobility array
    pub mobility: Option<Vec<f32>>,
}

impl RawSpectrum {
    /// Return a [`RawSpectrum`] with default values, but with the `file_id` field
    /// properly set
    pub fn default_with_file_id(file_id: usize) -> Self {
        Self {
            file_id,
            ..Default::default()
        }
    }
}

#[derive(Copy, Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Default)]
pub enum Representation {
    #[default]
    Profile,
    Centroid,
}

/// Binary search followed by linear search to select the most intense peak within `tolerance` window
/// * `offset` - this parameter allows for a static adjustment to the lower and upper bounds of the search window.
///
/// Sage subtracts a proton (and assumes z=1) for all experimental peaks, and stores all fragments as monoisotopic
/// masses. This simplifies downstream calculations at multiple charge states, but it also subtly changes tolerance
/// bounds. For most applications this is completely OK to ignore - however, for exact similarity of TMT reporter ion
/// measurements with ProteomeDiscoverer, FragPipe, etc, we need to account for this minor difference (which has an impact
/// perhaps 0.01% of the time)
pub fn select_most_intense_peak(
    masses: &[f32],
    intensities: &[f32],
    center: f32,
    tolerance: Tolerance,
    offset: Option<f32>,
) -> Option<usize> {
    debug_assert_eq!(masses.len(), intensities.len());
    let (lo, hi) = tolerance.bounds(center);
    let (lo, hi) = (
        lo + offset.unwrap_or_default(),
        hi + offset.unwrap_or_default(),
    );

    // Fragment tolerances cover very few peaks, so a second binary search for
    // the upper bound costs more than scanning forward from the lower bound.
    // Keep the exact inclusive bounds and later-peak tie breaking used here.
    let i = masses
        .partition_point(|mass| mass.total_cmp(&lo).is_lt())
        .saturating_sub(1);

    let mut best_peak = None;
    let mut max_int = 0.0;
    for idx in i..masses.len() {
        if masses[idx].total_cmp(&hi).is_gt() {
            break;
        }
        if masses[idx] < lo || masses[idx] > hi {
            continue;
        }
        if intensities[idx] >= max_int {
            max_int = intensities[idx];
            best_peak = Some(idx);
        }
    }
    best_peak
}

// pub fn find_spectrum_by_id(
//     spectra: &[ProcessedSpectrum],
//     scan_id: usize,
// ) -> Option<&ProcessedSpectrum> {
//     // First try indexing by scan
//     if let Some(first) = spectra.get(scan_id.saturating_sub(1)) {
//         if first.scan == scan_id {
//             return Some(first);
//         }
//     }
//     // Fall back to binary search
//     let idx = spectra
//         .binary_search_by(|spec| spec.scan.cmp(&scan_id))
//         .ok()?;
//     spectra.get(idx)
// }

#[derive(Debug)]
struct EnvelopeCandidate {
    indices: [usize; 4],
    len: usize,
    charge: u8,
    score: f32,
    mean_ppm_error: f32,
}

fn isotope_pattern_score(observed: &[f32], theoretical: &[f32; 4]) -> f32 {
    let observed_sum = observed.iter().sum::<f32>();
    let theoretical_coverage = theoretical[..observed.len()].iter().sum::<f32>();
    if observed_sum <= 0.0 || theoretical_coverage <= 0.0 {
        return 0.0;
    }

    let coefficient = observed
        .iter()
        .zip(theoretical)
        .map(|(&observed, &theoretical)| {
            ((observed / observed_sum) * (theoretical / theoretical_coverage)).sqrt()
        })
        .sum::<f32>();
    (coefficient * theoretical_coverage).clamp(0.0, 1.0)
}

/// Deisotope peaks with a bounded averagine-scored candidate search.
///
/// Every peak is considered as a monoisotopic seed for every allowed charge.
/// Candidates grow upward through at most four envelope peaks. Accepted
/// candidates claim peaks exclusively, which preserves total intensity and
/// makes charge assignment independent of input traversal side effects.
pub fn deisotope(
    mz: &[f32],
    int: &[f32],
    max_charge: u8,
    settings: DeisotopeSettings,
    min_mz: f32,
) -> Vec<Deisotoped> {
    debug_assert_eq!(mz.len(), int.len());
    debug_assert!(settings.validate().is_ok());
    let max_envelope_peaks = settings.max_envelope_peaks.clamp(2, 4);

    let sorted = mz
        .windows(2)
        .all(|pair| pair[0].total_cmp(&pair[1]).is_le());
    let (sorted_mz, sorted_int): (Cow<'_, [f32]>, Cow<'_, [f32]>) = if sorted {
        (Cow::Borrowed(mz), Cow::Borrowed(int))
    } else {
        let mut order = (0..mz.len()).collect::<Vec<_>>();
        order.sort_unstable_by(|&left, &right| mz[left].total_cmp(&mz[right]));
        (
            Cow::Owned(order.iter().map(|&idx| mz[idx]).collect()),
            Cow::Owned(order.iter().map(|&idx| int[idx]).collect()),
        )
    };
    let mut candidates = Vec::with_capacity(sorted_mz.len() * usize::from(max_charge));

    for seed in 0..sorted_mz.len() {
        if sorted_mz[seed] < min_mz {
            continue;
        }
        for charge in 1..=max_charge {
            let neutral_mass = (sorted_mz[seed] - PROTON) * f32::from(charge);
            if neutral_mass <= 0.0 {
                continue;
            }
            let theoretical = crate::isotopes::averagine_isotopes_cached(neutral_mass);
            let mut indices = [seed; 4];
            let mut len = 1usize;
            let mut ppm_error_sum = 0.0f32;

            for isotope in 1..max_envelope_peaks {
                let target = sorted_mz[seed] + NEUTRON * isotope as f32 / f32::from(charge);
                let tolerance = Tolerance::ppm_to_delta_mass(target, settings.ppm_tolerance);
                let lo = sorted_mz.partition_point(|&value| value < target - tolerance);
                let hi = sorted_mz.partition_point(|&value| value <= target + tolerance);
                let previous = indices[len - 1];
                let expected_ratio = theoretical[isotope] / theoretical[isotope - 1].max(1e-12);

                let best = (lo..hi)
                    .filter(|idx| !indices[..len].contains(idx))
                    .filter_map(|idx| {
                        let previous_intensity = sorted_int[previous];
                        let candidate_intensity = sorted_int[idx];
                        if previous_intensity <= 0.0 || candidate_intensity <= 0.0 {
                            return None;
                        }
                        let observed_ratio = candidate_intensity / previous_intensity;
                        let ratio_error = (observed_ratio / expected_ratio.max(1e-12)).log2().abs();
                        if ratio_error > settings.max_isotope_log2_ratio {
                            return None;
                        }
                        let ppm_error = (sorted_mz[idx] - target).abs() * 1e6 / target;
                        Some((idx, ppm_error, ratio_error))
                    })
                    .min_by(|left, right| {
                        left.1
                            .total_cmp(&right.1)
                            .then_with(|| left.2.total_cmp(&right.2))
                            .then_with(|| left.0.cmp(&right.0))
                    });

                let Some((idx, ppm_error, _)) = best else {
                    break;
                };
                indices[len] = idx;
                len += 1;
                ppm_error_sum += ppm_error;
            }

            if len < settings.min_envelope_peaks {
                continue;
            }
            let mut observed = [0.0f32; 4];
            for position in 0..len {
                observed[position] = sorted_int[indices[position]];
            }
            candidates.push(EnvelopeCandidate {
                indices,
                len,
                charge,
                score: isotope_pattern_score(&observed[..len], &theoretical),
                mean_ppm_error: ppm_error_sum / (len - 1) as f32,
            });
        }
    }

    candidates.sort_unstable_by(|left, right| {
        right
            .len
            .cmp(&left.len)
            .then_with(|| right.score.total_cmp(&left.score))
            .then_with(|| left.mean_ppm_error.total_cmp(&right.mean_ppm_error))
            .then_with(|| left.indices[0].cmp(&right.indices[0]))
            .then_with(|| left.charge.cmp(&right.charge))
    });

    let mut peaks = sorted_mz
        .iter()
        .zip(sorted_int.iter())
        .map(|(&mz, &intensity)| Deisotoped {
            mz,
            intensity,
            charge: None,
            envelope: None,
        })
        .collect::<Vec<_>>();
    let mut claimed = vec![false; peaks.len()];
    for candidate in candidates {
        if candidate.score < settings.min_score
            || candidate.indices[..candidate.len]
                .iter()
                .any(|&idx| claimed[idx])
        {
            continue;
        }

        let root = candidate.indices[0];
        let intensity = candidate.indices[..candidate.len]
            .iter()
            .map(|&idx| sorted_int[idx])
            .sum::<f32>();
        peaks[root].intensity = intensity;
        peaks[root].charge = Some(candidate.charge);
        claimed[root] = true;
        for &idx in &candidate.indices[1..candidate.len] {
            peaks[idx].charge = Some(candidate.charge);
            peaks[idx].envelope = Some(root);
            claimed[idx] = true;
        }
    }
    peaks
}

fn retain_top_n_by_intensity(
    indices: &mut Vec<usize>,
    masses: &[f32],
    intensities: &[f32],
    n: usize,
    prefer_low_mass: bool,
) {
    debug_assert_eq!(masses.len(), intensities.len());
    if n == 0 {
        indices.clear();
        return;
    }

    if indices.len() <= n {
        return;
    }

    let keep_from = indices.len() - n;
    indices.select_nth_unstable_by(keep_from, |&a, &b| {
        intensities[a].total_cmp(&intensities[b]).then_with(|| {
            if prefer_low_mass {
                masses[b].total_cmp(&masses[a])
            } else {
                masses[a].total_cmp(&masses[b])
            }
        })
    });
    indices.drain(..keep_from);
}

fn select_columns(
    indices: Vec<usize>,
    masses: &[f32],
    intensities: &[f32],
    charges: &[u8],
    charge_is_known: &[bool],
) -> (Vec<f32>, Vec<f32>, Vec<u8>, Vec<bool>) {
    debug_assert_eq!(masses.len(), intensities.len());
    debug_assert_eq!(masses.len(), charges.len());
    debug_assert_eq!(masses.len(), charge_is_known.len());

    let mut selected_masses = Vec::with_capacity(indices.len());
    let mut selected_intensities = Vec::with_capacity(indices.len());
    let mut selected_charges = Vec::with_capacity(indices.len());
    let mut selected_charge_is_known = Vec::with_capacity(indices.len());

    for idx in indices {
        selected_masses.push(masses[idx]);
        selected_intensities.push(intensities[idx]);
        selected_charges.push(charges[idx]);
        selected_charge_is_known.push(charge_is_known[idx]);
    }

    (
        selected_masses,
        selected_intensities,
        selected_charges,
        selected_charge_is_known,
    )
}

type SortedColumns = (Vec<f32>, Vec<f32>, Vec<u8>, Vec<bool>, Vec<f32>);

fn sort_columns_by_mass(
    masses: Vec<f32>,
    intensities: Vec<f32>,
    charges: Vec<u8>,
    charge_is_known: Vec<bool>,
    mobilities: Vec<f32>,
) -> SortedColumns {
    debug_assert_eq!(masses.len(), intensities.len());
    debug_assert_eq!(masses.len(), charges.len());
    debug_assert_eq!(masses.len(), charge_is_known.len());
    debug_assert!(mobilities.is_empty() || mobilities.len() == masses.len());

    if masses.len() <= 1 {
        return (masses, intensities, charges, charge_is_known, mobilities);
    }

    if masses.windows(2).all(|window| window[0] <= window[1]) {
        return (masses, intensities, charges, charge_is_known, mobilities);
    }

    let has_mobility = !mobilities.is_empty();
    let mut order = (0..masses.len()).collect::<Vec<_>>();
    order.sort_by(|&a, &b| masses[a].total_cmp(&masses[b]));

    let mut sorted_masses = Vec::with_capacity(masses.len());
    let mut sorted_intensities = Vec::with_capacity(intensities.len());
    let mut sorted_charges = Vec::with_capacity(charges.len());
    let mut sorted_charge_is_known = Vec::with_capacity(charge_is_known.len());
    let mut sorted_mobilities = Vec::with_capacity(mobilities.len());

    for idx in order {
        sorted_masses.push(masses[idx]);
        sorted_intensities.push(intensities[idx]);
        sorted_charges.push(charges[idx]);
        sorted_charge_is_known.push(charge_is_known[idx]);
        if has_mobility {
            sorted_mobilities.push(mobilities[idx]);
        }
    }

    (
        sorted_masses,
        sorted_intensities,
        sorted_charges,
        sorted_charge_is_known,
        sorted_mobilities,
    )
}

impl ProcessedSpectrum {
    pub fn len(&self) -> usize {
        self.masses.len()
    }

    pub fn is_empty(&self) -> bool {
        self.masses.is_empty()
    }

    pub fn peak_mz(&self, idx: usize) -> f32 {
        self.masses[idx] / self.charges[idx] as f32 + PROTON
    }

    pub fn has_known_charge(&self, idx: usize) -> bool {
        self.charge_is_known.get(idx).copied().unwrap_or(false)
    }

    pub fn extract_ms1_precursor(&self) -> Option<(f32, u8)> {
        let precursor = self.precursors.first()?;
        let charge = precursor.charge?;
        let mass = (precursor.mz - PROTON) * charge as f32;
        Some((mass, charge))
    }
    pub fn in_isolation_window(&self, mz: f32) -> Option<bool> {
        let precursor = self.precursors.first()?;
        let (lo, hi) = precursor.isolation_window?.bounds(precursor.mz - PROTON);
        Some(mz >= lo && mz <= hi)
    }
}

impl SpectrumProcessor {
    /// Create a new [`SpectrumProcessor`]
    ///
    /// # Arguments
    /// * `take_top_n`: Keep only the top N most intense peaks from the spectrum
    /// * `min_fragment_mz`: Keep only fragments >= this m/z
    /// * `max_fragment_mz`: Keep only fragments <= this m/z
    /// * `deisotope`: Perform deisotoping & charge state deconvolution
    pub fn new(take_top_n: usize, deisotope: bool, min_deisotope_mz: f32) -> Self {
        Self {
            take_top_n,
            min_deisotope_mz,
            deisotope: DeisotopeSettings::from_enabled(deisotope),
        }
    }

    pub fn with_deisotope_settings(
        take_top_n: usize,
        deisotope: DeisotopeSettings,
        min_deisotope_mz: f32,
    ) -> Self {
        Self {
            take_top_n,
            min_deisotope_mz,
            deisotope,
        }
    }

    fn process_ms2(
        &self,
        should_deisotope: bool,
        spectrum: &RawSpectrum,
    ) -> (Vec<f32>, Vec<f32>, Vec<u8>, Vec<bool>) {
        if spectrum.representation != Representation::Centroid {
            // Panic, because there's really nothing we can do with profile data
            panic!(
                "Scan {} contains profile data! Please convert to centroid",
                spectrum.id
            );
        }

        // If there is no precursor charge from the mzML file, then deisotope fragments up to z=3
        let charge = spectrum
            .precursors
            .first()
            .and_then(|p| p.charge)
            .unwrap_or(3);

        if should_deisotope {
            let charge = self
                .deisotope
                .max_charge
                .map(|configured| charge.min(configured))
                .unwrap_or(charge)
                .max(1);
            let peaks = deisotope(
                &spectrum.mz,
                &spectrum.intensity,
                charge,
                self.deisotope,
                self.min_deisotope_mz,
            );
            let mz = peaks.iter().map(|peak| peak.mz).collect::<Vec<_>>();
            let intensities = peaks.iter().map(|peak| peak.intensity).collect::<Vec<_>>();
            let charges = peaks
                .iter()
                .map(|peak| peak.charge.unwrap_or(1))
                .collect::<Vec<_>>();
            let charge_is_known = peaks
                .iter()
                .map(|peak| peak.charge.is_some())
                .collect::<Vec<_>>();
            let mut indices = peaks
                .iter()
                .enumerate()
                .filter_map(|(idx, peak)| peak.envelope.is_none().then_some(idx))
                .collect::<Vec<_>>();

            retain_top_n_by_intensity(&mut indices, &mz, &intensities, self.take_top_n, true);

            let mut masses = Vec::with_capacity(indices.len());
            let mut selected_intensities = Vec::with_capacity(indices.len());
            let mut selected_charges = Vec::with_capacity(indices.len());
            let mut selected_charge_is_known = Vec::with_capacity(indices.len());
            for idx in indices {
                let charge = charges[idx];
                masses.push((mz[idx] - PROTON) * charge as f32);
                selected_intensities.push(intensities[idx]);
                selected_charges.push(charge);
                selected_charge_is_known.push(charge_is_known[idx]);
            }

            (
                masses,
                selected_intensities,
                selected_charges,
                selected_charge_is_known,
            )
        } else {
            let masses = spectrum
                .mz
                .iter()
                .map(|mz| (mz - PROTON) * 1.0)
                .collect::<Vec<_>>();
            let intensities = spectrum.intensity.clone();
            let charges = vec![1; masses.len()];
            let charge_is_known = vec![false; masses.len()];
            let mut indices = (0..masses.len()).collect::<Vec<_>>();
            retain_top_n_by_intensity(&mut indices, &masses, &intensities, self.take_top_n, false);
            select_columns(indices, &masses, &intensities, &charges, &charge_is_known)
        }
    }

    pub fn process(&self, mut spectrum: RawSpectrum) -> ProcessedSpectrum {
        debug_assert_eq!(spectrum.mz.len(), spectrum.intensity.len());
        if let Some(mobilities) = spectrum.mobility.as_ref() {
            debug_assert_eq!(spectrum.mz.len(), mobilities.len());
        }

        let (masses, intensities, charges, charge_is_known, mobilities) = if spectrum.ms_level == 1
            && spectrum.mobility.is_some()
        {
            let mut masses = spectrum.mz;
            masses.iter_mut().for_each(|mass| *mass -= PROTON);
            let intensities = spectrum.intensity;
            let charges = vec![1; masses.len()];
            let charge_is_known = vec![false; masses.len()];
            let mobilities = spectrum.mobility.take().expect("checked above");
            sort_columns_by_mass(masses, intensities, charges, charge_is_known, mobilities)
        } else {
            let (masses, intensities, charges, charge_is_known) = match spectrum.ms_level {
                2 if self.deisotope.enabled => self.process_ms2(true, &spectrum),
                2 => {
                    if spectrum.representation != Representation::Centroid {
                        panic!(
                            "Scan {} contains profile data! Please convert to centroid",
                            spectrum.id
                        );
                    }
                    let mut masses = spectrum.mz;
                    masses.iter_mut().for_each(|mass| *mass -= PROTON);
                    let intensities = spectrum.intensity;
                    let charges = vec![1; masses.len()];
                    let charge_is_known = vec![false; masses.len()];
                    if masses.len() <= self.take_top_n {
                        (masses, intensities, charges, charge_is_known)
                    } else {
                        let mut indices = (0..masses.len()).collect::<Vec<_>>();
                        retain_top_n_by_intensity(
                            &mut indices,
                            &masses,
                            &intensities,
                            self.take_top_n,
                            false,
                        );
                        select_columns(indices, &masses, &intensities, &charges, &charge_is_known)
                    }
                }
                _ => {
                    let mut masses = spectrum.mz;
                    masses.iter_mut().for_each(|mass| *mass -= PROTON);
                    let intensities = spectrum.intensity;
                    let charges = vec![1; masses.len()];
                    let charge_is_known = vec![false; masses.len()];
                    (masses, intensities, charges, charge_is_known)
                }
            };
            sort_columns_by_mass(masses, intensities, charges, charge_is_known, Vec::new())
        };

        let total_ion_current = intensities.iter().sum::<f32>();

        ProcessedSpectrum {
            level: spectrum.ms_level,
            id: spectrum.id,
            file_id: spectrum.file_id,
            scan_start_time: spectrum.scan_start_time,
            ion_injection_time: spectrum.ion_injection_time,
            precursors: spectrum.precursors,
            masses,
            intensities,
            charges,
            charge_is_known,
            mobilities,
            total_ion_current,
        }
    }
}

#[cfg(test)]
#[path = "../tests/unit/spectrum.rs"]
mod test;
