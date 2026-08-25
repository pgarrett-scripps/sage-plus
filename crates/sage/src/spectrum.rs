use crate::database::binary_search_slice;
use crate::mass::{Tolerance, NEUTRON, PROTON};

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
    pub deisotope: bool,
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

    let (i, j) = binary_search_slice(masses, |mass, query| mass.total_cmp(query), lo, hi);

    let mut best_peak = None;
    let mut max_int = 0.0;
    for idx in (i..j).filter(|&idx| masses[idx] >= lo && masses[idx] <= hi) {
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

/// Deisotope a set of peaks by attempting to find C13 peaks under a given `ppm` tolerance
pub fn deisotope(
    mz: &[f32],
    int: &[f32],
    max_charge: u8,
    ppm: f32,
    min_mz: f32,
) -> Vec<Deisotoped> {
    let mut peaks = mz
        .iter()
        .zip(int.iter())
        .map(|(mz, int)| Deisotoped {
            mz: *mz,
            intensity: *int,
            envelope: None,
            charge: None,
        })
        .collect::<Vec<_>>();

    // Is the peak at index `i` an isotopic peak?
    for i in (0..mz.len()).rev() {
        // Two pointer approach, j is fast pointer
        let mut j = i.saturating_sub(1);
        while mz[i] - mz[j] <= NEUTRON + Tolerance::ppm_to_delta_mass(mz[i], ppm) && mz[j] >= min_mz
        {
            let delta = mz[i] - mz[j];
            let tol = Tolerance::ppm_to_delta_mass(mz[i], ppm);
            for charge in 1..=max_charge {
                let iso = NEUTRON / charge as f32;
                if (delta - iso).abs() <= tol && int[i] < int[j] {
                    // Make sure this peak isn't already part of an isotopic envelope
                    if let Some(existing) = peaks[i].charge {
                        if existing != charge {
                            continue;
                        }
                    }
                    peaks[j].intensity += peaks[i].intensity;
                    peaks[j].charge = Some(charge);
                    peaks[i].charge = Some(charge);
                    peaks[i].envelope = Some(j);
                }
            }
            j = j.saturating_sub(1);
            if j == 0 {
                break;
            }
        }
    }
    peaks
}

/// Path compression of isotopic envelope links
pub fn path_compression(peaks: &mut [Deisotoped]) {
    for idx in 0..peaks.len() {
        if let Some(parent) = peaks[idx].envelope {
            if let Some(upper) = peaks[parent].envelope {
                peaks[idx].envelope = Some(upper);
            }
            peaks[idx].intensity = 0.0;
        }
    }
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
) -> (Vec<f32>, Vec<f32>, Vec<u8>) {
    debug_assert_eq!(masses.len(), intensities.len());
    debug_assert_eq!(masses.len(), charges.len());

    let mut selected_masses = Vec::with_capacity(indices.len());
    let mut selected_intensities = Vec::with_capacity(indices.len());
    let mut selected_charges = Vec::with_capacity(indices.len());

    for idx in indices {
        selected_masses.push(masses[idx]);
        selected_intensities.push(intensities[idx]);
        selected_charges.push(charges[idx]);
    }

    (selected_masses, selected_intensities, selected_charges)
}

fn sort_columns_by_mass(
    masses: Vec<f32>,
    intensities: Vec<f32>,
    charges: Vec<u8>,
    mobilities: Vec<f32>,
) -> (Vec<f32>, Vec<f32>, Vec<u8>, Vec<f32>) {
    debug_assert_eq!(masses.len(), intensities.len());
    debug_assert_eq!(masses.len(), charges.len());
    debug_assert!(mobilities.is_empty() || mobilities.len() == masses.len());

    if masses.len() <= 1 {
        return (masses, intensities, charges, mobilities);
    }

    if masses.windows(2).all(|window| window[0] <= window[1]) {
        return (masses, intensities, charges, mobilities);
    }

    let has_mobility = !mobilities.is_empty();
    let mut order = (0..masses.len()).collect::<Vec<_>>();
    order.sort_by(|&a, &b| masses[a].total_cmp(&masses[b]));

    let mut sorted_masses = Vec::with_capacity(masses.len());
    let mut sorted_intensities = Vec::with_capacity(intensities.len());
    let mut sorted_charges = Vec::with_capacity(charges.len());
    let mut sorted_mobilities = Vec::with_capacity(mobilities.len());

    for idx in order {
        sorted_masses.push(masses[idx]);
        sorted_intensities.push(intensities[idx]);
        sorted_charges.push(charges[idx]);
        if has_mobility {
            sorted_mobilities.push(mobilities[idx]);
        }
    }

    (
        sorted_masses,
        sorted_intensities,
        sorted_charges,
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
            deisotope,
        }
    }

    fn process_ms2(
        &self,
        should_deisotope: bool,
        spectrum: &RawSpectrum,
    ) -> (Vec<f32>, Vec<f32>, Vec<u8>) {
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
            let peaks = deisotope(
                &spectrum.mz,
                &spectrum.intensity,
                charge,
                10.0,
                self.min_deisotope_mz,
            );
            let mz = peaks.iter().map(|peak| peak.mz).collect::<Vec<_>>();
            let intensities = peaks.iter().map(|peak| peak.intensity).collect::<Vec<_>>();
            let charges = peaks
                .iter()
                .map(|peak| peak.charge.unwrap_or(1))
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
            for idx in indices {
                let charge = charges[idx];
                masses.push((mz[idx] - PROTON) * charge as f32);
                selected_intensities.push(intensities[idx]);
                selected_charges.push(charge);
            }

            (masses, selected_intensities, selected_charges)
        } else {
            let masses = spectrum
                .mz
                .iter()
                .map(|mz| (mz - PROTON) * 1.0)
                .collect::<Vec<_>>();
            let intensities = spectrum.intensity.clone();
            let charges = vec![1; masses.len()];
            let mut indices = (0..masses.len()).collect::<Vec<_>>();
            retain_top_n_by_intensity(&mut indices, &masses, &intensities, self.take_top_n, false);
            select_columns(indices, &masses, &intensities, &charges)
        }
    }

    pub fn process(&self, mut spectrum: RawSpectrum) -> ProcessedSpectrum {
        debug_assert_eq!(spectrum.mz.len(), spectrum.intensity.len());
        if let Some(mobilities) = spectrum.mobility.as_ref() {
            debug_assert_eq!(spectrum.mz.len(), mobilities.len());
        }

        let (masses, intensities, charges, mobilities) =
            if spectrum.ms_level == 1 && spectrum.mobility.is_some() {
                let mut masses = spectrum.mz;
                masses.iter_mut().for_each(|mass| *mass -= PROTON);
                let intensities = spectrum.intensity;
                let charges = vec![1; masses.len()];
                let mobilities = spectrum.mobility.take().expect("checked above");
                sort_columns_by_mass(masses, intensities, charges, mobilities)
            } else {
                let (masses, intensities, charges) = match spectrum.ms_level {
                    2 if self.deisotope => self.process_ms2(true, &spectrum),
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
                        if masses.len() <= self.take_top_n {
                            (masses, intensities, charges)
                        } else {
                            let mut indices = (0..masses.len()).collect::<Vec<_>>();
                            retain_top_n_by_intensity(
                                &mut indices,
                                &masses,
                                &intensities,
                                self.take_top_n,
                                false,
                            );
                            select_columns(indices, &masses, &intensities, &charges)
                        }
                    }
                    _ => {
                        let mut masses = spectrum.mz;
                        masses.iter_mut().for_each(|mass| *mass -= PROTON);
                        let intensities = spectrum.intensity;
                        let charges = vec![1; masses.len()];
                        (masses, intensities, charges)
                    }
                };
                sort_columns_by_mass(masses, intensities, charges, Vec::new())
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
            mobilities,
            total_ion_current,
        }
    }
}

#[cfg(test)]
#[path = "../tests/unit/spectrum.rs"]
mod test;
