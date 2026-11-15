use std::io;
use std::path::Path;

use mzdata::io::{DetailLevel, MzMLbReader as MzMLbReaderImpl};
use mzdata::params::ParamValue;
use mzdata::prelude::{
    IonMobilityMeasure, MSDataFileMetadata, ParamDescribed, PrecursorSelection, SpectrumLike,
    SpectrumSource,
};
use mzdata::spectrum::{MultiLayerSpectrum, RawSpectrum as MzDataSpectrum, SignalContinuity};
use sage_core::mass::Tolerance;
use sage_core::spectrum::{
    AcquisitionGroup, Activation, MassAnalyzer, Precursor, RawSpectrum, Representation,
};

pub struct MzMLbReader {
    file_id: usize,
    /// Whether MS1 arrays are read from the HDF5 file, as in
    /// [`crate::mzml::MzMLReader::set_decode_ms1`].
    decode_ms1: bool,
}

/// PSI-MS CURIE for an accession number, e.g. `MS:1000484`.
fn psi_ms(accession: u32) -> String {
    format!("MS:{accession:07}")
}

fn normalize_fragment_charges(
    charges: &[i32],
    peak_count: usize,
    spectrum_id: &str,
) -> Option<Vec<u8>> {
    if charges.len() != peak_count {
        log::warn!(
            "ignoring fragment charge array with length different from m/z array for spectrum {}",
            spectrum_id
        );
        return None;
    }
    let converted = charges
        .iter()
        .map(|&charge| u8::try_from(charge).ok())
        .collect::<Option<Vec<_>>>();
    if converted.is_none() {
        log::warn!(
            "ignoring fragment charge array with values outside 0 through 255 for spectrum {}",
            spectrum_id
        );
    }
    converted
}

impl MzMLbReader {
    pub fn with_file_id(file_id: usize) -> Self {
        Self {
            file_id,
            decode_ms1: true,
        }
    }

    /// Read MS1 arrays (the default), or keep MS1 spectra without peaks, as
    /// [`crate::mzml::MzMLReader::set_decode_ms1`] does for mzML.
    pub fn set_decode_ms1(&mut self, decode_ms1: bool) -> &mut Self {
        self.decode_ms1 = decode_ms1;
        self
    }

    pub fn parse(&self, path: impl AsRef<Path>) -> io::Result<Vec<RawSpectrum>> {
        let mut reader = MzMLbReaderImpl::new(&path.as_ref().to_path_buf())?;
        // The analyzer of each instrument configuration, as in mzML.
        let analyzers = reader
            .instrument_configurations()
            .iter()
            .map(|(&id, configuration)| {
                let analyzer = configuration
                    .components
                    .iter()
                    .filter_map(|component| component.mass_analyzer())
                    .filter_map(|term| MassAnalyzer::from_psi_ms(&psi_ms(term.accession())))
                    .fold(MassAnalyzer::Unknown, MassAnalyzer::combine);
                (id, analyzer)
            })
            .collect::<std::collections::HashMap<_, _>>();
        let convert = |spectrum: MultiLayerSpectrum| {
            let ms_level = spectrum.ms_level();
            let id = spectrum.id().to_string();
            let scan_start_time = spectrum.start_time() as f32;
            let ion_injection_time = spectrum
                .acquisition()
                .first_scan()
                .map(|scan| scan.injection_time)
                .unwrap_or_default();
            let scan_mobility = spectrum.ion_mobility().map(|value| value as f32);
            let scan = spectrum.acquisition().first_scan();
            let activation = spectrum
                .precursor_iter()
                .flat_map(|precursor| precursor.activation.methods())
                .filter_map(|method| Activation::from_psi_ms(&psi_ms(method.accession())))
                .fold(Activation::Unknown, Activation::combine);
            let acquisition = AcquisitionGroup::resolve(
                scan.and_then(|scan| scan.filter_string()).as_deref(),
                scan.and_then(|scan| analyzers.get(&scan.instrument_configuration_id))
                    .copied()
                    .unwrap_or_default(),
                activation,
            );
            let representation = match spectrum.signal_continuity() {
                SignalContinuity::Centroid => Representation::Centroid,
                SignalContinuity::Profile | SignalContinuity::Unknown => Representation::Profile,
            };
            let total_ion_current = spectrum.peaks().tic();
            let precursors = spectrum
                .precursor_iter()
                .flat_map(|precursor| {
                    let window = precursor.isolation_window();
                    let isolation_window = (!window.is_empty()).then_some({
                        Tolerance::Da(
                            window.lower_bound - window.target,
                            window.upper_bound - window.target,
                        )
                    });
                    let spectrum_ref = precursor.precursor_id().cloned();
                    let ions = precursor
                        .iter()
                        .map(|ion| Precursor {
                            mz: ion.mz as f32,
                            intensity: Some(ion.intensity),
                            charge: ion.charge.and_then(|charge| u8::try_from(charge).ok()),
                            spectrum_ref: spectrum_ref.clone(),
                            isolation_window,
                            inverse_ion_mobility: ion
                                .ion_mobility()
                                .map(|value| value as f32)
                                .or(scan_mobility),
                        })
                        .collect::<Vec<_>>();
                    if ions.is_empty() && window.target > 0.0 {
                        vec![Precursor {
                            mz: window.target,
                            spectrum_ref,
                            isolation_window,
                            inverse_ion_mobility: scan_mobility,
                            ..Precursor::default()
                        }]
                    } else {
                        ions
                    }
                })
                .collect();
            let raw: MzDataSpectrum = spectrum.into();
            let fragment_charges = raw
                .arrays
                .charges()
                .ok()
                .and_then(|charges| normalize_fragment_charges(&charges, raw.mzs().len(), &id));

            RawSpectrum {
                file_id: self.file_id,
                ms_level,
                id,
                precursors,
                representation,
                scan_start_time,
                ion_injection_time,
                total_ion_current,
                mz: raw.mzs().iter().map(|value| *value as f32).collect(),
                intensity: raw.intensities().to_vec(),
                fragment_charges,
                mobility: None,
                acquisition,
                metadata_only: false,
            }
        };
        // mzdata chooses whether to read a spectrum's arrays before it parses
        // the spectrum, and the MS level is inside it, so a search that reads
        // only MS2 reads each spectrum's metadata first, through the index,
        // and the arrays only for a spectrum that is not MS1. Without an index
        // there is nothing to seek to, and every spectrum is read whole.
        if self.decode_ms1 || reader.get_index().is_empty() {
            return Ok(reader.map(convert).collect());
        }
        let mut spectra = Vec::with_capacity(reader.len());
        let mut undecoded_ms1 = 0usize;
        for index in 0..reader.len() {
            reader.set_detail_level(DetailLevel::MetadataOnly);
            let Some(metadata) = reader.get_spectrum_by_index(index) else {
                continue;
            };
            if metadata.ms_level() == 1 {
                spectra.push(self.undecoded(metadata, &analyzers));
                undecoded_ms1 += 1;
                continue;
            }
            reader.set_detail_level(DetailLevel::Full);
            if let Some(spectrum) = reader.get_spectrum_by_index(index) {
                spectra.push(convert(spectrum));
            }
        }
        if undecoded_ms1 > 0 {
            log::debug!(
                "{}: {} MS1 spectra were not decoded; this search reads only MS2",
                path.as_ref().display(),
                undecoded_ms1
            );
        }
        Ok(spectra)
    }

    /// An MS1 spectrum kept without its arrays: the metadata fields `parse`
    /// fills for a decoded MS1 spectrum, the total ion current from its
    /// `MS:1000285` parameter, and no peaks. It dates the file for quality
    /// control and measures nothing.
    fn undecoded(
        &self,
        spectrum: MultiLayerSpectrum,
        analyzers: &std::collections::HashMap<u32, MassAnalyzer>,
    ) -> RawSpectrum {
        let scan = spectrum.acquisition().first_scan();
        RawSpectrum {
            file_id: self.file_id,
            ms_level: spectrum.ms_level(),
            id: spectrum.id().to_string(),
            precursors: Vec::new(),
            representation: match spectrum.signal_continuity() {
                SignalContinuity::Centroid => Representation::Centroid,
                SignalContinuity::Profile | SignalContinuity::Unknown => Representation::Profile,
            },
            scan_start_time: spectrum.start_time() as f32,
            ion_injection_time: scan.map(|scan| scan.injection_time).unwrap_or_default(),
            total_ion_current: spectrum
                .get_param_by_curie(&mzdata::curie!(MS:1000285))
                .and_then(|param| param.to_f32().ok())
                .unwrap_or_default(),
            mz: Vec::new(),
            intensity: Vec::new(),
            fragment_charges: None,
            mobility: None,
            acquisition: AcquisitionGroup::resolve(
                scan.and_then(|scan| scan.filter_string()).as_deref(),
                scan.and_then(|scan| analyzers.get(&scan.instrument_configuration_id))
                    .copied()
                    .unwrap_or_default(),
                Activation::Unknown,
            ),
            metadata_only: true,
        }
    }
}

#[cfg(test)]
#[path = "../tests/unit/mzmlb.rs"]
mod test;
