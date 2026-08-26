use async_compression::tokio::bufread::ZlibDecoder;
use quick_xml::events::Event;
use quick_xml::Reader;
use sage_core::spectrum::{Precursor, Representation};
use sage_core::{mass::Tolerance, spectrum::RawSpectrum};
use std::collections::HashMap;
use tokio::io::{AsyncBufRead, AsyncReadExt};

#[derive(Copy, Clone, Debug, PartialEq, Eq, PartialOrd, Ord)]
/// Which tag are we inside?
enum State {
    Spectrum,
    Scan,
    BinaryDataArray,
    Binary,
    Precursor,
    SelectedIon,
}

#[derive(Copy, Clone, Debug)]
enum BinaryKind {
    Intensity,
    Mz,
    Noise,
}

#[derive(Copy, Clone, Debug)]
enum Dtype {
    F32,
    F64,
}

#[derive(Clone, Debug)]
struct CvParam {
    accession: Vec<u8>,
    value: Option<String>,
    unit_accession: Option<Vec<u8>>,
}

impl CvParam {
    fn from_event(ev: &quick_xml::events::BytesStart<'_>) -> Result<Self, MzMLError> {
        let accession = ev
            .try_get_attribute(b"accession")?
            .ok_or(MzMLError::Malformed)?
            .value
            .into_owned();
        let value = ev
            .try_get_attribute(b"value")?
            .map(|attribute| std::str::from_utf8(&attribute.value).map(str::to_owned))
            .transpose()?;
        let unit_accession = ev
            .try_get_attribute(b"unitAccession")?
            .map(|attribute| attribute.value.into_owned());
        Ok(Self {
            accession,
            value,
            unit_accession,
        })
    }

    fn parse_value<T>(&self) -> Result<T, MzMLError>
    where
        T: std::str::FromStr,
        MzMLError: From<T::Err>,
    {
        self.value
            .as_deref()
            .ok_or(MzMLError::Malformed)?
            .parse()
            .map_err(Into::into)
    }
}

// MUST supply only one of the following
const ZLIB_COMPRESSION: &[u8] = b"MS:1000574";
const NO_COMPRESSION: &[u8] = b"MS:1000576";

// MUST supply only one of the following
const INTENSITY_ARRAY: &[u8] = b"MS:1000515";
const MZ_ARRAY: &[u8] = b"MS:1000514";
const NOISE_ARRAY: &[u8] = b"MS:1002744";

// MUST supply only one of the following
const FLOAT_64: &[u8] = b"MS:1000523";
const FLOAT_32: &[u8] = b"MS:1000521";

const MS_LEVEL: &[u8] = b"MS:1000511";
const PROFILE: &[u8] = b"MS:1000128";
const CENTROID: &[u8] = b"MS:1000127";
const TOTAL_ION_CURRENT: &[u8] = b"MS:1000285";

const SCAN_START_TIME: &[u8] = b"MS:1000016";
const UNIT_SECONDS: &[u8] = b"UO:0000010";
const UNIT_MINUTES: &[u8] = b"UO:0000031";
const ION_INJECTION_TIME: &[u8] = b"MS:1000927";

const SELECTED_ION_MZ: &[u8] = b"MS:1000744";
const SELECTED_ION_INT: &[u8] = b"MS:1000042";
const SELECTED_ION_CHARGE: &[u8] = b"MS:1000041";

const ISO_WINDOW_TARGET: &[u8] = b"MS:1000827";
const ISO_WINDOW_LOWER: &[u8] = b"MS:1000828";
const ISO_WINDOW_UPPER: &[u8] = b"MS:1000829";

const INVERSE_ION_MOBILITY: &[u8] = b"MS:1002815";

pub struct MzMLReader {
    ms_level: Option<u8>,
    // If set to Some(level) and noise intensities are present in the MzML file,
    // divide intensities at this MS-level by noise to calculate S/N
    signal_to_noise: Option<u8>,

    file_id: usize,
}

impl MzMLReader {
    /// Create a new [`MzMlReader`] with a minimum MS level filter
    ///
    /// # Example
    ///
    /// A minimum level of 2 will not parse or return MS1 scans
    pub fn with_file_id_and_level_filter(file_id: usize, ms_level: u8) -> Self {
        Self {
            ms_level: Some(ms_level),
            file_id,
            signal_to_noise: None,
        }
    }

    pub fn with_file_id(file_id: usize) -> Self {
        Self {
            ms_level: None,
            signal_to_noise: None,
            file_id,
        }
    }

    pub fn set_file_id(&mut self, file_id: usize) -> &mut Self {
        self.file_id = file_id;
        self
    }

    pub fn set_signal_to_noise(&mut self, sn: Option<u8>) -> &mut Self {
        self.signal_to_noise = sn;
        self
    }

    /// Here be dragons -
    /// Seriously, this kinda sucks because it's a giant imperative, stateful loop.
    /// But I also don't want to spend any more time working on an mzML parser...
    pub async fn parse<B: AsyncBufRead + Unpin>(
        &self,
        b: B,
    ) -> Result<Vec<RawSpectrum>, MzMLError> {
        let mut reader = Reader::from_reader(b);
        let mut buf = Vec::new();

        let mut state = None;
        let mut compression = false;
        let mut output_buffer = Vec::with_capacity(4096);
        let mut binary_dtype = Dtype::F64;
        let mut binary_array = None;

        let mut spectrum = RawSpectrum::default_with_file_id(self.file_id);
        let mut precursor = Precursor::default();
        let mut iso_window_lo: Option<f32> = None;
        let mut iso_window_hi: Option<f32> = None;
        let mut spectra = Vec::new();

        let mut noise_array = Vec::new();
        let mut referenceable_params: HashMap<String, Vec<CvParam>> = HashMap::new();
        let mut current_referenceable_group: Option<String> = None;

        macro_rules! extract {
            ($ev:expr, $key:expr) => {
                $ev.try_get_attribute($key)?
                    .ok_or(MzMLError::Malformed)?
                    .value
            };
        }

        macro_rules! apply_cv_param {
            ($param:expr) => {{
                let param = $param;
                match state {
                    Some(State::BinaryDataArray) => match param.accession.as_slice() {
                        ZLIB_COMPRESSION => compression = true,
                        NO_COMPRESSION => compression = false,
                        FLOAT_64 => binary_dtype = Dtype::F64,
                        FLOAT_32 => binary_dtype = Dtype::F32,
                        INTENSITY_ARRAY => binary_array = Some(BinaryKind::Intensity),
                        MZ_ARRAY => binary_array = Some(BinaryKind::Mz),
                        NOISE_ARRAY => binary_array = Some(BinaryKind::Noise),
                        _ => {}
                    },
                    Some(State::Spectrum) => match param.accession.as_slice() {
                        MS_LEVEL => {
                            let level = param.parse_value()?;
                            if let Some(filter) = self.ms_level {
                                if level != filter {
                                    spectrum = RawSpectrum::default_with_file_id(self.file_id);
                                    state = None;
                                }
                            }
                            spectrum.ms_level = level;
                        }
                        PROFILE => spectrum.representation = Representation::Profile,
                        CENTROID => spectrum.representation = Representation::Centroid,
                        TOTAL_ION_CURRENT => {
                            let value = param.parse_value()?;
                            if value == 0.0 {
                                spectrum = RawSpectrum::default_with_file_id(self.file_id);
                                state = None;
                            } else {
                                spectrum.total_ion_current = value;
                            }
                        }
                        _ => {}
                    },
                    Some(State::Precursor) => match param.accession.as_slice() {
                        ISO_WINDOW_TARGET => {
                            if precursor.mz == 0.0 {
                                precursor.mz = param.parse_value()?;
                            }
                        }
                        ISO_WINDOW_LOWER => iso_window_lo = Some(param.parse_value()?),
                        ISO_WINDOW_UPPER => iso_window_hi = Some(param.parse_value()?),
                        _ => {}
                    },
                    Some(State::SelectedIon) => match param.accession.as_slice() {
                        SELECTED_ION_CHARGE => precursor.charge = Some(param.parse_value()?),
                        SELECTED_ION_MZ => {
                            let value = param.parse_value()?;
                            if value != 0.0 {
                                precursor.mz = value;
                            }
                        }
                        SELECTED_ION_INT => precursor.intensity = Some(param.parse_value()?),
                        INVERSE_ION_MOBILITY => {
                            precursor.inverse_ion_mobility = Some(param.parse_value()?);
                        }
                        _ => {}
                    },
                    Some(State::Scan) => match param.accession.as_slice() {
                        SCAN_START_TIME => {
                            let scan_start_time: f32 = param.parse_value()?;
                            spectrum.scan_start_time = match param.unit_accession.as_deref() {
                                Some(UNIT_SECONDS) => scan_start_time / 60.0,
                                Some(UNIT_MINUTES) => scan_start_time,
                                _ => return Err(MzMLError::Malformed),
                            };
                        }
                        ION_INJECTION_TIME => {
                            spectrum.ion_injection_time = param.parse_value()?;
                        }
                        INVERSE_ION_MOBILITY => {
                            precursor.inverse_ion_mobility = Some(param.parse_value()?);
                        }
                        _ => {}
                    },
                    _ => {}
                }
            }};
        }

        loop {
            match reader.read_event_into_async(&mut buf).await {
                Ok(Event::Start(ref ev)) => {
                    // State transition into child tag
                    state = match (ev.name().into_inner(), state) {
                        (b"spectrum", _) => Some(State::Spectrum),
                        (b"scan", Some(State::Spectrum)) => Some(State::Scan),
                        (b"binaryDataArray", Some(State::Spectrum)) => Some(State::BinaryDataArray),
                        (b"binary", Some(State::BinaryDataArray)) => Some(State::Binary),
                        (b"precursor", Some(State::Spectrum)) => Some(State::Precursor),
                        (b"selectedIon", Some(State::Precursor)) => Some(State::SelectedIon),
                        _ => state,
                    };
                    match ev.name().into_inner() {
                        b"spectrum" => {
                            let id = extract!(ev, b"id");
                            let id = std::str::from_utf8(&id)?;
                            spectrum.id = id.to_string();
                        }
                        b"precursor" => {
                            // Not all precursor fields have a spectrumRef
                            if let Some(scan) = ev.try_get_attribute(b"spectrumRef")? {
                                let scan = std::str::from_utf8(&scan.value)?;
                                precursor.spectrum_ref = Some(scan.to_string())
                            }
                        }
                        b"referenceableParamGroup" => {
                            let id = extract!(ev, b"id");
                            let id = std::str::from_utf8(&id)?.to_owned();
                            referenceable_params.entry(id.clone()).or_default();
                            current_referenceable_group = Some(id);
                        }
                        _ => {}
                    }
                }
                Ok(Event::Empty(ref ev)) => {
                    if ev.name().into_inner() == b"cvParam" {
                        let param = CvParam::from_event(ev)?;
                        if let Some(group) = current_referenceable_group.as_ref() {
                            referenceable_params
                                .get_mut(group)
                                .ok_or(MzMLError::Malformed)?
                                .push(param);
                        } else {
                            apply_cv_param!(&param);
                        }
                    } else if ev.name().into_inner() == b"referenceableParamGroupRef" {
                        let id = extract!(ev, b"ref");
                        let id = std::str::from_utf8(&id)?;
                        let params = referenceable_params
                            .get(id)
                            .ok_or_else(|| MzMLError::UnknownReferenceableParamGroup(id.into()))?
                            .clone();
                        for param in &params {
                            apply_cv_param!(param);
                        }
                    }
                }
                Ok(Event::Text(text)) => {
                    if let Some(State::Binary) = state {
                        if let Some(filter) = self.ms_level {
                            if spectrum.ms_level != filter {
                                continue;
                            }
                        }
                        let raw = text.unescape()?;
                        // There are occasionally empty binary data arrays, or unknown CVs
                        if raw.is_empty() || binary_array.is_none() {
                            continue;
                        }
                        let decoded = base64::decode(raw.as_bytes())?;
                        let bytes = match compression {
                            false => &decoded,
                            true => {
                                let mut r = ZlibDecoder::new(decoded.as_slice());
                                let n = r.read_to_end(&mut output_buffer).await?;
                                &output_buffer[..n]
                            }
                        };

                        let array = match binary_dtype {
                            Dtype::F32 => {
                                let mut buf: [u8; 4] = [0; 4];
                                bytes
                                    .chunks(4)
                                    .filter(|chunk| chunk.len() == 4)
                                    .map(|chunk| {
                                        buf.copy_from_slice(chunk);
                                        f32::from_le_bytes(buf)
                                    })
                                    .collect::<Vec<f32>>()
                            }
                            Dtype::F64 => {
                                let mut buf: [u8; 8] = [0; 8];
                                bytes
                                    .chunks(8)
                                    .map(|chunk| {
                                        buf.copy_from_slice(chunk);
                                        f64::from_le_bytes(buf) as f32
                                    })
                                    .collect::<Vec<f32>>()
                            }
                        };
                        output_buffer.clear();

                        match binary_array {
                            Some(BinaryKind::Intensity) => {
                                spectrum.intensity = array;
                            }
                            Some(BinaryKind::Mz) => {
                                spectrum.mz = array;
                            }
                            Some(BinaryKind::Noise) => {
                                noise_array = array;
                            }
                            None => {}
                        }

                        binary_array = None;
                    }
                }
                Ok(Event::End(ev)) => {
                    state = match (state, ev.name().into_inner()) {
                        (Some(State::Binary), b"binary") => Some(State::BinaryDataArray),
                        (Some(State::BinaryDataArray), b"binaryDataArray") => Some(State::Spectrum),
                        (Some(State::SelectedIon), b"selectedIon") => Some(State::Precursor),
                        (Some(State::Precursor), b"precursor") => {
                            if precursor.mz != 0.0 {
                                precursor.isolation_window = match (iso_window_lo, iso_window_hi) {
                                    (Some(lo), Some(hi)) => Some(Tolerance::Da(-lo, hi)),
                                    _ => None,
                                };
                                spectrum.precursors.push(precursor);
                                precursor = Precursor::default();
                            }
                            Some(State::Spectrum)
                        }
                        (Some(State::Scan), b"scan") => Some(State::Spectrum),
                        (_, b"referenceableParamGroup") => {
                            current_referenceable_group = None;
                            state
                        }
                        (_, b"spectrum") => {
                            let allow = self
                                .ms_level
                                .as_ref()
                                .map(|&level| level == spectrum.ms_level)
                                .unwrap_or(true);

                            match (allow, self.signal_to_noise) {
                                (true, Some(level))
                                    if level == spectrum.ms_level && !noise_array.is_empty() =>
                                {
                                    spectrum
                                        .intensity
                                        .iter_mut()
                                        .zip(noise_array.iter())
                                        .for_each(|(int, noise)| *int /= noise);
                                    noise_array.clear();
                                    spectra.push(spectrum);
                                }
                                (true, _) => {
                                    spectra.push(spectrum);
                                }
                                (false, _) => {}
                            }
                            spectrum = RawSpectrum::default_with_file_id(self.file_id);
                            None
                        }
                        _ => state,
                    };
                }
                Ok(Event::Eof) => break,
                Ok(_) => {}
                Err(err) => {
                    log::error!("unhandled XML error while parsing mzML: {}", err)
                }
            }
            buf.clear();
        }
        Ok(spectra)
    }
}

#[derive(thiserror::Error, Debug)]
pub enum MzMLError {
    #[error("malformed MzML")]
    Malformed,
    #[error("unsupported cvParam {0}")]
    UnsupportedCV(String),
    #[error("unknown referenceableParamGroup `{0}`")]
    UnknownReferenceableParamGroup(String),
    #[error("XML parsing error: {0}")]
    XMLError(#[from] quick_xml::Error),
    #[error("io error: {0}")]
    IOError(#[from] std::io::Error),
    #[error("utf8 error: {0}")]
    Utf8Error(#[from] std::str::Utf8Error),
    #[error("error parsing float: {0}")]
    FloatError(#[from] std::num::ParseFloatError),
    #[error("error parsing int: {0}")]
    IntError(#[from] std::num::ParseIntError),
    #[error("error decoding base64: {0}")]
    Base64Error(#[from] base64::DecodeError),
}

#[cfg(test)]
#[path = "../tests/unit/mzml.rs"]
mod test;
