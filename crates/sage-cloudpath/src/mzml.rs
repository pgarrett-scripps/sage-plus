use async_compression::tokio::bufread::ZlibDecoder;
use quick_xml::events::Event;
use quick_xml::Reader;
use sage_core::spectrum::{Precursor, Representation};
use sage_core::{mass::Tolerance, spectrum::RawSpectrum};
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

        macro_rules! extract {
            ($ev:expr, $key:expr) => {
                $ev.try_get_attribute($key)?
                    .ok_or(MzMLError::Malformed)?
                    .value
            };
        }

        macro_rules! extract_value {
            ($ev:expr) => {{
                let s = $ev
                    .try_get_attribute(b"value")?
                    .ok_or(MzMLError::Malformed)?
                    .value;
                std::str::from_utf8(&s)?.parse()?
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
                        _ => {}
                    }
                }
                Ok(Event::Empty(ref ev)) => match (state, ev.name().into_inner()) {
                    (Some(State::BinaryDataArray), b"cvParam") => {
                        let accession = extract!(ev, b"accession");
                        match accession.as_ref() {
                            ZLIB_COMPRESSION => compression = true,
                            NO_COMPRESSION => compression = false,
                            FLOAT_64 => binary_dtype = Dtype::F64,
                            FLOAT_32 => binary_dtype = Dtype::F32,
                            INTENSITY_ARRAY => binary_array = Some(BinaryKind::Intensity),
                            MZ_ARRAY => binary_array = Some(BinaryKind::Mz),
                            NOISE_ARRAY => binary_array = Some(BinaryKind::Noise),
                            _ => {
                                // Unknown CV - perhaps noise
                                binary_array = None;
                            }
                        }
                    }
                    (Some(State::Spectrum), b"cvParam") => {
                        let accession = extract!(ev, b"accession");
                        match accession.as_ref() {
                            MS_LEVEL => {
                                let level = extract_value!(ev);
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
                                let value = extract_value!(ev);
                                if value == 0.0 {
                                    // No ion current, break out of current state
                                    spectrum = RawSpectrum::default_with_file_id(self.file_id);
                                    state = None;
                                } else {
                                    spectrum.total_ion_current = value;
                                }
                            }
                            _ => {}
                        }
                    }
                    (Some(State::Precursor), b"cvParam") => {
                        let accession = extract!(ev, b"accession");
                        match accession.as_ref() {
                            ISO_WINDOW_TARGET => {
                                // use isolation window target for precursor m/z, e.g. to handle
                                // DIA setups where the mzML conversion software doesn't write
                                // a selection ion tag
                                if precursor.mz == 0.0 {
                                    precursor.mz = extract_value!(ev)
                                }
                            }
                            ISO_WINDOW_LOWER => iso_window_lo = Some(extract_value!(ev)),
                            ISO_WINDOW_UPPER => iso_window_hi = Some(extract_value!(ev)),
                            _ => {}
                        }
                    }
                    (Some(State::SelectedIon), b"cvParam") => {
                        let accession = extract!(ev, b"accession");
                        match accession.as_ref() {
                            SELECTED_ION_CHARGE => {
                                precursor.charge = Some(extract_value!(ev));
                            }
                            SELECTED_ION_MZ => {
                                let val = extract_value!(ev);
                                if val != 0.0 {
                                    precursor.mz = val;
                                }
                            }
                            SELECTED_ION_INT => {
                                precursor.intensity = Some(extract_value!(ev));
                            }
                            INVERSE_ION_MOBILITY => {
                                precursor.inverse_ion_mobility = Some(extract_value!(ev));
                            }
                            _ => {}
                        }
                    }
                    (Some(State::Scan), b"cvParam") => {
                        let accession = extract!(ev, b"accession");
                        match accession.as_ref() {
                            SCAN_START_TIME => {
                                let scan_start_time = extract_value!(ev);
                                let unit = extract!(ev, b"unitAccession");

                                spectrum.scan_start_time = match unit.as_ref() {
                                    UNIT_SECONDS => scan_start_time / 60.0,
                                    UNIT_MINUTES => scan_start_time,
                                    _ => return Err(MzMLError::Malformed),
                                };
                            }
                            ION_INJECTION_TIME => {
                                spectrum.ion_injection_time = extract_value!(ev);
                            }
                            INVERSE_ION_MOBILITY => {
                                precursor.inverse_ion_mobility = Some(extract_value!(ev));
                            }
                            _ => {}
                        }
                    }

                    _ => {}
                },
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
