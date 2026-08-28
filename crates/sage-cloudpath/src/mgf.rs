use regex::Regex;
use sage_core::mass::Tolerance;
use sage_core::spectrum::RawSpectrum;
use sage_core::spectrum::{Precursor, Representation};

#[derive(Clone)]
pub struct DefaultParams {
    is_query_start: bool,
    regex_for_charge: Regex,
    file_id: usize,
    tol: Option<f32>,
    tol_unit: Option<String>,
    charge_array: Option<Vec<u8>>,
}

impl Default for DefaultParams {
    fn default() -> Self {
        Self {
            is_query_start: false,
            file_id: 0,
            regex_for_charge: Regex::new(r"(\d)\+?").unwrap(),
            tol: None,
            tol_unit: None,
            charge_array: None,
        }
    }
}

impl DefaultParams {
    pub fn default_with_file_id(file_id: usize) -> Self {
        Self {
            file_id,
            ..Default::default()
        }
    }
}

#[derive(Default, Clone)]
pub struct QueryData {
    default_params: DefaultParams,
    spectra: Vec<RawSpectrum>,
    in_spectrum: bool,

    id: String,
    precursors: Vec<Precursor>,
    precursor_tol: Option<f32>,
    precursor_tol_unit: Option<String>,
    precursor_charge_array: Option<Vec<u8>>,
    rt_in_minutes: Option<f32>,
    inverse_ion_mobility: Option<f32>,
    ion_mz_array: Vec<f32>,
    ion_intensity_array: Vec<f32>,
}

impl QueryData {
    pub fn default_with_params(default_params: DefaultParams) -> Self {
        Self {
            default_params,
            in_spectrum: true,
            ..Default::default()
        }
    }
    pub fn init(&mut self) {
        self.in_spectrum = false;
        self.id = String::default();
        self.precursors = Vec::new();
        self.precursor_tol = self.default_params.tol;
        self.precursor_tol_unit = self.default_params.tol_unit.clone();
        self.precursor_charge_array = self.default_params.charge_array.clone();
        self.rt_in_minutes = None;
        self.inverse_ion_mobility = None;
        self.ion_mz_array = Vec::new();
        self.ion_intensity_array = Vec::new();
    }

    pub fn get_isolation_window(&mut self) -> Option<Tolerance> {
        if let Some(tol_value) = self.precursor_tol {
            if let Some(unit_str) = self.precursor_tol_unit.as_deref() {
                match unit_str {
                    "Da" => return Some(Tolerance::Da(-tol_value.abs(), tol_value.abs())),
                    "ppm" => return Some(Tolerance::Ppm(-tol_value.abs(), tol_value.abs())),
                    _ => return None,
                }
            }
        }
        None
    }

    //precursors with isolation window and charge
    pub fn get_precursors_with_charge(&mut self) -> Vec<Precursor> {
        let mut new_precursors = Vec::new();
        let isolation_window = self.get_isolation_window();

        for precursor in &mut self.precursors {
            precursor.isolation_window = isolation_window;
            precursor.inverse_ion_mobility = self.inverse_ion_mobility;

            if let Some(charge_array) = &self.precursor_charge_array {
                for &charge in charge_array.iter() {
                    let mut precursor_with_charge = precursor.clone();
                    precursor_with_charge.charge = Some(charge);
                    new_precursors.push(precursor_with_charge);
                }
            } else {
                new_precursors.push(precursor.clone());
            }
        }
        new_precursors
    }

    pub fn default_spectrum(&self, file_id: usize) -> RawSpectrum {
        RawSpectrum {
            file_id,
            ms_level: 2,
            representation: Representation::Centroid,
            ..Default::default()
        }
    }

    pub fn check_spectrum(&self, spectrum: &RawSpectrum) -> Result<bool, MgfError> {
        if spectrum.id.is_empty() {
            return Err(MgfError::Malformed {
                message: "spectrum is missing TITLE",
            });
        }
        if spectrum.precursors.is_empty() {
            return Err(MgfError::Malformed {
                message: "spectrum is missing PEPMASS",
            });
        }
        if spectrum
            .precursors
            .iter()
            .any(|precursor| !precursor.mz.is_finite() || precursor.mz <= 0.0)
        {
            return Err(MgfError::Malformed {
                message: "spectrum contains an invalid precursor mass",
            });
        }
        if spectrum.mz.is_empty() {
            return Err(MgfError::Malformed {
                message: "spectrum contains no peaks",
            });
        }
        if spectrum.mz.len() != spectrum.intensity.len() {
            return Err(MgfError::Malformed {
                message: "peak mass and intensity arrays have different lengths",
            });
        }
        if spectrum
            .mz
            .iter()
            .chain(&spectrum.intensity)
            .any(|value| !value.is_finite())
        {
            return Err(MgfError::Malformed {
                message: "spectrum contains a nonfinite peak value",
            });
        }
        Ok(true)
    }
}

pub struct DefaultParser;

type DefaultParserFn = fn(&str, &mut DefaultParams) -> Result<bool, MgfError>;

impl DefaultParser {
    pub fn get_parsers(&self) -> Vec<DefaultParserFn> {
        vec![
            Self::parse_begin,
            Self::parse_tol,
            Self::parse_tol_unit,
            Self::parse_charge,
        ]
    }
    pub fn parse_begin(line: &str, default_params: &mut DefaultParams) -> Result<bool, MgfError> {
        if line.starts_with("BEGIN IONS") {
            default_params.is_query_start = true;
            return Ok(true);
        }
        Ok(false)
    }
    pub fn parse_tol(line: &str, default_params: &mut DefaultParams) -> Result<bool, MgfError> {
        if let Some(tol_str) = line.strip_prefix("TOL=") {
            default_params.tol = Some(tol_str.parse().map_err(|_| MgfError::Malformed {
                message: "invalid TOL value",
            })?);
            return Ok(true);
        }
        Ok(false)
    }
    pub fn parse_tol_unit(
        line: &str,
        default_params: &mut DefaultParams,
    ) -> Result<bool, MgfError> {
        if let Some(tol_unit_str) = line.strip_prefix("TOLU=") {
            default_params.tol_unit = Some(tol_unit_str.to_string());
            return Ok(true);
        }
        Ok(false)
    }
    pub fn parse_charge(line: &str, default_params: &mut DefaultParams) -> Result<bool, MgfError> {
        let regex_for_charge = &default_params.regex_for_charge;

        if let Some(charge_str) = line.strip_prefix("CHARGE=") {
            let mut charge_array: Vec<u8> = Vec::new();
            for cap in regex_for_charge.captures_iter(charge_str) {
                if let Some(charge) = cap[0].chars().next().unwrap().to_digit(10) {
                    charge_array.push(charge as u8);
                }
            }
            default_params.charge_array = Some(charge_array);
            return Ok(true);
        }
        Ok(false)
    }
}

pub struct QueryParser;
type QueryParserFn = fn(&str, &mut QueryData) -> Result<bool, MgfError>;

impl QueryParser {
    pub fn get_parsers(&self) -> Vec<QueryParserFn> {
        vec![
            Self::parse_begin,
            Self::parse_mz,
            Self::parse_end,
            Self::parse_pepmass,
            Self::parse_title,
            Self::parse_charge,
            Self::parse_tol,
            Self::parse_tol_unit,
            Self::parse_rt,
            Self::parse_ion_mobility,
        ]
    }

    pub fn parse_begin(line: &str, query_data: &mut QueryData) -> Result<bool, MgfError> {
        if line.starts_with("BEGIN IONS") {
            if query_data.in_spectrum {
                return Err(MgfError::Malformed {
                    message: "nested BEGIN IONS marker",
                });
            }
            query_data.in_spectrum = true;
            return Ok(true);
        }
        Ok(false)
    }

    pub fn parse_pepmass(line: &str, query_data: &mut QueryData) -> Result<bool, MgfError> {
        if let Some(pepmass_str) = line.strip_prefix("PEPMASS=") {
            let mut precursor = Precursor::default();
            let mut pepmass = pepmass_str.split_ascii_whitespace();
            if let Some(mz_str) = pepmass.next() {
                match mz_str.parse::<f32>() {
                    Ok(mz) => precursor.mz = mz,
                    Err(_) => {
                        return Err(MgfError::Malformed {
                            message: "invalid PEPMASS value",
                        })
                    }
                }
            }
            if let Some(intensity_str) = pepmass.next() {
                precursor.intensity =
                    Some(intensity_str.parse().map_err(|_| MgfError::Malformed {
                        message: "invalid PEPMASS intensity",
                    })?);
            }
            query_data.precursors.push(precursor);
            return Ok(true);
        }
        Ok(false)
    }

    pub fn parse_charge(line: &str, query_data: &mut QueryData) -> Result<bool, MgfError> {
        let regex_for_charge = &query_data.default_params.regex_for_charge;

        if let Some(charge_str) = line.strip_prefix("CHARGE=") {
            let mut charge_array = Vec::new();
            for cap in regex_for_charge.captures_iter(charge_str) {
                if let Some(charge) = cap[0].chars().next().unwrap().to_digit(10) {
                    charge_array.push(charge as u8);
                }
            }
            query_data.precursor_charge_array = Some(charge_array);
            return Ok(true);
        }
        Ok(false)
    }

    pub fn parse_rt(line: &str, query_data: &mut QueryData) -> Result<bool, MgfError> {
        if let Some(rt_str) = line.strip_prefix("RTINSECONDS=") {
            let rt_in_seconds = rt_str.parse::<f32>().map_err(|_| MgfError::Malformed {
                message: "invalid RTINSECONDS value",
            })?;
            query_data.rt_in_minutes = Some(rt_in_seconds / 60.0);
            return Ok(true);
        }
        Ok(false)
    }

    pub fn parse_title(line: &str, query_data: &mut QueryData) -> Result<bool, MgfError> {
        if let Some(id_str) = line.strip_prefix("TITLE=") {
            query_data.id = id_str.to_string();
            if let Some(mobility) = parse_inverse_ion_mobility(id_str) {
                query_data.inverse_ion_mobility = Some(mobility);
            }
            return Ok(true);
        }
        Ok(false)
    }

    pub fn parse_ion_mobility(line: &str, query_data: &mut QueryData) -> Result<bool, MgfError> {
        if let Some(value) = line.strip_prefix("INVERSE_REDUCED_ION_MOBILITY=") {
            query_data.inverse_ion_mobility =
                Some(value.parse::<f32>().map_err(|_| MgfError::Malformed {
                    message: "invalid inverse ion mobility value",
                })?);
            return Ok(true);
        }
        Ok(false)
    }

    pub fn parse_tol(line: &str, query_data: &mut QueryData) -> Result<bool, MgfError> {
        if let Some(tol_str) = line.strip_prefix("TOL=") {
            query_data.precursor_tol = Some(tol_str.parse().map_err(|_| MgfError::Malformed {
                message: "invalid TOL value",
            })?);
            return Ok(true);
        }
        Ok(false)
    }

    pub fn parse_tol_unit(line: &str, query_data: &mut QueryData) -> Result<bool, MgfError> {
        if let Some(tol_unit_str) = line.strip_prefix("TOLU=") {
            query_data.precursor_tol_unit = Some(tol_unit_str.to_string());
            return Ok(true);
        }
        Ok(false)
    }

    pub fn parse_mz(line: &str, query_data: &mut QueryData) -> Result<bool, MgfError> {
        if line.chars().nth(0).unwrap_or_default().is_numeric() {
            let mut mz_intensity = line.split_ascii_whitespace();
            if let Some(mz_str) = mz_intensity.next() {
                match mz_str.parse::<f32>() {
                    Ok(mz) => query_data.ion_mz_array.push(mz),
                    Err(_) => {
                        return Err(MgfError::Malformed {
                            message: "invalid peak mass",
                        })
                    }
                }
            }
            if let Some(intensity_str) = mz_intensity.next() {
                query_data
                    .ion_intensity_array
                    .push(intensity_str.parse().map_err(|_| MgfError::Malformed {
                        message: "invalid peak intensity",
                    })?);
            } else {
                query_data.ion_intensity_array.push(1.0)
            }
            return Ok(true);
        }
        Ok(false)
    }

    pub fn parse_end(line: &str, query_data: &mut QueryData) -> Result<bool, MgfError> {
        if line.starts_with("END IONS") {
            let mut spectrum = query_data.default_spectrum(query_data.default_params.file_id);

            spectrum.id = query_data.id.to_string();
            spectrum.precursors = query_data.get_precursors_with_charge();
            spectrum.scan_start_time = query_data.rt_in_minutes.unwrap_or_default();
            spectrum.total_ion_current = query_data.ion_intensity_array.iter().sum();
            spectrum.mz = std::mem::take(&mut query_data.ion_mz_array);
            spectrum.intensity = std::mem::take(&mut query_data.ion_intensity_array);

            query_data.check_spectrum(&spectrum)?;
            query_data.spectra.push(spectrum);
            query_data.init();

            return Ok(true);
        }
        Ok(false)
    }
}

fn parse_inverse_ion_mobility(title: &str) -> Option<f32> {
    let value = title.split_once("1/K0=")?.1;
    value
        .split(|character: char| character.is_ascii_whitespace() || character == ',')
        .next()?
        .parse()
        .ok()
}

pub struct MgfReader {
    file_id: usize,
}

impl MgfReader {
    pub fn with_file_id(file_id: usize) -> Self {
        Self { file_id }
    }

    pub fn parse(&self, contents: String) -> Result<Vec<RawSpectrum>, MgfError> {
        let default_parsers = DefaultParser.get_parsers();
        let query_parsers = QueryParser.get_parsers();

        let mut default_params = DefaultParams::default_with_file_id(self.file_id);
        let mut lines = contents.as_str().lines().enumerate();

        // embedded parameters
        while !default_params.is_query_start {
            let (line_index, line) = lines.next().ok_or(MgfError::MissingBeginIons)?;
            let line = line.trim();
            for parser in &default_parsers {
                match parser(line, &mut default_params) {
                    Ok(true) => break,
                    Ok(false) => continue,
                    Err(err) => return Err(err.at_line(line_index + 1)),
                }
            }
        }

        let mut query_data = QueryData::default_with_params(default_params);

        // query
        for (line_index, line) in lines {
            if line.is_empty() {
                continue;
            }
            let line = line.trim();
            if !query_data.in_spectrum && !line.starts_with("BEGIN IONS") {
                continue;
            }
            for parser in &query_parsers {
                match parser(line, &mut query_data) {
                    Ok(true) => break,
                    Ok(false) => {}
                    Err(err) => return Err(err.at_line(line_index + 1)),
                }
            }
        }
        if query_data.in_spectrum {
            return Err(MgfError::UnterminatedSpectrum);
        }
        Ok(query_data.spectra)
    }
}

#[derive(thiserror::Error, Debug)]
pub enum MgfError {
    #[error("MGF does not contain a BEGIN IONS marker")]
    MissingBeginIons,
    #[error("MGF spectrum is missing END IONS")]
    UnterminatedSpectrum,
    #[error("malformed MGF: {message}")]
    Malformed { message: &'static str },
    #[error("malformed MGF at line {line}: {message}")]
    MalformedLine { line: usize, message: String },
    #[error("unsupported cvParam {0}")]
    UnsupportedCV(String),
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

impl MgfError {
    fn at_line(self, line: usize) -> Self {
        match self {
            Self::Malformed { message } => Self::MalformedLine {
                line,
                message: message.to_string(),
            },
            error => error,
        }
    }
}

#[cfg(test)]
#[path = "../tests/unit/mgf.rs"]
mod test;
