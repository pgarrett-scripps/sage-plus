//! Bruker timsTOF scan-to-1/K0 conversion.
//!
//! Uncalibrated conversions interpolate scan numbers between the acquisition
//! range limits `OneOverK0AcqRangeLower` and `OneOverK0AcqRangeUpper`. Bruker
//! software instead applies the per-frame calibration stored in the
//! `TimsCalibration` table of `analysis.tdf`. The two scales differ by up to
//! about 2%, so the calibrated scale is the default.
//!
//! ModelType 2 is the only supported model. For a scan number `s`,
//!
//! ```text
//! V    = C2 + (C2 - C3) / C1 * (C4 + C0 - s)
//! 1/K0 = V / (C6 * V + C7)
//! ```
//!
//! Outside `C8 <= V <= C9` the model continues as a straight line with the
//! slope at the limit. This reproduces `tims_scannum_to_oneoverk0` from the
//! Bruker SDK to within floating-point rounding, including fractional scans
//! and scans beyond both limits.

use sage_plus_tdf::{MetadataDatabase, MetadataRow, SqlValue};
use serde::{Deserialize, Serialize};
use std::collections::HashMap;
use std::path::{Path, PathBuf};

/// Which scan-to-1/K0 scale timsTOF mobilities are reported on.
#[derive(
    Deserialize, Serialize, Debug, Clone, Copy, Default, PartialEq, Eq, schemars::JsonSchema,
)]
#[serde(rename_all = "snake_case")]
pub enum BrukerMobilityScale {
    /// Per-frame `TimsCalibration` model, as used by Bruker software.
    #[default]
    Calibrated,
    /// timsrust interpolation between the acquisition range limits, as
    /// reported by Sage Plus Beta 6 and earlier.
    Linear,
}

impl BrukerMobilityScale {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::Calibrated => "calibrated",
            Self::Linear => "linear",
        }
    }

    /// The scale actually reported for the Bruker input at `path`. Inputs
    /// without an `analysis.tdf`, such as miniTDF `.ms2` directories, have no
    /// calibration table and always use the linear scale.
    pub fn effective_for(self, path: impl AsRef<Path>) -> Self {
        match self {
            Self::Calibrated if analysis_tdf(path).is_none() => Self::Linear,
            scale => scale,
        }
    }
}

/// The `analysis.tdf` of a Bruker input.
///
/// Sage accepts the `.d` directory or any path inside it, such as
/// `analysis.tdf_bin`, and walks up to the first directory containing both
/// `analysis.tdf` and `analysis.tdf_bin`. Returns `None` for inputs that are
/// not TDF acquisitions, such as miniTDF `.ms2` directories.
pub fn analysis_tdf(path: impl AsRef<Path>) -> Option<PathBuf> {
    path.as_ref()
        .ancestors()
        .filter(|directory| !directory.as_os_str().is_empty())
        .map(|directory| directory.join("analysis.tdf"))
        .find(|tdf| tdf.is_file() && tdf.with_file_name("analysis.tdf_bin").is_file())
}

/// [`analysis_tdf`], or the path it was expected at for error messages.
fn analysis_tdf_or_guess(path: &Path) -> PathBuf {
    analysis_tdf(path).unwrap_or_else(|| {
        if path.is_dir() {
            path.join("analysis.tdf")
        } else {
            path.to_path_buf()
        }
    })
}

#[derive(Debug, thiserror::Error)]
pub enum MobilityCalibrationError {
    #[error("failed to read ion mobility calibration from {path}: {source}")]
    Tdf {
        path: PathBuf,
        #[source]
        source: sage_plus_tdf::Error,
    },
    #[error(
        "{path}: TimsCalibration {id} uses unsupported ModelType {model}; set \
         `bruker_config.ion_mobility_scale` to \"linear\" to report the uncalibrated scale"
    )]
    UnsupportedModel { path: PathBuf, id: i64, model: i64 },
    #[error("{path}: frame {frame} references missing TimsCalibration {id}")]
    MissingCalibration { path: PathBuf, frame: i64, id: i64 },
    #[error("{path}: no TimsCalibration rows")]
    Empty { path: PathBuf },
    #[error("{path}: GlobalMetadata {key} is not a number: {value:?}")]
    InvalidMetadata {
        path: PathBuf,
        key: String,
        value: String,
    },
    #[error("{path}: {table}.{column} is missing or has the wrong type")]
    InvalidColumn {
        path: PathBuf,
        table: &'static str,
        column: &'static str,
    },
}

/// The metadata tables of one `analysis.tdf`, read without the SQLite C library.
struct Tables {
    path: PathBuf,
    database: MetadataDatabase,
}

impl Tables {
    fn open(path: &Path) -> Result<Self, MobilityCalibrationError> {
        let path = analysis_tdf_or_guess(path);
        let database =
            MetadataDatabase::open(&path).map_err(|source| MobilityCalibrationError::Tdf {
                path: path.clone(),
                source,
            })?;
        Ok(Self { path, database })
    }

    fn error(&self, source: sage_plus_tdf::Error) -> MobilityCalibrationError {
        MobilityCalibrationError::Tdf {
            path: self.path.clone(),
            source,
        }
    }

    fn table(&self, name: &str) -> Result<Vec<MetadataRow>, MobilityCalibrationError> {
        self.database.table(name).map_err(|error| self.error(error))
    }

    fn has_table(&self, name: &str) -> Result<bool, MobilityCalibrationError> {
        self.database
            .has_table(name)
            .map_err(|error| self.error(error))
    }

    fn column<T>(
        &self,
        row: &MetadataRow,
        table: &'static str,
        column: &'static str,
        read: impl Fn(&SqlValue) -> Option<T>,
    ) -> Result<T, MobilityCalibrationError> {
        row.get(column)
            .and_then(read)
            .ok_or_else(|| MobilityCalibrationError::InvalidColumn {
                path: self.path.clone(),
                table,
                column,
            })
    }

    fn metadata(&self, key: &str) -> Result<f64, MobilityCalibrationError> {
        let rows = self.table("GlobalMetadata")?;
        let value = rows
            .iter()
            .find(|row| row.get("Key").and_then(SqlValue::as_str) == Some(key))
            .and_then(|row| row.get("Value"))
            .ok_or_else(|| MobilityCalibrationError::InvalidMetadata {
                path: self.path.clone(),
                key: key.to_string(),
                value: String::new(),
            })?;
        match value {
            SqlValue::Text(text) => text.parse().ok(),
            value => value.as_f64(),
        }
        .ok_or_else(|| MobilityCalibrationError::InvalidMetadata {
            path: self.path.clone(),
            key: key.to_string(),
            value: format!("{value:?}"),
        })
    }
}

/// The uncalibrated scale reported by Sage Plus Beta 6 and earlier.
///
/// Scan numbers are interpolated between `OneOverK0AcqRangeUpper` (scan 0) and
/// `OneOverK0AcqRangeLower` (the largest `Frames.NumScans`), linearly in
/// `sqrt(1/K0)`. This is the conversion of timsrust 0.6.5, kept here with the
/// same floating-point operations so the scale does not change with the reader.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct LinearMobilityScale {
    intercept: f64,
    slope: f64,
}

impl LinearMobilityScale {
    pub fn new(lower: f64, upper: f64, scans: u64) -> Self {
        let intercept = upper.sqrt();
        Self {
            intercept,
            slope: (lower.sqrt() - intercept) / scans as u32 as f64,
        }
    }

    /// Read the acquisition range of a `.d` directory or a file inside it.
    pub fn from_path(path: impl AsRef<Path>) -> Result<Self, MobilityCalibrationError> {
        let tables = Tables::open(path.as_ref())?;
        let lower = tables.metadata("OneOverK0AcqRangeLower")?;
        let upper = tables.metadata("OneOverK0AcqRangeUpper")?;
        let mut scans = 0i64;
        for row in tables.table("Frames")? {
            scans = scans.max(tables.column(&row, "Frames", "NumScans", SqlValue::as_i64)?);
        }
        Ok(Self::new(lower, upper, scans as u64))
    }

    pub fn one_over_k0(&self, scan: u32) -> f64 {
        let value = self.intercept + self.slope * f64::from(scan);
        value * value
    }
}

/// One ModelType 2 `TimsCalibration` row, coefficients `C0` to `C9`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct MobilityModel {
    coefficients: [f64; 10],
}

impl MobilityModel {
    pub fn new(coefficients: [f64; 10]) -> Self {
        Self { coefficients }
    }

    /// Convert a possibly fractional scan number to 1/K0.
    pub fn one_over_k0(&self, scan: f64) -> f64 {
        let [c0, c1, c2, c3, c4, _, c6, c7, c8, c9] = self.coefficients;
        let v = c2 + (c2 - c3) / c1 * (c4 + c0 - scan);
        let model = |v: f64| v / (c6 * v + c7);
        let slope = |v: f64| c7 / (c6 * v + c7).powi(2);
        if v < c8 {
            model(c8) + slope(c8) * (v - c8)
        } else if v > c9 {
            model(c9) + slope(c9) * (v - c9)
        } else {
            model(v)
        }
    }
}

/// Calibration models of one acquisition and the model used by each frame.
#[derive(Debug, Clone)]
pub struct MobilityCalibration {
    models: Vec<MobilityModel>,
    frame_models: HashMap<usize, usize>,
    /// Model used by the most frames, for values without a frame.
    dominant: usize,
}

impl MobilityCalibration {
    /// Read the calibration of a `.d` directory or a file inside it.
    pub fn from_path(path: impl AsRef<Path>) -> Result<Self, MobilityCalibrationError> {
        let tables = Tables::open(path.as_ref())?;
        let tdf = tables.path.clone();

        let mut rows = Vec::new();
        for row in tables.table("TimsCalibration")? {
            let id = tables.column(&row, "TimsCalibration", "Id", SqlValue::as_i64)?;
            let model = tables.column(&row, "TimsCalibration", "ModelType", SqlValue::as_i64)?;
            const NAMES: [&str; 10] = ["C0", "C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8", "C9"];
            let mut coefficients = [0.0; 10];
            for (value, name) in coefficients.iter_mut().zip(NAMES) {
                *value = tables.column(&row, "TimsCalibration", name, SqlValue::as_f64)?;
            }
            rows.push((id, model, coefficients));
        }
        rows.sort_by_key(|&(id, _, _)| id);
        let mut ids = HashMap::new();
        let mut models = Vec::new();
        for (id, model, coefficients) in rows {
            if model != 2 {
                return Err(MobilityCalibrationError::UnsupportedModel {
                    path: tdf.clone(),
                    id,
                    model,
                });
            }
            ids.insert(id, models.len());
            models.push(MobilityModel::new(coefficients));
        }
        if models.is_empty() {
            return Err(MobilityCalibrationError::Empty { path: tdf });
        }

        let mut frame_models = HashMap::new();
        let mut counts = vec![0usize; models.len()];
        for row in tables.table("Frames")? {
            let frame = tables.column(&row, "Frames", "Id", SqlValue::as_i64)?;
            let id = tables.column(&row, "Frames", "TimsCalibration", SqlValue::as_i64)?;
            let model =
                *ids.get(&id)
                    .ok_or_else(|| MobilityCalibrationError::MissingCalibration {
                        path: tdf.clone(),
                        frame,
                        id,
                    })?;
            frame_models.insert(frame as usize, model);
            counts[model] += 1;
        }
        let dominant = (0..models.len())
            .max_by_key(|&model| (counts[model], std::cmp::Reverse(model)))
            .unwrap_or(0);
        Ok(Self {
            models,
            frame_models,
            dominant,
        })
    }

    /// The model of `frame` (the `Frames.Id` value). Frames without a
    /// calibration row fall back to the model used by the most frames.
    pub fn frame(&self, frame: usize) -> &MobilityModel {
        let model = self
            .frame_models
            .get(&frame)
            .copied()
            .unwrap_or(self.dominant);
        &self.models[model]
    }

    /// The model used by the most frames.
    pub fn dominant(&self) -> &MobilityModel {
        &self.models[self.dominant]
    }

    /// Precursor `Id` to `(Parent frame, fractional ScanNumber)` for DDA runs.
    /// Precursors without a `ScanNumber` or `Parent` are left out.
    pub fn dda_precursor_scans(
        path: impl AsRef<Path>,
    ) -> Result<HashMap<usize, (usize, f64)>, MobilityCalibrationError> {
        let tables = Tables::open(path.as_ref())?;
        if !tables.has_table("Precursors")? {
            return Ok(HashMap::new());
        }
        let mut scans = HashMap::new();
        for row in tables.table("Precursors")? {
            let (Some(scan), Some(parent)) = (
                row.get("ScanNumber").and_then(SqlValue::as_f64),
                row.get("Parent").and_then(SqlValue::as_i64),
            ) else {
                continue;
            };
            let id = tables.column(&row, "Precursors", "Id", SqlValue::as_i64)?;
            scans.insert(id as usize, (parent as usize, scan));
        }
        Ok(scans)
    }
}

#[cfg(test)]
#[path = "../tests/unit/tims_mobility.rs"]
mod test;
