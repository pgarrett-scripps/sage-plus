//! Bruker timsTOF scan-to-1/K0 conversion.
//!
//! timsrust interpolates scan numbers between the acquisition range limits
//! `OneOverK0AcqRangeLower` and `OneOverK0AcqRangeUpper`, linearly in
//! `sqrt(1/K0)`. Bruker
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
}

#[derive(Debug, thiserror::Error)]
pub enum MobilityCalibrationError {
    #[error("failed to read ion mobility calibration from {path}: {source}")]
    Sql {
        path: PathBuf,
        #[source]
        source: rusqlite::Error,
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
    /// Read the calibration of a `.d` directory or an `analysis.tdf` file.
    pub fn from_path(path: impl AsRef<Path>) -> Result<Self, MobilityCalibrationError> {
        let path = path.as_ref();
        let tdf = if path.is_dir() {
            path.join("analysis.tdf")
        } else {
            path.to_path_buf()
        };
        let sql = |source| MobilityCalibrationError::Sql {
            path: tdf.clone(),
            source,
        };
        let connection =
            rusqlite::Connection::open_with_flags(&tdf, rusqlite::OpenFlags::SQLITE_OPEN_READ_ONLY)
                .map_err(sql)?;

        let mut ids = HashMap::new();
        let mut models = Vec::new();
        let mut statement = connection
            .prepare(
                "SELECT Id, ModelType, C0, C1, C2, C3, C4, C5, C6, C7, C8, C9 \
                 FROM TimsCalibration ORDER BY Id",
            )
            .map_err(sql)?;
        let rows = statement
            .query_map([], |row| {
                let mut coefficients = [0.0; 10];
                for (index, value) in coefficients.iter_mut().enumerate() {
                    *value = row.get(index + 2)?;
                }
                Ok((row.get::<_, i64>(0)?, row.get::<_, i64>(1)?, coefficients))
            })
            .map_err(sql)?;
        for row in rows {
            let (id, model, coefficients) = row.map_err(sql)?;
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
        let mut statement = connection
            .prepare("SELECT Id, TimsCalibration FROM Frames")
            .map_err(sql)?;
        let rows = statement
            .query_map([], |row| Ok((row.get::<_, i64>(0)?, row.get::<_, i64>(1)?)))
            .map_err(sql)?;
        for row in rows {
            let (frame, id) = row.map_err(sql)?;
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
    pub fn dda_precursor_scans(
        path: impl AsRef<Path>,
    ) -> Result<HashMap<usize, (usize, f64)>, MobilityCalibrationError> {
        let path = path.as_ref();
        let tdf = if path.is_dir() {
            path.join("analysis.tdf")
        } else {
            path.to_path_buf()
        };
        let sql = |source| MobilityCalibrationError::Sql {
            path: tdf.clone(),
            source,
        };
        let connection =
            rusqlite::Connection::open_with_flags(&tdf, rusqlite::OpenFlags::SQLITE_OPEN_READ_ONLY)
                .map_err(sql)?;
        let exists: bool = connection
            .query_row(
                "SELECT EXISTS(SELECT 1 FROM sqlite_master WHERE type = 'table' AND name = 'Precursors')",
                [],
                |row| row.get(0),
            )
            .map_err(sql)?;
        if !exists {
            return Ok(HashMap::new());
        }
        let mut statement = connection
            .prepare("SELECT Id, Parent, ScanNumber FROM Precursors WHERE ScanNumber IS NOT NULL")
            .map_err(sql)?;
        let rows = statement
            .query_map([], |row| {
                Ok((
                    row.get::<_, i64>(0)? as usize,
                    (row.get::<_, i64>(1)? as usize, row.get::<_, f64>(2)?),
                ))
            })
            .map_err(sql)?;
        rows.collect::<Result<_, _>>().map_err(sql)
    }
}

#[cfg(test)]
#[path = "../tests/unit/tims_mobility.rs"]
mod test;
