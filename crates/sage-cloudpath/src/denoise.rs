//! Opt-in timsTOF MS1 denoising with [dnoise](https://github.com/pgarrett-scripps/dnoise).
//!
//! When `bruker_config.denoise.enabled` is set, every MS1 frame of a TDF input
//! passes through dnoise's per-frame pipeline (vertical ion-mobility streak
//! filter, horizontal halo filter and the acquisition's MS1 gate) before Sage
//! centroids it. The result matches searching a `.d` written by the dnoise CLI
//! with the same settings, without writing a denoised copy.
//!
//! Only MS1 frames are denoised, which is also the dnoise CLI default: MS/MS
//! spectra come from timsrust's spectrum reader, which merges frames per
//! precursor, and dnoise leaves MS/MS unchanged unless asked.
//!
//! The gates are built exactly as the dnoise CLI builds them (`src/writer.rs`
//! in dnoise 0.5.0): the ddaPASEF selection polygon from `GroupProperties`, the
//! diaPASEF isolation windows from `DiaFrameMsMsWindows`, both converted with the
//! run-level m/z calibration of `GlobalMetadata` and the run's `TimsCalibration`.

use crate::tims_mobility::BrukerMobilityScale;
use dnoise_core::{
    convert::{ConvertableDomain, Scan2ImConverter, Tof2MzConverter},
    mobility::TimsCalibrationModel,
    polygon::PolygonGate,
    windows::DiaMs1Box,
    Acquisition, Denoiser, DiaMs1WindowParams, FilterParams, FlatFrame, FrameMeta, HaloParams,
    Ms1PolygonParams, ScanToMobility, Stages,
};
use rusqlite::{Connection, OpenFlags, OptionalExtension};
use serde::{Deserialize, Serialize};
use std::collections::{BTreeSet, HashMap};
use std::path::{Path, PathBuf};

/// dnoise settings for timsTOF MS1 frames. Off by default; every other key
/// defaults to the dnoise 0.5.0 CLI default and uses the same name as
/// `dnoise.toml`.
#[derive(Deserialize, Serialize, Debug, Clone, Copy, PartialEq, schemars::JsonSchema)]
#[serde(default, deny_unknown_fields)]
pub struct BrukerDenoiseConfig {
    /// Denoise MS1 frames before centroiding. Needs an `analysis.tdf`.
    pub enabled: bool,
    /// Streak filter column half-width, in TOF indices.
    pub mz_half_width: u32,
    /// Minimum number of occupied mobility scans in a kept streak.
    pub min_feature_length: usize,
    /// Largest run of empty scans bridged inside a streak.
    pub max_internal_gap: usize,
    /// Per-scan summed-intensity floor for an occupied scan.
    pub min_window_intensity: u64,
    /// Total intensity a streak needs to be kept.
    pub min_feature_intensity: u64,
    /// Number of streak filter passes.
    pub iterations: usize,
    /// Run the horizontal halo filter after the streak filter.
    pub halo: bool,
    /// Drop a point below this fraction of its off-column neighbourhood maximum.
    pub halo_peak_fraction: f64,
    /// Halo reference box half-width, in TOF indices.
    pub halo_mz_idx_half_width: u32,
    /// Halo reference box half-width, in mobility scans.
    pub halo_scan_half_width: usize,
    /// ddaPASEF: drop MS1 signal outside the padded precursor selection polygon.
    pub ms1_polygon: bool,
    /// Keep a whole feature when any of its points is inside the polygon.
    pub ms1_polygon_overlap: bool,
    /// Polygon m/z padding, in Th.
    pub ms1_polygon_mz_pad: f64,
    /// Polygon mobility padding, in 1/K0.
    pub ms1_polygon_im_pad: f64,
    /// diaPASEF: drop MS1 signal outside every padded isolation window.
    pub dia_ms1_window: bool,
    /// Keep a whole feature when any of its points is inside a window.
    pub dia_ms1_overlap: bool,
    /// Isolation window m/z padding, in Th.
    pub dia_ms1_mz_pad: f64,
    /// Isolation window mobility padding, in 1/K0.
    pub dia_ms1_im_pad: f64,
    /// Scan-to-1/K0 scale the gates are placed on.
    pub mobility_scale: BrukerMobilityScale,
}

impl Default for BrukerDenoiseConfig {
    fn default() -> Self {
        let filter = FilterParams::default();
        let halo = HaloParams::default();
        let polygon = Ms1PolygonParams::default();
        let dia = DiaMs1WindowParams::default();
        Self {
            enabled: false,
            mz_half_width: filter.mz_half_width,
            min_feature_length: filter.min_feature_length,
            max_internal_gap: filter.max_internal_gap,
            min_window_intensity: filter.min_window_intensity,
            min_feature_intensity: filter.min_feature_intensity,
            iterations: filter.num_iterations,
            halo: true,
            halo_peak_fraction: halo.peak_fraction,
            halo_mz_idx_half_width: halo.mz_idx_half_width,
            halo_scan_half_width: halo.scan_half_width,
            ms1_polygon: true,
            ms1_polygon_overlap: polygon.overlap,
            ms1_polygon_mz_pad: polygon.mz_pad,
            ms1_polygon_im_pad: polygon.im_pad,
            dia_ms1_window: true,
            dia_ms1_overlap: dia.overlap,
            dia_ms1_mz_pad: dia.mz_pad,
            dia_ms1_im_pad: dia.im_pad,
            mobility_scale: BrukerMobilityScale::Calibrated,
        }
    }
}

#[derive(Debug, thiserror::Error)]
pub enum DenoiseError {
    #[error("invalid bruker_config.denoise: {0}")]
    Config(String),
    #[error("{path}: failed to read denoising metadata: {source}")]
    Sql {
        path: PathBuf,
        #[source]
        source: rusqlite::Error,
    },
    #[error("{path}: {message}")]
    Metadata { path: PathBuf, message: String },
    #[error("{path}: denoising failed: {source}")]
    Dnoise {
        path: PathBuf,
        #[source]
        source: dnoise_core::Error,
    },
}

/// The dnoise parameter structs, owned so a [`Denoiser`] can borrow them.
#[derive(Debug, Clone, Copy)]
pub struct DenoiseParams {
    filter: FilterParams,
    halo: Option<HaloParams>,
    polygon: Option<Ms1PolygonParams>,
    dia_ms1: Option<DiaMs1WindowParams>,
    scale: BrukerMobilityScale,
}

impl BrukerDenoiseConfig {
    pub fn params(&self) -> DenoiseParams {
        DenoiseParams {
            filter: FilterParams {
                mz_half_width: self.mz_half_width,
                min_feature_length: self.min_feature_length,
                max_internal_gap: self.max_internal_gap,
                min_window_intensity: self.min_window_intensity,
                min_feature_intensity: self.min_feature_intensity,
                num_iterations: self.iterations,
            },
            halo: self.halo.then_some(HaloParams {
                peak_fraction: self.halo_peak_fraction,
                mz_idx_half_width: self.halo_mz_idx_half_width,
                scan_half_width: self.halo_scan_half_width,
            }),
            polygon: self.ms1_polygon.then_some(Ms1PolygonParams {
                mz_pad: self.ms1_polygon_mz_pad,
                im_pad: self.ms1_polygon_im_pad,
                overlap: self.ms1_polygon_overlap,
            }),
            dia_ms1: self.dia_ms1_window.then_some(DiaMs1WindowParams {
                mz_pad: self.dia_ms1_mz_pad,
                im_pad: self.dia_ms1_im_pad,
                overlap: self.dia_ms1_overlap,
            }),
            scale: self.mobility_scale,
        }
    }

    /// Reject settings dnoise would reject, before any file is read.
    pub fn validate(&self) -> Result<(), DenoiseError> {
        let params = self.params();
        dnoise_core::pipeline::validate(&params.filter, &params.stages())
            .map_err(|error| DenoiseError::Config(error.to_string()))
    }
}

impl DenoiseParams {
    fn stages(&self) -> Stages<'_> {
        Stages {
            halo: self.halo.as_ref(),
            ms1_polygon: self.polygon.as_ref(),
            dia_ms1: self.dia_ms1.as_ref(),
            mobility_scale: match self.scale {
                BrukerMobilityScale::Calibrated => dnoise_core::MobilityScale::Calibrated,
                BrukerMobilityScale::Linear => dnoise_core::MobilityScale::Linear,
            },
            ..Stages::default()
        }
    }

    /// Build the run's denoiser from its `analysis.tdf`: frame table,
    /// acquisition scheme and the MS1 gate that scheme uses.
    pub fn denoiser(&self, tdf: &Path) -> Result<Ms1Denoiser<'_>, DenoiseError> {
        let db = Tdf::open(tdf)?;
        let meta = db.frame_meta()?;
        let acquisition = acquisition(&meta);
        let num_scans = meta.iter().map(|m| m.num_scans).max().unwrap_or(0);
        let mut denoiser = Denoiser::new(&self.filter, &self.stages(), acquisition, meta, false)
            .map_err(|source| db.dnoise(source))?;

        let stages = *denoiser.stages();
        if num_scans > 0 {
            if let Some(p) = stages.ms1_polygon {
                denoiser.set_polygon_gate(db.polygon_gate(p, self.scale, num_scans)?);
            }
            if let Some(p) = stages.dia_ms1 {
                denoiser.set_dia_ms1_gate(db.dia_ms1_gate(p, self.scale, num_scans)?);
            }
        }
        let gate = if denoiser.polygon_gate().is_some() {
            "selection polygon"
        } else if denoiser.dia_ms1_gate().is_some() {
            "diaPASEF isolation windows"
        } else {
            "none"
        };
        log::info!(
            "{}: denoising MS1 frames ({acquisition}; MS1 gate: {gate})",
            tdf.display()
        );
        Ok(Ms1Denoiser {
            denoiser,
            path: db.path,
        })
    }
}

/// A run's [`Denoiser`], applied to one timsrust MS1 frame at a time.
pub struct Ms1Denoiser<'a> {
    denoiser: Denoiser<'a>,
    path: PathBuf,
}

impl Ms1Denoiser<'_> {
    /// Surviving `(scan, tof, intensity)` points of MS1 frame `frame_id`, in
    /// scan order. `scan_offsets` is the frame's CSR row pointer; it may be
    /// shorter than `Frames.NumScans` when trailing scans are empty.
    pub fn denoise(
        &self,
        frame_id: usize,
        scan_offsets: &[usize],
        tof: Vec<u32>,
        intensity: Vec<u32>,
    ) -> Result<Vec<(u32, u32, u32)>, DenoiseError> {
        let meta = self.denoiser.meta();
        // Frames.Id is normally 1-based and contiguous.
        let index = match frame_id
            .checked_sub(1)
            .filter(|&i| meta.get(i).is_some_and(|m| m.id == frame_id))
        {
            Some(index) => index,
            None => meta
                .binary_search_by_key(&frame_id, |m| m.id)
                .map_err(|_| {
                    self.metadata(format!("frame {frame_id} is not in the Frames table"))
                })?,
        };
        let num_scans = meta[index].num_scans;
        let scans = scan_offsets.len().saturating_sub(1);
        if scans > num_scans {
            return Err(self.metadata(format!(
                "frame {frame_id} has {scans} scans but Frames.NumScans is {num_scans}"
            )));
        }
        let mut scan = Vec::with_capacity(tof.len());
        for (s, w) in scan_offsets.windows(2).enumerate() {
            scan.extend(std::iter::repeat_n(s as u32, w[1] - w[0]));
        }
        let frame = FlatFrame {
            frame_id,
            num_scans,
            scan,
            tof,
            intensity,
        };
        self.denoiser
            .denoise_ms1(index, &frame)
            .map(|decoded| decoded.survivors)
            .map_err(|source| DenoiseError::Dnoise {
                path: self.path.clone(),
                source,
            })
    }

    fn metadata(&self, message: String) -> DenoiseError {
        DenoiseError::Metadata {
            path: self.path.clone(),
            message,
        }
    }
}

/// Acquisition scheme from the non-MS1 `MsMsType` values, as dnoise detects it.
pub(crate) fn acquisition(meta: &[FrameMeta]) -> Acquisition {
    let kinds: BTreeSet<i64> = meta
        .iter()
        .filter(|m| !m.is_ms1())
        .map(|m| m.ms_ms_type)
        .collect();
    if kinds.len() > 1 {
        return Acquisition::Mixed;
    }
    match kinds.first() {
        None => Acquisition::Ms1Only,
        Some(8) => Acquisition::DdaPasef,
        Some(9) => Acquisition::DiaPasef,
        Some(10) => Acquisition::PrmPasef,
        Some(_) => Acquisition::Unknown,
    }
}

/// Polygon vertex columns: m/z values and the matching 1/K0 values.
type PolygonVertices = (Vec<f64>, Vec<f64>);

struct Tdf {
    path: PathBuf,
    connection: Connection,
}

impl Tdf {
    fn open(path: &Path) -> Result<Self, DenoiseError> {
        let connection = Connection::open_with_flags(path, OpenFlags::SQLITE_OPEN_READ_ONLY)
            .map_err(|source| DenoiseError::Sql {
                path: path.to_path_buf(),
                source,
            })?;
        Ok(Self {
            path: path.to_path_buf(),
            connection,
        })
    }

    fn sql(&self, source: rusqlite::Error) -> DenoiseError {
        DenoiseError::Sql {
            path: self.path.clone(),
            source,
        }
    }

    fn metadata(&self, message: impl Into<String>) -> DenoiseError {
        DenoiseError::Metadata {
            path: self.path.clone(),
            message: message.into(),
        }
    }

    fn dnoise(&self, source: dnoise_core::Error) -> DenoiseError {
        DenoiseError::Dnoise {
            path: self.path.clone(),
            source,
        }
    }

    fn has_table(&self, name: &str) -> Result<bool, DenoiseError> {
        self.connection
            .query_row(
                "SELECT EXISTS(SELECT 1 FROM sqlite_master WHERE type='table' AND name=?1)",
                [name],
                |row| row.get(0),
            )
            .map_err(|e| self.sql(e))
    }

    fn frame_meta(&self) -> Result<Vec<FrameMeta>, DenoiseError> {
        let mut statement = self
            .connection
            .prepare("SELECT Id, NumScans, NumPeaks, MsMsType, Time FROM Frames ORDER BY Id")
            .map_err(|e| self.sql(e))?;
        let rows = statement
            .query_map([], |r| {
                Ok(FrameMeta {
                    id: r.get::<_, i64>(0)? as usize,
                    num_scans: r.get::<_, i64>(1)? as usize,
                    num_peaks: r.get::<_, i64>(2)? as u64,
                    ms_ms_type: r.get(3)?,
                    rt: r.get(4)?,
                })
            })
            .map_err(|e| self.sql(e))?;
        rows.collect::<Result<_, _>>().map_err(|e| self.sql(e))
    }

    /// dnoise converts gate geometry once with the run-level calibration, which
    /// is exact only when every frame shares one calibration row.
    fn require_single_calibration(&self, gate: &str, key: &str) -> Result<(), DenoiseError> {
        let mut columns = BTreeSet::new();
        let mut statement = self
            .connection
            .prepare("PRAGMA table_info(Frames)")
            .map_err(|e| self.sql(e))?;
        let rows = statement
            .query_map([], |r| r.get::<_, String>(1))
            .map_err(|e| self.sql(e))?;
        for column in rows {
            columns.insert(column.map_err(|e| self.sql(e))?);
        }
        let count = |column: &str| -> Result<i64, DenoiseError> {
            if !columns.contains(column) {
                return Ok(1);
            }
            self.connection
                .query_row(
                    &format!("SELECT COUNT(DISTINCT {column}) FROM Frames"),
                    [],
                    |r| r.get(0),
                )
                .map_err(|e| self.sql(e))
        };
        let (mz, im) = (count("MzCalibration")?, count("TimsCalibration")?);
        if mz > 1 || im > 1 {
            return Err(self.metadata(format!(
                "the {gate} gate needs one calibration for the whole run, but frames reference \
                 {mz} m/z and {im} mobility calibrations; set `bruker_config.denoise.{key}` to false"
            )));
        }
        Ok(())
    }

    fn global_metadata(&self) -> Result<HashMap<String, String>, DenoiseError> {
        let mut statement = self
            .connection
            .prepare("SELECT Key, Value FROM GlobalMetadata")
            .map_err(|e| self.sql(e))?;
        let rows = statement
            .query_map([], |r| Ok((r.get::<_, String>(0)?, r.get::<_, String>(1)?)))
            .map_err(|e| self.sql(e))?;
        rows.collect::<Result<_, _>>().map_err(|e| self.sql(e))
    }

    /// The run-level converters dnoise uses: timsrust 0.4.2's m/z line (with the
    /// otofControl range widening) and the chosen scan-to-1/K0 scale.
    fn converters(
        &self,
        scale: BrukerMobilityScale,
        num_scans: usize,
    ) -> Result<(Tof2MzConverter, ScanToMobility), DenoiseError> {
        let metadata = self.global_metadata()?;
        let number = |key: &str| -> Result<f64, DenoiseError> {
            metadata
                .get(key)
                .and_then(|value| value.trim().parse().ok())
                .ok_or_else(|| {
                    self.metadata(format!("GlobalMetadata {key} is missing or not a number"))
                })
        };
        let (mut mz_min, mut mz_max) = (number("MzAcqRangeLower")?, number("MzAcqRangeUpper")?);
        if metadata.get("AcquisitionSoftware").map(String::as_str) == Some("Bruker otofControl") {
            mz_min -= 5.0;
            mz_max += 5.0;
        }
        let tof_max = number("DigitizerNumSamples")? as u32;
        let mz = Tof2MzConverter::from_boundaries(mz_min, mz_max, tof_max);
        let im = match scale {
            BrukerMobilityScale::Linear => {
                ScanToMobility::Linear(Scan2ImConverter::from_boundaries(
                    number("OneOverK0AcqRangeLower")?,
                    number("OneOverK0AcqRangeUpper")?,
                    num_scans as u32,
                ))
            }
            BrukerMobilityScale::Calibrated => {
                let ids: Vec<i64> = self
                    .connection
                    .prepare("SELECT DISTINCT TimsCalibration FROM Frames")
                    .and_then(|mut s| s.query_map([], |r| r.get(0))?.collect())
                    .map_err(|e| self.sql(e))?;
                let [id] = ids[..] else {
                    return Err(self.metadata(format!(
                        "frames reference {} TimsCalibration rows, expected one",
                        ids.len()
                    )));
                };
                let (model, c) = self
                    .connection
                    .query_row(
                        "SELECT ModelType, C0, C1, C2, C3, C4, C5, C6, C7, C8, C9 \
                         FROM TimsCalibration WHERE Id = ?1",
                        [id],
                        |r| {
                            let mut c = [0.0; 10];
                            for (i, v) in c.iter_mut().enumerate() {
                                *v = r.get(i + 1)?;
                            }
                            Ok((r.get::<_, i64>(0)?, c))
                        },
                    )
                    .map_err(|e| self.sql(e))?;
                ScanToMobility::Calibrated(TimsCalibrationModel::new(model, c).map_err(|e| {
                    self.metadata(format!(
                        "{e}; set `bruker_config.denoise.mobility_scale` to \"linear\""
                    ))
                })?)
            }
        };
        Ok((mz, im))
    }

    /// The ddaPASEF selection polygon as `(m/z, 1/K0)` vertex columns, or `None`.
    fn selection_polygon(&self) -> Result<Option<PolygonVertices>, DenoiseError> {
        if !self.has_table("PropertyDefinitions")? || !self.has_table("GroupProperties")? {
            return Ok(None);
        }
        let property = |name: &str| -> Result<Option<i64>, DenoiseError> {
            self.connection
                .query_row(
                    "SELECT Id FROM PropertyDefinitions WHERE PermanentName=?1",
                    [name],
                    |r| r.get(0),
                )
                .optional()
                .map_err(|e| self.sql(e))
        };
        let blob = |id: i64| -> Result<Option<Vec<f64>>, DenoiseError> {
            let blob: Option<Vec<u8>> = self
                .connection
                .query_row(
                    "SELECT Value FROM GroupProperties WHERE Property=?1 LIMIT 1",
                    [id],
                    |r| r.get(0),
                )
                .optional()
                .map_err(|e| self.sql(e))?;
            Ok(blob.map(|b| {
                b.chunks_exact(8)
                    .map(|c| f64::from_le_bytes(c.try_into().expect("8-byte chunk")))
                    .collect()
            }))
        };
        let (Some(mz_id), Some(im_id)) = (
            property("IMS_PolygonFilter_Mass")?,
            property("IMS_PolygonFilter_Mobility")?,
        ) else {
            return Ok(None);
        };
        let (Some(mz), Some(im)) = (blob(mz_id)?, blob(im_id)?) else {
            return Ok(None);
        };
        Ok((mz.len() >= 3 && mz.len() == im.len()).then_some((mz, im)))
    }

    /// True when some frame maps to a nonempty diaPASEF window group, as
    /// dnoise's `read_dia_windows` would report.
    fn has_dia_windows(&self) -> Result<bool, DenoiseError> {
        if !self.has_table("DiaFrameMsMsInfo")? || !self.has_table("DiaFrameMsMsWindows")? {
            return Ok(false);
        }
        self.connection
            .query_row(
                "SELECT EXISTS(SELECT 1 FROM DiaFrameMsMsInfo i JOIN DiaFrameMsMsWindows w \
                 ON w.WindowGroup = i.WindowGroup WHERE w.ScanNumEnd > w.ScanNumBegin)",
                [],
                |r| r.get(0),
            )
            .map_err(|e| self.sql(e))
    }

    fn dia_ms1_boxes(&self) -> Result<Vec<DiaMs1Box>, DenoiseError> {
        if !self.has_table("DiaFrameMsMsWindows")? {
            return Ok(Vec::new());
        }
        let mut statement = self
            .connection
            .prepare(
                "SELECT DISTINCT ScanNumBegin, ScanNumEnd, IsolationMz, IsolationWidth \
                 FROM DiaFrameMsMsWindows \
                 WHERE IsolationMz IS NOT NULL AND IsolationWidth IS NOT NULL",
            )
            .map_err(|e| self.sql(e))?;
        let rows = statement
            .query_map([], |r| {
                let mz: f64 = r.get(2)?;
                let width: f64 = r.get(3)?;
                Ok(DiaMs1Box {
                    scan_begin: r.get::<_, i64>(0)? as u32,
                    scan_end: r.get::<_, i64>(1)? as u32,
                    mz_lo: mz - width / 2.0,
                    mz_hi: mz + width / 2.0,
                })
            })
            .map_err(|e| self.sql(e))?;
        let boxes: Vec<DiaMs1Box> = rows.collect::<Result<_, _>>().map_err(|e| self.sql(e))?;
        Ok(boxes
            .into_iter()
            .filter(|b| b.scan_end > b.scan_begin)
            .collect())
    }

    /// dnoise's `build_polygon_gate`. Skipped on runs with diaPASEF windows,
    /// whose polygon property holds window anchors rather than one ring.
    fn polygon_gate(
        &self,
        p: &Ms1PolygonParams,
        scale: BrukerMobilityScale,
        num_scans: usize,
    ) -> Result<Option<PolygonGate>, DenoiseError> {
        if self.has_dia_windows()? {
            return Ok(None);
        }
        let Some((mz, im)) = self.selection_polygon()? else {
            return Ok(None);
        };
        self.require_single_calibration("MS1 polygon", "ms1_polygon")?;
        let (mz_converter, k0) = self.converters(scale, num_scans)?;
        let im_at_scan = |s: u32| k0.convert(s as f64);
        let mz_to_tof = |value: f64| mz_converter.invert(value);
        let Some(mut gate) = PolygonGate::build(
            &mz, &im, num_scans, im_at_scan, mz_to_tof, p.mz_pad, p.im_pad,
        ) else {
            return Ok(None);
        };
        gate.check_contains_unpadded(&mz, &im, im_at_scan, mz_to_tof)
            .map_err(|e| self.metadata(format!("MS1 polygon gate self-check failed: {e}")))?;
        gate.overlap = p.overlap;
        Ok(Some(gate))
    }

    /// dnoise's `build_dia_ms1_gate`.
    fn dia_ms1_gate(
        &self,
        p: &DiaMs1WindowParams,
        scale: BrukerMobilityScale,
        num_scans: usize,
    ) -> Result<Option<dnoise_core::dia_ms1::DiaMs1Gate>, DenoiseError> {
        let boxes = self.dia_ms1_boxes()?;
        if boxes.is_empty() {
            return Ok(None);
        }
        self.require_single_calibration("diaPASEF MS1 window", "dia_ms1_window")?;
        let (mz_converter, k0) = self.converters(scale, num_scans)?;
        Ok(dnoise_core::dia_ms1::DiaMs1Gate::from_windows(
            &boxes,
            p,
            |value: f64| mz_converter.invert(value),
            &k0,
            num_scans,
        ))
    }
}

#[cfg(test)]
#[path = "../tests/unit/denoise.rs"]
mod tests;
