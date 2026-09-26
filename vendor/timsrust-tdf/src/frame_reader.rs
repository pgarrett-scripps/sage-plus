// Legacy compression types remain reachable for upstream error API compatibility.
#![allow(dead_code)]
pub(crate) mod compression1;
pub(crate) mod compression2;
pub(crate) mod frame_info_reader;

use std::collections::HashMap;

use timsrust_core::{FrameIons, utils::reader::Reader};

pub(crate) use frame_info_reader::FrameInfoReader;

use crate::{
    Metadata, TDFPathError,
    frame_reader::{
        compression1::{TdfBlobReaderCompression1, TdfBlobReaderErrorCompression1},
        compression2::{TdfBlobReader, TdfBlobReaderError},
    },
};

pub use frame_info_reader::FrameReaderErrorInternal;

use super::{
    MetadataReaderError, QuadrupoleSettingsReaderError, TDFPathLike,
    file_readers::sql_reader::SqlReaderError,
};

/// Adapts a raw blob reader (indexed by binary file offset) to be indexed by
/// frame id, using the offset table from [`FrameInfoReader`].
#[derive(Debug)]
struct TdfOffsetIonReader<B> {
    blob_reader: B,
    offsets: HashMap<usize, usize>,
}

impl<B> TdfOffsetIonReader<B> {
    fn new(blob_reader: B, offsets: HashMap<usize, usize>) -> Self {
        Self {
            blob_reader,
            offsets,
        }
    }
}

impl<B: Reader<FrameIons>> Reader<FrameIons> for TdfOffsetIonReader<B>
where
    FrameReaderError: From<B::Error>,
{
    type Error = FrameReaderError;

    fn get(&self, index: usize) -> Result<FrameIons, Self::Error> {
        let offset = self
            .offsets
            .get(&index)
            .copied()
            .ok_or(FrameReaderError::IndexOutOfBounds)?;
        self.blob_reader.get(offset).map_err(FrameReaderError::from)
    }
}

/// Unifies the two TDF compression schemes into a single [`Reader<FrameIons>`].
#[allow(private_interfaces)]
#[derive(Debug)]
pub enum TdfIonReader {
    Compression1(TdfOffsetIonReader<TdfBlobReaderCompression1>),
    Native(NativeIonReader),
    Compression2(TdfOffsetIonReader<TdfBlobReader>),
}

impl Reader<FrameIons> for TdfIonReader {
    type Error = FrameReaderError;

    fn get(&self, index: usize) -> Result<FrameIons, FrameReaderError> {
        match self {
            Self::Compression1(r) => r.get(index),
            Self::Native(r) => r.get(index),
            Self::Compression2(r) => r.get(index),
        }
    }
}

/// A concrete frame reader for Bruker TDF files.
///
/// Thin newtype around
/// [`timsrust_core::FrameReader<TdfIonReader, FrameInfoReader>`].  All
/// [`FrameReader`] methods (`get_frame`, `get_info`, `iter_indices`,
/// `parallel_filter`, …) are available via [`Deref`].
#[derive(Debug)]
pub struct TdfFrameReader(timsrust_core::FrameReader<TdfIonReader, FrameInfoReader>);

impl std::ops::Deref for TdfFrameReader {
    type Target = timsrust_core::FrameReader<TdfIonReader, FrameInfoReader>;

    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl TdfFrameReader {
    pub fn new(path: impl TDFPathLike) -> Result<Self, FrameReaderError> {
        let metadata = Metadata::new(&path)?;
        Self::without_metadata(
            &path,
            metadata.compression_type(),
            metadata.max_peaks_per_scan(),
        )
    }

    pub fn without_metadata(
        path: impl TDFPathLike,
        compression_type: u8,
        max_peaks_per_scan: usize,
    ) -> Result<Self, FrameReaderError> {
        let _ = (compression_type, max_peaks_per_scan);
        let source = path.to_timstof_path()?;
        let native = timsrust_reader::TdfReader::open(source.tdf().as_ref())?;
        let info_reader = FrameInfoReader::new(&path)?;
        let ion_reader = TdfIonReader::Native(NativeIonReader(std::sync::Mutex::new(native)));
        Ok(Self(timsrust_core::FrameReader::new(
            ion_reader,
            info_reader,
        )))
    }

    /// Consume `self` and return the underlying generic
    /// [`timsrust_core::FrameReader`].
    pub fn into_inner(self) -> timsrust_core::FrameReader<TdfIonReader, FrameInfoReader> {
        self.0
    }

    /// Return the acquisition type detected from the frame metadata.
    pub fn get_acquisition(&self) -> timsrust_core::AcquisitionType {
        self.info_reader().get_acquisition()
    }
}

#[allow(private_interfaces)]
#[derive(Debug, thiserror::Error)]
pub enum FrameReaderError {
    #[error("{0}")]
    Native(#[from] timsrust_reader::Error),
    #[error("Timscompress error")]
    TimscompressError,
    #[error("{0}")]
    TdfBlobReaderError(#[from] TdfBlobReaderError),
    #[error("{0}")]
    TdfBlobReaderErrorCompression1(#[from] TdfBlobReaderErrorCompression1),
    #[error("{0}")]
    MetadataReaderError(#[from] MetadataReaderError),
    #[error("{0}")]
    FileNotFound(String),
    #[error("{0}")]
    SqlReaderError(#[from] SqlReaderError),
    #[error("Corrupt Frame")]
    CorruptFrame,
    #[error("{0}")]
    QuadrupoleSettingsReaderError(#[from] QuadrupoleSettingsReaderError),
    #[error("Index out of bounds")]
    IndexOutOfBounds,
    #[error("Compression type {0} not understood")]
    CompressionTypeError(u8),
    #[error("Failed to read path: {0}")]
    PathError(#[from] TDFPathError),
    #[error("Got unexpected TdfBlob type")]
    UnexpectedTdfBlobError,
    #[error("{0}")]
    FrameInfoReaderError(#[from] frame_info_reader::FrameReaderErrorInternal),
    #[error("{0}")]
    CoreFrameReaderError(#[from] timsrust_core::FrameReaderError),
}

/// Consumer adapter preserving the upstream spectrum processing contract.
pub struct NativeIonReader(std::sync::Mutex<timsrust_reader::TdfReader>);
impl std::fmt::Debug for NativeIonReader {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("NativeIonReader").finish_non_exhaustive()
    }
}
impl Reader<FrameIons> for NativeIonReader {
    type Error = FrameReaderError;
    fn get(&self, index: usize) -> Result<FrameIons, Self::Error> {
        let frame = self
            .0
            .lock()
            .map_err(|_| FrameReaderError::CorruptFrame)?
            .read_frame(index as u64)?;
        let tof = frame
            .tof_indices
            .into_iter()
            .map(timsrust_core::TofIndex::try_from)
            .collect::<Result<Vec<_>, _>>()
            .map_err(|_| FrameReaderError::CorruptFrame)?;
        let intensity = frame
            .intensities
            .into_iter()
            .map(timsrust_core::IntensityIndex::try_from)
            .collect::<Result<Vec<_>, _>>()
            .map_err(|_| FrameReaderError::CorruptFrame)?;
        Ok(FrameIons::new(frame.scan_offsets, tof, intensity))
    }
}
