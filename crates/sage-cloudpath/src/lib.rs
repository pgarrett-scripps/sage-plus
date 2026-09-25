use async_compression::tokio::bufread::GzipDecoder;
use async_compression::tokio::write::GzipEncoder;
#[cfg(feature = "cloud")]
use futures::TryStreamExt;
#[cfg(feature = "cloud")]
use object_store::{ObjectStore, ObjectStoreExt};
use std::io::Write;
#[cfg(feature = "cloud")]
use std::sync::Arc;
use tokio::io::{AsyncBufRead, AsyncRead, AsyncWriteExt, BufReader};

pub use url::Url;

pub mod denoise;
pub mod mgf;
pub mod mzml;
#[cfg(feature = "mzmlb")]
pub mod mzmlb;
pub mod tdf;
pub mod thermoraw;
pub mod tims_mobility;
pub mod util;
pub use util::FileFormat;

#[cfg(feature = "parquet")]
pub mod parquet;

/// Schemes recognized by `object_store::parse_url_opts`. Anything outside
/// this set — most importantly Windows drive letters like `C:` which parse
/// as single-letter URL schemes — is treated as a local path.
const OBJECT_STORE_SCHEMES: &[&str] = &[
    "file", "memory", "s3", "s3a", "gs", "az", "adl", "azure", "abfs", "abfss", "http", "https",
];

/// Parse `s` as a URL, but only accept schemes that `object_store` knows how
/// to handle. Returns `None` for local paths (including Windows paths like
/// `C:\foo` that would otherwise parse as a URL with scheme `c`).
pub fn try_parse_url(s: &str) -> Option<Url> {
    Url::parse(s)
        .ok()
        .filter(|u| OBJECT_STORE_SCHEMES.contains(&u.scheme()))
}

/// Convert a path string (local path or cloud URL) into a [`Url`].
pub fn to_url(s: &str) -> Result<Url, Error> {
    if let Some(url) = try_parse_url(s) {
        return Ok(url);
    }
    let path = std::path::Path::new(s);
    let canonical = path.canonicalize()?;
    Url::from_file_path(&canonical).map_err(|_| Error::InvalidUri)
}

/// Does the URL path end in "gz" or "gzip"?
fn gzip_heuristic(url: &Url) -> bool {
    let p = url.path();
    p.ends_with("gz") || p.ends_with("gzip")
}

/// Return the filename portion of a URL path. If the filename ends with `.tdf`,
/// return the parent directory name instead (Bruker `.d` convention).
pub fn filename(url: &Url) -> Option<&str> {
    let path = url.path();
    let name = path.rsplit('/').next().filter(|s| !s.is_empty());
    match name {
        Some(n) if n.ends_with("tdf") => {
            let mut iter = path.rsplit('/');
            iter.next();
            iter.next().filter(|s| !s.is_empty())
        }
        other => other,
    }
}

#[cfg(feature = "cloud")]
fn parse_url(url: &Url) -> Result<(Box<dyn ObjectStore>, object_store::path::Path), Error> {
    // AWS and Azure require lowercased config keys. By default, these aren't pulled from the env
    object_store::parse_url_opts(
        url,
        std::env::vars().map(|(k, v)| (k.to_ascii_lowercase(), v)),
    )
    .map_err(Error::ObjectStore)
}

/// A bounded-memory synchronous writer backed by object-store multipart uploads.
///
/// CSV generation in the CLI is synchronous. This adapter buffers writes and
/// drives the asynchronous object-store writer on an internal current-thread
/// runtime, switching to multipart upload once its 10 MiB buffer is full.
#[cfg(feature = "cloud")]
pub struct CloudWriter {
    runtime: tokio::runtime::Runtime,
    writer: object_store::buffered::BufWriter,
}

/// Without the `cloud` feature, a remote writer cannot be constructed.
#[cfg(not(feature = "cloud"))]
pub struct CloudWriter {
    _unconstructible: std::convert::Infallible,
}

#[cfg(not(feature = "cloud"))]
impl CloudWriter {
    pub fn new(url: &Url) -> Result<Self, Error> {
        Err(Error::CloudDisabled(url.scheme().to_string()))
    }

    pub fn finish(self) -> Result<(), Error> {
        match self._unconstructible {}
    }
}

#[cfg(not(feature = "cloud"))]
impl Write for CloudWriter {
    fn write(&mut self, _buf: &[u8]) -> std::io::Result<usize> {
        match self._unconstructible {}
    }

    fn flush(&mut self) -> std::io::Result<()> {
        match self._unconstructible {}
    }
}

#[cfg(feature = "cloud")]
impl CloudWriter {
    pub fn new(url: &Url) -> Result<Self, Error> {
        let (store, path) = parse_url(url)?;
        let runtime = tokio::runtime::Builder::new_current_thread()
            .enable_all()
            .build()?;
        let store: Arc<dyn ObjectStore> = Arc::from(store);
        Ok(Self {
            runtime,
            writer: object_store::buffered::BufWriter::new(store, path),
        })
    }

    pub fn finish(mut self) -> Result<(), Error> {
        self.runtime.block_on(self.writer.shutdown())?;
        Ok(())
    }
}

#[cfg(feature = "cloud")]
impl Write for CloudWriter {
    fn write(&mut self, buf: &[u8]) -> std::io::Result<usize> {
        self.runtime.block_on(self.writer.write(buf))
    }

    fn flush(&mut self) -> std::io::Result<()> {
        self.runtime.block_on(self.writer.flush())
    }
}

/// Open a streaming reader for the given URL.
async fn read_url(url: &Url) -> Result<Box<dyn AsyncBufRead + Unpin + Send>, Error> {
    let reader: BufReader<Box<dyn AsyncRead + Unpin + Send>> = BufReader::new(open_url(url).await?);

    if gzip_heuristic(url) {
        Ok(Box::new(BufReader::new(GzipDecoder::new(reader))))
    } else {
        Ok(Box::new(reader))
    }
}

#[cfg(feature = "cloud")]
async fn open_url(url: &Url) -> Result<Box<dyn AsyncRead + Unpin + Send>, Error> {
    let (store, obj_path) = parse_url(url)?;
    let result = store.get(&obj_path).await.map_err(Error::ObjectStore)?;
    let stream = result.into_stream().map_err(std::io::Error::other);
    Ok(Box::new(tokio_util::io::StreamReader::new(stream)))
}

#[cfg(not(feature = "cloud"))]
async fn open_url(url: &Url) -> Result<Box<dyn AsyncRead + Unpin + Send>, Error> {
    let path = local_path(url)?;
    Ok(Box::new(tokio::fs::File::open(path).await?))
}

/// Resolve a `file://` URL, or explain that other schemes need the `cloud` feature.
#[cfg(not(feature = "cloud"))]
fn local_path(url: &Url) -> Result<std::path::PathBuf, Error> {
    if url.scheme() != "file" {
        return Err(Error::CloudDisabled(url.scheme().to_string()));
    }
    url.to_file_path().map_err(|_| Error::InvalidUri)
}

/// Write bytes to the given URL. Gzip-compresses if the path ends in `.gz`.
pub async fn write_bytes_async(url: &Url, bytes: Vec<u8>) -> Result<(), Error> {
    let bytes: Vec<u8> = if gzip_heuristic(url) {
        let inner = Vec::with_capacity(bytes.len() / 2);
        let mut wtr = GzipEncoder::new(inner);
        wtr.write_all(&bytes).await?;
        wtr.shutdown().await?;
        wtr.into_inner()
    } else {
        bytes
    };

    // Ensure parent directories exist for local paths
    if let Ok(local) = url.to_file_path() {
        if let Some(parent) = local.parent() {
            std::fs::create_dir_all(parent)?;
        }
    }

    put_bytes(url, bytes).await
}

#[cfg(feature = "cloud")]
async fn put_bytes(url: &Url, bytes: Vec<u8>) -> Result<(), Error> {
    let (store, obj_path) = parse_url(url)?;
    store
        .put(&obj_path, bytes::Bytes::from(bytes).into())
        .await
        .map_err(Error::ObjectStore)?;
    Ok(())
}

#[cfg(not(feature = "cloud"))]
async fn put_bytes(url: &Url, bytes: Vec<u8>) -> Result<(), Error> {
    let path = local_path(url)?;
    // Write beside the target and rename, as object_store's local store does, so
    // readers never observe a partially written file.
    let file_name = path.file_name().ok_or(Error::InvalidUri)?;
    let mut staging = file_name.to_os_string();
    staging.push(format!(".{}.tmp", std::process::id()));
    let staging = path.with_file_name(staging);
    tokio::fs::write(&staging, bytes).await?;
    if let Err(error) = tokio::fs::rename(&staging, &path).await {
        let _ = tokio::fs::remove_file(&staging).await;
        return Err(error.into());
    }
    Ok(())
}

/// Synchronous wrapper around [`write_bytes_async`].
pub fn write_bytes_sync(url: &Url, bytes: Vec<u8>) -> Result<(), Error> {
    let rt = tokio::runtime::Builder::new_current_thread()
        .enable_all()
        .build()?;
    rt.block_on(write_bytes_async(url, bytes))
}

/// Open a reader for `path`, then execute `func` on it. Handles both local and
/// remote paths, with transparent gzip decompression.
pub fn read_and_execute<P, F, Fut, T>(path: P, func: F) -> Result<T, Error>
where
    P: AsRef<str>,
    Fut: futures::Future<Output = Result<T, Error>>,
    F: FnOnce(Box<dyn AsyncBufRead + Unpin>) -> Fut,
{
    let url = to_url(path.as_ref())?;

    // Reject remote URLs that have no object key
    if url.scheme() != "file" {
        let key = url.path().strip_prefix('/').unwrap_or(url.path());
        if key.is_empty() {
            return Err(Error::InvalidUri);
        }
    }

    let rt = tokio::runtime::Builder::new_current_thread()
        .enable_all()
        .build()?;

    rt.block_on(async {
        let reader = read_url(&url).await?;
        func(reader).await
    })
}

#[derive(thiserror::Error, Debug)]
pub enum Error {
    #[error("invalid uri")]
    InvalidUri,
    #[error("unsupported input option: {0}")]
    Unsupported(String),
    #[cfg(feature = "cloud")]
    #[error("object store error: {0}")]
    ObjectStore(#[from] object_store::Error),
    #[error(
        "`{0}://` paths need cloud storage support, which this build of Sage Plus omits. \
         Use a local path, or build with the `cloud` feature (enabled by default)"
    )]
    CloudDisabled(String),
    #[error(transparent)]
    IO(#[from] tokio::io::Error),
    #[error(transparent)]
    Json(#[from] serde_json::Error),
    #[error("MzML error: {0}")]
    MzML(#[from] mzml::MzMLError),
    #[error("TDF error: {0}")]
    TDF(#[from] timsrust::TimsRustError),
    #[error(transparent)]
    MobilityCalibration(#[from] tims_mobility::MobilityCalibrationError),
    #[error(transparent)]
    Denoise(#[from] denoise::DenoiseError),
    #[error("Thermo RAW error: {0}")]
    ThermoRaw(#[from] opentfraw::Error),
    #[error("MGF error: {0}")]
    MGF(#[from] mgf::MgfError),
    #[error("FASTA error: {0}")]
    Fasta(#[from] sage_core::fasta::FastaError),
}

#[cfg(test)]
#[path = "../tests/unit/lib.rs"]
mod test;
