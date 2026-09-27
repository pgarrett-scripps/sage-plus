//! Streaming content hashes of input files, for run provenance.

use crate::{read_and_execute, Error};
use sha2::{Digest, Sha256};
use std::io::Read;
use tokio::io::{AsyncBufReadExt, AsyncReadExt};
use url::Url;

/// Bytes read per step while hashing. Hashing never holds more than this.
const CHUNK_BYTES: usize = 1 << 20;

/// Size and SHA-256 of the bytes stored in one file.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct FileHash {
    pub size_bytes: u64,
    /// Lowercase hex, as printed by `sha256sum`.
    pub sha256: String,
}

/// Hash the bytes stored at `url`, as stored: gzip files are not decompressed,
/// so the result matches `sha256sum` on the file. Local files and remote
/// objects are both streamed in 1 MiB chunks and never loaded whole.
/// Directories (such as Bruker `.d` inputs) are an error.
pub fn sha256_url(url: &Url) -> Result<FileHash, Error> {
    if url.scheme() == "file" {
        let path = url.to_file_path().map_err(|_| Error::InvalidUri)?;
        return sha256_reader(std::fs::File::open(path)?);
    }
    let rt = tokio::runtime::Builder::new_current_thread()
        .enable_all()
        .build()?;
    rt.block_on(async {
        let mut reader = crate::open_url(url).await?;
        let mut hasher = Sha256::new();
        let mut buffer = vec![0; CHUNK_BYTES];
        let mut size_bytes = 0;
        loop {
            let read = reader.read(&mut buffer).await?;
            if read == 0 {
                break;
            }
            hasher.update(&buffer[..read]);
            size_bytes += read as u64;
        }
        Ok(FileHash {
            size_bytes,
            sha256: hex(&hasher.finalize()),
        })
    })
}

/// Hash everything `reader` yields, in 1 MiB chunks.
pub fn sha256_reader(mut reader: impl Read) -> Result<FileHash, Error> {
    let mut hasher = Sha256::new();
    let mut buffer = vec![0; CHUNK_BYTES];
    let mut size_bytes = 0;
    loop {
        let read = match reader.read(&mut buffer) {
            Ok(0) => break,
            Ok(read) => read,
            Err(error) if error.kind() == std::io::ErrorKind::Interrupted => continue,
            Err(error) => return Err(error.into()),
        };
        hasher.update(&buffer[..read]);
        size_bytes += read as u64;
    }
    Ok(FileHash {
        size_bytes,
        sha256: hex(&hasher.finalize()),
    })
}

/// Call `visit` on each line of the text at `path`, without its line ending.
/// Gzip files are decompressed; the text is streamed, not loaded whole.
pub fn for_each_line<P, F>(path: P, mut visit: F) -> Result<(), Error>
where
    P: AsRef<str>,
    F: FnMut(&str),
{
    read_and_execute(path, |mut reader| async move {
        let mut line = String::new();
        loop {
            line.clear();
            if reader.read_line(&mut line).await? == 0 {
                return Ok(());
            }
            visit(line.trim_end_matches(['\n', '\r']));
        }
    })
}

fn hex(bytes: &[u8]) -> String {
    bytes.iter().map(|byte| format!("{byte:02x}")).collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn empty_input_has_the_known_sha256() {
        let hash = sha256_reader(std::io::empty()).unwrap();
        assert_eq!(hash.size_bytes, 0);
        assert_eq!(
            hash.sha256,
            "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855"
        );
    }

    #[test]
    fn inputs_longer_than_one_chunk_hash_like_a_single_update() {
        let bytes = (0..(CHUNK_BYTES * 2 + 17))
            .map(|index| index as u8)
            .collect::<Vec<_>>();
        let streamed = sha256_reader(bytes.as_slice()).unwrap();
        assert_eq!(streamed.size_bytes, bytes.len() as u64);
        assert_eq!(streamed.sha256, hex(&Sha256::digest(&bytes)));
    }
}
