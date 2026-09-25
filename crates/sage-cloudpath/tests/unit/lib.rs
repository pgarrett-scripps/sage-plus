use super::*;

#[test]
fn filename_gcs() {
    let url = Url::parse("gs://my-bucket/path/to/file.mzML").unwrap();
    assert_eq!(filename(&url), Some("file.mzML"));
}

#[test]
fn filename_azure() {
    let url = Url::parse("az://my-container/path/to/file.mzML").unwrap();
    assert_eq!(filename(&url), Some("file.mzML"));
}

#[test]
fn invalid_remote_read() {
    assert!(read_and_execute("s3://my-bucket", |_| async move { Ok(()) }).is_err())
}

#[test]
fn windows_drive_letter_is_not_a_url() {
    // `Url::parse("C:\\...")` succeeds with `c` as a single-letter
    // scheme, which object_store later rejects with "Unable to recognise
    // URL". `to_url` must treat such inputs as local paths — here the
    // path doesn't exist on this machine, so we expect an IO error
    // rather than a bogus `Ok(Url { scheme: "c", ... })`.
    let backslash = to_url(r"C:\Users\nonexistent\bar.json");
    assert!(
        matches!(backslash, Err(Error::IO(_))),
        "expected IO error for Windows path with backslashes, got {:?}",
        backslash
    );

    let forwardslash = to_url("C:/Users/nonexistent/bar.json");
    assert!(
        matches!(forwardslash, Err(Error::IO(_))),
        "expected IO error for Windows path with forward slashes, got {:?}",
        forwardslash
    );
}

#[test]
fn cloud_urls_still_parse() {
    assert_eq!(to_url("s3://bucket/key").unwrap().scheme(), "s3");
    assert_eq!(to_url("gs://bucket/key").unwrap().scheme(), "gs");
    assert_eq!(to_url("az://container/key").unwrap().scheme(), "az");
    assert_eq!(to_url("https://example.com/key").unwrap().scheme(), "https");
}

#[test]
fn bruker_filenames() {
    let url = Url::parse("file:///data/20251005_sample_a.d/analysis.tdf").unwrap();
    assert_eq!(filename(&url), Some("20251005_sample_a.d"));

    let url = Url::parse("s3://bucket/baz/20251005_sample_a.d/analysis.tdf").unwrap();
    assert_eq!(filename(&url), Some("20251005_sample_a.d"));

    let url = Url::parse("file:///data/baz/20251005_sample_a.mzML").unwrap();
    assert_eq!(filename(&url), Some("20251005_sample_a.mzML"));
}

#[test]
fn gzip_detection() {
    assert!(gzip_heuristic(&Url::parse("file:///file.mzML.gz").unwrap()));
    assert!(gzip_heuristic(
        &Url::parse("s3://bucket/file.mzML.gzip").unwrap()
    ));
    assert!(!gzip_heuristic(&Url::parse("file:///file.mzML").unwrap()));
}

#[cfg(feature = "cloud")]
#[test]
fn cloud_writer_completes_multipart_upload() {
    let url = Url::parse("memory:///multipart-output.tsv").unwrap();
    let mut writer = CloudWriter::new(&url).unwrap();
    writer.write_all(&vec![b'x'; 11 * 1024 * 1024]).unwrap();
    writer.finish().unwrap();
}

#[test]
fn gzip_writer_finishes_complete_round_trips() {
    let directory = tempfile::tempdir().unwrap();
    for payload in [
        Vec::new(),
        b"round trip payload".to_vec(),
        (0..100_000).map(|i| (i % 251) as u8).collect(),
    ] {
        let path = directory.path().join("output.gz");
        let url = Url::from_file_path(&path).unwrap();
        write_bytes_sync(&url, payload.clone()).unwrap();
        let actual = read_and_execute(path.to_str().unwrap(), |mut reader| async move {
            let mut decoded = Vec::new();
            tokio::io::AsyncReadExt::read_to_end(&mut reader, &mut decoded).await?;
            Ok(decoded)
        })
        .unwrap();
        assert_eq!(actual, payload);
    }
}

#[cfg(not(feature = "cloud"))]
#[test]
fn cloud_urls_fail_clearly_without_the_cloud_feature() {
    let url = Url::parse("s3://bucket/results.tsv").unwrap();
    let message = CloudWriter::new(&url).err().unwrap().to_string();
    assert!(message.contains("`s3://`") && message.contains("`cloud` feature"));
    let error = write_bytes_sync(&url, b"x".to_vec()).unwrap_err();
    assert!(matches!(error, Error::CloudDisabled(scheme) if scheme == "s3"));
    let error = read_and_execute("gs://bucket/file.mzML", |_| async move { Ok(()) }).unwrap_err();
    assert!(matches!(error, Error::CloudDisabled(scheme) if scheme == "gs"));
}

#[test]
fn local_writes_replace_existing_files() {
    let directory = tempfile::tempdir().unwrap();
    let path = directory.path().join("nested").join("out.tsv");
    let url = Url::from_file_path(&path).unwrap();
    write_bytes_sync(&url, b"first".to_vec()).unwrap();
    write_bytes_sync(&url, b"second".to_vec()).unwrap();
    assert_eq!(std::fs::read(&path).unwrap(), b"second");
    assert_eq!(
        std::fs::read_dir(path.parent().unwrap()).unwrap().count(),
        1
    );
}
