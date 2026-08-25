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

#[test]
fn cloud_writer_completes_multipart_upload() {
    let url = Url::parse("memory:///multipart-output.tsv").unwrap();
    let mut writer = CloudWriter::new(&url).unwrap();
    writer.write_all(&vec![b'x'; 11 * 1024 * 1024]).unwrap();
    writer.finish().unwrap();
}
