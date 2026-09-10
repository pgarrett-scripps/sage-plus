use super::*;
use object_store::memory::InMemory;
use std::io::{BufRead, BufReader, Write};
use std::net::TcpListener;
use std::time::{Duration, Instant};

#[test]
fn storage_operations_preserve_bytes_and_paths() {
    let object = CloudObject {
        store: Arc::new(InMemory::new()),
        path: "sample.bin".into(),
        base: "s3://fixture".into(),
    };
    object.upload_bytes(b"0123456789".to_vec()).unwrap();
    assert_eq!(object.len().unwrap(), 10);
    assert_eq!(object.range(2..6).unwrap(), b"2345");
    assert_eq!(object.range(7..).unwrap(), b"789");
    let directory = CloudObject {
        store: object.store.clone(),
        path: "".into(),
        base: object.base.clone(),
    };
    assert_eq!(
        directory.children().unwrap(),
        vec!["s3://fixture/sample.bin"]
    );

    let dir = tempfile::tempdir().unwrap();
    let destination = dir.path().join("nested/download.bin");
    object.download_to(&destination).unwrap();
    assert_eq!(std::fs::read(&destination).unwrap(), b"0123456789");
    std::fs::write(&destination, b"replacement").unwrap();
    object.upload_from(&destination).unwrap();
    assert_eq!(object.range(..).unwrap(), b"replacement");
}

#[test]
fn s3_xml_listing_and_http_ranges_remain_compatible() {
    let listener = TcpListener::bind("127.0.0.1:0").unwrap();
    listener.set_nonblocking(true).unwrap();
    let endpoint = format!("http://{}", listener.local_addr().unwrap());
    let server = std::thread::spawn(move || {
        let deadline = Instant::now() + Duration::from_secs(20);
        let mut requests = Vec::new();
        while requests.len() < 4 {
            let (mut stream, _) = match listener.accept() {
                Ok(pair) => pair,
                Err(error) if error.kind() == std::io::ErrorKind::WouldBlock => {
                    assert!(Instant::now() < deadline, "missing S3 fixture request");
                    std::thread::sleep(Duration::from_millis(5));
                    continue;
                }
                Err(error) => panic!("S3 fixture accept failed: {error}"),
            };
            stream
                .set_read_timeout(Some(Duration::from_secs(5)))
                .unwrap();
            let mut reader = BufReader::new(stream.try_clone().unwrap());
            let mut request = String::new();
            loop {
                let mut line = String::new();
                assert!(reader.read_line(&mut line).unwrap() > 0);
                if line == "\r\n" {
                    break;
                }
                request.push_str(&line);
            }
            let listing = request.contains("list-type=2");
            let head = request.starts_with("HEAD ");
            let body = if listing {
                "<?xml version=\"1.0\"?><ListBucketResult xmlns=\"http://s3.amazonaws.com/doc/2006-03-01/\"><Name>fixture</Name><Prefix>sample.bin</Prefix><KeyCount>1</KeyCount><MaxKeys>1000</MaxKeys><IsTruncated>false</IsTruncated><Contents><Key>sample.bin</Key><LastModified>2026-09-10T00:00:00.000Z</LastModified><ETag>\"fixture\"</ETag><Size>10</Size><StorageClass>STANDARD</StorageClass></Contents></ListBucketResult>"
            } else if head {
                ""
            } else {
                assert!(request.to_lowercase().contains("range: bytes=2-5"));
                "2345"
            };
            let status = if !listing && !head {
                "206 Partial Content"
            } else {
                "200 OK"
            };
            let content_range = if !listing && !head {
                "Content-Range: bytes 2-5/10\r\n"
            } else {
                ""
            };
            let length = if head { 10 } else { body.len() };
            write!(stream, "HTTP/1.1 {status}\r\nContent-Length: {length}\r\nLast-Modified: Thu, 10 Sep 2026 00:00:00 GMT\r\nETag: \"fixture\"\r\n{content_range}Connection: close\r\n\r\n{body}").unwrap();
            requests.push(request);
        }
        requests
    });

    let store = object_store::aws::AmazonS3Builder::new()
        .with_bucket_name("fixture")
        .with_region("us-east-1")
        .with_endpoint(endpoint)
        .with_allow_http(true)
        .with_skip_signature(true)
        .build()
        .unwrap();
    let object = CloudObject {
        store: Arc::new(store),
        path: "sample.bin".into(),
        base: "s3://fixture".into(),
    };
    assert_eq!(object.len().unwrap(), 10);
    assert_eq!(object.range(2..6).unwrap(), b"2345");
    let directory = CloudObject {
        store: object.store.clone(),
        path: "".into(),
        base: object.base.clone(),
    };
    assert_eq!(
        directory.children().unwrap(),
        vec!["s3://fixture/sample.bin"]
    );
    let requests = server.join().unwrap();
    assert_eq!(
        requests.iter().filter(|r| r.starts_with("HEAD ")).count(),
        2
    );
}
