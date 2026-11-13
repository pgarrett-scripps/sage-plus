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
