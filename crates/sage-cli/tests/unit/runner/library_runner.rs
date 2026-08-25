use super::source_basename;

#[test]
fn source_basename_normalizes_platform_paths() {
    assert_eq!(
        source_basename("runs/One.mzML").as_deref(),
        Some("one.mzml")
    );
    assert_eq!(
        source_basename(r"C:\\runs\\Two.RAW").as_deref(),
        Some("two.raw")
    );
    assert_eq!(source_basename(""), None);
}
