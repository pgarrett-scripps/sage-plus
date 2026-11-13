use super::scan_number;

#[test]
fn scan_numbers_from_native_ids() {
    assert_eq!(
        scan_number("controllerType=0 controllerNumber=1 scan=35020"),
        35020
    );
    assert_eq!(scan_number("1234"), 1234);
    assert_eq!(scan_number("index=7"), 0);
}
