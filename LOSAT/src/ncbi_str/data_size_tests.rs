use super::*;

/// Values the source makes certain (decimal digits scaled by the suffix, 'K' = 1000, "KiB" = 1024).
#[test]
fn parsed_values() {
    let cases: &[(&str, u64)] = &[
        ("0", 0),
        ("300", 300),
        ("+5", 5),
        ("00010", 10),
        ("300B", 300),
        ("1KB", 1000),
        ("1K", 1000),
        ("1k", 1000),
        ("1 KB", 1000),
        ("1kb", 1000),
        ("1KiB", 1024),
        ("1kib", 1024),
        ("1 KiB", 1024),
        ("1MiB", 1048576),
        ("1MB", 1_000_000),
        ("3GB", 3_000_000_000),
        ("4GB", 4_000_000_000),
        ("4GiB", 4 << 30),
        ("1TB", 1_000_000_000_000),
        ("1TiB", 1 << 40),
        ("1PB", 1_000_000_000_000_000),
        ("1PiB", 1 << 50),
        ("1EB", 1_000_000_000_000_000_000),
        ("1EiB", 1 << 60),
        ("15EiB", 15 << 60),
        ("4294967295", 4294967295),
        ("18446744073709551615", u64::MAX),
        ("18446744073709551.615K", u64::MAX),
        // a fraction is rounded half up on its first dropped digit
        ("1.4", 1),
        ("1.5", 2),
        ("0.5", 1),
        ("0.4", 0),
        ("1.5KB", 1500),
        ("1.5KiB", 1536),
        ("1.5MiB", 1572864),
        ("2.5GiB", 2684354560),
        ("0.0005K", 1),
        ("0.0004999K", 0),
        ("1.2345K", 1235),
        ("1.2344K", 1234),
        ("10239.5", 10240),
        ("10239.4", 10239),
        ("18446744073709551615.4", u64::MAX),
    ];
    for (input, want) in cases {
        assert_eq!(
            string_to_uint8_data_size(input),
            Ok(*want),
            "input {:?}",
            input
        );
    }
}

#[test]
fn overflow_errors() {
    // ERANGE: m_Pos is the index of the digit that overflows
    let e = string_to_uint8_data_size("18446744073709551616").unwrap_err();
    assert_eq!(
        e.message,
        "Cannot convert string '18446744073709551616' to Uint8"
    );
    assert_eq!(e.pos, 19);
    assert_eq!(
        string_to_uint8_data_size("18446744073709551615.5")
            .unwrap_err()
            .pos,
        20
    );
    assert_eq!(string_to_uint8_data_size("16EiB").unwrap_err().pos, 19);
    assert_eq!(
        string_to_uint8_data_size("18446744073709552K")
            .unwrap_err()
            .pos,
        19
    );
}

#[test]
fn error_text_and_position() {
    let e = string_to_uint8_data_size("abc").unwrap_err();
    assert_eq!(e.message, "Cannot convert string 'abc' to Uint8");
    assert_eq!(e.pos, 0);
    // '-' gives m_Pos of the sign; trailing garbage gives m_Pos of the first unconsumed byte
    assert_eq!(string_to_uint8_data_size("-5").unwrap_err().pos, 0);
    assert_eq!(string_to_uint8_data_size("1KBB").unwrap_err().pos, 3);
    assert_eq!(string_to_uint8_data_size("").unwrap_err().pos, 0);
}
