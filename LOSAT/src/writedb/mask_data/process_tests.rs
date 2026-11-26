use super::*;

#[test]
fn spaces_and_commas_are_split_as_ncbi_does() {
    assert_eq!(truncate_spaces(b" \t a b \n"), b"a b");
    assert_eq!(truncate_spaces(" \u{a0}".as_bytes()), "\u{a0}".as_bytes());
    assert_eq!(truncate_spaces(b"  "), b"");
    assert!(split_commas(b"").is_empty());
    assert_eq!(split_commas(b","), [b"".as_slice(), b"".as_slice()]);
    assert_eq!(split_commas(b"a,"), [b"a".as_slice(), b"".as_slice()]);
    assert_eq!(
        split_commas(b"a,,b"),
        [b"a".as_slice(), b"".as_slice(), b"b".as_slice()]
    );
    assert_eq!(split_commas(b"a, b"), [b"a".as_slice(), b" b".as_slice()]);
}

/// The answers of makeblastdb 2.17.0 (`-mask_id` with one or two mask files).
#[test]
fn numeric_ids_are_those_that_stringtonumeric_converts() {
    for numeric in [
        "0",
        "5",
        "007",
        "+5",
        "-0",
        "-5",
        "2147483647",
        "-2147483648",
    ] {
        assert!(is_numeric_int4(numeric.as_bytes()), "{numeric}");
    }
    for not_numeric in [
        "",
        "+",
        "-",
        " 5",
        "5 ",
        "5a",
        "a5",
        "1e3",
        "0x10",
        "٣",
        "+-5",
        "2147483648",
        "-2147483649",
        "99999999999999999999",
        "9223372036854775808",
    ] {
        assert!(!is_numeric_int4(not_numeric.as_bytes()), "{not_numeric}");
    }
}
