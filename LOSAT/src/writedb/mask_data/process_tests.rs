use super::*;

#[test]
fn spaces_and_commas_are_split_as_ncbi_does() {
    assert_eq!(truncate_spaces_str(" \t a b \n"), "a b");
    assert_eq!(truncate_spaces_str(" \u{a0}"), "\u{a0}");
    assert_eq!(truncate_spaces_str("  "), "");
    assert!(split_commas("").is_empty());
    assert_eq!(split_commas(","), ["", ""]);
    assert_eq!(split_commas("a,"), ["a", ""]);
    assert_eq!(split_commas("a,,b"), ["a", "", "b"]);
    assert_eq!(split_commas("a, b"), ["a", " b"]);
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
        assert!(is_numeric_int4(numeric), "{numeric}");
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
        assert!(!is_numeric_int4(not_numeric), "{not_numeric}");
    }
}
