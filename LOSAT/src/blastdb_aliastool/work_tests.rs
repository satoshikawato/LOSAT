use super::*;

#[test]
fn split_names_keeps_empty_pieces_of_a_text_with_names() {
    // NCBI 2.17.0: `-dblist 'in/p4  in/q4'` writes `DBLIST "in/p4" "" "in/q4" ` (probe D08 of
    // row RG-D2a) and `-dblist ' in/p4'` writes `DBLIST "" "in/p4" ` (D17); `-vdblist ''` writes
    // nothing (V04).
    assert_eq!(split_names("in/p4  in/q4"), ["in/p4", "", "in/q4"]);
    assert_eq!(split_names(" in/p4"), ["", "in/p4"]);
    assert_eq!(split_names("a "), ["a", ""]);
    assert!(split_names("").is_empty());
}

#[test]
fn gi_ids_skip_comments_and_lines_that_do_not_convert() {
    // NCBI 2.17.0 (probe G03): ` 102`, `103 `, `abc`, `12x` and an Int8 overflow are skipped;
    // a negative, a zero and a 33-bit GI are kept, in the file's order.
    let text = b"101\n#102\n\n 102\n103 \nabc\n12x\n-5\n0\n4294967297\n99999999999999999999\n";
    assert_eq!(gi_ids(text), [101, -5, 0, 4_294_967_297]);
}

#[cfg(unix)]
#[test]
fn gi_ids_keep_a_carriage_return_on_unix() {
    // NCBI 2.17.0 on Linux (probe G06): a list with CRLF ends converts no GI.
    assert!(gi_ids(b"101\r\n102\r\n").is_empty());
}

#[test]
fn exclude_model_is_the_low_bit_of_the_int() {
    assert!(exclude_model(1));
    assert!(exclude_model(3));
    assert!(exclude_model(-1));
    assert!(!exclude_model(0));
    assert!(!exclude_model(2));
}
