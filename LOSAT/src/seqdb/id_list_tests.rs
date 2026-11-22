use super::*;

fn binary(words: &[u32]) -> Vec<u8> {
    words.iter().flat_map(|word| word.to_be_bytes()).collect()
}

#[test]
fn a_text_list_loses_its_last_number_without_a_separator_after_it() {
    assert_eq!(read_gi_list(b"11\n22").unwrap(), vec![11]);
    assert_eq!(read_gi_list(b"11\n22\n").unwrap(), vec![11, 22]);
    assert_eq!(read_gi_list(b"11 0 \r\n\n33\r").unwrap(), vec![11, 33]);
}

#[test]
fn a_text_gi_wraps_as_an_unsigned_32_bit_number() {
    assert_eq!(read_gi_list(b"4294967296\n").unwrap(), Vec::<u32>::new());
    assert_eq!(read_gi_list(b"4294967297\n").unwrap(), vec![1]);
}

#[test]
fn the_first_byte_decides_text_or_binary() {
    assert_eq!(read_gi_list(b"").unwrap_err(), "Specified file is empty.");
    assert_eq!(
        read_gi_list(b" 11\n").unwrap_err(),
        "Specified file is not a valid GI/TI list."
    );
    assert_eq!(
        read_gi_list(b"#c\n11\n").unwrap_err(),
        "Invalid byte in textGI list [35]."
    );
    assert_eq!(
        read_gi_list(b"11\t22\n").unwrap_err(),
        "Invalid byte in textGI list [9]."
    );
    assert_eq!(
        read_gi_list(&[0xFF, 0xFF, 0xFF, 0xFF, 0, 0, 0]).unwrap_err(),
        "Specified file is not a valid GI/TI list."
    );
}

#[test]
fn a_binary_gi_list_needs_its_marker_and_count() {
    assert_eq!(
        read_gi_list(&binary(&[0xFFFF_FFFF, 2, 7, 5])).unwrap(),
        vec![7, 5]
    );
    // A length that is not a whole number of words is rounded down.
    let mut odd = binary(&[0xFFFF_FFFF, 1, 9]);
    odd.push(1);
    assert_eq!(read_gi_list(&odd).unwrap(), vec![9]);
    for bad in [
        binary(&[0xFFFF_FFFF, 3, 7, 5]),
        binary(&[0xFFFF_FFFD, 1, 7]),
    ] {
        assert_eq!(
            read_gi_list(&bad).unwrap_err(),
            "Specified file is not a valid binary GI file."
        );
    }
}

#[test]
fn a_ti_list_is_signed_and_widened_past_int8() {
    assert_eq!(
        read_ti_list(b"9223372036854775807\n9223372036854775808\n5\n").unwrap(),
        vec![
            ListNumber::Value(i64::MAX),
            ListNumber::Past,
            ListNumber::Value(5)
        ]
    );
    assert_eq!(
        read_ti_list(&binary(&[0xFFFF_FFFD, 2, 0xFFFF_FFFF, 3])).unwrap(),
        vec![ListNumber::Value(-1), ListNumber::Value(3)]
    );
    assert_eq!(
        read_ti_list(&binary(&[0xFFFF_FFFC, 1, 1, 2])).unwrap(),
        vec![ListNumber::Value((1 << 32) | 2)]
    );
    for bad in [
        binary(&[0xFFFF_FFFF, 1, 3]),
        binary(&[0xFFFF_FFFD, 2, 3]),
        binary(&[0xFFFF_FFFC, 1, 1, 2, 3]),
    ] {
        assert_eq!(
            read_ti_list(&bad).unwrap_err(),
            "Specified file is not a valid binary GI or TI file."
        );
    }
    assert_eq!(
        read_ti_list(b"1a\n").unwrap_err(),
        "Invalid byte in textTI list [97]."
    );
}

#[test]
fn a_taxid_list_is_a_set_of_int4_numbers() {
    let set = read_tax_id_list(b"9606\n10090\n9606\n2147483648\n").unwrap();
    assert_eq!(
        set.into_iter().collect::<Vec<_>>(),
        vec![
            ListNumber::Value(9606),
            ListNumber::Value(10090),
            ListNumber::Past
        ]
    );
    assert_eq!(
        read_tax_id_list(&binary(&[0xFFFF_FFFF, 1, 562]))
            .unwrap()
            .into_iter()
            .collect::<Vec<_>>(),
        vec![ListNumber::Value(562)]
    );
    assert_eq!(
        read_tax_id_list(&binary(&[0xFFFF_FFFF, 2, 562])).unwrap_err(),
        "Specified file is not a valid binary Tax Id List file."
    );
    assert_eq!(
        read_tax_id_list(b"9606,\n").unwrap_err(),
        "Invalid byte in textTAXID list [44]."
    );
}
