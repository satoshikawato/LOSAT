use super::*;

#[test]
fn reverse_complement_keeps_case_and_ambiguity() {
    let mut letters = b"ACGTMRWSYKVHDBN-acgtu".to_vec();
    reverse_complement(&mut letters);
    assert_eq!(letters, b"aacgt-NVHDBMRSWYKACGT");
}

#[test]
fn crc32_is_zlib_s() {
    // The CRC-32 of "123456789" is the standard check value 0xCBF43926.
    assert_eq!(crc32_zip(b"123456789"), 0xCBF4_3926);
    assert_eq!(sequence_hash(b"1234\n56789"), !0xCBF4_3926);
    assert_eq!(hash_text(0), "0");
    assert_eq!(hash_text(0x09E5_E1FA), "0X9E5E1FA");
}

#[test]
fn masks_are_lower_case_inside_the_letters() {
    let masked = |length: usize, masks: &[(u32, u32)], range| {
        let mut letters = vec![b'A'; length];
        apply_masks(&mut letters, masks, range);
        String::from_utf8(letters).unwrap()
    };
    assert_eq!(masked(10, &[(2, 4), (8, 10)], None), "AAaaAAAAaa");
    assert_eq!(masked(6, &[(2, 4), (8, 10)], Some((3, 8))), "aAAAAa");
    assert_eq!(masked(5, &[(2, 4)], Some((4, 8))), "AAAAA");
    assert_eq!(masked(4, &[(2, 100), (3, 2), (1, 1)], None), "AAaa");
}
