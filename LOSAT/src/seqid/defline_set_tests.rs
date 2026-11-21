use super::*;

// `CObject_id::GetIdType` (Object_id.cpp:93-161): the strings that are integers, and those
// that stay strings (an error of `StringToInt8`, a `+`, leading zeros, `-0`). The conversion
// ends at a 0 byte, so `5\0x` is the integer 5 (the string's size is still 3).
#[test]
fn object_id_strings_are_typed_as_ncbis_get_id_type() {
    let number = |text: &[u8]| object_id_number(&ObjectId::Str(text.to_vec()));
    let integers: [(&[u8], i64); 7] = [
        (b"5", 5),
        (b"0", 0),
        (b"-12", -12),
        (b"9223372036854775807", i64::MAX),
        (b"-9223372036854775808", i64::MIN),
        (b"5\x00x", 5),
        (b"-7\x00", -7),
    ];
    for (text, value) in integers {
        assert_eq!(number(text), Some(value), "{text:?}");
    }
    let strings: [&[u8]; 14] = [
        b"",
        b"-",
        b"+5",
        b"05",
        b"-05",
        b"-0",
        b"00",
        b"0\x00",
        b" 5",
        b"5 ",
        b"5x",
        b"9223372036854775808",
        b"abc",
        b"\xff5",
    ];
    for text in strings {
        assert_eq!(number(text), None, "{text:?}");
    }
    assert_eq!(object_id_number(&ObjectId::Id(-3)), Some(-3));
}
