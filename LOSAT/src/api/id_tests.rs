use super::Id;
use std::collections::{BTreeMap, HashMap, HashSet};

#[test]
fn the_text_of_an_id_is_its_utf8_and_u_fffd_for_the_other_bytes() {
    let latin1 = Id::new(b"q\xE9");
    assert!(latin1.to_str().is_err());
    assert_eq!(latin1.to_str().unwrap_err().valid_up_to(), 1);
    assert_eq!(latin1.to_string_lossy(), "q\u{FFFD}");
    // One U+FFFD per maximal ill-formed subpart: a truncated `€` (E2 82) is one subpart, a
    // lone continuation byte (80) and FF are one each.
    let truncated = Id::new(b"a\xE2\x82b\x80\xFF");
    assert_eq!(truncated.to_string_lossy(), "a\u{FFFD}b\u{FFFD}\u{FFFD}");
    let utf8 = Id::new("qé€".as_bytes());
    assert_eq!(utf8.to_str(), Ok("qé€"));
    assert!(matches!(
        utf8.to_string_lossy(),
        std::borrow::Cow::Borrowed("qé€")
    ));
}

#[test]
fn the_debug_form_of_an_id_escapes_the_bytes_that_are_not_utf8() {
    assert_eq!(format!("{:?}", Id::new(b"q\xE9")), r#""q\xE9""#);
    assert_eq!(
        format!("{:?}", Id::new(b"a\"b\x80\xFF\tc")),
        r#""a\"b\x80\xFF\tc""#
    );
    // An id that is UTF-8 is written as its `str` (as the `Arc<str>` field before it was).
    for text in ["seq1", "qé", "a'b\"c\\d\n", "\u{301}x"] {
        assert_eq!(
            format!("{:?}", Id::new(text.as_bytes())),
            format!("{text:?}")
        );
    }
    // Inside an id that is not UTF-8, the UTF-8 parts are escaped as `str` escapes them (a
    // combining accent after the byte, too).
    let mixed = ["é\n".as_bytes(), b"\xC0", "\u{301}".as_bytes()].concat();
    assert_eq!(format!("{:?}", Id::new(&mixed)), r#""é\n\xC0\u{301}""#);
}

#[test]
fn ids_compare_by_their_bytes() {
    let id = Id::new(b"seq1");
    let (string, slice, vec): (String, &[u8], Vec<u8>) = ("seq1".into(), b"seq1", b"seq1".into());
    assert!(id == "seq1");
    assert!(id == *"seq1");
    assert!(id == string);
    assert!(id == b"seq1");
    assert!(id == *b"seq1");
    assert!(id == slice);
    assert!(id == *slice);
    assert!(id == vec);
    assert!("seq1" == id);
    assert!(*"seq1" == id);
    assert!(*slice == id);
    assert!(id != "seq2");
    assert!("seq2" != id);
    assert_eq!(id, id.clone());
    assert_eq!(id.as_ref(), b"seq1");
    // The text with U+FFFD is not the id.
    let latin1 = Id::new(b"q\xE9");
    assert!(latin1 != "q\u{FFFD}");
    assert!("q\u{FFFD}" != latin1);
    assert!(latin1 == b"q\xE9");
    assert!(latin1 != Id::new("q\u{FFFD}".as_bytes()));
    // Ordered as the bytes: 0xE9 after every ASCII byte; a prefix first.
    let mut ids = [Id::new(b"q\xE9"), Id::new(b"qz"), Id::new(b"q")];
    ids.sort();
    let sorted: Vec<&[u8]> = ids.iter().map(Id::as_bytes).collect();
    assert_eq!(sorted, [&b"q"[..], b"qz", b"q\xE9"]);
    let distinct: HashSet<Id> = [
        Id::new(b"q\xE9"),
        Id::new(b"q\xE9"),
        Id::new("q\u{FFFD}".as_bytes()),
    ]
    .into_iter()
    .collect();
    assert_eq!(distinct.len(), 2);
}

#[test]
fn maps_keyed_by_ids_are_looked_up_by_bytes() {
    let ids = [Id::new(b"q\xE9"), Id::new(b"seq1")];
    let hashed: HashMap<Id, usize> = ids.iter().cloned().zip(0..).collect();
    let ordered: BTreeMap<Id, usize> = ids.iter().cloned().zip(0..).collect();
    for (key, value) in [(&b"q\xE9"[..], 0), ("seq1".as_bytes(), 1)] {
        assert_eq!(hashed.get(key), Some(&value));
        assert_eq!(ordered.get(key), Some(&value));
    }
    assert_eq!(hashed.get("q\u{FFFD}".as_bytes()), None);
    assert_eq!(ordered.get("q\u{FFFD}".as_bytes()), None);
}
