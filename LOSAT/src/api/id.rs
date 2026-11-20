//! The id of a query or a subject in the HSP records (`HspRecord::query_id`, `subject_id`), as
//! the tabular reports write it: its bytes.
//!
//! NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:214-217
//! ```c++
//! void CBlastTabularInfo::x_PrintQueryAccessionVersion()
//! {
//!     m_Ostream << s_GetSeqIdListString(m_QueryId, eAccVersion);
//! }
//! ```
//! NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:252-255
//! ```c++
//! void CBlastTabularInfo::x_PrintSubjectAccessionVersion(void)
//! {
//!     m_Ostream << s_GetSeqIdListString(m_SubjectId, eAccVersion);
//! }
//! ```
//! NCBI's id is a `std::string` (`s_GetSeqIdListString`, tabular.cpp:165-195), whose bytes the
//! stream writes as they are: a byte of a title that is not UTF-8 reaches the report as itself.
//! `Id` keeps those bytes; its text (`to_str`, `to_string_lossy`) is a view of them.
#![warn(missing_docs, unused)]

use std::borrow::{Borrow, Cow};
use std::cmp::Ordering;
use std::fmt::{self, Write as _};
use std::hash::{Hash, Hasher};
use std::str::Utf8Error;
use std::sync::Arc;

/// The id of a query or a subject as the tabular reports write it: its bytes, exactly, also those
/// that are not UTF-8 (a title read from a Latin-1 file, for example). Today it holds the columns
/// `qaccver` and `saccver` (outfmt 7's `query acc.ver` and `subject acc.ver`) in
/// [`HspRecord::query_id`](super::HspRecord::query_id) and
/// [`subject_id`](super::HspRecord::subject_id); a subject for which NCBI writes no `saccver`
/// (its lookup fails) has the id that `subject_id` describes.
///
/// Cloning is cheap; it does not copy the bytes.
///
/// Two ids are equal, ordered and hashed as their bytes ([`as_bytes`](Self::as_bytes)). An id
/// equals a `str`, a `String`, a `[u8]`, a byte array or a `Vec<u8>` of the same bytes
/// (`hsp.query_id == "seq1"`, `"seq1" == hsp.query_id`, `hsp.query_id == b"seq1"`). It is
/// `Borrow<[u8]>`, so a `HashMap` or `BTreeMap` keyed by ids is looked up by a `&[u8]`
/// (`map.get(name.as_bytes())`). It will never be `Borrow<str>`: a `str` hashes otherwise than
/// its bytes as a `[u8]` do, so a lookup by a `&str` would miss.
///
/// As text, [`to_str`](Self::to_str) gives the id where it is UTF-8, and
/// [`to_string_lossy`](Self::to_string_lossy) gives one U+FFFD per maximal ill-formed subpart, as
/// `String::from_utf8_lossy` and Python's `decode("utf-8", "replace")` do. `Id` has no `Display`
/// on purpose: text with U+FFFD is not the id, so it is asked for by name. `Debug` writes the id
/// as `str` does where it is UTF-8, and otherwise as a quoted string with each byte that is not
/// UTF-8 as `\xNN` (`"q\xE9"`).
#[derive(Clone)]
pub struct Id {
    bytes: Arc<[u8]>,
}

impl Id {
    /// The id of these bytes.
    pub(crate) fn new(bytes: &[u8]) -> Self {
        Self {
            bytes: Arc::from(bytes),
        }
    }

    /// The bytes, as the reports write them.
    pub fn as_bytes(&self) -> &[u8] {
        &self.bytes
    }

    /// The id as text, or the error that says where its bytes stop being UTF-8 (as
    /// `CStr::to_str`).
    pub fn to_str(&self) -> Result<&str, Utf8Error> {
        std::str::from_utf8(self.as_bytes())
    }

    /// The id as text, with one U+FFFD per maximal ill-formed subpart (`String::from_utf8_lossy`);
    /// borrowed where every byte is UTF-8.
    pub fn to_string_lossy(&self) -> Cow<'_, str> {
        String::from_utf8_lossy(self.as_bytes())
    }
}

// Equality, order and hash are those of `as_bytes()`, written here rather than derived, so that
// they stay those of the bytes whatever holds them: `Borrow<[u8]>` needs it.
impl PartialEq for Id {
    fn eq(&self, other: &Id) -> bool {
        self.as_bytes() == other.as_bytes()
    }
}

impl Eq for Id {}

impl PartialOrd for Id {
    fn partial_cmp(&self, other: &Id) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for Id {
    fn cmp(&self, other: &Id) -> Ordering {
        self.as_bytes().cmp(other.as_bytes())
    }
}

impl Hash for Id {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.as_bytes().hash(state);
    }
}

impl Borrow<[u8]> for Id {
    fn borrow(&self) -> &[u8] {
        self.as_bytes()
    }
}

impl AsRef<[u8]> for Id {
    fn as_ref(&self) -> &[u8] {
        self.as_bytes()
    }
}

/// `Id == other` for each type: by the bytes of `other`.
macro_rules! id_equals {
    ($([$($generics:tt)*] $other:ty, $value:ident => $bytes:expr;)*) => {$(
        impl<$($generics)*> PartialEq<$other> for Id {
            fn eq(&self, $value: &$other) -> bool {
                self.as_bytes() == $bytes
            }
        }
    )*};
}

id_equals! {
    [] str, other => other.as_bytes();
    [] &str, other => other.as_bytes();
    [] String, other => other.as_bytes();
    [] [u8], other => other;
    [] &[u8], other => *other;
    [const N: usize] [u8; N], other => &other[..];
    [const N: usize] &[u8; N], other => &other[..];
    [] Vec<u8>, other => &other[..];
}

impl PartialEq<Id> for str {
    fn eq(&self, other: &Id) -> bool {
        self.as_bytes() == other.as_bytes()
    }
}

impl PartialEq<Id> for &str {
    fn eq(&self, other: &Id) -> bool {
        self.as_bytes() == other.as_bytes()
    }
}

impl PartialEq<Id> for [u8] {
    fn eq(&self, other: &Id) -> bool {
        self == other.as_bytes()
    }
}

/// A quoted string: the UTF-8 parts as `str`'s `Debug` writes them, each other byte as `\xNN`
/// (upper-case hex, as std writes the bytes of an `OsStr` on Unix that are not UTF-8).
impl fmt::Debug for Id {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if let Ok(text) = self.to_str() {
            return fmt::Debug::fmt(text, f);
        }
        f.write_char('"')?;
        for chunk in self.as_bytes().utf8_chunks() {
            let valid = format!("{:?}", chunk.valid());
            f.write_str(&valid[1..valid.len() - 1])?;
            for byte in chunk.invalid() {
                write!(f, "\\x{byte:02X}")?;
            }
        }
        f.write_char('"')
    }
}

#[cfg(test)]
mod tests {
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
        let (string, slice, vec): (String, &[u8], Vec<u8>) =
            ("seq1".into(), b"seq1", b"seq1".into());
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
}
