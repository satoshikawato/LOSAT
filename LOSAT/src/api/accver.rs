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
//! `AccVer` keeps those bytes; its text (`to_str`, `to_string_lossy`, `Display`) is a view of
//! them.
#![warn(missing_docs, unused)]

use std::borrow::Cow;
use std::fmt::{self, Write as _};
use std::str::Utf8Error;
use std::sync::Arc;

/// The id of a query or a subject as the tabular reports write it in the columns `qaccver` and
/// `saccver` (outfmt 7's `query acc.ver` and `subject acc.ver`): its bytes, exactly, also those
/// that are not UTF-8 (a title read from a Latin-1 file, for example).
///
/// A clone shares the bytes: the HSPs of one record share one copy. Two ids are equal, ordered
/// and hashed by their bytes, and an id equals a `str` or a `[u8]` of the same bytes
/// (`hsp.query_id == "seq1"`). As text, [`to_str`](Self::to_str) gives the id where it is UTF-8,
/// and [`to_string_lossy`](Self::to_string_lossy) and `Display` give U+FFFD for each sequence of
/// bytes that is not, as `String::from_utf8_lossy` does. `Debug` writes it as a quoted string
/// with those bytes as `\xNN` escapes (`"q\xE9"`), and as `str` writes it where it is UTF-8.
#[derive(Clone, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct AccVer {
    bytes: Arc<[u8]>,
}

impl AccVer {
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
        std::str::from_utf8(&self.bytes)
    }

    /// The id as text, with U+FFFD for each sequence of bytes that is not UTF-8
    /// (`String::from_utf8_lossy`); borrowed where every byte is UTF-8.
    pub fn to_string_lossy(&self) -> Cow<'_, str> {
        String::from_utf8_lossy(&self.bytes)
    }
}

impl PartialEq<str> for AccVer {
    fn eq(&self, other: &str) -> bool {
        *self.bytes == *other.as_bytes()
    }
}

impl PartialEq<&str> for AccVer {
    fn eq(&self, other: &&str) -> bool {
        *self.bytes == *other.as_bytes()
    }
}

impl PartialEq<[u8]> for AccVer {
    fn eq(&self, other: &[u8]) -> bool {
        *self.bytes == *other
    }
}

impl AsRef<[u8]> for AccVer {
    fn as_ref(&self) -> &[u8] {
        &self.bytes
    }
}

/// The text of `to_string_lossy`, padded as a `str` (width, fill, alignment, precision).
impl fmt::Display for AccVer {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.pad(&self.to_string_lossy())
    }
}

/// A quoted string: the UTF-8 parts as `str`'s `Debug` writes them, each other byte as `\xNN`
/// (upper-case hex, as std writes the bytes of an `OsStr` on Unix that are not UTF-8).
impl fmt::Debug for AccVer {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if let Ok(text) = self.to_str() {
            return fmt::Debug::fmt(text, f);
        }
        f.write_char('"')?;
        for chunk in self.bytes.utf8_chunks() {
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
    use super::AccVer;
    use std::collections::HashSet;

    #[test]
    fn the_text_of_an_id_is_its_utf8_and_u_fffd_for_the_other_bytes() {
        let latin1 = AccVer::new(b"q\xE9");
        assert!(latin1.to_str().is_err());
        assert_eq!(latin1.to_str().unwrap_err().valid_up_to(), 1);
        assert_eq!(latin1.to_string_lossy(), "q\u{FFFD}");
        assert_eq!(latin1.to_string(), "q\u{FFFD}");
        // One U+FFFD for each maximal sequence that is not UTF-8: a truncated `€` (E2 82) is one.
        let truncated = AccVer::new(b"a\xE2\x82b\x80\xFF");
        assert_eq!(truncated.to_string(), "a\u{FFFD}b\u{FFFD}\u{FFFD}");
        let utf8 = AccVer::new("qé€".as_bytes());
        assert_eq!(utf8.to_str(), Ok("qé€"));
        assert!(matches!(
            utf8.to_string_lossy(),
            std::borrow::Cow::Borrowed("qé€")
        ));
        assert_eq!(utf8.to_string(), "qé€");
        // Padded as a `str`.
        assert_eq!(format!("[{:>5}]", AccVer::new(b"ab")), "[   ab]");
        assert_eq!(format!("[{:<4}]", latin1), "[q\u{FFFD}  ]");
        assert_eq!(format!("[{:.1}]", utf8), "[q]");
    }

    #[test]
    fn the_debug_form_of_an_id_escapes_the_bytes_that_are_not_utf8() {
        assert_eq!(format!("{:?}", AccVer::new(b"q\xE9")), r#""q\xE9""#);
        assert_eq!(
            format!("{:?}", AccVer::new(b"a\"b\x80\xFF\tc")),
            r#""a\"b\x80\xFF\tc""#
        );
        // An id that is UTF-8 is written as its `str` (as the `Arc<str>` field before it was).
        for text in ["seq1", "qé", "a'b\"c\\d\n", "\u{301}x"] {
            assert_eq!(
                format!("{:?}", AccVer::new(text.as_bytes())),
                format!("{text:?}")
            );
        }
        // Inside an id that is not UTF-8, the UTF-8 parts are escaped as `str` escapes them (a
        // combining accent after the byte, too).
        let mixed = ["é\n".as_bytes(), b"\xC0", "\u{301}".as_bytes()].concat();
        assert_eq!(format!("{:?}", AccVer::new(&mixed)), r#""é\n\xC0\u{301}""#);
    }

    #[test]
    fn ids_compare_by_their_bytes() {
        let id = AccVer::new(b"seq1");
        assert!(id == "seq1");
        assert!(id == *"seq1");
        assert!(id == b"seq1"[..]);
        assert!(id != "seq2");
        assert_eq!(id, id.clone());
        assert_eq!(id.as_ref(), b"seq1");
        // The text with U+FFFD is not the id.
        let latin1 = AccVer::new(b"q\xE9");
        assert!(latin1 != "q\u{FFFD}");
        assert!(latin1 == b"q\xE9"[..]);
        assert!(latin1 != AccVer::new("q\u{FFFD}".as_bytes()));
        // Ordered as the bytes: 0xE9 after every ASCII byte; a prefix first.
        let mut ids = [AccVer::new(b"q\xE9"), AccVer::new(b"qz"), AccVer::new(b"q")];
        ids.sort();
        let sorted: Vec<&[u8]> = ids.iter().map(AccVer::as_bytes).collect();
        assert_eq!(sorted, [&b"q"[..], b"qz", b"q\xE9"]);
        let distinct: HashSet<AccVer> = [
            AccVer::new(b"q\xE9"),
            AccVer::new(b"q\xE9"),
            AccVer::new("q\u{FFFD}".as_bytes()),
        ]
        .into_iter()
        .collect();
        assert_eq!(distinct.len(), 2);
    }
}
