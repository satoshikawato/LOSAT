//! The id of a query or a subject in the HSP records (`HspRecord::query_id`, `subject_id`), as
//! the tabular reports write it: its bytes.
//!
//! NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:214-217
//! NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:252-255
//! NCBI's id is a `std::string` (`s_GetSeqIdListString`, c++/src/objtools/align_format/tabular.cpp:165-195), whose bytes the
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

impl fmt::Debug for Id {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        fmt::Debug::fmt(&DebugBytes(self.as_bytes()), f)
    }
}

/// Bytes as a quoted string: as `str`'s `Debug` writes them where they are UTF-8, otherwise the
/// UTF-8 parts as `str`'s `Debug` writes them and each other byte as `\xNN` (upper-case hex, as
/// std writes the bytes of an `OsStr` on Unix that are not UTF-8). The `Debug` of `Id` and of
/// `Error`'s message.
pub(super) struct DebugBytes<'a>(pub(super) &'a [u8]);

impl fmt::Debug for DebugBytes<'_> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if let Ok(text) = std::str::from_utf8(self.0) {
            return fmt::Debug::fmt(text, f);
        }
        f.write_char('"')?;
        for chunk in self.0.utf8_chunks() {
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
#[path = "id_tests.rs"]
mod tests;
