//! NCBI's repair of a VisibleString as its object streams read and write it: the method
//! `eFNP_ReplaceAndWarn` (the default of `SERIAL_WRONG_CHARS_READ` and of
//! `SERIAL_WRONG_CHARS_WRITE`) stores or writes `#` for each byte outside `' '..='~'` and posts
//! an error for it with the stream's stack trace, its position and the bytes being read or the
//! string being written (`ReplaceVisibleChar`). The readers of ASN.1 text
//! (`asn_text_reader::Reader::string`, which also cuts the string into NCBI's chunks) and of
//! ASN.1 binary (`fix_visible_chars`, the whole string), the writer of definition lines
//! (`writedb::header`) and the ASN.1 text and binary encoders of the serialized documents
//! (`asn_text`, `asn_ber`, `Encoder::take_bad_chars`) note the bad bytes of each string as one
//! [`BadChars`], which carries the stream's trace (`object_stack::ObjectStack`) and builds
//! NCBI's messages ([`BadChars::message`], the one producer of their text).
//!
//! UTF-8 is not let through (`SERIAL_ALLOW_UTF8_IN_VISIBLESTRING_ON_READING` is 0).
// NCBI reference (598d8ae6): c++/include/serial/objistr.hpp:1160-1163
// ```c++
// inline
// bool GoodVisibleChar(char c);
//
// #define SERIAL_ALLOW_UTF8_IN_VISIBLESTRING_ON_READING 0
// ```
// NCBI reference (598d8ae6): c++/src/serial/objistr.cpp:307-311
// ```c++
// NCBI_PARAM_ENUM_DECL(EFixNonPrint, SERIAL, WRONG_CHARS_READ);
// NCBI_PARAM_ENUM_DEF(EFixNonPrint, SERIAL, WRONG_CHARS_READ, eFNP_ReplaceAndWarn);
// ```
// NCBI reference (598d8ae6): c++/src/serial/objostr.cpp:277-278
// ```c++
// NCBI_PARAM_ENUM_DECL(EFixNonPrint, SERIAL, WRONG_CHARS_WRITE);
// NCBI_PARAM_ENUM_DEF(EFixNonPrint, SERIAL, WRONG_CHARS_WRITE, eFNP_ReplaceAndWarn);
// ```
// NCBI reference (598d8ae6): c++/include/serial/objistr.hpp:1164
// ```c++
// #define SERIAL_ALLOW_UTF8_IN_VISIBLESTRING_ON_WRITING 0
// ```
// NCBI reference (598d8ae6): c++/src/serial/objistr.cpp:1855-1898
// ```c++
// char ReplaceVisibleChar(char c, EFixNonPrint fix_method,
//     const CObjectStack* io, const CTempString& str, char subst)
// {
//     ...
//     if ( fix_method != eFNP_Replace ) {
//         string message;
//         if (io != NULL) {
//             message += io->GetStackTrace() + "\n";
//         }
//         message += "Bad char [0x" +
//                    NStr::NumericToString((unsigned char)c,0,16)+
//                    "] in string";
//         if (io != NULL) {
//             message += " at " + io->GetPosition();
//         }
//         if (!str.empty()) {
//             message += "\n" + str;
//         }
//         switch (fix_method) {
//         case eFNP_ReplaceAndWarn:
//             ...
//             ERR_POST_X(7, message);
//             break;
//     ...
//     return subst;
// }
// ```

/// `GoodVisibleChar`: a byte that a VisibleString keeps.
// NCBI reference (598d8ae6): c++/include/serial/impl/objistr.inl:455-459
// ```c++
// inline
// bool GoodVisibleChar(char c)
// {
//     return c >= ' ' && c <= '~';
// }
// ```
#[inline]
pub(crate) fn good_visible_char(c: u8) -> bool {
    (b' '..=b'~').contains(&c)
}

/// `m_NonPrintSubst`, the byte stored for a bad one (objistr.cpp:463).
pub(crate) const SUBST: u8 = b'#';

/// `GetPosition()` of the reader at a bad byte: the line of a text stream
/// (`CObjectIStreamAsn::Location`, objistrasn.cpp:111-115) or the bytes a binary stream has read
/// (`CObjectIStream::Location`, objistr.cpp:815-819).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum Position {
    Line(usize),
    Byte(usize),
}

/// How NCBI's message shows the string of a bad byte (`str` of `ReplaceVisibleChar`).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum Shown {
    /// The bytes as read or given: a text reader's chunk (`CTempString(data, count)`,
    /// objistrasn.cpp:893), a writer's `original` string (objostrasnb.cpp:983).
    AsIs,
    /// A binary reader's string as it is being fixed in place (objistrasnb.cpp:703-749): the bad
    /// bytes before this one already [`SUBST`].
    BeingFixed,
}

/// The bytes of one string (a text chunk, a binary string, a string written) that NCBI's stream
/// stores or writes as [`SUBST`], with what each message shows: the stream's trace at the string
/// (`GetStackTraceASN`), the stream's position at the byte and the string, which is kept once
/// (not per byte).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BadChars {
    trace: String,
    bytes: Vec<u8>,
    shown: Shown,
    /// Each bad byte: its place in `bytes` and the stream's position at it.
    bad: Vec<(usize, Position)>,
    /// For a writer that models NCBI's output buffer (`asn_text`), the bytes of the document
    /// that NCBI's stream has handed to its output when it posts each message (one per entry of
    /// `bad`); empty otherwise (the messages then come before the document).
    after: Vec<usize>,
}

impl BadChars {
    /// The bad bytes of `bytes` (`None` when it has none), all at `at`, for the stream whose
    /// trace is `trace`.
    pub(crate) fn find(trace: String, bytes: &[u8], shown: Shown, at: Position) -> Option<Self> {
        let bad: Vec<(usize, Position)> = bytes
            .iter()
            .enumerate()
            .filter(|&(_, &c)| !good_visible_char(c))
            .map(|(i, _)| (i, at))
            .collect();
        (!bad.is_empty()).then(|| BadChars {
            trace,
            bytes: bytes.to_vec(),
            shown,
            bad,
            after: Vec::new(),
        })
    }

    /// A string written by a stream whose trace is `trace`, before its bad bytes are noted
    /// (`note`), each at the position the writer is at.
    pub(crate) fn written(trace: String, original: &[u8]) -> Self {
        BadChars {
            trace,
            bytes: original.to_vec(),
            shown: Shown::AsIs,
            bad: Vec::new(),
            after: Vec::new(),
        }
    }

    /// The bad byte at `index` of the string, written at `at`.
    pub(crate) fn note(&mut self, index: usize, at: Position) {
        self.bad.push((index, at));
    }

    /// `note`, by a writer whose model of NCBI's buffer says that NCBI's stream has handed the
    /// first `after` bytes of the document to its output when it posts the message.
    pub(crate) fn note_after(&mut self, index: usize, at: Position, after: usize) {
        self.bad.push((index, at));
        self.after.push(after);
    }

    /// The number of bytes noted (one message each).
    pub(crate) fn len(&self) -> usize {
        self.bad.len()
    }

    /// The `k`-th bad byte and the position of the stream at it.
    pub(crate) fn byte(&self, k: usize) -> (u8, Position) {
        let (index, at) = self.bad[k];
        (self.bytes[index], at)
    }

    /// The trace of the stream at the string.
    #[cfg(test)]
    pub(crate) fn trace(&self) -> &str {
        &self.trace
    }

    /// The string as the message of the `k`-th bad byte shows it.
    pub(crate) fn shown(&self, k: usize) -> Vec<u8> {
        match self.shown {
            Shown::AsIs => self.bytes.clone(),
            Shown::BeingFixed => {
                let index = self.bad[k].0;
                let mut s = self.bytes.clone();
                for c in &mut s[..index] {
                    if !good_visible_char(*c) {
                        *c = SUBST;
                    }
                }
                s
            }
        }
    }

    /// The text that NCBI posts for the `k`-th bad byte (`ERR_POST_X(7, message)`, without the
    /// prefix of the diagnostic): the trace, "Bad char [0x..] in string at line N" or "at byte
    /// N", and the string as it shows it.
    pub(crate) fn message(&self, k: usize) -> Vec<u8> {
        let (byte, at) = self.byte(k);
        let (kind, n) = match at {
            Position::Line(n) => ("line", n),
            Position::Byte(n) => ("byte", n),
        };
        let mut message = format!(
            "{}\nBad char [0x{byte:X}] in string at {kind} {n}",
            self.trace
        )
        .into_bytes();
        if !self.bytes.is_empty() {
            message.push(b'\n');
            message.extend_from_slice(&self.shown(k));
        }
        message
    }

    /// The lines that NCBI's application `program` writes to its standard error for these bytes:
    /// `Error: (803.7) [program] <message>` and a newline, one per byte, appended to `out`; an
    /// application without a name in its diagnostics (`None`: the masking programs, such as
    /// convert2blastmask 2.17.0) writes `Error: (803.7) <message>`.
    pub(crate) fn post(&self, program: Option<&str>, out: &mut Vec<u8>) {
        for k in 0..self.len() {
            self.post_one(k, program, out);
        }
    }

    /// The post of the `k`-th bad byte (`post`).
    fn post_one(&self, k: usize, program: Option<&str>, out: &mut Vec<u8>) {
        out.extend_from_slice(b"Error: (803.7) ");
        if let Some(program) = program {
            out.extend_from_slice(format!("[{program}] ").as_bytes());
        }
        out.extend_from_slice(&self.message(k));
        out.push(b'\n');
    }

    /// `post` of every string of `bad`, in order.
    pub(crate) fn post_all(bad: &[BadChars], program: Option<&str>, out: &mut Vec<u8>) {
        for bad_chars in bad {
            bad_chars.post(program, out);
        }
    }
}

/// Each post of `BadChars::post_all` with the bytes of its document that come before it on a
/// stream that both go to (`BadChars::note_after`; 0 where the writer has no model).
pub(crate) fn posts_after(bad: &[BadChars], program: Option<&str>) -> Vec<(usize, Vec<u8>)> {
    let mut posts = Vec::new();
    for bad_chars in bad {
        for k in 0..bad_chars.len() {
            let mut text = Vec::new();
            bad_chars.post_one(k, program, &mut text);
            posts.push((bad_chars.after.get(k).copied().unwrap_or(0), text));
        }
    }
    posts
}

/// Each post of `BadChars::post_all` on its own (one message of NCBI's diagnostic handler each).
pub(crate) fn each_post(bad: &[BadChars], program: Option<&str>) -> Vec<Vec<u8>> {
    let mut posts = Vec::new();
    for bad_chars in bad {
        for k in 0..bad_chars.len() {
            let mut text = Vec::new();
            bad_chars.post_one(k, program, &mut text);
            posts.push(text);
        }
    }
    posts
}

/// `BadChars::post_all` written to `diagnostics` (nothing for no string).
pub(crate) fn write_posts(
    bad: &[BadChars],
    program: Option<&str>,
    diagnostics: &mut dyn std::io::Write,
) -> std::io::Result<()> {
    if bad.is_empty() {
        return Ok(());
    }
    let mut text = Vec::new();
    BadChars::post_all(bad, program, &mut text);
    diagnostics.write_all(&text)
}

/// `CObjectIStreamAsnBinary::FixVisibleChars(string&, eFNP_ReplaceAndWarn)`: each bad byte of
/// `s` becomes [`SUBST`] in turn; the bytes are noted at `at` for the stream whose trace is
/// `trace`, each message showing the string as it is then (the bytes before it already replaced).
// NCBI reference (598d8ae6): c++/src/serial/objistrasnb.cpp:703-749
// ```c++
// bool CObjectIStreamAsnBinary::FixVisibleChars(string& str, EFixNonPrint fix_method)
// {
//     ...
//     } else {
//         for(string::iterator i = str.begin(); i < to; ++i) {
//             if ( BadVisibleChar(*i) ) {
//                 ...
//                 {
//                     char c = ReplaceVisibleChar(*i, fix_method, this, str, subst);
//                     if (c != '\0') {
//                         *i = subst;
// ```
pub(crate) fn fix_visible_chars(s: &mut [u8], trace: String, at: Position) -> Option<BadChars> {
    let bad = BadChars::find(trace, s, Shown::BeingFixed, at)?;
    for c in s.iter_mut() {
        if !good_visible_char(*c) {
            *c = SUBST;
        }
    }
    Some(bad)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// makeblastdb 2.17.0 on a binary mask file whose `algo-options` holds 0x09 and 0xE9
    /// (`M/m4a/tmp/fb_probe/b1.asnb`): two messages, the second shows the first byte replaced.
    #[test]
    fn a_binary_string_is_fixed_byte_by_byte_as_ncbi_shows_it() {
        let mut s = b"\ti\xe9dow".to_vec();
        let trace = "Blast-db-mask-info.algo-options".to_string();
        let bad = fix_visible_chars(&mut s, trace, Position::Byte(49)).unwrap();
        assert_eq!(s, b"#i#dow");
        let mut posted = Vec::new();
        bad.post(Some("makeblastdb"), &mut posted);
        assert_eq!(
            posted,
            b"Error: (803.7) [makeblastdb] Blast-db-mask-info.algo-options\nBad char [0x9] in string at byte 49\n\ti\xe9dow\n\
              Error: (803.7) [makeblastdb] Blast-db-mask-info.algo-options\nBad char [0xE9] in string at byte 49\n#i\xe9dow\n"
                .to_vec()
        );
        let mut ok = b"window".to_vec();
        assert!(fix_visible_chars(&mut ok, String::new(), Position::Byte(1)).is_none());
    }

    /// A string written: each byte at its own position, the original string shown.
    #[test]
    fn a_written_string_shows_the_original() {
        let mut bad = BadChars::written("T.E.title".to_string(), b"t\xe9t\x01e");
        bad.note(1, Position::Byte(21));
        bad.note(3, Position::Byte(23));
        assert_eq!(bad.len(), 2);
        assert_eq!(
            bad.message(0),
            b"T.E.title\nBad char [0xE9] in string at byte 21\nt\xe9t\x01e"
        );
        assert_eq!(
            bad.message(1),
            b"T.E.title\nBad char [0x1] in string at byte 23\nt\xe9t\x01e"
        );
    }
}
