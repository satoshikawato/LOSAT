//! NCBI's repair of a VisibleString as its object streams read it: the method
//! `eFNP_ReplaceAndWarn` (the default of `SERIAL_WRONG_CHARS_READ`) stores `#` for each byte
//! outside `' '..='~'` and posts an error for it with the reader's stack trace, its position and
//! the bytes being read (`ReplaceVisibleChar`). The readers of ASN.1 text
//! (`asn_text_reader::Reader::string`, which also cuts the string into NCBI's chunks) and of
//! ASN.1 binary (`fix_visible_chars`, the whole string) note each such byte as a [`BadChar`];
//! the reader that knows its frames gives the trace ([`BadChar::message`]).
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

/// One byte that NCBI's reader replaced in a VisibleString, with what its message shows.
#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct BadChar {
    pub(crate) byte: u8,
    pub(crate) at: Position,
    /// `str`: the chunk of a text string being appended, as in the text; the whole binary string
    /// as it is being fixed (the bytes before this one already `#`); empty for a text string that
    /// is skipped (`SkipString` passes `kEmptyStr`).
    pub(crate) chunk: Vec<u8>,
    /// `GetStackTraceASN` of the text reader at the byte (`Reader::frames`: the type, then a
    /// name for each member, variant and `E` for each element); empty when the reader keeps no
    /// frames (the binary reader and the mask reader give the trace themselves).
    pub(crate) trace: String,
}

impl BadChar {
    /// `message` with the trace the text reader noted.
    pub(crate) fn traced_message(&self) -> Vec<u8> {
        self.message(&self.trace)
    }

    /// The text that NCBI posts (`ERR_POST_X(7, message)`, without the prefix of the
    /// diagnostic) for the reader whose stack trace is `trace` (`GetStackTraceASN`).
    pub(crate) fn message(&self, trace: &str) -> Vec<u8> {
        let (kind, n) = match self.at {
            Position::Line(n) => ("line", n),
            Position::Byte(n) => ("byte", n),
        };
        let mut message = format!(
            "{trace}\nBad char [0x{:X}] in string at {kind} {n}",
            self.byte
        )
        .into_bytes();
        if !self.chunk.is_empty() {
            message.push(b'\n');
            message.extend_from_slice(&self.chunk);
        }
        message
    }
}

/// `CObjectIStreamAsnBinary::FixVisibleChars(string&, eFNP_ReplaceAndWarn)`: each bad byte of
/// `s` becomes [`SUBST`] in turn, and is noted in `notes` at `at` with the string as it is then
/// (the bytes before it already replaced).
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
pub(crate) fn fix_visible_chars(s: &mut [u8], at: Position, notes: &mut Vec<BadChar>) {
    for i in 0..s.len() {
        let byte = s[i];
        if !good_visible_char(byte) {
            notes.push(BadChar {
                byte,
                at,
                chunk: s.to_vec(),
                trace: String::new(),
            });
            s[i] = SUBST;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// makeblastdb 2.17.0 on a binary mask file whose `algo-options` holds 0x09 and 0xE9
    /// (`M/m4a/tmp/fb_probe/b1.asnb`): two messages, the second shows the first byte replaced.
    #[test]
    fn a_binary_string_is_fixed_byte_by_byte_as_ncbi_shows_it() {
        let mut s = b"\ti\xe9dow".to_vec();
        let mut notes = Vec::new();
        fix_visible_chars(&mut s, Position::Byte(49), &mut notes);
        assert_eq!(s, b"#i#dow");
        let trace = "Blast-db-mask-info.algo-options";
        assert_eq!(
            notes
                .iter()
                .map(|note| note.message(trace))
                .collect::<Vec<_>>(),
            [
                b"Blast-db-mask-info.algo-options\nBad char [0x9] in string at byte 49\n\ti\xe9dow"
                    .to_vec(),
                b"Blast-db-mask-info.algo-options\nBad char [0xE9] in string at byte 49\n#i\xe9dow"
                    .to_vec(),
            ]
        );
        let skipped = BadChar {
            byte: 0x7f,
            at: Position::Line(3),
            chunk: Vec::new(),
            trace: String::new(),
        };
        assert_eq!(
            skipped.message("T.m"),
            b"T.m\nBad char [0x7F] in string at line 3"
        );
    }
}
