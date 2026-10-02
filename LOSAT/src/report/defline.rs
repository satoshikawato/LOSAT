//! The subject title of the pairwise report (outfmt 0) as NCBI makes it for a nucleotide
//! subject read from FASTA: `CDeflineGenerator::GenerateDefline` over the defline, which
//! is the title of a subject read without `-parse_deflines` (its ID included).
//!
//! NCBI reference: c++/src/objmgr/util/create_defline.cpp:3952-3960
//! ```c
//!     if (! m_Reconstruct) {
//!         // x_SetFlags set m_MainTitle from a suitable descriptor, if any;
//!         // now strip trailing periods, commas, semicolons, and spaces.
//!         size_t pos = m_MainTitle.find_last_not_of (".,;~ ");
//!         if (pos != NPOS) {
//!             m_MainTitle.erase (pos + 1);
//!         }
//! ```
//! NCBI reference: c++/src/objmgr/util/create_defline.cpp:4050-4095
//! ```c
//!     if (! (flags & fLeavePrefixSuffix)) {
//!         // remove TPA or TSA prefix, will rely on other data in record to set
//!         for (size_t i = 0; i < sizeof (s_tpaPrefixList) / sizeof (const char*); i++) {
//!             string str = s_tpaPrefixList [i];
//!             if (NStr::StartsWith (m_MainTitle, str, NStr::eNocase)) {
//!                 m_MainTitle.erase (0, str.length());
//!                 // strip leading spaces remaining after removal of old MAG before TPA or TSA prefixes
//!                 m_MainTitle.erase (0, m_MainTitle.find_first_not_of (' '));
//!             }
//!         }
//!     }
//!
//!     // strip leading spaces remaining after removal of old TPA or TSA prefixes
//!     m_MainTitle.erase (0, m_MainTitle.find_first_not_of (' '));
//!
//!     CStringUTF8 decoded = NStr::HtmlDecode (m_MainTitle);
//!
//!     // strip trailing commas, semicolons, and spaces (period may be an sp.
//!     // species)
//!     size_t pos = decoded.find_last_not_of (",;~ ");
//!     if (pos != NPOS) {
//!         decoded.erase (pos + 1);
//!     }
//!     ...
//!     // produce final result
//!     string penult = mag + prefix + decoded + suffix;
//!
//!     x_CleanAndCompress (final, penult, m_IsAA);
//! ```
//! A local FASTA subject has no MolInfo, source or metagenome data, so `mag`, `prefix` and
//! `suffix` are empty; its title is not empty (LOSAT rejects empty deflines), so it is not
//! capitalized. `NStr::HtmlDecode` changes only character references, which the callers
//! reject (`has_html_character_reference`, a superset of them). The alignment heading
//! removes the prefixes; the description table keeps them (`fLeavePrefixSuffix`,
//! showdefline.cpp:498). For some titles of commas, semicolons, tildes and spaces, NCBI's
//! `x_CleanAndCompress` lets its count of the remaining letters wrap and reads past the
//! end of the string (and crashes); `clean_and_compress` stops at the end of the string
//! there (approved exception 2 of PD-LOSAT-NCBI-DEFECTS).

/// NCBI reference: c++/src/objmgr/util/create_defline.cpp:3431-3446
/// ```c
/// static const char* s_tpaPrefixList [] = {
///   "MAG ",
///   "MAG:",
///   "MULTISPECIES:",
///   "TLS:",
///   "TPA:",
///   "TPA_exp:",
///   "TPA_inf:",
///   "TPA_reasm:",
///   "TPA_asm:",
///   "TSA:",
///   "UNVERIFIED_ORG:",
///   "UNVERIFIED_ASMBLY:",
///   "UNVERIFIED_CONTAM:",
///   "UNVERIFIED:"
/// };
/// ```
const TPA_PREFIXES: [&str; 14] = [
    "MAG ",
    "MAG:",
    "MULTISPECIES:",
    "TLS:",
    "TPA:",
    "TPA_exp:",
    "TPA_inf:",
    "TPA_reasm:",
    "TPA_asm:",
    "TSA:",
    "UNVERIFIED_ORG:",
    "UNVERIFIED_ASMBLY:",
    "UNVERIFIED_CONTAM:",
    "UNVERIFIED:",
];

/// The title of a nucleotide subject with the defline `defline` (the text after `>`),
/// for the alignment heading (`leave_prefix` false) or the description table (true).
pub fn ncbi_nucleotide_title(defline: &str, leave_prefix: bool) -> String {
    let mut title = defline.as_bytes().to_vec();
    trim_end_of(&mut title, b".,;~ ");
    if !leave_prefix {
        for prefix in TPA_PREFIXES {
            if title.len() >= prefix.len()
                && title[..prefix.len()].eq_ignore_ascii_case(prefix.as_bytes())
            {
                title.drain(..prefix.len());
                trim_leading_spaces(&mut title);
            }
        }
    }
    trim_leading_spaces(&mut title);
    trim_end_of(&mut title, b",;~ ");
    let cleaned = clean_and_compress(&title, false);
    String::from_utf8(cleaned).expect("ASCII deflines")
}

/// Whether the text may hold a character reference that `NStr::HtmlDecode` decodes: an `&`
/// followed by a letter, `#` and a digit, or `#x` and a hexadecimal digit, with a `;`
/// within the next 16 characters before another `&` or `#` (ncbistr.cpp:4545-4570). This
/// is a superset: it does not look the names up in NCBI's table (`&foo;` stays) and does
/// not follow NCBI's trims (a final `;` is trimmed before the decoding).
pub fn has_html_character_reference(text: &str) -> bool {
    let bytes = text.as_bytes();
    bytes.iter().enumerate().any(|(index, &byte)| {
        if byte != b'&' {
            return false;
        }
        let rest = &bytes[index + 1..];
        let start = match rest {
            [first, ..] if first.is_ascii_alphabetic() => 0,
            [b'#', digit, ..] if digit.is_ascii_digit() => 1,
            [b'#', b'x' | b'X', hex, ..] if hex.is_ascii_hexdigit() => 2,
            _ => return false,
        };
        rest[start..]
            .iter()
            .take(16)
            .take_while(|&&c| c != b'&' && c != b'#')
            .any(|&c| c == b';')
    })
}

/// `find_last_not_of(chars)` then `erase(pos + 1)`; a text made only of `chars` stays.
fn trim_end_of(text: &mut Vec<u8>, chars: &[u8]) {
    if let Some(pos) = text.iter().rposition(|byte| !chars.contains(byte)) {
        text.truncate(pos + 1);
    }
}

/// `erase(0, find_first_not_of(' '))`: a text of spaces only becomes empty.
fn trim_leading_spaces(text: &mut Vec<u8>) {
    let first = text
        .iter()
        .position(|&byte| byte != b' ')
        .unwrap_or(text.len());
    text.drain(..first);
}

/// NCBI's `x_CleanAndCompress`, byte for byte (the titles are ASCII), except that it stops
/// at the end of the string where NCBI's count of the remaining letters wraps (see below).
///
/// NCBI reference: c++/src/objmgr/util/create_defline.cpp:219-312
/// ```c
///     char curr = *in++; // initialize with first character
///     left--;
///
///     char next = 0;
///     Uint2 two_chars = curr; // this is two bytes storage where we see current and previous symbols
///
///     while (left > 0) {
///         next = *in++;
///
///         two_chars = Uint2((two_chars << 8) | next);
///
///         switch (two_chars)
///         {
///         case twocommas: // replace double commas with comma+space
///             *out++ = curr;
///             next = ' ';
///             break;
///         case twospaces: // skip multispaces (only print last one)
///             break;
///         case bracket_space: // skip space after bracket
///             next = curr;
///             two_chars = curr;
///             break;
///         case space_bracket: // skip space before bracket
///             break;
///         case space_comma:
///         case space_semicolon: // swap characters
///             *out++ = next;
///             next = curr;
///             two_chars = curr;
///             break;
///         case comma_space:
///             *out++ = curr;
///             *out++ = ' ';
///             while (next == ' ' || next == ',') {
///                 next = *in;
///                 in++;
///                 left--;
///             }
///             two_chars = next;
///             break;
///         case semicolon_space:
///     ...
///         default:
///             *out++ = curr;
///             break;
///         }
///
///         curr = next;
///         left--;
///     }
///
///     if (curr > 0 && curr != ' ') {
///         *out++ = curr;
///     }
/// ```
fn clean_and_compress(input: &[u8], is_protein: bool) -> Vec<u8> {
    let mut start = 0;
    let mut end = input.len();
    while start < end && input[start] == b' ' {
        start += 1;
    }
    while end > start && input[end - 1] == b' ' {
        end -= 1;
    }
    let text = &input[start..end];
    let mut out = Vec::with_capacity(text.len());
    if text.is_empty() {
        return out;
    }
    // Past the end reads the terminating NUL of the C++ string.
    let at = |index: usize| text.get(index).copied().unwrap_or(0);
    let mut index = 0;
    let mut left = text.len();
    let mut curr = at(index);
    index += 1;
    left -= 1;
    let mut two_chars = u16::from(curr);
    while left > 0 {
        let mut next = at(index);
        index += 1;
        two_chars = (two_chars << 8) | u16::from(next);
        match two_chars.to_be_bytes() {
            [b',', b','] => {
                out.push(curr);
                next = b' ';
            }
            [b' ', b' '] | [b' ', b')'] => {}
            [b'(', b' '] => {
                next = curr;
                two_chars = u16::from(curr);
            }
            [b' ', b','] | [b' ', b';'] => {
                out.push(next);
                next = curr;
                two_chars = u16::from(curr);
            }
            [separator @ (b',' | b';'), b' '] => {
                out.push(curr);
                out.push(b' ');
                while next == b' ' || next == separator {
                    next = at(index);
                    index += 1;
                    // Where NCBI's `left` (a size_t) wraps, for a run of spaces and
                    // separators that reaches the end of the string, its loop runs past the
                    // string (it reads the terminating NUL and on, and crashes). LOSAT stops
                    // at the end of the string (approved exception 2 of
                    // PD-LOSAT-NCBI-DEFECTS): `next` is the NUL and nothing more is written.
                    left = left.saturating_sub(1);
                }
                two_chars = u16::from(next);
            }
            _ => out.push(curr),
        }
        curr = next;
        left = left.saturating_sub(1);
    }
    if curr > 0 && curr != b' ' {
        out.push(curr);
    }
    if is_protein {
        let replaced = String::from_utf8_lossy(&out)
            .replace(". [", " [")
            .replace(", [", " [");
        return replaced.into_bytes();
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn titles_follow_ncbi_defline_generator() {
        // NCBI BLAST+ 2.17.0 outfmt 0 (alignment heading, description table).
        for (defline, heading, description) in [
            ("s1 abc.", "s1 abc", "s1 abc"),
            ("s1 abc ,def", "s1 abc, def", "s1 abc, def"),
            ("s1 a  b", "s1 a b", "s1 a b"),
            ("s1 (a )", "s1 (a)", "s1 (a)"),
            ("s1,,x", "s1, x", "s1, x"),
            ("TPA: s1 x", "s1 x", "TPA: s1 x"),
            ("s1 TPA: x", "s1 TPA: x", "s1 TPA: x"),
            ("s1 MAG x", "s1 MAG x", "s1 MAG x"),
            ("s1 x ; y", "s1 x; y", "s1 x; y"),
            ("s1 x,  y", "s1 x, y", "s1 x, y"),
            ("s1 x~", "s1 x", "s1 x"),
            ("s1 E. coli sp.", "s1 E. coli sp", "s1 E. coli sp"),
            ("s1 a, ,b", "s1 a, b", "s1 a, b"),
            ("s1 x;;y", "s1 x;;y", "s1 x;;y"),
            ("s1  lead", "s1 lead", "s1 lead"),
            ("s1 a ( b", "s1 a (b", "s1 a (b"),
            ("s1..x", "s1..x", "s1..x"),
            ("s1.x.", "s1.x", "s1.x"),
        ] {
            assert_eq!(
                ncbi_nucleotide_title(defline, false),
                heading,
                "{defline:?}"
            );
            assert_eq!(
                ncbi_nucleotide_title(defline, true),
                description,
                "{defline:?}"
            );
        }
        // NCBI BLAST+ 2.17.0 crashes on the outfmt 0 titles of the first group
        // (x_CleanAndCompress reads past the string); LOSAT stops at the end of the string
        // (approved exception 2 of PD-LOSAT-NCBI-DEFECTS). The second group is NCBI's.
        for (defline, title) in [
            (", ,", ", "),
            ("; ;", "; "),
            ("~, ,", "~, "),
            (",, ,", ",  "),
            (", ,,", ", "),
            (";  ;", "; "),
            (", , ,", ", "),
            (",~, ,", ",~, "),
            (", ;", ", ;"),
            (",,", ","),
            ("a, ,b", "a, b"),
            (", ,a", ", a"),
        ] {
            assert_eq!(ncbi_nucleotide_title(defline, false), title, "{defline:?}");
        }
    }

    #[test]
    fn html_character_references_are_found() {
        for text in ["a&amp;b", "a&#38;b", "a&#x26;b", "x &lt; y"] {
            assert!(has_html_character_reference(text), "{text}");
        }
        for text in ["a & b", "a&b", "a&;", "a&#;", "a&amp b", "a&x#y;"] {
            assert!(!has_html_character_reference(text), "{text}");
        }
    }
}
