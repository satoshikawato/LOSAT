//! HTML entity decoding, encoding detection, string cleaning, and byte helpers for
//! `CDeflineGenerator`.
//!
//! NCBI reference (598d8ae6): c++/src/corelib/ncbistr.cpp:4221-4618,6711-7172
//! NCBI reference (598d8ae6): c++/src/objmgr/util/create_defline.cpp:219-312

use crate::ncbi_str::replace_in_place;

use super::UnknownEncoding;

/// The entities of `NStr::HtmlDecode`, in NCBI's order.
///
/// NCBI reference (598d8ae6): c++/src/corelib/ncbistr.cpp:4221-4508
#[rustfmt::skip]
const NCBI_HTML_ENTITIES: [(u32, &str); 280] = [
    (9, "Tab"), (10, "NewLine"), (33, "excl"), (34, "quot"), (35, "num"), (36, "dollar"),
    (37, "percnt"), (38, "amp"), (39, "apos"), (40, "lpar"), (41, "rpar"), (42, "ast"),
    (43, "plus"), (44, "comma"), (46, "period"), (47, "sol"), (58, "colon"), (59, "semi"),
    (60, "lt"), (61, "equals"), (62, "gt"), (63, "quest"), (64, "commat"), (91, "lsqb"),
    (92, "bsol"), (93, "rsqb"), (94, "Hat"), (95, "lowbar"), (96, "grave"), (123, "lcub"),
    (124, "verbar"), (125, "rcub"), (160, "nbsp"), (161, "iexcl"), (162, "cent"),
    (163, "pound"), (164, "curren"), (165, "yen"), (166, "brvbar"), (167, "sect"),
    (168, "uml"), (169, "copy"), (170, "ordf"), (171, "laquo"), (172, "not"), (173, "shy"),
    (174, "reg"), (175, "macr"), (176, "deg"), (177, "plusmn"), (178, "sup2"), (179, "sup3"),
    (180, "acute"), (181, "micro"), (182, "para"), (183, "middot"), (184, "cedil"),
    (185, "sup1"), (186, "ordm"), (187, "raquo"), (188, "frac14"), (189, "frac12"),
    (190, "frac34"), (191, "iquest"), (192, "Agrave"), (193, "Aacute"), (194, "Acirc"),
    (195, "Atilde"), (196, "Auml"), (197, "Aring"), (198, "AElig"), (199, "Ccedil"),
    (200, "Egrave"), (201, "Eacute"), (202, "Ecirc"), (203, "Euml"), (204, "Igrave"),
    (205, "Iacute"), (206, "Icirc"), (207, "Iuml"), (208, "ETH"), (209, "Ntilde"),
    (210, "Ograve"), (211, "Oacute"), (212, "Ocirc"), (213, "Otilde"), (214, "Ouml"),
    (215, "times"), (216, "Oslash"), (217, "Ugrave"), (218, "Uacute"), (219, "Ucirc"),
    (220, "Uuml"), (221, "Yacute"), (222, "THORN"), (223, "szlig"), (224, "agrave"),
    (225, "aacute"), (226, "acirc"), (227, "atilde"), (228, "auml"), (229, "aring"),
    (230, "aelig"), (231, "ccedil"), (232, "egrave"), (233, "eacute"), (234, "ecirc"),
    (235, "euml"), (236, "igrave"), (237, "iacute"), (238, "icirc"), (239, "iuml"),
    (240, "eth"), (241, "ntilde"), (242, "ograve"), (243, "oacute"), (244, "ocirc"),
    (245, "otilde"), (246, "ouml"), (247, "divide"), (248, "oslash"), (249, "ugrave"),
    (250, "uacute"), (251, "ucirc"), (252, "uuml"), (253, "yacute"), (254, "thorn"),
    (255, "yuml"), (338, "OElig"), (339, "oelig"), (352, "Scaron"), (353, "scaron"),
    (376, "Yuml"), (402, "fnof"), (710, "circ"), (732, "tilde"), (913, "Alpha"),
    (914, "Beta"), (915, "Gamma"), (916, "Delta"), (917, "Epsilon"), (918, "Zeta"),
    (919, "Eta"), (920, "Theta"), (921, "Iota"), (922, "Kappa"), (923, "Lambda"), (924, "Mu"),
    (925, "Nu"), (926, "Xi"), (927, "Omicron"), (928, "Pi"), (929, "Rho"), (931, "Sigma"),
    (932, "Tau"), (933, "Upsilon"), (934, "Phi"), (935, "Chi"), (936, "Psi"), (937, "Omega"),
    (945, "alpha"), (946, "beta"), (947, "gamma"), (948, "delta"), (949, "epsilon"),
    (950, "zeta"), (951, "eta"), (952, "theta"), (953, "iota"), (954, "kappa"),
    (955, "lambda"), (956, "mu"), (957, "nu"), (958, "xi"), (959, "omicron"), (960, "pi"),
    (961, "rho"), (962, "sigmaf"), (963, "sigma"), (964, "tau"), (965, "upsilon"),
    (966, "phi"), (967, "chi"), (968, "psi"), (969, "omega"), (977, "thetasym"),
    (978, "upsih"), (982, "piv"), (8194, "ensp"), (8195, "emsp"), (8201, "thinsp"),
    (8204, "zwnj"), (8205, "zwj"), (8206, "lrm"), (8207, "rlm"), (8211, "ndash"),
    (8212, "mdash"), (8216, "lsquo"), (8217, "rsquo"), (8218, "sbquo"), (8220, "ldquo"),
    (8221, "rdquo"), (8222, "bdquo"), (8224, "dagger"), (8225, "Dagger"), (8226, "bull"),
    (8230, "hellip"), (8240, "permil"), (8242, "prime"), (8243, "Prime"), (8249, "lsaquo"),
    (8250, "rsaquo"), (8254, "oline"), (8260, "frasl"), (8364, "euro"), (8472, "weierp"),
    (8465, "image"), (8476, "real"), (8482, "trade"), (8501, "alefsym"), (8592, "larr"),
    (8593, "uarr"), (8594, "rarr"), (8595, "darr"), (8596, "harr"), (8629, "crarr"),
    (8656, "lArr"), (8657, "uArr"), (8658, "rArr"), (8659, "dArr"), (8660, "hArr"),
    (8704, "forall"), (8706, "part"), (8707, "exist"), (8709, "empty"), (8711, "nabla"),
    (8712, "isin"), (8713, "notin"), (8715, "ni"), (8719, "prod"), (8721, "sum"),
    (8722, "minus"), (8727, "lowast"), (8730, "radic"), (8733, "prop"), (8734, "infin"),
    (8736, "ang"), (8743, "and"), (8744, "or"), (8745, "cap"), (8746, "cup"), (8747, "int"),
    (8756, "there4"), (8764, "sim"), (8773, "cong"), (8776, "asymp"), (8800, "ne"),
    (8801, "equiv"), (8804, "le"), (8805, "ge"), (8834, "sub"), (8835, "sup"), (8836, "nsub"),
    (8838, "sube"), (8839, "supe"), (8853, "oplus"), (8855, "otimes"), (8869, "perp"),
    (8901, "sdot"), (8968, "lceil"), (8969, "rceil"), (8970, "lfloor"), (8971, "rfloor"),
    (9001, "lang"), (9002, "rang"), (9674, "loz"), (9824, "spades"), (9827, "clubs"),
    (9829, "hearts"), (9830, "diams"),
];

/// `NStr::HtmlDecode` of `text` with its default `eEncoding_Unknown`: the decoded bytes and
/// whether a character reference was decoded.
///
/// NCBI reference (598d8ae6): c++/src/corelib/ncbistr.cpp:4523-4618
/// NCBI reference (598d8ae6): c++/include/corelib/ncbistr.hpp:5729-5759
pub(super) fn html_decode(text: &[u8]) -> Result<(Vec<u8>, bool), UnknownEncoding> {
    let encoding = guess_encoding(text);
    if encoding == Encoding::Unknown {
        return Err(UnknownEncoding);
    }
    let mut decoded = Vec::with_capacity(text.len());
    let mut char_refs = false;
    let e = text.len();
    let mut i = 0;
    while i < e {
        let ch = text[i];
        i += 1;
        if i != e && ch == b'&' {
            let mut itmp = i;
            let mut end_of_entity = i;
            let mut ent = text[itmp].is_ascii_alphabetic();
            let mut dec = false;
            if !ent && text[itmp] == b'#' {
                itmp += 1;
                dec = itmp != e && text[itmp].is_ascii_digit();
            }
            let mut hex = false;
            if !dec && itmp != e && matches!(text[itmp], b'x' | b'X') {
                itmp += 1;
                hex = itmp != e && text[itmp].is_ascii_hexdigit();
            }
            let start_of_entity = itmp;
            if itmp != e && (ent || dec || hex) {
                let mut len = 0;
                while len < 16 && itmp != e {
                    let c = text[itmp];
                    if c == b'&' || c == b'#' {
                        break;
                    }
                    if c == b';' {
                        end_of_entity = itmp;
                        break;
                    }
                    ent = ent && c.is_ascii_alphanumeric();
                    dec = dec && c.is_ascii_digit();
                    hex = hex && c.is_ascii_hexdigit();
                    len += 1;
                    itmp += 1;
                }
                if end_of_entity != i && (ent || dec || hex) {
                    let mut uch: u32 = 0;
                    let mut parsed = false;
                    if ent {
                        let entity = &text[start_of_entity..end_of_entity];
                        if let Some(&(symbol, _)) = NCBI_HTML_ENTITIES
                            .iter()
                            .find(|(_, name)| name.as_bytes() == entity)
                        {
                            uch = symbol;
                            parsed = true;
                        }
                    } else {
                        parsed = true;
                        for &digit in &text[start_of_entity..end_of_entity] {
                            let mut ud = u32::from(digit);
                            if dec {
                                uch = uch.wrapping_mul(10).wrapping_add(ud.wrapping_sub(0x30));
                            } else if hex {
                                if digit.is_ascii_digit() {
                                    ud -= u32::from(b'0');
                                } else if (b'a'..=b'f').contains(&digit) {
                                    ud = ud - u32::from(b'a') + 10;
                                } else if (b'A'..=b'F').contains(&digit) {
                                    ud = ud - u32::from(b'A') + 10;
                                }
                                uch = uch.wrapping_mul(16).wrapping_add(ud);
                            }
                        }
                    }
                    if parsed {
                        char_refs = true;
                        append_char(&mut decoded, uch);
                        i = end_of_entity + 1;
                        continue;
                    }
                }
            }
        }
        match encoding {
            Encoding::Utf8 | Encoding::Ascii => decoded.push(ch),
            _ => append_as_utf8(&mut decoded, ch, encoding),
        }
    }
    Ok((decoded, char_refs))
}

/// NCBI's `EEncoding` values that `CUtf8::GuessEncoding` returns.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(super) enum Encoding {
    Unknown,
    Utf8,
    Ascii,
    Iso8859_1,
    Windows1252,
    Cesu8,
}

/// `CUtf8::GuessEncoding`.
///
/// NCBI reference (598d8ae6): c++/src/corelib/ncbistr.cpp:6711-6769
pub(super) fn guess_encoding(src: &[u8]) -> Encoding {
    let mut more = 0;
    let (mut cp1252, mut iso1, mut ascii, mut utf8, mut cesu8) = (true, true, true, true, false);
    for (index, &ch) in src.iter().enumerate() {
        let mut skip = false;
        if more != 0 {
            if eval_next(ch) {
                more -= 1;
                if more == 0 {
                    ascii = false;
                }
                skip = true;
            } else {
                more = 0;
                utf8 = false;
            }
        }
        if ch > 0x7F {
            ascii = false;
            if ch < 0xA0 {
                iso1 = false;
                if matches!(ch, 0x81 | 0x8D | 0x8F | 0x90 | 0x9D) {
                    cp1252 = false;
                }
            }
            if !skip && utf8 && !eval_first(ch, &mut more) {
                utf8 = false;
            }
            if utf8 && !cesu8 && ch == 0xED && src.len() - index > 5 {
                let c1 = src[index + 1];
                let c3 = src[index + 3];
                let c4 = src[index + 4];
                if (c1 & 0xA0) == 0xA0 && c3 == 0xED && (c4 & 0xB0) == 0xB0 {
                    cesu8 = true;
                }
            }
        }
    }
    if more != 0 {
        utf8 = false;
    }
    if ascii {
        Encoding::Ascii
    } else if utf8 {
        if cesu8 {
            Encoding::Cesu8
        } else {
            Encoding::Utf8
        }
    } else if cp1252 {
        if iso1 {
            Encoding::Iso8859_1
        } else {
            Encoding::Windows1252
        }
    } else {
        Encoding::Unknown
    }
}

/// NCBI reference (598d8ae6): c++/src/corelib/ncbistr.cpp:7143-7166
fn eval_first(ch: u8, more: &mut usize) -> bool {
    *more = 0;
    if ch & 0x80 != 0 {
        if ch & 0xE0 == 0xC0 {
            if ch & 0xFE == 0xC0 {
                return false;
            }
            *more = 1;
        } else if ch & 0xF0 == 0xE0 {
            *more = 2;
        } else if ch & 0xF8 == 0xF0 {
            if ch > 0xF4 {
                return false;
            }
            *more = 3;
        } else {
            return false;
        }
    }
    true
}

/// NCBI reference (598d8ae6): c++/src/corelib/ncbistr.cpp:7169-7172
fn eval_next(ch: u8) -> bool {
    ch & 0xC0 == 0x80
}

/// `CUtf8::AsUTF8` of a one-byte slice and `encoding`.
///
/// NCBI reference (598d8ae6): c++/src/corelib/ncbistr.cpp:7062-7115
/// NCBI reference (598d8ae6): c++/src/corelib/ncbistr.cpp:6851-6881
fn append_as_utf8(out: &mut Vec<u8>, ch: u8, encoding: Encoding) {
    const CP1252_TABLE: [u32; 32] = [
        0x20AC, 0x003F, 0x201A, 0x0192, 0x201E, 0x2026, 0x2020, 0x2021, 0x02C6, 0x2030, 0x0160,
        0x2039, 0x0152, 0x003F, 0x017D, 0x003F, 0x003F, 0x2018, 0x2019, 0x201C, 0x201D, 0x2022,
        0x2013, 0x2014, 0x02DC, 0x2122, 0x0161, 0x203A, 0x0153, 0x003F, 0x017E, 0x0178,
    ];
    match encoding {
        Encoding::Cesu8 => out.push(ch),
        Encoding::Windows1252 if (0x80..0xA0).contains(&ch) => {
            append_char(out, CP1252_TABLE[usize::from(ch - 0x80)]);
        }
        _ => append_char(out, u32::from(ch)),
    }
}

/// NCBI reference (598d8ae6): c++/src/corelib/ncbistr.cpp:7040-7060
fn append_char(out: &mut Vec<u8>, ch: u32) {
    if ch < 0x80 {
        out.push(ch as u8);
    } else if ch < 0x800 {
        out.push(((ch >> 6) | 0xC0) as u8);
        out.push(((ch & 0x3F) | 0x80) as u8);
    } else if ch < 0x10000 {
        out.push(((ch >> 12) | 0xE0) as u8);
        out.push((((ch >> 6) & 0x3F) | 0x80) as u8);
        out.push(((ch & 0x3F) | 0x80) as u8);
    } else {
        out.push(((ch >> 18) | 0xF0) as u8);
        out.push((((ch >> 12) & 0x3F) | 0x80) as u8);
        out.push((((ch >> 6) & 0x3F) | 0x80) as u8);
        out.push(((ch & 0x3F) | 0x80) as u8);
    }
}

/// Removes trailing bytes in `chars`, preserving a slice composed solely of `chars`.
pub(super) fn trim_end_of(text: &mut Vec<u8>, chars: &[u8]) {
    if let Some(pos) = text.iter().rposition(|byte| !chars.contains(byte)) {
        text.truncate(pos + 1);
    }
}

/// Removes leading ASCII space bytes; a slice of spaces only becomes empty.
pub(super) fn trim_leading_spaces(text: &mut Vec<u8>) {
    let first = text
        .iter()
        .position(|&byte| byte != b' ')
        .unwrap_or(text.len());
    text.drain(..first);
}

/// C-locale `isspace`: space and HT..CR.
#[inline]
pub(super) fn is_c_space(byte: u8) -> bool {
    matches!(byte, b' ' | b'\t' | b'\n' | 0x0b | 0x0c | b'\r')
}

/// C-locale `ispunct`: ASCII printable non-alphanumeric non-space.
#[inline]
pub(super) fn is_c_punct(byte: u8) -> bool {
    byte.is_ascii_punctuation()
}

/// `NStr::TruncateSpacesInPlace` with both-ends truncation under the C locale.
///
/// NCBI reference (598d8ae6): c++/src/corelib/ncbistr.cpp:3197-3234
pub(super) fn truncate_spaces_both(text: &mut Vec<u8>) {
    let start = text.iter().position(|&b| !is_c_space(b));
    let Some(start) = start else {
        text.clear();
        return;
    };
    let end = text.iter().rposition(|&b| !is_c_space(b)).unwrap_or(start);
    if start > 0 {
        text.copy_within(start..=end, 0);
        text.truncate(end - start + 1);
    } else {
        text.truncate(end + 1);
    }
}

/// Byte-slice exact substring test.
#[inline]
pub(super) fn contains_bytes(text: &[u8], pattern: &[u8]) -> bool {
    if pattern.is_empty() {
        return true;
    }
    text.windows(pattern.len()).any(|window| window == pattern)
}

/// Case-insensitive forward search (`NStr::FindNoCase`), returning `None` for empty pattern.
pub(super) fn find_nocase(text: &[u8], pattern: &[u8]) -> Option<usize> {
    if pattern.is_empty() || pattern.len() > text.len() {
        return None;
    }
    text.windows(pattern.len())
        .position(|w| w.eq_ignore_ascii_case(pattern))
}

/// Case-insensitive reverse search (`NStr::Find` with `eNocase` and `eReverseSearch`), returning
/// `None` for empty pattern.
pub(super) fn find_nocase_rev(text: &[u8], pattern: &[u8]) -> Option<usize> {
    if pattern.is_empty() || pattern.len() > text.len() {
        return None;
    }
    text.windows(pattern.len())
        .rposition(|w| w.eq_ignore_ascii_case(pattern))
}

/// Case-insensitive substring test (`NStr::FindNoCase` != NPOS).
#[inline]
pub(super) fn contains_nocase(text: &[u8], pattern: &[u8]) -> bool {
    find_nocase(text, pattern).is_some()
}

/// Case-sensitive forward search starting at `start`.
pub(super) fn find_bytes_from(text: &[u8], pattern: &[u8], start: usize) -> Option<usize> {
    if pattern.is_empty() || start + pattern.len() > text.len() {
        return None;
    }
    text[start..]
        .windows(pattern.len())
        .position(|w| w == pattern)
        .map(|rel| start + rel)
}

/// NCBI's `x_CleanAndCompress`, byte for byte, stopping at the end of the string where NCBI's
/// remaining-letter count wraps (approved exception 2 of PD-LOSAT-NCBI-DEFECTS).
///
/// NCBI reference (598d8ae6): c++/src/objmgr/util/create_defline.cpp:219-312
pub(super) fn clean_and_compress(input: &[u8], is_protein: bool) -> (Vec<u8>, bool) {
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
        return (out, false);
    }
    let at = |index: usize| text.get(index).copied().unwrap_or(0);
    let mut index = 0;
    let mut left = text.len();
    let mut wrapped = false;
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
                    wrapped |= left == 0;
                    left = left.saturating_sub(1);
                }
                two_chars = u16::from(next);
            }
            _ => out.push(curr),
        }
        curr = next;
        wrapped |= left == 0;
        left = left.saturating_sub(1);
    }
    if (1..0x80).contains(&curr) && curr != b' ' {
        out.push(curr);
    }
    if is_protein {
        replace_in_place(&mut out, b". [", b" [");
        replace_in_place(&mut out, b", [", b" [");
    }
    (out, wrapped)
}
