//! The HTML form of the pairwise and organism reports (`-html`, NCBI's `CBlastFormat::m_IsHTML`):
//! NCBI's `CHTMLHelper::HTMLEncode` and its `NStr::Wrap` with `fWrap_HTMLPre`, which the report
//! writers call only when `-html` is set, once per query or per subject (the plain path does not
//! reach this module).

use std::io::{self, Write};

// NCBI reference (2.17.0): c++/src/html/htmlhelper.cpp:146-149
// ```c++
// string CHTMLHelper::HTMLEncode(const string& str, THTMLEncodeFlags flags)
// {
//     return s_HTMLEncode(str, "\"&<>", flags);
// }
// ```
// NCBI reference (2.17.0): c++/src/html/htmlhelper.cpp:58-130 (`s_HTMLEncode`): `"` -> `&quot;`,
// `&` -> `&` + `amp;` (the default `fEncodeAll` sets no `fSkipEntities`), `<` -> `&lt;`,
// `>` -> `&gt;`; every other byte (`'` too) is copied.
/// `CHTMLHelper::HTMLEncode(text)` with its default flags (`fEncodeAll`), appended to `out`.
pub(crate) fn html_encode_into(text: &[u8], out: &mut Vec<u8>) {
    let mut last = 0;
    for (i, &c) in text.iter().enumerate() {
        let entity: &[u8] = match c {
            b'"' => b"&quot;",
            b'&' => b"&amp;",
            b'<' => b"&lt;",
            b'>' => b"&gt;",
            _ => continue,
        };
        out.extend_from_slice(&text[last..i]);
        out.extend_from_slice(entity);
        last = i + 1;
    }
    out.extend_from_slice(&text[last..]);
}

/// `CHTMLHelper::HTMLEncode(text)` (`html_encode_into`).
pub(crate) fn html_encode(text: &[u8]) -> Vec<u8> {
    let mut out = Vec::with_capacity(text.len() + 16);
    html_encode_into(text, &mut out);
    out
}

// NCBI reference (2.17.0): c++/src/corelib/ncbistr.cpp:5038-5052
// ```c++
// static SIZE_TYPE s_EndOfReference(const string& str, SIZE_TYPE start)
// {
//     SIZE_TYPE pos = str.find_first_not_of
//         ("#0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz",
//          start + 1);
//     if (pos != NPOS  &&  str[pos] == ';') {
//         return pos;
//     } else {
//         return start;
//     }
// }
// ```
/// The `;` that ends the character reference starting with the `&` at `start`, or `start` when it
/// is not one.
fn end_of_reference(text: &[u8], start: usize) -> usize {
    match text[start + 1..]
        .iter()
        .position(|&c| !(c == b'#' || c.is_ascii_alphanumeric()))
    {
        Some(off) if text[start + 1 + off] == b';' => start + 1 + off,
        _ => start,
    }
}

/// C's `isspace` in the "C" locale (NCBI's `_isspace`): `\v` too, which Rust's
/// `is_ascii_whitespace` leaves out.
fn c_isspace(c: u8) -> bool {
    matches!(c, b' ' | b'\t' | b'\n' | 0x0b | 0x0c | b'\r')
}

/// The scores of `NStr::WrapIt`'s break points, worst to best (ncbistr.cpp:5106-5112).
#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
enum Score {
    Forced,
    Punct,
    Comma,
    Space,
    Newline,
}

// NCBI reference (2.17.0): c++/src/objtools/align_format/align_format_util.cpp:279-294
// ```c++
// void CAlignFormatUtil::x_WrapOutputLine(string str, size_t line_len,
//                                         CNcbiOstream& out, bool html)
// {
//     list<string> string_l;
//     NStr::TWrapFlags flags = NStr::fWrap_FlatFile;
//     if (html) {
//         flags = NStr::fWrap_HTMLPre;
//         str = CHTMLHelper::HTMLEncode(str);
//     }
//     NStr::Wrap(str, line_len, string_l, flags);
//     list<string>::iterator iter = string_l.begin();
//     while(iter != string_l.end())
//     {
//        out << *iter;
//        out << "\n";
// ```
/// `x_WrapOutputLine(text, width, out, true)`: `text` HTML-encoded and wrapped by `NStr::Wrap`
/// with `fWrap_HTMLPre` (`wrap_html_pre`), each line followed by `\n`.
pub(crate) fn write_html_wrapped<W: Write + ?Sized>(
    writer: &mut W,
    text: &[u8],
    width: usize,
) -> io::Result<()> {
    let encoded = html_encode(text);
    let mut out = Vec::with_capacity(encoded.len() + 8);
    wrap_html_pre(&encoded, width, &mut out);
    crate::cli::insert(writer, &out)
}

// NCBI reference (2.17.0): c++/src/corelib/ncbistr.cpp:5127-5339 (`NStr::WrapIt`), the branches
// that `fWrap_HTMLPre` without `fWrap_FlatFile` and without `fWrap_Hyphenate`, with an empty
// prefix, takes:
// ```c++
//         for (SIZE_TYPE pos2 = pos0;  pos2 < len  &&  column <= width;
//             ++pos2, ++column) {
//             EScore    score = eForced;
//             SIZE_TYPE score_pos = pos2;
//             const char      c = str[pos2];
//             if (c == '\n') {
//                 best_pos = pos2;
//                 best_score = eNewline;
//                 best_link = latest_link;
//                 break;
//             }
//             else if (_isspace((unsigned char)c)) {
//                 if (!do_flat  &&  pos2 > 0 &&
//                     _isspace((unsigned char)str[pos2 - 1])) {
//                     if (pos2 < len - 1 && str[pos2 + 1] == '\b') {
//                         thisPartHasBackspace = true;
//                     }
//                     continue; // take the first space of a group
//                 }
//                 score = eSpace;
//             }
//             else if (is_html && c == '<') { ... tags as zero-width ... }
//             else if (is_html && c == '&') {
//                 // ...and references as single characters
//                 pos2 = s_EndOfReference(str, pos2);
//                 if (pos2 == NPOS) {
//                     break;
//                 }
//             }
//             else if (c == ','  && column < width && score_pos < len - 1) {
//                 score = eComma;
//                 ++score_pos;
//             }
//             else if (do_flat ? c == '-' : ispunct((unsigned char)c)) {
//                 switch (c) {
//                 case '(': case '[': case '{': case '<': case '`':
//                     score = ePunct;
//                     break;
//                 default:
//                     if (score_pos < len - 1 && column < width) {
//                         score = ePunct;
//                         ++score_pos;
//                     }
//                     break;
//                 }
//             }
//             if (score >= best_score  &&  score_pos > pos0) {
//                 best_pos = score_pos;
//                 best_score = score;
//                 best_link = latest_link;
//             }
//             while (pos2 < len - 1 && str[pos2 + 1] == '\b') {
//                 ++pos2;
//                 if (column > column0) {
//                     --column;
//                 }
//                 thisPartHasBackspace = true;
//             }
//         }
//         if (best_score != eNewline  &&  column <= width) {
//             if (best_pos != len) {
//                 best_pos = len;
//                 best_link = latest_link;
//                 thisPartHasBackspace = true;
//             }
//         }
//         ... (the part [pos, best_pos), backspaces eaten, :5274-5302) ...
//         pos = best_pos;
//         ...
//             if (best_score == eSpace || best_score == eNewline) {
//                 ++pos;
//             }
//         while (pos < len  &&  str[pos] == '\b') {
//             ++pos;
//         }
// ```
// The text is HTML-encoded (`x_WrapOutputLine`), so it holds no `<`: the tag and link branches
// (:5186-5216, :5142-5146, :5306-5308) never run. A reference (`&amp;`) is one column and no
// break point.
/// `NStr::Wrap(encoded, width, lines, fWrap_HTMLPre)`, each line followed by `\n`, appended to
/// `out`. `encoded` is HTML-encoded (no `<`).
pub(crate) fn wrap_html_pre(encoded: &[u8], width: usize, out: &mut Vec<u8>) {
    let str = encoded;
    let len = str.len();
    debug_assert!(
        !str.contains(&b'<'),
        "fWrap_HTMLPre text must be HTML-encoded"
    );
    let mut pos = 0;
    while pos < len {
        let mut column = 0usize;
        let column0 = column;
        let mut best_pos = usize::MAX;
        let mut best_score = Score::Forced;
        let mut has_backspace = false;
        let pos0 = pos;
        let mut pos2 = pos0;
        while pos2 < len && column <= width {
            let mut score = Score::Forced;
            let mut score_pos = pos2;
            let c = str[pos2];
            if c == b'\n' {
                best_pos = pos2;
                best_score = Score::Newline;
                break;
            } else if c_isspace(c) {
                if pos2 > 0 && c_isspace(str[pos2 - 1]) {
                    if pos2 + 1 < len && str[pos2 + 1] == 8 {
                        has_backspace = true;
                    }
                    pos2 += 1;
                    column += 1;
                    continue;
                }
                score = Score::Space;
            } else if c == b'&' {
                pos2 = end_of_reference(str, pos2);
            } else if c == b',' && column < width && score_pos + 1 < len {
                score = Score::Comma;
                score_pos += 1;
            } else if c.is_ascii_punctuation() {
                match c {
                    b'(' | b'[' | b'{' | b'<' | b'`' => score = Score::Punct,
                    _ => {
                        if score_pos + 1 < len && column < width {
                            score = Score::Punct;
                            score_pos += 1;
                        }
                    }
                }
            }
            if score >= best_score && score_pos > pos0 {
                best_pos = score_pos;
                best_score = score;
            }
            while pos2 + 1 < len && str[pos2 + 1] == 8 {
                pos2 += 1;
                column = column.saturating_sub(usize::from(column > column0));
                has_backspace = true;
            }
            pos2 += 1;
            column += 1;
        }
        if best_score != Score::Newline && column <= width && best_pos != len {
            best_pos = len;
            has_backspace = true;
        }
        let line_start = out.len();
        append_eating_backspaces(&str[pos..best_pos], has_backspace, line_start, out);
        out.push(b'\n');
        pos = best_pos;
        if matches!(best_score, Score::Space | Score::Newline) {
            pos += 1;
        }
        while pos < len && str[pos] == 8 {
            pos += 1;
        }
    }
}

// NCBI reference (2.17.0): c++/src/corelib/ncbistr.cpp:5277-5301
// ```c++
//             if (thisPartHasBackspace) {
//                 string::const_iterator bs; // position of next backspace
//                 while ((bs = find(begin, end, '\b')) != end) {
//                     if (bs != begin) {
//                         temp_back.append(begin, bs - 1);
//                     }
//                     else {
//                         SIZE_TYPE size = temp_back.size();
//                         if (size > prefix1->size()) { // current size > prefix size
//                             temp_back.resize(size - 1);
//                         }
//                     }
//                     begin = bs + 1;
//                 }
//             }
//             if (begin != end) {
//                 temp_back.append(begin, end);
//             }
// ```
/// The part of a line, each backspace eating the byte before it (`line_start` is where the line
/// starts in `out`, NCBI's empty prefix).
fn append_eating_backspaces(
    part: &[u8],
    has_backspace: bool,
    line_start: usize,
    out: &mut Vec<u8>,
) {
    let mut begin = 0;
    if has_backspace {
        while let Some(off) = part[begin..].iter().position(|&c| c == 8) {
            let bs = begin + off;
            if bs != begin {
                out.extend_from_slice(&part[begin..bs - 1]);
            } else if out.len() > line_start {
                out.pop();
            }
            begin = bs + 1;
        }
    }
    out.extend_from_slice(&part[begin..]);
}

#[cfg(test)]
#[path = "html_tests.rs"]
mod tests;
