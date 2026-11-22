//! The HTML form of the pairwise and organism reports (`-html`, NCBI's `CBlastFormat::m_IsHTML`):
//! NCBI's `CHTMLHelper::HTMLEncode` and its `NStr::Wrap` with `fWrap_HTMLPre`, which the report
//! writers call only when `-html` is set, once per query or per subject (the plain path does not
//! reach this module).

use std::io::{self, Write};

// NCBI reference (598d8ae6): c++/src/html/htmlhelper.cpp:146-149
// NCBI reference (598d8ae6): c++/src/html/htmlhelper.cpp:58-130 (`s_HTMLEncode`): `"` -> `&quot;`,
// `&` -> `&` + `amp;` (the default `fEncodeAll` sets no `fSkipEntities`), `<` -> `&lt;`,
// `>` -> `&gt;`; every other byte (`'` too) is copied.
/// `CHTMLHelper::HTMLEncode` of `text` with its default flags (`fEncodeAll`), appended to `out`.
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

/// `CHTMLHelper::HTMLEncode` of `text` (`html_encode_into`).
pub(crate) fn html_encode(text: &[u8]) -> Vec<u8> {
    let mut out = Vec::with_capacity(text.len() + 16);
    html_encode_into(text, &mut out);
    out
}

// NCBI reference (598d8ae6): c++/src/corelib/ncbistr.cpp:5038-5052
/// The `;` that ends the character reference starting with the `&` at `start` (the first byte
/// after it that is not a letter, a digit or `#`), or `start` when it is not one.
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

// NCBI reference (598d8ae6): c++/src/objtools/align_format/align_format_util.cpp:279-294
/// `x_WrapOutputLine` of `text` at `width` with HTML on: `text` HTML-encoded and wrapped by
/// `NStr::Wrap` with `fWrap_HTMLPre` (`wrap_html_pre`), each line followed by `\n`.
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

// NCBI reference (598d8ae6): c++/src/corelib/ncbistr.cpp:5127-5339 (`NStr::WrapIt`), the branches
// that `fWrap_HTMLPre` without `fWrap_FlatFile` and without `fWrap_Hyphenate`, with an empty
// prefix, takes. Each line looks for its break within `width` columns: a newline ends it; else
// the best-scored break point (`Score`) wins, the last one of that score: a space (not one after
// a space), a comma (the break after it), other punctuation (before an opening bracket or a
// backquote, after the rest); none gives a forced break at the width. A backspace takes back a
// column, and a line that reaches the end of the text within the width takes the rest. A space
// or newline at the break is dropped, as are backspaces at the start of the next line.
// The text is HTML-encoded (`x_WrapOutputLine`), so it holds no `<`: the tag and link branches
// (:5186-5216, :5142-5146, :5306-5308) never run. A reference (`&amp;`) is one column and no
// break point.
/// `NStr::Wrap` of `encoded` at `width` with the flag `fWrap_HTMLPre`, each line followed by `\n`,
/// appended to `out`. `encoded` is HTML-encoded (no `<`).
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

// NCBI reference (598d8ae6): c++/src/corelib/ncbistr.cpp:5277-5301
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
