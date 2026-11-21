//! LOSAT's own words where NCBI's messages name NCBI's product or show the internals of NCBI's
//! C++ Toolkit (`PD-LOSAT-CLI-NONSEARCH-DIFFERENCES` exception 12; the Owner, 2026-11-20).
//!
//! This module is the one producer of those words: the user-manual footer of the
//! `BLAST ... error:` messages ([`USER_MANUAL_FOOTER`]) and the text of an NCBI exception report
//! ([`losat_text`]: the message of each exception, without NCBI's `NCBI C++ Exception:` header,
//! its thread, the source file of NCBI's build tree and the line, the class and error code, and
//! NCBI's C++ function). The prefixes `BLAST query/options error:`, `BLAST engine error:`,
//! `BLAST Database error:` and `BLAST options error:`, NCBI's diagnostic codes such as `(106.15)`,
//! the message bodies, the exit statuses and every search output stay NCBI's.
//!
//! [`losat_text`] is also the one mapping from NCBI's frozen bytes to LOSAT's (AGENTS.md rule 8):
//! LOSAT's message is `losat_text` of NCBI's, so every comparison of LOSAT's messages with a
//! freeze of NCBI's applies it to the freeze (`losat_bytes`; the Python harnesses through
//! `LOSAT/tests/losat_wording.py`, which reads the constants below from this file).
//!
//! The rules, each applied where it matches, left to right:
//! - [`NCBI_USER_MANUAL_FOOTER`] becomes [`USER_MANUAL_FOOTER`].
//! - `NCBI C++ Exception:`, the white space after it and the place of the first exception
//!   (`T0 "<file>", line <n>: <severity>: (<class>::<code>) ` and NCBI's function
//!   `<function>() - ` when the report names one) are left out: the line keeps the message.
//! - The place of each further exception of the chain (after a space: its indentation, or the
//!   space of a message whose white space is collapsed, as in an archive's message) is left out
//!   the same way, so each further message keeps its indentation.
//! - Any other `(<class>::<code>) ` (and the function after it) is left out after a colon or at
//!   the start of a line; elsewhere it becomes `: `.
//! - `(m_Pos = <n>)` becomes `(position <n>)`.
//!
//! Each rule keeps its text on its line and maps a collapsed space as it maps a line break with
//! an indentation, so the mapping of an archive's message (white space collapsed,
//! `report::archive::error_lines`) is the collapsed mapping of the message.

use std::borrow::Cow;

/// The footer of LOSAT's `BLAST query/options error:` messages (NCBI's [`NCBI_USER_MANUAL_FOOTER`]).
pub const USER_MANUAL_FOOTER: &str = "Please refer to the LOSAT user manual.";

/// NCBI's footer of its query and options errors, which LOSAT replaces with [`USER_MANUAL_FOOTER`].
pub const NCBI_USER_MANUAL_FOOTER: &str = "Please refer to the BLAST+ user manual.";

/// The header of NCBI's exception report (`CException::ReportAll`).
pub const NCBI_EXCEPTION_HEADER: &str = "NCBI C++ Exception:";

/// NCBI's text of a `CStringException`'s position (`ReportExtra`), and LOSAT's.
pub const NCBI_POSITION: &str = "(m_Pos = ";
pub const POSITION: &str = "(position ";

/// The classes of NCBI's serial exceptions whose names do not end with `Exception`.
pub const SERIAL_EXCEPTION_CLASSES: &[&str] = &["CInvalidChoiceSelection", "CUnassignedMember"];

/// LOSAT's text of `text`, NCBI's (the module's rules).
pub fn losat_text(text: &str) -> Cow<'_, str> {
    match losat_bytes(text.as_bytes()) {
        Cow::Borrowed(_) => Cow::Borrowed(text),
        // The rules remove or replace ASCII runs only, so UTF-8 stays UTF-8.
        Cow::Owned(bytes) => Cow::Owned(String::from_utf8(bytes).expect("UTF-8 in, UTF-8 out")),
    }
}

/// [`losat_text`] of bytes.
pub fn losat_bytes(text: &[u8]) -> Cow<'_, [u8]> {
    let footer = NCBI_USER_MANUAL_FOOTER.as_bytes();
    let header = NCBI_EXCEPTION_HEADER.as_bytes();
    let position = NCBI_POSITION.as_bytes();
    if !(contains(text, footer)
        || contains(text, header)
        || contains(text, position)
        || contains(text, b"::"))
    {
        return Cow::Borrowed(text);
    }
    let mut out = Vec::with_capacity(text.len());
    let mut i = 0;
    let mut changed = false;
    while i < text.len() {
        let rest = &text[i..];
        if rest.starts_with(footer) {
            out.extend_from_slice(USER_MANUAL_FOOTER.as_bytes());
            i += footer.len();
            changed = true;
            continue;
        }
        if rest.starts_with(header) {
            let mut j = header.len();
            while matches!(rest.get(j), Some(b' ' | b'\n')) {
                j += 1;
            }
            if j > header.len() {
                if let Some(place) = place(&rest[j..]) {
                    // A report glued to the text before it (`Object type= BlastOutput2NCBI C++
                    // Exception:`) is set off by `: `, as a tag inside a line is; after a space
                    // (`[blastn] NCBI C++ Exception:`) nothing stands for it.
                    if out
                        .last()
                        .is_some_and(|&b| !matches!(b, b' ' | b'\n' | b':'))
                    {
                        glue(&mut out);
                    }
                    i += j + place;
                    changed = true;
                    continue;
                }
            }
        }
        if rest[0] == b'T' && out.last() == Some(&b' ') {
            if let Some(place) = place(rest) {
                i += place;
                changed = true;
                continue;
            }
        }
        if rest[0] == b'(' {
            if let Some(tag) = tag(rest) {
                glue(&mut out);
                i += tag;
                changed = true;
                continue;
            }
        }
        if rest.starts_with(position) {
            out.extend_from_slice(POSITION.as_bytes());
            i += position.len();
            changed = true;
            continue;
        }
        out.push(rest[0]);
        i += 1;
    }
    if changed {
        Cow::Owned(out)
    } else {
        Cow::Borrowed(text)
    }
}

/// The member of an ASN.1 text archive that holds a message (`Blast4-error.message`).
const ARCHIVE_MESSAGE: &[u8] = b"message \"";

/// [`losat_bytes`] of a whole output (standard output, an `-out` file, or both merged): the
/// messages of an ASN.1 text archive (`message "..."`, which the archive writer writes in LOSAT's
/// words) are read without the ends of line that wrap them, mapped, and written again as the
/// writer writes them (`report::serial::asn_text::visible_string_line`); the rest is mapped as
/// text.
pub fn losat_output_bytes(output: &[u8]) -> Cow<'_, [u8]> {
    if !contains(output, ARCHIVE_MESSAGE) {
        return losat_bytes(output);
    }
    let mut out = Vec::with_capacity(output.len());
    let mut done = 0;
    let mut search = 0;
    while let Some(found) = find(&output[search..], ARCHIVE_MESSAGE) {
        let at = search + found;
        let line_start = output[..at]
            .iter()
            .rposition(|&b| b == b'\n')
            .map_or(0, |i| i + 1);
        let quote = at + ARCHIVE_MESSAGE.len() - 1;
        search = quote + 1;
        // The member starts its line (after the indentation).
        if !output[line_start..at].iter().all(|&b| b == b' ') {
            continue;
        }
        let Some((value, end)) = asn_string(output, quote) else {
            continue;
        };
        search = end;
        let mapped = losat_bytes(&value);
        if let Cow::Owned(mapped) = mapped {
            out.extend_from_slice(&losat_bytes(&output[done..line_start]));
            out.extend_from_slice(&crate::report::serial::asn_text::visible_string_line(
                &output[line_start..quote],
                &mapped,
            ));
            done = end;
        }
    }
    if done == 0 {
        return losat_bytes(output);
    }
    out.extend_from_slice(&losat_bytes(&output[done..]));
    Cow::Owned(out)
}

/// The value of the ASN.1 text string whose opening quote is at `quote` (its doubled quotes
/// read as one, the ends of line of its wraps left out) and the end of its closing quote.
fn asn_string(text: &[u8], quote: usize) -> Option<(Vec<u8>, usize)> {
    let mut value = Vec::new();
    let mut i = quote + 1;
    loop {
        match *text.get(i)? {
            b'"' if text.get(i + 1) == Some(&b'"') => {
                value.push(b'"');
                i += 2;
            }
            b'"' => return Some((value, i + 1)),
            b'\n' => i += 1,
            byte => {
                value.push(byte);
                i += 1;
            }
        }
    }
}

fn find(text: &[u8], needle: &[u8]) -> Option<usize> {
    text.windows(needle.len())
        .position(|window| window == needle)
}

fn contains(text: &[u8], needle: &[u8]) -> bool {
    text.windows(needle.len()).any(|window| window == needle)
}

/// What stands for a removed report or tag after the text before it: nothing at the start, after
/// a line feed or after `:` (spaces after it kept); otherwise `: ` in place of a space before it.
fn glue(out: &mut Vec<u8>) {
    let kept = out.iter().rposition(|&b| b != b' ');
    if kept.is_some_and(|k| !matches!(out[k], b':' | b'\n')) {
        if out.last() == Some(&b' ') {
            out.pop();
        }
        out.extend_from_slice(b": ");
    }
}

/// The length of the place of an exception at the start of `text`:
/// `T<thread> "<file>", line <n>: <severity>: ` and its [`tag`].
fn place(text: &[u8]) -> Option<usize> {
    let mut i = literal(text, b"T")?;
    i += digits(&text[i..])?;
    i += literal(&text[i..], b" \"")?;
    i += text[i..].iter().position(|&b| b == b'"' || b == b'\n')?;
    i += literal(&text[i..], b"\", line ")?;
    i += digits(&text[i..])?;
    i += literal(&text[i..], b": ")?;
    let word = text[i..]
        .iter()
        .take_while(|b| b.is_ascii_alphabetic())
        .count();
    if word == 0 {
        return None;
    }
    i += word;
    i += literal(&text[i..], b": ")?;
    Some(i + tag(&text[i..])?)
}

/// The length of `(<class>::<code>) ` at the start of `text` and of the function after it
/// (`<function>() - `), where `<class>` is a class of NCBI's exceptions.
fn tag(text: &[u8]) -> Option<usize> {
    let mut i = literal(text, b"(")?;
    let class = text[i..]
        .iter()
        .take_while(|b| b.is_ascii_alphanumeric() || **b == b'_')
        .count();
    let name = &text[i..i + class];
    if !(name.ends_with(b"Exception")
        || SERIAL_EXCEPTION_CLASSES
            .iter()
            .any(|serial| serial.as_bytes() == name))
    {
        return None;
    }
    i += class;
    i += literal(&text[i..], b"::")?;
    i += text[i..].iter().position(|&b| b == b')' || b == b'\n')?;
    i += literal(&text[i..], b") ")?;
    let function = text[i..]
        .iter()
        .take_while(|b| b.is_ascii_alphanumeric() || matches!(b, b'_' | b':' | b'~'))
        .count();
    if function > 0 {
        if let Some(end) = literal(&text[i + function..], b"() - ") {
            i += function + end;
        }
    }
    Some(i)
}

fn digits(text: &[u8]) -> Option<usize> {
    let count = text.iter().take_while(|b| b.is_ascii_digit()).count();
    (count > 0).then_some(count)
}

fn literal(text: &[u8], expected: &[u8]) -> Option<usize> {
    text.starts_with(expected).then_some(expected.len())
}

#[cfg(test)]
#[path = "wording_tests.rs"]
pub(crate) mod tests;
