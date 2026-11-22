//! LOSAT's own words where NCBI's programs write NCBI's: the version statement, the prefixes of
//! the error reports, and NCBI's messages that name NCBI's product or show the internals of NCBI's
//! C++ Toolkit (`PD-LOSAT-CLI-NONSEARCH-DIFFERENCES` exceptions 12 and 14; the Owner, 2026-11-20
//! and 2026-11-21).
//!
//! This module is the one producer of those words:
//! - LOSAT's version ([`LOSAT_VERSION`], `LOSAT 0.2.0`, from the crate's version) where NCBI's
//!   programs write `2.17.0+` alone (an archive's `ident`, SAM's `@PG VN:`, `-version`), and
//!   [`REPORT_VERSION`] (`(LOSAT 0.2.0)`) where NCBI's follows the program's name: a report's
//!   header `BLASTP (LOSAT 0.2.0)` (outfmt 0 and 7, `BlastOutput_version`, the JSON and XML2
//!   `version`), NCBI's `BLASTP 2.17.0+`. The program's name stays first, so readers that take the
//!   first word as the program keep reading it.
//! - The prefixes of the error reports ([`ERROR_PREFIXES`]: `LOSAT engine error:` and the others,
//!   NCBI's `BLAST engine error:` ...).
//! - The user-manual footer of the query and options errors ([`USER_MANUAL_FOOTER`]) and the text
//!   of an NCBI exception report ([`losat_text`]: the message of each exception, without NCBI's
//!   `NCBI C++ Exception:` header, its thread, the source file of NCBI's build tree and the line,
//!   the class and error code, and NCBI's C++ function).
//!
//! NCBI's diagnostic codes such as `(106.15)`, the message bodies, the exit statuses and every
//! other byte of the search outputs stay NCBI's.
//!
//! [`losat_text`] is also the one mapping from NCBI's frozen bytes to LOSAT's (AGENTS.md rule 8):
//! LOSAT's text is `losat_text` of NCBI's, so every comparison of LOSAT's output or messages with
//! a freeze of NCBI's applies it to the freeze (`losat_output_bytes`, `losat_bytes`; the Python
//! harnesses through `LOSAT/tests/losat_wording.py`, which reads the constants below from this
//! file and LOSAT's version from `LOSAT/Cargo.toml`; on the command line
//! `python3 -I LOSAT/tests/losat_wording.py map-output <file>`).
//!
//! The rules, each applied where it matches, left to right:
//! - NCBI's version [`NCBI_VERSION`] (`2.17.0+`) as the VisibleString of a binary (BER) archive's
//!   `ident` (its tag and length byte before it) becomes [`LOSAT_VERSION`] with LOSAT's length byte.
//!   After a start of line, a tab, `"`, `>`, `:` or a
//!   space it becomes [`REPORT_VERSION`] when a program's name in capitals (`BLASTP`, `PSITBLASTN`,
//!   `BLAST`) and a space come before it, and [`LOSAT_VERSION`] otherwise (`blastn: 2.17.0+`,
//!   `VN:2.17.0+`, `ident "2.17.0+"`, `version 2.17.0+`); rpsbproc's banner
//!   [`RPSBPROC_TITLE`]` v`[`NCBI_RPSBPROC_VERSION`] ends with ` `[`REPORT_VERSION`] instead.
//! - The version of `-xmlhelp`'s program (`<version>2.17.0</version>` after `<program type="...">`
//!   and its `<name>`; dustmasker's and segmasker's `1.0.0`, rpsbproc's `0.5.1`) becomes
//!   [`VERSION_NUMBER`], the number that `-version-full-xml` gives too.
//! - Each NCBI prefix of [`ERROR_PREFIXES`] becomes LOSAT's.
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
//! [`losat_output_bytes`], the mapping of a whole frozen output or standard error, also moves the
//! position of each byte's post of a binary archive (`<trace>` `Blast4-archive...` or
//! `Blast4-request...`, then `Bad char [0x<hex>] in string at byte <n>`) by the bytes that LOSAT's
//! `ident` adds ([`IDENT_GROWTH`]), since the `ident` comes before every string of an archive. This
//! rule moves a number, so it is not one of [`losat_text`]'s, whose text of LOSAT's maps to itself.
//!
//! Each rule keeps its text on its line and maps a collapsed space as it maps a line break with
//! an indentation, so the mapping of an archive's message (white space collapsed,
//! `report::archive::error_lines`) is the collapsed mapping of the message. LOSAT's prefixes have
//! the length of NCBI's, so an archive's message keeps NCBI's wraps.

use std::borrow::Cow;

// `losat_version!()`: LOSAT's version where NCBI's programs write theirs alone (`LOSAT 0.2.0`).
// A macro, so that `concat!` builds `&'static str` constants with it.
macro_rules! losat_version {
    () => {
        concat!("LOSAT ", env!("CARGO_PKG_VERSION"))
    };
}
// `report_version!()`: LOSAT's version after a program's name (`BLASTP (LOSAT 0.2.0)`).
macro_rules! report_version {
    () => {
        concat!("(", $crate::wording::losat_version!(), ")")
    };
}
pub(crate) use losat_version;
pub(crate) use report_version;

/// LOSAT's version where NCBI's programs write `2.17.0+` (or their own version) alone: an
/// archive's or a search strategy's `ident`, SAM's `@PG VN:`, the first line of `-version`
/// (`blastn: LOSAT 0.2.0`).
pub const LOSAT_VERSION: &str = losat_version!();

/// LOSAT's version after a program's name: a report's header `BLASTP (LOSAT 0.2.0)` (NCBI's
/// `BLASTP 2.17.0+`), which the writers join as `<PROGRAM> <version>`.
pub const REPORT_VERSION: &str = report_version!();

/// LOSAT's version number (`0.2.0`) where NCBI's `-xmlhelp` writes its program's
/// (`<version>2.17.0</version>`): the element holds a number, as `-version-full-xml`'s
/// `version_info` does.
pub const VERSION_NUMBER: &str = env!("CARGO_PKG_VERSION");

/// What NCBI's `CBlastVersion::Print()` writes (`2.17.0+`), which LOSAT's text replaces with
/// [`LOSAT_VERSION`] or [`REPORT_VERSION`].
pub const NCBI_VERSION: &str = crate::NCBI_BLAST_VERSION;

/// The bytes that LOSAT's `ident` ([`LOSAT_VERSION`]) adds to an archive, NCBI's being
/// [`NCBI_VERSION`]: the positions of the bytes after it in a binary archive move by as many.
pub const IDENT_GROWTH: usize = LOSAT_VERSION.len() - NCBI_VERSION.len();

/// The name of rpsbproc in its banner (NCBI's and LOSAT's), and NCBI's version of it (`VERDEF`):
/// NCBI's banner is `<title> v0.5.1`, LOSAT's `<title> (LOSAT 0.2.0)` ([`REPORT_VERSION`]).
pub const RPSBPROC_TITLE: &str = "Post-RPSBLAST Processing Utility";
pub const NCBI_RPSBPROC_VERSION: &str = "0.5.1";

/// LOSAT's prefixes of the error reports (NCBI's `blast_app_util.hpp` `CATCH_ALL` and the
/// readers' reports), each after NCBI's in [`ERROR_PREFIXES`].
pub const QUERY_OPTIONS_ERROR: &str = "LOSAT query/options error:";
pub const ENGINE_ERROR: &str = "LOSAT engine error:";
pub const DATABASE_ERROR: &str = "LOSAT database error:";
pub const DATABASE_CREATION_ERROR: &str = "LOSAT database creation error:";
pub const OPTIONS_ERROR: &str = "LOSAT options error:";
pub const QUERY_ERROR: &str = "LOSAT query error:";
pub const NCBI_QUERY_OPTIONS_ERROR: &str = "BLAST query/options error:";
pub const NCBI_ENGINE_ERROR: &str = "BLAST engine error:";
pub const NCBI_DATABASE_ERROR: &str = "BLAST Database error:";
pub const NCBI_DATABASE_CREATION_ERROR: &str = "BLAST Database creation error:";
pub const NCBI_OPTIONS_ERROR: &str = "BLAST options error:";
pub const NCBI_QUERY_ERROR: &str = "BLAST query error:";

/// NCBI's prefixes and LOSAT's (`BLAST engine error:` -> `LOSAT engine error:`).
pub const ERROR_PREFIXES: [(&str, &str); 6] = [
    (NCBI_QUERY_OPTIONS_ERROR, QUERY_OPTIONS_ERROR),
    (NCBI_ENGINE_ERROR, ENGINE_ERROR),
    (NCBI_DATABASE_ERROR, DATABASE_ERROR),
    (NCBI_DATABASE_CREATION_ERROR, DATABASE_CREATION_ERROR),
    (NCBI_OPTIONS_ERROR, OPTIONS_ERROR),
    (NCBI_QUERY_ERROR, QUERY_ERROR),
];

/// The footer of LOSAT's `LOSAT query/options error:` messages (NCBI's [`NCBI_USER_MANUAL_FOOTER`]).
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
    let version = NCBI_VERSION.as_bytes();
    let title = RPSBPROC_TITLE.as_bytes();
    if !(contains(text, footer)
        || contains(text, header)
        || contains(text, position)
        || contains(text, b"::")
        || contains(text, version)
        || contains(text, b" error:")
        || contains(text, title)
        || contains(text, XMLHELP_PROGRAM))
    {
        return Cow::Borrowed(text);
    }
    let mut out = Vec::with_capacity(text.len());
    let mut i = 0;
    let mut changed = false;
    while i < text.len() {
        let rest = &text[i..];
        // A binary (BER) archive's `ident`: a VisibleString (tag 26) of NCBI's length, whose length
        // byte becomes LOSAT's.
        if rest.starts_with(version) && out.ends_with(&[BER_VISIBLE_STRING, version.len() as u8]) {
            out.pop();
            out.push(LOSAT_VERSION.len() as u8);
            out.extend_from_slice(LOSAT_VERSION.as_bytes());
            i += version.len();
            changed = true;
            continue;
        }
        if rest.starts_with(version)
            && matches!(
                out.last(),
                None | Some(b'\n' | b'\t' | b'"' | b'>' | b':' | b' ')
            )
        {
            let losat = if after_program_name(&out) {
                REPORT_VERSION
            } else {
                LOSAT_VERSION
            };
            out.extend_from_slice(losat.as_bytes());
            i += version.len();
            changed = true;
            continue;
        }
        if let Some((ncbi, losat)) = ERROR_PREFIXES
            .iter()
            .find(|(ncbi, _)| rest[0] == b'B' && rest.starts_with(ncbi.as_bytes()))
        {
            out.extend_from_slice(losat.as_bytes());
            i += ncbi.len();
            changed = true;
            continue;
        }
        if let Some(banner) = rpsbproc_banner(rest) {
            out.extend_from_slice(title);
            out.push(b' ');
            out.extend_from_slice(REPORT_VERSION.as_bytes());
            i += banner;
            changed = true;
            continue;
        }
        if rest[0] == b'<' {
            if let Some(number) = xmlhelp_version(&out, rest) {
                out.extend_from_slice(b"<version>");
                out.extend_from_slice(VERSION_NUMBER.as_bytes());
                out.extend_from_slice(b"</version>");
                changed |= &rest[9..number - 10] != VERSION_NUMBER.as_bytes();
                i += number;
                continue;
            }
        }
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
    let words = losat_output_words(output);
    match archive_positions(&words) {
        Some(moved) => Cow::Owned(moved),
        None => words,
    }
}

/// [`losat_output_bytes`] without the positions of the byte posts: the mapping of a frozen stream
/// whose bytes LOSAT's run keeps at NCBI's offsets (a test that balances the longer version with a
/// shorter title before the posts).
pub fn losat_output_words(output: &[u8]) -> Cow<'_, [u8]> {
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

/// `output` with the position of each byte's post of a binary archive moved by [`IDENT_GROWTH`],
/// or `None` when it holds no such post.
fn archive_positions(output: &[u8]) -> Option<Vec<u8>> {
    let mut out = Vec::new();
    let mut done = 0;
    let mut search = 0;
    while let Some(found) = find(&output[search..], BYTE_POSITION) {
        let at = search + found;
        let number_at = at + BYTE_POSITION.len();
        search = number_at;
        if !archive_byte_post(&output[..at]) {
            continue;
        }
        let Some(n) = digits(&output[number_at..]) else {
            continue;
        };
        let Some(number) = std::str::from_utf8(&output[number_at..number_at + n])
            .ok()
            .and_then(|text| text.parse::<u64>().ok())
        else {
            continue;
        };
        out.extend_from_slice(&output[done..number_at]);
        out.extend_from_slice((number + IDENT_GROWTH as u64).to_string().as_bytes());
        done = number_at + n;
    }
    if done == 0 {
        return None;
    }
    out.extend_from_slice(&output[done..]);
    Some(out)
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

/// The BER tag of a VisibleString (universal 26), before its length byte in a binary archive.
const BER_VISIBLE_STRING: u8 = 0x1a;

/// The words before a position in a binary stream's post of a byte (`report::serial::encoder::ReplacedChar`).
const BYTE_POSITION: &[u8] = b" in string at byte ";

/// Whether `out` ends with `Bad char [0x<hex>]` on the line after the trace of an archive
/// (`... Blast4-archive.request.body...` or `... Blast4-request...`), whose `ident` comes first.
fn archive_byte_post(out: &[u8]) -> bool {
    let Some(line_end) = out.iter().rposition(|&b| b == b'\n') else {
        return false;
    };
    if !out[line_end + 1..].starts_with(b"Bad char [0x") {
        return false;
    }
    let trace = &out[..line_end];
    let word = trace
        .iter()
        .rev()
        .take_while(|b| b.is_ascii_alphanumeric() || matches!(b, b'.' | b'-' | b'_'))
        .count();
    let word = &trace[trace.len() - word..];
    word.starts_with(b"Blast4-archive") || word.starts_with(b"Blast4-request")
}

/// Whether `out` ends with a program's name in capitals and a space (`BLASTP `, `# PSITBLASTN `,
/// `Specific BLAST `): the word before the space has a capital and no small letter.
fn after_program_name(out: &[u8]) -> bool {
    let Some(before) = out.strip_suffix(b" ") else {
        return false;
    };
    let word = before
        .iter()
        .rev()
        .take_while(|b| b.is_ascii_alphanumeric() || matches!(b, b'_' | b'-'))
        .count();
    let word = &before[before.len() - word..];
    word.iter().any(u8::is_ascii_uppercase) && !word.iter().any(u8::is_ascii_lowercase)
}

/// The start of the program's element in `-xmlhelp` (`ncbi_xmlhelp`).
const XMLHELP_PROGRAM: &[u8] = b"<program type=\"";

/// The length of `<version>V</version>` (V digits and dots) at the start of `text` when it is
/// the version of `-xmlhelp`'s program: `out` ends with the lines `<program type="...">` and
/// `<name>...</name>` (each line feed a space in a collapsed message).
fn xmlhelp_version(out: &[u8], text: &[u8]) -> Option<usize> {
    let line_end = |text: &[u8]| matches!(text.last(), Some(b'\n' | b' '));
    if !line_end(out) {
        return None;
    }
    let before = out[..out.len() - 1].strip_suffix(b"</name>")?;
    let name = rfind(before, b"<name>")?;
    if before[name..].contains(&b'\n')
        || before[name + 1..].contains(&b'<')
        || !line_end(&before[..name])
    {
        return None;
    }
    let before = before[..name - 1].strip_suffix(b">")?;
    let program = rfind(before, XMLHELP_PROGRAM)?;
    if (program > 0 && !line_end(&before[..program]))
        || before[program + XMLHELP_PROGRAM.len()..]
            .iter()
            .any(|b| matches!(b, b'\n' | b'<' | b'>'))
    {
        return None;
    }
    let mut i = literal(text, b"<version>")?;
    let number = text[i..]
        .iter()
        .take_while(|b| b.is_ascii_digit() || **b == b'.')
        .count();
    if number == 0 {
        return None;
    }
    i += number;
    Some(i + literal(&text[i..], b"</version>")?)
}

/// The length of NCBI's rpsbproc banner `<title> v0.5.1` at the start of `text`.
fn rpsbproc_banner(text: &[u8]) -> Option<usize> {
    let mut i = literal(text, RPSBPROC_TITLE.as_bytes())?;
    i += literal(&text[i..], b" v")?;
    Some(i + literal(&text[i..], NCBI_RPSBPROC_VERSION.as_bytes())?)
}

fn find(text: &[u8], needle: &[u8]) -> Option<usize> {
    text.windows(needle.len())
        .position(|window| window == needle)
}

fn rfind(text: &[u8], needle: &[u8]) -> Option<usize> {
    text.windows(needle.len())
        .rposition(|window| window == needle)
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
