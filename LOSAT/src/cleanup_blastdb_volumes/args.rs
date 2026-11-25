//! Command-line parser for `cleanup-blastdb-volumes.py`, reproducing Python 3.13 `argparse`'s
//! options, prefix abbreviations, terminal-width-wrapped usage/help text, error messages, and exit
//! status `2`.
//!
//! NCBI reference (598d8ae6): c++/src/app/blast/cleanup-blastdb-volumes.py:40-41,144-156

use std::ffi::OsString;

use crate::utils::os_bytes::os_str_bytes;

use super::config::py_trim;

/// Program name written in usage, help, version, and error lines.
pub(crate) const PROG: &str = "cleanup-blastdb-volumes.py";

/// Fallback terminal column count used when `COLUMNS` is unset/invalid and stdout is not a
/// terminal (Python's default 80-column terminal width).
pub(crate) const DEFAULT_COLUMNS: usize = 80;

const USAGE_PARTS: [&str; 6] = [
    "[-h]",
    "-db DB",
    "-dbtype {prot,nucl}",
    "[-dry-run]",
    "[-version]",
    "[-verbose]",
];

const DESCRIPTION_WORDS: [&str; 5] = ["Remove", "needless", "BLAST", "database", "volumes."];

const OPTION_HELP: [(&str, &[&str]); 6] = [
    (
        "-h, --help",
        &["show", "this", "help", "message", "and", "exit"],
    ),
    ("-db DB", &["BLAST", "database", "name"]),
    ("-dbtype {prot,nucl}", &["Molecule", "type"]),
    (
        "-dry-run",
        &[
            "Do", "not", "delete", "any", "files,", "just", "list", "them",
        ],
    ),
    (
        "-version",
        &["show", "program's", "version", "number", "and", "exit"],
    ),
    ("-verbose", &["Increase", "output", "verbosity"]),
];

/// Parses Python's `COLUMNS` environment value (`int(s)` when `> 0`, returning `None` on empty,
/// non-integer, or non-positive values; huge positive integers saturate to `usize::MAX`).
fn parse_columns_env(raw: &str) -> Option<usize> {
    let trimmed = py_trim(raw);
    if trimmed.is_empty() {
        return None;
    }
    let digits = if let Some(rest) = trimmed.strip_prefix('+') {
        rest
    } else if trimmed.starts_with('-') {
        return None;
    } else {
        trimmed
    };
    if digits.is_empty() || !digits.bytes().all(|b| b.is_ascii_digit()) {
        return None;
    }
    let trimmed_zeros = digits.trim_start_matches('0');
    if trimmed_zeros.is_empty() {
        return None;
    }
    Some(trimmed_zeros.parse::<usize>().unwrap_or(usize::MAX))
}

/// Resolves the terminal column width following Python's terminal-size precedence:
/// positive `COLUMNS`, else positive stdout terminal width, else `80`.
pub(crate) fn resolve_columns(columns_env: Option<&str>, stdout_columns: Option<usize>) -> usize {
    if let Some(cols) = columns_env.and_then(parse_columns_env) {
        return cols;
    }
    if let Some(cols) = stdout_columns.filter(|&c| c > 0) {
        return cols;
    }
    DEFAULT_COLUMNS
}

/// Greedily packs `words` into lines of at most `max_width` columns (each line takes at least one
/// item, matching `argparse.HelpFormatter`'s usage and `textwrap` wrapping for our short words).
fn wrap_items(words: &[&str], max_width: isize) -> Vec<String> {
    let mut lines: Vec<String> = Vec::new();
    let mut current = String::new();
    for &word in words {
        if current.is_empty() {
            current.push_str(word);
        } else if (current.len() + 1 + word.len()) as isize <= max_width {
            current.push(' ');
            current.push_str(word);
        } else {
            lines.push(std::mem::take(&mut current));
            current.push_str(word);
        }
    }
    if !current.is_empty() {
        lines.push(current);
    }
    lines
}

/// Computes `argparse.HelpFormatter`'s inner wrapping width (`columns - 2`) as signed `isize`.
fn formatter_width(columns: usize) -> isize {
    (columns.min(isize::MAX as usize) as isize) - 2
}

// NCBI reference (598d8ae6): c++/src/app/blast/cleanup-blastdb-volumes.py:144-156
/// Formats the usage banner wrapped for a terminal of `columns` columns (`width = columns - 2`),
/// matching Python 3.13 `argparse.HelpFormatter`'s usage wrapping.
pub(crate) fn usage_text(columns: usize) -> String {
    let width = formatter_width(columns);
    let prefix = format!("usage: {PROG}");
    let joined = USAGE_PARTS.join(" ");
    if (prefix.len() + 1 + joined.len()) as isize <= width {
        return format!("{prefix} {joined}\n");
    }

    let mut out = String::new();
    // `argparse.HelpFormatter` places options on the next line (indented past `"usage: "`) when
    // the prefix exceeds 75% of `width` (exact integer check: `4 * prefix.len() > 3 * width`).
    if width <= 0 || 4 * (prefix.len() as isize) > 3 * width {
        let indent = "usage: ".len();
        out.push_str(&prefix);
        out.push('\n');
        let pad = " ".repeat(indent);
        for line in wrap_items(&USAGE_PARTS, width - indent as isize) {
            out.push_str(&pad);
            out.push_str(&line);
            out.push('\n');
        }
    } else {
        let indent = prefix.len() + 1;
        let lines = wrap_items(&USAGE_PARTS, width - indent as isize);
        let pad = " ".repeat(indent);
        for (idx, line) in lines.into_iter().enumerate() {
            if idx == 0 {
                out.push_str(&prefix);
                out.push(' ');
            } else {
                out.push_str(&pad);
            }
            out.push_str(&line);
            out.push('\n');
        }
    }
    out
}

// NCBI reference (598d8ae6): c++/src/app/blast/cleanup-blastdb-volumes.py:41,144-156
/// Formats the complete `-h` / `--help` text wrapped for a terminal of `columns` columns,
/// matching Python 3.13 `argparse.HelpFormatter`'s full help layout.
pub(crate) fn help_text_for_columns(columns: usize) -> String {
    let width = formatter_width(columns);
    let mut out = usage_text(columns);
    out.push('\n');

    let desc_width = width.max(11);
    for line in wrap_items(&DESCRIPTION_WORDS, desc_width) {
        out.push_str(&line);
        out.push('\n');
    }

    out.push_str("\noptions:\n");
    let indent_inc: isize = 2;
    let max_help_position = 24_isize.min((width - 20).max(indent_inc * 2)) as usize;
    let max_action_header = OPTION_HELP
        .iter()
        .map(|(invoc, _)| 2 + invoc.len())
        .max()
        .unwrap_or(0);
    let help_position = (max_action_header + 2).min(max_help_position);
    let help_width = (width - help_position as isize).max(11);
    let help_pad = " ".repeat(help_position);

    for (invoc, help_words) in OPTION_HELP {
        let hlines = wrap_items(help_words, help_width);
        if 2 + invoc.len() + 2 <= help_position {
            let gap = help_position - 2 - invoc.len();
            out.push_str("  ");
            out.push_str(invoc);
            for _ in 0..gap {
                out.push(' ');
            }
            out.push_str(&hlines[0]);
            out.push('\n');
            for line in &hlines[1..] {
                out.push_str(&help_pad);
                out.push_str(line);
                out.push('\n');
            }
        } else {
            out.push_str("  ");
            out.push_str(invoc);
            out.push('\n');
            for line in &hlines {
                out.push_str(&help_pad);
                out.push_str(line);
                out.push('\n');
            }
        }
    }
    out
}

/// Formats the default 80-column `-h` / `--help` text.
#[cfg(test)]
pub(crate) fn help_text() -> String {
    help_text_for_columns(DEFAULT_COLUMNS)
}

/// Molecule type selected by `-dbtype`.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum DbType {
    Prot,
    Nucl,
}

impl DbType {
    /// Leading letter of BLAST database volume extensions (`p` or `n`).
    pub(crate) fn ext_letter(self) -> u8 {
        match self {
            Self::Prot => b'p',
            Self::Nucl => b'n',
        }
    }

    /// Molecule type word written in the missing-database error message.
    pub(crate) fn label(self) -> &'static str {
        match self {
            Self::Prot => "protein",
            Self::Nucl => "nucleotide",
        }
    }

    /// Whether this molecule type is protein.
    pub(crate) fn is_prot(self) -> bool {
        matches!(self, Self::Prot)
    }
}

/// Validated arguments of `cleanup-blastdb-volumes.py`.
#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) struct ParsedArgs {
    /// Raw bytes of the `-db` database name.
    pub(crate) db: Vec<u8>,
    /// Molecule type (`prot` or `nucl`).
    pub(crate) dbtype: DbType,
    /// Whether `-dry-run` was given.
    pub(crate) dry_run: bool,
    /// Number of times `-verbose` was given.
    pub(crate) verbose: usize,
}

/// Outcome of parsing the command line.
#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) enum ParseOutcome {
    Run(ParsedArgs),
    Help(String),
    Version(String),
    Error(String),
}

/// Single-dash multi-character options in parser registration order (matches NCBI's `ambiguous
/// option` listing order).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum LongOpt {
    Db,
    DbType,
    DryRun,
    Version,
    Verbose,
}

impl LongOpt {
    const ALL: [Self; 5] = [
        Self::Db,
        Self::DbType,
        Self::DryRun,
        Self::Version,
        Self::Verbose,
    ];

    fn flag(self) -> &'static str {
        match self {
            Self::Db => "-db",
            Self::DbType => "-dbtype",
            Self::DryRun => "-dry-run",
            Self::Version => "-version",
            Self::Verbose => "-verbose",
        }
    }
}

/// Version line printed by `-version`: LOSAT's version in place of NCBI's script version
/// (`PD-LOSAT-CLI-NONSEARCH-DIFFERENCES` exception 14).
pub(crate) fn version_text() -> String {
    format!("{PROG} {}\n", crate::wording::LOSAT_VERSION)
}

/// Formats an argument error with the wrapped usage banner and program error line.
fn error_text(columns: usize, message: &str) -> String {
    format!("{}{PROG}: error: {message}\n", usage_text(columns))
}

/// Builds a `ParseOutcome::Error` from `message`.
fn err_outcome(columns: usize, message: &str) -> ParseOutcome {
    ParseOutcome::Error(error_text(columns, message))
}

/// Builds a `ParseOutcome::Error` for a flag that does not accept an explicit `=value`.
fn ignored_arg(columns: usize, flag: &str, val: &[u8]) -> ParseOutcome {
    err_outcome(
        columns,
        &format!(
            "argument {flag}: ignored explicit argument {}",
            py_repr(val)
        ),
    )
}

/// Splits `token` at the first `=` into `(option_part, explicit_value)`.
fn split_eq(token: &[u8]) -> (&[u8], Option<&[u8]>) {
    match token.iter().position(|&b| b == b'=') {
        Some(eq) => (&token[..eq], Some(&token[eq + 1..])),
        None => (token, None),
    }
}

/// Whether `(opt_part, explicit_val)` matches `--help` (including `--h`, `--he`, `--hel`, `--=<v>`).
fn is_double_dash_help(opt_part: &[u8], explicit_val: Option<&[u8]>) -> bool {
    (opt_part.len() >= 3 || explicit_val.is_some()) && b"--help".starts_with(opt_part)
}

/// Consumes the required value of `-db` or `-dbtype` from `explicit_val` or the following token.
fn take_value(
    columns: usize,
    flag: &str,
    explicit_val: Option<&[u8]>,
    raw: &[Vec<u8>],
    i: &mut usize,
) -> Result<Vec<u8>, ParseOutcome> {
    if let Some(val) = explicit_val {
        *i += 1;
        return Ok(val.to_vec());
    }
    if *i + 1 < raw.len() && !is_optional_token(&raw[*i + 1]) {
        let val = raw[*i + 1].clone();
        *i += 2;
        return Ok(val);
    }
    Err(err_outcome(
        columns,
        &format!("argument {flag}: expected one argument"),
    ))
}

/// Formats `bytes` with Python `%r` string quoting (single quotes unless the string contains `'`
/// and no `"`).
fn py_repr(bytes: &[u8]) -> String {
    let text = String::from_utf8_lossy(bytes);
    let use_double = text.contains('\'') && !text.contains('"');
    let quote = if use_double { '"' } else { '\'' };
    let mut out = String::with_capacity(text.len() + 2);
    out.push(quote);
    for ch in text.chars() {
        match ch {
            '\\' => out.push_str("\\\\"),
            '\t' => out.push_str("\\t"),
            '\n' => out.push_str("\\n"),
            '\r' => out.push_str("\\r"),
            '\'' if !use_double => out.push_str("\\'"),
            c if (c as u32) < 0x20 || c as u32 == 0x7f => {
                out.push_str(&format!("\\x{:02x}", c as u32));
            }
            c => out.push(c),
        }
    }
    out.push(quote);
    out
}

/// Whether `token` matches a negative decimal number (`-1`, `-3.14`, `-.5`), which `argparse`
/// treats as a positional value rather than an option flag.
fn is_negative_number(token: &[u8]) -> bool {
    let Some(rest) = token.strip_prefix(b"-") else {
        return false;
    };
    if rest.is_empty() {
        return false;
    }
    let digits = |s: &[u8]| !s.is_empty() && s.iter().all(u8::is_ascii_digit);
    if digits(rest) {
        return true;
    }
    if let Some(dot) = rest.iter().position(|&b| b == b'.') {
        let head = &rest[..dot];
        let tail = &rest[dot + 1..];
        return (head.is_empty() || digits(head)) && digits(tail);
    }
    false
}

/// Matches a single-dash option prefix against `LongOpt::ALL`.
enum SingleDashMatch {
    Exact(LongOpt),
    Unique(LongOpt),
    Ambiguous(Vec<&'static str>),
    None,
}

fn match_single_dash(opt: &[u8]) -> SingleDashMatch {
    for candidate in LongOpt::ALL {
        if opt == candidate.flag().as_bytes() {
            return SingleDashMatch::Exact(candidate);
        }
    }
    let mut matches = Vec::new();
    for candidate in LongOpt::ALL {
        if candidate.flag().as_bytes().starts_with(opt) {
            matches.push(candidate);
        }
    }
    match matches.as_slice() {
        [] => SingleDashMatch::None,
        [only] => SingleDashMatch::Unique(*only),
        many => SingleDashMatch::Ambiguous(many.iter().map(|c| c.flag()).collect()),
    }
}

/// Builds `ParseOutcome::Help` unless an explicit value was supplied (`-h=1`, `--help=1`).
fn help_or_ignored(columns: usize, val: Option<&[u8]>) -> ParseOutcome {
    match val {
        Some(v) => ignored_arg(columns, "-h/--help", v),
        None => ParseOutcome::Help(help_text_for_columns(columns)),
    }
}

/// Whether `token` is treated as an option by Python `argparse._parse_optional` (a leading-`-`
/// token containing an ASCII space is treated as a positional value unless it matches `--help`,
/// `-h*`, or a registered option/prefix).
fn is_optional_token(token: &[u8]) -> bool {
    if token.len() < 2 || token[0] != b'-' || is_negative_number(token) {
        return false;
    }
    if token == b"--" {
        return true;
    }
    let (opt, val) = split_eq(token);
    if opt.starts_with(b"--") {
        if is_double_dash_help(opt, val) {
            return true;
        }
        return !token.contains(&b' ');
    }
    if opt == b"-h" || token.starts_with(b"-h") {
        return true;
    }
    if !matches!(match_single_dash(opt), SingleDashMatch::None) {
        return true;
    }
    !token.contains(&b' ')
}

/// Parses `words` with the default 80-column width.
#[cfg(test)]
pub(crate) fn parse(words: &[OsString]) -> ParseOutcome {
    parse_with_columns(words, DEFAULT_COLUMNS)
}

// NCBI reference (598d8ae6): c++/src/app/blast/cleanup-blastdb-volumes.py:106-107,144-156
/// Parses `words` (arguments after the program name) as NCBI's `argparse` parser parses them,
/// wrapping usage and help output to `columns`.
pub(crate) fn parse_with_columns(words: &[OsString], columns: usize) -> ParseOutcome {
    let raw: Vec<Vec<u8>> = words
        .iter()
        .map(|w| os_str_bytes(w.as_os_str()).into_owned())
        .collect();

    let mut db: Option<Vec<u8>> = None;
    let mut dbtype: Option<DbType> = None;
    let mut dry_run = false;
    let mut verbose = 0usize;
    let mut extras: Vec<String> = Vec::new();
    let mut after_double_dash = false;
    let mut push_extra = |tok: &[u8]| extras.push(String::from_utf8_lossy(tok).into_owned());

    let mut i = 0usize;
    while i < raw.len() {
        let token = &raw[i];
        if after_double_dash {
            push_extra(token);
            i += 1;
            continue;
        }
        if token == b"--" {
            after_double_dash = true;
            push_extra(b"--");
            i += 1;
            continue;
        }
        if !is_optional_token(token) {
            push_extra(token);
            i += 1;
            continue;
        }

        let (opt_part, explicit_val) = split_eq(token);

        // Exact `-h` or double-dash `--help` (and prefixes `--h`, `--he`, `--hel`, plus `--=<val>`).
        if opt_part == b"-h"
            || (opt_part.starts_with(b"--") && is_double_dash_help(opt_part, explicit_val))
        {
            return help_or_ignored(columns, explicit_val);
        }
        if opt_part.starts_with(b"--") {
            push_extra(token);
            i += 1;
            continue;
        }

        // Single-dash multi-character options and abbreviations.
        let matched = match match_single_dash(opt_part) {
            SingleDashMatch::Exact(opt) | SingleDashMatch::Unique(opt) => opt,
            SingleDashMatch::Ambiguous(candidates) => {
                return err_outcome(
                    columns,
                    &format!(
                        "ambiguous option: {} could match {}",
                        String::from_utf8_lossy(token),
                        candidates.join(", ")
                    ),
                );
            }
            SingleDashMatch::None => {
                // Short-option bundling: `-h<rest>` errors when `<rest>` starts with `-`, and
                // otherwise dispatches `-h` immediately (`-help`, `-he`, `-hx`, `-h bar`).
                if let Some(rest) = token.strip_prefix(b"-h") {
                    if rest.starts_with(b"-") {
                        return ignored_arg(columns, "-h/--help", rest);
                    }
                    return ParseOutcome::Help(help_text_for_columns(columns));
                }
                push_extra(token);
                i += 1;
                continue;
            }
        };

        match matched {
            LongOpt::Version => {
                if let Some(val) = explicit_val {
                    return ignored_arg(columns, "-version", val);
                }
                return ParseOutcome::Version(version_text());
            }
            LongOpt::DryRun => {
                if let Some(val) = explicit_val {
                    return ignored_arg(columns, "-dry-run", val);
                }
                dry_run = true;
                i += 1;
            }
            LongOpt::Verbose => {
                if let Some(val) = explicit_val {
                    return ignored_arg(columns, "-verbose", val);
                }
                verbose += 1;
                i += 1;
            }
            LongOpt::Db => match take_value(columns, "-db", explicit_val, &raw, &mut i) {
                Ok(value) => db = Some(value),
                Err(outcome) => return outcome,
            },
            LongOpt::DbType => {
                let value = match take_value(columns, "-dbtype", explicit_val, &raw, &mut i) {
                    Ok(value) => value,
                    Err(outcome) => return outcome,
                };
                match value.as_slice() {
                    b"prot" => dbtype = Some(DbType::Prot),
                    b"nucl" => dbtype = Some(DbType::Nucl),
                    other => {
                        return err_outcome(
                            columns,
                            &format!(
                                "argument -dbtype: invalid choice: {} (choose from prot, nucl)",
                                py_repr(other)
                            ),
                        );
                    }
                }
            }
        }
    }

    let mut missing = Vec::new();
    if db.is_none() {
        missing.push("-db");
    }
    if dbtype.is_none() {
        missing.push("-dbtype");
    }
    if !missing.is_empty() {
        return err_outcome(
            columns,
            &format!(
                "the following arguments are required: {}",
                missing.join(", ")
            ),
        );
    }
    if !extras.is_empty() {
        return err_outcome(
            columns,
            &format!("unrecognized arguments: {}", extras.join(" ")),
        );
    }

    ParseOutcome::Run(ParsedArgs {
        db: db.expect("checked above"),
        dbtype: dbtype.expect("checked above"),
        dry_run,
        verbose,
    })
}

#[cfg(test)]
#[path = "args_tests.rs"]
mod tests;
