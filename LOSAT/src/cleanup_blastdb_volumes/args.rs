//! Command-line parser for `cleanup-blastdb-volumes.py`, reproducing Python `argparse`'s options,
//! prefix abbreviations, error messages, help text, and exit status `2`.
//!
//! NCBI reference (598d8ae6): c++/src/app/blast/cleanup-blastdb-volumes.py:40-41,144-156

use std::ffi::OsString;

use crate::utils::os_bytes::os_str_bytes;

/// Program name written in usage, help, version, and error lines.
pub(crate) const PROG: &str = "cleanup-blastdb-volumes.py";

// NCBI reference (598d8ae6): c++/src/app/blast/cleanup-blastdb-volumes.py:144-156
/// Usage line printed by `-h`, `--help`, and argument errors (observed from NCBI BLAST+ 2.17.0).
pub(crate) const USAGE: &str = "\
usage: cleanup-blastdb-volumes.py [-h] -db DB -dbtype {prot,nucl} [-dry-run]
                                  [-version] [-verbose]
";

// NCBI reference (598d8ae6): c++/src/app/blast/cleanup-blastdb-volumes.py:41,144-156
/// Full help text printed by `-h` and `--help` (observed from NCBI BLAST+ 2.17.0).
const HELP_BODY: &str = "\n\
Remove needless BLAST database volumes.

options:
  -h, --help           show this help message and exit
  -db DB               BLAST database name
  -dbtype {prot,nucl}  Molecule type
  -dry-run             Do not delete any files, just list them
  -version             show program's version number and exit
  -verbose             Increase output verbosity
";

/// Formats the complete `-h` / `--help` text.
pub(crate) fn help_text() -> String {
    format!("{USAGE}{HELP_BODY}")
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

/// Formats an argument error with the usage banner and program error line.
fn error_text(message: &str) -> String {
    format!("{USAGE}{PROG}: error: {message}\n")
}

/// Builds a `ParseOutcome::Error` from `message`.
fn err_outcome(message: &str) -> ParseOutcome {
    ParseOutcome::Error(error_text(message))
}

/// Builds a `ParseOutcome::Error` for a flag that does not accept an explicit `=value`.
fn ignored_arg(flag: &str, val: &[u8]) -> ParseOutcome {
    err_outcome(&format!(
        "argument {flag}: ignored explicit argument {}",
        py_repr(val)
    ))
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
    Err(err_outcome(&format!(
        "argument {flag}: expected one argument"
    )))
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
fn help_or_ignored(val: Option<&[u8]>) -> ParseOutcome {
    match val {
        Some(v) => ignored_arg("-h/--help", v),
        None => ParseOutcome::Help(help_text()),
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

// NCBI reference (598d8ae6): c++/src/app/blast/cleanup-blastdb-volumes.py:106-107,144-156
/// Parses `words` (arguments after the program name) as NCBI's `argparse` parser parses them.
pub(crate) fn parse(words: &[OsString]) -> ParseOutcome {
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

        // Double-dash option (`--help` and prefixes `--h`, `--he`, `--hel`, plus `--=<val>`).
        if opt_part.starts_with(b"--") {
            if is_double_dash_help(opt_part, explicit_val) {
                return help_or_ignored(explicit_val);
            }
            push_extra(token);
            i += 1;
            continue;
        }

        // Exact `-h` (handles `-h` and `-h=<val>`).
        if opt_part == b"-h" {
            return help_or_ignored(explicit_val);
        }

        // Single-dash multi-character options and abbreviations.
        let matched = match match_single_dash(opt_part) {
            SingleDashMatch::Exact(opt) | SingleDashMatch::Unique(opt) => opt,
            SingleDashMatch::Ambiguous(candidates) => {
                return err_outcome(&format!(
                    "ambiguous option: {} could match {}",
                    String::from_utf8_lossy(token),
                    candidates.join(", ")
                ));
            }
            SingleDashMatch::None => {
                // Short-option bundling: `-h<rest>` errors when `<rest>` starts with `-`, and
                // otherwise dispatches `-h` immediately (`-help`, `-he`, `-hx`, `-h bar`).
                if let Some(rest) = token.strip_prefix(b"-h") {
                    if rest.starts_with(b"-") {
                        return ignored_arg("-h/--help", rest);
                    }
                    return ParseOutcome::Help(help_text());
                }
                push_extra(token);
                i += 1;
                continue;
            }
        };

        match matched {
            LongOpt::Version => {
                if let Some(val) = explicit_val {
                    return ignored_arg("-version", val);
                }
                return ParseOutcome::Version(version_text());
            }
            LongOpt::DryRun => {
                if let Some(val) = explicit_val {
                    return ignored_arg("-dry-run", val);
                }
                dry_run = true;
                i += 1;
            }
            LongOpt::Verbose => {
                if let Some(val) = explicit_val {
                    return ignored_arg("-verbose", val);
                }
                verbose += 1;
                i += 1;
            }
            LongOpt::Db => match take_value("-db", explicit_val, &raw, &mut i) {
                Ok(value) => db = Some(value),
                Err(outcome) => return outcome,
            },
            LongOpt::DbType => {
                let value = match take_value("-dbtype", explicit_val, &raw, &mut i) {
                    Ok(value) => value,
                    Err(outcome) => return outcome,
                };
                match value.as_slice() {
                    b"prot" => dbtype = Some(DbType::Prot),
                    b"nucl" => dbtype = Some(DbType::Nucl),
                    other => {
                        return err_outcome(&format!(
                            "argument -dbtype: invalid choice: {} (choose from prot, nucl)",
                            py_repr(other)
                        ));
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
        return err_outcome(&format!(
            "the following arguments are required: {}",
            missing.join(", ")
        ));
    }
    if !extras.is_empty() {
        return err_outcome(&format!("unrecognized arguments: {}", extras.join(" ")));
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
