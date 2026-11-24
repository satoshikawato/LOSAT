//! What LOSAT's blast_formatter refuses, by name: the one list of its refusals. Each refusal
//! stops the run with LOSAT's message instead of writing something else; the manual's "What
//! LOSAT refuses" names each one (`cli::manual_refusals_tests`). No NCBI counterpart: NCBI's
//! blast_formatter writes every output format from binary, text and XML archives, fetches the
//! results of `-rid` from NCBI's servers and asks them for the sequences of a database it does
//! not open (LOSAT's search programs refuse `-remote` alike, line RG DECISIONS 7 and F3-5), and
//! on some archives that its own programs do not write it reads memory the archive does not
//! decide or aborts; on a stored message of severity Fatal it stops with a stack trace of its own
//! binary.

use crate::ncbi_app::Diag;

/// The output formats that blast_formatter writes; it refuses the others by name
/// (`Refusal::OutputFormat`).
pub(crate) const WRITTEN_FORMATS: [i32; 7] = [6, 7, 8, 9, 10, 12, 20];

/// The output formats of NCBI's blast_formatter (0 to 18 and 20; `-help`, probe P2.201): the
/// tests' list of the refusals (`documented`).
#[cfg(test)]
const NCBI_FORMATS: [i32; 20] = [
    0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 20,
];

/// A run that LOSAT's blast_formatter refuses.
#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) enum Refusal {
    /// `-rid`: the results of a search that NCBI's servers ran, which NCBI's blast_formatter
    /// fetches from them; refused where NCBI asks for them (the row of the one list of refused
    /// arguments, `ncbi_app::args::refusal`, read by `cli::run_refusal`).
    Rid,
    /// An output format of NCBI's blast_formatter that LOSAT's does not write
    /// (`WRITTEN_FORMATS`).
    OutputFormat(i32),
    /// An archive in ASN.1 binary (`-outfmt 11` with `ARCHIVE_FORMAT=bin`).
    BinaryArchive,
    /// An archive in XML (`-outfmt 11` with `ARCHIVE_FORMAT=xml`).
    XmlArchive,
    /// The database of an archive (its name) that NCBI's database reader does not open with the
    /// archive's taxid list: NCBI's blast_formatter then asks NCBI's servers for its sequences.
    Database(String),
    /// An archive on which NCBI's code reads memory that the archive does not decide, follows a
    /// null pointer or aborts (`report::archive::read::Stop`, the format guess's undecided
    /// answer): what NCBI meets.
    Undecided(&'static str),
    /// A form of ASN.1 text that LOSAT's reader does not read (its text).
    UnreadText(String),
    /// A stored message of the severity `Fatal`, whose post ends NCBI's program (`abort`).
    FatalMessage,
    /// `PRE_FETCH_SEQS_LIMIT=2147483647` with an output format that LOSAT writes (8, 9, 12):
    /// NCBI's blast_formatter then loads each subject of the results through its data loaders
    /// (its BLAST database loader, or NCBI's servers where the database or the subjects lack
    /// one), which may post their own texts (decision C3-1 of line RG).
    Prefetch,
    /// An archive of the mapper (its task, `mapr2g`), whose name NCBI's formatter maps back to
    /// the mapper; NCBI's programs write no such archive (PHI-BLAST's archives end with NCBI's
    /// engine error instead, `search_program::blast_type`).
    Program(&'static str),
    /// A value in an archive's results that LOSAT's encoder of the output format does not write
    /// (`report::serial::SerialError`), what it is and the format: a NaN or an infinite real
    /// number in `-outfmt 12`, at which NCBI stops with an error of its JSON stream (agent C6 of
    /// line RG item F3, 2026-11-23, `remaining-gaps-20261112/F/f3/tmp/c6/pd`; in 9 NCBI writes
    /// its one special octet, as LOSAT's ASN.1 binary encoder does), and the values the encoders
    /// leave unwritten (in 12 a UTF8String, a NULL value and an OCTET STRING; a member read into
    /// a delay buffer). The bytes of the Seq-annot that went to the output before it stay
    /// written.
    Unwritable(&'static str, i32),
    /// A form of a BLAST archive that the tabular formats (6, 7, 10 and 20) do not write yet: the
    /// translated programs, queries given as locations, and the forms whose subjects or rows
    /// NCBI's scope may find otherwise than LOSAT (decision F3-13 of line RG). What it is.
    Tabular(String),
}

impl Refusal {
    /// The refusal's message (after `Error: `).
    pub(crate) fn reason(&self) -> String {
        match self {
            Refusal::Rid => crate::cli::run_refusal(super::PROGRAM, "rid"),
            Refusal::OutputFormat(number) => format!(
                "output format {number} is not supported by LOSAT's {} (it writes {})",
                super::PROGRAM,
                written_formats()
            ),
            Refusal::BinaryArchive => format!(
                "a BLAST archive in ASN.1 binary is not supported by LOSAT's {}: it reads the \
                 archives in ASN.1 text that -outfmt 11 writes by default",
                super::PROGRAM
            ),
            Refusal::XmlArchive => format!(
                "a BLAST archive in XML is not supported by LOSAT's {}: it reads the archives in \
                 ASN.1 text that -outfmt 11 writes by default",
                super::PROGRAM
            ),
            Refusal::Database(name) => format!(
                "the database '{name}' of the BLAST archive cannot be opened here; NCBI BLAST+ \
                 then asks NCBI's servers for its sequences, which is not supported by LOSAT's {}",
                super::PROGRAM
            ),
            Refusal::Undecided(what) => format!(
                "a BLAST archive on which NCBI BLAST+ has no defined result ({what}) is not \
                 supported by LOSAT's {}",
                super::PROGRAM
            ),
            Refusal::UnreadText(text) => format!(
                "a BLAST archive in a form of ASN.1 text that LOSAT does not read is not \
                 supported by LOSAT's {} ({text})",
                super::PROGRAM
            ),
            Refusal::FatalMessage => format!(
                "a BLAST archive with a stored message of severity Fatal, on which NCBI BLAST+ \
                 stops with exit status 255 and a stack trace of its own binary, is not supported \
                 by LOSAT's {}",
                super::PROGRAM
            ),
            Refusal::Prefetch => format!(
                "PRE_FETCH_SEQS_LIMIT=2147483647, with which NCBI BLAST+ loads the subject \
                 sequences of output formats 8, 9 and 12 ahead (from NCBI's servers where the \
                 database lacks one), is not supported by LOSAT's {}",
                super::PROGRAM
            ),
            Refusal::Program(task) => format!(
                "a BLAST archive of {task} is not supported by LOSAT's {}",
                super::PROGRAM
            ),
            Refusal::Unwritable(what, format) => format!(
                "output format {format} of a BLAST archive with {what} in its results is not \
                 supported by LOSAT's {}",
                super::PROGRAM
            ),
            Refusal::Tabular(what) => format!(
                "the tabular output formats (6, 7, 10 and 20) of a BLAST archive with {what} are \
                 not supported by LOSAT's {}",
                super::PROGRAM
            ),
        }
    }

    /// Stops the run with the refusal, as LOSAT's refusals of NCBI's applications
    /// (`Diag::refuse`, exit status 1).
    pub(crate) fn refuse(&self, diag: &mut Diag) -> i32 {
        diag.refuse(&self.reason())
    }

    /// The words of the manual's section on blast_formatter (`docs/manual/src/refusals.md`) that
    /// name the refusal (the tests' check of the manual).
    #[cfg(test)]
    pub(crate) fn manual_words(&self) -> String {
        match self {
            Refusal::Rid => "`-rid`".to_string(),
            Refusal::OutputFormat(number) => format!("`{number}`"),
            Refusal::BinaryArchive => "ASN.1 binary".to_string(),
            Refusal::XmlArchive => "XML".to_string(),
            Refusal::Database(_) => "cannot be opened".to_string(),
            Refusal::Undecided(_) => "no defined result".to_string(),
            Refusal::UnreadText(_) => "form of ASN.1 text".to_string(),
            Refusal::FatalMessage => "severity Fatal".to_string(),
            Refusal::Prefetch => "`PRE_FETCH_SEQS_LIMIT=2147483647`".to_string(),
            Refusal::Program(_) => "mapper".to_string(),
            Refusal::Unwritable(..) => "an archive's results".to_string(),
            Refusal::Tabular(_) => "tabular output formats".to_string(),
        }
    }
}

impl std::fmt::Display for Refusal {
    /// The refusal's message (`Refusal::reason`): what `super::catch_all` writes after `Error: `.
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(&self.reason())
    }
}

impl std::error::Error for Refusal {}

/// The written formats as a list in words: `8, 9 and 12`.
fn written_formats() -> String {
    let names: Vec<String> = WRITTEN_FORMATS.iter().map(i32::to_string).collect();
    match names.split_last() {
        Some((last, [])) => last.clone(),
        Some((last, rest)) => format!("{} and {last}", rest.join(", ")),
        None => String::new(),
    }
}

/// The refusal of the output format `number`, which `extract_formatting` read: `None` for a
/// format that blast_formatter writes.
pub(crate) fn output_format(number: i32) -> Option<Refusal> {
    (!WRITTEN_FORMATS.contains(&number)).then_some(Refusal::OutputFormat(number))
}

/// One refusal of each kind, the output formats one by one: what the manual documents
/// (`cli::manual_refusals_tests`).
#[cfg(test)]
pub(crate) fn documented() -> Vec<Refusal> {
    let mut refusals = vec![Refusal::Rid];
    refusals.extend(NCBI_FORMATS.into_iter().filter_map(output_format));
    refusals.extend([
        Refusal::BinaryArchive,
        Refusal::XmlArchive,
        Refusal::Database(String::new()),
        Refusal::Undecided(""),
        Refusal::UnreadText(String::new()),
        Refusal::FatalMessage,
        Refusal::Prefetch,
        Refusal::Program("mapr2g"),
        Refusal::Unwritable("", 12),
        Refusal::Tabular(String::new()),
    ]);
    refusals
}

#[cfg(test)]
#[path = "refusals_tests.rs"]
mod refusals_tests;
