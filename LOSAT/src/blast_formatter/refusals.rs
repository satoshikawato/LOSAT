//! What LOSAT's blast_formatter refuses, by name: the one list of its refusals. Each refusal
//! stops the run with LOSAT's message instead of writing something else; README.md's "What LOSAT
//! refuses" names each one (`cli::readme_refusals_tests`). No NCBI counterpart: NCBI's
//! blast_formatter writes every output format from binary, text and XML archives, fetches the
//! results of `-rid` from NCBI's servers and asks them for the sequences of a database it does
//! not find (LOSAT's search programs refuse `-remote` alike, line RG DECISIONS 7 and F3-5).

use crate::ncbi_app::Diag;

/// The output formats that blast_formatter writes; it refuses the others by name
/// (`Refusal::OutputFormat`).
pub(crate) const WRITTEN_FORMATS: [i32; 3] = [8, 9, 12];

/// The output formats of NCBI's blast_formatter (0 to 18 and 20; `-help`, probe P2.201).
const NCBI_FORMATS: [i32; 20] = [
    0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 20,
];

/// A run that LOSAT's blast_formatter refuses.
#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) enum Refusal {
    /// `-rid`: the results of a search that NCBI's servers ran, which NCBI's blast_formatter
    /// fetches from them; refused where NCBI asks for them, as the search programs refuse
    /// `-remote` (`cli::named_option_refusal`).
    Rid,
    /// An output format of NCBI's blast_formatter that LOSAT's does not write
    /// (`WRITTEN_FORMATS`).
    OutputFormat(i32),
    /// An archive in ASN.1 binary (`-outfmt 11` with `ARCHIVE_FORMAT=bin`).
    BinaryArchive,
    /// An archive in XML (`-outfmt 11` with `ARCHIVE_FORMAT=xml`).
    XmlArchive,
    /// The database of an archive (its name) that is not found as NCBI finds a database: NCBI's
    /// blast_formatter then asks NCBI's servers for its sequences.
    Database(String),
    /// Every run that reaches the archive, in a build without the archive's reader (removed when
    /// the reader joins: line RG item F3, step C).
    ArchiveReader,
}

impl Refusal {
    /// The refusal's message (after `Error: `; the clap form `error: ...` for `-rid`).
    pub(crate) fn reason(&self) -> String {
        match self {
            Refusal::Rid => crate::cli::command_line_message(
                &crate::cli::named_option_refusal("rid", super::PROGRAM),
                false,
            ),
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
                "the database '{name}' of the BLAST archive is not found here; NCBI BLAST+ then \
                 asks NCBI's servers for its sequences, which is not supported by LOSAT's {}",
                super::PROGRAM
            ),
            Refusal::ArchiveReader => format!(
                "reading a BLAST archive is not supported by this build of LOSAT's {} yet",
                super::PROGRAM
            ),
        }
    }

    /// Stops the run with the refusal: `-rid` as the search programs' `-remote` (clap's error on
    /// the standard error stream, exit status 2), the others as LOSAT's refusals of NCBI's
    /// applications (`Diag::refuse`, exit status 1).
    pub(crate) fn refuse(&self, diag: &mut Diag) -> i32 {
        match self {
            Refusal::Rid => {
                let error = crate::cli::named_option_refusal("rid", super::PROGRAM);
                diag.cerr(crate::cli::command_line_message(&error, false).as_bytes());
                error.exit_code()
            }
            _ => diag.refuse(&self.reason()),
        }
    }

    /// The words of README.md's paragraph on blast_formatter that name the refusal; `None` for
    /// `ArchiveReader`, which no release has.
    pub(crate) fn readme_words(&self) -> Option<String> {
        match self {
            Refusal::Rid => Some("`-rid`".to_string()),
            Refusal::OutputFormat(number) => Some(format!("`{number}`")),
            Refusal::BinaryArchive => Some("ASN.1 binary".to_string()),
            Refusal::XmlArchive => Some("XML".to_string()),
            Refusal::Database(_) => Some("database".to_string()),
            Refusal::ArchiveReader => None,
        }
    }
}

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

/// One refusal of each kind that a release has, the output formats one by one: what README.md
/// documents (`cli::readme_refusals_tests`).
pub(crate) fn documented() -> Vec<Refusal> {
    let mut refusals = vec![Refusal::Rid];
    refusals.extend(NCBI_FORMATS.into_iter().filter_map(output_format));
    refusals.extend([
        Refusal::BinaryArchive,
        Refusal::XmlArchive,
        Refusal::Database(String::new()),
    ]);
    refusals
}

#[cfg(test)]
#[path = "refusals_tests.rs"]
mod refusals_tests;
