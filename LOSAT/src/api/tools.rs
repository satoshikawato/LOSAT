//! The tools of the command line besides the searches: `makeblastdb` builds a BLAST database,
//! `blastdbcmd` reads one, `dustmasker` masks the low-complexity regions of nucleotide
//! sequences and `segmasker` those of protein sequences, `windowmasker` masks the over-represented
//! regions of nucleotide sequences by a counts file, `convert2blastmask` turns the lower case
//! letters of FASTA input into masking information, `makembindex` writes the megablast index
//! of a BLAST database or of FASTA input, and `blast_formatter` writes the reports of a BLAST
//! archive (`-outfmt 11`) again. They take the command line's words,
//! read the caller's reader where the command line reads its standard input, and write what it
//! writes, to the caller's writers; none touches the process's standard streams.
//!
//! NCBI reference (598d8ae6): c++/src/app/blastdb/makeblastdb.cpp:1360-1363
//! NCBI reference (598d8ae6): c++/src/app/blastdb/blastdbcmd.cpp:1295-1298
//! NCBI reference (598d8ae6): c++/src/app/dustmasker/main.cpp:39-42
//! NCBI reference (598d8ae6): c++/src/app/segmasker/segmasker.cpp:265-272
//! NCBI reference (598d8ae6): c++/src/app/winmasker/main.cpp:39-42
//! NCBI reference (598d8ae6): c++/src/app/blastdb/convert2blastmask.cpp:259-263
//! NCBI reference (598d8ae6): c++/src/algo/blast/dbindex/makeindex/main.cpp:51-54
//! NCBI reference (598d8ae6): c++/src/app/blast/blast_formatter.cpp:588-592
#![warn(missing_docs, unused)]

use std::ffi::OsString;
use std::io::{self, Read, Write};

use super::error::Error;
use crate::blastinput::ncbi_environment::{preparse, Preparsed};
use crate::cli::Commands;
use crate::writedb::makeblastdb;
use crate::writedb::time::LocalTime;
use crate::writedb::Platform;

/// The caller's stream of a tool's diagnostics, which keeps what the tool writes there: a tool
/// writes its failure there, and that text without its last newlines is the error's message, also
/// where the caller's writer fails.
struct Recorded<'a> {
    writer: &'a mut dyn Write,
    bytes: Vec<u8>,
}

impl Write for Recorded<'_> {
    /// Keeps the bytes that the caller's writer takes, and the whole `buffer` where it takes none
    /// (an error, or nothing written): the command line's standard error would hold them, and
    /// `write_all` gives them up there. An interrupted write keeps nothing, as `write_all` writes
    /// the buffer again.
    fn write(&mut self, buffer: &[u8]) -> io::Result<usize> {
        let result = self.writer.write(buffer);
        let kept = match &result {
            Ok(written) if *written > 0 => &buffer[..*written],
            Err(error) if error.kind() == io::ErrorKind::Interrupted => &[],
            _ => buffer,
        };
        self.bytes.extend_from_slice(kept);
        result
    }

    fn flush(&mut self) -> io::Result<()> {
        self.writer.flush()
    }
}

/// The result of a tool that ended with `status` after writing `diagnostics`: their bytes, as the
/// command line writes them to its standard error (`Error::written`).
fn outcome(status: i32, diagnostics: Vec<u8>) -> Result<(), Error> {
    if status == 0 {
        return Ok(());
    }
    Err(Error::written(diagnostics))
}

/// LOSAT's refusal `reason` in a tool's run, written to its diagnostics as the command line's
/// `main` writes it (`Error: ` and the reason; exit status 1 there), and the error of what the
/// run wrote (`outcome`), as the tools of `ncbi_app` refuse.
fn refused(mut recorded: Recorded<'_>, reason: &[u8]) -> Error {
    let _ = recorded.write_all(&[b"Error: ", reason, b"\n"].concat());
    match outcome(1, recorded.bytes) {
        Err(error) => error,
        Ok(()) => unreachable!("a status of 1 is a failure"),
    }
}

/// LOSAT's refusal `reason` of an argument of a tool's argv (makeblastdb's, after NCBI's reader:
/// `cli::reader::refusal_after_read`), written to its diagnostics as the command line writes it
/// (`Error: ` and the reason; exit status 1 there); the error is the refusal's text, of the argv's
/// kind (`Error::Argv`), as the searches' refusals of an argument.
fn refused_argument(mut recorded: Recorded<'_>, reason: String) -> Error {
    let _ = recorded.write_all(format!("Error: {reason}\n").as_bytes());
    Error::Argv(reason.into_bytes())
}

/// Runs `blastdbcmd` with the words after the program name, as on the command line. `input` is what
/// the command line reads as standard input (`-entry_batch -`, `-ipg_batch -`, `-taxidlist -`), so
/// a batch of ids can be given from memory. What the command line writes to standard output goes to
/// `output` (`-out` is a file, as there), what it writes to standard error to `diagnostics`. A
/// failure is an error whose message is that text without its last newlines (an argv that NCBI's
/// argument reader refuses too). The call runs NCBI's application start, as the command line does
/// (`-logfile`, `-conffile`, `-dryrun`, the versions); `-xmlhelp`, whose XML only the command line
/// writes, is refused.
pub fn run_blastdbcmd(
    argv: &[&str],
    input: &mut dyn Read,
    output: &mut dyn Write,
    diagnostics: &mut dyn Write,
) -> Result<(), Error> {
    let words: Vec<OsString> = argv.iter().map(OsString::from).collect();
    let mut recorded = Recorded {
        writer: diagnostics,
        bytes: Vec::new(),
    };
    let status = crate::blastdbcmd::run_with(&words, input, output, &mut recorded);
    outcome(status, recorded.bytes)
}

/// Runs `dustmasker` with the words after the program name, as on the command line. `input` is what
/// the command line reads as standard input (`-in -`, and without `-in`), so sequences in memory
/// can be masked. What the command line writes to standard output goes to `output` (`-out` is a
/// file, as there), what it writes to standard error to `diagnostics`. A failure is an error whose
/// message is that text without its last newlines (an argv that NCBI's argument reader refuses too:
/// its usage and report, as the command line prints them). The call runs NCBI's application start,
/// as the command line does (`-logfile`, `-conffile`, `-dryrun`, the versions); `-xmlhelp`, whose
/// XML only the command line writes, is refused.
pub fn run_dustmasker(
    argv: &[&str],
    input: &mut dyn Read,
    output: &mut dyn Write,
    diagnostics: &mut dyn Write,
) -> Result<(), Error> {
    let words: Vec<OsString> = argv.iter().map(OsString::from).collect();
    let mut recorded = Recorded {
        writer: diagnostics,
        bytes: Vec::new(),
    };
    let status = crate::dustmasker::run_with(&words, input, output, &mut recorded);
    outcome(status, recorded.bytes)
}

/// Runs `segmasker` with the words after the program name, as on the command line: as
/// [`run_dustmasker`], for protein sequences. `input` is what the command line reads as standard
/// input (`-in -`, and without `-in`), what it writes to standard output goes to `output` (`-out`
/// is a file, as there), what it writes to standard error to `diagnostics`, and a failure is an
/// error whose message is that text without its last newlines. The call runs NCBI's application
/// start, as the command line does (`-logfile`, `-conffile`, `-dryrun`, the versions); `-xmlhelp`,
/// whose XML only the command line writes, is refused.
pub fn run_segmasker(
    argv: &[&str],
    input: &mut dyn Read,
    output: &mut dyn Write,
    diagnostics: &mut dyn Write,
) -> Result<(), Error> {
    let words: Vec<OsString> = argv.iter().map(OsString::from).collect();
    let mut recorded = Recorded {
        writer: diagnostics,
        bytes: Vec::new(),
    };
    let status = crate::segmasker::run_with(&words, input, output, &mut recorded);
    outcome(status, recorded.bytes)
}

/// Runs `windowmasker` with the words after the program name, as on the command line: as
/// [`run_dustmasker`], in its three modes: making unit counts (`-mk_counts`), converting them to
/// another format (`-convert`), and masking with them (`-ustat`). `input` is what the command line
/// reads as standard input (`-in -`, and without `-in`), what it writes to standard output goes to
/// `output` (`-out` is a file, as there), what it writes to standard error to `diagnostics`, and a
/// failure is an error whose message is that text without its last newlines. The call runs NCBI's application
/// start, as the command line does (`-logfile`, `-conffile`, `-dryrun`, the versions); `-xmlhelp`,
/// whose XML only the command line writes, is refused.
pub fn run_windowmasker(
    argv: &[&str],
    input: &mut dyn Read,
    output: &mut dyn Write,
    diagnostics: &mut dyn Write,
) -> Result<(), Error> {
    let words: Vec<OsString> = argv.iter().map(OsString::from).collect();
    let mut recorded = Recorded {
        writer: diagnostics,
        bytes: Vec::new(),
    };
    let status = crate::windowmasker::run_with(&words, input, output, &mut recorded);
    outcome(status, recorded.bytes)
}

/// Runs `convert2blastmask` with the words after the program name, as on the command line: as
/// [`run_dustmasker`], for the lower case letters of FASTA input. `input` is what the command
/// line reads as standard input (`-in -`, and without `-in`), what it writes to standard output
/// goes to `output` (`-out` is a file, as there), what it writes to standard error to
/// `diagnostics`, and a failure is an error whose message is that text without its last newlines.
pub fn run_convert2blastmask(
    argv: &[&str],
    input: &mut dyn Read,
    output: &mut dyn Write,
    diagnostics: &mut dyn Write,
) -> Result<(), Error> {
    let words: Vec<OsString> = argv.iter().map(OsString::from).collect();
    let mut recorded = Recorded {
        writer: diagnostics,
        bytes: Vec::new(),
    };
    let status = crate::convert2blastmask::run_with(&words, input, output, &mut recorded);
    outcome(status, recorded.bytes)
}

/// Runs `makembindex` with the words after the program name, as on the command line: as
/// [`run_dustmasker`], for the index of a BLAST database (`-input`) or of FASTA input
/// (`-old_style_index true -iformat fasta`). `input` is what the command line reads as standard
/// input (FASTA input without `-input`), what it writes to standard output (`-show_filters`) goes
/// to `output`, what it writes to standard error (the `creating ...` lines too) to `diagnostics`,
/// and a failure is an error whose message is that text without its last newlines. The index files
/// are written where the command line writes them.
pub fn run_makembindex(
    argv: &[&str],
    input: &mut dyn Read,
    output: &mut dyn Write,
    diagnostics: &mut dyn Write,
) -> Result<(), Error> {
    let words: Vec<OsString> = argv.iter().map(OsString::from).collect();
    let mut recorded = Recorded {
        writer: diagnostics,
        bytes: Vec::new(),
    };
    let status = crate::makembindex::run_with(&words, input, output, &mut recorded);
    outcome(status, recorded.bytes)
}

/// Runs `blast_formatter` with the words after the program name, as on the command line: as
/// [`run_dustmasker`], for the reports of a BLAST archive (`-archive`, the file of `-outfmt 11`).
/// `input` is what the command line reads as standard input (`-archive -`), what it writes to
/// standard output goes to `output` (`-out` is a file, as there), what it writes to standard error
/// to `diagnostics`, and a failure is an error whose message is that text.
pub fn run_blast_formatter(
    argv: &[&str],
    input: &mut dyn Read,
    output: &mut dyn Write,
    diagnostics: &mut dyn Write,
) -> Result<(), Error> {
    let words: Vec<OsString> = argv.iter().map(OsString::from).collect();
    let mut recorded = Recorded {
        writer: diagnostics,
        bytes: Vec::new(),
    };
    let status = crate::blast_formatter::run_with(&words, input, output, &mut recorded);
    outcome(status, recorded.bytes)
}

// NCBI reference (598d8ae6): c++/src/corelib/ncbiapp.cpp:1031-1041
// NCBI reference (598d8ae6): c++/src/app/blastdb/makeblastdb.cpp:1134-1136
/// Runs `makeblastdb` with the words after the program name, as on the command line, which
/// writes the database's files as there. `input` is what the command line reads as standard input
/// (`-in -`, and without `-in`), so a database can be built from FASTA in memory. What the command
/// line writes to standard output goes to `output`, what it writes to standard error to
/// `diagnostics`. A failure is an error whose message is that text without its last newlines; an
/// argv that NCBI's argument reader refuses is `Error::Argv`, with the reader's report in the
/// diagnostics too (with `-logfile`, the report's posts go to its file and the diagnostics say
/// where they are, as on the command line), and so is LOSAT's refusal of an NCBI option that it
/// does not run. A help word (`-h`, `-help`, `-help-full`) writes the help to `output` and returns
/// `Ok(())`, as the command line prints it and exits with status 0, and as the library's other
/// tools do. The call runs NCBI's application start, as the command line does:
/// the file of `-logfile` gets its messages, the file of `-conffile` is loaded, the versions are
/// printed to `output`, and `-dryrun` ends the call after the arguments and the settings.
pub fn run_makeblastdb(
    argv: &[&str],
    input: &mut dyn Read,
    output: &mut dyn Write,
    diagnostics: &mut dyn Write,
) -> Result<(), Error> {
    let mut recorded = Recorded {
        writer: diagnostics,
        bytes: Vec::new(),
    };
    // NCBI's application start (`cli::main_entry`): the usage report's messages, the `-logfile`
    // file, the versions (`ncbi_version::print`, line RG item E3), then the configuration load.
    let mut opened = match super::start::open(
        "makeblastdb",
        argv,
        crate::report::log_inputs::Supplied::default(),
        &mut recorded,
    ) {
        Ok(opened) => opened,
        Err(error) => return Err(refused(recorded, error.message())),
    };
    let words = || ["losat", "makeblastdb"].iter().chain(argv).copied();
    if let Some(text) = crate::cli::ncbi_version_request(words()) {
        let _ = output.write_all(text.as_bytes());
        return Ok(());
    }
    match opened
        .start
        .write_config_load(&mut opened.target, &mut recorded)
    {
        Ok(None) => {}
        Ok(Some(status)) => return outcome(status, recorded.bytes),
        Err(reason) => return Err(refused(recorded, reason.as_bytes())),
    }
    // NCBI's reader of the words that the start leaves to it (`-dryrun` taken), as on the command
    // line (`cli::main_entry`): the help on `output`, an argument error's usage on the
    // diagnostics and its posts where the diagnostics go (the file of `-logfile`); `-xmlhelp`,
    // whose XML only the command line writes, is refused.
    let tool = crate::cli::reader::tool("makeblastdb").expect("makeblastdb is on the reader");
    let (words, dryrun) = match preparse(argv) {
        Preparsed::Words { words, dryrun, .. } => (words, dryrun),
        Preparsed::Version(_) => unreachable!("the versions are printed above"),
    };
    let words: Vec<OsString> = words.iter().map(OsString::from).collect();
    let mut diag = crate::ncbi_app::Diag::new(&opened.target, &mut recorded, tool);
    let read = crate::ncbi_app::read(
        tool,
        &words,
        output,
        &mut diag,
        None,
        &crate::report::run_files::OpenFiles::disk(),
    );
    drop(diag);
    let args = match read {
        Ok(args) => args,
        Err(crate::ncbi_app::ReadEnd::Exit(0)) => return Ok(()),
        Err(crate::ncbi_app::ReadEnd::Exit(_)) => return Err(Error::written_argv(recorded.bytes)),
        Err(crate::ncbi_app::ReadEnd::Refused(reason)) => {
            return Err(refused_argument(recorded, reason))
        }
    };
    // `-dryrun` ends here, after the settings that LOSAT refuses, as on the command line
    // (`cli::main_entry`).
    if dryrun {
        return match opened.start.settings() {
            Ok(_) => Ok(()),
            Err(reason) => Err(refused(recorded, reason.as_bytes())),
        };
    }
    let converted = match crate::cli::reader::refusal_after_read(tool, &args) {
        Some(reason) => Err(reason),
        None => crate::cli::reader::options(tool, &args),
    };
    let args = match converted {
        Ok(Commands::Makedb(args)) => args,
        Ok(_) => unreachable!("makeblastdb's arguments are makeblastdb's"),
        Err(reason) => return Err(refused_argument(recorded, reason)),
    };
    // NCBI's application layer reads the environment and the registry files before `Run`, and
    // LOSAT rejects the settings that change the output (`cli/main_entry.rs`).
    let settings = match opened.start.settings() {
        Ok(settings) => settings,
        Err(reason) => return Err(refused(recorded, reason.as_bytes())),
    };
    let checked = makeblastdb::check_registry_settings(
        &args,
        settings
            .max_seqid_length
            .as_ref()
            .map(Option::as_deref)
            .map_err(String::as_str),
    )
    .and_then(|()| makeblastdb::check_environment(&args, &|name: &str| std::env::var_os(name)));
    if let Err(reason) = checked {
        return Err(refused(recorded, reason.as_bytes()));
    }
    let now = LocalTime::now();
    let status = makeblastdb::run_with(
        &args,
        &settings,
        now.as_ref(),
        Platform::host(),
        input,
        output,
        &mut recorded,
    );
    outcome(status, recorded.bytes)
}
