//! The tools of the command line besides the searches: `makeblastdb` builds a BLAST database,
//! `blastdbcmd` reads one, `dustmasker` masks the low-complexity regions of nucleotide
//! sequences and `segmasker` those of protein sequences, and `convert2blastmask` turns the lower
//! case letters of FASTA input into masking information. They take the command line's words,
//! read the caller's reader where the command line reads its standard input, and write what it
//! writes, to the caller's writers; none touches the process's standard streams.
//!
//! NCBI reference (598d8ae6): c++/src/app/blastdb/makeblastdb.cpp:1360-1363
//! ```c++
//! int main(int argc, const char* argv[] /*, const char* envp[]*/)
//! {
//!     return CMakeBlastDBApp().AppMain(argc, argv);
//! }
//! ```
//! NCBI reference (598d8ae6): c++/src/app/blastdb/blastdbcmd.cpp:1295-1298
//! ```c++
//! int main(int argc, const char* argv[] /*, const char* envp[]*/)
//! {
//!     return CBlastDBCmdApp().AppMain(argc, argv);
//! }
//! ```
//! NCBI reference (598d8ae6): c++/src/app/dustmasker/main.cpp:39-42
//! ```c++
//! int main( int argc, char * argv[] )
//! {
//!     return CDustMaskApplication().AppMain(argc, argv);
//! }
//! ```
//! NCBI reference (598d8ae6): c++/src/app/segmasker/segmasker.cpp:265-272
//! ```c++
//! int main(int argc, const char* argv[])
//! {
//!     // Execute main application function
//!     return SegMaskerApplication().AppMain(argc, argv);
//! }
//! ```
//! NCBI reference (598d8ae6): c++/src/app/blastdb/convert2blastmask.cpp:259-263
//! ```c++
//! int main(int argc, const char* argv[])
//! {
//!     // Execute main application function
//!     return CConvert2BlastMaskApplication().AppMain(argc, argv);
//! }
//! ```
#![warn(missing_docs, unused)]

use std::ffi::OsString;
use std::io::{self, Read, Write};

use super::error::Error;
use crate::blastinput::ncbi_environment::check_ncbi_application_settings;
use crate::cli::{Cli, Commands};
use crate::writedb::makeblastdb;
use crate::writedb::time::LocalTime;
use crate::writedb::Platform;

/// The caller's stream of a tool's diagnostics, which keeps what it passes on: a tool writes its
/// failure there, and that text is the error's message.
struct Recorded<'a> {
    writer: &'a mut dyn Write,
    bytes: Vec<u8>,
}

impl Write for Recorded<'_> {
    fn write(&mut self, buffer: &[u8]) -> io::Result<usize> {
        let written = self.writer.write(buffer)?;
        self.bytes.extend_from_slice(&buffer[..written]);
        Ok(written)
    }

    fn flush(&mut self) -> io::Result<()> {
        self.writer.flush()
    }
}

/// The result of a tool that ended with `status` after writing `diagnostics`.
fn outcome(status: i32, diagnostics: Vec<u8>) -> Result<(), Error> {
    if status == 0 {
        return Ok(());
    }
    Err(Error::Search(
        String::from_utf8_lossy(&diagnostics)
            .trim_end_matches('\n')
            .to_string(),
    ))
}

/// Runs `blastdbcmd` with the words after the program name, as on the command line. `input` is
/// what the command line reads as standard input (`-entry_batch -`, `-ipg_batch -`,
/// `-taxidlist -`), so a batch of ids can be given from memory. What the command line writes to
/// standard output goes to `output` (`-out` is a file, as there), what it writes to standard
/// error to `diagnostics`. A failure is an error whose message is that text (an argv that NCBI's
/// argument reader refuses too). `-conffile` and `-dryrun`, which the command line's NCBI
/// application start reads, are refused.
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

/// Runs `dustmasker` with the words after the program name, as on the command line. `input` is
/// what the command line reads as standard input (`-in -`, and without `-in`), so sequences in
/// memory can be masked. What the command line writes to standard output goes to `output`
/// (`-out` is a file, as there), what it writes to standard error to `diagnostics`. A failure
/// is an error whose message is that text (an argv that NCBI's argument reader refuses too: its
/// usage and report, as the command line prints them). `-conffile` and `-dryrun`, which the
/// command line's NCBI application start reads, are refused.
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
/// error whose message is that text. `-conffile` and `-dryrun`, which the command line's NCBI
/// application start reads, are refused.
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

/// Runs `convert2blastmask` with the words after the program name, as on the command line: as
/// [`run_dustmasker`], for the lower case letters of FASTA input. `input` is what the command
/// line reads as standard input (`-in -`, and without `-in`), what it writes to standard output
/// goes to `output` (`-out` is a file, as there), what it writes to standard error to
/// `diagnostics`, and a failure is an error whose message is that text.
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

// NCBI reference (598d8ae6): c++/src/corelib/ncbiapp.cpp:1031-1041
// ```c++
//     // Setup some debugging features from environment variables.
//     if ( !m_Environ->Get(DIAG_TRACE).empty() ) {
//         SetDiagTrace(eDT_Enable, eDT_Enable);
//     }
//     string post_level = m_Environ->Get(DIAG_POST_LEVEL);
//     if ( !post_level.empty() ) {
//         EDiagSev sev;
//         if (CNcbiDiag::StrToSeverityLevel(post_level.c_str(), sev)) {
//             SetDiagFixedPostLevel(sev);
//         }
//     }
// ```
// NCBI reference (598d8ae6): c++/src/app/blastdb/makeblastdb.cpp:1134-1136
// ```c++
//     const string dont_scan_bioseq =
//         GetEnvironment().Get("NCBI_MAKEBLASTDB_DONT_SCAN_BIOSEQ_FOR_CFASTAREADER_USER_OBJECTS");
//     const bool scan_bioseq_4_cfastareader_usrobj = static_cast<bool>(dont_scan_bioseq.empty());
// ```
/// Runs `makeblastdb` with the words after the program name, as on the command line, which
/// writes the database's files as there. `input` is what the command line reads as standard input
/// (`-in -`, and without `-in`), so a database can be built from FASTA in memory. What the command
/// line writes to standard output goes to `output`, what it writes to standard error to
/// `diagnostics`. A failure is an error whose message is that text; an argv that the parser
/// refuses is `Error::Argv`.
pub fn run_makeblastdb(
    argv: &[&str],
    input: &mut dyn Read,
    output: &mut dyn Write,
    diagnostics: &mut dyn Write,
) -> Result<(), Error> {
    let words = ["losat", "makeblastdb"].iter().chain(argv).copied();
    let cli: Cli = crate::cli::try_parse_from(words)
        .map_err(|error| Error::Argv(crate::cli::render_message(&error)))?;
    let Commands::Makedb(args) = cli.command else {
        unreachable!("the first word names makeblastdb");
    };
    // NCBI's application layer reads the environment and the registry files before `Run`, and
    // LOSAT rejects the settings that change the output (`cli/main_entry.rs`).
    let settings = check_ncbi_application_settings("makeblastdb").map_err(Error::Search)?;
    makeblastdb::check_registry_settings(
        &args,
        settings
            .max_seqid_length
            .as_ref()
            .map(Option::as_deref)
            .map_err(String::as_str),
    )
    .map_err(Error::Search)?;
    makeblastdb::check_environment(&|name: &str| std::env::var_os(name)).map_err(Error::Search)?;
    let now = LocalTime::now();
    let mut recorded = Recorded {
        writer: diagnostics,
        bytes: Vec::new(),
    };
    let status = makeblastdb::run_with(
        &args,
        now.as_ref(),
        Platform::host(),
        input,
        output,
        &mut recorded,
    );
    outcome(status, recorded.bytes)
}
