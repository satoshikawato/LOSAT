//! The names that run a program without its subcommand word (`losat --install-shims`).

use std::ffi::{OsStr, OsString};

/// The names that run a program without its subcommand word (`entry_argv`), as `losat
/// --install-shims` links them. No NCBI counterpart: NCBI has one executable per program.
pub const ENTRY_NAMES: [&str; 12] = [
    "blastn",
    "blastp",
    "blastx",
    "tblastn",
    "tblastx",
    "psiblast",
    "rpsblast",
    "makeblastdb",
    "blastdbcmd",
    "dustmasker",
    "segmasker",
    "rpsbproc",
];

/// The program of an executable called `argv0`: its file name without `.exe`, when that is one
/// of `ENTRY_NAMES`.
fn entry_program(argv0: &OsStr) -> Option<&'static str> {
    let path = std::path::Path::new(argv0);
    if !path
        .extension()
        .is_none_or(|extension| extension.eq_ignore_ascii_case("exe"))
    {
        return None;
    }
    let stem = path.file_stem()?.to_str()?;
    // Windows finds `blastn.exe` for the word `BLASTN`, which is then `argv[0]`.
    #[cfg(windows)]
    let stem = stem.to_ascii_lowercase();
    ENTRY_NAMES.into_iter().find(|name| *name == stem)
}

/// The words of the command line as `losat PROGRAM ...` reads them: an executable called by the
/// name of a program (`ENTRY_NAMES`) runs that program, with the texts of `losat PROGRAM`. Any
/// other name is `losat`, and its words stay as they are.
pub fn entry_argv<I, S>(argv: I) -> Vec<OsString>
where
    I: IntoIterator<Item = S>,
    S: Into<OsString>,
{
    let mut words: Vec<OsString> = argv.into_iter().map(Into::into).collect();
    if let Some(program) = words.first().and_then(|argv0| entry_program(argv0)) {
        words[0] = format!("losat{}", std::env::consts::EXE_SUFFIX).into();
        words.insert(1, program.into());
    }
    words
}

#[cfg(test)]
#[path = "entry_names_tests.rs"]
mod entry_names_tests;
