//! `losat cleanup-blastdb-volumes.py`: pure-Rust port of NCBI BLAST+ 2.17.0's
//! `cleanup-blastdb-volumes.py` (`c++/src/app/blast/cleanup-blastdb-volumes.py`), which removes
//! orphaned BLAST database volume files not referenced by a database's `.pal` or `.nal` alias
//! file.
//!
//! Because NCBI's utility is a standalone Python script rather than an NCBI C++ Toolkit
//! `CNcbiApplication`, it does not read `-logfile`, `-conffile`, `-dryrun`, `-xmlhelp`, or
//! `BLAST_USAGE_REPORT`. Its `-version` prints LOSAT's version (`crate::wording::LOSAT_VERSION`)
//! in place of NCBI's script version (`PD-LOSAT-CLI-NONSEARCH-DIFFERENCES` exception 14).

mod args;
mod config;
mod glob;
mod work;

use std::ffi::OsString;
use std::io::Write;

pub(crate) use work::Env;

// NCBI reference (598d8ae6): c++/src/app/blast/cleanup-blastdb-volumes.py:104-161
/// Runs `cleanup-blastdb-volumes.py` with `words` (the command line after
/// `losat cleanup-blastdb-volumes.py`) on the process's environment and standard streams; returns
/// the exit status.
pub(crate) fn run(words: &[OsString]) -> i32 {
    let env = Env::process();
    let mut out = crate::cli::standard_output();
    let mut err = std::io::stderr();
    run_with(words, &env, &mut out, &mut err)
}

// NCBI reference (598d8ae6): c++/src/app/blast/cleanup-blastdb-volumes.py:104-141
/// Runs `cleanup-blastdb-volumes.py` with `words`, `env`, `out`, and `err`; returns the exit
/// status (`0` on success, `-h`, or `-version`; `1` when the database or alias file is missing or
/// an I/O or config error occurs; `2` on argument errors).
pub(crate) fn run_with(
    words: &[OsString],
    env: &Env,
    out: &mut dyn Write,
    err: &mut dyn Write,
) -> i32 {
    match args::parse(words) {
        args::ParseOutcome::Help(text) | args::ParseOutcome::Version(text) => {
            let _ = out.write_all(text.as_bytes());
            let _ = out.flush();
            0
        }
        args::ParseOutcome::Error(text) => {
            let _ = err.write_all(text.as_bytes());
            let _ = err.flush();
            2
        }
        args::ParseOutcome::Run(parsed) => work::run(&parsed, env, out, err),
    }
}
