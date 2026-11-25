//! `losat cleanup-blastdb-volumes.py`: pure-Rust port of NCBI BLAST+ 2.17.0's
//! `cleanup-blastdb-volumes.py` (`c++/src/app/blast/cleanup-blastdb-volumes.py`), which removes
//! orphaned BLAST database volume files not referenced by a database's `.pal` or `.nal` alias
//! file.
//!
//! Because NCBI's utility is a standalone Python script (matched to Python 3.13) rather than an
//! NCBI C++ Toolkit `CNcbiApplication`, it does not read `-logfile`, `-conffile`, `-dryrun`,
//! `-xmlhelp`, or `BLAST_USAGE_REPORT`. Its `-version` prints LOSAT's version
//! (`crate::wording::LOSAT_VERSION`) in place of NCBI's script version
//! (`PD-LOSAT-CLI-NONSEARCH-DIFFERENCES` exception 14).

mod args;
mod config;
mod glob;
mod work;

use std::ffi::OsString;
use std::io::Write;

pub(crate) use work::Env;

/// Buffer threshold matching Python 3.13's `TextIOWrapper` / `BufferedWriter` default block size
/// (8 KiB) for non-terminal standard output (pipes and regular files).
const PY_STDOUT_CHUNK_SIZE: usize = 8192;

/// Python 3.13-compatible stdout buffer: 8 KiB block-buffered on pipes/files, line-buffered on a
/// terminal. Mid-run flush failures stop immediately with exit status `1`, while exit-time flush
/// failures print Python's `Exception ignored on flushing sys.stdout:` message and exit with
/// status `120` (observed from NCBI BLAST+ 2.17.0 under Python 3.13.3).
pub(crate) struct PyStdout<'a> {
    inner: &'a mut dyn Write,
    buf: Vec<u8>,
    line_buffered: bool,
    failed_mid_run: bool,
}

impl<'a> PyStdout<'a> {
    pub(crate) fn new(inner: &'a mut dyn Write, line_buffered: bool) -> Self {
        Self {
            inner,
            buf: Vec::with_capacity(PY_STDOUT_CHUNK_SIZE),
            line_buffered,
            failed_mid_run: false,
        }
    }

    fn flush_now(&mut self) -> std::io::Result<()> {
        let res = self
            .inner
            .write_all(&self.buf)
            .and_then(|()| self.inner.flush());
        self.buf.clear();
        if res.is_err() {
            self.failed_mid_run = true;
        }
        res
    }

    fn maybe_flush(&mut self, has_newline: bool) -> std::io::Result<()> {
        if self.buf.len() >= PY_STDOUT_CHUNK_SIZE || (self.line_buffered && has_newline) {
            self.flush_now()?;
        }
        Ok(())
    }

    /// Writes `text` (such as `-h` or `-version`) as a single stdout write call.
    pub(crate) fn write_all_text(&mut self, text: &[u8]) -> std::io::Result<()> {
        self.buf.extend_from_slice(text);
        self.maybe_flush(text.contains(&b'\n'))
    }

    /// Writes `prefix` + `path` followed by `\n`, mirroring Python's `print` (writing the line
    /// text and then the trailing newline as two stdout writes).
    pub(crate) fn print_line(&mut self, prefix: &[u8], path: &[u8]) -> std::io::Result<()> {
        self.buf.extend_from_slice(prefix);
        self.buf.extend_from_slice(path);
        self.maybe_flush(prefix.contains(&b'\n') || path.contains(&b'\n'))?;
        self.buf.push(b'\n');
        self.maybe_flush(true)
    }

    /// Flushes any remaining buffered stdout bytes at program exit. Matches Python 3.13's exit
    /// flush: if flushing fails, prints `Exception ignored on flushing sys.stdout:` + the OS error
    /// line to `err` and overrides the exit status to `120`.
    pub(crate) fn finish_flush(&mut self, current_rc: i32, err: &mut dyn Write) -> i32 {
        if self.failed_mid_run || self.buf.is_empty() {
            return current_rc;
        }
        let res = self
            .inner
            .write_all(&self.buf)
            .and_then(|()| self.inner.flush());
        self.buf.clear();
        match res {
            Ok(()) => current_rc,
            Err(e) => {
                let _ = writeln!(
                    err,
                    "Exception ignored on flushing sys.stdout:\n{}",
                    format_py_os_error(&e)
                );
                let _ = err.flush();
                120
            }
        }
    }
}

/// Formats `e` like Python 3.13's `BrokenPipeError` / `PermissionError` / `OSError` line
/// (`BrokenPipeError: [Errno 32] Broken pipe`, `OSError: [Errno 28] No space left on device`).
fn format_py_os_error(e: &std::io::Error) -> String {
    let exc = match e.kind() {
        std::io::ErrorKind::BrokenPipe => "BrokenPipeError",
        std::io::ErrorKind::PermissionDenied => "PermissionError",
        _ => "OSError",
    };
    if let Some(code) = e.raw_os_error() {
        let raw = e.to_string();
        let suffix = format!(" (os error {code})");
        let desc = raw.strip_suffix(&suffix).unwrap_or(&raw);
        format!("{exc}: [Errno {code}] {desc}")
    } else {
        format!("{exc}: {e}")
    }
}

/// Reports a mid-run output write failure on `err` and returns exit status `1`.
pub(crate) fn write_failed(err: &mut dyn Write, e: &std::io::Error) -> i32 {
    let _ = writeln!(err, "Error: failed to write output: {e}");
    let _ = err.flush();
    1
}

// NCBI reference (598d8ae6): c++/src/app/blast/cleanup-blastdb-volumes.py:104-161
/// Runs `cleanup-blastdb-volumes.py` with `words` (the command line after
/// `losat cleanup-blastdb-volumes.py`) on the process's environment and standard streams; returns
/// the exit status.
pub(crate) fn run(words: &[OsString]) -> i32 {
    let env = Env::process();
    let mut out = crate::cli::standard_output();
    let mut err = std::io::stderr();
    let line_buffered = std::io::IsTerminal::is_terminal(&std::io::stdout());
    run_with_buffering(words, &env, &mut out, line_buffered, &mut err)
}

// NCBI reference (598d8ae6): c++/src/app/blast/cleanup-blastdb-volumes.py:104-141
/// Runs `cleanup-blastdb-volumes.py` with `words`, `env`, block-buffered `out`, and `err`; returns
/// the exit status (`0` on success, `-h`, or `-version`; `1` when the database or alias file is
/// missing or a mid-run I/O or config error occurs; `2` on argument errors; `120` when exit-time
/// stdout flushing fails).
#[cfg(test)]
pub(crate) fn run_with(
    words: &[OsString],
    env: &Env,
    out: &mut dyn Write,
    err: &mut dyn Write,
) -> i32 {
    run_with_buffering(words, env, out, false, err)
}

// NCBI reference (598d8ae6): c++/src/app/blast/cleanup-blastdb-volumes.py:104-141
/// Runs `cleanup-blastdb-volumes.py` with explicit stdout line-buffering control (true for a
/// terminal, false for a pipe or file).
pub(crate) fn run_with_buffering(
    words: &[OsString],
    env: &Env,
    out: &mut dyn Write,
    line_buffered: bool,
    err: &mut dyn Write,
) -> i32 {
    let columns = args::resolve_columns(env.columns.as_deref(), env.stdout_columns);
    let mut stdout = PyStdout::new(out, line_buffered);
    let rc = match args::parse_with_columns(words, columns) {
        args::ParseOutcome::Help(text) | args::ParseOutcome::Version(text) => {
            match stdout.write_all_text(text.as_bytes()) {
                Ok(()) => 0,
                Err(e) => write_failed(err, &e),
            }
        }
        args::ParseOutcome::Error(text) => {
            let _ = err.write_all(text.as_bytes());
            let _ = err.flush();
            2
        }
        args::ParseOutcome::Run(parsed) => work::run(&parsed, env, &mut stdout, err),
    };
    stdout.finish_flush(rc, err)
}
