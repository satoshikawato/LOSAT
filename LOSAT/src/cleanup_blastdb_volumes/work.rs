//! Core database lookup and orphaned volume removal for `cleanup-blastdb-volumes.py`.
//!
//! NCBI reference (598d8ae6): c++/src/app/blast/cleanup-blastdb-volumes.py:67-142

use std::collections::VecDeque;
use std::io::{Read, Write};

use crate::utils::os_bytes::{os_str_bytes, path_bytes, path_of_bytes};

use super::args::ParsedArgs;
use super::config::{get_blastdb_from_ncbi_config, py_split_whitespace};
use super::glob::{glob, glob_sorted};
use super::{write_failed, PyStdout};

/// Host path separator byte (`/` on Unix/WASI, `\` on Windows).
const SEP: u8 = if cfg!(windows) { b'\\' } else { b'/' };

/// Chunk size used by Python 3.13's `TextIOWrapper` when decoding the alias file (`8 KiB`).
const ALIAS_CHUNK_SIZE: usize = 8192;

/// Queries the column width of standard output when attached to a terminal, or `None` otherwise
/// (matching Python's terminal-size lookup on stdout).
fn detect_stdout_columns() -> Option<usize> {
    #[cfg(unix)]
    {
        rustix::termios::tcgetwinsize(std::io::stdout())
            .ok()
            .map(|ws| usize::from(ws.ws_col))
            .filter(|&cols| cols > 0)
    }
    #[cfg(windows)]
    {
        use windows_sys::Win32::System::Console::{
            GetConsoleScreenBufferInfo, GetStdHandle, CONSOLE_SCREEN_BUFFER_INFO, STD_OUTPUT_HANDLE,
        };
        // SAFETY: `GetStdHandle` returns the process's stdout console handle (or null/invalid),
        // and `GetConsoleScreenBufferInfo` writes only into our initialized stack `info` struct
        // and returns 0 when stdout is not a Windows console.
        unsafe {
            let handle = GetStdHandle(STD_OUTPUT_HANDLE);
            if handle.is_null() {
                return None;
            }
            let mut info = CONSOLE_SCREEN_BUFFER_INFO::default();
            if GetConsoleScreenBufferInfo(handle, &mut info) == 0 {
                return None;
            }
            let cols = i32::from(info.srWindow.Right) - i32::from(info.srWindow.Left) + 1;
            usize::try_from(cols).ok().filter(|&c| c > 0)
        }
    }
    #[cfg(not(any(unix, windows)))]
    {
        None
    }
}

/// Environment inputs needed by `cleanup-blastdb-volumes.py`, injected so unit tests do not mutate
/// process state.
#[derive(Clone, Debug, Default)]
pub(crate) struct Env {
    /// Current working directory (`Some(Err(_))` when `current_dir()` fails, e.g. removed cwd).
    pub(crate) cwd: Option<Result<Vec<u8>, String>>,
    pub(crate) home: Option<Vec<u8>>,
    pub(crate) blastdb: Option<Vec<u8>>,
    pub(crate) ncbi: Option<Vec<u8>>,
    /// Value of the `COLUMNS` environment variable, if set.
    pub(crate) columns: Option<String>,
    /// Terminal column width of standard output, if attached to a terminal.
    pub(crate) stdout_columns: Option<usize>,
}

impl Env {
    /// Captures the current process's working directory, home directory, `BLASTDB`, `NCBI`,
    /// `COLUMNS`, and stdout terminal width.
    pub(crate) fn process() -> Self {
        #[cfg(unix)]
        let home = match std::env::var_os("HOME") {
            Some(v) if v.is_empty() => Some(b"/".to_vec()),
            Some(v) => Some(os_str_bytes(&v).into_owned()),
            None => std::env::home_dir().map(|p| path_bytes(&p).into_owned()),
        };
        #[cfg(not(unix))]
        let home = std::env::home_dir().map(|p| path_bytes(&p).into_owned());
        let cwd = Some(
            std::env::current_dir()
                .map(|p| path_bytes(&p).into_owned())
                .map_err(|e| format!("cannot determine current working directory: {e}")),
        );
        Self {
            cwd,
            home,
            blastdb: std::env::var_os("BLASTDB").map(|v| os_str_bytes(&v).into_owned()),
            ncbi: std::env::var_os("NCBI").map(|v| os_str_bytes(&v).into_owned()),
            columns: std::env::var_os("COLUMNS").map(|v| v.to_string_lossy().into_owned()),
            stdout_columns: detect_stdout_columns(),
        }
    }
}

/// Returns whether `b` is a path separator on the host OS.
fn is_sep(b: u8) -> bool {
    b == b'/' || (cfg!(windows) && b == b'\\')
}

/// Returns the length of a Windows drive prefix (`C:`), or 0 otherwise.
fn drive_prefix_len(path: &[u8]) -> usize {
    if cfg!(windows) && path.len() >= 2 && path[0].is_ascii_alphabetic() && path[1] == b':' {
        2
    } else {
        0
    }
}

/// Returns whether `path` is absolute (matching Python's `os.path.isabs`).
fn is_abs(path: &[u8]) -> bool {
    if cfg!(windows) {
        let drive = drive_prefix_len(path);
        path.get(drive).is_some_and(|&b| is_sep(b))
    } else {
        path.starts_with(b"/")
    }
}

/// Joins `dir` and `name` following Python's `os.path.join`.
pub(crate) fn path_join(dir: &[u8], name: &[u8]) -> Vec<u8> {
    if is_abs(name) || dir.is_empty() {
        return name.to_vec();
    }
    let name_drive = drive_prefix_len(name);
    if name_drive > 0 {
        let dir_drive = drive_prefix_len(dir);
        if dir_drive == 0 || !dir[0].eq_ignore_ascii_case(&name[0]) {
            return name.to_vec();
        }
        return path_join(dir, &name[name_drive..]);
    }
    let mut out = Vec::with_capacity(dir.len() + 1 + name.len());
    out.extend_from_slice(dir);
    let drive_only = drive_prefix_len(dir) == dir.len();
    if !drive_only && !dir.last().is_some_and(|&b| is_sep(b)) {
        out.push(SEP);
    }
    out.extend_from_slice(name);
    out
}

/// Returns the final path component following Python's `os.path.basename`.
fn basename(path: &[u8]) -> &[u8] {
    let rest = &path[drive_prefix_len(path)..];
    match rest.iter().rposition(|&b| is_sep(b)) {
        Some(pos) => &rest[pos + 1..],
        None => rest,
    }
}

/// Whether `path` exists, following symbolic links (matching Python's `os.path.exists`).
fn exists(path: &[u8]) -> bool {
    !path.is_empty() && std::fs::metadata(path_of_bytes(path)).is_ok()
}

/// Builds `<base>.<p|n><tail>` from raw path bytes.
fn with_ext(base: &[u8], letter: u8, tail: &[u8]) -> Vec<u8> {
    let mut out = Vec::with_capacity(base.len() + 2 + tail.len());
    out.extend_from_slice(base);
    out.push(b'.');
    out.push(letter);
    out.extend_from_slice(tail);
    out
}

// NCBI reference (598d8ae6): c++/src/app/blast/cleanup-blastdb-volumes.py:67-94
/// Locates the BLAST database `name` by checking:
/// 1. Direct path (`{name}.{p|n}al` or `{name}.{p|n}in`);
/// 2. `$BLASTDB` joined with `name`;
/// 3. `[BLAST] BLASTDB` in `.ncbirc` and `ncbi.ini` under `cwd`, `home`, and `$NCBI`.
pub(crate) fn find_blastdb(
    name: &[u8],
    is_prot: bool,
    env: &Env,
) -> Result<Option<Vec<u8>>, String> {
    let letter = if is_prot { b'p' } else { b'n' };
    let alias_file = with_ext(name, letter, b"al");
    let index_file = with_ext(name, letter, b"in");
    if exists(&alias_file) || exists(&index_file) {
        return Ok(Some(name.to_vec()));
    }

    if let Some(blastdb) = &env.blastdb {
        let alf = path_join(blastdb, &alias_file);
        let idxf = path_join(blastdb, &index_file);
        if exists(&alf) || exists(&idxf) {
            return Ok(Some(path_join(blastdb, name)));
        }
    }

    let mut search_dirs: Vec<&[u8]> = Vec::with_capacity(3);
    if let Some(cwd_res) = &env.cwd {
        let cwd = cwd_res.as_ref().map_err(Clone::clone)?;
        search_dirs.push(cwd);
    }
    if let Some(home) = &env.home {
        search_dirs.push(home);
    }
    if let Some(ncbi) = &env.ncbi {
        search_dirs.push(ncbi);
    }

    for dir in search_dirs {
        for fname in [b".ncbirc".as_slice(), b"ncbi.ini".as_slice()] {
            let ncbirc = path_join(dir, fname);
            if exists(&ncbirc) {
                let path = path_of_bytes(&ncbirc);
                if let Some(blastdb) = get_blastdb_from_ncbi_config(&path)? {
                    let alf = path_join(&blastdb, &alias_file);
                    let idxf = path_join(&blastdb, &index_file);
                    if exists(&alf) || exists(&idxf) {
                        return Ok(Some(path_join(&blastdb, name)));
                    }
                }
            }
        }
    }

    Ok(None)
}

/// Error raised while streaming an alias file in 8 KiB UTF-8 chunks.
enum AliasReadError {
    Io(std::io::Error),
    Utf8(String),
}

/// Streams lines from an alias file in 8 KiB UTF-8 chunks with universal newlines (`\r\n`, `\n`,
/// `\r`), matching Python 3.13's `open(alias_file, "rt")` so `DBLIST` lines in chunks preceding an
/// invalid UTF-8 byte are yielded and processed before failing.
struct AliasLineReader<R> {
    reader: R,
    pending_bytes: Vec<u8>,
    decoded_lines: VecDeque<String>,
    partial_line: String,
    pending_cr: bool,
    eof: bool,
}

impl<R: Read> AliasLineReader<R> {
    fn new(reader: R) -> Self {
        Self {
            reader,
            pending_bytes: Vec::new(),
            decoded_lines: VecDeque::new(),
            partial_line: String::new(),
            pending_cr: false,
            eof: false,
        }
    }

    fn read_chunk(&mut self, buf: &mut [u8; ALIAS_CHUNK_SIZE]) -> std::io::Result<usize> {
        let mut filled = 0;
        while filled < buf.len() {
            match self.reader.read(&mut buf[filled..]) {
                Ok(0) => break,
                Ok(n) => filled += n,
                Err(e) if e.kind() == std::io::ErrorKind::Interrupted => continue,
                Err(e) => return Err(e),
            }
        }
        Ok(filled)
    }

    fn emit_partial(&mut self) {
        let finished = std::mem::take(&mut self.partial_line);
        self.decoded_lines.push_back(finished);
    }

    fn push_decoded_chunk(&mut self, text: &str) {
        let mut slice = text;
        if self.pending_cr {
            self.pending_cr = false;
            self.emit_partial();
            if let Some(rest) = slice.strip_prefix('\n') {
                slice = rest;
            }
        }
        let bytes = slice.as_bytes();
        let mut start = 0usize;
        let mut i = 0usize;
        while i < bytes.len() {
            if bytes[i] == b'\n' || bytes[i] == b'\r' {
                self.partial_line.push_str(&slice[start..i]);
                if bytes[i] == b'\n' {
                    self.emit_partial();
                    i += 1;
                } else if i + 1 == bytes.len() {
                    // Wait for the next chunk (or EOF) to see whether `\n` follows this `\r`.
                    self.pending_cr = true;
                    i += 1;
                } else {
                    self.emit_partial();
                    i += if bytes[i + 1] == b'\n' { 2 } else { 1 };
                }
                start = i;
            } else {
                i += 1;
            }
        }
        if start < bytes.len() {
            self.partial_line.push_str(&slice[start..]);
        }
    }

    fn next_line(&mut self) -> Result<Option<String>, AliasReadError> {
        if let Some(line) = self.decoded_lines.pop_front() {
            return Ok(Some(line));
        }
        let mut raw = [0u8; ALIAS_CHUNK_SIZE];
        while self.decoded_lines.is_empty() && !self.eof {
            let n = self.read_chunk(&mut raw).map_err(AliasReadError::Io)?;
            if n == 0 {
                self.eof = true;
                if !self.pending_bytes.is_empty() {
                    return Err(AliasReadError::Utf8(
                        "unexpected end of UTF-8 sequence".to_string(),
                    ));
                }
                if self.pending_cr || !self.partial_line.is_empty() {
                    self.pending_cr = false;
                    return Ok(Some(std::mem::take(&mut self.partial_line)));
                }
                return Ok(None);
            }

            let mut combined = Vec::with_capacity(self.pending_bytes.len() + n);
            combined.extend_from_slice(&self.pending_bytes);
            combined.extend_from_slice(&raw[..n]);
            self.pending_bytes.clear();

            match std::str::from_utf8(&combined) {
                Ok(valid) => self.push_decoded_chunk(valid),
                Err(err) => {
                    if err.error_len().is_none() {
                        let good = err.valid_up_to();
                        let valid =
                            std::str::from_utf8(&combined[..good]).expect("valid UTF-8 prefix");
                        self.pending_bytes.extend_from_slice(&combined[good..]);
                        self.push_decoded_chunk(valid);
                    } else {
                        return Err(AliasReadError::Utf8(err.to_string()));
                    }
                }
            }
        }
        Ok(self.decoded_lines.pop_front())
    }
}

/// Processes a single `DBLIST` line from the alias file; returns `Some(rc)` to exit early.
fn process_dblist_line(
    line: &str,
    vol_glob: &[u8],
    args: &ParsedArgs,
    out: &mut PyStdout<'_>,
    err: &mut dyn Write,
) -> Option<i32> {
    let vols: Vec<String> = py_split_whitespace(line)
        .skip(1)
        .map(|token| token.replace('"', ""))
        .collect();

    for existing_vol in glob_sorted(vol_glob) {
        let base = basename(&existing_vol);
        if base.len() < 4 || existing_vol.len() < 4 {
            continue;
        }
        let vol_name = &base[..base.len() - 4];
        if vols.iter().any(|v| v.as_bytes() == vol_name) {
            continue;
        }
        let vol_stem = &existing_vol[..existing_vol.len() - 4];
        if args.dry_run {
            if let Err(e) = out.print_line(b"Will remove extra volume ", vol_stem) {
                return Some(write_failed(err, &e));
            }
        }

        let mut ext_glob = Vec::with_capacity(existing_vol.len());
        ext_glob.extend_from_slice(&existing_vol[..existing_vol.len() - 2]);
        ext_glob.extend_from_slice(b"??");
        let mut md5_glob = Vec::with_capacity(vol_stem.len() + 11);
        md5_glob.extend_from_slice(vol_stem);
        md5_glob.extend_from_slice(b".tar.gz.md5");

        // Keep `read_dir` order within each extra volume before `.tar.gz.md5`, matching Python's
        // `glob.glob` (`os.scandir`).
        let mut to_rm = glob(&ext_glob);
        to_rm.extend(glob(&md5_glob));

        for f in to_rm {
            if !args.dry_run {
                if let Err(e) = std::fs::remove_file(path_of_bytes(&f)) {
                    let _ = writeln!(
                        err,
                        "Error: cannot remove {}: {e}",
                        String::from_utf8_lossy(&f)
                    );
                    let _ = err.flush();
                    return Some(1);
                }
                if let Err(e) = out.print_line(b"Removed ", &f) {
                    return Some(write_failed(err, &e));
                }
            } else if args.verbose > 0 {
                if let Err(e) = out.print_line(b"Will remove ", &f) {
                    return Some(write_failed(err, &e));
                }
            }
        }
    }
    None
}

// NCBI reference (598d8ae6): c++/src/app/blast/cleanup-blastdb-volumes.py:104-142
/// Executes the cleanup logic after argument parsing.
pub(crate) fn run(
    args: &ParsedArgs,
    env: &Env,
    out: &mut PyStdout<'_>,
    err: &mut dyn Write,
) -> i32 {
    let ext = args.dbtype.ext_letter();
    let db = match find_blastdb(&args.db, args.dbtype.is_prot(), env) {
        Ok(Some(db)) => db,
        Ok(None) => {
            let _ = err.write_all(b"Cannot find ");
            let _ = err.write_all(args.dbtype.label().as_bytes());
            let _ = err.write_all(b" ");
            let _ = err.write_all(&args.db);
            let _ = err.write_all(b" BLAST database\n");
            let _ = err.flush();
            return 1;
        }
        Err(msg) => {
            let _ = writeln!(err, "Error: {msg}");
            let _ = err.flush();
            return 1;
        }
    };

    let alias_file = with_ext(&db, ext, b"al");
    if !exists(&alias_file) {
        return 1;
    }

    let file = match std::fs::File::open(path_of_bytes(&alias_file)) {
        Ok(file) => file,
        Err(e) => {
            let _ = writeln!(
                err,
                "Error: cannot read {}: {e}",
                String::from_utf8_lossy(&alias_file)
            );
            let _ = err.flush();
            return 1;
        }
    };

    let mut vol_glob = Vec::with_capacity(db.len() + 6);
    vol_glob.extend_from_slice(&db);
    vol_glob.extend_from_slice(b".*.");
    vol_glob.push(ext);
    vol_glob.extend_from_slice(b"in");

    let mut reader = AliasLineReader::new(file);
    loop {
        let line = match reader.next_line() {
            Ok(Some(line)) => line,
            Ok(None) => break,
            Err(AliasReadError::Io(e)) => {
                let _ = writeln!(
                    err,
                    "Error: cannot read {}: {e}",
                    String::from_utf8_lossy(&alias_file)
                );
                let _ = err.flush();
                return 1;
            }
            Err(AliasReadError::Utf8(e)) => {
                let _ = writeln!(
                    err,
                    "Error: {}: invalid UTF-8: {e}",
                    String::from_utf8_lossy(&alias_file)
                );
                let _ = err.flush();
                return 1;
            }
        };
        if !line.starts_with("DBLIST") {
            continue;
        }
        if let Some(rc) = process_dblist_line(&line, &vol_glob, args, out, err) {
            return rc;
        }
    }

    0
}

#[cfg(test)]
#[path = "work_tests.rs"]
mod tests;
