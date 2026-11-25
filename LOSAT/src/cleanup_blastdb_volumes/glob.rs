//! Pure-Rust byte-path glob matcher reproducing Python's `glob.glob` behavior for
//! `cleanup-blastdb-volumes.py` (`{db}.*.{ext}in`, `{vol}.p??`, `{vol}.tar.gz.md5`).
//!
//! `glob` preserves `std::fs::read_dir` order (matching Python `glob.glob`), while `glob_sorted`
//! sorts using Python's `surrogateescape` code-point order (invalid Unix byte `b` mapped to
//! `0xDC00 + b`) for the outer `sorted(glob(...))` loop.
//!
//! NCBI reference (598d8ae6): c++/src/app/blast/cleanup-blastdb-volumes.py:126,132-133

use crate::utils::os_bytes::{os_str_bytes, path_of_bytes};

/// Returns whether `b` is a path separator for the given OS mode.
fn is_sep_os(b: u8, windows: bool) -> bool {
    b == b'/' || (windows && b == b'\\')
}

/// Returns whether `pattern` contains any glob metacharacter (`*`, `?`, `[`).
fn has_magic(pattern: &[u8]) -> bool {
    pattern.iter().any(|&b| matches!(b, b'*' | b'?' | b'['))
}

/// Whether `path` exists as a directory entry (including broken symlinks, matching Python's
/// `os.path.lexists`).
fn lexists(path: &[u8]) -> bool {
    !path.is_empty() && std::fs::symlink_metadata(path_of_bytes(path)).is_ok()
}

/// Whether `path` is an existing directory (following symlinks, matching Python's `os.path.isdir`).
fn is_dir(path: &[u8]) -> bool {
    !path.is_empty() && std::fs::metadata(path_of_bytes(path)).is_ok_and(|m| m.is_dir())
}

/// Joins a directory prefix (with its separator, if any) and a child component.
pub(crate) fn join_with_sep(dir: &[u8], sep: Option<u8>, name: &[u8]) -> Vec<u8> {
    let mut out = Vec::with_capacity(dir.len() + usize::from(sep.is_some()) + name.len());
    out.extend_from_slice(dir);
    if let Some(s) = sep {
        out.push(s);
    }
    out.extend_from_slice(name);
    out
}

/// Returns the length of a Windows drive (`C:`) or UNC (`\\server\share`) prefix, or 0 otherwise.
fn windows_drive_or_unc_len(path: &[u8]) -> usize {
    if path.len() >= 2 && path[0].is_ascii_alphabetic() && path[1] == b':' {
        return 2;
    }
    if path.len() >= 2
        && is_sep_os(path[0], true)
        && is_sep_os(path[1], true)
        && !path.get(2).is_some_and(|&b| is_sep_os(b, true))
    {
        let rest = &path[2..];
        if let Some(s1) = rest.iter().position(|&b| is_sep_os(b, true)) {
            let after_host = &rest[s1 + 1..];
            if !after_host.is_empty() && !is_sep_os(after_host[0], true) {
                let s2 = after_host
                    .iter()
                    .position(|&b| is_sep_os(b, true))
                    .unwrap_or(after_host.len());
                return 2 + s1 + 1 + s2;
            }
        }
    }
    0
}

/// Splits `pattern` into `(dir_part, join_separator, base_part)` matching Python's `os.path.split`
/// + `os.path.join` (`posixpath` on Unix, `ntpath` on Windows: Windows joins non-root prefixes
/// with `\` and collapses multiple root separators after a drive/UNC prefix to one).
pub(crate) fn split_dir_base_os(pattern: &[u8], windows: bool) -> (&[u8], Option<u8>, &[u8]) {
    let drive = if windows {
        windows_drive_or_unc_len(pattern)
    } else {
        0
    };
    let rest = &pattern[drive..];
    let Some(rel_pos) = rest.iter().rposition(|&b| is_sep_os(b, windows)) else {
        if drive > 0 {
            return (&pattern[..drive], None, rest);
        }
        return (b"", None, pattern);
    };
    let pos = drive + rel_pos;
    let after_drive = &pattern[drive..=pos];
    if after_drive.iter().all(|&b| is_sep_os(b, windows)) {
        if windows {
            // `ntpath.split`: keep the drive/UNC prefix plus its single root separator (e.g.
            // `C://p5` -> `("C:/", "p5")`, `/p5` -> `("/", "p5")`, `\\srv\sh\p5` -> `("\\\\srv\\sh\\", "p5")`),
            // except a bare `//foo` (no share) which `ntpath.split` treats as UNC root with empty
            // tail.
            if drive == 0
                && after_drive.len() >= 2
                && pattern[pos + 1..].iter().all(|&b| !is_sep_os(b, true))
            {
                return (pattern, Some(b'\\'), b"");
            }
            return (&pattern[..=drive], None, &pattern[pos + 1..]);
        }
        // `posixpath.split`: root prefix such as `/`, `//`, `///` stays intact with no extra `/`.
        return (&pattern[..=pos], None, &pattern[pos + 1..]);
    }
    let mut dir_end = pos;
    while dir_end > drive && is_sep_os(pattern[dir_end - 1], windows) {
        dir_end -= 1;
    }
    let join_sep = if windows { b'\\' } else { b'/' };
    (&pattern[..dir_end], Some(join_sep), &pattern[pos + 1..])
}

/// One element of a filename or glob pattern: a decoded Unicode scalar value, or a raw non-UTF-8
/// byte on Unix (mapped to Python's `surrogateescape` range `U+DC80`–`U+DCFF`).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Atom {
    Char(char),
    Byte(u8),
}

impl Atom {
    /// Returns Python's `ord()` code point for `self`, placing invalid Unix bytes `b` in the
    /// `surrogateescape` range `0xDC00 + b` (`U+DC80`–`U+DCFF`).
    fn as_u32(self) -> u32 {
        match self {
            Self::Char(c) => c as u32,
            Self::Byte(b) => 0xDC00 + u32::from(b),
        }
    }
}

/// Decodes `bytes` into `Atom`s (valid UTF-8 scalars, with raw `Atom::Byte` for invalid bytes).
/// When `fold_case` is true (Windows), valid UTF-8 chunks are Unicode-lowercased via
/// `str::to_lowercase`, matching Python's `ntpath.normcase`.
fn to_atoms(bytes: &[u8], fold_case: bool) -> Vec<Atom> {
    let mut out = Vec::with_capacity(bytes.len());
    let mut i = 0;
    let push_valid = |s: &str, out: &mut Vec<Atom>| {
        if fold_case {
            out.extend(s.to_lowercase().chars().map(Atom::Char));
        } else {
            out.extend(s.chars().map(Atom::Char));
        }
    };
    while i < bytes.len() {
        match std::str::from_utf8(&bytes[i..]) {
            Ok(valid) => {
                push_valid(valid, &mut out);
                break;
            }
            Err(err) => {
                let good = err.valid_up_to();
                if good > 0 {
                    let prefix = std::str::from_utf8(&bytes[i..i + good]).expect("valid prefix");
                    push_valid(prefix, &mut out);
                    i += good;
                }
                let bad_len = err.error_len().unwrap_or(1);
                for &b in &bytes[i..i + bad_len] {
                    out.push(Atom::Byte(b));
                }
                i += bad_len;
            }
        }
    }
    out
}

/// Matches a bracket class `[...]` starting at `pat[p]` (where `pat[p] == Atom::Char('[')`).
/// Returns `Some((matched, next_p))` if a closing `]` exists, or `None` if unclosed (treat `[` as
/// literal, matching Python `fnmatch`).
fn match_bracket(pat: &[Atom], p: usize, ch: Atom) -> Option<(bool, usize)> {
    let mut j = p + 1;
    let negated = pat.get(j) == Some(&Atom::Char('!'));
    if negated {
        j += 1;
    }
    // A `]` immediately after `[` or `[!` is a literal member of the class.
    if pat.get(j) == Some(&Atom::Char(']')) {
        j += 1;
    }
    while j < pat.len() && pat[j] != Atom::Char(']') {
        j += 1;
    }
    if j >= pat.len() {
        return None;
    }
    let start = p + 1 + usize::from(negated);
    let body = &pat[start..j];
    let val = ch.as_u32();
    let mut found = false;
    let mut k = 0;
    while k < body.len() {
        if k + 2 < body.len() && body[k + 1] == Atom::Char('-') {
            let lo = body[k].as_u32();
            let hi = body[k + 2].as_u32();
            if lo <= val && val <= hi {
                found = true;
            }
            k += 3;
        } else {
            if body[k] == ch {
                found = true;
            }
            k += 1;
        }
    }
    Some((if negated { !found } else { found }, j + 1))
}

/// Matches a single path component `name` against glob `pattern` with explicit case-folding mode
/// (`true` on Windows for `ntpath.normcase`, `false` on Unix).
pub(crate) fn fnmatch_bytes_case(pattern: &[u8], name: &[u8], fold_case: bool) -> bool {
    let pat = to_atoms(pattern, fold_case);
    let txt = to_atoms(name, fold_case);
    let mut memo = vec![None; (pat.len() + 1) * (txt.len() + 1)];
    fn dp(pat: &[Atom], txt: &[Atom], p: usize, t: usize, memo: &mut [Option<bool>]) -> bool {
        let key = p * (txt.len() + 1) + t;
        if let Some(cached) = memo[key] {
            return cached;
        }
        let ans = if p == pat.len() {
            t == txt.len()
        } else if pat[p] == Atom::Char('*') {
            dp(pat, txt, p + 1, t, memo) || (t < txt.len() && dp(pat, txt, p, t + 1, memo))
        } else if t == txt.len() {
            false
        } else if pat[p] == Atom::Char('?') {
            dp(pat, txt, p + 1, t + 1, memo)
        } else if let Some((matched, next_p)) = (pat[p] == Atom::Char('['))
            .then(|| match_bracket(pat, p, txt[t]))
            .flatten()
        {
            matched && dp(pat, txt, next_p, t + 1, memo)
        } else {
            pat[p] == txt[t] && dp(pat, txt, p + 1, t + 1, memo)
        };
        memo[key] = Some(ans);
        ans
    }
    dp(&pat, &txt, 0, 0, &mut memo)
}

/// Matches a single path component `name` against glob `pattern` using host OS case sensitivity.
pub(crate) fn fnmatch_bytes(pattern: &[u8], name: &[u8]) -> bool {
    fnmatch_bytes_case(pattern, name, cfg!(windows))
}

/// Compares two raw path byte slices in Python's `surrogateescape` string order (`U+DC80`–`U+DCFF`
/// for invalid UTF-8 bytes).
pub(crate) fn cmp_py_path(a: &[u8], b: &[u8]) -> std::cmp::Ordering {
    to_atoms(a, false)
        .into_iter()
        .map(Atom::as_u32)
        .cmp(to_atoms(b, false).into_iter().map(Atom::as_u32))
}

/// Expands `pattern` as Python's `sorted(glob.glob(pattern))` does (sorted in Python's
/// `surrogateescape` code-point order).
pub(crate) fn glob_sorted(pattern: &[u8]) -> Vec<Vec<u8>> {
    let mut matches = glob(pattern);
    matches.sort_by(|a, b| cmp_py_path(a, b));
    matches
}

/// Expands `pattern` as Python's `glob.glob(pattern)` does, preserving filesystem directory
/// (`read_dir`) iteration order.
pub(crate) fn glob(pattern: &[u8]) -> Vec<Vec<u8>> {
    if pattern.is_empty() {
        return Vec::new();
    }
    if !has_magic(pattern) {
        return if lexists(pattern) {
            vec![pattern.to_vec()]
        } else {
            Vec::new()
        };
    }

    let (dir_part, sep, base_part) = split_dir_base_os(pattern, cfg!(windows));
    let dirs: Vec<Vec<u8>> = if dir_part.is_empty() {
        vec![Vec::new()]
    } else if has_magic(dir_part) {
        glob(dir_part).into_iter().filter(|d| is_dir(d)).collect()
    } else {
        vec![dir_part.to_vec()]
    };

    let mut out = Vec::new();
    for dir in dirs {
        if base_part.is_empty() {
            let candidate = join_with_sep(&dir, sep, b"");
            if is_dir(&candidate) {
                out.push(candidate);
            }
        } else if !has_magic(base_part) {
            let candidate = join_with_sep(&dir, sep, base_part);
            if lexists(&candidate) {
                out.push(candidate);
            }
        } else {
            let open_dir = if dir.is_empty() {
                std::fs::read_dir(".")
            } else {
                std::fs::read_dir(path_of_bytes(&dir))
            };
            let Ok(entries) = open_dir else {
                continue;
            };
            let allow_hidden = base_part.starts_with(b".");
            for entry in entries.flatten() {
                let name = os_str_bytes(&entry.file_name()).into_owned();
                if name.starts_with(b".") && !allow_hidden {
                    continue;
                }
                if fnmatch_bytes(base_part, &name) {
                    out.push(join_with_sep(&dir, sep, &name));
                }
            }
        }
    }
    out
}

#[cfg(test)]
#[path = "glob_tests.rs"]
mod tests;
