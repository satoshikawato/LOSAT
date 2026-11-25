//! Pure-Rust byte-path glob matcher reproducing Python's `glob.glob` behavior for
//! `cleanup-blastdb-volumes.py` (`{db}.*.{ext}in`, `{vol}.p??`, `{vol}.tar.gz.md5`).
//!
//! Results are sorted lexicographically by byte path so that volume and file output order is
//! deterministic across operating systems and filesystems.
//!
//! NCBI reference (598d8ae6): c++/src/app/blast/cleanup-blastdb-volumes.py:126,132-133

use crate::utils::os_bytes::{os_str_bytes, path_of_bytes};

/// Returns whether `b` is a path separator on the host OS.
fn is_sep(b: u8) -> bool {
    b == b'/' || (cfg!(windows) && b == b'\\')
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

/// Joins a directory prefix (with its original separator, if any) and a child component.
fn join_with_sep(dir: &[u8], sep: Option<u8>, name: &[u8]) -> Vec<u8> {
    let mut out = Vec::with_capacity(dir.len() + usize::from(sep.is_some()) + name.len());
    out.extend_from_slice(dir);
    if let Some(s) = sep {
        out.push(s);
    }
    out.extend_from_slice(name);
    out
}

/// Returns the length of a Windows drive prefix (`C:`), or 0 otherwise.
fn drive_prefix_len(path: &[u8]) -> usize {
    if cfg!(windows) && path.len() >= 2 && path[0].is_ascii_alphabetic() && path[1] == b':' {
        2
    } else {
        0
    }
}

/// Splits `pattern` into `(dir_part, separator, base_part)` at the last separator following
/// Python's `os.path.split` (collapsing redundant trailing separators on non-root directory
/// prefixes, while keeping root/drive-root prefixes intact).
fn split_dir_base(pattern: &[u8]) -> (&[u8], Option<u8>, &[u8]) {
    let drive = drive_prefix_len(pattern);
    let rest = &pattern[drive..];
    let Some(rel_pos) = rest.iter().rposition(|&b| is_sep(b)) else {
        if drive > 0 {
            return (&pattern[..drive], None, rest);
        }
        return (b"", None, pattern);
    };
    let pos = drive + rel_pos;
    let after_drive = &pattern[drive..=pos];
    if after_drive.iter().all(|&b| is_sep(b)) {
        // Root prefix such as `/`, `//`, `C:\`, `C:\\`: keep intact with no extra separator.
        (&pattern[..=pos], None, &pattern[pos + 1..])
    } else {
        let mut dir_end = pos;
        while dir_end > drive && is_sep(pattern[dir_end - 1]) {
            dir_end -= 1;
        }
        (&pattern[..dir_end], Some(pattern[pos]), &pattern[pos + 1..])
    }
}

/// One element of a filename or glob pattern: a decoded Unicode scalar value, or a raw non-UTF-8
/// byte on Unix.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Atom {
    Char(char),
    Byte(u8),
}

impl Atom {
    /// Normalizes ASCII case on Windows (matching Python `fnmatch` / `os.path.normcase`).
    fn norm(self) -> Self {
        if cfg!(windows) {
            match self {
                Self::Char(c) => Self::Char(c.to_ascii_lowercase()),
                Self::Byte(b) => Self::Byte(b.to_ascii_lowercase()),
            }
        } else {
            self
        }
    }

    fn as_u32(self) -> u32 {
        match self {
            Self::Char(c) => c as u32,
            Self::Byte(b) => 0x0011_0000 + u32::from(b),
        }
    }
}

/// Decodes `bytes` into `Atom`s (valid UTF-8 scalars, with raw `Atom::Byte` for invalid bytes).
fn to_atoms(bytes: &[u8]) -> Vec<Atom> {
    let mut out = Vec::with_capacity(bytes.len());
    let mut i = 0;
    while i < bytes.len() {
        match std::str::from_utf8(&bytes[i..]) {
            Ok(valid) => {
                out.extend(valid.chars().map(Atom::Char));
                break;
            }
            Err(err) => {
                let good = err.valid_up_to();
                if good > 0 {
                    let prefix = std::str::from_utf8(&bytes[i..i + good]).expect("valid prefix");
                    out.extend(prefix.chars().map(Atom::Char));
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
    let target = ch.norm();
    let mut found = false;
    let mut k = 0;
    while k < body.len() {
        if k + 2 < body.len() && body[k + 1] == Atom::Char('-') {
            let lo = body[k].norm().as_u32();
            let hi = body[k + 2].norm().as_u32();
            let val = target.as_u32();
            if lo <= val && val <= hi {
                found = true;
            }
            k += 3;
        } else {
            if body[k].norm() == target {
                found = true;
            }
            k += 1;
        }
    }
    Some((if negated { !found } else { found }, j + 1))
}

/// Matches a single path component `name` against glob `pattern` (Python `fnmatch` semantics).
pub(crate) fn fnmatch_bytes(pattern: &[u8], name: &[u8]) -> bool {
    let pat = to_atoms(pattern);
    let txt = to_atoms(name);
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
            pat[p].norm() == txt[t].norm() && dp(pat, txt, p + 1, t + 1, memo)
        };
        memo[key] = Some(ans);
        ans
    }
    dp(&pat, &txt, 0, 0, &mut memo)
}

/// Expands `pattern` as Python's `glob.glob(pattern)` does and sorts the matching byte paths.
pub(crate) fn glob_sorted(pattern: &[u8]) -> Vec<Vec<u8>> {
    let mut matches = glob_raw(pattern);
    matches.sort();
    matches
}

fn glob_raw(pattern: &[u8]) -> Vec<Vec<u8>> {
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

    let (dir_part, sep, base_part) = split_dir_base(pattern);
    let dirs: Vec<Vec<u8>> = if dir_part.is_empty() {
        vec![Vec::new()]
    } else if has_magic(dir_part) {
        glob_raw(dir_part)
            .into_iter()
            .filter(|d| is_dir(d))
            .collect()
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
