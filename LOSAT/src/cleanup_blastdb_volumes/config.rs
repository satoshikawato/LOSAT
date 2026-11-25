//! Reader for `[BLAST] BLASTDB` in `.ncbirc` / `ncbi.ini`, following Python's
//! `configparser.ConfigParser` semantics used by `cleanup-blastdb-volumes.py`.
//!
//! NCBI reference (598d8ae6): c++/src/app/blast/cleanup-blastdb-volumes.py:96-101

use std::collections::BTreeMap;
use std::path::Path;

/// Whether `ch` is whitespace according to Python's `str.isspace()` (Unicode `White_Space` plus
/// ASCII information separators U+001C–U+001F).
#[inline]
pub(crate) fn is_py_space(ch: char) -> bool {
    ch.is_whitespace() || matches!(ch, '\u{001c}'..='\u{001f}')
}

/// Trims leading and trailing Python whitespace (`str.strip()`).
#[inline]
pub(crate) fn py_trim(s: &str) -> &str {
    s.trim_matches(is_py_space)
}

/// Trims trailing Python whitespace (`str.rstrip()`).
#[inline]
pub(crate) fn py_trim_end(s: &str) -> &str {
    s.trim_end_matches(is_py_space)
}

/// Splits `s` on runs of Python whitespace (`str.split()`).
#[inline]
pub(crate) fn py_split_whitespace(s: &str) -> impl Iterator<Item = &str> {
    s.split(is_py_space).filter(|part| !part.is_empty())
}

/// Splits `text` on `\r\n`, `\n`, or standalone `\r` (Python universal newlines).
pub(crate) fn universal_lines(text: &str) -> Vec<&str> {
    let mut out = Vec::new();
    let bytes = text.as_bytes();
    let mut start = 0usize;
    let mut i = 0usize;
    while i < bytes.len() {
        if bytes[i] == b'\r' {
            out.push(&text[start..i]);
            i += if bytes.get(i + 1) == Some(&b'\n') {
                2
            } else {
                1
            };
            start = i;
        } else if bytes[i] == b'\n' {
            out.push(&text[start..i]);
            i += 1;
            start = i;
        } else {
            i += 1;
        }
    }
    if start < bytes.len() {
        out.push(&text[start..]);
    }
    out
}

fn bad_interpolation(raw: &str) -> String {
    format!("invalid interpolation syntax in '{raw}'")
}

/// Returns the mutable option map for `sec_name` (`DEFAULT` or a named section).
fn section_map_mut<'a>(
    sec_name: &str,
    defaults: &'a mut BTreeMap<String, String>,
    sections: &'a mut BTreeMap<String, BTreeMap<String, String>>,
) -> &'a mut BTreeMap<String, String> {
    if sec_name == "DEFAULT" {
        defaults
    } else {
        sections.get_mut(sec_name).expect("section inserted")
    }
}

/// Expands Python `BasicInterpolation` (`%%` -> `%`, `%(key)s` -> option value) up to 10 levels,
/// returning `Err` on malformed `%` sequences or missing keys.
fn interpolate(
    raw: &str,
    section: &BTreeMap<String, String>,
    defaults: &BTreeMap<String, String>,
    depth: usize,
) -> Result<String, String> {
    if depth > 10 {
        return Err("interpolation depth exceeded".to_string());
    }
    let mut out = String::with_capacity(raw.len());
    let mut rest = raw;
    while let Some(pos) = rest.find('%') {
        out.push_str(&rest[..pos]);
        let after = &rest[pos + 1..];
        if let Some(tail) = after.strip_prefix('%') {
            out.push('%');
            rest = tail;
        } else if let Some(after_lparen) = after.strip_prefix('(') {
            let Some(close) = after_lparen.find(')') else {
                return Err(bad_interpolation(raw));
            };
            if close == 0 || !after_lparen[close..].starts_with(")s") {
                return Err(bad_interpolation(raw));
            }
            let key = after_lparen[..close].to_lowercase();
            let val = section
                .get(&key)
                .or_else(|| defaults.get(&key))
                .ok_or_else(|| format!("missing interpolation option '{key}'"))?;
            if val.contains('%') {
                let expanded = interpolate(val, section, defaults, depth + 1)?;
                out.push_str(&expanded);
            } else {
                out.push_str(val);
            }
            rest = &after_lparen[close + 2..];
        } else {
            return Err(bad_interpolation(raw));
        }
    }
    out.push_str(rest);
    Ok(out)
}

/// Parses INI `text` with strict Python ConfigParser rules and returns `[BLAST] BLASTDB`
/// if present.
pub(crate) fn parse_blastdb_from_ini(
    text: &str,
    source_name: &str,
) -> Result<Option<Vec<u8>>, String> {
    let mut defaults: BTreeMap<String, String> = BTreeMap::new();
    let mut sections: BTreeMap<String, BTreeMap<String, String>> = BTreeMap::new();
    let mut current_section: Option<String> = None;
    let mut last_option: Option<String> = None;
    let mut last_indent: Option<usize> = None;
    let mut pending_empty_lines = 0usize;

    for (idx, line) in universal_lines(text).into_iter().enumerate() {
        let lineno = idx + 1;
        let trimmed = py_trim(line);
        if trimmed.is_empty() {
            if current_section.is_some() && last_option.is_some() {
                pending_empty_lines += 1;
            }
            continue;
        }
        if trimmed.starts_with('#') || trimmed.starts_with(';') {
            continue;
        }

        // Python ConfigParser compares leading whitespace counts in Unicode code points.
        let indent = line.chars().take_while(|&ch| is_py_space(ch)).count();
        if indent > 0
            && last_indent.is_some_and(|prev| indent > prev)
            && current_section.is_some()
            && last_option.is_some()
        {
            let sec_name = current_section.as_deref().unwrap();
            let opt_name = last_option.as_deref().unwrap();
            let target = section_map_mut(sec_name, &mut defaults, &mut sections);
            if let Some(val) = target.get_mut(opt_name) {
                for _ in 0..=pending_empty_lines {
                    val.push('\n');
                }
                val.push_str(trimmed);
            }
            pending_empty_lines = 0;
            continue;
        }

        pending_empty_lines = 0;

        if let Some(rest) = trimmed.strip_prefix('[') {
            if let Some(end) = rest.rfind(']') {
                if end > 0 {
                    let inside = &rest[..end];
                    if inside != "DEFAULT" && sections.contains_key(inside) {
                        return Err(format!(
                            "{source_name}:{lineno}: duplicate section '{inside}'"
                        ));
                    }
                    if inside != "DEFAULT" {
                        sections.insert(inside.to_string(), BTreeMap::new());
                    }
                    current_section = Some(inside.to_string());
                    last_option = None;
                    last_indent = Some(indent);
                    continue;
                }
            }
        }

        let Some(sec_name) = current_section.as_deref() else {
            return Err(format!(
                "{source_name}:{lineno}: file contains no section headers"
            ));
        };

        let sep_pos = line
            .bytes()
            .position(|b| b == b'=' || b == b':')
            .ok_or_else(|| format!("{source_name}:{lineno}: source contains parsing errors"))?;
        let key = py_trim(&line[..sep_pos]).to_lowercase();
        if key.is_empty() {
            return Err(format!("{source_name}:{lineno}: empty option name"));
        }
        let val = py_trim(&line[sep_pos + 1..]).to_string();
        let map = section_map_mut(sec_name, &mut defaults, &mut sections);
        if map.contains_key(&key) {
            return Err(format!(
                "{source_name}:{lineno}: duplicate option '{key}' in section '{sec_name}'"
            ));
        }
        map.insert(key.clone(), val);
        last_option = Some(key);
        last_indent = Some(indent);
    }

    // Strip trailing blank-line padding in values (Python ConfigParser._join_multiline_values
    // preserves internal blank lines and strips trailing whitespace on final join).
    let Some(blast_sec) = sections.get("BLAST") else {
        return Ok(None);
    };
    let Some(raw_val) = blast_sec.get("blastdb").or_else(|| defaults.get("blastdb")) else {
        return Ok(None);
    };
    let value = interpolate(py_trim_end(raw_val), blast_sec, &defaults, 1)?;
    Ok(Some(value.into_bytes()))
}

// NCBI reference (598d8ae6): c++/src/app/blast/cleanup-blastdb-volumes.py:96-101
/// Reads `BLASTDB` from the `[BLAST]` section of `config_file`, or `Ok(None)` when the file
/// cannot be opened or does not set `[BLAST] BLASTDB`.
pub(crate) fn get_blastdb_from_ncbi_config(config_file: &Path) -> Result<Option<Vec<u8>>, String> {
    // Python's `ConfigParser.read` ignores unreadable / missing files (OSError) and returns empty.
    let Ok(bytes) = std::fs::read(config_file) else {
        return Ok(None);
    };
    let text = String::from_utf8(bytes)
        .map_err(|e| format!("{}: invalid UTF-8: {e}", config_file.display()))?;
    parse_blastdb_from_ini(&text, &config_file.display().to_string())
}

#[cfg(test)]
#[path = "config_tests.rs"]
mod tests;
