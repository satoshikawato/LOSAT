//! A minimal JSON writer for the ABI v2 responses and HSP records (docs/web/abi_v2.md
//! §8-9). The adapter writes only objects of strings, integers, numbers, booleans,
//! arrays and null, so no JSON library is needed.

use std::fmt::Write;

/// Appends `text` as a JSON string.
pub fn string(out: &mut String, text: &str) {
    out.push('"');
    for character in text.chars() {
        match character {
            '"' => out.push_str("\\\""),
            '\\' => out.push_str("\\\\"),
            '\n' => out.push_str("\\n"),
            '\r' => out.push_str("\\r"),
            '\t' => out.push_str("\\t"),
            c if (c as u32) < 0x20 => {
                let _ = write!(out, "\\u{:04x}", c as u32);
            }
            c => out.push(c),
        }
    }
    out.push('"');
}

/// Appends a number with the shortest text that reads back as the same `f64`, or null
/// when it is not finite (JSON has no infinities).
pub fn number(out: &mut String, value: f64) {
    if value.is_finite() {
        let _ = write!(out, "{value:?}");
    } else {
        out.push_str("null");
    }
}

/// Appends `[start, end]` or null.
pub fn range(out: &mut String, range: Option<(u64, u64)>) {
    match range {
        Some((start, end)) => {
            let _ = write!(out, "[{start},{end}]");
        }
        None => out.push_str("null"),
    }
}

/// Appends an optional integer or null.
pub fn optional<T: std::fmt::Display>(out: &mut String, value: Option<T>) {
    match value {
        Some(value) => {
            let _ = write!(out, "{value}");
        }
        None => out.push_str("null"),
    }
}

/// Appends an optional string or null.
pub fn optional_string(out: &mut String, value: Option<&str>) {
    match value {
        Some(value) => string(out, value),
        None => out.push_str("null"),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn strings_and_numbers_are_valid_json() {
        let mut out = String::new();
        string(&mut out, "a\"b\\c\n\u{1}é");
        assert_eq!(out, "\"a\\\"b\\\\c\\n\\u0001é\"");
        let mut out = String::new();
        number(&mut out, 1e-100);
        out.push(' ');
        number(&mut out, 0.5);
        out.push(' ');
        number(&mut out, f64::INFINITY);
        assert_eq!(out, "1e-100 0.5 null");
    }
}
