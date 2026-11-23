//! LOSAT's top level (`losat` and the words before a program's name, `top_level`) writes the bytes
//! frozen in `tests/fixtures/top_level` from the top level of LOSAT 0.2.0 before line RG item
//! K4i-3c (clap's): its help, its version, and its errors, for every form that it took or refused.

use std::ffi::OsString;

use super::{top_level, TopLevel};

/// The fixture: one line per form, `name<TAB>status<TAB>stream<TAB>text<TAB>words`; the text is
/// `@name` for the file `name.txt` beside it, `@version` for `losat <version>` and a line end, or
/// the text with `\n` for its line ends; the words are separated by spaces, `''` an empty word,
/// `\xff` a byte.
const FORMS: &str = include_str!("../../tests/fixtures/top_level/forms.tsv");

/// The text that the fixture names (`@help`): the top level's help, which `losat` alone writes to
/// the standard error.
const HELP: &str = include_str!("../../tests/fixtures/top_level/help.txt");

/// The word `text` of the fixture (`\xff` a byte that is not UTF-8).
#[cfg(unix)]
fn word(text: &str) -> Option<OsString> {
    use std::os::unix::ffi::OsStringExt;
    let parts: Vec<&[u8]> = text.split("\\xff").map(str::as_bytes).collect();
    Some(OsString::from_vec(parts.join(&0xff_u8)))
}

/// The word `text` of the fixture, `None` for one that is not UTF-8.
#[cfg(not(unix))]
fn word(text: &str) -> Option<OsString> {
    (!text.contains("\\xff")).then(|| OsString::from(text))
}

/// The forms of the fixture: each one's argv and what the top level writes.
fn forms() -> Vec<(String, Vec<OsString>, TopLevel)> {
    FORMS
        .lines()
        .filter(|line| !line.starts_with('#'))
        .filter_map(|line| {
            let fields: Vec<&str> = line.split('\t').collect();
            let (name, status, stream, text) = (fields[0], fields[1], fields[2], fields[3]);
            let words = fields.get(4).copied().unwrap_or("");
            let words: Option<Vec<OsString>> = match words.is_empty() {
                true => Some(Vec::new()),
                false => words
                    .split(' ')
                    .map(|text| word(if text == "''" { "" } else { text }))
                    .collect(),
            };
            let argv: Option<Vec<OsString>> = words.map(|words| {
                std::iter::once(OsString::from("losat"))
                    .chain(words)
                    .collect()
            });
            let text = match text.strip_prefix('@') {
                Some("version") => format!("losat {}\n", env!("CARGO_PKG_VERSION")),
                Some("help") => HELP.to_string(),
                Some(other) => panic!("no text {other}"),
                None => text.replace("\\n", "\n"),
            };
            let expected = TopLevel {
                text,
                error: stream == "stderr",
                status: status.parse().unwrap(),
            };
            argv.map(|argv| (name.to_string(), argv, expected))
        })
        .collect()
}

#[test]
fn the_top_level_writes_its_frozen_bytes() {
    let forms = forms();
    assert!(forms.len() >= 20);
    for (name, argv, expected) in forms {
        assert_eq!(top_level(&argv), expected, "{name} {argv:?}");
    }
}
