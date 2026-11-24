//! LOSAT's top level (`losat` and the words before a program's name, `top_level`) writes the bytes
//! frozen in `tests/fixtures/top_level` from the top level of LOSAT 0.2.0 before line RG item
//! K4i-3c (clap's top level then): its help, whose usage names the executable as `argv[0]` called
//! it, its version, and its errors, for every form that it took or refused.

use std::ffi::OsString;

use super::{top_level, TopLevel};

/// The fixture: one line per form,
/// `name<TAB>os<TAB>status<TAB>stream<TAB>text<TAB>argv0<TAB>words`; `os` is `-` for every
/// platform, `unix` or `windows`; the text is `@help NAME` for the file `help.txt` beside it with
/// its usage naming `NAME`, `@version` for `losat <version>` and a line end, or the text with `\n`
/// for its line ends; `argv0` is `(none)` for an empty argv; the words are separated by spaces;
/// `''` is an empty word and `\xff` a byte.
const FORMS: &str = include_str!("../../tests/fixtures/top_level/forms.tsv");

/// The text that the fixture names (`@help`): the top level's help, which `losat` alone writes to
/// the standard error, with its usage naming `losat`.
const HELP: &str = include_str!("../../tests/fixtures/top_level/help.txt");

/// The usage line of `HELP`.
const USAGE: &str = "USAGE\n  losat <COMMAND>\n";

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

/// Whether a row of the fixture for `os` is one of this platform's.
fn on_this_platform(os: &str) -> bool {
    match os {
        "-" => true,
        "unix" => cfg!(unix),
        "windows" => cfg!(windows),
        other => panic!("no platform {other}"),
    }
}

/// The help whose usage names `name`.
fn help(name: &str) -> String {
    assert!(HELP.starts_with(USAGE));
    HELP.replacen(USAGE, &format!("USAGE\n  {name} <COMMAND>\n"), 1)
}

/// The forms of the fixture on this platform: each one's argv and what the top level writes.
fn forms() -> Vec<(String, Vec<OsString>, TopLevel)> {
    FORMS
        .lines()
        .filter(|line| !line.starts_with('#'))
        .filter_map(|line| {
            let fields: Vec<&str> = line.split('\t').collect();
            assert_eq!(fields.len(), 7, "{line}");
            let (name, os, status, stream, text, argv0, words) = (
                fields[0], fields[1], fields[2], fields[3], fields[4], fields[5], fields[6],
            );
            if !on_this_platform(os) {
                return None;
            }
            let word = |text: &str| word(if text == "''" { "" } else { text });
            let argv0: Option<Vec<OsString>> = match argv0 {
                "(none)" => Some(Vec::new()),
                argv0 => word(argv0).map(|argv0| vec![argv0]),
            };
            let words: Option<Vec<OsString>> = match words.is_empty() {
                true => Some(Vec::new()),
                false => words.split(' ').map(word).collect(),
            };
            let text = match text.strip_prefix('@') {
                Some("version") => format!("losat {}\n", env!("CARGO_PKG_VERSION")),
                Some(other) => match other.strip_prefix("help ") {
                    Some(name) => help(name),
                    None => panic!("no text {other}"),
                },
                None => text.replace("\\n", "\n"),
            };
            let expected = TopLevel {
                text,
                error: stream == "stderr",
                status: status.parse().unwrap(),
            };
            let argv = [argv0?, words?].concat();
            Some((name.to_string(), argv, expected))
        })
        .collect()
}

#[test]
fn the_top_level_writes_its_frozen_bytes() {
    let forms = forms();
    assert!(forms.len() >= 90);
    // The fixture's helps name the executable as `losat`, a copy (`losat.exe`) and a link
    // (`mylosat`) were called.
    for name in ["losat", "losat.exe", "mylosat"] {
        let named = help(name);
        assert!(
            forms.iter().any(|(_, _, form)| form.text == named),
            "{name}"
        );
    }
    for (name, argv, expected) in forms {
        assert_eq!(top_level(&argv), expected, "{name} {argv:?}");
    }
}
