use super::*;
use std::path::{Path, PathBuf};

/// The samples that `LOSAT/tests/losat_wording.py` checks too: NCBI's text, then LOSAT's.
pub(crate) const SAMPLES: &str = include_str!("../tests/fixtures/wording/samples.tsv");

/// The samples of [`losat_output_bytes`] alone (`LOSAT/tests/losat_wording.py` checks them too).
const OUTPUT_SAMPLES: &str = include_str!("../tests/fixtures/wording/output_samples.tsv");

/// A field of the samples: `\\n` and `\\t` are a line feed and a tab, `\\xNN` the byte NN
/// (below 0x80), `{version}` is the crate's version, `{verlen}` the byte of the length of
/// `LOSAT <version>` and `{at<n>}` the number n + [`IDENT_GROWTH`].
fn unescape(field: &str) -> String {
    let field = field
        .replace("\\n", "\n")
        .replace("\\t", "\t")
        .replace("{version}", env!("CARGO_PKG_VERSION"))
        .replace(
            "{verlen}",
            &char::from(LOSAT_VERSION.len() as u8).to_string(),
        );
    let mut field = field;
    while let Some(at) = field.find("{at") {
        let end = at + field[at..].find('}').expect("{at<n>}");
        let n: usize = field[at + 3..end].parse().expect("{at<n>}");
        field.replace_range(at..=end, &(n + IDENT_GROWTH).to_string());
    }
    let mut out = String::with_capacity(field.len());
    let mut rest = field.as_str();
    while let Some(at) = rest.find("\\x") {
        out.push_str(&rest[..at]);
        let byte = u8::from_str_radix(&rest[at + 2..at + 4], 16).expect("\\xNN");
        out.push(char::from(byte));
        rest = &rest[at + 4..];
    }
    out.push_str(rest);
    out
}

/// The `.rs` files below `dir`.
pub(crate) fn rust_files(dir: &Path, out: &mut Vec<PathBuf>) {
    for entry in std::fs::read_dir(dir).unwrap() {
        let path = entry.unwrap().path();
        if path.is_dir() {
            rust_files(&path, out);
        } else if path.extension().is_some_and(|ext| ext == "rs") {
            out.push(path);
        }
    }
}

/// The lines of `text` that are code outside the tests: comment lines and the top-level
/// `#[cfg(test)] mod name { ... }` blocks are left out (rustfmt puts the closing brace of a
/// top-level module at the start of a line).
pub(crate) fn code_lines(text: &str) -> Vec<(usize, &str)> {
    let mut lines = Vec::new();
    let mut in_tests = false;
    let mut after_cfg_test = false;
    for (number, line) in text.lines().enumerate() {
        if in_tests {
            in_tests = line != "}";
            continue;
        }
        if line == "#[cfg(test)]" {
            after_cfg_test = true;
            continue;
        }
        if after_cfg_test {
            if line.starts_with("#[") {
                continue;
            }
            after_cfg_test = false;
            if line.starts_with("mod ") && line.ends_with('{') {
                in_tests = true;
                continue;
            }
        }
        // A top-level test module under another attribute (`#[cfg(all(test, ...))]`).
        if line == "mod tests {" {
            in_tests = true;
            continue;
        }
        if !line.trim_start().starts_with("//") {
            lines.push((number + 1, line));
        }
    }
    lines
}

#[test]
fn code_lines_leave_out_comments_and_test_modules() {
    let text =
        "fn a() {}\n// \"NCBI C++ Exception:\\n\"\n#[cfg(test)]\n#[path = \"x.rs\"]\nmod x;\n\
                    fn b() {}\n#[cfg(test)]\nmod tests {\n    fn c() {}\n}\nfn d() {}\n";
    assert_eq!(
        code_lines(text),
        [
            (1, "fn a() {}"),
            (5, "mod x;"),
            (6, "fn b() {}"),
            (11, "fn d() {}")
        ]
    );
}

/// Whether `path` (relative to the repository) is a test file, which may hold NCBI's words as
/// the expected values of NCBI's side.
fn test_file(path: &Path) -> bool {
    let name = path.file_name().unwrap().to_string_lossy();
    let stem = path.file_stem().unwrap().to_string_lossy();
    stem == "tests"
        || stem.ends_with("_tests")
        || stem.starts_with("test_")
        || stem.starts_with("tests_")
        || name.contains(".test.")
        || name.contains(".spec.")
        || path.components().any(|part| {
            let part = part.as_os_str().to_string_lossy();
            part == "tests" || part.ends_with("_tests")
        })
}

/// The files of `dir` with one of `extensions`, below it (`node_modules` and build outputs
/// left out).
fn files(dir: &Path, extensions: &[&str], out: &mut Vec<PathBuf>) {
    let Ok(entries) = std::fs::read_dir(dir) else {
        return;
    };
    for entry in entries {
        let path = entry.unwrap().path();
        let name = path.file_name().unwrap().to_string_lossy().into_owned();
        if path.is_dir() {
            if !matches!(
                name.as_str(),
                "node_modules" | "dist" | "target" | "__pycache__"
            ) {
                files(&path, extensions, out);
            }
        } else if path
            .extension()
            .is_some_and(|ext| extensions.iter().any(|e| ext == *e))
        {
            out.push(path);
        }
    }
}

/// Whether `line` names NCBI's exception classes in a tag, `(<class>::<code>)`, as NCBI's
/// reports do: a tag written by hand outside this module (the reports of
/// `ncbi_exception::NcbiException` and its posts pass through [`losat_bytes`]).
fn has_tag(line: &str) -> bool {
    line.match_indices('(').any(|(at, _)| {
        let rest = &line[at + 1..];
        let class: String = rest
            .chars()
            .take_while(|c| c.is_ascii_alphanumeric() || *c == '_')
            .collect();
        let code = rest[class.len()..].strip_prefix("::").unwrap_or("(");
        (class.ends_with("Exception") || SERIAL_EXCEPTION_CLASSES.contains(&class.as_str()))
            && code
                .find(')')
                .is_some_and(|close| !code[..close].contains('('))
    })
}

// NCBI's words in a user-facing string outside this module fail here (lines NW and NW2, the
// Owner's 2026-11-20 and 2026-11-21 instructions): NCBI's user-manual footer, NCBI's
// exception header, the name of NCBI's C++ Toolkit, a hand-written exception tag in the
// engine's Rust code, BLAST+ named as the program that runs (in LOSAT Web's texts, where
// `NCBI BLAST+` names NCBI's), NCBI's version `2.17.0+` stated as the program's (outside
// NCBI's package names `ncbi-blast-2.17.0+-...`) and NCBI's `BLAST ... error:` prefixes.
// Comment lines, test files and test modules are not read; `ncbi_exception.rs` builds NCBI's
// report and maps it (its own test keeps it the one owner of that text).
#[test]
fn ncbi_words_are_written_only_here() {
    let repo = Path::new(env!("CARGO_MANIFEST_DIR")).parent().unwrap();
    let mut paths = Vec::new();
    for dir in ["LOSAT/src", "web/adapter/src", "python/src"] {
        rust_files(&repo.join(dir), &mut paths);
    }
    let rust = paths.len();
    assert!(rust > 100, "{rust} Rust files");
    for dir in ["web/app/src", "web/packages"] {
        files(&repo.join(dir), &["ts", "vue", "js"], &mut paths);
    }
    files(&repo.join("python/pylosat"), &["py"], &mut paths);
    assert!(paths.len() > rust + 20, "{} files", paths.len());
    // NCBI's package names (`ncbi-blast-2.17.0+-x64-linux.tar.gz`) are not version statements.
    let package = format!("ncbi-blast-{NCBI_VERSION}");
    let mut offenders = Vec::new();
    for path in paths {
        let relative = path.strip_prefix(repo).unwrap();
        if test_file(relative)
            || relative == Path::new("LOSAT/src/wording.rs")
            || relative == Path::new("LOSAT/src/ncbi_exception.rs")
        {
            continue;
        }
        let text = std::fs::read_to_string(&path).unwrap();
        let is_rust = relative.extension().is_some_and(|ext| ext == "rs");
        let lines: Vec<(usize, &str)> = if is_rust {
            code_lines(&text)
        } else {
            text.lines()
                .enumerate()
                .map(|(number, line)| (number + 1, line))
                .filter(|(_, line)| {
                    let line = line.trim_start();
                    !["//", "/*", "*", "#", "<!--"]
                        .iter()
                        .any(|comment| line.starts_with(comment))
                })
                .collect()
        };
        for (number, line) in lines {
            let bare_blast_plus =
                !is_rust && line.matches("BLAST+").count() > line.matches("NCBI BLAST+").count();
            if [
                NCBI_USER_MANUAL_FOOTER,
                "NCBI C++ Exception",
                "NCBI C++ Toolkit",
            ]
            .iter()
            .any(|needle| line.contains(needle))
                || (is_rust && has_tag(line))
                || bare_blast_plus
                || line.matches(NCBI_VERSION).count() > line.matches(&package).count()
                || ERROR_PREFIXES.iter().any(|(ncbi, _)| line.contains(ncbi))
            {
                offenders.push(format!("{}:{number}: {}", relative.display(), line.trim()));
            }
        }
    }
    assert!(
        offenders.is_empty(),
        "write LOSAT's words with crate::wording (PD-LOSAT-CLI-NONSEARCH-DIFFERENCES \
             exceptions 12 and 14):\n{}",
        offenders.join("\n")
    );
}

/// `text` with each `{db_dir}` token written as the path of `tests/db_search/rps/db_dir_len.txt`
/// characters that it stands for in the frozen rows.
fn expand_db_dir(text: &[u8]) -> Vec<u8> {
    const TOKEN: &[u8] = b"{db_dir}";
    let length: usize = include_str!("../tests/db_search/rps/db_dir_len.txt")
        .trim()
        .parse()
        .unwrap();
    let mut out = Vec::with_capacity(text.len() + length);
    let mut i = 0;
    while i < text.len() {
        if text[i..].starts_with(TOKEN) {
            out.resize(out.len() + length, b'd');
            i += TOKEN.len();
        } else {
            out.push(text[i]);
            i += 1;
        }
    }
    out
}

// `asn_text::visible_string_line` writes each message of the frozen archives as NCBI wrote
// it: the wraps of a message mapped by `losat_output_bytes` are the writer's.
#[test]
fn archive_messages_are_wrapped_again_as_ncbi_wraps_them() {
    let fixtures = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures");
    let mut paths = Vec::new();
    files(&fixtures, &["out", "asn", "txt"], &mut paths);
    let mut count = 0;
    for path in paths {
        let data = std::fs::read(&path).unwrap();
        let mut search = 0;
        while let Some(found) = find(&data[search..], ARCHIVE_MESSAGE) {
            let at = search + found;
            let line_start = data[..at]
                .iter()
                .rposition(|&b| b == b'\n')
                .map_or(0, |i| i + 1);
            let quote = at + ARCHIVE_MESSAGE.len() - 1;
            search = quote + 1;
            if !data[line_start..at].iter().all(|&b| b == b' ') {
                continue;
            }
            let Some((value, end)) = asn_string(&data, quote) else {
                continue;
            };
            search = end;
            // A message with a byte outside 0x20..0x7E was not written by an archive writer
            // (NCBI's and LOSAT's write `#` there): a hand-made input of the rpsbproc fixtures
            // (`root/vis_all.asn`), which says nothing of the wraps.
            if !value.iter().all(|b| (b' '..=b'~').contains(b)) {
                continue;
            }
            // A frozen `{db_dir}` stands for a path of the length the rows were frozen with
            // (`db_search_fixtures.py` `inflate_rps`): the wraps of NCBI's bytes are those of
            // the whole path.
            assert_eq!(
                crate::report::serial::asn_text::visible_string_line(
                    &data[line_start..quote],
                    &expand_db_dir(&value)
                ),
                expand_db_dir(&data[line_start..end]),
                "{} byte {at}",
                path.display()
            );
            count += 1;
        }
    }
    assert!(count > 300, "{count} messages");
}

#[test]
fn the_samples_map_to_losats_text() {
    let mut count = 0;
    for line in SAMPLES.lines().filter(|line| !line.starts_with('#')) {
        let (ncbi, losat) = line.split_once('\t').expect("two fields");
        let (ncbi, losat) = (unescape(ncbi), unescape(losat));
        assert_eq!(losat_text(&ncbi), losat, "{ncbi:?}");
        // Mapped text maps to itself.
        assert_eq!(losat_text(&losat), losat, "{losat:?}");
        count += 1;
    }
    assert!(count >= 10, "{count} samples");
}

#[test]
fn the_output_samples_map_to_losats_bytes() {
    let mut count = 0;
    for line in OUTPUT_SAMPLES.lines().filter(|line| !line.starts_with('#')) {
        let (ncbi, losat) = line.split_once('\t').expect("two fields");
        // The fields are read as Latin-1 (`\xNN` any byte), as `losat_wording.py` reads them.
        let latin1 = |field: &str| -> Vec<u8> {
            unescape(field)
                .chars()
                .map(|c| u8::try_from(c).expect("a Latin-1 sample"))
                .collect()
        };
        let (ncbi, losat) = (latin1(ncbi), latin1(losat));
        assert_eq!(
            losat_output_bytes(&ncbi).as_ref(),
            losat.as_slice(),
            "{:?}",
            String::from_utf8_lossy(&ncbi)
        );
        count += 1;
    }
    assert!(count >= 7, "{count} samples");
}

#[test]
fn text_without_ncbis_words_is_borrowed() {
    for text in [
        "LOSAT query/options error: No hits are being saved\n",
        "Warning: (1.1) [blastn] x::y\n",
        "ncbi-blast-2.17.0+-x64-linux.tar.gz, NCBI BLAST+ 2.17.0 aborts",
        "",
    ] {
        assert!(matches!(losat_text(text), Cow::Borrowed(_)), "{text:?}");
    }
}

// A collapsed message (an archive's) maps to the collapsed mapping of the message.
#[test]
fn the_mapping_of_a_collapsed_message_is_the_collapsed_mapping() {
    let collapse = |text: &str| text.split_whitespace().collect::<Vec<_>>().join(" ");
    for line in SAMPLES.lines().filter(|line| !line.starts_with('#')) {
        let ncbi = unescape(line.split_once('\t').unwrap().0);
        assert_eq!(
            collapse(&losat_text(&collapse(&ncbi))),
            collapse(&losat_text(&ncbi)),
            "{ncbi:?}"
        );
    }
}
