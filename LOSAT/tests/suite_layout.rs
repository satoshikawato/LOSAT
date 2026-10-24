//! Checks how the files of `LOSAT/tests/` are built into test binaries. `Cargo.toml` turns
//! automatic test discovery off and lists eight binaries: six suites (`tests/suites/<name>.rs`,
//! one `#[path]` module per file) and two files that change the process environment and so stay
//! alone. A file that no binary lists would silently not build, a file in two binaries would run
//! twice, and a suite member that changes the environment would disturb the other tests of its
//! process. Each test below names every offender and says how to fix it. Standard library only.

use std::collections::BTreeMap;
use std::fs;
use std::path::{Path, PathBuf};

/// The programs, each also the name of a suite.
const PROGRAMS: [&str; 5] = ["tblastx", "tblastn", "blastx", "blastp", "blastn"];

/// The calls that change the process environment.
const ENV_CALLS: [&str; 3] = ["set_var", "remove_var", "set_current_dir"];

/// The two files that change the process environment and therefore stay alone. The walk of their
/// module trees must find an environment change (`walk_finds_the_lone_files_changes`); when it
/// stops finding one, the file can join a suite.
const LONE_FILES: [(&str, &str); 2] = [
    (
        "blastn_report_stream_api",
        "tests/blastn_report_stream_api.rs",
    ),
    ("unit_tests", "tests/unit_tests.rs"),
];

fn manifest_dir() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
}

fn read(path: &Path) -> String {
    fs::read_to_string(path)
        .unwrap_or_else(|error| panic!("cannot read {}: {error}", path.display()))
        .replace('\r', "")
}

/// A path as text with `/` separators, relative to the crate directory when it is below it.
fn shown(path: &Path) -> String {
    let base = manifest_dir();
    path.strip_prefix(&base)
        .unwrap_or(path)
        .to_string_lossy()
        .replace('\\', "/")
}

// ---------------------------------------------------------------------------------------------
// Cargo.toml
// ---------------------------------------------------------------------------------------------

struct Target {
    name: String,
    path: String,
}

fn quoted_value(line: &str, key: &str) -> Option<String> {
    let rest = line.strip_prefix(key)?.trim_start().strip_prefix('=')?;
    let rest = rest.trim().strip_prefix('"')?;
    Some(rest[..rest.find('"')?].to_string())
}

/// The `[[test]]` blocks of `Cargo.toml`, read by lines.
fn cargo_targets() -> (bool, Vec<Target>) {
    let text = read(&manifest_dir().join("Cargo.toml"));
    let mut autotests_off = false;
    let mut targets: Vec<Target> = Vec::new();
    let mut in_test = false;
    for line in text.lines().map(str::trim) {
        if line.starts_with('[') {
            in_test = line == "[[test]]";
            if in_test {
                targets.push(Target {
                    name: String::new(),
                    path: String::new(),
                });
            }
            continue;
        }
        if line.replace(' ', "") == "autotests=false" {
            autotests_off = true;
        }
        if in_test {
            let target = targets.last_mut().expect("a [[test]] block");
            if let Some(name) = quoted_value(line, "name") {
                target.name = name;
            }
            if let Some(path) = quoted_value(line, "path") {
                target.path = path;
            }
        }
    }
    (autotests_off, targets)
}

fn top_level_test_files() -> Vec<String> {
    let mut files: Vec<String> = fs::read_dir(manifest_dir().join("tests"))
        .expect("read tests/")
        .map(|entry| entry.expect("directory entry").path())
        .filter(|path| path.is_file() && path.extension().is_some_and(|ext| ext == "rs"))
        .map(|path| path.file_name().unwrap().to_string_lossy().into_owned())
        .collect();
    files.sort();
    files
}

fn suite_roots() -> Vec<String> {
    let mut roots: Vec<String> = fs::read_dir(manifest_dir().join("tests/suites"))
        .expect("read tests/suites/")
        .map(|entry| entry.expect("directory entry").path())
        .filter(|path| path.extension().is_some_and(|ext| ext == "rs"))
        .map(|path| {
            format!(
                "tests/suites/{}",
                path.file_name().unwrap().to_string_lossy()
            )
        })
        .collect();
    roots.sort();
    roots
}

/// The suite a new file belongs to, from its name.
fn suite_for(stem: &str) -> &'static str {
    let program_word = stem.strip_prefix("run_local_").unwrap_or(stem);
    for program in PROGRAMS {
        if program_word.starts_with(program) {
            return program;
        }
    }
    if stem.starts_with("megablast") {
        "blastn"
    } else {
        "common"
    }
}

// ---------------------------------------------------------------------------------------------
// A small Rust tokenizer: identifiers, punctuation and string contents; comments, character
// literals and the insides of raw strings never reach the checks.
// ---------------------------------------------------------------------------------------------

#[derive(Debug, Clone, PartialEq)]
enum Token {
    Ident(String),
    Str(String),
    Punct(char),
}

fn tokenize(source: &str) -> Vec<Token> {
    let chars: Vec<char> = source.chars().collect();
    let mut tokens = Vec::new();
    let mut i = 0;
    while i < chars.len() {
        let c = chars[i];
        if c.is_whitespace() {
            i += 1;
        } else if c == '/' && chars.get(i + 1) == Some(&'/') {
            while i < chars.len() && chars[i] != '\n' {
                i += 1;
            }
        } else if c == '/' && chars.get(i + 1) == Some(&'*') {
            let mut depth = 0;
            while i < chars.len() {
                if chars[i] == '/' && chars.get(i + 1) == Some(&'*') {
                    depth += 1;
                    i += 2;
                } else if chars[i] == '*' && chars.get(i + 1) == Some(&'/') {
                    depth -= 1;
                    i += 2;
                    if depth == 0 {
                        break;
                    }
                } else {
                    i += 1;
                }
            }
        } else if c == '"' {
            let (text, next) = read_string(&chars, i + 1);
            tokens.push(Token::Str(text));
            i = next;
        } else if c == '\'' {
            i = skip_char_or_lifetime(&chars, i);
        } else if c.is_alphanumeric() || c == '_' {
            let start = i;
            while i < chars.len() && (chars[i].is_alphanumeric() || chars[i] == '_') {
                i += 1;
            }
            let word: String = chars[start..i].iter().collect();
            match (word.as_str(), chars.get(i)) {
                ("r" | "br" | "cr", Some('#' | '"')) => {
                    let mut j = i;
                    while chars.get(j) == Some(&'#') {
                        j += 1;
                    }
                    if chars.get(j) == Some(&'"') {
                        i = skip_raw_string(&chars, i, j - i);
                        tokens.push(Token::Str(String::new()));
                    } else {
                        tokens.push(Token::Ident(word));
                    }
                }
                ("b" | "c", Some('"')) => {
                    let (text, next) = read_string(&chars, i + 1);
                    tokens.push(Token::Str(text));
                    i = next;
                }
                ("b", Some('\'')) => {
                    i = skip_char_or_lifetime(&chars, i);
                }
                _ => tokens.push(Token::Ident(word)),
            }
        } else {
            tokens.push(Token::Punct(c));
            i += 1;
        }
    }
    tokens
}

/// Reads a string whose opening quote is before `start`; returns the text and the index after the
/// closing quote.
fn read_string(chars: &[char], start: usize) -> (String, usize) {
    let mut text = String::new();
    let mut i = start;
    while i < chars.len() {
        match chars[i] {
            '\\' => {
                text.push('\\');
                if let Some(&next) = chars.get(i + 1) {
                    text.push(next);
                }
                i += 2;
            }
            '"' => return (text, i + 1),
            other => {
                text.push(other);
                i += 1;
            }
        }
    }
    (text, i)
}

/// `hashes` hash marks start at `from`, then the opening quote; returns the index after the closing
/// quote and its hash marks.
fn skip_raw_string(chars: &[char], from: usize, hashes: usize) -> usize {
    let mut i = from + hashes + 1;
    while i < chars.len() {
        if chars[i] == '"' && (0..hashes).all(|k| chars.get(i + 1 + k) == Some(&'#')) {
            return i + 1 + hashes;
        }
        i += 1;
    }
    i
}

/// At a `'`: skips a character literal (`'{'`, `'"'`, `'\n'`, `'\''`); a lifetime leaves only the
/// quote behind.
fn skip_char_or_lifetime(chars: &[char], at: usize) -> usize {
    if chars.get(at + 1) == Some(&'\\') {
        let mut i = at + 2;
        // The escaped character itself, then up to the closing quote (`'\u{1F600}'` has more).
        i += 1;
        while i < chars.len() && chars[i] != '\'' && chars[i] != '\n' {
            i += 1;
        }
        i + 1
    } else if chars.get(at + 2) == Some(&'\'') {
        at + 3
    } else {
        at + 1
    }
}

// ---------------------------------------------------------------------------------------------
// The module tree of a test file
// ---------------------------------------------------------------------------------------------

#[derive(Default)]
struct Walk {
    visited: Vec<PathBuf>,
    /// (file, call) for every environment change found.
    changes: Vec<(PathBuf, String)>,
    /// Declarations that resolve to no file or to two.
    unresolved: Vec<String>,
}

/// Walks `file` and every file it pulls in with `mod x;`. `own_dir` says whether the file looks for
/// its child modules in its own directory (a root file, a `mod.rs`, a `#[path]` file) or in the
/// directory named after it (`foo.rs` looks in `foo/`).
fn walk_module_tree(file: &Path, own_dir: bool, walk: &mut Walk) {
    if walk.visited.iter().any(|seen| seen == file) {
        return;
    }
    walk.visited.push(file.to_path_buf());
    let dir = file.parent().expect("a file has a directory").to_path_buf();
    let children_dir = if own_dir {
        dir.clone()
    } else {
        dir.join(file.file_stem().expect("a file stem"))
    };
    let tokens = tokenize(&read(file));

    // Open braces: Some(name) for `mod name {`, None for any other.
    let mut braces: Vec<Option<String>> = Vec::new();
    let mut pending_path: Option<String> = None;
    let mut i = 0;
    while i < tokens.len() {
        match &tokens[i] {
            Token::Punct('#') => {
                // An attribute: `#[...]` or `#![...]`; remember `path = "..."`.
                let mut j = i + 1;
                if tokens.get(j) == Some(&Token::Punct('!')) {
                    j += 1;
                }
                if tokens.get(j) == Some(&Token::Punct('[')) {
                    let mut depth = 0;
                    let start = j;
                    while j < tokens.len() {
                        match tokens[j] {
                            Token::Punct('[') => depth += 1,
                            Token::Punct(']') => {
                                depth -= 1;
                                if depth == 0 {
                                    break;
                                }
                            }
                            _ => {}
                        }
                        j += 1;
                    }
                    if let (
                        Some(Token::Ident(key)),
                        Some(Token::Punct('=')),
                        Some(Token::Str(value)),
                    ) = (
                        tokens.get(start + 1),
                        tokens.get(start + 2),
                        tokens.get(start + 3),
                    ) {
                        if key == "path" {
                            pending_path = Some(value.clone());
                        }
                    }
                    i = j + 1;
                    continue;
                }
                i += 1;
            }
            Token::Ident(word) if word == "mod" => {
                let name = match tokens.get(i + 1) {
                    Some(Token::Ident(name)) => name.clone(),
                    _ => {
                        i += 1;
                        continue;
                    }
                };
                match tokens.get(i + 2) {
                    Some(Token::Punct('{')) => {
                        braces.push(Some(name));
                        pending_path = None;
                        i += 3;
                    }
                    Some(Token::Punct(';')) => {
                        let inline: Vec<&str> =
                            braces.iter().flatten().map(String::as_str).collect();
                        let mut here = children_dir.clone();
                        for part in &inline {
                            here = here.join(part);
                        }
                        let path_base = if inline.is_empty() {
                            dir.clone()
                        } else {
                            here.clone()
                        };
                        match pending_path.take() {
                            Some(relative) => {
                                let target = path_base.join(relative);
                                if target.is_file() {
                                    walk_module_tree(&target, true, walk);
                                } else {
                                    walk.unresolved.push(format!(
                                        "`mod {name};` with #[path] in {} points to {}, which does not exist",
                                        shown(file),
                                        shown(&target)
                                    ));
                                }
                            }
                            None => {
                                let plain = here.join(format!("{name}.rs"));
                                let nested = here.join(&name).join("mod.rs");
                                match (plain.is_file(), nested.is_file()) {
                                    (true, false) => walk_module_tree(&plain, false, walk),
                                    (false, true) => walk_module_tree(&nested, true, walk),
                                    (false, false) => walk.unresolved.push(format!(
                                        "`mod {name};` in {} resolves to neither {} nor {}",
                                        shown(file),
                                        shown(&plain),
                                        shown(&nested)
                                    )),
                                    (true, true) => walk.unresolved.push(format!(
                                        "`mod {name};` in {} resolves to both {} and {}",
                                        shown(file),
                                        shown(&plain),
                                        shown(&nested)
                                    )),
                                }
                            }
                        }
                        i += 3;
                    }
                    _ => i += 1,
                }
            }
            Token::Ident(word) if ENV_CALLS.contains(&word.as_str()) => {
                if tokens.get(i + 1) == Some(&Token::Punct('(')) {
                    walk.changes.push((file.to_path_buf(), format!("{word}(")));
                }
                i += 1;
            }
            Token::Punct('{') => {
                braces.push(None);
                pending_path = None;
                i += 1;
            }
            Token::Punct('}') => {
                braces.pop();
                pending_path = None;
                i += 1;
            }
            Token::Punct(';') => {
                pending_path = None;
                i += 1;
            }
            _ => i += 1,
        }
    }
}

fn walk_from(root_file: &Path) -> Walk {
    let mut walk = Walk::default();
    walk_module_tree(root_file, true, &mut walk);
    walk
}

// ---------------------------------------------------------------------------------------------
// Suite membership
// ---------------------------------------------------------------------------------------------

struct Member {
    suite_root: String,
    /// The path attribute as written, e.g. `../blastn_strand.rs`.
    path_attr: Option<String>,
    mod_name: String,
}

/// The `mod x;` declarations at the top level of a suite root, with their `#[path]`.
fn suite_members(root: &str) -> Vec<Member> {
    let tokens = tokenize(&read(&manifest_dir().join(root)));
    let mut members = Vec::new();
    let mut pending_path: Option<String> = None;
    let mut depth = 0usize;
    let mut i = 0;
    while i < tokens.len() {
        match &tokens[i] {
            Token::Punct('#') => {
                let mut j = i + 1;
                if tokens.get(j) == Some(&Token::Punct('!')) {
                    j += 1;
                }
                if tokens.get(j) == Some(&Token::Punct('[')) {
                    if let (
                        Some(Token::Ident(key)),
                        Some(Token::Punct('=')),
                        Some(Token::Str(value)),
                    ) = (tokens.get(j + 1), tokens.get(j + 2), tokens.get(j + 3))
                    {
                        if key == "path" && depth == 0 {
                            pending_path = Some(value.clone());
                        }
                    }
                }
                i += 1;
            }
            Token::Punct('{') => {
                depth += 1;
                i += 1;
            }
            Token::Punct('}') => {
                depth = depth.saturating_sub(1);
                i += 1;
            }
            Token::Ident(word) if word == "mod" && depth == 0 => {
                if let (Some(Token::Ident(name)), Some(Token::Punct(';'))) =
                    (tokens.get(i + 1), tokens.get(i + 2))
                {
                    members.push(Member {
                        suite_root: root.to_string(),
                        path_attr: pending_path.take(),
                        mod_name: name.clone(),
                    });
                }
                i += 1;
            }
            _ => i += 1,
        }
    }
    members
}

/// The file stem a member declares through `#[path = "../<file>.rs"]`, if it has that form.
fn member_file(member: &Member) -> Option<String> {
    let file = member.path_attr.as_deref()?.strip_prefix("../")?;
    if file.contains('/') || file.contains('\\') {
        return None;
    }
    file.strip_suffix(".rs").map(str::to_string)
}

fn all_members() -> Vec<Member> {
    suite_roots()
        .iter()
        .flat_map(|root| suite_members(root))
        .collect()
}

fn lines(items: &[String]) -> String {
    items
        .iter()
        .map(|item| format!("  {item}"))
        .collect::<Vec<_>>()
        .join("\n")
}

// ---------------------------------------------------------------------------------------------
// The checks
// ---------------------------------------------------------------------------------------------

#[test]
fn cargo_toml_lists_the_binaries() {
    let (autotests_off, targets) = cargo_targets();
    let mut problems = Vec::new();
    if !autotests_off {
        problems.push(
            "LOSAT/Cargo.toml must set `autotests = false` in [package]: with automatic discovery on, \
             every file of tests/ would also build as a binary of its own"
                .to_string(),
        );
    }
    for target in &targets {
        if target.path.is_empty() || !manifest_dir().join(&target.path).is_file() {
            problems.push(format!(
                "[[test]] `{}` has path `{}`, which is not a file: fix the path in LOSAT/Cargo.toml",
                target.name, target.path
            ));
        }
    }
    for root in suite_roots() {
        if !targets.iter().any(|target| target.path == root) {
            problems.push(format!(
                "{root} is not the path of a [[test]]: add `[[test]]` with name = \"{}\" and path = \"{root}\" to LOSAT/Cargo.toml",
                root.trim_start_matches("tests/suites/").trim_end_matches(".rs")
            ));
        }
    }
    assert!(problems.is_empty(), "\n{}\n", lines(&problems));
}

#[test]
fn every_test_file_is_in_exactly_one_binary() {
    let (_, targets) = cargo_targets();
    let mut homes: BTreeMap<String, Vec<String>> = BTreeMap::new();
    for member in all_members() {
        if let Some(file) = member_file(&member) {
            homes
                .entry(format!("{file}.rs"))
                .or_default()
                .push(member.suite_root.clone());
        }
    }
    for target in &targets {
        if let Some(file) = target.path.strip_prefix("tests/") {
            if !file.contains('/') {
                homes
                    .entry(file.to_string())
                    .or_default()
                    .push(format!("its own [[test]] `{}`", target.name));
            }
        }
    }
    let mut problems = Vec::new();
    for file in top_level_test_files() {
        let stem = file.trim_end_matches(".rs");
        match homes.get(&file).map(Vec::as_slice) {
            None | Some([]) => {
                let suite = suite_for(stem);
                problems.push(format!(
                    "tests/{file} is in no test binary: add `#[path = \"../{file}\"]` and `mod {stem};` to tests/suites/{suite}.rs"
                ));
            }
            Some([_]) => {}
            Some(several) => problems.push(format!(
                "tests/{file} is in {} test binaries ({}): keep it in one",
                several.len(),
                several.join(" and ")
            )),
        }
    }
    assert!(problems.is_empty(), "\n{}\n", lines(&problems));
}

#[test]
fn suite_members_are_top_level_files_named_like_their_module() {
    let mut problems = Vec::new();
    for member in all_members() {
        let root = &member.suite_root;
        match member_file(&member) {
            None => problems.push(format!(
                "{root}: `mod {};` needs `#[path = \"../<file>.rs\"]` naming a file directly in tests/ \
                 (found {:?})",
                member.mod_name, member.path_attr
            )),
            Some(file) => {
                if file != member.mod_name {
                    problems.push(format!(
                        "{root}: module `{}` loads tests/{file}.rs: name the module `{file}`, because test names are \
                         `<file>::<test>` and nextest filters and `child_test_name` rely on it",
                        member.mod_name
                    ));
                }
                if !manifest_dir().join("tests").join(format!("{file}.rs")).is_file() {
                    problems.push(format!("{root}: tests/{file}.rs does not exist"));
                }
            }
        }
    }
    assert!(problems.is_empty(), "\n{}\n", lines(&problems));
}

#[test]
fn no_test_directory_has_a_main_file() {
    let mut problems = Vec::new();
    for entry in fs::read_dir(manifest_dir().join("tests")).expect("read tests/") {
        let path = entry.expect("directory entry").path();
        if path.is_dir() && path.join("main.rs").is_file() {
            problems.push(format!(
                "{} exists: Cargo built tests/<dir>/main.rs as a test binary by itself, and with autotests off \
                 nothing builds it; rename it or add a [[test]] for it in LOSAT/Cargo.toml",
                shown(&path.join("main.rs"))
            ));
        }
    }
    assert!(problems.is_empty(), "\n{}\n", lines(&problems));
}

#[test]
fn suite_members_do_not_change_the_process_environment() {
    let mut problems = Vec::new();
    for member in all_members() {
        let Some(file) = member_file(&member) else {
            continue;
        };
        let root_file = manifest_dir().join("tests").join(format!("{file}.rs"));
        if !root_file.is_file() {
            continue;
        }
        let walk = walk_from(&root_file);
        for unresolved in &walk.unresolved {
            problems.push(format!(
                "{} (member {file}): {unresolved}; the walk cannot follow it, so fix the declaration",
                member.suite_root
            ));
        }
        for (changing, call) in &walk.changes {
            problems.push(format!(
                "member {file}: {} calls {call}...); the tests of a suite share one process: give the file its own \
                 `[[test]]` in LOSAT/Cargo.toml, as blastn_report_stream_api",
                shown(changing)
            ));
        }
    }
    assert!(problems.is_empty(), "\n{}\n", lines(&problems));
}

/// The walk must see the environment changes of the two lone files, the second one through
/// `tests/unit/common/diagnostics.rs` inside `mod unit { pub mod common; }`. When this stops being
/// true, the file no longer needs its own binary and can join a suite.
#[test]
fn walk_finds_the_lone_files_changes() {
    let mut problems = Vec::new();
    for (name, path) in LONE_FILES {
        let walk = walk_from(&manifest_dir().join(path));
        if !walk.unresolved.is_empty() {
            problems.push(format!("{name}: {}", walk.unresolved.join("; ")));
        }
        if walk.changes.is_empty() {
            problems.push(format!(
                "{name} ({path}) no longer changes the process environment, or the module walk misses it: \
                 if the change is gone, the file can join a suite; otherwise fix the walk in tests/suite_layout.rs"
            ));
        }
        if name == "unit_tests"
            && !walk
                .changes
                .iter()
                .any(|(file, _)| shown(file) == "tests/unit/common/diagnostics.rs")
        {
            problems.push(
                "the walk of tests/unit_tests.rs does not reach tests/unit/common/diagnostics.rs through \
                 `mod unit { pub mod common; }`: fix the inline-module handling of the walk"
                    .to_string(),
            );
        }
    }
    assert!(problems.is_empty(), "\n{}\n", lines(&problems));
}
