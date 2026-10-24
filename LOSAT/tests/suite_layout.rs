//! Checks how the files of `LOSAT/tests/` are built into test binaries. `Cargo.toml` turns
//! automatic test discovery off and lists eight binaries: six suites (`tests/suites/<name>.rs`,
//! one `#[path]` module per file) and two files that change the process environment and so stay
//! alone. A file that no binary lists would silently not build, a file in two binaries would run
//! twice, and a suite member that changes the process environment or process-wide state (the panic
//! hook, rayon's global pool) would disturb the other tests of its process. Each test below names
//! every offender and says how to fix it. Standard library only.

use std::collections::BTreeMap;
use std::fs;
use std::path::{Path, PathBuf};

/// The programs, each also the name of a suite.
const PROGRAMS: [&str; 5] = ["tblastx", "tblastn", "blastx", "blastp", "blastn"];

/// The calls that change the process environment or other process-wide state: environment
/// variables, the working directory, the panic hook and rayon's global thread pool.
const PROCESS_CALLS: [&str; 6] = [
    "set_var",
    "remove_var",
    "set_current_dir",
    "set_hook",
    "take_hook",
    "build_global",
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
    /// Every key of the block, in file order.
    keys: Vec<String>,
}

impl Target {
    /// The file name of `path` when it is directly below `tests/` (not in `tests/suites/`).
    fn lone_file(&self) -> Option<&str> {
        self.path
            .strip_prefix("tests/")
            .filter(|file| !file.contains('/'))
    }
}

fn quoted_value(line: &str, key: &str) -> Option<String> {
    let rest = line.strip_prefix(key)?.trim_start().strip_prefix('=')?;
    let rest = rest.trim().strip_prefix('"')?;
    Some(rest[..rest.find('"')?].to_string())
}

/// `autotests = false` under `[package]`, the `[[test]]` blocks of `Cargo.toml`, and every
/// `autotests` line found under another header; read by lines without comments.
fn cargo_targets() -> (bool, Vec<Target>, Vec<String>) {
    let text = read(&manifest_dir().join("Cargo.toml"));
    let mut autotests_off = false;
    let mut targets: Vec<Target> = Vec::new();
    let mut in_test = false;
    let mut section = String::new();
    let mut misplaced: Vec<String> = Vec::new();
    for line in text
        .lines()
        .map(|line| line.split('#').next().unwrap_or("").trim())
    {
        if line.starts_with('[') {
            section = line.to_string();
            in_test = line == "[[test]]";
            if in_test {
                targets.push(Target {
                    name: String::new(),
                    path: String::new(),
                    keys: Vec::new(),
                });
            }
        } else if line.replace(' ', "").starts_with("autotests=") {
            if section != "[package]" {
                misplaced.push(format!("`{line}` is in {section}"));
            } else if line.replace(' ', "") == "autotests=false" {
                autotests_off = true;
            }
        } else if let (true, Some((key, _))) = (in_test, line.split_once('=')) {
            let target = targets.last_mut().expect("a [[test]] block");
            target.keys.push(key.trim().to_string());
            if let Some(name) = quoted_value(line, "name") {
                target.name = name;
            }
            if let Some(path) = quoted_value(line, "path") {
                target.path = path;
            }
        }
    }
    (autotests_off, targets, misplaced)
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
// A small Rust tokenizer: identifiers, punctuation and string contents (only `#[path]` values are
// read, so escapes are dropped); comments and character literals never reach the checks.
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
                        let (text, next) = read_raw_string(&chars, j + 1, j - i);
                        tokens.push(Token::Str(text));
                        i = next;
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

/// Reads a string whose opening quote is before `start`; returns the text (an escape adds nothing)
/// and the index after the closing quote.
fn read_string(chars: &[char], start: usize) -> (String, usize) {
    let mut text = String::new();
    let mut i = start;
    while i < chars.len() && chars[i] != '"' {
        if chars[i] == '\\' {
            i += 1;
        } else {
            text.push(chars[i]);
        }
        i += 1;
    }
    (text, i + 1)
}

/// Reads a raw string whose opening quote is before `start`, closed by a quote and `hashes` hash
/// marks; returns the text and the index after the closing marks.
fn read_raw_string(chars: &[char], start: usize, hashes: usize) -> (String, usize) {
    let mut i = start;
    while i < chars.len() {
        if chars[i] == '"' && (0..hashes).all(|k| chars.get(i + 1 + k) == Some(&'#')) {
            return (chars[start..i].iter().collect(), i + 1 + hashes);
        }
        i += 1;
    }
    (chars[start..].iter().collect(), i)
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
// Module declarations, and the module tree of a test file
// ---------------------------------------------------------------------------------------------

/// A `mod name;` declaration.
struct Decl {
    name: String,
    /// The value of a `#[path = "..."]` attribute on it.
    path: Option<String>,
    /// The names of the inline `mod a { ... }` blocks around it, outermost first.
    inline: Vec<String>,
}

/// Every `mod x;` declaration of a token stream (with or without `pub` and attributes).
fn declarations(tokens: &[Token]) -> Vec<Decl> {
    let mut decls = Vec::new();
    // Open braces: Some(name) for `mod name {`, None for any other.
    let mut braces: Vec<Option<String>> = Vec::new();
    let mut path: Option<String> = None;
    let mut i = 0;
    while i < tokens.len() {
        match &tokens[i] {
            Token::Punct('#') => {
                // An attribute, `#[...]` or `#![...]`: remember `path = "..."`, skip the rest.
                let start = i + 1 + usize::from(tokens.get(i + 1) == Some(&Token::Punct('!')));
                if tokens.get(start) != Some(&Token::Punct('[')) {
                    i += 1;
                    continue;
                }
                if let Some([Token::Ident(key), Token::Punct('='), Token::Str(value)]) =
                    tokens.get(start + 1..start + 4)
                {
                    if key == "path" {
                        path = Some(value.clone());
                    }
                }
                let mut depth = 0;
                let mut end = start;
                while end < tokens.len() {
                    match tokens[end] {
                        Token::Punct('[') => depth += 1,
                        Token::Punct(']') => depth -= 1,
                        _ => {}
                    }
                    if depth == 0 {
                        break;
                    }
                    end += 1;
                }
                i = end + 1;
            }
            Token::Ident(word) if word == "mod" => match (tokens.get(i + 1), tokens.get(i + 2)) {
                (Some(Token::Ident(name)), Some(Token::Punct('{'))) => {
                    braces.push(Some(name.clone()));
                    path = None;
                    i += 3;
                }
                (Some(Token::Ident(name)), Some(Token::Punct(';'))) => {
                    decls.push(Decl {
                        name: name.clone(),
                        path: path.take(),
                        inline: braces.iter().flatten().cloned().collect(),
                    });
                    i += 3;
                }
                _ => i += 1,
            },
            Token::Punct('{') => {
                braces.push(None);
                path = None;
                i += 1;
            }
            Token::Punct('}') => {
                braces.pop();
                path = None;
                i += 1;
            }
            Token::Punct(';') => {
                path = None;
                i += 1;
            }
            _ => i += 1,
        }
    }
    decls
}

/// The calls of `PROCESS_CALLS` in a token stream, written `name(`.
fn process_calls(tokens: &[Token]) -> Vec<String> {
    let mut calls = Vec::new();
    for pair in tokens.windows(2) {
        if let [Token::Ident(word), Token::Punct('(')] = pair {
            if PROCESS_CALLS.contains(&word.as_str()) {
                calls.push(format!("{word}("));
            }
        }
    }
    calls
}

#[derive(Default)]
struct Walk {
    visited: Vec<PathBuf>,
    /// (file, call) for every environment or process-wide change found.
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
    for call in process_calls(&tokens) {
        walk.changes.push((file.to_path_buf(), call));
    }
    for decl in declarations(&tokens) {
        let name = &decl.name;
        // An inline block `mod a { mod x; }` looks in `a/`; so does a `#[path]` inside it.
        let mut here = children_dir.clone();
        for part in &decl.inline {
            here = here.join(part);
        }
        if let Some(relative) = decl.path {
            let target = if decl.inline.is_empty() { &dir } else { &here }.join(relative);
            if target.is_file() {
                walk_module_tree(&target, true, walk);
            } else {
                walk.unresolved.push(format!(
                    "`mod {name};` with #[path] in {} points to {}, which does not exist",
                    shown(file),
                    shown(&target)
                ));
            }
            continue;
        }
        let plain = here.join(format!("{name}.rs"));
        let nested = here.join(name).join("mod.rs");
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

/// The declarations of a suite root outside any inline block.
fn suite_members(root: &str) -> Vec<Member> {
    declarations(&tokenize(&read(&manifest_dir().join(root))))
        .into_iter()
        .filter(|decl| decl.inline.is_empty())
        .map(|decl| Member {
            suite_root: root.to_string(),
            path_attr: decl.path,
            mod_name: decl.name,
        })
        .collect()
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
    let (autotests_off, targets, misplaced) = cargo_targets();
    let mut problems = Vec::new();
    for line in &misplaced {
        problems.push(format!(
            "LOSAT/Cargo.toml: {line}, but `autotests` must be in [package]: Cargo ignores it anywhere else \
             and discovers every file of tests/ as a binary of its own"
        ));
    }
    if !autotests_off {
        problems.push(
            "LOSAT/Cargo.toml must set `autotests = false` in [package]: with automatic discovery on, \
             every file of tests/ would also build as a binary of its own"
                .to_string(),
        );
    }
    for target in &targets {
        let mut keys = target.keys.clone();
        keys.sort();
        if keys != ["name", "path"] {
            problems.push(format!(
                "[[test]] `{}` has the keys {:?}: it must have exactly `name` and `path` (`test`, `harness`, \
                 `required-features` and the like can stop a suite from running or run it without the features \
                 of the gates)",
                target.name, target.keys
            ));
        } else if !manifest_dir().join(&target.path).is_file() {
            problems.push(format!(
                "[[test]] `{}` has path `{}`, which is not a file: fix the path in LOSAT/Cargo.toml",
                target.name, target.path
            ));
        } else if Path::new(&target.path)
            .file_stem()
            .and_then(|stem| stem.to_str())
            != Some(target.name.as_str())
        {
            problems.push(format!(
                "[[test]] `{}` has path `{}`: the name must be the file stem (`--test <suite>` and the \
                 nextest binary id `LOSAT::<suite>` rely on it)",
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
    let (_, targets, _) = cargo_targets();
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
        if let Some(file) = target.lone_file() {
            homes
                .entry(file.to_string())
                .or_default()
                .push(format!("its own [[test]] `{}`", target.name));
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
    for root in suite_roots() {
        let tokens = tokenize(&read(&manifest_dir().join(&root)));
        for decl in declarations(&tokens)
            .iter()
            .filter(|decl| !decl.inline.is_empty())
        {
            problems.push(format!(
                "{root}: `mod {};` is inside `mod {} {{ ... }}`: list members only at the top level",
                decl.name,
                decl.inline.join("::")
            ));
        }
        for call in process_calls(&tokens) {
            problems.push(format!(
                "{root}: calls {call}...) in the suite root itself, which changes the process environment or \
                 process-wide state for every test of the suite: a suite root holds only member declarations"
            ));
        }
    }
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
            let fix = if *changing == root_file {
                "give the file its own `[[test]]` in LOSAT/Cargo.toml, as blastn_report_stream_api"
                    .to_string()
            } else {
                format!(
                    "move the change out of {}, or give member {file} its own `[[test]]` in LOSAT/Cargo.toml, \
                     as blastn_report_stream_api",
                    shown(changing)
                )
            };
            problems.push(format!(
                "member {file}: {} calls {call}...), which changes the process environment or process-wide \
                 state, and the tests of a suite share one process: {fix}",
                shown(changing)
            ));
        }
    }
    assert!(problems.is_empty(), "\n{}\n", lines(&problems));
}

/// The walk must see the environment changes of the lone files, the one of `unit_tests` through
/// `tests/unit/common/diagnostics.rs` inside `mod unit { pub mod common; }`. When this stops being
/// true, the file no longer needs its own binary and can join a suite.
#[test]
fn walk_finds_the_lone_files_changes() {
    let (_, targets, _) = cargo_targets();
    let mut problems = Vec::new();
    for target in targets.iter().filter(|target| target.lone_file().is_some()) {
        let (name, path) = (&target.name, &target.path);
        let walk = walk_from(&manifest_dir().join(path));
        if !walk.unresolved.is_empty() {
            problems.push(format!("{name}: {}", walk.unresolved.join("; ")));
        }
        if walk.changes.is_empty() {
            problems.push(format!(
                "{name} ({path}) no longer changes the process environment or process-wide state, or the module walk misses it: \
                 if the change is gone, the file can join a suite; otherwise fix the walk in tests/suite_layout.rs"
            ));
        }
    }
    // The inline-module handling, checked on the file whose change sits inside `mod unit { ... }`,
    // whether or not the file still has a binary of its own.
    let unit_walk = walk_from(&manifest_dir().join("tests/unit_tests.rs"));
    if !unit_walk
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
    assert!(problems.is_empty(), "\n{}\n", lines(&problems));
}
