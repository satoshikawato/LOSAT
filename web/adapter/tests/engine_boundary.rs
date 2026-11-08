//! The adapter's boundary with the engine crate (`LOSAT/src/web_adapter.rs`; line LIB,
//! DECISIONS 13): the adapter names the engine through the public `api` and the hidden façade
//! `web_adapter` only, the façade names each item on its own, and every item of the façade is
//! one that the adapter uses. A new engine item for the adapter is then one line of the façade
//! in the pull request's diff.

use std::collections::BTreeSet;
use std::path::{Path, PathBuf};

/// The engine crate's name, kept out of a literal path so that this file passes its own scan.
const ENGINE: &str = "LOSAT";
const FACADE: &str = "web_adapter";

fn adapter_dir() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
}

fn rust_files(dir: &Path, out: &mut Vec<PathBuf>) {
    let mut entries: Vec<PathBuf> = std::fs::read_dir(dir)
        .unwrap_or_else(|error| panic!("{}: {error}", dir.display()))
        .map(|entry| entry.expect("a directory entry").path())
        .collect();
    entries.sort();
    for path in entries {
        if path.is_dir() {
            rust_files(&path, out);
        } else if path.extension().is_some_and(|ext| ext == "rs") {
            out.push(path);
        }
    }
}

/// The Rust sources of the adapter crate (`src/`, `tests/`), with their text.
fn adapter_sources() -> Vec<(PathBuf, String)> {
    let mut paths = Vec::new();
    for dir in ["src", "tests"] {
        rust_files(&adapter_dir().join(dir), &mut paths);
    }
    assert!(paths.len() > 3, "the adapter's sources: {paths:?}");
    paths
        .into_iter()
        .map(|path| {
            let text = std::fs::read_to_string(&path).expect("a source file");
            (path, text)
        })
        .collect()
}

fn facade_source() -> String {
    let path = adapter_dir().join("../../LOSAT/src/web_adapter.rs");
    std::fs::read_to_string(&path).unwrap_or_else(|error| panic!("{}: {error}", path.display()))
}

/// `text` without its `//` comments (the sources here hold `//` in no string).
fn code(text: &str) -> String {
    text.lines()
        .map(|line| line.find("//").map_or(line, |at| &line[..at]))
        .collect::<Vec<_>>()
        .join("\n")
}

fn is_ident(byte: u8) -> bool {
    byte.is_ascii_alphanumeric() || byte == b'_'
}

/// The identifier at the start of `text`.
fn ident(text: &str) -> &str {
    let end = text
        .bytes()
        .position(|byte| !is_ident(byte))
        .unwrap_or(text.len());
    &text[..end]
}

/// The starts of `word` in `text` that are not inside a longer identifier.
fn word_starts<'a>(text: &'a str, word: &'a str) -> impl Iterator<Item = usize> + 'a {
    text.match_indices(word)
        .map(|(at, _)| at)
        .filter(move |&at| {
            let before = at.checked_sub(1).map(|i| text.as_bytes()[i]);
            let after = text.as_bytes().get(at + word.len()).copied();
            !before.is_some_and(is_ident) && !after.is_some_and(is_ident)
        })
}

/// The names of a use tree after its prefix: `name`, `name as alias` or `{a, b as c}` (one
/// level). `take_alias` takes the alias of `as` (what a `pub use` exports) instead of the name.
fn use_tree_names(tree: &str, take_alias: bool) -> Vec<String> {
    let tree = tree.trim();
    let items: Vec<&str> = match tree.strip_prefix('{') {
        Some(group) => {
            let group = group
                .strip_suffix('}')
                .unwrap_or_else(|| panic!("a group: {tree}"));
            assert!(!group.contains('{'), "one level of braces: {tree}");
            group
                .split(',')
                .map(str::trim)
                .filter(|item| !item.is_empty())
                .collect()
        }
        None => vec![tree],
    };
    items
        .into_iter()
        .map(|item| {
            assert!(!item.contains('*'), "a glob: {item}");
            let (name, alias) = item
                .split_once(" as ")
                .map_or((item, None), |(n, a)| (n, Some(a)));
            let name = if take_alias {
                alias.unwrap_or(name)
            } else {
                name
            };
            assert!(!name.contains("::"), "one name per item: {item}");
            name.trim().to_string()
        })
        .collect()
}

#[test]
fn the_adapter_names_the_engine_through_the_api_and_the_facade_only() {
    let mut wrong = Vec::new();
    for (path, text) in adapter_sources() {
        for (line_number, line) in text.lines().enumerate() {
            for at in word_starts(line, ENGINE) {
                let rest = &line[at + ENGINE.len()..];
                let Some(path_rest) = rest.strip_prefix("::") else {
                    // `LOSAT` as a word of prose, or an alias of the crate.
                    if rest.trim_start().starts_with("as ") || line.trim_start().starts_with("use ")
                    {
                        wrong.push(format!("{}:{}: {line}", path.display(), line_number + 1));
                    }
                    continue;
                };
                let module = ident(path_rest);
                let ok = match module {
                    // The public surface; its hidden child `local_blast` is the façade's.
                    "api" => !path_rest.starts_with("api::local_blast"),
                    FACADE => {
                        // One item after the façade; a further `::` only for an associated item
                        // of a type (`ScoringMatrix::ALL`), not for a module re-exported whole.
                        let after = &path_rest[FACADE.len()..];
                        match after.strip_prefix("::") {
                            Some(item) if item.starts_with('{') => true,
                            Some(item) => {
                                let name = ident(item);
                                !item[name.len()..].starts_with("::")
                                    || name.starts_with(|c: char| c.is_ascii_uppercase())
                            }
                            None => false,
                        }
                    }
                    _ => false,
                };
                if !ok {
                    wrong.push(format!("{}:{}: {line}", path.display(), line_number + 1));
                }
            }
        }
    }
    assert!(
        wrong.is_empty(),
        "the adapter names the engine through `api` and `web_adapter` only (add the item to \
         LOSAT/src/web_adapter.rs):\n{}",
        wrong.join("\n")
    );
}

/// The names that the façade exports: its `pub use` items (their aliases) and its `pub fn`s.
fn facade_names() -> BTreeSet<String> {
    let source = code(&facade_source());
    let mut names = BTreeSet::new();
    let mut rest = source.as_str();
    while let Some(at) = rest.find("pub ") {
        let item = &rest[at + "pub ".len()..];
        if let Some(tree) = item.strip_prefix("use ") {
            let end = tree.find(';').expect("a `pub use` ends with `;`");
            let tree = tree[..end].split_whitespace().collect::<Vec<_>>().join(" ");
            let tree = tree
                .replace("{ ", "{")
                .replace(" }", "}")
                .replace(", }", "}");
            let path = tree
                .strip_prefix("crate::")
                .expect("a `pub use` of this crate");
            // The prefix ends at the last `::` before the group or the item.
            let split = match path.find('{') {
                Some(brace) => brace,
                None => path.rfind("::").map(|at| at + 2).expect("a path"),
            };
            names.extend(use_tree_names(&path[split..], true));
            rest = tree_rest(item);
        } else if let Some(function) = item.strip_prefix("fn ") {
            names.insert(ident(function).to_string());
            rest = function;
        } else {
            panic!(
                "the façade holds `pub use` and `pub fn` items only: pub {}",
                ident(item)
            );
        }
    }
    names
}

/// `item` after the end of its `pub use` statement.
fn tree_rest(item: &str) -> &str {
    &item[item.find(';').expect("a `pub use` ends with `;`") + 1..]
}

/// The façade's names that the adapter's sources use.
fn used_names() -> BTreeSet<String> {
    let prefix = format!("{ENGINE}::{FACADE}::");
    let mut names = BTreeSet::new();
    for (_, text) in adapter_sources() {
        let text = code(&text);
        for (at, _) in text.match_indices(&prefix) {
            let rest = &text[at + prefix.len()..];
            if rest.starts_with('{') {
                let end = rest.find('}').expect("a closed group");
                let group = rest[..=end]
                    .split_whitespace()
                    .collect::<Vec<_>>()
                    .join(" ");
                names.extend(use_tree_names(&group, false));
            } else {
                names.insert(ident(rest).to_string());
            }
        }
    }
    names
}

#[test]
fn every_item_of_the_facade_is_one_the_adapter_uses() {
    let exported = facade_names();
    let used = used_names();
    assert!(exported.len() > 40, "the façade's names: {exported:?}");
    let unused: Vec<&String> = exported.difference(&used).collect();
    let missing: Vec<&String> = used.difference(&exported).collect();
    assert!(
        unused.is_empty() && missing.is_empty(),
        "LOSAT/src/web_adapter.rs exports what the adapter does not use: {unused:?}; the \
         adapter uses what it does not export: {missing:?}"
    );
}
