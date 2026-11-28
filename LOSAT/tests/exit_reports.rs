//! The reports that the end of a run writes (`src/utils/xexit.rs`, `write`) run only where a
//! variable of its `NAMES` is set, or, in a diagnostic build, one of `xdiag::NAMES`: a run that
//! asks for none runs none of their code (line RG-KB item 9a-4).
//!
//! The first test reads the source. From each call of `write` it finds the function that runs
//! (by its module's file and through `pub use .. as` and `use`), then the functions that function
//! calls, and so on (the report's functions), and fails:
//! - where a report's functions quote a `LOSAT_*` name that is not in `NAMES`, in `xdiag::NAMES`
//!   or a switch of `utils/xswitch.rs`, or read a variable through a name that is not quoted (a
//!   constant, a computed name), unless the name is the function's own parameter (a reader such
//!   as `xswitch::on`, whose callers quote it);
//! - where a report tests no variable of `NAMES` or of `xdiag` on its way (a report with an empty
//!   body writes nothing), or a name of `NAMES` is read by no report;
//! - where a line is written (`eprintln!`, `println!`, `write!`, `fs::write`, ..) that no test of
//!   such a variable guards: the condition of a block around it (`if`, `if let`, `match`, `for`,
//!   `while`; not `else`) or an earlier block of its level whose condition tests and which
//!   returns (`if !on() { return; }`) tests one, by quoting it or by calling a function that does,
//!   or the function was called where a test guards the call. A counter or a switch alone is not
//!   a test: a counter counts whatever its own code counts, and a switch is on unless `0`, so a
//!   line under either alone would be written, or its report's code run, in every run;
//! - where `cli/main_entry.rs` does not end the function that calls `write` with its `match` on
//!   the command, then `crate::utils::xexit::write(); Ok(())`, or calls a report of `write`
//!   elsewhere: a report called there, before or after `write`, would run in every run again.
//!
//! The finder knows a function by its name and the file stem of its module (`a/b.rs` or
//! `a/b/mod.rs` holds module `b`), the nearest directory first where several files have that
//! stem, as LOSAT's modules are laid out; a call it cannot place (a method, a closure, `std`) is
//! not followed, and a test it reads must be in the text of a condition, not in a local
//! variable (`let on = ..; if on`).
//!
//! The second test runs one search with no such variable and with each variable alone: without one
//! no report line is written, and each variable still gives its report.
//!
//! No NCBI counterpart: LOSAT's own diagnostics.

use std::collections::HashMap;
use std::path::Path;
use std::process::{Command, Output};

use super::function_size_limit::{
    closers, fn_item, ident, literal_end, mask, read_rust_files, word_before,
};

/// A file below `src` (not a `_tests.rs` file): its path, its text with comments blanked
/// (literals kept), the text with comments and literals blanked (`mask`), its functions (name,
/// offset of `fn`, offset of the closing brace) and its `use` items (the name each binds, the
/// path).
struct Source {
    rel: String,
    code: Vec<u8>,
    masked: Vec<u8>,
    fns: Vec<(String, usize, usize)>,
    uses: Vec<(String, Vec<String>)>,
}

fn src_dir() -> std::path::PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("src")
}

fn sources() -> Vec<Source> {
    let mut files = Vec::new();
    read_rust_files(&src_dir(), "", &mut files);
    files.retain(|(rel, _)| !rel.ends_with("_tests.rs"));
    files
        .into_iter()
        .map(|(rel, text)| {
            let masked = mask(&text);
            let close = closers(&masked);
            let fns = (0..masked.len())
                .filter(|&i| !word_before(&masked, i) && masked[i..].starts_with(b"fn"))
                .filter_map(|i| {
                    let (name, end) = fn_item(&masked, &close, i)?;
                    Some((String::from_utf8_lossy(name).into_owned(), i, end))
                })
                .collect();
            let uses = uses(&masked);
            Source {
                rel,
                code: strip_comments(&text),
                masked,
                fns,
                uses,
            }
        })
        .collect()
}

/// `text` with its comments turned into spaces (newlines kept); its literals stay.
fn strip_comments(text: &[u8]) -> Vec<u8> {
    let (mut code, mut i) = (text.to_vec(), 0);
    while i < text.len() {
        let Some(end) = literal_end(text, i) else {
            i += 1;
            continue;
        };
        if text[i] == b'/' {
            code[i..end]
                .iter_mut()
                .filter(|b| **b != b'\n')
                .for_each(|b| *b = b' ');
        }
        i = end;
    }
    code
}

/// The `use` items of a masked file: (the name bound, the path), groups expanded.
fn uses(m: &[u8]) -> Vec<(String, Vec<String>)> {
    let mut out = Vec::new();
    for i in (0..m.len()).filter(|&i| !word_before(m, i) && ident(m, i) == b"use") {
        let Some(end) = m[i..].iter().position(|&b| b == b';') else {
            continue;
        };
        let words = String::from_utf8_lossy(&m[i + 3..i + end]).into_owned();
        let tree = words.split_whitespace().collect::<Vec<_>>().join(" ");
        let tree = ["::", "{", "}", ","].iter().fold(tree, |t, p| {
            t.replace(&format!(" {p}"), p).replace(&format!("{p} "), p)
        });
        use_tree(&tree, &[], &mut out);
    }
    out
}

fn use_tree(tree: &str, prefix: &[String], out: &mut Vec<(String, Vec<String>)>) {
    if tree.is_empty() {
        return;
    }
    let words = |text: &str| -> Vec<String> {
        let parts = text.split("::").filter(|part| !part.is_empty());
        prefix
            .iter()
            .cloned()
            .chain(parts.map(String::from))
            .collect()
    };
    if let (Some(open), Some(close)) = (tree.find('{'), tree.rfind('}')) {
        let head = words(tree[..open].trim_end_matches(':'));
        let (mut depth, mut start) = (0, open + 1);
        for (k, b) in tree.bytes().enumerate().take(close + 1).skip(open + 1) {
            match b {
                b'{' => depth += 1,
                b'}' if depth > 0 => depth -= 1,
                b',' | b'}' if depth == 0 => {
                    use_tree(&tree[start..k], &head, out);
                    start = k + 1;
                }
                _ => {}
            }
        }
        return;
    }
    let (path, alias) = match tree.split_once(" as ") {
        Some((path, alias)) => (path, Some(alias.to_string())),
        None => (tree, None),
    };
    let mut path = words(path);
    if path.last().is_some_and(|last| last == "self") {
        path.pop();
    }
    if let Some(last) = path.last().filter(|last| last.as_str() != "*") {
        out.push((alias.unwrap_or_else(|| last.clone()), path));
    }
}

/// The module a file holds, by its name: `a/b.rs` and `a/b/mod.rs` hold `b`.
fn stem(rel: &str) -> &str {
    let mut parts = rel.rsplit('/');
    let file = parts.next().unwrap_or_default();
    match file {
        "mod.rs" => parts.next().unwrap_or_default(),
        _ => file.trim_end_matches(".rs"),
    }
}

/// The directories two files share, from the top.
fn shared_dirs(a: &str, b: &str) -> usize {
    let dirs = |rel: &str| rel.rsplit_once('/').map_or("", |d| d.0).to_string();
    let (a, b) = (dirs(a), dirs(b));
    a.split('/')
        .zip(b.split('/'))
        .take_while(|(x, y)| x == y)
        .count()
}

/// A function found: its file, the offset of its `fn` and of its closing brace.
type Found = (usize, usize, usize);

/// The functions that `path`, called in file `k`, names (see the module doc); `depth` bounds the
/// `use` items followed.
fn resolve(src: &[Source], k: usize, path: &[String], depth: usize) -> Vec<Found> {
    let segs: Vec<&String> = path
        .iter()
        .filter(|s| !matches!(s.as_str(), "crate" | "self" | "super"))
        .collect();
    // `std::..` and `core::..` are the standard library's; `crate::core::..` is LOSAT's `core`.
    if depth > 8
        || path
            .first()
            .is_some_and(|s| matches!(s.as_str(), "std" | "core" | "alloc"))
    {
        return Vec::new();
    }
    let (module, name) = match segs.as_slice() {
        [] => return Vec::new(),
        [name] => return named(src, k, name, depth),
        [.., module, name] => (module.as_str(), name.as_str()),
    };
    // A module this file imports (`use super::x_kmerge;`) is found by its own name below; one
    // imported under another name is replaced by its path.
    if segs.len() == 2 {
        if let Some((_, p)) = src[k]
            .uses
            .iter()
            .find(|(n, p)| n == module && p.last().is_some_and(|l| l != module))
        {
            let full: Vec<String> = p.iter().cloned().chain([name.to_string()]).collect();
            return resolve(src, k, &full, depth + 1);
        }
    }
    let mut best: (usize, Vec<Found>) = (0, Vec::new());
    for (j, file) in src
        .iter()
        .enumerate()
        .filter(|(_, f)| stem(&f.rel) == module)
    {
        let found = named(src, j, name, depth);
        let near = shared_dirs(&src[k].rel, &file.rel);
        if found.is_empty() || (near < best.0 && !best.1.is_empty()) {
            continue;
        }
        if near > best.0 || best.1.is_empty() {
            best = (near, found);
        } else {
            best.1.extend(found);
        }
    }
    best.1
}

/// The functions named `name` in file `k`: its own, or the ones a `use` item binds to the name.
fn named(src: &[Source], k: usize, name: &str, depth: usize) -> Vec<Found> {
    let file = &src[k];
    let mut found: Vec<Found> = file
        .fns
        .iter()
        .filter(|f| f.0 == name)
        .map(|f| (k, f.1, f.2))
        .collect();
    for (bound, path) in &file.uses {
        if bound == name
            && path
                .last()
                .is_some_and(|last| last != name || path.len() > 1)
        {
            found.extend(resolve(src, k, path, depth + 1));
        }
    }
    found
}

/// The paths called in `m[start..end]` (`a::b::f(`), outside method calls, with their offsets.
fn calls_at(m: &[u8], start: usize, end: usize) -> Vec<(usize, Vec<String>)> {
    let (mut out, mut i) = (Vec::new(), start);
    while i < end {
        let first = ident(m, i);
        if first.is_empty() || word_before(m, i) {
            i += 1;
            continue;
        }
        let method = m[..i].iter().rev().find(|b| !b.is_ascii_whitespace()) == Some(&b'.');
        let (mut path, mut j) = (Vec::new(), i);
        loop {
            let word = ident(m, j);
            path.push(String::from_utf8_lossy(word).into_owned());
            j += word.len();
            if m[j..].starts_with(b"::") && !ident(m, j + 2).is_empty() {
                j += 2;
            } else {
                break;
            }
        }
        let next = m[j..].iter().find(|b| !b.is_ascii_whitespace());
        if !method && next == Some(&b'(') {
            out.push((i, path));
        }
        i = j;
    }
    out
}

/// The variables that `code` reads from the environment by a quoted name: (through `xdiag`, the
/// name).
fn reads(code: &[u8]) -> Vec<(bool, String)> {
    let squeezed: String = String::from_utf8_lossy(code).split_whitespace().collect();
    let mut out = Vec::new();
    for pattern in ["var_os(\"", "var(\""] {
        for (at, _) in squeezed.match_indices(pattern) {
            let before = squeezed[..at].chars().next_back();
            if before.is_some_and(|c| c.is_ascii_alphanumeric() || c == '_') {
                continue;
            }
            let rest = &squeezed[at + pattern.len()..];
            let name = &rest[..rest.find('"').unwrap_or(0)];
            if !name.is_empty()
                && name
                    .bytes()
                    .all(|b| b.is_ascii_uppercase() || b.is_ascii_digit() || b == b'_')
            {
                out.push((squeezed[..at].ends_with("xdiag::"), name.to_string()));
            }
        }
    }
    out
}

/// The arguments of the environment reads of `code` that are not a quoted name (`var_os(name)`).
fn unquoted_reads(code: &[u8]) -> Vec<String> {
    let squeezed: String = String::from_utf8_lossy(code).split_whitespace().collect();
    let mut out = Vec::new();
    for pattern in ["var_os(", "var("] {
        for (at, _) in squeezed.match_indices(pattern) {
            let before = squeezed[..at].chars().next_back();
            if before.is_some_and(|c| c.is_ascii_alphanumeric() || c == '_' || c == '.') {
                continue;
            }
            let rest = &squeezed[at + pattern.len()..];
            let arg = &rest[..rest.find([')', ',']).unwrap_or(0)];
            if !arg.is_empty() && !arg.starts_with('"') {
                out.push(arg.trim_start_matches('&').to_string());
            }
        }
    }
    out
}

/// Whether `arg` is a parameter of the function whose text is `code` (in its signature, or of a
/// closure in it: `|name|`): a reader that gets the name from its caller, who quotes it.
fn bound(code: &[u8], arg: &str) -> bool {
    let squeezed: String = String::from_utf8_lossy(code).split_whitespace().collect();
    let signature = &squeezed[..squeezed.find('{').unwrap_or(0)];
    let word = |text: &str| {
        text.split(|c: char| !(c.is_ascii_alphanumeric() || c == '_'))
            .any(|w| w == arg)
    };
    word(signature)
        || [format!("|{arg}|"), format!("|&{arg}|")]
            .iter()
            .any(|c| squeezed.contains(c))
}

/// The string literals of `code` that are a whole `LOSAT_*` name.
fn name_literals(code: &[u8]) -> Vec<String> {
    let (mut out, mut i) = (Vec::new(), 0);
    while i < code.len() {
        let Some(end) = literal_end(code, i) else {
            i += 1;
            continue;
        };
        let literal = &code[i..end];
        if literal.len() > 2 && literal[0] == b'"' {
            let name = &literal[1..literal.len() - 1];
            if name.starts_with(b"LOSAT_")
                && name
                    .iter()
                    .all(|b| b.is_ascii_uppercase() || b.is_ascii_digit() || *b == b'_')
            {
                out.push(String::from_utf8_lossy(name).into_owned());
            }
        }
        i = end;
    }
    out
}

/// The offsets in `m[start..end]` where a line is written: `eprintln!`, `eprint!`, `println!`,
/// `print!`, `writeln!`, `write!`, `fs::write(`.
fn writes(m: &[u8], start: usize, end: usize) -> Vec<usize> {
    let tokens: [&[u8]; 7] = [
        b"eprintln!",
        b"eprint!",
        b"println!",
        b"print!",
        b"writeln!",
        b"write!",
        b"fs::write(",
    ];
    (start..end)
        .filter(|&i| !word_before(m, i) && tokens.iter().any(|t| m[i..].starts_with(t)))
        .collect()
}

/// The offset of the `{` that the `}` at `close` closes (`m` masked).
fn opener(m: &[u8], close: usize) -> usize {
    let mut depth = 0;
    for i in (0..=close).rev() {
        match m[i] {
            b'}' => depth += 1,
            b'{' => {
                depth -= 1;
                if depth == 0 {
                    return i;
                }
            }
            _ => {}
        }
    }
    0
}

/// Which functions test a variable of `NAMES` or of `xdiag` (in their own text or in a function
/// they call), remembered per function, and where a test guards a function's lines.
struct Tests<'a> {
    src: &'a [Source],
    names: &'a [String],
    known: HashMap<(usize, usize), bool>,
}

impl Tests<'_> {
    /// Whether the function at `s..=e` of `src[j]` tests a variable.
    fn function(&mut self, j: usize, s: usize, e: usize) -> bool {
        if let Some(&known) = self.known.get(&(j, s)) {
            return known;
        }
        // A cycle reads as untested until its first function is known.
        self.known.insert((j, s), false);
        let tests = self.range(j, s, e + 1);
        self.known.insert((j, s), tests);
        tests
    }

    /// Whether `src[j]` tests a variable in `from..to`: quotes one in a read, or calls a function
    /// that tests one.
    fn range(&mut self, j: usize, from: usize, to: usize) -> bool {
        let src = self.src;
        if reads(&src[j].code[from..to])
            .iter()
            .any(|(xdiag, name)| *xdiag || self.names.contains(name))
        {
            return true;
        }
        for (_, path) in calls_at(&src[j].masked, from, to) {
            for (a, b, c) in resolve(src, j, &path, 0) {
                if self.function(a, b, c) {
                    return true;
                }
            }
        }
        false
    }

    /// For each offset of the function at `s..=e` of `src[j]`: whether a test guards it (see the
    /// module doc); `entry`: the function is called where a test guards the call.
    fn guards(&mut self, j: usize, s: usize, e: usize, entry: bool) -> Vec<bool> {
        let src = self.src;
        let m = &src[j].masked;
        let open = opener(m, e);
        let mut out = vec![entry; e + 1 - s];
        // (guarded, its condition tests, its `{`) for each block around the offset.
        let mut blocks = vec![(entry, false, open)];
        let mut statement = open + 1;
        for i in open + 1..e {
            out[i - s] = blocks.last().map_or(entry, |b| b.0);
            match m[i] {
                b'{' => {
                    let head: String = String::from_utf8_lossy(&m[statement..i])
                        .split_whitespace()
                        .collect();
                    // An `else` block runs where its `if` did not test true.
                    let plain_else = head.starts_with("else") && !head.starts_with("elseif");
                    let tests = !plain_else && self.range(j, statement, i);
                    let guarded = blocks.last().map_or(entry, |b| b.0) || tests;
                    blocks.push((guarded, tests, i));
                    statement = i + 1;
                }
                b'}' if blocks.len() > 1 => {
                    let (_, tests, block) = blocks.pop().expect("a block");
                    let returns =
                        (block..i).any(|k| !word_before(m, k) && ident(m, k) == b"return");
                    if tests && returns {
                        blocks.last_mut().expect("the function's block").0 = true;
                    }
                    statement = i + 1;
                }
                b';' => statement = i + 1,
                _ => {}
            }
        }
        out
    }
}

/// The line of offset `at` in `text`.
fn line(text: &[u8], at: usize) -> usize {
    1 + text[..at].iter().filter(|&&b| b == b'\n').count()
}

/// The quoted names after `anchor` up to the next `]` in `text`.
fn quoted_after(text: &str, anchor: &str) -> Vec<String> {
    let start = text
        .find(anchor)
        .unwrap_or_else(|| panic!("`{anchor}` not found"))
        + anchor.len();
    let block = &text[start..start + text[start..].find(']').expect("the end of the list")];
    block
        .split('"')
        .skip(1)
        .step_by(2)
        .map(String::from)
        .collect()
}

/// The names the source passes to `xswitch::on` and `xswitch::value`: the strict set's switches.
fn switches(src: &[Source]) -> Vec<String> {
    let mut names = Vec::new();
    for file in src {
        let squeezed: String = String::from_utf8_lossy(&file.code)
            .split_whitespace()
            .collect();
        for pattern in ["xswitch::on(\"", "xswitch::value(\""] {
            for (at, _) in squeezed.match_indices(pattern) {
                let rest = &squeezed[at + pattern.len()..];
                names.push(rest[..rest.find('"').unwrap_or(0)].to_string());
            }
        }
    }
    names
}

/// `cli/main_entry.rs`: the function that calls `xexit::write` ends with its `match` on the
/// command, then `crate::utils::xexit::write(); Ok(())`, and calls no report of `write` (`reports`)
/// elsewhere.
fn entry_ends_with_write(src: &[Source], reports: &[Vec<String>]) -> Vec<String> {
    let entry = src
        .iter()
        .find(|f| f.rel == "cli/main_entry.rs")
        .expect("cli/main_entry.rs");
    let text = String::from_utf8_lossy(&entry.masked);
    let rule = "cli/main_entry.rs: the function that calls crate::utils::xexit::write() must end \
                with its `match` on the command, then `crate::utils::xexit::write(); Ok(())`: a \
                report of the run's end belongs in xexit::write, where it runs only when asked for";
    let calls: Vec<usize> = text
        .match_indices("crate::utils::xexit::write()")
        .map(|(at, _)| at)
        .collect();
    let mut failures: Vec<String> = reports
        .iter()
        .map(|report| report.join("::"))
        .filter(|report| text.contains(report.as_str()))
        .map(|report| format!("cli/main_entry.rs calls {report} outside xexit::write"))
        .collect();
    let [at] = calls[..] else {
        failures.push(format!("{rule} ({} calls of it)", calls.len()));
        return failures;
    };
    let ends = entry
        .fns
        .iter()
        .filter(|f| f.1 < at && at < f.2)
        .max_by_key(|f| f.1)
        .is_some_and(|f| {
            let tail: String = text[at..=f.2].split_whitespace().collect();
            let before = text[..at].trim_end();
            before.ends_with('}') && tail == "crate::utils::xexit::write();Ok(())}" && {
                let open = opener(&entry.masked, before.len() - 1);
                let statement = text[..open].rfind([';', '{', '}']).map_or(0, |b| b + 1);
                text[statement..open].trim_start().starts_with("match ")
            }
        });
    if !ends {
        failures.push(rule.to_string());
    }
    failures
}

#[test]
fn every_report_of_the_runs_end_tests_a_listed_variable() {
    let src = sources();
    let file = |rel: &str| {
        src.iter()
            .position(|f| f.rel == rel)
            .unwrap_or_else(|| panic!("src/{rel}"))
    };
    let (k, x) = (file("utils/xexit.rs"), file("utils/xdiag.rs"));
    let list = |j: usize| {
        quoted_after(
            &String::from_utf8_lossy(&src[j].code),
            "const NAMES: &[&str] = &[",
        )
    };
    let (names, xdiag_names) = (list(k), list(x));
    assert!(
        names.windows(2).all(|w| w[0] < w[1]),
        "NAMES is not sorted: {names:?}"
    );
    let switches = switches(&src);
    let write = src[k]
        .fns
        .iter()
        .find(|f| f.0 == "write")
        .expect("xexit::write");
    let reports: Vec<Vec<String>> = calls_at(&src[k].masked, write.1, write.2)
        .into_iter()
        .map(|(_, path)| path)
        .filter(|path| path[0] == "crate")
        .collect();
    assert!(reports.len() > 30, "the calls of xexit::write: {reports:?}");
    let mut tests = Tests {
        src: &src,
        names: &names,
        known: HashMap::new(),
    };
    let (mut failures, mut used) = (Vec::new(), vec![false; names.len()]);
    for report in &reports {
        let shown = report.join("::");
        let found = resolve(&src, k, report, 0);
        if found.is_empty() {
            failures.push(format!(
                "{shown}: its function was not found (see the module doc)"
            ));
            continue;
        }
        let empty = found.iter().all(|&(j, _, e)| {
            let body = &src[j].masked[opener(&src[j].masked, e)..=e];
            body[1..body.len() - 1].iter().all(u8::is_ascii_whitespace)
        });
        if !empty && !found.iter().any(|&(j, s, e)| tests.function(j, s, e)) {
            failures.push(format!(
                "{shown}: no variable of xexit::NAMES (or of xdiag) is tested on its way to writing"
            ));
        }
        let mut queue: Vec<(usize, usize, usize, bool)> =
            found.iter().map(|&(j, s, e)| (j, s, e, false)).collect();
        let mut seen = Vec::new();
        while let Some((j, s, e, entry)) = queue.pop() {
            if seen.contains(&(j, s, entry)) {
                continue;
            }
            seen.push((j, s, entry));
            let (rel, code) = (&src[j].rel, &src[j].code[s..=e]);
            for (_, name) in reads(code) {
                if let Some(index) = names.iter().position(|n| *n == name) {
                    used[index] = true;
                }
            }
            for name in name_literals(code) {
                if !names.contains(&name)
                    && !xdiag_names.contains(&name)
                    && !switches.contains(&name)
                {
                    failures.push(format!(
                        "{shown}: src/{rel} names {name}, which is not in xexit::NAMES (a report \
                         written under it must be asked for there)"
                    ));
                }
            }
            for arg in unquoted_reads(code) {
                if !bound(code, &arg) {
                    failures.push(format!(
                        "{shown}: src/{rel} reads a variable through `{arg}`, not a quoted name, \
                         so which one is not known here (quote the name in the call)"
                    ));
                }
            }
            let guards = tests.guards(j, s, e, entry);
            for at in writes(&src[j].masked, s, e) {
                if !guards[at - s] {
                    failures.push(format!(
                        "{shown}: src/{rel}:{} writes a line that no test of a variable of \
                         xexit::NAMES (or of xdiag) guards; a counter or a switch alone is not one",
                        line(&src[j].code, at)
                    ));
                }
            }
            for (at, path) in calls_at(&src[j].masked, s, e) {
                for (a, b, c) in resolve(&src, j, &path, 0) {
                    queue.push((a, b, c, guards[at - s]));
                }
            }
        }
    }
    for (name, used) in names.iter().zip(&used) {
        if !used {
            failures.push(format!("{name} is in xexit::NAMES but no report reads it"));
        }
    }
    failures.extend(entry_ends_with_write(&src, &reports));
    failures.dedup();
    assert!(failures.is_empty(), "{failures:#?}");
}

fn fasta(name: &str) -> std::path::PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/fasta")
        .join(name)
}

/// A blastx search (30,000 bases of AP027131 against the first proteins of AP027133, one thread,
/// `-outfmt 6`), with `set` as its only `LOSAT_*` variable besides `LOSAT_BLASTX_WORKER_LOG`
/// (which a `--all-features` binary reads and which selects nothing).
fn search(dir: &Path, set: Option<(&str, &str)>) -> Output {
    let mut command = Command::new(env!("CARGO_BIN_EXE_losat"));
    for (name, _) in std::env::vars_os() {
        if name.to_string_lossy().starts_with("LOSAT_") && name != "LOSAT_BLASTX_WORKER_LOG" {
            command.env_remove(name);
        }
    }
    command.envs(set);
    command
        .current_dir(dir)
        .arg("blastx")
        .arg("-query")
        .arg(dir.join("query.fa"))
        .arg("-subject")
        .arg(dir.join("subject.faa"))
        .args([
            "-query_gencode",
            "4",
            "-outfmt",
            "6",
            "-evalue",
            "10",
            "-num_threads",
            "1",
        ])
        .output()
        .expect("run LOSAT")
}

/// Without a variable of `NAMES` no report line is written; with each alone the search reports the
/// same and its report is written (`LOSAT_X_REACH` with a value that is not an absolute path says
/// so, `utils/xreach.rs`).
#[test]
fn each_listed_variable_alone_still_writes_its_report() {
    let dir = std::env::temp_dir().join(format!("losat-exit-reports-{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).expect("directory");
    let query = std::fs::read_to_string(fasta("AP027131.fasta")).expect("AP027131.fasta");
    let mut lines = query.lines();
    let header = lines.next().expect("a header");
    let bases: String = lines.collect::<String>().chars().take(30_000).collect();
    std::fs::write(dir.join("query.fa"), format!("{header}\n{bases}\n")).expect("query");
    let proteins = std::fs::read_to_string(fasta("AP027133.faa")).expect("AP027133.faa");
    let subject: String = proteins
        .lines()
        .take(150)
        .map(|l| format!("{l}\n"))
        .collect();
    std::fs::write(dir.join("subject.faa"), subject).expect("subject");
    let plain = search(&dir, None);
    let plain_err = String::from_utf8_lossy(&plain.stderr).into_owned();
    assert!(plain.status.success(), "{plain_err}");
    assert!(!plain.stdout.is_empty(), "no hits");
    assert!(
        !plain_err
            .lines()
            .any(|l| l.starts_with("[X_") || l.starts_with("[TIMING]")),
        "a report without its variable: {plain_err}"
    );
    for (name, value, line) in [
        ("LOSAT_TIMING", "1", "[TIMING] common: since_start "),
        (
            "LOSAT_X_BNTBRANGELOG",
            "1",
            "[X_BNTBRANGE] sequences in ranges ",
        ),
        (
            "LOSAT_X_BXDIAGSHADOW",
            "1",
            "[X_SHADOW] BXDIAG subjects_compared=",
        ),
        (
            "LOSAT_X_BXSUBJFLOWSHADOW",
            "1",
            "[X_SHADOW] BXSUBJFLOW tasks_compared=",
        ),
        (
            "LOSAT_X_BXUNGAPSENTSHADOW",
            "1",
            "[X_SHADOW] BXUNGAPSENT subjects=",
        ),
        (
            "LOSAT_X_REACH",
            "0",
            "[X_REACH] ignored: LOSAT_X_REACH=0 is not an absolute path",
        ),
        ("LOSAT_X_STATS", "1", "[X_STATS] DP_SO_CALLS="),
    ] {
        let run = search(&dir, Some((name, value)));
        let stderr = String::from_utf8_lossy(&run.stderr);
        assert!(run.status.success(), "{name}: {stderr}");
        assert!(run.stdout == plain.stdout, "{name}: the report differs");
        assert!(
            stderr.lines().any(|l| l.starts_with(line)),
            "{name}: no `{line}` line: {stderr}"
        );
    }
    let _ = std::fs::remove_dir_all(&dir);
}
