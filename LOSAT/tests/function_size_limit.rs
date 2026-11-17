//! Caps the length of LOSAT's longest functions. Every function of `LOSAT/src` outside test code is
//! measured from the line of its `fn` keyword to the line of its closing brace; none may be longer
//! than `LIMIT` lines but those listed in `CAPS`, and those may not grow past their caps.
//!
//! The source is read as text. Comments, string literals (plain, byte, C, raw) and character and
//! byte literals are blanked first, newlines kept, so their braces do not count; a lifetime or a
//! label (`'a`) is not a literal. Then each `fn` item's parameter list and body braces are matched;
//! a declaration without a body has no length. Closures and nested functions count in the function
//! around them, and a nested function is also measured on its own. Test code is left out:
//! `#[cfg(test)]` modules, the files they declare (with or without `#[path]`, and the modules those
//! files declare), files under a `tests` directory, files with `#![cfg(test)]`, and items with a
//! `#[test]` attribute or a `#[cfg(..)]` that holds only under `test` (`all(test, ..)`; not
//! `any(test, ..)` or `not(test)`). These are the rules, and the counts, of the scanner of line
//! SW's audit (`rust_fns.py`, outside the repository).
//!
//! No NCBI counterpart: the test reads LOSAT's own source; NCBI BLAST+ has no such check.
//! Standard library only.

use std::collections::BTreeMap;
use std::fs;
use std::path::Path;

/// The functions over 1,000 lines, each capped at its length when this test was added:
/// (path below `LOSAT/src`, function name, cap in lines). The Owner decided on 2026-11-17 (answer a
/// to line SW's question Q12) that these functions may not grow and that no other function may
/// pass 1,000 lines; line SW's pull requests lower a cap when they shorten its function. When the
/// test fails, move the new code into a function of its own (the preferred way: one owner per
/// concept, AGENTS.md), or, when the growth is intended, raise that function's number here in the
/// same pull request. A new function over 1,000 lines is split, or gets a row here with a comment
/// that gives the reason. A name that occurs twice in one file has the cap on each.
#[rustfmt::skip]
const CAPS: &[(&str, &str, usize)] = &[
    ("algorithm/blastn/blast_engine/run.rs", "run_in_pool", 1516),
    ("algorithm/blastp/blast_engine.rs", "run_call", 1794),
    ("algorithm/tblastn/tblastn_app.rs", "search", 1451),
    ("algorithm/tblastx/sum_stats_linking/linking.rs", "link_hsp_group_ncbi", 1422),
];

/// The most lines a function without a row in `CAPS` may have.
const LIMIT: usize = 1000;

#[test]
fn functions_stay_within_their_caps() {
    let src = Path::new(concat!(env!("CARGO_MANIFEST_DIR"), "/src"));
    let offenders = offenders(&scan_tree(src), CAPS);
    assert!(
        offenders.is_empty(),
        "functions longer than allowed:\n{}\nTo fix: move the new code into a function of its own \
         (preferred), or, when the growth is intended, raise the function's cap in CAPS of \
         LOSAT/tests/function_size_limit.rs in the same pull request (a new function over the \
         limit: split it, or add a row with a reason)",
        offenders.join("\n")
    );
}

/// Each function outside test code that is longer than its cap, or than `LIMIT` without one.
fn offenders(functions: &[Function], caps: &[(&str, &str, usize)]) -> Vec<String> {
    let mut found = Vec::new();
    for f in functions.iter().filter(|f| !f.test) {
        let (most, bound) = match caps.iter().find(|c| c.0 == f.file && c.1 == f.name) {
            Some(&(_, _, cap)) => (cap, format!("cap {cap}")),
            None => (LIMIT, format!("limit {LIMIT}")),
        };
        if f.lines > most {
            let (file, line, name, lines) = (&f.file, f.line, &f.name, f.lines);
            found.push(format!("src/{file}:{line} {name} {lines} lines ({bound})"));
        }
    }
    found
}

/// A function with a body: its file below `src` (with `/`), its name, the line of its `fn`
/// keyword, its length in lines, and whether it is test code.
struct Function {
    file: String,
    name: String,
    line: usize,
    lines: usize,
    test: bool,
}

fn scan_tree(src: &Path) -> Vec<Function> {
    let mut files = Vec::new();
    read_rust_files(src, "", &mut files);
    files.sort();
    scan_sources(&files)
}

/// The `.rs` files below `dir`: (path below `src`, text).
fn read_rust_files(dir: &Path, rel: &str, out: &mut Vec<(String, Vec<u8>)>) {
    let entries = fs::read_dir(dir).unwrap_or_else(|e| panic!("{}: {e}", dir.display()));
    for entry in entries.map(Result::unwrap) {
        let (path, rel) = (
            entry.path(),
            format!("{rel}{}", entry.file_name().to_string_lossy()),
        );
        if entry.file_type().unwrap().is_dir() {
            read_rust_files(&path, &format!("{rel}/"), out);
        } else if rel.ends_with(".rs") {
            out.push((rel, fs::read(&path).unwrap()));
        }
    }
}

/// The functions of the given files (path below `src`, text), with test code marked.
fn scan_sources(files: &[(String, Vec<u8>)]) -> Vec<Function> {
    let scans: Vec<FileScan> = files.iter().map(|(_, text)| scan_file(text)).collect();
    // Test files: under a `tests` directory, with `#![cfg(test)]`, declared with a test attribute
    // or inside a test item, and every file that a test file declares.
    let mut queue = Vec::new();
    for (k, scan) in scans.iter().enumerate() {
        if scan.whole_test
            || files[k]
                .0
                .split('/')
                .rev()
                .skip(1)
                .any(|dir| dir == "tests")
        {
            queue.push(k);
        }
        for decl in scan.decls.values() {
            if decl.test || in_spans(&scan.spans, decl.pos) {
                queue.extend(resolve(files, k, decl));
            }
        }
    }
    let mut test_file = vec![false; files.len()];
    while let Some(k) = queue.pop() {
        if !std::mem::replace(&mut test_file[k], true) {
            for decl in scans[k].decls.values() {
                queue.extend(resolve(files, k, decl));
            }
        }
    }
    let mut functions = Vec::new();
    for (k, scan) in scans.into_iter().enumerate() {
        for (pos, name, line, lines) in scan.fns {
            let (file, test) = (
                files[k].0.clone(),
                test_file[k] || in_spans(&scan.spans, pos),
            );
            functions.push(Function {
                file,
                name,
                line,
                lines,
                test,
            });
        }
    }
    functions
}

/// One file: its functions with a body (offset of `fn`, name, line, lines), its test items (from
/// the test attribute to the end of the item), its `mod name;` declarations by the offset after
/// their `;`, and whether it has `#![cfg(test)]`.
#[derive(Default)]
struct FileScan {
    fns: Vec<(usize, String, usize, usize)>,
    spans: Vec<(usize, usize)>,
    decls: BTreeMap<usize, Decl>,
    whole_test: bool,
}

/// A `mod name;` declaration: the offset of `mod`, the name, its `#[path]`, and whether a test
/// attribute is on it.
struct Decl {
    pos: usize,
    name: String,
    path: Option<String>,
    test: bool,
}

fn in_spans(spans: &[(usize, usize)], pos: usize) -> bool {
    spans.iter().any(|&(a, b)| a <= pos && pos < b)
}

fn scan_file(src: &[u8]) -> FileScan {
    let m = mask(src);
    let close = closers(&m);
    let newlines: Vec<usize> = (0..src.len()).filter(|&i| src[i] == b'\n').collect();
    let line = |pos: usize| newlines.partition_point(|&nl| nl < pos) + 1;
    let mut scan = FileScan::default();
    for i in (0..m.len()).filter(|&i| !word_before(&m, i)) {
        if m[i..].starts_with(b"fn") {
            if let Some((name, end)) = fn_item(&m, &close, i) {
                scan.fns
                    .push((i, string(name), line(i), line(end) - line(i) + 1));
            }
        } else if m[i..].starts_with(b"mod") {
            let k = skip_ws(&m, i + 3);
            let name = string(ident(&m, k));
            let semi = skip_ws(&m, k + name.len());
            if k > i + 3 && !name.is_empty() && at(&m, semi) == b';' {
                let (pos, path, test) = (i, None, false);
                scan.decls.insert(
                    semi + 1,
                    Decl {
                        pos,
                        name,
                        path,
                        test,
                    },
                );
            }
        }
    }
    for a in 0..m.len() {
        let Some((mut end, inner, attr)) = attribute(src, &m, &close, a) else {
            continue;
        };
        if inner {
            scan.whole_test |= cfg_test(&attr);
            continue;
        }
        // The attributes that follow this one (the last `#[path]` counts), then the item.
        let mut path = path_of(&attr);
        while let Some((next_end, _, next)) = attribute(src, &m, &close, skip_ws(&m, end)) {
            (path, end) = (path_of(&next).or(path), next_end);
        }
        let test = cfg_test(&attr) || attr == b"test";
        let Some(end) = item_extent(&m, &close, skip_ws(&m, end)) else {
            continue;
        };
        if let Some(decl) = scan.decls.get_mut(&end) {
            decl.test |= test;
            decl.path = decl.path.take().or(path);
        } else if test {
            scan.spans.push((a, end));
        }
    }
    scan
}

/// The files that a `mod` declaration of file `k` loads.
fn resolve(files: &[(String, Vec<u8>)], k: usize, decl: &Decl) -> Vec<usize> {
    let rel = &files[k].0;
    let dir = &rel[..rel.rfind('/').map_or(0, |p| p + 1)];
    let own = match &rel[dir.len()..] {
        "mod.rs" | "lib.rs" | "main.rs" => dir.to_string(),
        base => format!("{dir}{}/", base.strip_suffix(".rs").unwrap_or(base)),
    };
    let name = &decl.name;
    let candidates = match &decl.path {
        Some(path) => vec![normalize(&format!("{dir}{path}"))],
        None => vec![format!("{own}{name}.rs"), format!("{own}{name}/mod.rs")],
    };
    let index = |c: &String| files.iter().position(|f| f.0 == *c);
    candidates.iter().filter_map(index).collect()
}

fn normalize(path: &str) -> String {
    let mut parts: Vec<&str> = Vec::new();
    for part in path.split('/').filter(|&p| p != "." && !p.is_empty()) {
        if part == ".." && parts.last().is_some_and(|&p| p != "..") {
            parts.pop();
        } else {
            parts.push(part);
        }
    }
    parts.join("/")
}

/// `src` with every comment and literal turned into spaces, newlines kept.
fn mask(src: &[u8]) -> Vec<u8> {
    let (mut m, mut i) = (src.to_vec(), 0);
    while i < src.len() {
        let Some(end) = literal_end(src, i) else {
            i += 1;
            continue;
        };
        m[i..end]
            .iter_mut()
            .filter(|b| **b != b'\n')
            .for_each(|b| *b = b' ');
        i = end;
    }
    m
}

/// The end of the comment or literal that starts at `i`, tried in this order: line comment, block
/// comment (nested), raw string (`r`, `br`, `cr`, any number of `#`), string (plain, `b`, `c`),
/// character or byte literal. A literal's prefix must not end an identifier.
fn literal_end(s: &[u8], i: usize) -> Option<usize> {
    let (rest, n) = (&s[i..], s.len());
    let prefix_ok = i == 0 || !(s[i - 1].is_ascii_alphanumeric() || s[i - 1] == b'_');
    if rest.starts_with(b"//") {
        return Some(find(s, i, b"\n").unwrap_or(n));
    }
    if rest.starts_with(b"/*") {
        let (mut depth, mut k) = (1, i + 2);
        while depth > 0 && k < n {
            let Some(c) = find(s, k, b"*/") else {
                return Some(n);
            };
            (depth, k) = match find(&s[..c + 1], k, b"/*") {
                Some(o) => (depth + 1, o + 2),
                None => (depth - 1, c + 2),
            };
        }
        return Some(k);
    }
    let raw = 1 + usize::from(rest.starts_with(b"br") || rest.starts_with(b"cr"));
    if prefix_ok && (raw == 2 || rest[0] == b'r') {
        let q = i + raw + rest[raw..].iter().take_while(|&&b| b == b'#').count();
        if at(s, q) == b'"' {
            let mut term = vec![b'"'];
            term.resize(q - i - raw + 1, b'#');
            return Some(find(s, q + 1, &term).map_or(n, |e| e + term.len()));
        }
    }
    let prefix = |p: u8, quote: u8| prefix_ok && rest[0] == p && at(s, i + 1) == quote;
    let mut k = i + usize::from(prefix(b'b', b'"') || prefix(b'c', b'"'));
    if at(s, k) == b'"' {
        while k + 1 < n {
            k += 1;
            match s[k] {
                b'"' => return Some(k + 1),
                b'\\' => k += 1,
                _ => {}
            }
        }
    }
    // One character, or an escape (`\'`, `\x41`, `\u{7d}`), between single quotes.
    let q = i + usize::from(prefix(b'b', b'\''));
    match (at(s, q), at(s, q + 1)) {
        (b'\'', b'\\') => find(s, q + 3, b"'").map(|e| e + 1),
        (b'\'', b'\'' | b'\n' | 0) => None,
        (b'\'', c) => {
            let e = q + 1 + c.leading_ones().max(1) as usize;
            (at(s, e) == b'\'').then_some(e + 1)
        }
        _ => None,
    }
}

fn find(s: &[u8], from: usize, pattern: &[u8]) -> Option<usize> {
    let mut windows = s.get(from..)?.windows(pattern.len());
    Some(from + windows.position(|w| w == pattern)?)
}

/// For each `{`, `(` and `[` of `m`, the offset of the one that closes it (`u32::MAX`: none).
fn closers(m: &[u8]) -> Vec<u32> {
    let (mut close, mut open) = (vec![u32::MAX; m.len()], [vec![], vec![], vec![]]);
    for (i, b) in m.iter().enumerate() {
        if let Some(kind) = b"{([".iter().position(|o| o == b) {
            open[kind].push(i);
        } else if let Some(kind) = b"})]".iter().position(|c| c == b) {
            if let Some(o) = open[kind].pop() {
                close[o] = i as u32;
            }
        }
    }
    close
}

fn closer(close: &[u32], i: usize) -> Option<usize> {
    let c = *close.get(i)?;
    (c != u32::MAX).then_some(c as usize)
}

fn at(m: &[u8], i: usize) -> u8 {
    m.get(i).copied().unwrap_or(0)
}

fn is_space(b: u8) -> bool {
    matches!(b, b' ' | b'\t'..=b'\r')
}

fn skip_ws(m: &[u8], mut i: usize) -> usize {
    while is_space(at(m, i)) {
        i += 1;
    }
    i
}

fn word_before(m: &[u8], i: usize) -> bool {
    i > 0 && (m[i - 1].is_ascii_alphanumeric() || m[i - 1] == b'_' || m[i - 1] >= 0x80)
}

/// The identifier at `i` (empty when none starts there).
fn ident(m: &[u8], i: usize) -> &[u8] {
    let rest = m.get(i..).unwrap_or(&[]);
    if rest.first().is_some_and(u8::is_ascii_digit) {
        return &[];
    }
    let word = |b: &u8| b.is_ascii_alphanumeric() || *b == b'_';
    &rest[..rest.iter().take_while(|b| word(b)).count()]
}

fn string(bytes: &[u8]) -> String {
    String::from_utf8_lossy(bytes).into_owned()
}

/// The `fn` item whose keyword is at `i`: its name and the offset of its body's closing brace.
fn fn_item<'m>(m: &'m [u8], close: &[u32], i: usize) -> Option<(&'m [u8], usize)> {
    let j = skip_ws(m, i + 2);
    let name = ident(m, j);
    if j == i + 2 || name.is_empty() {
        return None;
    }
    let mut k = skip_ws(m, j + name.len());
    if at(m, k) == b'<' {
        k = skip_ws(m, skip_angle(m, k)?);
    }
    if at(m, k) != b'(' {
        return None;
    }
    // The body: the first `{` outside parentheses, brackets and angle brackets; a `;` before it
    // (outside parentheses and brackets) or a `}` means a declaration without a body.
    let (mut depth, mut angle, mut p) = (0, 0, closer(close, k)? + 1);
    while p < m.len() {
        match m[p] {
            b'-' | b'=' if at(m, p + 1) == b'>' => p += 1,
            b'(' | b'[' => depth += 1,
            b')' | b']' => depth -= 1,
            b'<' => angle += 1,
            b'>' => angle -= 1,
            b';' if depth <= 0 => return None,
            b'{' if depth <= 0 && angle <= 0 => return Some((name, closer(close, p)?)),
            b'}' => return None,
            _ => {}
        }
        p += 1;
    }
    None
}

/// The offset after the `>` that closes the generic parameters opened at `i` (looked for in the
/// next 20,000 bytes, before any `{`, `}` or `;`).
fn skip_angle(m: &[u8], mut i: usize) -> Option<usize> {
    let (end, mut depth) = (m.len().min(i + 20_000), 0);
    while i < end {
        match m[i] {
            b'-' | b'=' if at(m, i + 1) == b'>' => i += 1,
            b'<' => depth += 1,
            b'>' if depth == 1 => return Some(i + 1),
            b'>' => depth -= 1,
            b'{' | b'}' | b';' => return None,
            _ => {}
        }
        i += 1;
    }
    None
}

/// The attribute whose `#` is at `a`: the offset after its `]`, whether it is inner (`#!`), and
/// its text inside the brackets without white space.
fn attribute(src: &[u8], m: &[u8], close: &[u32], a: usize) -> Option<(usize, bool, Vec<u8>)> {
    if at(m, a) != b'#' {
        return None;
    }
    let mut k = skip_ws(m, a + 1);
    let inner = at(m, k) == b'!';
    k = if inner { skip_ws(m, k + 1) } else { k };
    let rb = closer(close, k).filter(|_| at(m, k) == b'[')?;
    let mut text = src[k + 1..rb].to_vec();
    text.retain(|&b| !is_space(b));
    Some((rb + 1, inner, text))
}

/// Whether the attribute is a `cfg(..)` that holds only when `test` is set.
fn cfg_test(attr: &[u8]) -> bool {
    let Some(predicate) = attr.strip_prefix(b"cfg(") else {
        return false;
    };
    predicate
        .strip_suffix(b")")
        .is_some_and(|e| requires_test(e, &mut 0))
}

/// Whether the cfg predicate at `e[*i..]` needs `test`: `test`, `all(..)` with one such part,
/// `any(..)` with only such parts; not `not(..)`, `name = ".."` or another name. Moves `*i` past it.
fn requires_test(e: &[u8], i: &mut usize) -> bool {
    let name = ident(e, *i);
    *i += name.len();
    if name.is_empty() || at(e, *i) == b'=' {
        while *i < e.len() && !matches!(e[*i], b',' | b')') {
            *i += 1;
        }
        return false;
    }
    if at(e, *i) != b'(' {
        return name == b"test";
    }
    let mut parts = Vec::new();
    *i += 1;
    while *i < e.len() && e[*i] != b')' {
        parts.push(requires_test(e, i));
        *i += usize::from(at(e, *i) == b',');
    }
    *i += 1;
    match name {
        b"all" => parts.contains(&true),
        b"any" => !parts.is_empty() && !parts.contains(&false),
        _ => false,
    }
}

fn path_of(attr: &[u8]) -> Option<String> {
    let path = attr.strip_prefix(b"path=\"")?.strip_suffix(b"\"")?;
    (!path.contains(&b'"')).then(|| string(path))
}

const QUALIFIERS: &str = "pub default unsafe async extern const";
const KINDS: &str = "mod fn impl trait struct enum union use static const type macro_rules crate";

fn listed(list: &str, word: &[u8]) -> bool {
    !word.is_empty() && list.split(' ').any(|w| w.as_bytes() == word)
}

/// The end of the item that starts at `start` (after its attributes): after its `;`, or after the
/// brace that closes its first brace outside parentheses and brackets (`use`, `static`, `const`,
/// `type` items: after the first `;` outside all brackets). `None` where no item starts: a
/// statement, a field, a match arm.
fn item_extent(m: &[u8], close: &[u32], start: usize) -> Option<usize> {
    // Qualifiers (`pub(..)` first), then the kind of item; `const` is either.
    let (mut i, mut kind) = (start, None);
    loop {
        let word = ident(m, i);
        let mut end = i + word.len();
        if word == b"pub" && at(m, skip_ws(m, end)) == b'(' {
            end = closer(close, skip_ws(m, end))? + 1;
        }
        if listed(KINDS, word) {
            kind = Some((word, end));
        }
        let next = skip_ws(m, end);
        if !listed(QUALIFIERS, word) || next == end {
            break;
        }
        i = next;
    }
    let (word, bang) = (ident(m, start), skip_ws(m, start + ident(m, start).len()));
    let opens = b"{([".contains(&at(m, skip_ws(m, bang + 1)));
    let (kind, from) = match kind {
        Some(kind) => kind,
        None if !word.is_empty() && at(m, bang) == b'!' && opens => (&b"macro"[..], start),
        None => return None,
    };
    let semi = listed("use static const type crate", kind);
    let (first, mut depth) = (if semi { start } else { from }, 0);
    for k in first..m.len() {
        match m[k] {
            b'(' | b'[' => depth += 1,
            b')' | b']' => depth -= 1,
            b'{' if semi => depth += 1,
            b'}' if semi => depth -= 1,
            b';' if depth <= 0 => return Some(k + 1),
            b'{' if depth <= 0 => return closer(close, k).map(|e| e + 1),
            _ => {}
        }
    }
    None
}

/// `file:line name lines`, and ` test` for test code, of each function of the given files.
fn found(files: &[(&str, &str)]) -> Vec<String> {
    let own = |&(rel, text): &(&str, &str)| (rel.to_string(), text.as_bytes().to_vec());
    let files: Vec<_> = files.iter().map(own).collect();
    let mut rows = Vec::new();
    for f in scan_sources(&files) {
        let (file, line, name, lines) = (f.file, f.line, f.name, f.lines);
        let test = if f.test { " test" } else { "" };
        rows.push(format!("{file}:{line} {name} {lines}{test}"));
    }
    rows
}

const LITERALS: &str = r##"fn a<'a>(x: &'a str) -> &'a str {
    // }
    /* } /* } */ } */
    let s = "}\"}";
    let r = r#"}"}"#;
    let c = ('}', b'{', '\'', b"{", '\u{7d}');
    'b: { break 'b; }
    x
}
fn b() {}"##;

#[test]
fn braces_in_comments_and_literals_do_not_count() {
    let expected = ["lib.rs:1 a 9", "lib.rs:10 b 1"];
    assert_eq!(found(&[("lib.rs", LITERALS)]), expected);
}

const ITEMS: &str = r#"trait T {
    fn declared(&self);
    fn provided(&self) {}
}
pub(crate) const unsafe fn outer() {
    fn inner() {
    }
}
#[test]
#[ignore = "slow"]
fn unit() {}
#[cfg(all(test, unix))]
fn helper() {}
#[cfg(any(test, unix))]
fn shared() {}
#[cfg(test)]
mod tests {
    fn in_module() {}
}
#[cfg(test)]
#[path = "../x_tests.rs"]
mod x;"#;

#[test]
fn nested_functions_and_test_code() {
    let files = [
        ("a/lib.rs", ITEMS),
        ("x_tests.rs", "mod deep;\nfn in_file() {}"),
        ("x_tests/deep.rs", "fn deeper() {}"),
        ("writedb/tests/b.rs", "fn in_dir() {}"),
        ("plain.rs", "#[cfg(not(test))]\nfn plain() {}"),
    ];
    let expected = [
        "a/lib.rs:3 provided 1",
        "a/lib.rs:5 outer 4",
        "a/lib.rs:6 inner 2",
        "a/lib.rs:11 unit 1 test",
        "a/lib.rs:13 helper 1 test",
        "a/lib.rs:15 shared 1",
        "a/lib.rs:18 in_module 1 test",
        "x_tests.rs:2 in_file 1 test",
        "x_tests/deep.rs:1 deeper 1 test",
        "writedb/tests/b.rs:1 in_dir 1 test",
        "plain.rs:2 plain 1",
    ];
    assert_eq!(found(&files), expected);
}

#[test]
fn offenders_name_their_cap_or_the_limit() {
    let body = "\n".repeat(LIMIT);
    let text = format!("fn long() {{{body}}}\nfn capped() {{\n}}\n#[test]\nfn unit() {{{body}}}");
    let functions = scan_sources(&[("a.rs".to_string(), text.into_bytes())]);
    let expected = [
        "src/a.rs:1 long 1001 lines (limit 1000)",
        "src/a.rs:1002 capped 2 lines (cap 1)",
    ];
    assert_eq!(offenders(&functions, &[("a.rs", "capped", 1)]), expected);
    assert!(offenders(&functions, &[("a.rs", "long", 1001), ("a.rs", "capped", 2)]).is_empty());
}
