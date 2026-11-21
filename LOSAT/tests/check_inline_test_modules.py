#!/usr/bin/env python3
"""Keep the unit tests out of the published `losat` crate: no inline `#[cfg(test)] mod x { ... }`.

No NCBI counterpart: this checks LOSAT's own packaging, not ported behavior.

crates.io refuses a `.crate` file over 10 MiB (the "Crate size" step of ci.yml). `exclude` in
LOSAT/Cargo.toml leaves the files named `*_tests.rs`, `tests.rs`, `tests_*.rs` out of the package
(`check_crate_test_files.py` guards that), but a test module written inline in a source file
(`#[cfg(test)] mod tests { ... }`) is packaged with it: on 2026-11-21 the 391 such modules were 1.29
MiB of the 6.20 MiB crate. A unit test lives in a sibling file instead (AGENTS.md):

    #[cfg(test)]
    #[path = "<stem>_tests.rs"]
    mod tests;

The module keeps its name and its place in the module tree, so the test names, `use super::*` and the
`include_str!` paths (relative to the file, whose directory does not change) stay as they were.

Check: fails, naming `file:line`, for every inline module of a packaged `LOSAT/src` file whose `cfg`
holds only in a test build (`#[cfg(test)]`, `#[cfg(all(test, ..))]`, or an inner `#![cfg(test)]`).
`cfg(any(test, ..))` and `cfg(not(test))` are not test-only and stay. A test-only module inside
another test-only one moves with it.

`--fix` (a branch that merges a dev with such modules runs it; it is idempotent) moves each such
module's body, unchanged but for one indentation level (the module's indentation and 4 spaces, taken
off the lines that do not start inside a string, raw string, byte string or block comment), into
`<stem>_tests.rs` beside the file (`<stem>_<name>_tests.rs` for a module `name`, `<stem>_<name>.rs`
when `name` ends in `_tests`; a `mod.rs` has no stem), and leaves the module's other attributes and
doc comments where they are, with `#[path = ".."]` added. The files are found by the patterns of
`exclude`. It reads the source with a lexer (nested block comments, raw strings, byte strings, char
literals against lifetimes), not with regular expressions, and refuses, changing nothing, when
- the file name is taken (an existing file, or another module of the run) and so is the one fallback,
  `<stem>_inline_tests.rs` for `<stem>_tests.rs` (2026-11-21: `args.rs` and `translation.rs` have a
  declared `args_tests.rs` as well as an inline `mod tests`);
- the module sits inside something that is not a `mod` block (a function, an `impl`, a macro);
- the body declares `mod x;` or `#[path]` (their files would be looked for elsewhere) or uses
  `file!`, `line!` or `column!` (their values would change);
- the module is not laid out as rustfmt writes a module (`{` ending its line, `}` alone on a line at
  the module's indentation), or a line of the body that starts outside a literal is indented less than
  the module's body, or the file is not read back to the text it came from (the self-check);
- a path of a `#[path]` is not decided by the Rust reference's rules (a `mod.rs`, a crate root, or a
  file that a `#[path]` loads, owns its directory; a module of another file inside an inline
  module is looked for under the file's name and the inline modules' names).
Run `cargo fmt` after `--fix` (the dedented lines may fit one more level of rustfmt's width), and
`check_foreign_quotes.py --fix`, which reads the moved files too. `--fix` also rewrites the
`losat_test` paths of docs/evidence/ncbi_unit_cases/LEDGER.tsv that name a test it moved.

Usage: python3 LOSAT/tests/check_inline_test_modules.py --root . [--fix]
"""

from __future__ import annotations

import argparse
import posixpath
import re
import sys
import tomllib
from dataclasses import dataclass
from pathlib import Path
from typing import Callable, Sequence

sys.path.insert(0, str(Path(__file__).resolve().parent))
import check_crate_test_files as crate  # noqa: E402  (the one owner of `exclude`, cfg and module-path rules)

LEDGER = "docs/evidence/ncbi_unit_cases/LEDGER.tsv"
INDENT = 4  # rustfmt's indentation level


# ---------------------------------------------------------------- the lexer
@dataclass(frozen=True)
class Token:
    kind: str  # ws, line, block, ident, life, str, char, num, punct
    start: int
    end: int


class LexError(ValueError):
    pass


RAW_STRING_RE = re.compile(r'(?:br|cr|r)(#*)"')
RAW_IDENT_RE = re.compile(r"r#[A-Za-z_]")


def ident_end(text: str, i: int) -> int:
    while i < len(text) and (text[i].isalnum() or text[i] == "_"):
        i += 1
    return i


def string_end(text: str, i: int) -> int:
    """The end of the string literal whose opening quote is at `i` (a backslash takes the next char)."""
    j = i + 1
    while j < len(text):
        if text[j] == "\\":
            j += 2
        elif text[j] == '"':
            return j + 1
        else:
            j += 1
    raise LexError(f"unterminated string at offset {i}")


def char_or_lifetime(text: str, i: int) -> Token:
    """`'x'`, `'\\n'` and `b'x'` are chars; `'a` and `'static` are lifetimes (and labels)."""
    if text.startswith("\\", i + 1):
        close = text.find("'", i + 3)
        if close < 0:
            raise LexError(f"unterminated char at offset {i}")
        return Token("char", i, close + 1)
    if i + 2 < len(text) and text[i + 2] == "'":
        return Token("char", i, i + 3)
    return Token("life", i, ident_end(text, i + 1))


def lex(text: str) -> list[Token]:
    """The tokens of Rust source, whitespace and comments included; they cover `text` without gaps."""
    tokens: list[Token] = []
    i, n = 0, len(text)
    while i < n:
        c = text[i]
        if c.isspace():
            j = i
            while j < n and text[j].isspace():
                j += 1
            tokens.append(Token("ws", i, j))
        elif text.startswith("//", i):
            j = text.find("\n", i)
            j = n if j < 0 else j
            tokens.append(Token("line", i, j))
        elif text.startswith("/*", i):
            depth, j = 1, i + 2
            while depth:
                inner, close = text.find("/*", j), text.find("*/", j)
                if close < 0:
                    raise LexError(f"unterminated block comment at offset {i}")
                if 0 <= inner < close:
                    depth, j = depth + 1, inner + 2
                else:
                    depth, j = depth - 1, close + 2
            tokens.append(Token("block", i, j))
        elif c == '"':
            tokens.append(Token("str", i, string_end(text, i)))
        elif c == "'":
            tokens.append(char_or_lifetime(text, i))
        elif c.isalpha() or c == "_":
            raw = RAW_STRING_RE.match(text, i)
            if raw:  # r"..", r#".."#, br"..", cr".."
                close = text.find('"' + raw.group(1), raw.end())
                if close < 0:
                    raise LexError(f"unterminated raw string at offset {i}")
                tokens.append(Token("str", i, close + 1 + len(raw.group(1))))
            elif text.startswith(("b\"", "c\""), i):
                tokens.append(Token("str", i, string_end(text, i + 1)))
            elif text.startswith("b'", i):
                close = char_or_lifetime(text, i + 1)
                if close.kind != "char":
                    raise LexError(f"bad byte literal at offset {i}")
                tokens.append(Token("char", i, close.end))
            elif RAW_IDENT_RE.match(text, i):
                tokens.append(Token("ident", i, ident_end(text, i + 2)))
            else:
                tokens.append(Token("ident", i, ident_end(text, i)))
        elif c.isdigit():
            j = i
            while j < n and (text[j].isalnum() or text[j] == "_"):
                j += 1
            tokens.append(Token("num", i, j))
        else:
            tokens.append(Token("punct", i, i + 1))
        i = tokens[-1].end
    return tokens


def literal_lines(text: str, tokens: Sequence[Token]) -> set[int]:
    """The 0-based numbers of the lines that start inside a string, byte string, raw string or block comment."""
    starts = [0] + [m.end() for m in re.finditer("\n", text)]
    inside: set[int] = set()
    for t in tokens:
        if t.kind in ("str", "block") and "\n" in text[t.start:t.end]:
            first = text.count("\n", 0, t.start) + 1  # the first line that starts after the token's start
            last = text.count("\n", 0, t.end - 1)  # the line of the token's last character
            inside.update(range(first, last + 1))
    return inside & set(range(len(starts)))


# ---------------------------------------------------------------- the items
@dataclass
class Found:
    name: str
    enclosing: tuple[str, ...]  # names of the inline modules around it
    nested_in_other: str | None  # what is around it when that is not a `mod` block (a refusal)
    item_start: int  # first token of `[pub(..)] mod`
    open: int  # offset of `{`
    close: int  # offset of `}`
    line: int  # 1-based line of `mod`


def significant(tokens: Sequence[Token]) -> list[int]:
    return [k for k, t in enumerate(tokens) if t.kind not in ("ws", "line", "block")]


def matching(tokens: Sequence[Token], text: str, sig: Sequence[int], at: int, step: int) -> int:
    """The index in `sig` of the bracket that closes (step 1) or opens (step -1) the one at `sig[at]`."""
    pairs = {"(": ")", "[": "]", "{": "}"}
    reverse = {v: k for k, v in pairs.items()}
    depth = 0
    k = at
    while 0 <= k < len(sig):
        ch = text[tokens[sig[k]].start]
        if tokens[sig[k]].kind == "punct":
            if (step == 1 and ch in pairs) or (step == -1 and ch in reverse):
                depth += 1
            elif (step == 1 and ch in reverse) or (step == -1 and ch in pairs):
                depth -= 1
                if depth == 0:
                    return k
        k += step
    raise LexError("unbalanced brackets")


def cfg_predicates(text: str, tokens: Sequence[Token], sig: Sequence[int], lo: int, hi: int) -> list[str]:
    """The `cfg(..)` predicates of the attribute whose tokens are `sig[lo:hi]` (`[` .. `]` included)."""
    inner = sig[lo + 1:hi - 1]
    if len(inner) < 3 or text[tokens[inner[0]].start:tokens[inner[0]].end] != "cfg":
        return []
    if text[tokens[inner[1]].start] != "(":
        return []
    return ["".join(text[tokens[k].start:tokens[k].end] for k in inner[2:-1])]


def find_modules(text: str) -> list[Found]:
    """The outermost test-only inline modules of a source file, in order."""
    tokens = lex(text)
    sig = significant(tokens)
    txt = lambda k: text[tokens[sig[k]].start:tokens[sig[k]].end]  # noqa: E731
    stack: list[dict] = []
    found: list[Found] = []
    for k in range(len(sig)):
        t = tokens[sig[k]]
        if t.kind != "punct":
            continue
        ch = text[t.start]
        if ch == "}":
            frame = stack.pop()
            if frame["found"] is not None:
                frame["found"].close = t.start
            continue
        if ch != "{":
            continue
        is_mod = k >= 2 and txt(k - 2) == "mod" and tokens[sig[k - 1]].kind == "ident"
        if not is_mod:
            first = k
            while first > 0 and txt(first - 1) not in (";", "{", "}"):
                first -= 1
            what = " ".join(txt(i) for i in range(first, k))
            stack.append({"kind": "other", "test": False, "found": None, "what": what[:40]})
            continue
        name = txt(k - 1)
        mod_at = k - 2
        # `pub`, `pub(crate)`, `pub(in path)` before `mod`
        item_at = mod_at
        if mod_at >= 1 and txt(mod_at - 1) == "pub":
            item_at = mod_at - 1
        elif mod_at >= 1 and txt(mod_at - 1) == ")":
            opener = matching(tokens, text, sig, mod_at - 1, -1)
            if opener >= 1 and txt(opener - 1) == "pub":
                item_at = opener - 1
        # outer attributes (and comments, which are not significant) before the item
        test = False
        at = item_at - 1
        while at >= 1 and txt(at) == "]":
            opener = matching(tokens, text, sig, at, -1)
            if opener < 1 or txt(opener - 1) != "#":
                break
            test = test or any(crate.test_only(p) for p in cfg_predicates(text, tokens, sig, opener, at + 1))
            at = opener - 2
        # an inner `#![cfg(test)]` first in the body
        if k + 3 < len(sig) and txt(k + 1) == "#" and txt(k + 2) == "!" and txt(k + 3) == "[":
            closer = matching(tokens, text, sig, k + 3, 1)
            test = test or any(crate.test_only(p) for p in cfg_predicates(text, tokens, sig, k + 3, closer + 1))
        within_test = any(f["test"] for f in stack)
        frame = {"kind": "mod", "test": test, "found": None, "name": name}
        if test and not within_test:
            others = [f["what"] for f in stack if f["kind"] == "other"]
            item = Found(
                name=name,
                enclosing=tuple(f["name"] for f in stack if f["kind"] == "mod"),
                nested_in_other=others[-1] if others else None,
                item_start=tokens[sig[item_at]].start,
                open=t.start,
                close=-1,
                line=text.count("\n", 0, tokens[sig[mod_at]].start) + 1,
            )
            frame["found"] = item
            found.append(item)
        stack.append(frame)
    if stack:
        raise LexError("unbalanced braces")
    return found


# ---------------------------------------------------------------- moving a module
class Refusal(Exception):
    pass


def sibling_name(stem: str, found: Found) -> str:
    parts = [p for p in (stem, *found.enclosing) if p]
    if found.name == "tests":
        parts.append("tests")
        return "_".join(parts) + ".rs"
    parts.append(found.name)
    return "_".join(parts) + (".rs" if found.name.endswith("_tests") else "_tests.rs")


def dedent(content: str, n: int) -> str:
    tokens = lex(content)
    inside = literal_lines(content, tokens)
    out = []
    for number, line in enumerate(content.splitlines(keepends=True)):
        if number not in inside and line.startswith(" " * n):
            line = line[n:]
        out.append(line)
    return "".join(out)


def reindent(content: str, n: int) -> str:
    tokens = lex(content)
    inside = literal_lines(content, tokens)
    out = []
    for number, line in enumerate(content.splitlines(keepends=True)):
        if number not in inside and line.strip("\n") != "":
            line = " " * n + line
        out.append(line)
    return "".join(out)


def check_body(content: str) -> None:
    """Refuse a body whose meaning depends on the file it is in."""
    tokens = lex(content)
    sig = significant(tokens)
    words = [content[tokens[k].start:tokens[k].end] for k in sig]
    for i, word in enumerate(words):
        nxt = words[i + 1] if i + 1 < len(words) else ""
        after = words[i + 2] if i + 2 < len(words) else ""
        if word == "mod" and i + 1 < len(words) and tokens[sig[i + 1]].kind == "ident" and after == ";":
            raise Refusal(f"its body declares `mod {nxt};`, whose file would be looked for elsewhere")
        if word == "path" and i >= 2 and words[i - 1] == "[" and words[i - 2] == "#" and nxt == "=":
            raise Refusal("its body has a `#[path]` attribute, which would be read from another directory")
        if word in ("file", "line", "column") and nxt == "!" and after in ("(", "[", "{"):
            raise Refusal(f"its body uses `{word}!`, whose value would change")


@dataclass
class Move:
    parent: str  # path below LOSAT/
    file: str  # path below LOSAT/ of the sibling
    attr: str  # the #[path] text
    found: Found
    content: str  # the sibling's text
    indent: str


def mod_rs_files(crate_dir: Path, files: set[str], roots: Sequence[str]) -> set[str]:
    """The files that own their directory: crate roots, `mod.rs`, and what a `#[path]` loads."""
    owners = set(roots) | {f for f in files if posixpath.basename(f) == "mod.rs"}
    parsed = {}
    changed = True
    while changed:
        changed = False
        for f in sorted(files):
            if not (f.startswith("src/") and f.endswith(".rs")):
                continue
            if f not in parsed:
                parsed[f] = crate.edges((crate_dir / f).read_text(encoding="utf-8", errors="replace"))[0]
            for edge in parsed[f]:
                if edge.kind == "mod" and edge.path_attr is not None:
                    target = crate.module_candidates(f, f in owners, edge)[0]
                    if target not in owners:
                        owners.add(target)
                        changed = True
    return owners


def plan_file(parent: str, text: str, found: Sequence[Found], owners: Callable[[], set[str]],
              occupied: set[str]) -> list[Move]:
    stem = posixpath.basename(parent)[:-3]
    if stem == "mod":
        stem = ""
    moves = []
    for f in found:
        where = f"LOSAT/{parent}:{f.line}: mod {f.name}"
        try:
            if f.nested_in_other is not None:
                raise Refusal(f"it sits inside `{f.nested_in_other}`, which is not a `mod` block")
            line_start = text.rfind("\n", 0, f.item_start) + 1
            before = text[line_start:f.item_start]
            indent = before if before.strip() == "" else re.match(r"[ \t]*", text[line_start:]).group(0)
            if "\t" in indent:
                raise Refusal("its line is indented with a tab")
            if not text.startswith("\n", f.open + 1):
                raise Refusal("its `{` does not end its line")
            body = text[f.open + 1:f.close]
            tail = "\n" + indent
            if not body.endswith(tail):
                raise Refusal("its `}` is not alone on a line at the module's indentation")
            content = body[1:len(body) - len(indent)]
            check_body(content)
            n = len(indent) + INDENT
            sibling = dedent(content, n)
            if reindent(sibling, n) != content:
                raise Refusal(f"a line of its body that starts outside a literal is indented less than {n} spaces")
            name = sibling_name(stem, f)
            directory = posixpath.dirname(parent)
            target = posixpath.join(directory, name)
            if target in occupied:  # the one fallback: `x_tests.rs` -> `x_inline_tests.rs`
                name = name[:-len("_tests.rs")] + "_inline_tests.rs"
                target = posixpath.join(directory, name)
                if target in occupied:
                    raise Refusal(f"the files LOSAT/{target} and the one before it are taken")
            occupied.add(target)
            if f.enclosing:
                base = posixpath.join(directory if parent in owners() else posixpath.join(directory, stem), *f.enclosing)
                rel = posixpath.relpath(target, base)
            else:
                rel = name
            moves.append(Move(parent, target, f'#[path = "{rel}"]', f, sibling, indent))
        except Refusal as why:
            raise Refusal(f"{where}: {why}") from None
    return moves


def apply_moves(text: str, moves: Sequence[Move]) -> str:
    for move in reversed(moves):
        f = move.found
        line_start = text.rfind("\n", 0, f.item_start) + 1
        own_line = text[line_start:f.item_start].strip() == ""
        head = text[f.item_start:f.open].rstrip()
        glue = "\n" + move.indent if own_line else " "
        text = text[:f.item_start] + move.attr + glue + head + ";" + text[f.close + 1:]
    return text


def restore(text: str, sibling: Callable[[str, str], str | None]) -> str:
    """The parent as it was before `apply_moves`, from its new text and the sibling files (the self-check)."""
    tokens = lex(text)
    sig = significant(tokens)
    word = lambda k: text[tokens[sig[k]].start:tokens[sig[k]].end]  # noqa: E731
    edits = []
    for k in range(2, len(sig)):
        if word(k) != ";" or word(k - 2) != "mod" or tokens[sig[k - 1]].kind != "ident":
            continue
        j = k - 2  # the first token of `[pub(..)] mod name;`
        if j >= 1 and word(j - 1) == "pub":
            j -= 1
        elif j >= 1 and word(j - 1) == ")":
            opener = matching(tokens, text, sig, j - 1, -1)
            if opener >= 1 and word(opener - 1) == "pub":
                j = opener - 1
        if j < 1 or word(j - 1) != "]":
            continue
        opener = matching(tokens, text, sig, j - 1, -1)
        attr = text[tokens[sig[opener - 1]].start:tokens[sig[j - 1]].end]  # `#` .. `]`
        m = re.fullmatch(r'#\[path = "([^"]*)"\]', attr)
        content = sibling(m.group(1), word(k - 1)) if m else None
        if content is None:
            continue
        item_start = tokens[sig[j]].start
        line_start = text.rfind("\n", 0, item_start) + 1
        before = text[line_start:item_start]
        indent = before if before.strip() == "" else re.match(r"[ \t]*", text[line_start:]).group(0)
        head = text[item_start:tokens[sig[k - 1]].end]
        body = "{\n" + reindent(content, len(indent) + INDENT) + indent + "}"
        edits.append((tokens[sig[opener - 1]].start, tokens[sig[k]].end, head + " " + body))
    out = text
    for start, end, replacement in reversed(edits):
        out = out[:start] + replacement + out[end:]
    return out


# ---------------------------------------------------------------- the run
def sibling_loader(moves: Sequence[Move]) -> Callable[[str, str], str | None]:
    by_attr = {m.attr[len('#[path = "'):-2]: m for m in moves}
    return lambda rel, name: by_attr[rel].content if rel in by_attr and by_attr[rel].found.name == name else None


def packaged_sources(root: Path) -> tuple[list[str], set[str], list[str], dict]:
    manifest = tomllib.loads((root / crate.CRATE / "Cargo.toml").read_text(encoding="utf-8"))
    files = crate.git_files(root)
    packaged, _ = crate.package_files(files, crate.Exclude(manifest.get("package", {}).get("exclude", [])))
    roots = crate.crate_roots(manifest, set(files))
    return [f for f in packaged if f.startswith("src/") and f.endswith(".rs")], set(files), roots, manifest


def ledger_fix(root: Path, moved: dict[str, list[Move]]) -> int:
    """Rewrite the `losat_test` paths of the ledger that name a test of a moved module; the number of rows changed."""
    path = root / LEDGER
    if not path.is_file() or not moved:
        return 0
    lines = path.read_text(encoding="utf-8").split("\n")
    changed = 0
    for number, line in enumerate(lines):
        fields = line.split("\t")
        if len(fields) < 4 or "::" not in fields[3]:
            continue
        refs = fields[3].split(";")
        new = []
        for ref in refs:
            rel, sep, name = ref.strip().partition("::")
            parent = rel[len(crate.CRATE) + 1:] if rel.startswith(crate.CRATE + "/") else None
            if sep and parent in moved:
                homes = [m.file for m in moved[parent] if re.search(rf"\bfn\s+{re.escape(name)}\b", m.content)]
                if len(homes) > 1:
                    raise Refusal(f"{LEDGER}:{number + 1}: {ref.strip()} is a test of several moved modules: {homes}")
                if homes:
                    ref = f"{crate.CRATE}/{homes[0]}::{name}"
            new.append(ref)
        if new != refs:
            fields[3] = ";".join(new)
            lines[number] = "\t".join(fields)
            changed += 1
    if changed:
        path.write_text("\n".join(lines), encoding="utf-8")
    return changed


def run(root: Path, fix: bool) -> tuple[list[str], str]:
    """The findings and a summary line; with `fix`, the files are rewritten first."""
    sources, files, roots, _ = packaged_sources(root)
    crate_dir = root / crate.CRATE
    findings: list[str] = []
    plans: dict[str, list[Move]] = {}
    texts: dict[str, str] = {}
    owners_cache: list[set[str]] = []

    def owners() -> set[str]:
        if not owners_cache:
            owners_cache.append(mod_rs_files(crate_dir, files, roots))
        return owners_cache[0]

    occupied = set(files)
    total = 0
    for rel in sources:
        text = (crate_dir / rel).read_text(encoding="utf-8")
        try:
            found = find_modules(text)
        except LexError as error:
            findings.append(f"LOSAT/{rel}: cannot be read ({error})")
            continue
        if not found:
            continue
        total += len(found)
        if not fix:
            findings += [f"LOSAT/{rel}:{f.line}\tinline test-only module `{f.name}`: move it to a sibling file "
                         f"(python3 LOSAT/tests/check_inline_test_modules.py --fix)" for f in found]
            continue
        try:
            moves = plan_file(rel, text, found, owners, occupied)
            new_text = apply_moves(text, moves)
            if restore(new_text, sibling_loader(moves)) != text:
                raise Refusal(f"LOSAT/{rel}: reading the file back from its parent and siblings gives another text")
        except Refusal as why:
            findings.append(str(why))
            continue
        plans[rel] = moves
        texts[rel] = new_text
    if not fix or findings:
        if fix and findings:
            findings.append("nothing was changed")
        return findings, f"inline test modules: {total} in {len(plans) or 0} files"
    for rel, moves in plans.items():
        for m in moves:
            (crate_dir / m.file).write_text(m.content, encoding="utf-8")
        (crate_dir / rel).write_text(texts[rel], encoding="utf-8")
    rows = ledger_fix(root, plans)
    return [], (f"moved {sum(len(m) for m in plans.values())} inline test modules of {len(plans)} files into "
                f"sibling files; {rows} ledger rows rewritten; run cargo fmt")


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--root", type=Path, default=Path("."), help="the repository root (holds LOSAT/)")
    parser.add_argument("--fix", action="store_true", help="move the inline test modules to sibling files")
    args = parser.parse_args(argv)
    try:
        findings, summary = run(args.root.resolve(), args.fix)
    except (OSError, ValueError, Refusal, crate.subprocess.CalledProcessError, tomllib.TOMLDecodeError) as error:
        print(f"inline test-module check error: {error}", file=sys.stderr)
        return 2
    for finding in findings:
        print(finding)
    if findings:
        print(f"inline test-module check: {len(findings)} findings", file=sys.stderr)
        return 1
    print(summary)
    return 0


if __name__ == "__main__":
    sys.exit(main())
