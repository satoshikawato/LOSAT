#!/usr/bin/env python3
"""Keep the references between the programs under `LOSAT/src/algorithm/` from growing in number.

No NCBI counterpart: this checks LOSAT's own module structure, not ported behavior.

The programs are the directories right under `LOSAT/src/algorithm/` (blastn, blastp, blastx,
common, psiblast, pssm_search, tblastn and tblastx at 2026-11-07). In the SWEBOK audit of 2026-11-04
(KA3, S9) each of them reached every other through cycles (common -> blastn -> common, tblastx ->
blastx -> tblastx, ...), so a change to any program could break any other. They are layered,
bottom to top,

    common < blastn < tblastx < blastp < tblastn < blastx

and a program uses only the programs below it; `common` uses none. psiblast and pssm_search are
not layered yet (see LAYERS).

An edge A -> B is a reference from a module under `algorithm::A` to one under `algorithm::B`,
A != B. A reference is every `crate::...` path and `super::` chain, resolved to an absolute module
path through the real module tree (from `lib.rs` through `mod x;`, `#[path = "..."]` and inline
`mod x { ... }`, which add a level for `super::`), with `use` groups expanded (`{a, b::{c, d}}`
gives one site per leaf). Comments, string literals and char literals are blanked first; `$crate`
and the `pub(in ...)` visibility are not references. A site is test code when it is inside an
item, statement or field under a `cfg` that requires `test` (`all(test, x)`; `any(test, x)` is
code), inside a `#[test]` fn, or in a file reached through a `#[cfg(test)] mod x;` (transitively)
or under an inner `#![cfg(test)]`. An edge's scope is `code` when it has a site outside test code,
else `test`.

`LOSAT/tests/program_edges.tsv` lists the reviewed edges (`src<TAB>dst<TAB>scope`) and only
shrinks. The ratchet reports, one tab-separated line each with the first site and the referenced
path:
- UNEXPECTED: an edge that is not listed, or a listed `test` edge that now has code sites;
- STALE: a listed edge that is gone, or a listed `code` edge whose sites are all test code now;
- UPWARD: an edge up the layering, listed or not (it is never also reported as UNEXPECTED).

Usage: python3 LOSAT/tests/check_program_edges.py --root .
       python3 LOSAT/tests/check_program_edges.py --root . --list
       (`--list` prints every edge; `--list | cut -f1-3` regenerates the list's rows)
Exit status: 0 clean, 1 findings, 2 a source file or the list could not be read or is malformed.
"""

from __future__ import annotations

import argparse
import bisect
import os
import posixpath
import re
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import Sequence

# Bottom to top: a program may use only the ones before it, and `common` (the first) uses none.
# psiblast and pssm_search are left out until their cycle with blastp (blastp -> psiblast,
# blastp -> pssm_search and back) is cut; that step inserts them after blastp
# (`blastp < pssm_search < psiblast`). Until then only edges between two programs listed here are
# checked, except that nothing may leave `common`.
LAYERS = ("common", "blastn", "tblastx", "blastp", "tblastn", "blastx")
SRC_DIR = ("LOSAT", "src")
ALGORITHM = "algorithm"
EDGES_HEADER = ("src", "dst", "scope")
SCOPES = ("code", "test")

# Where a comment, string or char literal may start, found left to right (each branch starts with
# a literal, so the regex engine skips the code between them quickly). Prefixes (`b`, `c`, `r#`,
# `br`, `cr`) are read backwards from the quote. A block comment nests; a raw string ends at `"`
# plus as many `#` as it opened with; a char literal needs its closing quote, so a lifetime (`'a`)
# is never one.
CANDIDATE_RE = re.compile(r"/[/*]|\"|'")
CHAR_RE = re.compile(r"'(?:[^'\\\n]|\\(?:x[0-9A-Fa-f]{2}|u\{[0-9A-Fa-f_]+\}|[^\n]))'")
STRING_REST_RE = re.compile(r'[^"\\]*(?:\\.[^"\\]*)*"', re.S)
BLOCK_COMMENT_RE = re.compile(r"/\*|\*/")
BRACKET_RE = re.compile(r"[()\[\]{}]")
LINE_RE = re.compile(r"\n")

ATTR_RE = re.compile(r"#\s*(!?)\s*\[")
PATH_ATTR_RE = re.compile(r'\s*path\s*=\s*"([^"]*)"\s*')
CFG_ATTR_RE = re.compile(r"cfg\((.*)\)")
WORD_RE = re.compile(r"[A-Za-z_][A-Za-z0-9_]*")
MOD_RE = re.compile(r"\bmod\s+(?:r#)?([A-Za-z_][A-Za-z0-9_]*)\s*([;{])")
VISIBILITY_BEFORE_RE = re.compile(r"\bpub\s*(?:\([^()]*\))?\s*$")

# What an attribute is attached to, by its first tokens (visibility and qualifiers skipped).
ITEM_RE = re.compile(
    r"(?:pub\s*(?:\([^()]*\))?\s*)?"
    r"(?:(?:default|const|async|unsafe|safe|extern)\s+)*"
    r"(fn|mod|impl|trait|struct|enum|union|extern|use|static|const|type|let|"
    r"if|match|while|for|loop|unsafe)\b"
)
BLOCK_ITEMS = frozenset({"fn", "mod", "impl", "trait", "struct", "enum", "union", "extern"})
FLOW_ITEMS = frozenset({"if", "match", "while", "for", "loop", "unsafe"})
MACRO_CALL_RE = re.compile(
    r"(?:[A-Za-z_][A-Za-z0-9_]*\s*::\s*)*[A-Za-z_][A-Za-z0-9_]*\s*!\s*"
    r"(?:[A-Za-z_][A-Za-z0-9_]*\s*)?([(\[{])"  # `name!(`, `name! {`, `macro_rules! name {`
)
FIELD_RE = re.compile(r"(?:pub\s*(?:\([^()]*\))?\s*)?[A-Za-z_][A-Za-z0-9_]*\s*:(?!:)")
EXTERN_CRATE_RE = re.compile(r"extern\s+crate\b")
ELSE_RE = re.compile(r"\s*else\b\s*")
SCAN_RE = re.compile(r"[;,()\[\]{}]")
FIELD_SCAN_RE = re.compile(r"->|[;,()\[\]{}<>]")

# A path that starts at the crate root or a parent module. `$crate` (macros) and `a::crate` are not.
REF_RE = re.compile(r"(?<![A-Za-z0-9_$:])(?:crate|super)\s*::")
VISIBILITY_IN_RE = re.compile(r"\(\s*in\s+$")
SEGMENT_RE = re.compile(r"(?:r#)?([A-Za-z_][A-Za-z0-9_]*)")
RENAME_RE = re.compile(r"\s*as\s+(?:r#)?[A-Za-z_][A-Za-z0-9_]*")
SPACE_RE = re.compile(r"\s*")


def blank(text: str) -> str:
    """`text` as spaces, keeping its newlines (so offsets and line numbers stay)."""
    if "\n" not in text:
        return " " * len(text)
    return "\n".join(" " * len(part) for part in text.split("\n"))


def is_word(text: str, i: int) -> bool:
    return i >= 0 and (text[i].isalnum() or text[i] == "_")


def literal_prefix(text: str, quote: int, letters: str) -> int:
    """Start of the literal whose quote is at `quote`: before a `b`/`c` (or `r`, `br`, `cr` when
    `letters` is "r") that is not the end of an identifier."""
    i = quote
    if letters == "r":
        if i > 0 and text[i - 1] == "r":
            i -= 1
            if i > 0 and text[i - 1] in "bc":
                i -= 1
            return quote if is_word(text, i - 1) else i
        return quote
    if i > 0 and text[i - 1] in letters and not is_word(text, i - 2):
        return i - 1
    return quote


def mask(text: str) -> str:
    """`text` with comments (doc and nested block comments too), strings and chars blanked."""
    out = []
    last = pos = 0
    n = len(text)
    while True:
        m = CANDIDATE_RE.search(text, pos)
        if m is None:
            break
        start = m.start()
        c = text[start]
        if c == "/" and text[start + 1] == "/":
            end = text.find("\n", start)
            end = n if end < 0 else end
        elif c == "/":
            depth = 1
            end = n
            for token in BLOCK_COMMENT_RE.finditer(text, start + 2):
                depth += 1 if token.group() == "/*" else -1
                if depth == 0:
                    end = token.end()
                    break
        elif c == '"':
            hashes = start
            while hashes > last and text[hashes - 1] == "#":
                hashes -= 1
            raw = literal_prefix(text, hashes, "r")
            if raw < hashes:  # r"..", r#".."#, br"..", cr".."
                close = '"' + text[hashes:start]
                found = text.find(close, start + 1)
                end = n if found < 0 else found + len(close)
                start = raw
            else:
                rest = STRING_REST_RE.match(text, start + 1)
                end = n if rest is None else rest.end()
                start = literal_prefix(text, start, "bc")
        else:
            char = CHAR_RE.match(text, start)
            if char is None:  # a lifetime or a label
                pos = start + 1
                continue
            end = char.end()
            start = literal_prefix(text, start, "b")
        out.append(text[last:start])
        out.append(blank(text[start:end]))
        last = pos = end
    out.append(text[last:])
    return "".join(out)


def match_brackets(masked: str) -> dict[int, int]:
    """Offset of each `(`, `[` and `{` -> offset of its closer."""
    pairs: dict[int, int] = {}
    stack: list[int] = []
    for m in BRACKET_RE.finditer(masked):
        if m.group() in "([{":
            stack.append(m.start())
        elif stack:
            pairs[stack.pop()] = m.start()
    return pairs


def cfg_requires_test(predicate: str) -> bool:
    """Whether the cfg predicate (the inside of `cfg(...)`, no spaces, string values blanked) can
    hold only when `test` does: `test`, `all(test, x)`, `any(test, all(test, y))`; not
    `any(test, x)`."""

    def term(i: int) -> tuple[bool, int]:
        m = WORD_RE.match(predicate, i)
        if m is None:
            while i < len(predicate) and predicate[i] not in ",)":
                i += 1
            return False, i
        name, i = m.group(), m.end()
        if predicate.startswith("(", i):
            kids = []
            i += 1
            while i < len(predicate) and predicate[i] != ")":
                value, i = term(i)
                kids.append(value)
                if predicate.startswith(",", i):
                    i += 1
            i += 1
            if name == "all":
                return any(kids), i
            if name == "any":
                return bool(kids) and all(kids), i
            return False, i
        if predicate.startswith("=", i):
            while i < len(predicate) and predicate[i] not in ",)":
                i += 1
            return False, i
        return name == "test", i

    return term(0)[0]


@dataclass
class Attr:
    start: int
    end: int  # after the `]`
    inner: bool
    body: str  # masked, without spaces
    raw: str  # the source between the brackets


@dataclass
class ModDecl:
    name: str
    at: int  # offset of `mod`
    inline: bool
    path: str | None  # `#[path = "..."]`
    open: int = -1  # inline: the braces
    close: int = -1


@dataclass
class SourceFile:
    rel: str  # under LOSAT/src, with `/`
    text: str
    masked: str = ""
    pairs: dict[int, int] = field(default_factory=dict)
    line_starts: list[int] = field(default_factory=list)
    mods: list[ModDecl] = field(default_factory=list)
    inline_mods: list[ModDecl] = field(default_factory=list)
    test_spans: list[tuple[int, int]] = field(default_factory=list)
    whole_file_test: bool = False

    def __post_init__(self) -> None:
        self.masked = mask(self.text)
        self.pairs = match_brackets(self.masked)
        self.line_starts = [0] + [m.end() for m in LINE_RE.finditer(self.text)]
        attrs = self._attributes()
        self._find_mods(attrs)
        self._find_test_spans(attrs)

    def line(self, offset: int) -> int:
        return bisect.bisect_right(self.line_starts, offset)

    def _attributes(self) -> list[Attr]:
        attrs = []
        for m in ATTR_RE.finditer(self.masked):
            open_ = m.end() - 1
            close = self.pairs.get(open_)
            if close is None:
                continue
            body = re.sub(r"\s+", "", self.masked[open_ + 1 : close])
            raw = self.text[open_ + 1 : close]
            attrs.append(Attr(m.start(), close + 1, bool(m.group(1)), body, raw))
        return attrs

    def _find_mods(self, attrs: list[Attr]) -> None:
        by_end = {attr.end: attr for attr in attrs if not attr.inner}
        for m in MOD_RE.finditer(self.masked):
            at = m.start()
            vis = VISIBILITY_BEFORE_RE.search(self.masked, max(0, at - 160), at)
            q = vis.start() if vis else at
            path = None
            while True:  # the outer attributes right before the item, last to first
                j = q
                while j > 0 and self.masked[j - 1].isspace():
                    j -= 1
                attr = by_end.get(j)
                if attr is None:
                    break
                found = PATH_ATTR_RE.fullmatch(attr.raw)
                if found and path is None:
                    path = found.group(1)
                q = attr.start
            decl = ModDecl(m.group(1), at, m.group(2) == "{", path)
            if decl.inline:
                decl.open = m.end() - 1
                decl.close = self.pairs.get(decl.open, len(self.masked))
                self.inline_mods.append(decl)
            else:
                self.mods.append(decl)

    def _find_test_spans(self, attrs: list[Attr]) -> None:
        index = 0
        while index < len(attrs):
            attr = attrs[index]
            index += 1
            cfg = CFG_ATTR_RE.fullmatch(attr.body)
            is_test = attr.body == "test" and not attr.inner
            if not is_test and not (cfg and cfg_requires_test(cfg.group(1))):
                continue
            if attr.inner:
                enclosing = [(o, c) for o, c in self.pairs.items()
                             if o < attr.start < c and self.masked[o] == "{"]
                if enclosing:
                    self.test_spans.append(max(enclosing))
                else:
                    self.whole_file_test = True
                continue
            item = attr.end
            while True:  # skip the attributes after this one
                item = SPACE_RE.match(self.masked, item).end()
                following = ATTR_RE.match(self.masked, item)
                if following is None or following.group(1):
                    break
                close = self.pairs.get(following.end() - 1)
                if close is None:
                    break
                item = close + 1
            self.test_spans.append((attr.start, self._extent(item)))

    def _scan(self, i: int, stop_at_block: bool, stop_at_comma: bool) -> int:
        """End of a run of tokens from `i`: after the first `;` (or `,`, or `{...}` when asked)
        outside brackets, or at the closer of the enclosing group."""
        masked, pairs = self.masked, self.pairs
        while True:
            m = SCAN_RE.search(masked, i)
            if m is None:
                return len(masked)
            c, p = m.group(), m.start()
            if c in "([{":
                close = pairs.get(p)
                if close is None:
                    return len(masked)
                if c == "{" and stop_at_block:
                    return close + 1
                i = close + 1
            elif c == ";" or (c == "," and stop_at_comma):
                return p + 1
            elif c == ",":
                i = p + 1
            else:
                return p

    def _field_end(self, i: int) -> int:
        """A field `name: Type<A, B>,` ends at a `,` outside `<...>` too."""
        masked, pairs = self.masked, self.pairs
        angle = 0
        while True:
            m = FIELD_SCAN_RE.search(masked, i)
            if m is None:
                return len(masked)
            c, p = m.group(), m.start()
            if c in "([{":
                close = pairs.get(p)
                if close is None:
                    return len(masked)
                i = close + 1
                continue
            i = m.end()
            if c == "<":
                angle += 1
            elif c == ">":
                angle = max(0, angle - 1)
            elif c == ";" or (c == "," and angle == 0):
                return p + 1
            elif c in ")]}":
                return p

    def _extent(self, i: int) -> int:
        """End of the item, statement, field, variant or arm that starts at `i`."""
        masked, pairs = self.masked, self.pairs
        m = ITEM_RE.match(masked, i)
        if m:
            keyword = m.group(1)
            extern_crate = keyword == "extern" and EXTERN_CRATE_RE.match(masked, m.start(1))
            if keyword in BLOCK_ITEMS and not extern_crate:
                return self._scan(m.end(), stop_at_block=True, stop_at_comma=False)
            if keyword in FLOW_ITEMS:
                end = self._scan(m.end(), stop_at_block=True, stop_at_comma=False)
                while keyword == "if":
                    other = ELSE_RE.match(masked, end)
                    if other is None:
                        break
                    end = self._scan(other.end(), stop_at_block=True, stop_at_comma=False)
                return end
            return self._scan(m.end(), stop_at_block=False, stop_at_comma=False)
        m = MACRO_CALL_RE.match(masked, i)
        if m:
            close = pairs.get(m.start(1))
            if close is None:
                return len(masked)
            end = close + 1
            if m.group(1) != "{":
                after = SPACE_RE.match(masked, end).end()
                if masked.startswith(";", after):
                    end = after + 1
            return end
        if masked.startswith("{", i):
            return pairs.get(i, len(masked) - 1) + 1
        if FIELD_RE.match(masked, i):
            return self._field_end(i)
        return self._scan(i, stop_at_block=False, stop_at_comma=True)

    def in_test(self, offset: int) -> bool:
        return self.whole_file_test or any(a <= offset < b for a, b in self.test_spans)

    def inline_chain(self, offset: int) -> list[ModDecl]:
        """The inline modules around `offset`, outermost first."""
        return [decl for decl in self.inline_mods if decl.open < offset < decl.close]


@dataclass
class Context:
    """Where a file sits in the module tree, and how rustc resolves its `mod x;` (rustc_expand
    module.rs: a `foo.rs` looks for `foo/x.rs`; a `mod.rs`, the crate root and a `#[path]` file for
    `x.rs`)."""

    modpath: tuple[str, ...]
    dir: str
    relative: str | None
    test: bool


@dataclass
class Site:
    rel: str  # under the repository root
    line: int
    target: tuple[str, ...]
    test: bool

    def where(self) -> str:
        return f"{self.rel}:{self.line}"

    def reference(self) -> str:
        return "crate::" + "::".join(self.target)


@dataclass
class Edge:
    src: str
    dst: str
    code: list[Site] = field(default_factory=list)
    test: list[Site] = field(default_factory=list)

    @property
    def scope(self) -> str:
        return "code" if self.code else "test"

    def first(self) -> Site:
        return (self.code or self.test)[0]


def join(directory: str, path: str) -> str:
    return posixpath.normpath(posixpath.join(directory, path) if directory else path)


def child_file(
    ctx: Context, chain: list[ModDecl], decl: ModDecl, existing: set[str]
) -> tuple[str, Context] | None:
    """The file and context of `mod <decl>;` in the inline modules `chain` of a file at `ctx`."""
    directory, relative = ctx.dir, ctx.relative
    for inline in chain:
        if inline.path is not None:
            directory = join(directory, inline.path)
        else:
            if relative:
                directory = join(directory, relative)
            directory = join(directory, inline.name)
        relative = None
    modpath = ctx.modpath + tuple(inline.name for inline in chain) + (decl.name,)
    if decl.path is not None:
        rel = join(directory, decl.path)
        if rel not in existing:
            return None
        return rel, Context(modpath, posixpath.dirname(rel), None, ctx.test)
    base = join(directory, relative) if relative else directory
    flat = join(base, decl.name + ".rs")
    if flat in existing:
        return flat, Context(modpath, base, decl.name, ctx.test)
    nested = join(base, decl.name + "/mod.rs")
    if nested in existing:
        return nested, Context(modpath, posixpath.dirname(nested), None, ctx.test)
    return None


def parse_use_tree(masked: str, i: int) -> tuple[list[list[str]], int]:
    """The paths named by the path or `use` tree at `i` (one per leaf of a group), and its end."""
    i = SPACE_RE.match(masked, i).end()
    if masked.startswith("{", i):
        return parse_group(masked, i)
    segments: list[str] = []
    while True:
        i = SPACE_RE.match(masked, i).end()
        m = SEGMENT_RE.match(masked, i)
        if m:
            segments.append(m.group(1))
            i = m.end()
        elif masked.startswith("*", i):
            segments.append("*")
            return [segments], i + 1
        else:
            break
        j = SPACE_RE.match(masked, i).end()
        if not masked.startswith("::", j):
            break
        k = SPACE_RE.match(masked, j + 2).end()
        if masked.startswith("{", k):
            leaves, end = parse_group(masked, k)
            return [segments + leaf for leaf in leaves], end
        i = j + 2
    return [segments], i


def parse_group(masked: str, i: int) -> tuple[list[list[str]], int]:
    leaves: list[list[str]] = []
    i += 1
    while i < len(masked):
        i = SPACE_RE.match(masked, i).end()
        if masked.startswith("}", i):
            return leaves, i + 1
        if masked.startswith(",", i):
            i += 1
            continue
        sub, end = parse_use_tree(masked, i)
        if end <= i:
            i += 1
            continue
        leaves.extend(leaf for leaf in sub if leaf)
        i = end
        rename = RENAME_RE.match(masked, i)
        if rename:
            i = rename.end()
    return leaves, i


def resolve(leaf: list[str], here: tuple[str, ...]) -> tuple[str, ...] | None:
    """The absolute module path of `leaf` (`crate::...` or `super::...`) written in `here`."""
    if leaf[0] == "crate":
        base, rest = (), leaf[1:]
    else:
        up = 0
        while up < len(leaf) and leaf[up] == "super":
            up += 1
        if up > len(here):
            return None
        base, rest = here[: len(here) - up], leaf[up:]
    return base + tuple(segment for segment in rest if segment not in ("self", "*"))


def program_of(modpath: tuple[str, ...], programs: frozenset[str]) -> str | None:
    if len(modpath) >= 2 and modpath[0] == ALGORITHM and modpath[1] in programs:
        return modpath[1]
    return None


def find_edges(root: Path) -> dict[tuple[str, str], Edge]:
    """Every edge between two programs, with its sites in file and line order."""
    src_root = root.joinpath(*SRC_DIR)
    algorithm_dir = src_root / ALGORITHM
    if not (src_root / "lib.rs").is_file() or not algorithm_dir.is_dir():
        raise ValueError(f"{src_root}: no lib.rs or {ALGORITHM}/ (is --root the repository root?)")
    programs = frozenset(
        entry.name for entry in algorithm_dir.iterdir()
        if entry.is_dir() and any(entry.rglob("*.rs"))
    )
    existing = set()
    for directory, subdirs, files in os.walk(src_root):
        subdirs.sort()
        for name in files:
            if name.endswith(".rs"):
                existing.add(Path(directory, name).relative_to(src_root).as_posix())

    sources: dict[str, SourceFile] = {}
    contexts: dict[str, Context] = {}

    def load(rel: str) -> SourceFile:
        if rel not in sources:
            text = src_root.joinpath(rel).read_text(encoding="utf-8")
            sources[rel] = SourceFile(rel, text)
        return sources[rel]

    # Walk the module tree from lib.rs into `algorithm` (only modules under it can be a source). A
    # file reached twice keeps its first module path; it is test code only if every way to it is.
    stack = ["lib.rs"]
    contexts["lib.rs"] = Context((), "", None, False)
    while stack:
        rel = stack.pop()
        ctx = contexts[rel]
        source = load(rel)
        for decl in source.mods:
            chain = source.inline_chain(decl.at)
            child = child_file(ctx, chain, decl, existing)
            if child is None:
                continue
            child_rel, child_ctx = child
            if child_ctx.modpath[:1] != (ALGORITHM,):
                continue
            child_ctx.test = ctx.test or source.whole_file_test or source.in_test(decl.at)
            seen = contexts.get(child_rel)
            if seen is None:
                contexts[child_rel] = child_ctx
                stack.append(child_rel)
            elif seen.test and not child_ctx.test:
                seen.test = False
                stack.append(child_rel)

    edges: dict[tuple[str, str], Edge] = {}
    prefix = "/".join(SRC_DIR) + "/"
    for rel in sorted(contexts):
        ctx = contexts[rel]
        source = sources[rel]
        if len(ctx.modpath) >= 2 and program_of(ctx.modpath, programs) is None:
            continue
        masked = source.masked
        for m in REF_RE.finditer(masked):
            at = m.start()
            if VISIBILITY_IN_RE.search(masked, max(0, at - 16), at):
                continue  # `pub(in crate::a::b)` names an ancestor, never another program
            here = ctx.modpath + tuple(decl.name for decl in source.inline_chain(at))
            src = program_of(here, programs)
            if src is None:
                continue
            leaves, _ = parse_use_tree(masked, at)
            test = None
            for leaf in leaves:
                target = resolve(leaf, here) if leaf else None
                dst = program_of(target, programs) if target else None
                if dst is None or dst == src:
                    continue
                if test is None:
                    test = ctx.test or source.in_test(at)
                edge = edges.setdefault((src, dst), Edge(src, dst))
                site = Site(prefix + rel, source.line(at), target, test)
                (edge.test if test else edge.code).append(site)
    return dict(sorted(edges.items()))


def load_listed(path: Path) -> dict[tuple[str, str], str]:
    """The rows of the edge list: (src, dst) -> scope."""
    rows: list[tuple[int, list[str]]] = []
    for number, line in enumerate(path.read_text(encoding="utf-8").splitlines(), start=1):
        if not line.strip() or line.startswith("#"):
            continue
        rows.append((number, line.split("\t")))
    if not rows or tuple(rows[0][1]) != EDGES_HEADER:
        raise ValueError(f"{path}: the first row must be the header `{'<TAB>'.join(EDGES_HEADER)}`")
    listed: dict[tuple[str, str], str] = {}
    for number, cells in rows[1:]:
        if len(cells) != 3 or not all(cells) or cells[2] not in SCOPES or cells[0] == cells[1]:
            raise ValueError(
                f"{path}:{number}: want `src<TAB>dst<TAB>scope` with src != dst and scope "
                f"{' or '.join(SCOPES)}"
            )
        key = (cells[0], cells[1])
        if key in listed:
            raise ValueError(f"{path}:{number}: {cells[0]} -> {cells[1]} is listed twice")
        listed[key] = cells[2]
    return listed


def is_upward(src: str, dst: str) -> bool:
    if src == LAYERS[0]:
        return True
    return src in LAYERS and dst in LAYERS and LAYERS.index(src) < LAYERS.index(dst)


def evaluate(
    edges: dict[tuple[str, str], Edge], listed: dict[tuple[str, str], str]
) -> tuple[list[Edge], list[tuple[tuple[str, str], str, Edge | None]], list[Edge]]:
    """(unexpected, stale, upward); a stale entry is (listed edge, its scope, edge now or None)."""
    unexpected, stale, upward = [], [], []
    for key, edge in edges.items():
        if is_upward(*key):
            upward.append(edge)
        elif key not in listed or (listed[key] == "test" and edge.scope == "code"):
            unexpected.append(edge)
    for key, scope in listed.items():
        edge = edges.get(key)
        if is_upward(*key) and edge is not None:
            continue  # reported as UPWARD
        if edge is None or (scope == "code" and edge.scope == "test"):
            stale.append((key, scope, edge))
    return unexpected, stale, upward


def print_list(edges: dict[tuple[str, str], Edge]) -> None:
    print("\t".join(EDGES_HEADER + ("code_sites", "test_sites", "first_site", "reference")))
    for edge in edges.values():
        first = edge.first()
        print(f"{edge.src}\t{edge.dst}\t{edge.scope}\t{len(edge.code)}\t{len(edge.test)}\t"
              f"{first.where()}\t{first.reference()}")


def parse_args(argv: Sequence[str] | None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--root", type=Path, default=Path("."), help="repository root")
    parser.add_argument("--edges", type=Path, help="default: <root>/LOSAT/tests/program_edges.tsv")
    parser.add_argument("--list", action="store_true",
                        help="print every edge with its code and test sites and first site")
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    root = args.root.resolve()
    edges_path = args.edges or root / "LOSAT" / "tests" / "program_edges.tsv"
    try:
        edges = find_edges(root)
        if args.list:
            print_list(edges)
            return 0
        listed = load_listed(edges_path)
    except (OSError, ValueError) as error:
        print(f"program-edge ratchet error: {error}", file=sys.stderr)
        return 2
    unexpected, stale, upward = evaluate(edges, listed)
    layering = " < ".join(LAYERS)
    bottom = f"{ALGORITHM}/{LAYERS[0]}"
    for edge in upward:
        first = edge.first()
        if edge.src == LAYERS[0]:
            fix = f"{bottom} uses no program: move the item into {bottom}"
        else:
            fix = (f"goes up the layering {layering}: move the item into {ALGORITHM}/{edge.src}, "
                   f"a program below it or {bottom}")
        print(f"UPWARD\t{edge.src} -> {edge.dst}\t{edge.scope}\t{first.where()}\t"
              f"{first.reference()}\t{fix}; a row in {edges_path.name} does not allow it")
    for edge in unexpected:
        first = edge.first()
        was = " (listed as test only)" if (edge.src, edge.dst) in listed else ""
        print(f"UNEXPECTED\t{edge.src} -> {edge.dst}\t{edge.scope}{was}\t{first.where()}\t"
              f"{first.reference()}\tnew reference between programs: take the item from a program "
              f"below {edge.src} or from {bottom}, or move it there; a new downward edge needs a "
              f"reviewed row in {edges_path.name}")
    for (src, dst), scope, edge in stale:
        if edge is None:
            print(f"STALE\t{src} -> {dst}\t{scope}\t-\t-\tlisted but gone: remove the row from "
                  f"{edges_path.name} (fewer edges are welcome)")
        else:
            first = edge.first()
            print(f"STALE\t{src} -> {dst}\t{scope}\t{first.where()}\t{first.reference()}\t"
                  f"listed as code, now only in test code: change the row's scope to test")
    print(f"program-edge ratchet\tedges={len(edges)}\tlisted={len(listed)}\t"
          f"unexpected={len(unexpected)}\tstale={len(stale)}\tupward={len(upward)}")
    return 0 if not unexpected and not stale and not upward else 1


if __name__ == "__main__":
    raise SystemExit(main())
