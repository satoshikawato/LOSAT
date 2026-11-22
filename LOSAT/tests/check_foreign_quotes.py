#!/usr/bin/env python3
"""Fail when a comment holds source text, or when a citation of NCBI's source is not pinned.

No NCBI counterpart: this checks LOSAT's own comments, not ported behavior.

Comments cite source code by path and lines and do not copy it, NCBI's included (AGENTS.md
rule 4; the maintainer's decisions of 2026-11-20): `NCBI reference (598d8ae6): c++/<path>:<a>-<b>`
for NCBI BLAST+ 2.17.0 (`NCBI reference (RpsbProc 0.5.1, 2025-05-13): RpsbProc/<path>:<lines>`
for rpsbproc), `<project> <version>: <path>:<lines>` for another project, also when NCBI's tree
carries a copy of it (`LMDB 0.9.24 (NCBI 598d8ae6 copy): c++/src/util/lmdb/mdb.c:...`).
`scripts/ncbi_ref.py` prints the cited lines from a checkout of the pinned tree.

Checks, over the comments of the repository's source files (`COMMENT_STYLES`: line comments of
Rust, C, C++, JavaScript and TypeScript; `#` comments of TOML, YAML, shell and Python, and
Python's docstrings; every line of Markdown and text files):
- QUOTE: source text in a comment.
  - A fenced block (```` ``` ````) whose language is one of `SOURCE_LANGUAGES` (C, C++, ASN.1,
    CMake), whatever it cites.
  - In a source file, a fenced block right under a citation: the comment lines above it, back
    to an empty comment line, cite a source (a labelled citation, `NCBI reference ...`, a path
    with a source extension or of NCBI's tree), unless it is a Rust example (`RUST_FENCES`).
  - A fenced block whose citation is another project's: the last of the comment lines right
    above it that starts with a citation label, else each of them that names a source; without
    one, in a source file, the last citation above it in the file (a quote inside a function
    cited above the function). Another project: it names a source of `FOREIGN_SOURCES` or a
    path of `NON_NCBI_PATHS`, or its label is not one of `NCBI_LABEL_RE`.
  - A comment line that is a line of C or C++ (`CODE_LINE_RE`) under a labelled citation, in
    the same comment: a quote without a fence.
  - In a source file, a C or C++ statement in backticks inside a comment's prose
    (`inline_statement`), cited or not.
  - In a Rust, Python, TOML, YAML or shell file, a C or C++ expression in backticks inside a
    comment's prose (`inline_expression`): a span holds names, values, paths, command lines,
    messages and LOSAT's own code, not an assignment, an operation or a call with arguments of
    C or C++. One with a mark of C or C++ is a quote; in Rust, so is one that uses a name
    LOSAT's code nowhere uses.
  - In a Rust, Python, TOML, YAML or shell file, an unbackticked C or C++ statement or
    expression outside backticks (`unbackticked_c_code`) in a line comment, docstring or block
    comment, cited or not (valid 1-line citations, Rust pseudo-code, plain math equations,
    output labels and observed runtime messages pass).
- LABEL: a comment line labelled as NCBI's (`NCBI reference ...`, `NCBI <path>`) whose label or
  path names another project's source (the code is not NCBI's).
- PIN (source files): a citation of NCBI's source that does not name the pinned tree: an
  `NCBI reference` label other than `PINNED_LABELS`, or a path through a local checkout
  (`CHECKOUT_PATH_RE`: a directory named for NCBI's tree, or `GitHub`, before `c++/`).
- REF (source files): a path of the pinned trees (`c++/...`, `RpsbProc/...`) that is not a file
  or directory of them, or lines (`:<a>-<b>`) outside the file. The trees' files and their line
  counts are in `ncbi_source_lines.tsv` (`--write-table` regenerates it from the checkouts).
A fenced block with no citation above it is not a quote unless its language is a source
language (Rust examples, commands and output stay). In Markdown and text files only the
paragraph right above a fence is read for its citation (their fences are mostly commands and
examples), and PIN and REF do not apply (dated records keep the paths they were written with).

`--fix` removes what QUOTE finds and pins what PIN finds (a branch that merges a dev with
quotes removed runs it; it changes comment lines only), then checks again; REF and statements
and expressions in backticks are fixed by hand (name the function or member and say in words
what it does).

The tables below are the one place that says which sources are NCBI's. There is no
allowlist of exceptions: cite the file and lines instead of quoting them.

Usage: python3 LOSAT/tests/check_foreign_quotes.py --root . [--fix]
       python3 LOSAT/tests/check_foreign_quotes.py --write-table --ncbi-src $NCBI_SRC \\
           --rpsbproc-src <dir holding RpsbProc/>
"""

from __future__ import annotations

import argparse
import ast
import io
import os
import re
import sys
import tokenize
from dataclasses import dataclass
from pathlib import Path
from typing import Callable, Iterable, Sequence

# Labels of NCBI's own source: `NCBI reference (598d8ae6): c++/...`, `NCBI unit test ...`,
# `NCBI 2.17.0: ...`, `Reference: c++/...` (PIN then asks for the pinned label).
NCBI_LABEL_RE = re.compile(r"^(?:NCBI(?!['\u2019]s)\b|Reference\b)")
# A line that says it cites NCBI's source.
NCBI_CITATION_RE = re.compile(r"^(?:NCBI reference|NCBI unit test|Reference)\b")

# The languages of a fenced block that is source text, whatever it cites.
SOURCE_LANGUAGES = frozenset({"c", "c++", "cpp", "cc", "cxx", "h", "hpp", "asn", "asn1", "cmake"})
# The fences of Rust examples (rustdoc's attributes), which are LOSAT's code, not a quote.
RUST_FENCES = frozenset({"rust", "no_run", "ignore", "should_panic", "compile_fail", "edition2021"})
# The labels of a citation of NCBI's source that name the pinned trees.
PINNED_LABELS = (
    "NCBI reference (598d8ae6)",
    "NCBI references (598d8ae6)",
    "NCBI reference (RpsbProc 0.5.1, 2025-05-13)",
    "NCBI unit test (598d8ae6)",
)
# `NCBI reference` used as a label: followed by `:` or by a parenthesis.
NCBI_LABEL_USE_RE = re.compile(r"\bNCBI (?:reference|unit test)s?(?=\s*[:(])(?:\s*\([^)]*\))?")
# A path of a local checkout of NCBI's tree.
CHECKOUT_PATH_RE = re.compile(r"(?:ncbi-blast[\w.+-]*|GitHub)/c\+\+/")
# A path of the pinned trees, with its lines: `c++/<path>:12-20,31`.
TREE_PATH_RE = re.compile(
    r"(?<![\w.+/-])(?P<path>(?:c\+\+|RpsbProc)/[\w.+/-]*[\w+/])"
    r"(?::(?P<ranges>\d+(?:-\d+)?(?:,\s?\d+(?:-\d+)?)*))?"
)
TABLE_PATH = Path(__file__).with_name("ncbi_source_lines.tsv")
# What follows the last cited path and lines of a line (`...c:410-420: opt.SetX(true);`).
CITATION_TAIL_RE = re.compile(r"\s*[:,;]?\s*(?P<tail>[^`].*)$")
CODE_AFTER_CITATION = "code after the citation"

# Files in NCBI's tree that are other projects' code: they carry no NCBI PUBLIC DOMAIN NOTICE.
# A path is matched when it contains the entry (after `c++/`; `boost_erf.c` by its name alone).
NON_NCBI_PATHS: tuple[tuple[str, str], ...] = (
    ("LMDB", "src/util/lmdb/"),
    ("LMDB", "include/util/lmdb/"),
    ("lmdb++", "include/util/lmdbxx/"),
    ("Boost", "boost_erf.c"),
    ("BitMagic", "include/util/bitset/"),
    ("bzip2", "src/util/compress/bzip2/"),
    ("bzip2", "include/util/compress/bzip2/"),
    ("zlib", "src/util/compress/zlib"),
    ("zlib", "include/util/compress/zlib/"),
    ("zlib", "include/util/compress/zlib_cloudflare/"),
    ("PCRE", "src/util/regexp/"),
    ("PCRE", "include/util/regexp/"),
    ("CityHash", "src/util/checksum/cityhash/"),
    ("FarmHash", "src/util/checksum/farmhash/"),
    ("MurmurHash", "src/util/checksum/murmurhash/"),
    ("Mbed TLS", "src/connect/mbedtls/"),
    ("RapidJSON", "include/misc/jsonwrapp/rapidjson"),
)

# Other projects named in a citation line, by their name or by a path that only they have.
FOREIGN_SOURCES: tuple[tuple[str, re.Pattern[str]], ...] = (
    (
        "libstdc++",
        re.compile(
            r"libstdc\+\+|\bGCC\b|/usr/include/|\binclude/c\+\+/|\bc\+\+/\d+/"
            r"|(?:^|[\s(`])(?:bits|ext)/\w+\.(?:h|tcc)\b|\bstl_(?!msvc)\w+\.h\b"
            r"|\bstd/(?:fstream|[io]?stream|streambuf)\b"
        ),
    ),
    ("glibc", re.compile(r"\bglibc\b|\blibio/|\bsysdeps/|\bGNU C Library\b")),
    (
        "Microsoft C runtime",
        re.compile(r"\bUCRT\b|\bucrt/|CRTSource|Windows SDK|\bMSVC\b|\bmsvcrt\b"),
    ),
    ("musl", re.compile(r"\bmusl\b")),
    ("LLVM", re.compile(r"\bLLVM\b|\blibc\+\+(?:abi)?\b")),
    ("Wine/ReactOS/MinGW", re.compile(r"\bWine\b|\bReactOS\b|\bmingw", re.I)),
    ("SQLite", re.compile(r"\bSQLite\b[^:]*:\s*\S*\.c\b|\bsrc/btree\.c\b")),
    ("bzip2 crate", re.compile(r"\bbzip2 \d")),
    ("ruzstd", re.compile(r"\bruzstd\b")),
    ("Arm Optimized Routines", re.compile(r"\bOptimized Routines\b|\bAOR\b")),
    ("fdlibm", re.compile(r"\bfdlibm\b|\be_\w+\.c\b")),
)

COMMENT_RE = re.compile(r"^\s*//[/!]?(?P<text>.*)$")
HASH_COMMENT_RE = re.compile(r"^\s*#(?![!\[])(?P<text>.*)$")
# A Python docstring's opening line: the quotes start the line.
DOCSTRING_RE = re.compile(r"^\s*[rRuU]?(?P<quotes>\"\"\"|\'\'\')(?P<text>.*)$")
# The comment syntax of each kind of source file, by its extension.
COMMENT_STYLES: dict[str, str] = {
    ".rs": "slash",
    ".c": "slash",
    ".h": "slash",
    ".cc": "slash",
    ".cpp": "slash",
    ".hpp": "slash",
    ".js": "slash",
    ".mjs": "slash",
    ".cjs": "slash",
    ".ts": "slash",
    ".tsx": "slash",
    ".toml": "hash",
    ".yml": "hash",
    ".yaml": "hash",
    ".sh": "hash",
    ".py": "python",
    ".md": "text",
    ".txt": "text",
}
FENCE_RE = re.compile(r"^\s*```")
# A cited source: a path with a source extension, a libstdc++ header path, or `<path>:<line>`.
# (Each path is matched from the start of its token only: a long token, such as a line of a
# sequence in a text file, is then read once, not once from each of its letters.)
CITE_RE = re.compile(
    r"(?<![\w.+/-])[\w.+-]*(?:/[\w.+-]+)*\.(?:c|cc|cpp|cxx|h|hh|hpp|hxx|inl|tcc|asn|def|prt)\b"
    r"|\b(?:bits|ext|std)/\w+"
    r"|(?<![\w.+/-])[\w.+-]+(?:/[\w.+-]+)+:\d"
)
# A labelled citation: the line starts with `<label>: <path>:<line>` or `<label>: <path> (`.
LABELLED_RE = re.compile(
    r"^\s*(?P<label>[A-Za-z][^:`]{0,100}?)\s*:\s+(?P<path>[\w.+/-]+)(?::\d|\s*\(|\s*$|,)"
)
# A line that starts with `<project> [<version>] reference`.
LABEL_START_RE = re.compile(r"^[A-Za-z][\w.+ -]{0,40}\breference\b")
# A line labelled as NCBI's by a path: `NCBI c++/src/...:<lines>`, `NCBI 2.17.0 boost_erf.c:...`,
# the path in backticks or not.
NCBI_PATH_LABEL_RE = re.compile(
    r"^NCBI(?:\s+[^\s:`]+){0,4}?\s+`?(?P<path>[\w.+/-]+\.(?:c|cc|cpp|cxx|h|hh|hpp|hxx|inl|tcc))\b"
)
# A comment line that is a line of C or C++: a statement, a brace, a definition, a `#define`.
CODE_LINE_RE = re.compile(
    r"(?:=|\breturn\b|\w\()[^.`]*;\s*$|^\s*\}\s*$|\)\s*\{\s*$|^\s*#\s*define\b"
    r"|^\s*(?:static\s+)?(?:const\s+)?(?:unsigned\s+|signed\s+)?"
    r"(?:int|double|float|char|void|size_t|long|bool|Int4|Int8|Uint4|Uint8)\s+\**\w+\s*[(=;]"
    # C statements: an `if` with a `continue`, an increment, `break`, `} else {`, `case X:`, a call
    # with the semicolon after each
    r"|^\s*(?:\}\s*)?(?:if|for|while|switch)\s*\(.*\)\s*(?:\{|[^`]*;)?\s*$"
    r"|^\s*(?:break|continue)\s*;\s*$|\w(?:\+\+|--)\s*;\s*$|^\s*(?:\+\+|--)\w+\s*;\s*$"
    r"|^\s*\}\s*else\b|^\s*else\s*\{\s*$|^\s*case\s+[\w:]+\s*:\s*|^\s*default\s*:\s*$"
    r"|^\s*[\w*&>.\[\]-]*\w\([^`]*\)\s*;\s*$|^\s*[\w*&>.\[\]-]*\w\([^()`]*,\s*$"
    r"|^\s*[\w.>-]*\w\([^`]*\)\s*,\s*$"
)
# Four lowercase words in a row: prose, which ends a quoted call's continuation lines.
PROSE_RE = re.compile(r"\b[a-z]+ [a-z]+ [a-z]+ [a-z]+\b")
# The indentation of a comment line's text (after its marker; a docstring's line as it is):
# an indented block under a source label is quoted code even after a blank line.
MARKER_RE = re.compile(r"^\s*(?://[/!]?|#!?)")
SKIPPED_DIRS = frozenset(
    {".git", ".worktrees", "target", "node_modules", ".venv", "venv", "__pycache__"}
)


@dataclass(frozen=True)
class Finding:
    path: str
    line: int
    kind: str
    source: str
    text: str


def source_files(root: Path) -> Iterable[Path]:
    """Every source file under `root` whose comments are read (`COMMENT_STYLES`), skipping build
    output, worktrees, environments and package caches."""
    for directory, subdirs, files in os.walk(root):
        subdirs[:] = sorted(d for d in subdirs if d not in SKIPPED_DIRS)
        for name in sorted(files):
            if Path(name).suffix in COMMENT_STYLES:
                yield Path(directory) / name


def python_texts(lines: Sequence[str]) -> list[str | None] | None:
    """The comment text of each line of a Python file: full-line `#` comments (tokenize) and the
    lines of string statements, such as docstrings (ast); None when the file does not parse."""
    source = "\n".join(lines) + "\n"
    try:
        tree = ast.parse(source)
        tokens = list(tokenize.generate_tokens(io.StringIO(source).readline))
    except (SyntaxError, tokenize.TokenError, IndentationError, ValueError):
        return None
    texts: list[str | None] = [None] * len(lines)
    for node in ast.walk(tree):
        if (
            isinstance(node, ast.Expr)
            and isinstance(node.value, ast.Constant)
            and isinstance(node.value.value, str)
            and node.end_lineno is not None
        ):
            for number in range(node.lineno, node.end_lineno + 1):
                text = lines[number - 1].strip()
                if number == node.lineno:
                    opening = DOCSTRING_RE.match(lines[number - 1])
                    text = opening.group("text").strip() if opening is not None else text
                if number == node.end_lineno and number != node.lineno:
                    text = re.sub(r"(?:\"\"\"|\'\'\')\s*$", "", text).rstrip()
                texts[number - 1] = text
    for token in tokens:
        if token.type == tokenize.COMMENT:
            row, column = token.start
            line = lines[row - 1]
            if not line[:column].strip() and not line.lstrip().startswith(("#!", "#[")):
                texts[row - 1] = token.string[1:].strip()
    return texts


def comment_texts(style: str, lines: Sequence[str]) -> list[str | None]:
    """The text of each line that is a comment (`style` of `COMMENT_STYLES`), None for the others.

    A Python string statement (a docstring) is a comment from its first line to its last; every
    line of a Markdown or text file is text."""
    if style == "python":
        texts_of_python = python_texts(lines)
        if texts_of_python is not None:
            return texts_of_python
    texts: list[str | None] = []
    docstring: str | None = None
    for line in lines:
        if docstring is not None:
            if docstring in line:
                docstring = None
            texts.append(line.strip())
            continue
        if style == "text":
            texts.append(line.strip())
            continue
        if style == "slash":
            comment = COMMENT_RE.match(line)
        else:
            comment = HASH_COMMENT_RE.match(line)
            if comment is None and style == "python":
                opening = DOCSTRING_RE.match(line)
                if opening is not None:
                    if opening.group("quotes") not in opening.group("text"):
                        docstring = opening.group("quotes")
                    texts.append(opening.group("text").strip())
                    continue
        texts.append(None if comment is None else comment.group("text").strip())
    return texts


def named_source(text: str) -> str | None:
    """The other project that `text` names by a path of NCBI's tree or by a name or path."""
    for name, part in NON_NCBI_PATHS:
        if part in text:
            return name
    for name, pattern in FOREIGN_SOURCES:
        if pattern.search(text):
            return name
    return None


def ncbi_path_label(text: str) -> str | None:
    """The other project of a line labelled as NCBI's by a path (`NCBI <path>`), or None."""
    labelled = NCBI_PATH_LABEL_RE.match(text)
    return None if labelled is None else named_source(labelled.group("path"))


def foreign_source(text: str) -> str | None:
    """The other project whose source a citation line cites, or None for NCBI's or LOSAT's.

    A labelled citation (`<label>: <path>...`) is judged by its label and path; any other
    citation line by everything it names."""
    if NCBI_PATH_LABEL_RE.match(text):
        return ncbi_path_label(text)
    labelled = LABELLED_RE.match(text)
    if labelled is None:
        return named_source(text)
    label, path = labelled.group("label").strip(), labelled.group("path")
    source = named_source(label) or named_source(path)
    if source is not None:
        return source
    if NCBI_LABEL_RE.match(label) or path.endswith(".rs") or not CITE_RE.search(path + ":0"):
        return None
    return label


def is_citation(text: str) -> bool:
    """Whether a comment line cites a source: by a path, by an NCBI label, or by naming a
    project of `FOREIGN_SOURCES` or `NON_NCBI_PATHS`."""
    return (
        CITE_RE.search(text) is not None
        or NCBI_CITATION_RE.match(text) is not None
        or named_source(text) is not None
    )


def cites_source(text: str) -> bool:
    """Whether a citation line cites a source file: an NCBI label, a path with a source extension,
    or a path of the pinned trees (not a label over prose, `Result: zero.`)."""
    return (
        NCBI_CITATION_RE.match(text) is not None
        or NCBI_PATH_LABEL_RE.match(text) is not None
        or CITE_RE.search(text) is not None
        or TREE_PATH_RE.search(text) is not None
    )


def is_label_line(text: str) -> bool:
    """Whether a comment line starts with a citation label: `<label>: <path>...`,
    `NCBI reference ...` or `<project> [<version>] reference ...`."""
    return (
        LABELLED_RE.match(text) is not None
        or NCBI_CITATION_RE.match(text) is not None
        or NCBI_PATH_LABEL_RE.match(text) is not None
        or LABEL_START_RE.match(text) is not None
    )


def cited_source(cited: Sequence[tuple[int, str]]) -> tuple[str, str] | None:
    """The other project of the first of `cited` that cites one, with where it is cited."""
    for cited_at, cited_text in cited:
        source = foreign_source(cited_text)
        if source is not None:
            return source, f"cited at :{cited_at}: {cited_text}"
    return None


def load_table(path: Path = TABLE_PATH) -> dict[str, int]:
    """The files of the pinned trees with their line counts (`ncbi_source_lines.tsv`): a line that
    ends with `/` names a directory, the lines after it `<name>\t<lines>` its files."""
    table: dict[str, int] = {}
    directory = ""
    for line in path.read_text(encoding="utf-8").splitlines():
        if not line or line.startswith("#"):
            continue
        if line.endswith("/"):
            directory = line
            continue
        name, _, count = line.partition("\t")
        table[directory + name] = int(count)
    return table


def directories(table: dict[str, int]) -> frozenset[str]:
    """The directories of the table's files (`c++/src/`) and their names without the last
    extension (`.../data/user-column-db` of a BLAST database's files)."""
    out = set()
    for path in table:
        parts = path.split("/")
        for end in range(1, len(parts)):
            out.add("/".join(parts[:end]) + "/")
        stem, dot, _ = path.rpartition(".")
        if dot and "/" in stem:
            out.add(stem)
    return frozenset(out)


def ref_problem(path: str, ranges: str | None, table: dict[str, int], dirs: frozenset[str]) -> str | None:
    """Why a cited path of the pinned trees is wrong, or None."""
    if path.endswith("/"):
        return None if path in dirs else "no such directory"
    if path not in table:
        return None if path + "/" in dirs or path in dirs else "no such file"
    count = table[path]
    for part in (ranges or "").split(","):
        part = part.strip()
        if not part:
            continue
        first, _, last = part.partition("-")
        a, b = int(first), int(last or first)
        if not 1 <= a <= b <= count:
            return f"lines {part} outside the file's {count}"
    return None


def pin_problem(text: str) -> str | None:
    """Why a comment line's citation of NCBI's source does not name the pinned tree, or None."""
    for match in NCBI_LABEL_USE_RE.finditer(text):
        label = " ".join(match.group(0).split())
        if "(" not in label:
            return f"label `{label}` without the pinned tree"
        if label not in PINNED_LABELS:
            return f"label `{label}`"
    if CHECKOUT_PATH_RE.search(text):
        return "path of a local checkout"
    return None


# A comment after code on its line (`let x = 6; // NCBI reference ...`): read for PIN and REF.
TRAILING_COMMENT_RE = re.compile(r"[\w;,{}()\[\]]\s+//[/!]?\s(?P<text>.*)$")


def code_after_citation(text: str) -> int | None:
    """Where code that follows the last cited path of NCBI's tree starts in `text`, or None."""
    last = None
    for last in TREE_PATH_RE.finditer(text):
        pass
    if last is None:
        return None
    tail = CITATION_TAIL_RE.match(text, last.end())
    if tail is None or not CODE_LINE_RE.search(tail.group("tail")):
        return None
    return last.end()


# A C or C++ statement in backticks inside prose (an `if` with its condition, body and
# semicolon): a quote without a fence. Each statement of the span (they end with `;`) has a C
# form; a Rust statement (LOSAT's own code: `let`, `self`, `Some(`, `=>` ...) is not one.
INLINE_STATEMENT = "statement in backticks"
SPAN_RE = re.compile(r"`([^`]+)`")
RUST_STATEMENT_RE = re.compile(
    r"^(?:let|mod|use|pub|fn|impl|struct|enum|type|#\[)\b|=>|::<|&mut\b|\bas [iuf]\d"
    r"|\bSome\(|\bOk\(|\bErr\(|\bself\b|\.unwrap\(|!\(|\?$"
)
C_STATEMENT_RE = re.compile(
    r"^(?:if|while|for|switch)\s*\(|^(?:case\s+[\w:]+\s*:|return\b|break$|continue$|goto\b)"
    r"|->|\+\+|--$|\bm_\w"
    r"|^(?:static\s+|const\s+|unsigned\s+|long\s+)*(?:Int[1248]|Uint[1248]|int|long|short"
    r"|char|double|float|Boolean|size_t|void|vector<[^>]*>|string|[A-Z]\w*::\w+)\b[\s*&]+\w"
    r"|^\w+(?:\[[^\]]*\])*\s*[-+*/%|&^]?=\s*[^=\s]"
    r"|^\w+(?:\.\w+)*\.\w+\s*\([^;]*\)$|^[A-Za-z_]\w*\s*\([^;]*\)$"
)


def inline_statement(span: str) -> bool:
    """Whether the text of a backtick span is one or more C or C++ statements."""
    text = " ".join(span.split())
    if not text.endswith(";"):
        return False
    parts = [part.strip() for part in text.split(";") if part.strip()]
    if not parts or re.fullmatch(r"&[#\w]+", parts[0]):  # an HTML entity
        return False
    if any(RUST_STATEMENT_RE.search(part) for part in parts):
        return False
    return all(C_STATEMENT_RE.search(part) for part in parts)


# A C or C++ expression in backticks inside prose: a quote too. A span may hold a name (an
# identifier, `NStr::StringToInt`, a place such as `a->b[c].d`, a call without arguments such as
# `GetMsg()`), a value, a path, a setting, a command line, a message the binary prints (rule 5:
# observed; `(m_Pos = <n>)`), or LOSAT's own Rust. An expression is an assignment, a comparison,
# arithmetic, a stream `<<`, a logical or bitwise operation, a cast, `++`, `new`, a declaration,
# or a call with arguments. It is C or C++ when it has a mark of C or C++ (`C_EXPRESSION_MARK_RE`)
# and none of Rust (`RUST_MARK_RE`). Without either, in Rust, it is LOSAT's when its names are
# LOSAT's: a name that LOSAT's code nowhere uses (a test of NCBI's `fff` for the end of a string)
# makes it the cited code's, and so does, for a statement under a citation, a name it acts on that
# the file's code does not have (an assignment to `buffer` above the Rust that reads into `text`).
# Names of one or two letters are the symbols of a formula, not names of code.
INLINE_EXPRESSION = "expression in backticks"
# The files whose comments are read for expressions: their own language is not of C's family.
NON_C_FAMILY = frozenset({".rs", ".py", ".toml", ".yml", ".yaml", ".sh"})
STRING_LITERAL_RE = re.compile(r'"(?:[^"\\]|\\.)*"|\'(?:[^\'\\]|\\.)\'')
_SEGMENT = r"~?[A-Za-z_]\w*(?:<[\w:, *]*>)?"
_INDEX = r"\[[^\[\]]*(?:\[[^\[\]]*\][^\[\]]*)*\]"
NAME_SPAN_RE = re.compile(
    rf"^[&*]?\(?{_SEGMENT}(?:(?:::|\.|->){_SEGMENT}|{_INDEX}|\(\))*\)?"
    rf"(?:(?:::|\.|->){_SEGMENT}|{_INDEX}|\(\))*$"
)
VALUE_SPAN_RE = re.compile(r"^-?(?:0x[0-9a-fA-F]+|\d[\d_]*(?:\.\d+)?(?:[eE][-+]?\d+)?)[a-zA-Z]*$|^\"\"$")
PATH_SPAN_RE = re.compile(r"^[\w.+<>*-]*(?:/[\w.+<>*|-]+)+/?(?::[\d,-]+)?$")
# A setting or a row of them (`-evalue=10`, `a=1 b=2`): names and values joined by `=`; a value holds
# no space or `=`, and spaces after it come before the next name or at the end. The same strings as
# `^(?:[\w.\[\]-]*=[^\s=]*\s*)+$`, written so that no two parts take the same letter: there a value and
# the next name could split `-` between them, so a span that does not match took time exponential in
# its `=` (CodeQL py/redos, alert 49).
SETTING_SPAN_RE = re.compile(r"^[\w.\[\]-]*=(?:[^\s=]*(?:\s+[\w.\[\]-]*)?=)*[^\s=]*\s*$")
COMMAND_SPAN_RE = re.compile(
    r"(?:^|\s)--?[a-z][\w-]*(?:\s|$|=)|^(?:cargo|python3?|git|gh|ssh|losat|blast\w*|makeblastdb|npm"
    r"|bash|sh|cd|ls|env|export|make|rustc|perf|time|valgrind)\b(?!\s*\()|^[A-Z][A-Z0-9_]*=\S*(?:\s|$)"
    r"|^set\s+[A-Z_]\w*=|^\S+\.(?:sh|py)(?:\s|$)"
)
PLACEHOLDER_RE = re.compile(r"(?<![\w:])<[a-z][\w +*/-]*>")
# Text, not code: a sentence or a labelled line (`Database: ...`), a diagnostic line
# (`[X_STATS] ...`), a markup tag or a FASTA line,
# a word in quotes (`'abc'`), a label of output (`Length=`, `Frame type= ...`), a power (`2^k`).
# (A word between two operands is text too: `'x' to Int8 (m_Pos = 0)`, `(which is N)`.)
_C_WORDS = (
    "return|new|delete|sizeof|case|const|static|unsigned|signed|long|short|int|char|double|float|void"
    "|bool|struct|enum|class|union|else|goto|throw|typename|template|virtual|inline|extern|volatile"
    "|register|auto|string|operator|if|while|for|switch|do|catch|try|default|public|private|using"
    "|namespace|typedef"
)
TEXT_SPAN_RE = re.compile(
    r"^[A-Z][a-z]+(?: [a-z]+){2}\b|^[A-Z][a-z]+: |^\[[\w-]+\]|^[<>][A-Za-z]|'[^']{2,}'|=$|\w= |\d\^[\w(]"
    r"|\w\^\d"
    # an ASN.1 definition, an alignment row, a Seq-id in FASTA, a manual page (`unlink(2)`), a
    # setting of NCBI's registry (`DATA_LOADERS = none`)
    r"|\s::=\s|^-*[A-Z][A-Z*-]*$|^(?:lcl|gi|ref|gb|emb|dbj|pdb|sp|tr|gnl)\||^[a-z_]\w*\([1-8]\)$"
    r"|^[A-Z][A-Z0-9_]+ = \w+$"
    # a formula: an absolute value, a subscript (`max_{j<i}`, `sum_i`), a half-open interval, `~`
    r"|(?:^|\s)\|[^|\s][^|]*[^|\s]\|(?:\s|$)|\w_\{|\b(?:sum|max|min|prod)_[a-z]\s+[\w(]"
    r"|^\[[^\[\]]*,[^\[\]]*\)$|^\([^()]*,[^()]*\]$|\s~\s"
    rf"|[\w\"')\]]\s+(?!(?:{_C_WORDS})\b)[A-Za-z][a-z]+\s+(?=[\w\"'(])"
)
RUST_ATTRIBUTE_RE = re.compile(
    r"^#?\[?(?:inline|cfg|cfg_attr|derive|allow|deny|warn|expect|doc|target_feature|repr|must_use"
    r"|track_caller|cold|test|ignore)\b"
)
EXPRESSION_RE = re.compile(
    r"(?<![=!<>+\-*/%&|^])=(?![=>])|[-+*/%&|^]=|<<=|>>=|==|!=|<=|>=|\s[<>]\s|<<|>>|&&|\|\|"
    r"|\s\?\s.*\s:\s|\s[-+*/%]\s|\s[&|^]\s|(?:^|[\s(])[!~]\s?[\w(]"
    r"|[\w>]\s*\(\s*[^)\s]|\)\(\s*[^)\s]"  # a call with arguments
    r"|(?<=\w)(?<!\b[Cc])(?:\+\+|--)(?!\w)|(?<![\w+-])(?:\+\+|--)(?=[A-Za-z_])"
    r"|\((?:const\s+)?(?:unsigned\s+|signed\s+)?[A-Za-z_]\w*\s*\**\s*\)\s*[\w(&*-]"  # a cast
    r"|(?:^|[\s(&|^=])~\s*[\w(]|(?<=[\w)\]])[<>](?=-?\d)"
    r"|\bnew\s+[A-Z]|\bsizeof\b"
    # a declaration
    r"|^(?:static\s+|const\s+|unsigned\s+|signed\s+|long\s+)*(?:Int[1248]|Uint[1248]|int|long|short"
    r"|char|double|float|Boolean|size_t|void|string|bool|TSeqPos|TGi|[A-Z]\w*<[^>]*>)"
    r"(?:\s*[*&]+\s*|\s+)[*&]*(?!(?:int|long|short|char|double|float|unsigned|signed|const)\b)[A-Za-z_]\w*\b"
)
# What a statement without a mark acts on: the variable an assignment writes (`buffer` of
# `buffer[length] = 0`), and a call's function and receiver (`reserve` and `str` of
# `str.reserve(count)`); a function of arithmetic (`ceil(x / k)`) is not a call.
ASSIGNMENT_RE = re.compile(r"^(?P<place>[A-Za-z_][\w.]*(?:\[[^\]=]*\])*)\s*(?:[-+*/%&|^]|<<|>>)?=(?![=>])")
CALL_RE = re.compile(r"(?<![\w.:])(?P<receiver>(?:[A-Za-z_]\w*\.)*)(?P<function>[A-Za-z_]\w*)\s*\(\s*[^)\s]")
ARITHMETIC_FUNCTIONS = frozenset(
    "ceil floor exp expl exp2 log log2 log10 ln sqrt max min abs fabs pow round trunc sat erf erfc "
    "sum acos asin atan cos sin tan".split()
)
# Arithmetic without spaces counts in a span with a mark of C or C++ (`3 - i%4`; not the hyphen
# of `Site-annot-set`).
TIGHT_ARITHMETIC_RE = re.compile(r"[\w)\]][+*/%&|^][\w(]|[\d)\]]-[\w(]|[\w)\]]-[\d(]")
C_EXPRESSION_MARK_RE = re.compile(
    r"(?<=[\w)\]])->(?=[\w(*])"  # a member through a pointer (Rust's `->` has spaces)
    r"|(?<![\w.])m_[A-Za-z]\w*"  # a member
    r"|(?<![\w.])[kef][A-Z][a-z]\w*"  # a constant, an enum value, a flag
    r"|(?:\.|->)[A-Z]\w*\("  # a method in CamelCase
    r"|(?<![\w:.!])[A-Z][A-Z0-9]*(?:_[A-Z0-9]+)+\s*\(|(?<![\w:.!])[A-Z][A-Z0-9]{2,}\("  # a macro
    r"|(?<![\w:.!])(?:MIN|MAX|ABS|ASSERT|_ASSERT|ITERATE|NON_CONST_ITERATE)\s*\("
    r"|(?<!\w)[A-Za-z][a-z0-9]*_[A-Z][A-Za-z0-9]*[a-z]\w*\s*\("  # Prefix_Name(, x_Name(
    r"|(?<![\w:.])[A-Z][a-z0-9]+[A-Z]\w*\s*\("  # a function in CamelCase
    r"|\b[CNS][A-Z][a-z]\w*(?:<[^>]*>)?::[A-Z]\w*\s*\(|\bC[A-Z][a-z]\w*\s*\("  # a class's
    r"|<<\s*(?:\"|'|\b(?:endl|NcbiEndl|Warning|Error|Info|Fatal|Critical|Trace|Note|MSerial\w*|setw"
    r"|setprecision|flush|fixed|scientific)\b)|\b(?:cout|cerr|ostr|os|out|\w*[Ss]tream\w*)\s*<<"
    r"|\((?:const\s+)?(?:unsigned\s+|signed\s+)?(?:Int[1248]|Uint[1248]|int|long|short|char|double"
    r"|float|Boolean|size_t|void|TSeqPos|\w+_t)\s*\**\s*\)\s*[\w(&*-]"  # a cast
    r"|^(?:static\s+|const\s+|unsigned\s+|signed\s+)*(?:Int[1248]|Uint[1248]|int|long|short|char"
    r"|double|float|Boolean|size_t|void|string|bool|\w+_t)(?:\s*[*&]+\s*|\s+)[*&]*(?!(?:int|long|short|char|double|float|unsigned|signed|const)\b)[A-Za-z_]\w*\b"
    r"|\b(?:sizeof|malloc|calloc|realloc|memcpy|memset|memmove|strlen|strcmp|strncmp|strcpy|sprintf"
    r"|snprintf|fprintf|printf|static_cast|reinterpret_cast|const_cast|dynamic_cast|NULL|nullptr"
    r"|TRUE|FALSE|INT\d*_M(?:IN|AX)|UINT\d*_MAX|CHAR_BIT|DBL_\w+|FLT_\w+)\b"
    r"|\bnew\s+[A-Z]\w*|\bconst\s+\w+\s*[&*]|!!|\bios(?:_base)?::|::npos\b"
    r"|\.(?:c_str|push_back|emplace_back|substr|empty|find_first_of|find_last_of|find_first_not_of"
    r"|tellg|seekg|tellp|seekp|rdbuf)\(|(?<=\w)(?<!\b[Cc])\+\+(?!\w)|(?<![\w+-])\+\+(?=[A-Za-z_])"
    r"|\bstd::(?:string|vector|map|set|pair|ostream|istream|cout|cerr|endl)\b"
    r"|^(?:\}\s*)?(?:if|while|for|switch|catch)\s*\(|^\}?\s*else\b|\btry\s*\{"
    r"|(?<=[\w)\]])\s+\?\s+\S.*?\s:\s+\S"  # a ternary
    r"|(?<![\w.])(?:0[xX][0-9A-Fa-f]+|\d+)(?:[uU](?:ll|LL|l|L)?|(?:ll|LL|l|L)[uU]?)\b"  # 1U, 8191L
    r"|(?<=\w)--(?![\w-])"  # a decrement
    # a call of `operator()`, and the call of a template given a type of C (`min` of `int`)
    r"|\boperator\s*\(\)\s*\(|<\s*(?:unsigned\s+)?(?:int|long|short|char|double|float|size_t"
    r"|Int[1248]|Uint[1248]|TSeqPos)\b[^<>]*>\s*\("
)
# Rust has `!` for a bitwise not, so `~` there is C's (Python has `~` too).
BITWISE_NOT_RE = re.compile(r"(?:^|[\s(&|^=])~\s*(?:\(|\w+\b(?!:))")
# The arguments of a template (`map<CSeq_id_Handle, CConstRef<CSeq_loc>>` is a name).
TEMPLATE_ARGS_RE = re.compile(r"(?<=\w)<[\w:, *&]*>")
# A call of a name in camel case or in capitals: C or C++'s, unless it has arguments and LOSAT's
# Rust has the name (a tuple struct or an enum variant).
CAMEL_CASE_CALL_RE = re.compile(r"[A-Z][a-z0-9]+[A-Z]\w*\s*\(|[A-Z][A-Z0-9]{2,}\(")
# In Rust and Python a variable is in snake case: a name in camel case (`alnLen`) is C or C++'s,
# and a statement that uses one is theirs (arithmetic on one may be reasoning about its value).
CAMEL_CASE_VARIABLE_RE = re.compile(r"(?<![\w.:'\"])[a-z][a-z0-9]*[A-Z]\w*(?![\w(])")
RUST_MARK_RE = re.compile(
    r"^(?:let|mod|use|pub|fn|impl|struct|enum|type|match|loop|unsafe|#\[)\b|=>|::<|&mut\b"
    r"|\bas (?:[iuf]\d+|usize|isize)\b|\bSome\(|\bOk\(|\bErr\(|\bNone\b|\bself\b|\bSelf\b"
    r"|\.(?:unwrap|expect)\(|\w!\(|\w!\[|\?(?:$|[.;)])|\w\.\.=?\w|\.\.=\w|\|\w*\||\bfn\b"
    r"|\b(?:Vec|Option|Result|Box)<|&\[|\b[iuf](?:8|16|32|64|128|size)\b"
    r"|\.(?:len|is_empty|iter|into_iter|map|collect|as_ref|as_slice|to_vec|to_string|clone"
    r"|saturating_\w+|wrapping_\w+|checked_\w+)\(|'_\b|[<&]'[a-z]\w*\b|\b(?:Fn|FnMut|FnOnce)\("
    r"|\b(?:impl|dyn)\s|\b(?:core|std|alloc)::[a-z_]+::"
)
# Words of Rust and its types: not the names of a file's code.
RUST_WORDS = frozenset(
    "as break const continue crate else enum extern false fn for if impl in let loop match mod move "
    "mut pub ref return self Self static struct super trait true type unsafe use where while async "
    "await dyn u8 u16 u32 u64 u128 usize i8 i16 i32 i64 i128 isize f32 f64 bool char str String Vec "
    "Option Some None Ok Err Result Box".split()
)
IDENTIFIER_RE = re.compile(r"[A-Za-z_]\w*")
# A name in an expression (not the digits and suffix of a number: `0x11`, `1e-5`).
NAME_RE = re.compile(r"(?<!\w)[A-Za-z_]\w*")
# Words that are no name of code: the keywords of Rust and Python, functions and constants of
# arithmetic (`Lambda` of a score's statistics).
NOT_NAMES = RUST_WORDS | frozenset(
    "true false None True False and or not in is if else lambda Lambda NaN nan inf Inf INF".split()
)


def inline_expression(
    span: str,
    names: Callable[[], frozenset[str]],
    cited: bool,
    snake_case: bool = True,
    losat: Callable[[], frozenset[str]] = frozenset,
    rust: bool = False,
) -> bool:
    """Whether the text of a backtick span is a C or C++ expression. `names` gives the names of
    the file's code and the functions LOSAT defines, `losat` every name of LOSAT's code; `cited`
    says whether the span's comment cites a source; `snake_case`, that the file's language names
    variables in snake case (Rust, Python); `rust`, that it is Rust."""
    text = " ".join(span.split()).rstrip(";").strip()
    bare = STRING_LITERAL_RE.sub('""', text)
    plain = bare
    while TEMPLATE_ARGS_RE.search(plain):
        plain = TEMPLATE_ARGS_RE.sub("", plain)
    if (
        not bare
        or VALUE_SPAN_RE.match(bare)
        or NAME_SPAN_RE.match(bare)
        or NAME_SPAN_RE.match(plain)
        or PATH_SPAN_RE.match(bare)
        or SETTING_SPAN_RE.match(text)
        or RUST_ATTRIBUTE_RE.match(bare)
        or COMMAND_SPAN_RE.search(text)
        or PLACEHOLDER_RE.search(bare)
        or TEXT_SPAN_RE.search(bare)
        or re.match(r"^operator\W+$", bare)
        or cites_source(text)
        or PROSE_RE.search(bare)
    ):
        return False
    marks = [mark.group(0) for mark in C_EXPRESSION_MARK_RE.finditer(bare)]
    calls = [mark for mark in marks if CAMEL_CASE_CALL_RE.fullmatch(mark)]
    marked = (
        len(calls) < len(marks)
        or (rust and BITWISE_NOT_RE.search(bare) is not None)
        or any(
            not re.search(re.escape(call) + r"\s*[^)\s]", bare) or call.rstrip("( ") not in names()
            for call in calls
        )
    )
    if not (EXPRESSION_RE.search(bare) or (marked and TIGHT_ARITHMETIC_RE.search(bare))):
        return False
    if RUST_MARK_RE.search(bare):
        return False  # LOSAT's Rust
    if marked:
        return True
    # Without a mark: a statement that uses a variable in camel case, or that acts on a name
    # foreign to the file's code under a citation, is the cited code's.
    acted = acted_names(bare)
    if acted and snake_case and CAMEL_CASE_VARIABLE_RE.search(bare):
        return True
    if acted and cited and any(name not in names() for name in acted if len(name) > 1):
        return True
    # In Rust, any expression with a name that LOSAT's code nowhere uses is the cited code's, and
    # one with LOSAT's names is LOSAT's code or reasoning about its values (`last_hit <= diag_offset`).
    if not rust:
        return False
    foreign = {
        name
        for name in NAME_RE.findall(bare)
        if len(name) > 2 and name not in NOT_NAMES and name not in ARITHMETIC_FUNCTIONS
    }
    if not foreign:
        return False
    known = names()
    return any(name not in known and name not in losat() for name in foreign)


def acted_names(bare: str) -> list[str]:
    """The names that a statement without a mark acts on (`ASSIGNMENT_RE`, `CALL_RE`); none for
    a comparison or arithmetic."""
    acted = []
    assignment = ASSIGNMENT_RE.match(bare)
    if assignment is not None:
        acted.append(IDENTIFIER_RE.match(assignment.group("place")).group(0))
    for call in CALL_RE.finditer(bare):
        if call.group("function") in ARITHMETIC_FUNCTIONS:
            continue
        acted.append(call.group("function"))
        if call.group("receiver"):
            acted.append(IDENTIFIER_RE.match(call.group("receiver")).group(0))
    return acted


def code_names(style: str, lines: Sequence[str], texts: Sequence[str | None]) -> frozenset[str]:
    """The names of a file's code: its identifiers outside comments, docstrings and strings."""
    if style == "python":
        try:
            tokens = list(tokenize.generate_tokens(io.StringIO("\n".join(lines) + "\n").readline))
        except (SyntaxError, tokenize.TokenError, IndentationError, ValueError):
            tokens = []
        return frozenset(t.string for t in tokens if t.type == tokenize.NAME)
    return frozenset(IDENTIFIER_RE.findall(code_lines(lines, texts))) - RUST_WORDS


# The names of LOSAT's Rust anywhere in the tree: the functions it defines (`fn name`) and the
# names in camel case of its code (types, tuple structs, enum variants). A call of one is LOSAT's.
FN_DEFINITION_RE = re.compile(r"\bfn\s+(?P<name>[a-z_]\w*)")
CAMEL_CASE_NAME_RE = re.compile(r"\b[A-Z][a-z0-9]+[A-Z]\w*\b")


def code_lines(lines: Sequence[str], texts: Sequence[str | None]) -> str:
    """The code of a file without its comments, docstrings and strings (line comments after
    code too), as one text."""
    out = []
    for line, text in zip(lines, texts):
        if text is None:
            line = STRING_LITERAL_RE.sub(" ", line)
            out.append(line.split("//", 1)[0] if "//" in line else line)
    return "\n".join(out)


def rust_names(root: Path) -> tuple[frozenset[str], frozenset[str]]:
    """The names of LOSAT's Rust under `root`: those a call may use (`FN_DEFINITION_RE`,
    `CAMEL_CASE_NAME_RE`), and every identifier of its code."""
    functions: set[str] = set()
    identifiers: set[str] = set()
    for path in source_files(root):
        if path.suffix == ".rs":
            lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
            code = code_lines(lines, comment_texts("slash", lines))
            functions.update(match.group("name") for match in FN_DEFINITION_RE.finditer(code))
            functions.update(CAMEL_CASE_NAME_RE.findall(code))
            identifiers.update(IDENTIFIER_RE.findall(code))
    return frozenset(functions), frozenset(identifiers)


def inline_statement_findings(
    rel: str,
    texts: Sequence[str | None],
    lines: Sequence[str] = (),
    style: str = "slash",
    functions: frozenset[str] = frozenset(),
    identifiers: frozenset[str] = frozenset(),
) -> list[Finding]:
    """QUOTE findings for statements and expressions in backticks; a span may go on to the next
    comment line, so each comment run (up to a code line, an empty comment line or a fence) is
    read as one text. Fenced lines are left to the fence checks."""
    findings: list[Finding] = []
    run: list[tuple[int, str]] = []
    in_fence = False
    known: list[frozenset[str]] = []
    # In C, C++, JavaScript and TypeScript, an expression in backticks is the file's own syntax.
    non_c_family = Path(rel).suffix in NON_C_FAMILY
    snake_case = Path(rel).suffix in (".rs", ".py")
    rust = Path(rel).suffix == ".rs"

    def names() -> frozenset[str]:
        if not known:
            known.append(code_names(style, lines, texts) | functions)
        return known[0]

    def flush() -> None:
        joined, starts = "", []
        for number, text in run:
            starts.append((len(joined), number))
            joined += text.strip() + " "
        cited: bool | None = None
        for match in SPAN_RE.finditer(joined):
            span = match.group(1)
            source = INLINE_STATEMENT if inline_statement(span) else None
            if source is None:
                if cited is None:
                    cited = any(cites_source(text) for _, text in run)
                if non_c_family and inline_expression(
                    span, names, cited, snake_case, lambda: identifiers, rust
                ):
                    source = INLINE_EXPRESSION
            if source is not None:
                number = [n for offset, n in starts if offset <= match.start()][-1]
                findings.append(Finding(rel, number, "QUOTE", source, match.group(0)))
        run.clear()

    for number, text in enumerate(texts, start=1):
        if text is None:  # a code line ends the run, and a fence with its comment
            flush()
            in_fence = False
            continue
        if FENCE_RE.match(text):
            flush()
            in_fence = not in_fence
            continue
        if in_fence:
            continue
        if not text:
            flush()
            continue
        run.append((number, text))
    flush()
    return findings


# Unbackticked C or C++ statements and expressions in comment text outside backticks (line
# comments, docstrings and `/* ... */` block comments in `NON_C_FAMILY` files), cited or not.
# Every pattern below is linear-time (no adjacent or nested quantifiers competing for the same
# characters; balanced parentheses are scanned in Python by `_closing_paren`).
UNBACKTICKED_QUOTE = "unbackticked C/C++ code"
DOUBLE_QUOTED_RE = re.compile(r'"[^"\\]*(?:\\.[^"\\]*)*"')
QUOTE_PREFIX_RE = re.compile(
    r"^(?:(?:NCBI(?:\s+(?:BLAST|lines?\s+\d+(?:-\d+)?))?|C(?:\+\+)?)\s*:\s*|\[C(?:\+\+)?\]\s*)"
)
SQL_DDL_RE = re.compile(r"^\s*CREATE\s+(?:UNIQUE\s+)?(?:INDEX|TABLE)\b", re.I)
OBSERVED_MPOS_MSG_RE = re.compile(r"'[^']+'\s+to\s+\w+\s+\(m_Pos\s*=")
CONTROL_START_RE = re.compile(r"(?:^|[;{}]\s*|\belse\s+)(?P<kw>if|while|for|switch|catch)\s*\(")
CONTROL_BODY_RE = re.compile(
    r"\s*(?:\{|;|return\b|break\b|continue\b|goto\b|else\b|\+\+|--"
    r"|\*?[A-Za-z_]\w*(?:\[[^\[\]]*\])*(?:->|\.[A-Za-z_]\w*|\s*(?:[-+*/%&|^]?=|\+\+|--|\()))"
)
CONTROL_COND_MARK_RE = re.compile(r"&&|\|\||==|!=|<=|>=|(?<!\w)[A-Za-z_]\w*->|\bm_[A-Za-z]\w*")
UNBACKTICKED_STRUCT_RE = re.compile(
    r"(?:^|[;}]\s*)else\s*\{"
    r"|(?:^|[;{}]\s*)typedef\s+(?:struct|enum|union)\b"
    r"|(?:^|[;{}]\s*)(?:static\s+|const\s+)?struct\s*(?:[A-Za-z_]\w*\s+)?\{[^}]*\}\s*[A-Za-z_]\w*\s*;"
)
SEMI_BRACE_END_RE = re.compile(
    r";\s*\}\s*(?://[^\n]*|/\*[^*]*(?:\*(?!/)[^*]*)*\*/|\([^()]*\))?$"
)
UNBACKTICKED_STMT_RE = re.compile(
    r"(?:^|[;{}]\s*)(?:break|continue)\s*;"
    r"|(?:^|[;{}]\s*)return\s*;\s*(?:$|/\*|//|\})"
    r"|(?:^|[;{}]\s*)return\s+[^\s;.][^;.]*;\s*(?:$|/\*|//|\})"
    r"|(?:^|[;{}]\s*)(?:\+\+|--)[A-Za-z_]\w*\s*;"
    r"|(?:^|[;{}]\s*)[A-Za-z_]\w*(?:\[[^\[\]]*\])*(?<!\b[Cc])(?:\+\+|--)\s*;"
)
UNBACKTICKED_DECL_RE = re.compile(
    r"(?:^|[;{}]\s*)(?:static\s+|const\s+|unsigned\s+|signed\s+)*"
    r"(?:Int[1248]|Uint[1248]|Boolean|TSeqPos|TGi|size_t|Blast[A-Z]\w*)"
    r"(?:\s*[*&]+\s*[A-Za-z_]\w*|\s+[A-Za-z_]\w*\s*[=,;\[])"
)
# Standalone assignment or call statement with `;` (e.g. a variable assignment or function call).
STANDALONE_ASSIGN_SEMI_RE = re.compile(
    r"^(?:\*|\(\s*\*?[A-Za-z_]\w*\s*\)|[A-Za-z_]\w*(?:(?:\.|->)[A-Za-z_]\w*|\[[^\[\]]*\])*)"
    r"\s*(?:[-+*/%&|^]|<<|>>)?=\s*[^=\s][^;]*;\s*(?:$|/\*|//|\})"
)
STANDALONE_CALL_HEAD_RE = re.compile(r"^[A-Za-z_]\w*(?:(?:::|\.|->)[A-Za-z_]\w*)*\(")
TRAILING_SEMI_END_RE = re.compile(r"^\s*;\s*(?:$|/\*|//|\})")
ARROW_MEMBER_RE = re.compile(r"(?<!\w)[A-Za-z_]\w+->[A-Za-z_]\w+")
ASSIGN_HEAD_RE = re.compile(
    r"(?:^|[;{}]\s*)(?:\*?[A-Za-z_]\w*(?:(?:\.|->)[A-Za-z_]\w*|\[[^\[\]]*\])*|\(\s*\*?[A-Za-z_]\w*\s*\))"
    r"\s*(?:[-+*/%&|^]|<<|>>)?=(?![=>])\s*(?P<rhs>\S.*)$"
)
UNBACKTICKED_RUST_RE = re.compile(
    r"^(?:let|mod|use|pub|fn|impl|enum|type|match|loop|unsafe|#\[)\b"
    r"|^struct\s+[A-Z]\w*(?:\s*[({<;]|\s*$)"
    r"|=>|::<|&mut\b|\bas (?:[iuf]\d+|usize|isize)\b"
    r"|\bSome\(|\bOk\(|\bErr\(|\bself\b|\bSelf\b"
    r"|\.(?:unwrap|expect)\(|\w!\(|\w!\[|\w\.\.=?\w|\.\.=\w|\bfn\b"
    r"|\b(?:Vec|Option|Result|Box)<|&\[|\b[iuf](?:8|16|32|64|128|size)\b"
    r"|\.(?:len|is_empty|iter|into_iter|map|collect|as_ref|as_slice|to_vec|to_string|clone"
    r"|saturating_\w+|wrapping_\w+|checked_\w+)\(|'_\b|[<&]'[a-z]\w*\b|\b(?:Fn|FnMut|FnOnce)\("
    r"|\b(?:impl|dyn)\s|\b(?:core|std|alloc)::[a-z_]+::"
)
UNBACKTICKED_EXPR_RE = re.compile(
    # pointer member in assignment / comparison / inc / dec / call without space
    r"(?<!\w)[A-Za-z_]\w*(?:\[[^\[\]]*\])*->[A-Za-z_]\w+(?:\[[^\[\]]*\])*\s*(?:[-+*/%&|^]?=|==|!=|<=|>=|\+\+|--)"
    r"|(?<!\w)[A-Za-z_]\w*(?:\[[^\[\]]*\])*->[A-Za-z_]\w+\([^)\s][^)]*\)"
    # assignment whose right-hand side reads through a C pointer `->`
    r"|(?:\w|\])\s*(?:[-+*/%&|^]|<<|>>)?=(?![=>])\s*\*?[A-Za-z_]\w*(?:\[[^\[\]]*\])*->[A-Za-z_]\w+"
    # m_* member in expression (observed `(m_Pos = <n>)` / `'x' to T (m_Pos = 0)` are filtered)
    r"|(?<![\w.])m_[A-Za-z]\w*(?:\[[^\[\]]*\]|\.[A-Za-z_]\w*\([^)]*\)|\s*(?:[-+*/%&|^]?=|==|!=|<=|>=|<<|\+\+|--))"
    # C++ stream insertion `<<`
    r"|\b(?:m_[A-Za-z]\w*|cout|cerr|ostr|os|out|(?:o|i|io|f|of|if|string|ostring|istring|out_|err_|log_)?[Ss]tream[A-Za-z0-9_]*)\s*<<\s*(?:\"\"|'|[A-Za-z_])"
    r"|<<\s*(?:\"\"|\b(?:endl|NcbiEndl|Warning|Error|Info|Fatal|Critical)\b)"
    # C++ namespace::Call(args) or C++ STL method calls (no space before `(`)
    r"|\b(?:NStr|C[A-Z][a-z]\w*|std)::[A-Za-z_]\w*(?:<[^<>]*>)?\([^)\s][^)]*\)"
    r"|\.(?:c_str|push_back|emplace_back|substr|find_first_of|find_last_of|find_first_not_of|tellg|seekg|tellp|seekp|rdbuf)\([^)]*\)"
    # C cast attached to an expression (not English prose like `Int4` before "of" or "cast")
    r"|\((?:const\s+)?(?:unsigned\s+|signed\s+)?(?:Int[1248]|Uint[1248]|size_t|TSeqPos|double|float|int|long)\s*\*?\)"
    r"(?:\(|[A-Za-z_*]\w*(?:->|\.|\[|\()|\s+(?:ceil|floor|[A-Za-z_]\w*->|\*[A-Za-z_]))"
    # prefix increment or postfix increment/decrement at statement boundary (not `C++` or `--flag`)
    r"|(?<![\w+-])\+\+[A-Za-z_]\w*"
    r"|(?:^|[;{}]\s*)[A-Za-z_]\w*(?:\[[^\[\]]*\])*(?<!\b[Cc])(?:\+\+|--)(?=\s*[,;)}]|$)"
    # array equality/inequality: `a[i] == b[i]`
    r"|(?<!\w)[A-Za-z_]\w*\[[^\[\]]+\]\s*(?:==|!=)\s*[A-Za-z_]\w*\[[^\[\]]+\]"
    # C/NCBI function or macro call with arguments (no space before `(`), or `sizeof` operator
    r"|(?<![\w:.])(?:[sx]_[A-Z]\w*|BLAST_[A-Z]\w*|Blast_[A-Z]\w*|Blast[A-Z]\w*|LOG_POST|ERR_POST\w*|ITERATE|NON_CONST_ITERATE)\([^)\s][^)]*\)"
    r"|\bsizeof\s*\([^)\s][^)]*\)"
)
UNBACKTICKED_METHOD_CALL_RE = re.compile(r"(?:\.|->)[A-Z]\w*\([^)\s][^)]*\)")
UNBACKTICKED_CAMEL_CALL_RE = re.compile(r"(?<![\w:.])(?P<name>[A-Z][a-z0-9]+[A-Z]\w*)\([^)\s][^)]*\)")
PREFIXED_ASSIGN_RE = re.compile(
    r"^(?:\*?[A-Za-z_]\w*(?:(?:\.|->)[A-Za-z_]\w*|\[[^\[\]]*\])*|\(\s*\*?[A-Za-z_]\w*\s*\))"
    r"\s*(?:[-+*/%&|^]|<<|>>)?=(?![=>])\s*\S"
)
FOR_INIT_CODE_RE = re.compile(r"^\s*$|==|!=|<=|>=|=|\+\+|--|->|<|>")


def _is_c_for_header(inside: str) -> bool:
    """Whether `inside` (the text after `for (`) starts with a C `for` initializer before `;`."""
    if ";" not in inside:
        return False
    init = inside.split(";", 1)[0]
    return FOR_INIT_CODE_RE.search(init) is not None


def _closing_paren(text: str, open_idx: int) -> int | None:
    """Index of the `)` matching `text[open_idx] == '('`, or None if unclosed."""
    depth = 0
    for idx in range(open_idx, len(text)):
        ch = text[idx]
        if ch == "(":
            depth += 1
        elif ch == ")":
            depth -= 1
            if depth == 0:
                return idx
    return None


def _has_c_control_statement(clean: str) -> bool:
    """Whether `clean` (double-quoted strings masked) contains an unbackticked C/C++ control
    statement (`if`, `while`, `for`, `switch`, `catch`)."""
    for match in CONTROL_START_RE.finditer(clean):
        kw = match.group("kw")
        open_idx = match.end() - 1
        close_idx = _closing_paren(clean, open_idx)
        if close_idx is None:
            inside = clean[open_idx + 1 :]
            if kw == "for" and _is_c_for_header(inside):
                return True
            if kw in ("if", "while", "catch") and CONTROL_COND_MARK_RE.search(inside):
                return True
            continue
        inside = clean[open_idx + 1 : close_idx]
        rest = clean[close_idx + 1 :]
        if kw == "for":
            if _is_c_for_header(inside):
                return True
            continue
        if kw == "switch":
            if rest.lstrip().startswith("{"):
                return True
            continue
        if CONTROL_BODY_RE.match(rest) or CONTROL_COND_MARK_RE.search(inside):
            return True
    return False


def _is_standalone_call_semi(clean: str) -> bool:
    """Whether `clean` is a standalone C/C++ call ending with `;` (supporting nested parens)."""
    head = STANDALONE_CALL_HEAD_RE.match(clean)
    if head is None:
        return False
    close_idx = _closing_paren(clean, head.end() - 1)
    return close_idx is not None and TRAILING_SEMI_END_RE.match(clean[close_idx + 1 :]) is not None


def unbackticked_c_code(
    outside: str,
    raw: str,
    *,
    style: str,
    names: Callable[[], frozenset[str]],
) -> bool:
    """Whether comment text outside backticks (`outside`, from comment line `raw`) holds an
    unbackticked C or C++ statement or expression."""
    if not outside:
        return False
    if is_citation(raw) and code_after_citation(raw) is None:
        return False
    if style == "python" and outside.startswith(">>>"):
        return False
    if PLACEHOLDER_RE.search(outside) or SQL_DDL_RE.match(outside):
        return False
    had_prefix = False
    body = outside
    while True:
        pref = QUOTE_PREFIX_RE.match(body)
        if pref is None:
            break
        had_prefix = True
        body = body[pref.end() :].strip()
    if not body:
        return False
    clean = DOUBLE_QUOTED_RE.sub('""', body) if '"' in body else body
    if OBSERVED_MPOS_MSG_RE.search(clean):
        return False
    if UNBACKTICKED_STRUCT_RE.search(clean):
        return True
    if UNBACKTICKED_RUST_RE.search(body):
        return False
    if _has_c_control_statement(clean):
        return True
    if SEMI_BRACE_END_RE.search(clean):
        return True
    if not PROSE_RE.search(clean) and not clean.startswith("- "):
        if (
            UNBACKTICKED_DECL_RE.search(clean)
            or STANDALONE_ASSIGN_SEMI_RE.match(clean)
            or _is_standalone_call_semi(clean)
            or UNBACKTICKED_STMT_RE.search(clean)
        ):
            return True
    if UNBACKTICKED_EXPR_RE.search(clean):
        return True
    assign = ASSIGN_HEAD_RE.search(clean)
    if assign is not None and ARROW_MEMBER_RE.search(assign.group("rhs")):
        return True
    if had_prefix and PREFIXED_ASSIGN_RE.match(clean) and not PROSE_RE.search(clean):
        return True
    if style != "python":
        if UNBACKTICKED_METHOD_CALL_RE.search(clean):
            return True
        for match in UNBACKTICKED_CAMEL_CALL_RE.finditer(clean):
            if match.group("name") not in names():
                return True
    return False


def strip_backticks_in_run(run: Sequence[tuple[int, str]]) -> list[tuple[int, str, str]]:
    """For each `(line_no, raw_text)` in a comment run, return `(line_no, raw_text, outside)`
    where text inside backtick spans (including spans that cross comment lines) is replaced by
    a neutral token so adjacent words never merge across a removed span."""
    out: list[tuple[int, str, str]] = []
    in_bt = False
    for number, text in run:
        if not in_bt and "`" not in text:
            out.append((number, text, text.strip()))
            continue
        parts: list[str] = []
        for ch in text:
            if ch == "`":
                in_bt = not in_bt
                parts.append(" _BT_ " if not in_bt else " ")
            elif not in_bt:
                parts.append(ch)
        out.append((number, text, "".join(parts).strip()))
    return out


def block_comments(lines: Sequence[str]) -> list[tuple[int, str]]:
    """Every non-empty line inside `/* ... */` block comments of a slash-style file (skipping
    `//` comments, string/raw-string literals and char literals, supporting nested `/* ... */`),
    as `(1-based line_no, text)` with leading `*` decoration stripped."""
    text = "\n".join(lines)
    if "/*" not in text:
        return []
    n = len(text)
    i = 0
    line_no = 1
    out: list[tuple[int, str]] = []
    while i < n:
        ch = text[i]
        if ch == "\n":
            line_no += 1
            i += 1
            continue
        if text.startswith("//", i):
            j = text.find("\n", i)
            i = n if j < 0 else j
            continue
        if text.startswith("/*", i):
            start_line = line_no
            depth = 1
            i += 2
            buf_start = i
            while i < n and depth > 0:
                if text.startswith("/*", i):
                    depth += 1
                    i += 2
                elif text.startswith("*/", i):
                    depth -= 1
                    if depth == 0:
                        body = text[buf_start:i]
                        i += 2
                        break
                    i += 2
                else:
                    if text[i] == "\n":
                        line_no += 1
                    i += 1
            else:
                body = text[buf_start:i]
            for offset, raw_line in enumerate(body.split("\n")):
                cleaned = raw_line.strip()
                if cleaned.startswith("*"):
                    cleaned = cleaned.lstrip("*").strip()
                if cleaned:
                    out.append((start_line + offset, cleaned))
            continue
        raw = re.match(r'[bc]?r(#*)"', text[i : i + 260])
        if raw and (i == 0 or not (text[i - 1].isalnum() or text[i - 1] == "_")):
            close = '"' + raw.group(1)
            j = text.find(close, i + raw.end())
            end = n if j < 0 else j + len(close)
            line_no += text.count("\n", i, end)
            i = end
            continue
        if ch == '"':
            j = i + 1
            while j < n and text[j] != '"':
                if text[j] == "\n":
                    line_no += 1
                j += 2 if text[j] == "\\" else 1
            i = min(j + 1, n)
            continue
        char = re.match(r"'(?:\\.[^'\n]*|[^'\\\n])'", text[i : i + 16])
        if char:
            i += char.end()
            continue
        i += 1
    return out


def unbackticked_findings(
    rel: str,
    texts: Sequence[str | None],
    lines: Sequence[str],
    style: str,
    functions: frozenset[str],
    identifiers: frozenset[str],
    already_quoted: frozenset[int],
) -> list[Finding]:
    """QUOTE findings for unbackticked C/C++ statements and expressions in `NON_C_FAMILY` files
    (line comments, docstrings, trailing comments and `/* ... */` block comments)."""
    if Path(rel).suffix not in NON_C_FAMILY:
        return []
    findings: list[Finding] = []
    seen_lines: set[int] = set(already_quoted)
    known: list[frozenset[str]] = []

    def names() -> frozenset[str]:
        if not known:
            known.append(code_names(style, lines, texts) | functions)
        return known[0]

    def record(number: int, raw: str, outside: str) -> None:
        if number in seen_lines:
            return
        if unbackticked_c_code(outside, raw, style=style, names=names):
            seen_lines.add(number)
            findings.append(Finding(rel, number, "QUOTE", UNBACKTICKED_QUOTE, raw))

    def check_single(number: int, raw: str) -> None:
        _, _, outside = strip_backticks_in_run([(number, raw)])[0]
        record(number, raw, outside)

    run: list[tuple[int, str]] = []
    in_fence = False

    def flush() -> None:
        for number, raw, outside in strip_backticks_in_run(run):
            record(number, raw, outside)
        run.clear()

    for number, text in enumerate(texts, start=1):
        if text is None:
            flush()
            in_fence = False
            if style == "slash" and number not in seen_lines:
                no_str = DOUBLE_QUOTED_RE.sub('""', lines[number - 1])
                trailing = TRAILING_COMMENT_RE.search(no_str)
                if trailing is not None:
                    check_single(number, trailing.group("text").strip())
            continue
        if FENCE_RE.match(text):
            flush()
            in_fence = not in_fence
            continue
        if in_fence:
            continue
        if not text:
            flush()
            continue
        run.append((number, text))
    flush()

    if style == "slash":
        bc_lines = block_comments(lines)
        if bc_lines:
            snake_case = Path(rel).suffix in (".rs", ".py")
            rust = Path(rel).suffix == ".rs"
            for number, raw in bc_lines:
                if number in seen_lines:
                    continue
                for match in SPAN_RE.finditer(raw):
                    span = match.group(1)
                    if inline_statement(span):
                        seen_lines.add(number)
                        findings.append(Finding(rel, number, "QUOTE", INLINE_STATEMENT, match.group(0)))
                        break
                    if inline_expression(span, names, False, snake_case, lambda: identifiers, rust):
                        seen_lines.add(number)
                        findings.append(Finding(rel, number, "QUOTE", INLINE_EXPRESSION, match.group(0)))
                        break
                check_single(number, raw)

    return findings


def citation_findings(
    rel: str, number: int, text: str, table: dict[str, int], dirs: frozenset[str]
) -> list[Finding]:
    """PIN and REF of one comment's text."""
    findings: list[Finding] = []
    problem = pin_problem(text)
    if problem is not None:
        findings.append(Finding(rel, number, "PIN", problem, text))
    for match in TREE_PATH_RE.finditer(text):
        problem = ref_problem(match.group("path"), match.group("ranges"), table, dirs)
        if problem is not None:
            findings.append(Finding(rel, number, "REF", problem, match.group(0)))
    return findings


def text_indent(line: str) -> int:
    marker = MARKER_RE.match(line)
    rest = line[marker.end() :] if marker else line
    return len(rest) - len(rest.lstrip(" "))


def check_text(
    rel: str,
    lines: Sequence[str],
    style: str = "slash",
    table: dict[str, int] | None = None,
    dirs: frozenset[str] | None = None,
    functions: frozenset[str] = frozenset(),
    identifiers: frozenset[str] = frozenset(),
) -> list[Finding]:
    """The findings of one file, whose lines are `lines` (comments of `style`). PIN and REF are
    checked in source files when `table` is given; `functions` are the functions and types LOSAT
    defines, `identifiers` every name of LOSAT's code (`rust_names`)."""
    findings: list[Finding] = []
    texts = comment_texts(style, lines)
    if style != "text":
        findings.extend(inline_statement_findings(rel, texts, lines, style, functions, identifiers))
    last_citation: tuple[int, str] | None = None
    # The citations among the comment lines above the current line, back to a code line, an
    # empty comment line or a fence (the run): the last labelled one decides; without one, every
    # citation line of the run.
    run_label: tuple[int, str] | None = None
    run_citations: list[tuple[int, str]] = []
    in_fence = False
    # The last label line of the comment that cites source, with its indentation: a code line
    # indented 4 more under it is quoted, blank lines between or not.
    block_label: tuple[int, str, int] | None = None
    # A quoted call that goes on to the next line (`f(a,`): its lines up to the `;` are quoted.
    open_call: tuple[str, str] | None = None
    source_file = style != "text"
    if table is not None and dirs is None:
        dirs = directories(table)
    for number, text in enumerate(texts, start=1):
        if text is None:
            trailing = TRAILING_COMMENT_RE.search(lines[number - 1]) if style == "slash" else None
            if trailing is not None and table is not None:
                findings.extend(
                    citation_findings(rel, number, trailing.group("text"), table, dirs or frozenset())
                )
            in_fence = False  # a fence ends with its comment run
            run_label, run_citations = None, []
            block_label = None
            open_call = None
            continue
        cited = [run_label] if run_label is not None else list(run_citations)
        if FENCE_RE.match(text):
            if not in_fence:
                language = text.strip()[3:].strip().split(" ")[0].split(",")[0].lower()
                in_run = [entry for entry in cited if cites_source(entry[1])]
                # Without a citation in the run, the last citation above in the file (in source
                # files: a quote inside a function cited above the function).
                if not cited and last_citation is not None and source_file:
                    cited = [last_citation]
                found = cited_source(cited)
                if found is not None:
                    findings.append(Finding(rel, number, "QUOTE", found[0], found[1]))
                elif language in SOURCE_LANGUAGES:
                    where = f"cited at :{cited[0][0]}: {cited[0][1]}" if cited else "no citation"
                    findings.append(Finding(rel, number, "QUOTE", "source text", where))
                elif in_run and source_file and language not in RUST_FENCES:
                    findings.append(
                        Finding(rel, number, "QUOTE", "source text", f"cited at :{in_run[0][0]}: {in_run[0][1]}")
                    )
            in_fence = not in_fence
            run_label, run_citations = None, []
            continue
        if in_fence:
            continue
        if not text:
            run_label, run_citations = None, []
            continue
        citation = is_citation(text)
        if open_call is not None:
            if text and not citation and not text.endswith(".") and not PROSE_RE.search(text):
                findings.append(Finding(rel, number, "QUOTE", open_call[0], open_call[1]))
                if not text.endswith(","):
                    open_call = None
                continue
            open_call = None
        label = run_label
        if label is None and block_label is not None and text_indent(lines[number - 1]) >= block_label[2] + 4:
            label = (block_label[0], block_label[1])
        if label is not None and not citation and CODE_LINE_RE.search(text):
            found = cited_source([label])
            if found is None and cites_source(label[1]):
                found = ("source text", f"cited at :{label[0]}: {label[1]}")
            if found is not None:
                findings.append(Finding(rel, number, "QUOTE", found[0], found[1]))
                if text.endswith(",") and text.count("(") > text.count(")"):
                    open_call = found
        if source_file and code_after_citation(text) is not None:
            findings.append(Finding(rel, number, "QUOTE", CODE_AFTER_CITATION, text))
        if is_label_line(text):
            run_label = (number, text)
            if cites_source(text) and style != "text":
                block_label = (number, text, text_indent(lines[number - 1]))
        elif citation:
            run_citations.append((number, text))
        labelled = LABELLED_RE.match(text)
        if (
            NCBI_CITATION_RE.match(text)
            or NCBI_PATH_LABEL_RE.match(text)
            or (labelled is not None and NCBI_LABEL_RE.match(labelled.group("label")))
        ):
            source = foreign_source(text)
            if source is not None:
                findings.append(Finding(rel, number, "LABEL", source, text))
        if source_file and table is not None:
            findings.extend(citation_findings(rel, number, text, table, dirs or frozenset()))
        if citation:
            last_citation = (number, text)
    if style != "text":
        already_quoted = frozenset(f.line for f in findings if f.kind == "QUOTE")
        findings.extend(
            unbackticked_findings(rel, texts, lines, style, functions, identifiers, already_quoted)
        )
    return findings


def find(root: Path, table: dict[str, int] | None = None) -> list[Finding]:
    findings: list[Finding] = []
    if table is None:
        table = load_table()
    dirs = directories(table)
    functions, identifiers = rust_names(root)
    for path in source_files(root):
        text = path.read_text(encoding="utf-8", errors="replace")
        style = COMMENT_STYLES[path.suffix]
        rel = path.relative_to(root).as_posix()
        findings.extend(check_text(rel, text.splitlines(), style, table, dirs, functions, identifiers))
    return findings


# The unpinned forms of a citation of NCBI's tree that `--fix` rewrites to the pinned label.
CHECKOUT_PREFIX = r"(?:/mnt/c/Users/[^/\s]+/GitHub/|[\w./-]*/)?(?:ncbi-blast[\w.+-]*/)?"
PIN_FIXES: tuple[tuple[re.Pattern[str], str], ...] = (
    (
        re.compile(r"\bNCBI reference(?: \((?!598d8ae6\)|RpsbProc)[^)]*\))?:\s*" + CHECKOUT_PREFIX + r"(?=c\+\+/)"),
        "NCBI reference (598d8ae6): ",
    ),
    (re.compile(r"(?<![\w`])References?:\s*" + CHECKOUT_PREFIX + r"(?=c\+\+/)"), "NCBI reference (598d8ae6): "),
    (re.compile(r"(?<![\w`])(?:/mnt/c/Users/[^/\s]+/GitHub/)?ncbi-blast[\w.+-]*/(?=c\+\+/)"), ""),
)
# A label that names a file of NCBI's tree by its base name; `--fix` writes the file's path when
# the tree has one file of that name.
BARE_LABEL_RE = re.compile(
    r"\bNCBI (?P<kind>reference|unit test)(?P<s>s?)(?: \((?!598d8ae6\)|RpsbProc)[^)]*\))?:\s*"
    r"(?P<name>[\w+-]+\.(?:c|h|cc|cpp|cxx|hpp|inl|asn))\b"
)
BLANK_COMMENT_RE = re.compile(r"^\s*(?://[/!]?|#)\s*$")


def base_names(table: dict[str, int]) -> dict[str, str]:
    """The base names that one file of NCBI's tree has, with that file's path."""
    seen: dict[str, list[str]] = {}
    for path in table:
        if path.startswith("c++/"):
            seen.setdefault(path.rsplit("/", 1)[-1], []).append(path)
    return {name: paths[0] for name, paths in seen.items() if len(paths) == 1}


def fix_lines(
    lines: list[str],
    style: str,
    findings: Sequence[Finding],
    names: dict[str, str] | None = None,
) -> list[str]:
    """`lines` without the source text that `findings` name (QUOTE: a fence from its opening to
    its closing line, or one code line) and with PIN's labels pinned where the label names a
    path of NCBI's tree, or a base name that one file of the tree has (`names`, from
    `base_names`). A blank comment line left doubled, or left between a comment and code, goes
    too. Only comment lines change. REF is left to the writer: cite the right file and lines."""
    texts = comment_texts(style, lines)
    drop: set[int] = set()
    for finding in findings:
        index = finding.line - 1
        if finding.kind == "QUOTE":
            text = texts[index] or ""
            if finding.source in (INLINE_STATEMENT, INLINE_EXPRESSION, UNBACKTICKED_QUOTE):
                continue  # inside prose or unbackticked comment: the writer says it in words
            if finding.source == CODE_AFTER_CITATION:
                cut = code_after_citation(lines[index])
                if cut is not None:
                    lines[index] = lines[index][:cut].rstrip()
                continue
            if FENCE_RE.match(text):
                end = index + 1
                while end < len(lines) and texts[end] is not None and not FENCE_RE.match(texts[end] or ""):
                    end += 1
                closed = end < len(lines) and texts[end] is not None
                if closed or style != "text":
                    # A fence in a comment ends with the comment; in a Markdown file an
                    # unclosed one is left to the writer.
                    drop.update(range(index, end + 1 if closed else end))
            else:
                drop.add(index)
        elif finding.kind == "PIN":
            line = lines[index]
            for pattern, replacement in PIN_FIXES:
                line = pattern.sub(replacement, line)
            if names:
                line = BARE_LABEL_RE.sub(
                    lambda m: (
                        f"NCBI {m['kind']}{m['s']} (598d8ae6): {names[m['name']]}"
                        if m["name"] in names
                        else m[0]
                    ),
                    line,
                )
            lines[index] = line
    drop = {index for index in drop if texts[index] is not None}
    if style == "python":
        # A docstring's quotes stay: the text before closing quotes goes, and a fence that
        # shares the opening quotes' line loses only its backticks.
        for index in [i for i in drop if '"""' in lines[i] or "'''" in lines[i]]:
            drop.discard(index)
            closing = re.match(r"^(\s*)\S.*?(\"\"\"|''')\s*$", lines[index])
            if closing is not None and not lines[index].lstrip().startswith(closing.group(2)):
                lines[index] = closing.group(1) + closing.group(2)
            else:
                lines[index] = lines[index].replace("```", "", 1)

    def blank(line: str) -> bool:
        if style in ("text", "python"):
            return not line.strip() or BLANK_COMMENT_RE.match(line) is not None
        return BLANK_COMMENT_RE.match(line) is not None

    out: list[str] = []
    index = 0
    while index < len(lines):
        if index not in drop:
            out.append(lines[index])
            index += 1
            continue
        end = index
        while end in drop:
            end += 1
        before = out[-1] if out else None
        after = lines[end] if end < len(lines) else None
        blank_before = before is not None and blank(before)
        blank_after = after is not None and blank(after)
        # In Markdown and docstrings a blank line is a paragraph break: only a doubled one goes.
        comment_after = style in ("text", "python") or (
            after is not None and re.match(r"^\s*(?://|#)", after) is not None
        )
        if blank_before and blank_after:
            end += 1
        elif blank_before and not comment_after:
            out.pop()
        index = end
    return out


def rust_code(text: str) -> str:
    """The tokens of Rust source without its comments, one space between them where there was
    space: what a change of comments only leaves the same."""
    out: list[str] = []
    i, n = 0, len(text)
    while i < n:
        if text.startswith("//", i):
            j = text.find("\n", i)
            i = n if j < 0 else j
            continue
        if text.startswith("/*", i):
            depth = 0
            while i < n:
                if text.startswith("/*", i):
                    depth, i = depth + 1, i + 2
                elif text.startswith("*/", i):
                    depth, i = depth - 1, i + 2
                else:
                    i += 1
                if depth == 0:
                    break
            out.append(" ")
            continue
        raw = re.match(r'b?r(#*)"', text[i : i + 260])
        if raw and (i == 0 or not (text[i - 1].isalnum() or text[i - 1] == "_")):
            close = '"' + raw.group(1)
            j = text.find(close, i + raw.end())
            j = n if j < 0 else j + len(close)
            out.append(text[i:j])
            i = j
            continue
        if text[i] == '"':
            j = i + 1
            while j < n and text[j] != '"':
                j += 2 if text[j] == "\\" else 1
            out.append(text[i : j + 1])
            i = j + 1
            continue
        char = re.match(r"'(?:\\.[^']*|[^'\\])'", text[i : i + 16])
        if char:
            out.append(char.group(0))
            i += char.end()
            continue
        out.append(" " if text[i].isspace() else text[i])
        i += 1
    return " ".join("".join(out).split())


def python_code(text: str) -> str | None:
    """The AST of Python source without docstrings, or None when it does not parse."""
    try:
        tree = ast.parse(text)
    except (SyntaxError, ValueError):
        return None
    for node in ast.walk(tree):
        body = getattr(node, "body", None)
        if (
            isinstance(body, list)
            and body
            and isinstance(body[0], ast.Expr)
            and isinstance(body[0].value, ast.Constant)
            and isinstance(body[0].value.value, str)
        ):
            body[0] = ast.Pass()
    return ast.dump(tree)


def same_code(suffix: str, old: str, new: str) -> bool:
    """A fix of a Rust or Python file left its code as it was (strings included)."""
    if suffix == ".rs":
        return rust_code(old) == rust_code(new)
    if suffix == ".py":
        return python_code(old) == python_code(new)
    return True


def fix(root: Path) -> tuple[int, list[str]]:
    """Rewrite the files with findings (`fix_lines`): the number of files changed, and the files
    left as they were because the fix would have changed their code (a comment marker inside a
    string), which the writer fixes by hand. Line endings are kept."""
    table = load_table()
    dirs = directories(table)
    names = base_names(table)
    changed = 0
    refused: list[str] = []
    for path in source_files(root):
        try:
            with path.open(encoding="utf-8", newline="") as handle:
                text = handle.read()
        except UnicodeDecodeError:
            continue
        style = COMMENT_STYLES[path.suffix]
        lines = text.split("\n")
        rel = path.relative_to(root).as_posix()
        # Removing a line can bring the next under a label: fix until nothing changes.
        new = text
        for _ in range(5):
            findings = [
                f for f in check_text(rel, lines, style, table, dirs) if f.kind in ("QUOTE", "PIN")
            ]
            fixed = fix_lines(list(lines), style, findings, names) if findings else lines
            if fixed == lines:
                break
            lines = fixed
            new = "\n".join(lines)
        if new == text:
            continue
        if not same_code(path.suffix, text, new):
            refused.append(rel)
            continue
        with path.open("w", encoding="utf-8", newline="") as handle:
            handle.write(new)
        changed += 1
    return changed, refused


def write_table(ncbi_src: Path, rpsbproc_src: Path, out: Path = TABLE_PATH) -> int:
    """Write `ncbi_source_lines.tsv` from a checkout of NCBI's tree at 598d8ae6 (its git files)
    and the directory that holds rpsbproc 0.5.1's `RpsbProc/`."""
    import subprocess

    head = subprocess.run(
        ["git", "-C", str(ncbi_src), "rev-parse", "HEAD"], capture_output=True, text=True, check=True
    ).stdout.strip()
    if not head.startswith("598d8ae6"):
        raise SystemExit(f"{ncbi_src} is at {head}, not 598d8ae6")
    names = subprocess.run(
        ["git", "-C", str(ncbi_src), "ls-files"], capture_output=True, text=True, check=True
    ).stdout.split("\n")
    files = [(name, ncbi_src / name) for name in names if name]
    for directory, subdirs, entries in os.walk(rpsbproc_src / "RpsbProc"):
        subdirs[:] = [d for d in subdirs if d != ".svn"]
        for entry in entries:
            path = Path(directory) / entry
            files.append((path.relative_to(rpsbproc_src).as_posix(), path))
    rows = [
        "# The files of the trees that LOSAT's comments cite, with their line counts, read by",
        "# check_foreign_quotes.py (REF): c++/ = NCBI BLAST+ 2.17.0 (satoshikawato/ncbi-blast",
        "# 598d8ae6a72b923127ba2fbfaffd48e4c83bfbf4), RpsbProc/ = rpsbproc 0.5.1 (RpsbProc-src.tar.gz",
        "# of 2025-05-13, ftp.ncbi.nlm.nih.gov/pub/mmdb/cdd/rpsbproc/). A line that ends with / is a",
        "# directory; the lines after it are <name><TAB><lines> of its files. Written by",
        "# `check_foreign_quotes.py --write-table`.",
    ]
    directory = None
    for name, path in sorted(files):
        parent, _, base = name.rpartition("/")
        if parent + "/" != directory:
            directory = parent + "/"
            rows.append(directory)
        data = path.read_bytes()
        count = data.count(b"\n") + (0 if not data or data.endswith(b"\n") else 1)
        rows.append(f"{base}\t{count}")
    out.write_text("\n".join(rows) + "\n", encoding="utf-8")
    return len(files)


def parse_args(argv: Sequence[str] | None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--root", type=Path, default=Path("."), help="repository root")
    parser.add_argument("--write-table", action="store_true", help="regenerate the tree table")
    parser.add_argument(
        "--fix", action="store_true", help="remove the quoted source text and pin the labels"
    )
    parser.add_argument("--ncbi-src", type=Path, help="checkout of NCBI's tree at 598d8ae6")
    parser.add_argument("--rpsbproc-src", type=Path, help="directory that holds RpsbProc/")
    return parser.parse_args(argv)


ADVICE = {
    "QUOTE": "source text in a comment: keep or add the citation (path:lines), remove the text "
    "(a statement or expression in or outside backticks: name the function or member, say in "
    "words what it does)",
    "LABEL": "labels another project's code as NCBI's: name the project and version",
    "PIN": "write `NCBI reference (598d8ae6): c++/<path>:<lines>` (rpsbproc: "
    "`NCBI reference (RpsbProc 0.5.1, 2025-05-13): RpsbProc/<path>:<lines>`)",
    "REF": "cite a file and lines of the pinned tree (scripts/ncbi_ref.py prints them)",
}


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    if args.write_table:
        if args.ncbi_src is None or args.rpsbproc_src is None:
            raise SystemExit("--write-table needs --ncbi-src and --rpsbproc-src")
        count = write_table(args.ncbi_src, args.rpsbproc_src)
        print(f"wrote {TABLE_PATH.name}: {count} files")
        return 0
    root = args.root.resolve()
    if args.fix:
        changed, refused = fix(root)
        print(f"fixed {changed} files")
        for rel in refused:
            print(f"{rel}: not fixed, the fix would change code; fix it by hand")
    findings = find(root)
    for f in findings:
        print(f"{f.kind}\t{f.path}:{f.line}\t{f.source}\t{ADVICE[f.kind]}\t{f.text}")
    print(f"source quotes\tfindings={len(findings)}")
    return 0 if not findings else 1


if __name__ == "__main__":
    raise SystemExit(main())
