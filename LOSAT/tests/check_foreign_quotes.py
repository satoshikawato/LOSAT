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

The tables below are the one place that says which sources are NCBI's. There is no
allowlist of exceptions: cite the file and lines instead of quoting them.

Usage: python3 LOSAT/tests/check_foreign_quotes.py --root .
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
from typing import Iterable, Sequence

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
    "NCBI reference (RpsbProc 0.5.1, 2025-05-13)",
    "NCBI unit test (598d8ae6)",
)
# `NCBI reference` used as a label: followed by `:` or by a parenthesis.
NCBI_LABEL_USE_RE = re.compile(r"\bNCBI (?:reference|unit test)s?(?=\s*[:(])(?:\s*\([^)]*\))?")
# A path of a local checkout of NCBI's tree.
CHECKOUT_PATH_RE = re.compile(r"(?:ncbi-blast[\w.+-]*|GitHub)/c\+\+/")
# A path of the pinned trees, with its lines: `c++/src/x.c:12-20,31`.
TREE_PATH_RE = re.compile(
    r"(?<![\w.+/-])(?P<path>(?:c\+\+|RpsbProc)/[\w.+/-]*[\w+/])"
    r"(?::(?P<ranges>\d+(?:-\d+)?(?:,\s?\d+(?:-\d+)?)*))?"
)
TABLE_PATH = Path(__file__).with_name("ncbi_source_lines.tsv")

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
# A line labelled as NCBI's by a path: `NCBI c++/src/...:<lines>`, `NCBI 2.17.0 boost_erf.c:...`.
NCBI_PATH_LABEL_RE = re.compile(
    r"^NCBI(?:\s+[^\s:]+){0,4}?\s+(?P<path>[\w.+/-]+\.(?:c|cc|cpp|cxx|h|hh|hpp|hxx|inl|tcc))\b"
)
# A comment line that is a line of C or C++: a statement, a brace, a definition, a `#define`.
CODE_LINE_RE = re.compile(
    r"(?:=|\breturn\b|\w\()[^.`]*;\s*$|^\s*\}\s*$|\)\s*\{\s*$|^\s*#\s*define\b"
    r"|^\s*(?:static\s+)?(?:const\s+)?(?:unsigned\s+|signed\s+)?"
    r"(?:int|double|float|char|void|size_t|long|bool|Int4|Int8|Uint4|Uint8)\s+\**\w+\s*[(=;]"
)
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
    out = set()
    for path in table:
        parts = path.split("/")
        for end in range(1, len(parts)):
            out.add("/".join(parts[:end]) + "/")
    return frozenset(out)


def ref_problem(path: str, ranges: str | None, table: dict[str, int], dirs: frozenset[str]) -> str | None:
    """Why a cited path of the pinned trees is wrong, or None."""
    if path.endswith("/"):
        return None if path in dirs else "no such directory"
    if path not in table:
        return None if path + "/" in dirs else "no such file"
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


def check_text(
    rel: str,
    lines: Sequence[str],
    style: str = "slash",
    table: dict[str, int] | None = None,
    dirs: frozenset[str] | None = None,
) -> list[Finding]:
    """The findings of one file, whose lines are `lines` (comments of `style`). PIN and REF are
    checked in source files when `table` is given."""
    findings: list[Finding] = []
    last_citation: tuple[int, str] | None = None
    # The citations among the comment lines above the current line, back to a code line, an
    # empty comment line or a fence (the run): the last labelled one decides; without one, every
    # citation line of the run.
    run_label: tuple[int, str] | None = None
    run_citations: list[tuple[int, str]] = []
    in_fence = False
    source_file = style != "text"
    if table is not None and dirs is None:
        dirs = directories(table)
    for number, text in enumerate(comment_texts(style, lines), start=1):
        if text is None:
            in_fence = False  # a fence ends with its comment run
            run_label, run_citations = None, []
            continue
        cited = [run_label] if run_label is not None else list(run_citations)
        if FENCE_RE.match(text):
            if not in_fence:
                language = text.strip()[3:].strip().split(" ")[0].split(",")[0].lower()
                in_run = bool(cited)
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
                        Finding(rel, number, "QUOTE", "source text", f"cited at :{cited[0][0]}: {cited[0][1]}")
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
        if run_label is not None and not citation and CODE_LINE_RE.search(text):
            found = cited_source([run_label])
            if found is None:
                found = ("source text", f"cited at :{run_label[0]}: {run_label[1]}")
            findings.append(Finding(rel, number, "QUOTE", found[0], found[1]))
        if is_label_line(text):
            run_label = (number, text)
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
            problem = pin_problem(text)
            if problem is not None:
                findings.append(Finding(rel, number, "PIN", problem, text))
            for match in TREE_PATH_RE.finditer(text):
                problem = ref_problem(match.group("path"), match.group("ranges"), table, dirs or frozenset())
                if problem is not None:
                    findings.append(Finding(rel, number, "REF", problem, match.group(0)))
        if citation:
            last_citation = (number, text)
    return findings


def find(root: Path, table: dict[str, int] | None = None) -> list[Finding]:
    findings: list[Finding] = []
    if table is None:
        table = load_table()
    dirs = directories(table)
    for path in source_files(root):
        text = path.read_text(encoding="utf-8", errors="replace")
        style = COMMENT_STYLES[path.suffix]
        rel = path.relative_to(root).as_posix()
        findings.extend(check_text(rel, text.splitlines(), style, table, dirs))
    return findings


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
    parser.add_argument("--ncbi-src", type=Path, help="checkout of NCBI's tree at 598d8ae6")
    parser.add_argument("--rpsbproc-src", type=Path, help="directory that holds RpsbProc/")
    return parser.parse_args(argv)


ADVICE = {
    "QUOTE": "source text in a comment: keep or add the citation (path:lines), remove the text",
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
    findings = find(root)
    for f in findings:
        print(f"{f.kind}\t{f.path}:{f.line}\t{f.source}\t{ADVICE[f.kind]}\t{f.text}")
    print(f"source quotes\tfindings={len(findings)}")
    return 0 if not findings else 1


if __name__ == "__main__":
    raise SystemExit(main())
