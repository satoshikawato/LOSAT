#!/usr/bin/env python3
"""Fail when a comment quotes the source text of a project other than NCBI.

No NCBI counterpart: this checks LOSAT's own comments, not ported behavior.

Comments may quote only NCBI's own source, which carries NCBI's PUBLIC DOMAIN NOTICE (the
maintainer's decision of 2026-11-20). Code from any other project is cited by
`<project> <version>: <path>:<lines>` and not quoted, also when NCBI's tree carries a copy of
it (`LMDB 0.9.24 (NCBI 598d8ae6 copy): c++/src/util/lmdb/mdb.c:...`). The 2026-11-20 audit of
the published sources found quoted text of Microsoft's C runtime, libstdc++ and glibc (GPL and
LGPL code, some of it labelled "NCBI reference"), and of LMDB, Boost, lmdb++, SQLite and the
Rust crates bzip2 and ruzstd; those quotes were removed and their citations kept.

Two checks, over the comments of the repository's source files (`COMMENT_STYLES`: line
comments of Rust, C, C++, JavaScript and TypeScript; `#` comments of TOML, YAML, shell and
Python, and Python's docstrings):
- QUOTE: a fenced block (```` ``` ````) in a comment is a quote of the source that the
  comment lines right above it cite: the last of them that starts with a citation label
  (`<label>: <path>...`, `NCBI reference ...`, `NCBI <path>`, `<project> <version>
  reference ...`), else each of them that names a source; without one, the last citation
  above it in the file (a quote inside a function cited above the function). It fails when
  that citation is of another project: it names a source of `FOREIGN_SOURCES` or a path of
  `NON_NCBI_PATHS`, or its label is not one of `NCBI_LABEL_RE`. A comment line that is a
  line of C or C++ (`CODE_LINE_RE`) right under such a citation, in the same comment, is a
  quote without a fence and fails the same way.
- LABEL: a comment line labelled as NCBI's (`NCBI reference ...`, `NCBI <path>`) fails when
  its label or path names another project's source (the code is not NCBI's).
A fenced block with no citation above it in its file is not a quote and is not checked.
Markdown files are not read (their fences are commands and examples, not quotes).

The tables below are the one place that says which sources are NCBI's. There is no
allowlist of exceptions: cite the file and lines instead of quoting them.

Usage: python3 LOSAT/tests/check_foreign_quotes.py --root .
"""

from __future__ import annotations

import argparse
import os
import re
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Sequence

# Labels of NCBI's own source: `NCBI reference (598d8ae6): c++/...`, `NCBI unit test ...`,
# `NCBI 2.17.0: ...`, `Reference: ncbi-blast/c++/...`.
NCBI_LABEL_RE = re.compile(r"^(?:NCBI(?!['\u2019]s)\b|Reference\b)")
# A line that says it cites NCBI's source.
NCBI_CITATION_RE = re.compile(r"^(?:NCBI reference|NCBI unit test|Reference)\b")

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
}
FENCE_RE = re.compile(r"^\s*```")
# A cited source: a path with a source extension, a libstdc++ header path, or `<path>:<line>`.
CITE_RE = re.compile(
    r"[\w.+-]*(?:/[\w.+-]+)*\.(?:c|cc|cpp|cxx|h|hh|hpp|hxx|inl|tcc|asn|def|prt)\b"
    r"|\b(?:bits|ext|std)/\w+"
    r"|[\w.+-]+(?:/[\w.+-]+)+:\d"
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


def comment_texts(style: str, lines: Sequence[str]) -> list[str | None]:
    """The text of each line that is a comment (`style` of `COMMENT_STYLES`), None for the others.

    A Python docstring is a comment from its opening line to its closing one."""
    texts: list[str | None] = []
    docstring: str | None = None
    for line in lines:
        if docstring is not None:
            if docstring in line:
                docstring = None
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


def check_text(rel: str, lines: Sequence[str], style: str = "slash") -> list[Finding]:
    """The findings of one file, whose lines are `lines` (comments of `style`)."""
    findings: list[Finding] = []
    last_citation: tuple[int, str] | None = None
    # The comment lines above the current line, back to a code line, an empty comment line or
    # a fence: the citation of a fence that opens next is among them.
    run: list[tuple[int, str]] = []
    in_fence = False
    for number, text in enumerate(comment_texts(style, lines), start=1):
        if text is None:
            in_fence = False  # a fence ends with its comment run
            run = []
            continue
        # The last labelled citation of the run decides; without one, every citation line of
        # the run.
        labels = [(n, t) for n, t in run if is_label_line(t)]
        cited = labels[-1:] or [(n, t) for n, t in run if is_citation(t)]
        if FENCE_RE.match(text):
            if not in_fence:
                # Without a citation in the run, the last citation above in the file.
                if not cited and last_citation is not None:
                    cited = [last_citation]
                found = cited_source(cited)
                if found is not None:
                    findings.append(Finding(rel, number, "QUOTE", found[0], found[1]))
            in_fence = not in_fence
            run = []
            continue
        if in_fence:
            continue
        if not text:
            run = []
            continue
        if CODE_LINE_RE.search(text) and not is_citation(text):
            found = cited_source(cited)
            if found is not None:
                findings.append(Finding(rel, number, "QUOTE", found[0], found[1]))
        run.append((number, text))
        labelled = LABELLED_RE.match(text)
        if (
            NCBI_CITATION_RE.match(text)
            or NCBI_PATH_LABEL_RE.match(text)
            or (labelled is not None and NCBI_LABEL_RE.match(labelled.group("label")))
        ):
            source = foreign_source(text)
            if source is not None:
                findings.append(Finding(rel, number, "LABEL", source, text))
        if is_citation(text):
            last_citation = (number, text)
    return findings


def find(root: Path) -> list[Finding]:
    findings: list[Finding] = []
    for path in source_files(root):
        text = path.read_text(encoding="utf-8", errors="replace")
        style = COMMENT_STYLES[path.suffix]
        findings.extend(check_text(path.relative_to(root).as_posix(), text.splitlines(), style))
    return findings


def parse_args(argv: Sequence[str] | None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--root", type=Path, default=Path("."), help="repository root")
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    root = args.root.resolve()
    findings = find(root)
    for f in findings:
        if f.kind == "QUOTE":
            advice = "quotes another project's source: keep the citation, remove the quote"
        else:
            advice = "labels another project's code as NCBI's: name the project and version"
        print(f"{f.kind}\t{f.path}:{f.line}\t{f.source}\t{advice}\t{f.text}")
    print(f"foreign quotes\tfindings={len(findings)}")
    return 0 if not findings else 1


if __name__ == "__main__":
    raise SystemExit(main())
