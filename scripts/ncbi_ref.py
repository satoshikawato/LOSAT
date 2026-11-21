#!/usr/bin/env python3
"""Print the NCBI source lines that LOSAT's comments cite.

No NCBI counterpart: a reading aid for LOSAT's comments, not ported behavior.

LOSAT's comments cite NCBI's C/C++ source by path and lines in the pinned tree
(`NCBI reference (598d8ae6): c++/<path>:<a>-<b>`) and do not copy it (AGENTS.md rule 4). This
prints the cited lines from a checkout of that tree (`$NCBI_SRC`, or `--ncbi-src`); citations of
RpsbProc (`RpsbProc/<path>:<a>-<b>`, version 0.5.1) are read from the directory that holds its
`RpsbProc/` (`$RPSBPROC_SRC`, or `--rpsbproc-src`):

  scripts/ncbi_ref.py LOSAT/src/algorithm/blastn/hsp.rs:120   the citations nearest above line 120
  scripts/ncbi_ref.py LOSAT/src/algorithm/blastn/hsp.rs       every citation of the file
  scripts/ncbi_ref.py c++/src/algo/blast/core/blast_hits.c:993-1001

With a line, the citations are those of the comment that holds the line, else of the nearest
comment above it that cites NCBI's source. `--context N` adds N lines around each range.
"""

from __future__ import annotations

import argparse
import os
import re
import subprocess
import sys
from pathlib import Path
from typing import Iterable, Sequence

PINNED = "598d8ae6a72b923127ba2fbfaffd48e4c83bfbf4"
# A path of NCBI's tree (or of RpsbProc) with its lines: `c++/<path>:12-20,31`.
CITE_RE = re.compile(
    r"(?<![\w.+/-])(?P<path>(?:c\+\+|RpsbProc)/[\w.+/-]*[\w+])"
    r"(?::(?P<ranges>\d+(?:-\d+)?(?:\s*,\s*\d+(?:-\d+)?)*))?"
)
COMMENT_RE = re.compile(r"^\s*(?://|#)")


def ranges_of(text: str | None) -> list[tuple[int, int]]:
    out: list[tuple[int, int]] = []
    for part in (text or "").split(","):
        part = part.strip()
        if part:
            first, _, last = part.partition("-")
            out.append((int(first), int(last or first)))
    return out


def citations(text: str) -> list[tuple[str, list[tuple[int, int]]]]:
    return [(m.group("path"), ranges_of(m.group("ranges"))) for m in CITE_RE.finditer(text)]


def nearest(lines: Sequence[str], line: int) -> list[str]:
    """The comment lines with citations of the comment that holds `line` (1-based) or, when it
    has none, of the nearest comment above it that has some."""
    index = min(max(line, 1), len(lines)) - 1
    while index >= 0:
        if COMMENT_RE.match(lines[index]):
            start = index
            while start > 0 and COMMENT_RE.match(lines[start - 1]):
                start -= 1
            end = index
            while end + 1 < len(lines) and COMMENT_RE.match(lines[end + 1]):
                end += 1
            cited = [text for text in lines[start : end + 1] if CITE_RE.search(text)]
            if cited:
                return cited
            index = start - 1
        else:
            index -= 1
    return []


def show(ncbi: Path, path: str, ranges: list[tuple[int, int]], context: int, out=None) -> bool:
    out = out or sys.stdout
    source = ncbi / path
    if not source.is_file():
        print(f"{path}: not in {ncbi}", file=out)
        return False
    text = source.read_text(encoding="utf-8", errors="replace").split("\n")
    if text and text[-1] == "":
        text.pop()
    if not ranges:
        print(f"== {path} ({len(text)} lines; the citation names no lines)", file=out)
        return True
    ok = True
    for first, last in ranges:
        span = f"{first}" if first == last else f"{first}-{last}"
        if not 1 <= first <= last <= len(text):
            print(f"== {path}:{span}: outside the file's {len(text)} lines", file=out)
            ok = False
            continue
        print(f"== {path}:{span}", file=out)
        for number in range(max(1, first - context), min(len(text), last + context) + 1):
            mark = " " if first <= number <= last else "."
            print(f"{number:6}{mark} {text[number - 1]}", file=out)
    return ok


def resolve_ncbi(option: str | None) -> Path | None:
    value = option or os.environ.get("NCBI_SRC")
    if not value:
        return None
    root = Path(value)
    head = subprocess.run(
        ["git", "-C", str(root), "rev-parse", "HEAD"], capture_output=True, text=True
    ).stdout.strip()
    if head and head != PINNED:
        print(f"ncbi_ref.py: warning: {root} is at {head[:12]}, not {PINNED[:8]}", file=sys.stderr)
    return root


def targets(arg: str) -> Iterable[str]:
    if arg.startswith(("c++/", "RpsbProc/")):
        yield arg
        return
    file, _, line = arg.rpartition(":")
    if not file or not line.isdigit():
        file, line = arg, ""
    lines = Path(file).read_text(encoding="utf-8", errors="replace").split("\n")
    if line:
        yield from nearest(lines, int(line))
    else:
        yield from (text for text in lines if COMMENT_RE.match(text) and CITE_RE.search(text))


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("target", help="<file>[:<line>] of LOSAT, or c++/<path>:<lines>")
    parser.add_argument("--ncbi-src", help="checkout of NCBI's tree (default: $NCBI_SRC)")
    parser.add_argument(
        "--rpsbproc-src", help="directory that holds RpsbProc/ (default: $RPSBPROC_SRC)"
    )
    parser.add_argument("--context", type=int, default=0, help="lines around each range")
    args = parser.parse_args(argv)
    ncbi = resolve_ncbi(args.ncbi_src)
    rpsbproc_value = args.rpsbproc_src or os.environ.get("RPSBPROC_SRC")
    rpsbproc = Path(rpsbproc_value) if rpsbproc_value else None
    ok, any_cited = True, False
    for text in targets(args.target):
        for path, ranges in citations(text):
            any_cited = True
            if path.startswith("RpsbProc/"):
                root, need = rpsbproc, "$RPSBPROC_SRC (or --rpsbproc-src) to the directory that holds RpsbProc/"
            else:
                root, need = ncbi, f"$NCBI_SRC (or --ncbi-src) to a checkout of NCBI's tree at {PINNED}"
            if root is None:
                print(f"ncbi_ref.py: {path}: set {need}", file=sys.stderr)
                ok = False
                continue
            ok = show(root, path, ranges, args.context) and ok
    if not any_cited:
        print(f"ncbi_ref.py: no citation of NCBI's tree at {args.target}", file=sys.stderr)
        return 1
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
