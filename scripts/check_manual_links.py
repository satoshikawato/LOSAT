#!/usr/bin/env python3
"""Checks the user manual's Markdown sources (docs/manual/src) before mdBook builds them.

- Every page that SUMMARY.md lists exists, and every page in the tree is listed.
- Every relative link and image points at a file in the tree, and a `#fragment` names a
  heading of that page (mdBook's heading ids: lower case, punctuation dropped, spaces as `-`).
- Text outside code says none of "verification", "verify", "verified", "certified",
  "certification": the Owner's rule for text users read (2026-11-04), held by this check.

Usage: check_manual_links.py SRC_DIR. Exit status 1 with one line per problem.
"""

from __future__ import annotations

import re
import sys
import unicodedata
from pathlib import Path

LINK = re.compile(r"!?\[(?:[^\]\[]|\[[^\]]*\])*\]\(([^)\s]+)(?:\s+\"[^\"]*\")?\)")
HTML_SRC = re.compile(r"""<img\s[^>]*src=["']([^"']+)["']""", re.IGNORECASE)
HEADING = re.compile(r"^(#{1,6})\s+(.*?)\s*#*\s*$")
FENCE = re.compile(r"^\s*(```|~~~)")
INLINE_CODE = re.compile(r"`+[^`]*`+")
FORBIDDEN = re.compile(r"\b(verif(?:y|ied|ies|ication|ying)|certif(?:ied|ication|y))\b", re.IGNORECASE)


def heading_id(text: str) -> str:
    """mdBook's `normalize_id`: keep letters, digits, `_` and `-`, spaces become `-`, lower case."""
    text = re.sub(r"<[^>]+>", "", text)
    text = INLINE_CODE.sub(lambda m: m.group(0).strip("`"), text)
    text = re.sub(r"\[([^\]]*)\]\([^)]*\)", r"\1", text)
    out = []
    for ch in text.strip():
        if ch.isalnum() or ch in "_-":
            out.append(ch.lower())
        elif ch.isspace():
            out.append("-")
    return unicodedata.normalize("NFC", "".join(out))


def prose_lines(text: str):
    """(line number, line) of the lines outside fenced code blocks."""
    in_fence = False
    for number, line in enumerate(text.splitlines(), 1):
        if FENCE.match(line):
            in_fence = not in_fence
            continue
        if not in_fence:
            yield number, line


def anchors(path: Path) -> set[str]:
    ids: dict[str, int] = {}
    found = set()
    for _, line in prose_lines(path.read_text(encoding="utf-8")):
        match = HEADING.match(line)
        if not match:
            continue
        base = heading_id(match.group(2))
        count = ids.get(base, 0)
        ids[base] = count + 1
        found.add(base if count == 0 else f"{base}-{count}")
    return found


def check(src: Path) -> list[str]:
    problems: list[str] = []
    pages = {p.relative_to(src).as_posix() for p in src.rglob("*.md")}
    summary = src / "SUMMARY.md"
    if not summary.is_file():
        return [f"{summary}: missing"]
    listed = set()
    for number, line in prose_lines(summary.read_text(encoding="utf-8")):
        for target in LINK.findall(line):
            listed.add(target)
            if target not in pages:
                problems.append(f"SUMMARY.md:{number}: lists {target}, which does not exist")
    for page in sorted(pages - listed - {"SUMMARY.md"}):
        problems.append(f"{page}: not listed in SUMMARY.md")

    anchor_cache: dict[Path, set[str]] = {}
    for page in sorted(pages - {"SUMMARY.md"}):
        path = src / page
        for number, line in prose_lines(path.read_text(encoding="utf-8")):
            prose = INLINE_CODE.sub("", line)
            for word in FORBIDDEN.findall(re.sub(r"\]\([^)]*\)", "]", prose)):
                problems.append(f"{page}:{number}: says '{word}' (use 'check', 'Supported' or 'Not yet supported')")
            for target in LINK.findall(prose) + HTML_SRC.findall(prose):
                if re.match(r"^[a-z][a-z0-9+.-]*:", target, re.IGNORECASE):
                    continue
                file_part, _, fragment = target.partition("#")
                dest = (path.parent / file_part).resolve() if file_part else path.resolve()
                try:
                    dest.relative_to(src.resolve())
                except ValueError:
                    problems.append(f"{page}:{number}: {target} leaves the manual (use a full https:// link)")
                    continue
                if not dest.exists():
                    problems.append(f"{page}:{number}: {target} does not exist")
                    continue
                if fragment and dest.suffix == ".md":
                    if dest not in anchor_cache:
                        anchor_cache[dest] = anchors(dest)
                    if fragment not in anchor_cache[dest]:
                        problems.append(f"{page}:{number}: {target} names no heading of {dest.name}")
    return problems


def main(argv: list[str]) -> int:
    if len(argv) != 2:
        print(__doc__.strip().splitlines()[-1], file=sys.stderr)
        return 2
    problems = check(Path(argv[1]))
    for problem in problems:
        print(problem)
    return 1 if problems else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
