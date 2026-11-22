#!/usr/bin/env python3
"""Fail when a Product Decision gives one version or one exception number twice.

No NCBI counterpart: this checks LOSAT's own decision records in docs/product_decisions/.

Several work lines add entries to the same record on the same day, and each takes the next
version and the next exception number by hand. On 2026-11-21, two pull requests into dev both
wrote `PD-LOSAT-NCBI-DEFECTS` 1.68 = "approved exception 20" (#1451, #1460), after two earlier
races that day (exception 18 and 1.63, then 1.66 to 1.68). Each time Git stopped the second pull
request with a conflict in the version list. A conflict resolved by keeping both entries, or an
entry added where it does not conflict, would leave two entries with one number. This check
fails then.

For each `PD-*.md` in docs/product_decisions/, it reads the header's `- Date:` field (the
version list) and its `- Version:` field, and fails when:
- a version appears twice in the version list;
- the versions of the list do not increase from first to last;
- the last version of the list is not the `- Version:` field's;
- the parentheses of the version list do not balance (the entries after an unbalanced one
  cannot be told from the text inside an entry);
- one exception number is introduced twice. An entry introduces an exception as
  `exception N,` (`approved exception 20, TBLASTX's linking ...`). Other mentions
  (`exception 6 also for ...`, `exception 8: ...`, `exceptions 3 and 5`) refer to it.

A version entry is a number such as `1.68` followed by a date (`1.68 2026-11-21 (...)`,
`1.1, 2026-11-01 (...)`) or by `the same day`, outside every parenthesis of the field.
Versions named inside an entry's parentheses ("numbered after 1.67", "first 1.54, then
1.56") are not entries. A record whose Date field has no entry (only `- Version:`) has
nothing to compare.

Usage: check_pd_numbers.py [--root DIR]
Exit status 0 when every record passes, 1 otherwise (one `::error` line per finding).
"""

from __future__ import annotations

import argparse
import re
import sys
from dataclasses import dataclass
from pathlib import Path

PD_DIR = Path("docs") / "product_decisions"

FIELD_RE = re.compile(r"^- ([A-Z][A-Za-z ]*):")
VERSION_FIELD_RE = re.compile(r"^- Version:\s*(\d+(?:\.\d+)+)")
ENTRY_RE = re.compile(r"(?<![\w.])(\d+\.\d+),?\s+(?:\d{4}-\d{2}-\d{2}\b|the same day\b)")
EXCEPTION_RE = re.compile(r"\bexception (\d+),")


@dataclass(frozen=True)
class Finding:
    path: Path
    line: int
    message: str


def version_key(text: str) -> tuple[int, ...]:
    return tuple(int(part) for part in text.split("."))


def header_field(lines: list[str], name: str) -> tuple[int, list[str]] | None:
    """The 1-based first line and the lines of the header field `- <name>:`."""
    for index, line in enumerate(lines):
        match = FIELD_RE.match(line)
        if match is None or match.group(1) != name:
            continue
        field = [line]
        for following in lines[index + 1 :]:
            if not following.startswith("  "):
                break
            field.append(following)
        return index + 1, field
    return None


def depth_zero_entries(field: list[str], first_line: int) -> tuple[list[tuple[str, int]], int | None]:
    """Version entries outside every parenthesis, with their 1-based line numbers, and the line
    of the first unbalanced parenthesis (None when they balance).

    An entry may break between its number and its date (`1.3,` at the end of one line,
    `2026-11-01 (` at the start of the next). An unbalanced parenthesis hides the entries after
    it, so it is a finding of its own.
    """
    text = "\n".join(field)
    entries: list[tuple[str, int]] = []
    depth = 0
    unbalanced: int | None = None
    opened_at: list[int] = []
    position = 0
    while position < len(text):
        char = text[position]
        if char == "(":
            depth += 1
            opened_at.append(position)
        elif char == ")":
            if depth == 0:
                if unbalanced is None:
                    unbalanced = first_line + text.count("\n", 0, position)
            else:
                depth -= 1
                opened_at.pop()
        elif depth == 0 and (position == 0 or not text[position - 1].isalnum()):
            match = ENTRY_RE.match(text, position)
            if match is not None:
                entries.append((match.group(1), first_line + text.count("\n", 0, position)))
                position = match.end(1)
                continue
        position += 1
    if unbalanced is None and opened_at:
        unbalanced = first_line + text.count("\n", 0, opened_at[0])
    return entries, unbalanced


def check_record(path: Path, display: Path) -> list[Finding]:
    lines = path.read_text(encoding="utf-8").splitlines()
    findings: list[Finding] = []

    version_line = 0
    version = None
    for index, line in enumerate(lines):
        match = VERSION_FIELD_RE.match(line)
        if match is not None:
            version_line, version = index + 1, match.group(1)
            break

    date = header_field(lines, "Date")
    if date is None:
        return findings
    first_line, field = date

    entries, unbalanced = depth_zero_entries(field, first_line)
    if unbalanced is not None:
        findings.append(
            Finding(
                display,
                unbalanced,
                "unbalanced parenthesis in the version list; the entries after it cannot be read",
            )
        )
    seen: dict[str, int] = {}
    for number, line_no in entries:
        if number in seen:
            findings.append(
                Finding(
                    display,
                    line_no,
                    f"version {number} appears twice in the version list (also line {seen[number]}); "
                    "renumber the later entry after the last one on dev",
                )
            )
        else:
            seen[number] = line_no
    for (before, _), (after, line_no) in zip(entries, entries[1:]):
        if version_key(after) <= version_key(before) and after != before:
            findings.append(
                Finding(
                    display,
                    line_no,
                    f"version {after} follows {before} in the version list; versions must increase",
                )
            )
    if entries and version is not None and entries[-1][0] != version:
        findings.append(
            Finding(
                display,
                version_line,
                f"the Version field says {version} but the version list ends at {entries[-1][0]} "
                f"(line {entries[-1][1]})",
            )
        )

    introduced: dict[str, int] = {}
    for offset, line in enumerate(field):
        for match in EXCEPTION_RE.finditer(line):
            number = match.group(1)
            line_no = first_line + offset
            if number in introduced:
                findings.append(
                    Finding(
                        display,
                        line_no,
                        f"exception {number} is introduced twice in the version list (also line "
                        f"{introduced[number]}); renumber the later one after the last one on dev",
                    )
                )
            else:
                introduced[number] = line_no
    return findings


def check_root(root: Path) -> list[Finding]:
    findings: list[Finding] = []
    for path in sorted((root / PD_DIR).glob("PD-*.md")):
        findings.extend(check_record(path, path.relative_to(root)))
    return findings


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[2])
    args = parser.parse_args(argv)
    findings = check_root(args.root)
    for finding in findings:
        print(f"::error file={finding.path},line={finding.line}::{finding.message}")
    if findings:
        print(f"{len(findings)} numbering finding(s) in {PD_DIR}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
