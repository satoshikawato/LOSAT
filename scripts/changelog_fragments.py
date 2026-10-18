#!/usr/bin/env python3
"""Pending CHANGELOG entries, one file per pull request (changelog.d/README.md).

Every pull request that used to add a line under `## [Unreleased]` in CHANGELOG.md inserted it at
the same place, so any two open pull requests conflicted there. Each one now adds its own file
`changelog.d/<branch, "/" as "-">.md` instead, and a release folds the files into CHANGELOG.md.

Usage:
  changelog_fragments.py check            every fragment is well formed (CI)
  changelog_fragments.py collect [--dry-run]
                                          move every fragment under `## [Unreleased]` in
                                          CHANGELOG.md, newest first, and delete the files
"""
from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
FRAGMENTS = ROOT / "changelog.d"
CHANGELOG = ROOT / "CHANGELOG.md"
HEADING = "## [Unreleased]"


def fragments() -> list[Path]:
    return sorted(p for p in FRAGMENTS.glob("*.md") if p.name != "README.md")


def problems(path: Path) -> list[str]:
    """A fragment is one or more entries: a line starting with "- ", continued by lines that start
    with two spaces; blank lines between entries are allowed; nothing else (no heading)."""
    try:
        lines = path.read_text(encoding="utf-8").splitlines()
    except UnicodeDecodeError:
        return ["not UTF-8"]
    found = []
    entries = 0
    for number, line in enumerate(lines, 1):
        if line.startswith("- "):
            entries += 1
        elif line.startswith("  ") and entries and line.strip():
            continue
        elif not line.strip():
            continue
        else:
            found.append(f"line {number}: not an entry ('- ...') or its continuation ('  ...')")
    if entries == 0:
        found.append("no entry ('- ...')")
    return found


def added_at(path: Path) -> int:
    """Commit time of the commit that added the file (0 when it is not committed: it goes first)."""
    out = subprocess.run(
        ["git", "log", "--diff-filter=A", "--format=%ct", "-1", "--", str(path.relative_to(ROOT))],
        cwd=ROOT, capture_output=True, text=True, check=False,
    ).stdout.strip()
    return int(out) if out else 0


def check() -> int:
    bad = 0
    for path in fragments():
        for problem in problems(path):
            print(f"{path.relative_to(ROOT)}: {problem}")
            bad += 1
    print(f"{len(fragments())} changelog fragments, {bad} problems")
    return 1 if bad else 0


def collect(dry_run: bool) -> int:
    if check():
        return 1
    files = fragments()
    if not files:
        print("nothing to collect")
        return 0
    # newest first, as the entries under [Unreleased] are ordered; uncommitted files first
    files.sort(key=lambda p: (-added_at(p) if added_at(p) else -(1 << 62), p.name))
    entries = []
    for path in files:
        entries.extend(line for line in path.read_text(encoding="utf-8").splitlines() if line.strip())
    text = CHANGELOG.read_text(encoding="utf-8")
    at = text.find(HEADING + "\n")
    if at < 0:
        print(f"{CHANGELOG.name}: no '{HEADING}' heading")
        return 1
    head = at + len(HEADING) + 1
    rest = text[head:].lstrip("\n")
    new = text[:head] + "\n" + "\n".join(entries) + "\n" + ("\n" if rest.startswith("## ") else "") + rest
    for path in files:
        print(f"{'would move' if dry_run else 'moved'} {path.relative_to(ROOT)}")
    if dry_run:
        return 0
    CHANGELOG.write_text(new, encoding="utf-8")
    for path in files:
        path.unlink()
    print(f"{len(files)} fragments, {len(entries)} lines -> {CHANGELOG.name}; commit both")
    return 0


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = parser.add_subparsers(dest="mode", required=True)
    sub.add_parser("check")
    c = sub.add_parser("collect")
    c.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()
    return check() if args.mode == "check" else collect(args.dry_run)


if __name__ == "__main__":
    sys.exit(main())
