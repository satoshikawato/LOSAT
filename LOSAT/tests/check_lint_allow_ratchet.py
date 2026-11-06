#!/usr/bin/env python3
"""Keep the files with a blanket `#![allow(warnings ...)]` from growing in number.

No NCBI counterpart: this checks LOSAT's own lint settings, not ported behavior.

`#![allow(warnings, clippy::all)]` at the top of `LOSAT/src/lib.rs` and `main.rs`
(2026-06-23) and of most integration-test files silences every rustc and clippy
warning in their crate or module, so `cargo clippy -- -D warnings` in CI cannot
fail there. New
test files copied the line with the rest of their header (SWEBOK audit
2026-11-04, KA4: 22 -> 23 test files). This ratchet fails when a Rust file that
is not in `lint_allow_warnings.tsv` carries such an attribute (also an outer
`#[allow(warnings)]`, or `clippy::all`), and when a listed
file no longer does (remove it from the list: a shorter list is welcome, and
keeps the attribute from coming back there).

Usage: python3 LOSAT/tests/check_lint_allow_ratchet.py --root .
"""

from __future__ import annotations

import argparse
import os
import re
import sys
from pathlib import Path
from typing import Iterable, Sequence

# An attribute whose `allow(...)` or `expect(...)` names `warnings` or `clippy::all`, inner or outer,
# also inside `cfg_attr(...)` and across lines (rustfmt breaks a long list):
# `#![allow(warnings, clippy::all)]`, `#![cfg_attr(test, allow(warnings))]`, `#[allow(warnings)] mod x;`.
ALLOW_WARNINGS_RE = re.compile(
    r"#!?\s*\[[^\]]*\b(?:allow|expect)\s*\([^)]*\b(?:warnings|clippy\s*::\s*all)\b"
)
COMMENT_RE = re.compile(r"/\*.*?\*/|//[^\n]*", re.S)
SKIPPED_DIRS = frozenset({".git", ".worktrees", "target", "node_modules"})
ALLOWLIST_HEADER = "path"


def rust_files(root: Path) -> Iterable[Path]:
    """Every `.rs` file under `root`, skipping build output, worktrees and package caches."""
    for directory, subdirs, files in os.walk(root):
        subdirs[:] = sorted(d for d in subdirs if d not in SKIPPED_DIRS)
        for name in sorted(files):
            if name.endswith(".rs"):
                yield Path(directory) / name


def has_allow_warnings(text: str) -> bool:
    """Whether `text`, without its comments, holds the attribute."""
    return ALLOW_WARNINGS_RE.search(COMMENT_RE.sub(" ", text)) is not None


def find_files(root: Path) -> list[str]:
    """The repository-relative paths (with `/`) of the files that carry the attribute."""
    found = []
    for path in rust_files(root):
        text = path.read_text(encoding="utf-8", errors="replace")
        if has_allow_warnings(text):
            found.append(path.relative_to(root).as_posix())
    return sorted(found)


def load_allowlist(path: Path) -> list[str]:
    lines = [line.strip() for line in path.read_text(encoding="utf-8").splitlines()]
    rows = [line for line in lines if line and not line.startswith("#")]
    if not rows or rows[0] != ALLOWLIST_HEADER:
        raise ValueError(f"{path}: the first row must be the header `{ALLOWLIST_HEADER}`")
    entries = rows[1:]
    if len(set(entries)) != len(entries):
        raise ValueError(f"{path}: a path is listed twice")
    return entries


def evaluate(found: Sequence[str], allowed: Sequence[str]) -> tuple[list[str], list[str]]:
    """(unexpected, stale): files with the attribute that are not listed, listed files without it."""
    unexpected = sorted(set(found) - set(allowed))
    stale = sorted(set(allowed) - set(found))
    return unexpected, stale


def parse_args(argv: Sequence[str] | None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--root", type=Path, default=Path("."), help="repository root")
    parser.add_argument(
        "--allowlist",
        type=Path,
        help="default: <root>/LOSAT/tests/lint_allow_warnings.tsv",
    )
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    root = args.root.resolve()
    allowlist_path = args.allowlist or root / "LOSAT" / "tests" / "lint_allow_warnings.tsv"
    try:
        allowed = load_allowlist(allowlist_path)
        found = find_files(root)
    except (OSError, ValueError) as error:
        print(f"lint-allow ratchet error: {error}", file=sys.stderr)
        return 2
    unexpected, stale = evaluate(found, allowed)
    for path in unexpected:
        print(
            f"UNEXPECTED\t{path}\tnew `#![allow(warnings ...)]`: remove it and allow only the lints "
            "the file needs, by name",
        )
    for path in stale:
        print(f"STALE\t{path}\tno longer has the attribute: remove it from {allowlist_path.name}")
    print(f"lint-allow ratchet\tfiles={len(found)}\tlisted={len(allowed)}\t"
          f"unexpected={len(unexpected)}\tstale={len(stale)}")
    return 0 if not unexpected and not stale else 1


if __name__ == "__main__":
    raise SystemExit(main())
