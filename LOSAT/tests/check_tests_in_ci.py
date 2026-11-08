#!/usr/bin/env python3
"""Fail when a Python or Node test file of the repository is run by no workflow that runs on its own.

No NCBI counterpart: this checks LOSAT's own CI wiring, not ported behavior.

The SWEBOK audit of 2026-11-04 (KA5) found five `test_*.py` files with 35 tests
that no workflow ran (`test_audit_blastp_v010.py`, `test_audit_tblastx_v010.py`,
`test_benchmark_wasm_threading.py`, `test_certify_blastn_v010.py`,
`test_compare_blastn_parity.py`), so a change to the tools they test could break
them unseen. On 2026-11-08 the same audit found `LOSAT/tests/test_benchmark_wasi_reuse.js`
(the test of the Wasm benchmark runner) unrun; `wasm-threading.yml`, which the nightly
workflow calls, now runs it, and this check reads `test_*.js` files too. The workflows name
the test files they run one by one (`unittest discover -s LOSAT/tests -p 'test_x.py'`,
`unittest LOSAT/tests/test_x.py`, `node LOSAT/tests/test_x.js`), so this check fails when
the file name of a `test_*.py` or `test_*.js` file appears in the
`run:` commands of no workflow that runs without a person: one started by
push, pull_request(_target), merge_group or schedule, or one such a workflow
calls (`uses: ./.github/workflows/<file>`). A name only in a `paths:` filter, in
a comment, or in a workflow that runs only when dispatched (release readiness,
platform certification) does not count. A file name shared by two directories
counts as named for both.

Only files named `test_*.py` and `test_*.js` are test files here. web/app's vitest files
(`*.test.ts`) are run by web.yml as one suite and are out of scope, as is any other
`*.test.js` file.

`python/tests/` is pylosat's pytest suite, which `python.yml` runs as a directory
(`pytest tests`, working directory `python`): its files are run while that
workflow runs that command; without it, each file needs its name like any other.

The workflow files are read line by line (the Python of the runner has no YAML
parser), with these limits:
- `on:` is read inline (`on: push`), as a flow list, or as a block; a flow
  mapping (`on: {push: ...}`) is not, so its workflow does not count.
- A ` #` inside a quoted `run:` value is taken as the start of a comment.
- An action input named `run` under `with:` is read as a command, and a job's
  `if:` is not read: a job that runs only when dispatched, in a workflow that
  push or pull_request starts, counts.
The first two can only report a file as not run. The search walks the
directory, so a local untracked `test_*.py` file is reported too; CI checks a
clean checkout.

Usage: python3 LOSAT/tests/check_tests_in_ci.py --root .
"""

from __future__ import annotations

import argparse
import os
import re
import sys
from pathlib import Path
from typing import Iterable, Sequence

SKIPPED_DIRS = frozenset({".git", ".worktrees", "target", "node_modules", ".venv", "__pycache__"})
TEST_FILE_RE = re.compile(r"^test_.*\.(py|js)$")
# A YAML comment: `#` at the start of a line or after a space (`'#'` inside a value is not one).
YAML_COMMENT_RE = re.compile(r"(^|\s)#.*$")
# Events that start a workflow without a person.
AUTOMATIC_EVENTS = frozenset({"push", "pull_request", "pull_request_target", "merge_group", "schedule"})
ON_RE = re.compile(r"""^(?:on|"on"|'on'|true):\s*(.*)$""")
RUN_RE = re.compile(r"^(\s*(?:-\s+)?)run:\s*(.*)$")
CALL_RE = re.compile(r"^\s*(?:-\s+)?uses:\s*\./\.github/workflows/([^\s@]+)\s*$")
# Directories whose test files a workflow runs as a whole: directory -> (workflow file, command).
DIRECTORY_SUITES = {"python/tests": ("python.yml", re.compile(r"\bpytest tests(?:\s|$)"))}


def test_files(root: Path) -> Iterable[str]:
    """The repository-relative paths (with `/`) of the `test_*.py` and `test_*.js` files under `root`."""
    for directory, subdirs, files in os.walk(root):
        subdirs[:] = sorted(d for d in subdirs if d not in SKIPPED_DIRS)
        for name in sorted(files):
            if TEST_FILE_RE.match(name):
                yield (Path(directory) / name).relative_to(root).as_posix()


def indent(line: str) -> int:
    return len(line) - len(line.lstrip(" "))


def triggers(lines: list[str]) -> set[str]:
    """The events of the top-level `on:` key (`on: push`, `on: [push, schedule]`, or a block)."""
    for i, line in enumerate(lines):
        match = ON_RE.match(line)
        if not match:
            continue
        inline = match.group(1).strip()
        if inline:
            return {word for word in re.split(r"[\s,\[\]]+", inline) if word}
        events = set()
        child = None
        for below in lines[i + 1:]:
            if not below.strip():
                continue
            if indent(below) == 0:
                break
            child = indent(below) if child is None else child
            if indent(below) == child:
                events.add(below.strip().lstrip("-").strip().split(":")[0].strip())
        return events
    return set()


def run_text(lines: list[str]) -> str:
    """The values of every `run:` key, inline or block (`|`, `>`), one after another."""
    out = []
    i = 0
    while i < len(lines):
        match = RUN_RE.match(lines[i])
        i += 1
        if not match:
            continue
        value = match.group(2).strip()
        if value[:1] not in ("|", ">"):
            out.append(value)
            continue
        key_column = len(match.group(1))
        while i < len(lines) and (not lines[i].strip() or indent(lines[i]) > key_column):
            out.append(lines[i].strip())
            i += 1
    return "\n".join(out)


def automatic_run_texts(root: Path) -> dict[str, str]:
    """File name -> the `run:` text, for the workflows that run without a person."""
    directory = root / ".github" / "workflows"
    lines = {}
    for path in sorted(directory.glob("*.y*ml")):
        text = path.read_text(encoding="utf-8")
        lines[path.name] = [YAML_COMMENT_RE.sub("", line) for line in text.splitlines()]
    if not lines:
        raise ValueError(f"{directory}: no workflow files")
    running = [name for name, body in lines.items() if triggers(body) & AUTOMATIC_EVENTS]
    seen = set(running)
    while running:
        name = running.pop()
        for line in lines[name]:
            match = CALL_RE.match(line)
            if match and match.group(1) in lines and match.group(1) not in seen:
                seen.add(match.group(1))
                running.append(match.group(1))
    return {name: run_text(lines[name]) for name in sorted(seen)}


def is_run(path: str, workflows: dict[str, str]) -> bool:
    """Whether a running workflow names the file, or runs the directory suite that holds it."""
    for directory, (workflow, command) in DIRECTORY_SUITES.items():
        if path.startswith(directory + "/") and command.search(workflows.get(workflow, "")):
            return True
    name = path.rsplit("/", 1)[-1]
    return any(name in text for text in workflows.values())


def not_run(root: Path) -> tuple[list[str], list[str]]:
    """(files, the files no running workflow runs)."""
    workflows = automatic_run_texts(root)
    files = list(test_files(root))
    return files, [path for path in files if not is_run(path, workflows)]


def parse_args(argv: Sequence[str] | None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--root", type=Path, default=Path("."), help="repository root")
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    root = parse_args(argv).root.resolve()
    try:
        files, missing = not_run(root)
    except (OSError, ValueError) as error:
        print(f"tests-in-CI error: {error}", file=sys.stderr)
        return 2
    for path in missing:
        print(
            f"NOT RUN\t{path}\tno workflow that runs on push, pull request or schedule runs it: add a "
            "step (ci.yml `rust lint` for a quick test of a tool, wasm-threading.yml for one that needs the Wasm artifacts), or delete the test with the tool "
            "it tests",
        )
    print(f"tests in CI\tfiles={len(files)}\tnot_run={len(missing)}")
    return 0 if not missing else 1


if __name__ == "__main__":
    raise SystemExit(main())
