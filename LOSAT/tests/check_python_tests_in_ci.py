#!/usr/bin/env python3
"""Fail when a Python test file of the repository is run by no workflow.

No NCBI counterpart: this checks LOSAT's own CI wiring, not ported behavior.

The SWEBOK audit of 2026-11-04 (KA5) found five `test_*.py` files with 35 tests
that no workflow ran (`test_audit_blastp_v010.py`, `test_audit_tblastx_v010.py`,
`test_benchmark_wasm_threading.py`, `test_certify_blastn_v010.py`,
`test_compare_blastn_parity.py`), so a change to the tools they test could break
them unseen. Every workflow names the test files it runs one by one
(`unittest discover -s LOSAT/tests -p 'test_x.py'`, `unittest LOSAT/tests/test_x.py`),
so this check fails when the file name of a `test_*.py` file appears in no
`.github/workflows/*.yml` outside a comment. A file name shared by two
directories counts as named for both.

`python/tests/` is pylosat's pytest suite, which `python.yml` runs as a directory
(`pytest tests`, working directory `python`): its files are run while that
workflow holds that command; without it, each file needs its name in a workflow
like any other.

Usage: python3 LOSAT/tests/check_python_tests_in_ci.py --root .
"""

from __future__ import annotations

import argparse
import os
import re
import sys
from pathlib import Path
from typing import Iterable, Sequence

SKIPPED_DIRS = frozenset({".git", ".worktrees", "target", "node_modules", ".venv", "__pycache__"})
TEST_FILE_RE = re.compile(r"^test_.*\.py$")
# A YAML comment: `#` at the start of a line or after a space (`'#'` inside a value is not one).
YAML_COMMENT_RE = re.compile(r"(^|\s)#.*$", re.M)
# Directories whose test files a workflow runs as a whole: directory -> (workflow file, command).
DIRECTORY_SUITES = {"python/tests": ("python.yml", "pytest tests")}


def test_files(root: Path) -> Iterable[str]:
    """The repository-relative paths (with `/`) of the `test_*.py` files under `root`."""
    for directory, subdirs, files in os.walk(root):
        subdirs[:] = sorted(d for d in subdirs if d not in SKIPPED_DIRS)
        for name in sorted(files):
            if TEST_FILE_RE.match(name):
                yield (Path(directory) / name).relative_to(root).as_posix()


def workflow_texts(root: Path) -> dict[str, str]:
    """File name -> text without comments, for every workflow file."""
    directory = root / ".github" / "workflows"
    texts = {}
    for path in sorted(directory.glob("*.y*ml")):
        texts[path.name] = YAML_COMMENT_RE.sub("", path.read_text(encoding="utf-8"))
    if not texts:
        raise ValueError(f"{directory}: no workflow files")
    return texts


def is_run(path: str, workflows: dict[str, str]) -> bool:
    """Whether a workflow names the file, or runs the directory suite that holds it."""
    for directory, (workflow, command) in DIRECTORY_SUITES.items():
        if path.startswith(directory + "/") and command in workflows.get(workflow, ""):
            return True
    name = path.rsplit("/", 1)[-1]
    return any(name in text for text in workflows.values())


def not_run(root: Path) -> tuple[list[str], list[str]]:
    """(files, the files no workflow runs)."""
    workflows = workflow_texts(root)
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
        print(f"python-tests-in-CI error: {error}", file=sys.stderr)
        return 2
    for path in missing:
        print(
            f"NOT RUN\t{path}\tno workflow names it: add a step that runs it (ci.yml `rust lint` "
            "for a quick test of a tool), or delete the test with the tool it tests",
        )
    print(f"python tests in CI\tfiles={len(files)}\tnot_run={len(missing)}")
    return 0 if not missing else 1


if __name__ == "__main__":
    raise SystemExit(main())
