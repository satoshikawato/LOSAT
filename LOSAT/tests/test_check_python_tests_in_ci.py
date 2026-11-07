#!/usr/bin/env python3
"""Tests for the Python-tests-in-CI check (`check_python_tests_in_ci.py`).

No NCBI counterpart: the checker reads LOSAT's own workflow files.
"""

from __future__ import annotations

import contextlib
import importlib.util
import io
import sys
import tempfile
import unittest
from pathlib import Path


CHECKER_PATH = Path(__file__).with_name("check_python_tests_in_ci.py")
SPEC = importlib.util.spec_from_file_location("python_tests_in_ci", CHECKER_PATH)
assert SPEC is not None and SPEC.loader is not None
checker = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = checker
SPEC.loader.exec_module(checker)

CI = """jobs:
  lint:
    steps:
      - name: Tool tests
        run: python3 -m unittest discover -s LOSAT/tests -p 'test_tool.py'
"""
PYTHON_YML = """jobs:
  pylosat:
    defaults:
      run:
        working-directory: python
    steps:
      - run: .venv/bin/python -m pytest tests --timeout=300
"""


def write(root: Path, rel: str, text: str = "") -> None:
    path = root / rel
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")


def run(root: Path) -> tuple[int, str]:
    out = io.StringIO()
    with contextlib.redirect_stdout(out), contextlib.redirect_stderr(out):
        code = checker.main(["--root", str(root)])
    return code, out.getvalue()


class PythonTestsInCiTests(unittest.TestCase):
    def test_a_named_file_passes_and_an_unnamed_one_fails(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            write(root, ".github/workflows/ci.yml", CI)
            write(root, "LOSAT/tests/test_tool.py")
            write(root, "LOSAT/tests/tool.py")
            code, out = run(root)
            self.assertEqual(code, 0, out)
            self.assertIn("files=1\tnot_run=0", out)
            write(root, "LOSAT/tests/test_other.py")
            write(root, "scripts/test_script.py")
            code, out = run(root)
            self.assertEqual(code, 1)
            self.assertIn("NOT RUN\tLOSAT/tests/test_other.py\t", out)
            self.assertIn("NOT RUN\tscripts/test_script.py\t", out)
            self.assertNotIn("test_tool.py\t", out)

    def test_a_name_in_a_comment_does_not_count(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            write(
                root,
                ".github/workflows/ci.yml",
                CI + "      # - run: python3 -m unittest LOSAT/tests/test_old.py\n"
                "      - run: echo hi  # test_trailing.py\n",
            )
            write(root, "LOSAT/tests/test_old.py")
            write(root, "LOSAT/tests/test_trailing.py")
            code, out = run(root)
            self.assertEqual(code, 1)
            self.assertIn("NOT RUN\tLOSAT/tests/test_old.py", out)
            self.assertIn("NOT RUN\tLOSAT/tests/test_trailing.py", out)

    def test_a_directory_suite_counts_while_its_workflow_runs_it(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            write(root, ".github/workflows/ci.yml", CI)
            write(root, ".github/workflows/python.yml", PYTHON_YML)
            write(root, "python/tests/test_search.py")
            self.assertEqual(run(root)[0], 0)
            # The command in another workflow does not run python/tests.
            write(root, ".github/workflows/python.yml", "jobs: {}\n")
            write(root, ".github/workflows/ci.yml", CI + "      - run: pytest tests\n")
            code, out = run(root)
            self.assertEqual(code, 1)
            self.assertIn("NOT RUN\tpython/tests/test_search.py", out)
            # Without the directory command, a step that names the file runs it.
            write(root, ".github/workflows/ci.yml", CI + "      - run: pytest python/tests/test_search.py\n")
            self.assertEqual(run(root)[0], 0)

    def test_build_output_worktrees_and_packages_are_not_read(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            write(root, ".github/workflows/ci.yml", CI)
            write(root, "LOSAT/target/debug/build/x/test_gen.py")
            write(root, ".worktrees/other/LOSAT/tests/test_new.py")
            write(root, "web/app/node_modules/pkg/test_pkg.py")
            write(root, "python/.venv/lib/site-packages/test_dep.py")
            self.assertEqual(run(root)[0], 0)

    def test_no_workflow_files_is_an_error(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            write(root, "LOSAT/tests/test_tool.py")
            code, out = run(root)
            self.assertEqual(code, 2)
            self.assertIn("no workflow files", out)


if __name__ == "__main__":
    unittest.main()
