#!/usr/bin/env python3
"""Tests for the tests-in-CI check (`check_tests_in_ci.py`).

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


CHECKER_PATH = Path(__file__).with_name("check_tests_in_ci.py")
SPEC = importlib.util.spec_from_file_location("tests_in_ci", CHECKER_PATH)
assert SPEC is not None and SPEC.loader is not None
checker = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = checker
SPEC.loader.exec_module(checker)

CI = """on:
  push:
    branches: [main]
  pull_request:
jobs:
  lint:
    steps:
      - name: Tool tests
        run: python3 -m unittest discover -s LOSAT/tests -p 'test_tool.py'
"""
PYTHON_YML = """on:
  pull_request:
    paths:
      - "python/**"
jobs:
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


def workflow(on: str, run_lines: str) -> str:
    return f"on: {on}\njobs:\n  a:\n    steps:\n{run_lines}"


class TestsInCiTests(unittest.TestCase):
    def test_trigger_forms(self) -> None:
        self.assertEqual(checker.triggers(["on: push"]), {"push"})
        self.assertEqual(checker.triggers(["on: [push, schedule]"]), {"push", "schedule"})
        block = ["name: x", "on:", "  workflow_call:", "    inputs:", "      push:", "  schedule:",
                 "    - cron: '1 2 * * *'", "jobs:", "  push:"]
        self.assertEqual(checker.triggers(block), {"workflow_call", "schedule"})
        self.assertEqual(checker.triggers(["on:", "  - push", "  - pull_request"]), {"push", "pull_request"})
        self.assertEqual(checker.triggers(["jobs:", "  a:", "    on: push"]), set())

    def test_run_text_reads_inline_and_block_values_only(self) -> None:
        lines = [
            "    paths: ['tests/test_paths.py']",
            "      - run: python3 tests/test_inline.py",
            "      - name: block",
            "        run: |",
            "          python3 tests/test_block_1.py",
            "",
            "          python3 tests/test_block_2.py",
            "        env:",
            "          X: tests/test_env.py",
            "        run: >-",
            "          python3 tests/test_folded.py",
        ]
        text = checker.run_text(lines)
        for name in ("test_inline.py", "test_block_1.py", "test_block_2.py", "test_folded.py"):
            self.assertIn(name, text)
        for name in ("test_paths.py", "test_env.py"):
            self.assertNotIn(name, text)

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

    def test_a_node_test_file_counts_when_a_run_command_names_it(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            write(root, ".github/workflows/ci.yml", CI + "      - run: node LOSAT/tests/test_node_tool.js x.wasm\n")
            write(root, "LOSAT/tests/test_tool.py")
            write(root, "LOSAT/tests/test_node_tool.js")
            code, out = run(root)
            self.assertEqual(code, 0, out)
            self.assertIn("files=2\tnot_run=0", out)
            write(root, "LOSAT/tests/test_unnamed.js")
            code, out = run(root)
            self.assertEqual(code, 1)
            self.assertIn("NOT RUN\tLOSAT/tests/test_unnamed.js\t", out)
            self.assertNotIn("test_node_tool.js\t", out)

    def test_other_js_and_ts_names_are_not_test_files_for_this_check(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            write(root, ".github/workflows/ci.yml", CI)
            write(root, "LOSAT/tests/test_tool.py")
            for rel in ("LOSAT/tests/tool.test.js", "web/app/src/x.test.ts", "web/app/src/test_x.ts",
                        "LOSAT/tests/helper_test.js", "LOSAT/tests/test_data.json"):
                write(root, rel)
            code, out = run(root)
            self.assertEqual(code, 0, out)
            self.assertIn("files=1\tnot_run=0", out)

    def test_a_name_outside_a_run_command_does_not_count(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            write(
                root,
                ".github/workflows/ci.yml",
                CI.replace("  pull_request:\n", "  pull_request:\n    paths: ['LOSAT/tests/test_paths.py']\n")
                + "      # - run: python3 -m unittest LOSAT/tests/test_old.py\n"
                "      - run: echo hi  # test_trailing.py\n"
                "        env:\n          FILE: LOSAT/tests/test_env.py\n",
            )
            for name in ("test_tool", "test_paths", "test_old", "test_trailing", "test_env"):
                write(root, f"LOSAT/tests/{name}.py")
            code, out = run(root)
            self.assertEqual(code, 1)
            for name in ("test_paths", "test_old", "test_trailing", "test_env"):
                self.assertIn(f"NOT RUN\tLOSAT/tests/{name}.py", out)
            self.assertIn("not_run=4", out)

    def test_only_workflows_that_run_without_a_person_count(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            step = "      - run: python3 -m unittest discover -s LOSAT/tests -p 'test_{}.py'\n"
            write(root, ".github/workflows/ci.yml", CI)
            write(root, ".github/workflows/release.yml",
                  workflow("workflow_dispatch", step.format("release")
                           + "  b:\n    uses: ./.github/workflows/called_by_release.yml\n"))
            write(root, ".github/workflows/called_by_release.yml",
                  workflow("[workflow_call, workflow_dispatch]", step.format("called_by_release")))
            write(root, ".github/workflows/nightly.yml",
                  workflow("\n  schedule:\n    - cron: '17 18 * * *'\n  workflow_dispatch:", step.format("nightly"))
                  + "  wasm:\n    uses: ./.github/workflows/called.yml\n")
            write(root, ".github/workflows/called.yml", workflow("workflow_call", step.format("called")))
            for name in ("tool", "release", "called_by_release", "nightly", "called"):
                write(root, f"LOSAT/tests/test_{name}.py")
            code, out = run(root)
            self.assertEqual(code, 1)
            self.assertIn("NOT RUN\tLOSAT/tests/test_release.py", out)
            self.assertIn("NOT RUN\tLOSAT/tests/test_called_by_release.py", out)
            self.assertIn("not_run=2", out)

    def test_a_directory_suite_counts_while_its_workflow_runs_it(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            write(root, ".github/workflows/ci.yml", CI)
            write(root, ".github/workflows/python.yml", PYTHON_YML)
            write(root, "python/tests/test_search.py")
            write(root, "python/tests/test_other.py")
            self.assertEqual(run(root)[0], 0)
            # One file of the directory is not the directory.
            write(root, ".github/workflows/python.yml", PYTHON_YML.replace("pytest tests", "pytest tests/test_search.py"))
            code, out = run(root)
            self.assertEqual(code, 1)
            self.assertIn("NOT RUN\tpython/tests/test_other.py", out)
            self.assertIn("not_run=1", out)
            # The command in another workflow does not run python/tests.
            write(root, ".github/workflows/python.yml", "on: pull_request\njobs: {}\n")
            write(root, ".github/workflows/ci.yml", CI + "      - run: pytest tests\n")
            code, out = run(root)
            self.assertEqual(code, 1)
            self.assertIn("NOT RUN\tpython/tests/test_search.py", out)
            # Without the directory command, a step that names the file runs it.
            write(root, ".github/workflows/ci.yml",
                  CI + "      - run: pytest python/tests/test_search.py python/tests/test_other.py\n")
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
