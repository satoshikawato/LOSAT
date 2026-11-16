#!/usr/bin/env python3
"""Tests for the lint-allow ratchet (`check_lint_allow_ratchet.py`).

No NCBI counterpart: the checker reads LOSAT's own lint attributes.
"""

from __future__ import annotations

import contextlib
import importlib.util
import io
import sys
import tempfile
import unittest
from pathlib import Path


CHECKER_PATH = Path(__file__).with_name("check_lint_allow_ratchet.py")
SPEC = importlib.util.spec_from_file_location("lint_allow_ratchet", CHECKER_PATH)
assert SPEC is not None and SPEC.loader is not None
ratchet = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = ratchet
SPEC.loader.exec_module(ratchet)

BLANKET = "#![allow(warnings, clippy::all)]\n\nfn main() {}\n"


def write(root: Path, rel: str, text: str) -> None:
    path = root / rel
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")


def run(root: Path, listed: list[str]) -> tuple[int, str]:
    write(root, "LOSAT/tests/lint_allow_warnings.tsv", "# comment\npath\n" + "".join(f"{p}\n" for p in listed))
    out = io.StringIO()
    with contextlib.redirect_stdout(out), contextlib.redirect_stderr(out):
        code = ratchet.main(["--root", str(root)])
    return code, out.getvalue()


class LintAllowRatchetTests(unittest.TestCase):
    def test_the_attribute_forms_it_counts(self) -> None:
        self.assertTrue(ratchet.has_allow_warnings("#![allow(warnings)]\n"))
        self.assertTrue(ratchet.has_allow_warnings("  #![allow(dead_code, warnings)]\n"))
        self.assertTrue(ratchet.has_allow_warnings("#![cfg_attr(test, allow(warnings))]\n"))
        self.assertTrue(ratchet.has_allow_warnings("#![allow(\n    dead_code,\n    warnings\n)]\n"))
        self.assertTrue(ratchet.has_allow_warnings("#[allow(warnings)]\nmod x;\n"))
        self.assertTrue(ratchet.has_allow_warnings("#![allow(clippy::all)]\n"))
        self.assertTrue(ratchet.has_allow_warnings("#![expect(warnings)]\n"))
        self.assertFalse(ratchet.has_allow_warnings("#![allow(dead_code)]\n"))
        self.assertFalse(ratchet.has_allow_warnings("#![allow(clippy::too_many_arguments)]\n"))
        self.assertFalse(ratchet.has_allow_warnings("// #![allow(warnings)]\n"))
        self.assertFalse(ratchet.has_allow_warnings("/* #![allow(warnings)] */\nfn f() {}\n"))

    def test_listed_files_pass_and_a_new_one_fails(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            write(root, "LOSAT/src/lib.rs", BLANKET)
            write(root, "LOSAT/src/clean.rs", "fn f() {}\n")
            self.assertEqual(run(root, ["LOSAT/src/lib.rs"])[0], 0)
            write(root, "LOSAT/tests/new_test.rs", BLANKET)
            code, out = run(root, ["LOSAT/src/lib.rs"])
            self.assertEqual(code, 1)
            self.assertIn("UNEXPECTED\tLOSAT/tests/new_test.rs", out)

    def test_a_listed_file_without_the_attribute_must_leave_the_list(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            write(root, "LOSAT/src/lib.rs", "fn f() {}\n")
            code, out = run(root, ["LOSAT/src/lib.rs"])
            self.assertEqual(code, 1)
            self.assertIn("STALE\tLOSAT/src/lib.rs", out)

    def test_build_output_and_worktrees_are_not_read(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            write(root, "LOSAT/target/debug/build/x/out.rs", BLANKET)
            write(root, ".worktrees/other/LOSAT/src/lib.rs", BLANKET)
            self.assertEqual(run(root, [])[0], 0)


if __name__ == "__main__":
    unittest.main()
