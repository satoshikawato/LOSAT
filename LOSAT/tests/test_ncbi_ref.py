#!/usr/bin/env python3
"""Tests for `scripts/ncbi_ref.py`, which prints the NCBI lines that LOSAT's comments cite.

No NCBI counterpart: a reading aid for LOSAT's comments.
"""

from __future__ import annotations

import contextlib
import importlib.util
import io
import sys
import tempfile
import unittest
from pathlib import Path

SCRIPT = Path(__file__).resolve().parents[2] / "scripts" / "ncbi_ref.py"
SPEC = importlib.util.spec_from_file_location("ncbi_ref", SCRIPT)
assert SPEC is not None and SPEC.loader is not None
ncbi_ref = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = ncbi_ref
SPEC.loader.exec_module(ncbi_ref)

RUST = """\
// NCBI reference (598d8ae6): c++/src/a.c:2-3
// The second and third lines.
fn f() {
    let x = 1;
}

/// NCBI reference (598d8ae6): c++/src/a.c:5, c++/src/b.c:1
fn g() {}
"""


class NcbiRefTest(unittest.TestCase):
    def setUp(self) -> None:
        self.tmp = tempfile.TemporaryDirectory()
        root = Path(self.tmp.name)
        self.ncbi = root / "ncbi"
        (self.ncbi / "c++/src").mkdir(parents=True)
        (self.ncbi / "c++/src/a.c").write_text("one\ntwo\nthree\nfour\nfive\n", encoding="utf-8")
        (self.ncbi / "c++/src/b.c").write_text("bee\n", encoding="utf-8")
        self.rust = root / "x.rs"
        self.rust.write_text(RUST, encoding="utf-8")

    def tearDown(self) -> None:
        self.tmp.cleanup()

    def run_main(self, *args: str) -> tuple[int, str]:
        out = io.StringIO()
        with contextlib.redirect_stdout(out), contextlib.redirect_stderr(io.StringIO()):
            code = ncbi_ref.main([*args, "--ncbi-src", str(self.ncbi)])
        return code, out.getvalue()

    def test_a_line_takes_the_citations_of_the_comment_above_it(self) -> None:
        code, out = self.run_main(f"{self.rust}:4")
        self.assertEqual(code, 0)
        self.assertEqual(out, "== c++/src/a.c:2-3\n     2  two\n     3  three\n")

    def test_a_line_inside_a_comment_takes_that_comment(self) -> None:
        code, out = self.run_main(f"{self.rust}:7")
        self.assertEqual(code, 0)
        self.assertIn("== c++/src/a.c:5\n     5  five\n", out)
        self.assertIn("== c++/src/b.c:1\n     1  bee\n", out)

    def test_a_file_prints_every_citation_and_context_marks_the_lines_around(self) -> None:
        code, out = self.run_main(str(self.rust))
        self.assertEqual(code, 0)
        self.assertEqual(out.count("== "), 3)
        code, out = self.run_main("c++/src/a.c:3", "--context", "1")
        self.assertEqual(out, "== c++/src/a.c:3\n     2. two\n     3  three\n     4. four\n")

    def test_rpsbproc_citations_read_their_own_tree(self) -> None:
        rpsb = Path(self.tmp.name) / "rpsb"
        (rpsb / "RpsbProc/src").mkdir(parents=True)
        (rpsb / "RpsbProc/src/main.cpp").write_text("m1\nm2\n", encoding="utf-8")
        out = io.StringIO()
        with contextlib.redirect_stdout(out), contextlib.redirect_stderr(io.StringIO()):
            code = ncbi_ref.main(["RpsbProc/src/main.cpp:2", "--rpsbproc-src", str(rpsb)])
        self.assertEqual((code, out.getvalue()), (0, "== RpsbProc/src/main.cpp:2\n     2  m2\n"))
        code, _ = self.run_main("RpsbProc/src/main.cpp:2")
        self.assertEqual(code, 1)

    def test_lines_outside_the_file_fail(self) -> None:
        code, out = self.run_main("c++/src/a.c:5-6")
        self.assertEqual(code, 1)
        self.assertIn("outside the file's 5 lines", out)
        code, _ = self.run_main("c++/src/none.c:1")
        self.assertEqual(code, 1)


if __name__ == "__main__":
    unittest.main()
