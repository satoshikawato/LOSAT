#!/usr/bin/env python3
"""Tests for the check of quoted foreign source (`check_foreign_quotes.py`).

No NCBI counterpart: the checker reads LOSAT's own comments.
"""

from __future__ import annotations

import contextlib
import importlib.util
import io
import sys
import tempfile
import unittest
from pathlib import Path


CHECKER_PATH = Path(__file__).with_name("check_foreign_quotes.py")
SPEC = importlib.util.spec_from_file_location("foreign_quotes", CHECKER_PATH)
assert SPEC is not None and SPEC.loader is not None
checker = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = checker
SPEC.loader.exec_module(checker)


def findings(text: str) -> list[tuple[str, int, str]]:
    return [(f.kind, f.line, f.source) for f in checker.check_text("x.rs", text.splitlines())]


NCBI_QUOTE = """\
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:1379-1381
// ```c
//     qsort(hsp_list->hsp_array, hsp_list->hspcnt, sizeof(BlastHSP*), ScoreCompareHSPs);
// ```
fn f() {}
"""


class CheckTextTest(unittest.TestCase):
    def test_ncbi_quote_passes(self) -> None:
        self.assertEqual(findings(NCBI_QUOTE), [])

    def test_quote_after_another_projects_label_fails(self) -> None:
        text = (
            "// Microsoft.Windows.SDK.CRTSource 10.0.22621.3 (MIT): ucrt/time/timeset.cpp:67-69\n"
            "// ```c++\n//     _timezone.initialize(8 * 3600L);\n// ```\nconst X: i32 = 0;\n"
        )
        self.assertEqual(findings(text), [("QUOTE", 2, "Microsoft C runtime")])

    def test_quote_after_an_unknown_projects_label_fails(self) -> None:
        text = "// Foo 1.2: src/foo.c:10-12\n// ```c\n// int x;\n// ```\n"
        self.assertEqual(findings(text), [("QUOTE", 2, "Foo 1.2")])

    def test_quote_after_a_prose_citation_of_libstdcxx_fails(self) -> None:
        text = (
            "    // libstdc++ reference (bits/stl_algobase.h, `std::__lower_bound`)\n"
            "    // ```c++\n    //   __len = std::distance(__first, __last);\n    // ```\n"
        )
        self.assertEqual(findings(text), [("QUOTE", 2, "libstdc++")])

    def test_ncbi_label_of_glibc_fails_twice(self) -> None:
        text = (
            "// NCBI reference (glibc 2.17 to 2.39): libio/fileops.c `_IO_new_file_overflow`\n"
            "// ```c\n//   return EOF;\n// ```\n"
        )
        self.assertEqual(findings(text), [("LABEL", 1, "glibc"), ("QUOTE", 2, "glibc")])

    def test_ncbi_label_without_a_path_fails(self) -> None:
        text = "// NCBI reference: libstdc++'s `basic_filebuf::xsgetn` (GCC 14)\n"
        self.assertEqual(findings(text), [("LABEL", 1, "libstdc++")])

    def test_copies_in_ncbis_tree_are_not_ncbis(self) -> None:
        lmdb = "c++/src/util/lmdb/mdb.c:6540-7066"
        self.assertEqual(
            findings(f"// NCBI reference (598d8ae6): {lmdb}\nfn f() {{}}\n"), [("LABEL", 1, "LMDB")]
        )
        copy = f"// LMDB 0.9.24 (NCBI 598d8ae6 copy): {lmdb}\n"
        self.assertEqual(findings(copy + "fn f() {}\n"), [])
        quoted = copy + "// ```c\n// rc = 0;\n// ```\n"
        self.assertEqual(findings(quoted), [("QUOTE", 2, "LMDB")])
        boost = "// Boost (NCBI 598d8ae6 copy): c++/src/algo/blast/core/boost_erf.c:252-255\n"
        self.assertEqual(findings(boost + "// ```c\n// x\n// ```\n"), [("QUOTE", 2, "Boost")])

    def test_a_citation_over_two_lines_is_read_from_its_label(self) -> None:
        text = (
            "// NCBI reference (598d8ae6): /usr/include/c++/8/bits/stl_heap.h (__make_heap,\n"
            "// __push_heap), and stl_algo.h __heap_select\n"
            "// ```c++\n// __heap_select(__first, __middle, __last, __comp) {\n// ```\n"
        )
        self.assertEqual(findings(text), [("LABEL", 1, "libstdc++"), ("QUOTE", 3, "libstdc++")])

    def test_a_quote_inside_a_function_takes_the_citation_above_the_function(self) -> None:
        text = (
            "    // LMDB 0.9.24 (NCBI 598d8ae6 copy): c++/src/util/lmdb/mdb.c:6540-7066\n"
            "    fn cursor_put(&mut self) -> i32 {\n"
            "        let rc = 0;\n"
            "        // ```c\n        // \tnospill = flags & MDB_NOSPILL;\n        // ```\n"
            "        rc\n    }\n"
        )
        self.assertEqual(findings(text), [("QUOTE", 4, "LMDB")])

    def test_the_last_label_above_a_quote_decides(self) -> None:
        text = (
            "// libstdc++ (GCC 13): bits/list.tcc:614-662\n" + NCBI_QUOTE.replace("fn f() {}\n", "")
        )
        self.assertEqual(findings(text), [])
        prose = (
            "/// libstdc++'s insertion sort moves the second one (glibc's qsort, GCC, LLVM).\n"
            "/// NCBI reference: c++/src/algo/blast/api/setup_factory.cpp:170-184\n"
            "/// ```c\n///     Blast_Message2TSearchMessages(blast_msg.Get(), query_info, msgs);\n"
            "/// ```\n"
        )
        self.assertEqual(findings(prose), [])

    def test_a_fence_without_a_citation_is_not_checked(self) -> None:
        text = "/// Example:\n/// ```text\n/// losat blastn -query q.fa\n/// ```\nfn f() {}\n"
        self.assertEqual(findings(text), [])

    def test_ncbis_prose_naming_glibc_is_not_a_label(self) -> None:
        text = "// NCBI's `CArg_Double` reads glibc's hexadecimal numbers and a final `e`\n"
        self.assertEqual(findings(text), [])


class MainTest(unittest.TestCase):
    def run_main(self, files: dict[str, str]) -> tuple[int, str]:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            for rel, text in files.items():
                path = root / rel
                path.parent.mkdir(parents=True, exist_ok=True)
                path.write_text(text, encoding="utf-8")
            out = io.StringIO()
            with contextlib.redirect_stdout(out):
                code = checker.main(["--root", str(root)])
            return code, out.getvalue()

    def test_clean_tree_passes(self) -> None:
        code, out = self.run_main({"LOSAT/src/a.rs": NCBI_QUOTE})
        self.assertEqual(code, 0)
        self.assertIn("findings=0", out)

    def test_foreign_quote_fails_with_its_place(self) -> None:
        bad = "// SQLite 3.51.1: src/btree.c:404-411\n// ```c\n// x\n// ```\n"
        code, out = self.run_main({"LOSAT/src/a.rs": NCBI_QUOTE, "LOSAT/src/b.rs": bad})
        self.assertEqual(code, 1)
        self.assertIn("QUOTE\tLOSAT/src/b.rs:2\tSQLite", out)
        self.assertIn("findings=1", out)

    def test_build_output_is_skipped(self) -> None:
        bad = "// UCRT reference (10.0.26100.0): ucrt/time/tzset.cpp:1-2\n// ```c++\n// x\n// ```\n"
        code, _ = self.run_main({"target/debug/build/a.rs": bad, ".worktrees/w/b.rs": bad})
        self.assertEqual(code, 0)


if __name__ == "__main__":
    unittest.main()
